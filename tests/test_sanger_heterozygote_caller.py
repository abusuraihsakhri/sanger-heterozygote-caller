"""
Automated Pytest Test Suite for Sanger Heterozygote Caller.
Domain: Clinical & Biomedical AI
Standard: CAP / CLSI / ISO Standards
"""
import os
import sys
import warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from agents.base import PHIGuard, AuditLogger, SecurityException, AuditTrail
from agents.models import SystemTaskPayload, UrgencyLevel, SystemIntegrityStatus
from agents.workers import InvariantQCWorker, SafetyEscalationWorker, ProtocolConformanceWorker
from agents.supervisor import SystemSupervisor
from cli import main


def test_phi_guard_enforcement():
    with pytest.raises(SecurityException):
        PHIGuard.assert_no_phi("Patient MRN-994827 blood culture positive for Staphylococcus")

    # Clean text passes
    PHIGuard.assert_no_phi("Analytical assay specimen KEY-001 optimal")


def test_specialized_workers():
    # Worker 1: QC Invariant
    p1 = SystemTaskPayload(task_id="T1", target_identifier="KEY-01", primary_metric=35.0)
    alerts1 = InvariantQCWorker.evaluate(p1)
    assert len(alerts1) == 1
    assert alerts1[0].urgency == UrgencyLevel.ELEVATED

    # Worker 2: Safety
    p2 = SystemTaskPayload(task_id="T2", target_identifier="KEY-02", primary_metric=10.0, is_critical_flag=True)
    alerts2 = SafetyEscalationWorker.evaluate(p2)
    assert len(alerts2) == 1
    assert alerts2[0].urgency == UrgencyLevel.CRITICAL_STAT

    # Worker 3: Protocol Conformance
    p3 = SystemTaskPayload(task_id="T3", target_identifier="KEY-03", primary_metric=10.0, status_descriptor="DISCORDANT_ANOMALY")
    alerts3 = ProtocolConformanceWorker.evaluate(p3)
    assert len(alerts3) == 1


def test_supervisor_consensus_and_audit():
    supervisor = SystemSupervisor(model_provider="mock")
    payload = SystemTaskPayload(
        task_id="TASK-PROD-01",
        target_identifier="KEY-PROD-01",
        primary_metric=12.0,
        secondary_metric=4.0,
        status_descriptor="NOMINAL"
    )
    dossier = supervisor.process_task(payload)
    assert dossier.overall_urgency == UrgencyLevel.ROUTINE
    assert dossier.integrity_status == SystemIntegrityStatus.VALIDATED
    assert dossier.audit_hash != ""

    # Verify cryptographic audit trail
    assert AuditLogger.verify_integrity() is True

    # CLI tests
    assert main(["audit", "--task-id", "CLI-TEST-01"]) == 0
    assert main(["chat", "Explain", "specifications"]) == 0
    assert main(["verify-audit"]) == 0


# ---------------------------------------------------------------------------
# Security & reliability tests (added by improvement pass)
# ---------------------------------------------------------------------------

def test_audit_trail_uses_ephemeral_key_without_env():
    """When AUDIT_SECRET_KEY is unset, a random key is used (not a hardcoded default)."""
    saved = os.environ.pop("AUDIT_SECRET_KEY", None)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            trail = AuditTrail()
        assert len(trail.secret_key) > 0
        # Key should differ between instances (random)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            trail2 = AuditTrail()
        assert trail.secret_key != trail2.secret_key
    finally:
        if saved is not None:
            os.environ["AUDIT_SECRET_KEY"] = saved


def test_audit_trail_uses_env_key_when_set():
    """When AUDIT_SECRET_KEY is set, it is used directly."""
    os.environ["AUDIT_SECRET_KEY"] = "test-key-12345"
    try:
        trail = AuditTrail()
        assert trail.secret_key == b"test-key-12345"
    finally:
        os.environ.pop("AUDIT_SECRET_KEY", None)


def test_audit_trail_explicit_key_overrides_env():
    """An explicit secret_key argument takes precedence over the env var."""
    os.environ["AUDIT_SECRET_KEY"] = "env-key"
    try:
        trail = AuditTrail(secret_key="explicit-key")
        assert trail.secret_key == b"explicit-key"
    finally:
        os.environ.pop("AUDIT_SECRET_KEY", None)


def test_phi_guard_redacts_multiple_patterns():
    """PHIGuard.redact_phi should replace all PHI matches with a redaction marker."""
    text = "Patient John Doe, MRN-123456, SSN 123-45-6789, email test@example.com"
    redacted = PHIGuard.redact_phi(text)
    assert "John Doe" not in redacted
    assert "MRN" not in redacted
    assert "123-45-6789" not in redacted
    assert "test@example.com" not in redacted
    assert "[REDACTED_IDENTIFIER]" in redacted


def test_phi_guard_empty_and_none_safe():
    """Empty strings and None-like input should not raise."""
    PHIGuard.assert_no_phi("")
    PHIGuard.assert_no_phi("   ")
    PHIGuard.assert_no_phi("ACGTACGT")  # pure sequence
