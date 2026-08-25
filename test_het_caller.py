"""Tests for the Sanger heterozygote peak-ratio calling algorithm.

Run with: python test_het_caller.py
"""

import unittest

import numpy as np

from het_caller.core import (
    TraceData,
    call_peaks,
    iupac_code,
    map_flagged_to_reference,
    parse_ab1,
    trace_quality_summary,
)


class FakeRecord:
    """Minimal stand-in for a Bio.SeqIO 'abi' record, exercising parse_ab1
    without needing a real binary .ab1 file on disk."""

    def __init__(self, abif_raw, sample_id="sample1"):
        self.annotations = {"abif_raw": abif_raw}
        self.id = sample_id


def make_trace(sequence, primary_heights, secondary=None, n_points=None, spacing=10, baseline=20.0):
    """Build a synthetic TraceData for a given base-call string.

    `primary_heights`: list of peak heights for the called base at each position.
    `secondary`: optional dict {position: (other_base, height)} injecting a
    second, overlapping peak at that position (simulating a het/mosaic call).
    """
    n = len(sequence)
    peak_locations = [spacing * (i + 1) for i in range(n)]
    if n_points is None:
        n_points = (peak_locations[-1] + spacing) if peak_locations else spacing

    channels = {b: np.full(n_points, baseline, dtype=np.float64) for b in "ACGT"}
    for i, (base, loc, height) in enumerate(zip(sequence, peak_locations, primary_heights)):
        channels[base][loc] = height

    secondary = secondary or {}
    for pos, (other_base, height) in secondary.items():
        channels[other_base][peak_locations[pos]] = height

    return TraceData(
        channels=channels,
        peak_locations=peak_locations,
        called_bases=sequence,
        qualities=[40] * n,
        sample_name="synthetic",
    )


class TestIupacCode(unittest.TestCase):
    def test_single_base(self):
        self.assertEqual(iupac_code(["A"]), "A")

    def test_pair_codes(self):
        self.assertEqual(iupac_code(["A", "G"]), "R")
        self.assertEqual(iupac_code(["C", "T"]), "Y")
        self.assertEqual(iupac_code(["G", "C"]), "S")
        # order-independence
        self.assertEqual(iupac_code(["T", "C"]), "Y")

    def test_unknown_combo_raises(self):
        with self.assertRaises(ValueError):
            iupac_code(["A", "Z"])


class TestCallPeaks(unittest.TestCase):
    def test_clean_homozygous_trace_no_flags(self):
        # Every position has one dominant peak far above baseline noise.
        trace = make_trace("ACGTAC", primary_heights=[500, 480, 510, 495, 505, 490])
        calls = call_peaks(trace, threshold=0.25, window=2)

        self.assertEqual(len(calls), 6)
        self.assertTrue(all(not c.flagged for c in calls))
        for c, base in zip(calls, "ACGTAC"):
            self.assertEqual(c.primary_base, base)
            self.assertEqual(c.iupac, base)

    def test_heterozygous_snp_is_flagged_with_correct_iupac(self):
        # Position 2 (0-based) gets a near-equal-height overlapping peak: A/G -> R
        seq = "ACGTAC"
        trace = make_trace(
            seq,
            primary_heights=[500, 480, 510, 495, 505, 490],
            secondary={2: ("A", 400)},  # G primary (510) + A secondary (400) at pos 2
        )
        calls = call_peaks(trace, threshold=0.25, window=2)

        het_call = calls[2]
        self.assertTrue(het_call.flagged)
        self.assertEqual(het_call.primary_base, "G")
        self.assertEqual(het_call.secondary_base, "A")
        self.assertAlmostEqual(het_call.ratio, 400 / 510, places=4)
        self.assertEqual(het_call.iupac, "R")

        # every other position should remain unflagged
        for i, c in enumerate(calls):
            if i != 2:
                self.assertFalse(c.flagged)

    def test_low_level_mosaicism_threshold_sensitivity(self):
        # Secondary peak at ~30% of primary height: caught by a lenient
        # threshold (mosaicism-sensitive) but not by a strict het-only one.
        seq = "ACGTAC"
        trace = make_trace(
            seq,
            primary_heights=[500, 480, 510, 495, 505, 490],
            secondary={4: ("T", 150)},  # A primary (505) + T secondary (150) -> ratio ~0.297
        )

        lenient_calls = call_peaks(trace, threshold=0.25, window=2)
        strict_calls = call_peaks(trace, threshold=0.5, window=2)

        self.assertTrue(lenient_calls[4].flagged)
        self.assertEqual(lenient_calls[4].iupac, iupac_code(["A", "T"]))
        self.assertFalse(strict_calls[4].flagged)
        self.assertEqual(strict_calls[4].iupac, "A")

    def test_invalid_threshold_raises(self):
        trace = make_trace("AC", primary_heights=[500, 480])
        with self.assertRaises(ValueError):
            call_peaks(trace, threshold=0.0)
        with self.assertRaises(ValueError):
            call_peaks(trace, threshold=1.5)


class TestTraceQualitySummary(unittest.TestCase):
    def test_summary_counts_and_snr(self):
        seq = "ACGTAC"
        trace = make_trace(
            seq,
            primary_heights=[500, 480, 510, 495, 505, 490],
            secondary={2: ("A", 400)},
        )
        calls = call_peaks(trace, threshold=0.25, window=2)
        summary = trace_quality_summary(trace, calls)

        self.assertEqual(summary["num_bases"], 6)
        self.assertEqual(summary["num_flagged"], 1)
        self.assertGreater(summary["signal_to_noise_ratio"], 1.0)
        self.assertAlmostEqual(
            summary["mean_primary_peak_height"],
            sum([500, 480, 510, 495, 505, 490]) / 6,
            places=1,
        )

    def test_summary_with_no_calls(self):
        trace = make_trace("", primary_heights=[])
        summary = trace_quality_summary(trace, [])
        self.assertEqual(summary["num_bases"], 0)
        self.assertEqual(summary["num_flagged"], 0)
        self.assertEqual(summary["mean_primary_peak_height"], 0.0)


class TestParseAb1(unittest.TestCase):
    def test_parse_ab1_extracts_channels_and_calls(self):
        n_points = 100
        data_a = [20.0] * n_points
        data_c = [20.0] * n_points
        data_g = [20.0] * n_points
        data_t = [20.0] * n_points
        data_a[10] = 500.0
        data_c[20] = 480.0

        raw = {
            "FWO_1": b"ACGT",
            "DATA9": data_a,
            "DATA10": data_c,
            "DATA11": data_g,
            "DATA12": data_t,
            "PLOC2": [10, 20],
            "PBAS2": b"AC",
            "PCON2": bytes([40, 38]),
        }
        record = FakeRecord(raw, sample_id="test_sample")

        trace = parse_ab1(record)

        self.assertEqual(trace.sample_name, "test_sample")
        self.assertEqual(trace.called_bases, "AC")
        self.assertEqual(trace.peak_locations, [10, 20])
        self.assertEqual(trace.qualities, [40, 38])
        self.assertEqual(trace.channels["A"][10], 500.0)
        self.assertEqual(trace.channels["C"][20], 480.0)
        self.assertEqual(set(trace.channels.keys()), {"A", "C", "G", "T"})

        calls = call_peaks(trace, threshold=0.25, window=1)
        self.assertEqual(len(calls), 2)
        self.assertFalse(calls[0].flagged)
        self.assertFalse(calls[1].flagged)

    def test_parse_ab1_falls_back_to_ploc1_pbas1(self):
        n_points = 50
        raw = {
            "FWO_1": "ACGT",
            "DATA9": [20.0] * n_points,
            "DATA10": [20.0] * n_points,
            "DATA11": [20.0] * n_points,
            "DATA12": [20.0] * n_points,
            "PLOC1": [5],
            "PBAS1": "G",
        }
        raw["DATA11"][5] = 600.0
        record = FakeRecord(raw)

        trace = parse_ab1(record)
        self.assertEqual(trace.called_bases, "G")
        self.assertEqual(trace.peak_locations, [5])
        self.assertEqual(trace.qualities, [0])  # no PCON tag present -> defaults to zeros


class TestReferenceMapping(unittest.TestCase):
    def test_map_flagged_to_reference_reports_position_and_context(self):
        seq = "ACGTAC"
        reference = "ACGTAC"  # identical, so alignment is 1:1
        trace = make_trace(
            seq,
            primary_heights=[500, 480, 510, 495, 505, 490],
            secondary={2: ("A", 400)},
        )
        calls = call_peaks(trace, threshold=0.25, window=2)

        mapping = map_flagged_to_reference(calls, reference, seq, context=2)

        self.assertEqual(len(mapping), 1)
        entry = mapping[0]
        self.assertEqual(entry["query_position"], 2)
        self.assertEqual(entry["reference_position"], 3)  # 1-based
        self.assertEqual(entry["iupac"], "R")
        self.assertIn("[G]", entry["reference_context"])

    def test_map_flagged_to_reference_no_flags_returns_empty(self):
        seq = "ACGTAC"
        trace = make_trace(seq, primary_heights=[500, 480, 510, 495, 505, 490])
        calls = call_peaks(trace, threshold=0.25, window=2)
        mapping = map_flagged_to_reference(calls, seq, seq, context=2)
        self.assertEqual(mapping, [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
