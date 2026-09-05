# Sanger Heterozygote Caller

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Standards:** CAP / CLSI / ISO Quality Frameworks

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Sanger Heterozygote Caller flags heterozygous and low-level mosaic positions in Sanger `.ab1` chromatograms using a secondary-to-primary fluorescence peak-height ratio heuristic (similar to Mutation Surveyor / Poly Peak Parser). It also provides an enterprise multi-agent supervision layer with HMAC-SHA256 audit trails and zero-PHI outbound guards.

### Two Entry Points

1. **Core algorithm** — call heterozygote positions from an `.ab1` file:
   ```bash
   python -m het_caller.cli <file.ab1> [options]
   ```

2. **Enterprise CLI** — multi-agent task evaluation, batch processing, and REST API:
   ```bash
   python cli.py <command> [options]
   ```

---

## ⚙️ Core Algorithm (`het_caller/`)

Implements the peak-ratio heuristic described by Hill et al., *BioTechniques* 2014:

- **Four-channel trace parsing** from ABIF `DATA9-12` tags.
- **Peak detection** with a configurable window around each called base.
- **Heterozygote flagging** when a secondary channel exceeds `threshold` × primary height.
- **IUPAC ambiguity code** assignment (e.g. `R` = A/G, `Y` = C/T).
- **Reference alignment** (BioPython `PairwiseAligner`) for reporting variants in reference coordinates.
- **Quality summary** with estimated signal-to-noise ratio.

### Parameters

| Flag | Description | Default |
|:-----|:------------|:--------|
| `ab1_file` | Path to input `.ab1` chromatogram | (required) |
| `--threshold` | Min secondary/primary ratio to flag | 0.25 |
| `--window` | +/- samples around peak for local max | 2 |
| `--min-quality` | Ignore bases below this Phred score | 0 (off) |
| `--reference` | FASTA file or raw sequence to align against | None |
| `--context` | Flanking bases shown around variants | 5 |
| `--output` | Write JSON report to file | stdout text |

### Example

```bash
python -m het_caller.cli sample.ab1 --threshold 0.25 --window 2 --reference ref.fasta --output report.json
```

---

## ⚙️ Enterprise CLI (`cli.py`)

### Commands

| Command | Description |
|:--------|:------------|
| `audit` | Run single task evaluation with multi-agent consensus |
| `chat` | Query the supervisory chat interface |
| `batch` | Batch process CSV records |
| `verify-audit` | Verify HMAC audit trail integrity |
| `serve` | Launch FastAPI REST server |

### Examples

```bash
# Single task evaluation
python cli.py audit --task-id TASK-001 --primary 28.5 --secondary 14.2 --critical

# Supervisory chat
python cli.py chat "Explain the current protocol status"

# Batch processing
python cli.py batch -i input.csv -o results.csv

# Verify audit integrity
python cli.py verify-audit

# Start REST server
python cli.py serve --host 127.0.0.1 --port 8000
```

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active regex inspection blocking SSNs, MRNs, phone numbers, emails, DOBs, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation.
* **Configurable Audit Secret:** Set `AUDIT_SECRET_KEY` environment variable in production; ephemeral key generated at runtime with a warning if unset.
* **Path Traversal Protection:** Reference file paths are resolved to absolute paths before access.

---

## 🧪 Testing & Verification

Run the full test suite:

```bash
pytest -v
```

This executes 18 tests covering:
- IUPAC code assignment
- Peak calling and heterozygote flagging
- Trace quality summary
- AB1 parsing with fallbacks
- Reference alignment and variant mapping
- PHI guard enforcement
- Multi-agent worker evaluation
- Supervisor consensus and audit integrity
- Enrichment suite execution

---

## 🐳 Container Deployment

```bash
docker build -t sanger-heterozygote-caller .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=<your-secret> sanger-heterozygote-caller
```

Or with Docker Compose:

```bash
AUDIT_SECRET_KEY=<your-secret> docker-compose up
```

---

## 📦 Dependencies

- `biopython>=1.80` — ABIF trace parsing, sequence alignment
- `numpy>=1.24` — numerical array operations
- `pydantic>=2.0` — data models (enterprise CLI)
- `fastapi`, `uvicorn` — REST API server
- `pytest` — test suite

---

## 📁 Project Structure

```
sanger-heterozygote-caller/
├── het_caller/           # Core algorithm package
│   ├── core.py           # Peak-ratio calling, AB1 parsing, alignment
│   └── cli.py            # Core algorithm CLI
├── agents/               # Enterprise multi-agent framework
│   ├── base.py           # PHI guard, HMAC audit trail
│   ├── models.py         # Pydantic data models
│   ├── workers.py        # Domain worker agents
│   ├── supervisor.py     # Consensus orchestrator
│   ├── api.py            # FastAPI endpoints
│   └── ...
├── cli.py                # Enterprise CLI entry point
├── enrichment.py         # Batch screening & enrichment engines
├── simulator.py          # High-throughput simulation
├── tests/                # Pytest test suite
├── test_het_caller.py    # Core algorithm tests
├── web/                  # Operations console (HTML)
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

---

## 📄 License

MIT License — see [LICENSE](LICENSE).
