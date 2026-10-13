# Sanger Heterozygote Caller

Screen candidate mixed-base positions in Sanger AB1 chromatograms using analyzed four-channel peak heights. Provides a **local browser interface** and a **Python command-line tool** with optional reference alignment.

The browser interface is [web/index.html](web/index.html); the repository root includes a GitHub Pages entry page. **Live deployment has not yet been independently verified.**

## Features

- Parse ABIF tags DATA9–DATA12 (analyzed trace), FWO_1 (dye order), PLOC2/PLOC1 (peak locations), PBAS2/PBAS1 (base calls), and PCON2/PCON1 (quality).
- Evaluate a configurable secondary-to-primary peak-height threshold near each called base.
- Apply optional instrument-provided Phred quality filtering.
- Flag candidate A/C/G/T mixtures using two-base IUPAC codes.
- Plot four local fluorescence channels and download flagged calls as CSV or JSON in the browser.
- Optionally align Python CLI calls to a single-sequence FASTA reference or a DNA reference string.

**Scientific limitation:** This is a peak-ratio *screen*, not a validated genotyping assay. Secondary signals can arise from background, artifacts, dye effects, or biological mixture. A ratio alone does not establish heterozygosity, mosaic fraction, or clinical significance. Review chromatograms and obtain orthogonal confirmation where appropriate. The Python API's reported SNR is a heuristic ratio to the bottom-quartile baseline and is not a calibrated analytical SNR.

## Browser usage

Open [web/index.html](web/index.html) from a local HTTP server or correctly configured GitHub Pages deployment. Select an AB1 chromatogram, adjust the ratio threshold, peak radius and minimum quality, and press **Analyze trace**. Select **Load example** for explicitly synthetic data. Inspect candidate rows and plots, then export JSON or CSV.

**Privacy:** The browser processes the file locally in JavaScript. It does not upload the file or make API calls. Following the repository hyperlink navigates to GitHub. The Python CLI also reads local files. The optional **separate** FastAPI server processes data submitted to its endpoints and should not be exposed publicly without appropriate security controls.

The browser uses a compact JavaScript ABIF parser rather than Pyodide, avoiding a heavyweight Python/WebAssembly runtime. The Python implementation remains available for CLI analysis and reference mapping. Check consistency against real chromatogram fixtures before using either workflow in a validated laboratory process.

## Python CLI

Python 3.10 or later is recommended.

    python -m pip install -r requirements.txt
    python -m het_caller.cli sample.ab1 --threshold 0.25 --window 2 --min-quality 20 --output results.json
    python -m het_caller.cli sample.ab1 --reference reference.fasta --context 5 --output aligned.json

| Option | Purpose | Default |
| --- | --- | --- |
| --threshold | Minimum second-to-first ratio (greater than 0, up to 1) | 0.25 |
| --window | Peak search radius in trace samples | 2 |
| --min-quality | Exclude calls below instrument quality | 0 (off) |
| --reference | Single FASTA sequence path or DNA string | Not set |
| --context | Reference flanking bases | 5 |
| --output | Write JSON rather than console text | Standard output |

Python flagged_positions are **0-based**; reference coordinates and the browser table/export are **1-based**.

## Optional API and supervisory CLI

    python cli.py audit --task-id TEST-01 --primary 28.5 --secondary 14.2
    python cli.py batch --input tasks.csv --output results.csv
    python cli.py verify-audit
    python cli.py serve --host 127.0.0.1 --port 8000

FastAPI routes include /health, /metrics, /api/audit, /api/chat, and /api/audit/logs. **The API cannot run on GitHub Pages.** This subsystem evaluates generic numerical threshold rules and is separate from chromatogram analysis. Provider selection currently resolves to deterministic mock responses; it does not conduct remote model inference.

The PHI regex guard is only a limited pattern detector, **not** comprehensive de-identification or a compliance guarantee. The audit entries are HMAC-signed and chained but stored in memory. Set AUDIT_SECRET_KEY to use a stable signing key; persistence and key management are needed for production auditability.

## Development and tests

    python -m pip install -r requirements.txt pytest
    python -m pytest -q
    python -m compileall -q het_caller agents cli.py enrichment.py simulator.py
    node --check web/abif.js
    node --check web/app.js
    node --test tests/web.test.cjs
    python -m http.server 8000

Then visit http://localhost:8000/web/. CI covers Python 3.10–3.12 and browser binary parsing with Node 22. Real instrument files and multiple mobile browsers should still be tested.

## Technology and license

Python, NumPy, Biopython, FastAPI, Pydantic, plain JavaScript, HTML and CSS. Modern browsers need File.arrayBuffer, DataView and Canvas. No build step is required to serve the static browser UI.

MIT license; see [LICENSE](LICENSE).
