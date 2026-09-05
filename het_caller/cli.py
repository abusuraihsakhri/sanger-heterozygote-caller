"""Command-line interface for the Sanger Heterozygote Caller."""

from __future__ import annotations

import argparse
import json
import os
import sys

from Bio import SeqIO

from .core import call_peaks, map_flagged_to_reference, parse_ab1, trace_quality_summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sanger-het-caller",
        description=(
            "Flag heterozygous / low-level mosaic positions in a Sanger .ab1 "
            "chromatogram via secondary-to-primary fluorescence peak-height ratio."
        ),
    )
    parser.add_argument("ab1_file", help="Path to the input .ab1 chromatogram file")
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.25,
        help="Minimum secondary/primary peak-height ratio to flag a position (default: 0.25)",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=2,
        help="Search +/- this many trace samples around each called peak for the true "
        "local maximum per channel (default: 2)",
    )
    parser.add_argument(
        "--min-quality",
        type=int,
        default=0,
        help="Ignore called bases below this quality score, 0 disables filtering (default: 0)",
    )
    parser.add_argument(
        "--reference",
        default=None,
        help="Path to a FASTA file, or a raw sequence string, to align against for "
        "reporting flagged-variant reference position and context",
    )
    parser.add_argument(
        "--context",
        type=int,
        default=5,
        help="Flanking bases shown around each flagged variant in reference context (default: 5)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Write the full report as JSON to this path (default: print a text summary to stdout)",
    )
    return parser


def _load_reference(spec: str) -> str:
    """Load a reference sequence from a FASTA/plain-text file path, or treat `spec` as raw sequence."""
    if os.path.isfile(spec):
        # Resolve to absolute, normalized path to prevent traversal surprises.
        abs_path = os.path.realpath(spec)
        with open(abs_path) as handle:
            text = handle.read()
        if text.lstrip().startswith(">"):
            lines = text.splitlines()[1:]
            return "".join(line.strip() for line in lines)
        return "".join(text.split())
    return "".join(spec.split())


def _build_report(trace, calls, summary, args) -> dict:
    report = {
        "sample": trace.sample_name,
        "parameters": {
            "threshold": args.threshold,
            "window": args.window,
            "min_quality": args.min_quality,
        },
        "summary": summary,
        "flagged_positions": [
            {
                "position": c.position,
                "trace_index": c.trace_index,
                "called_base": c.called_base,
                "primary_base": c.primary_base,
                "primary_height": c.primary_height,
                "secondary_base": c.secondary_base,
                "secondary_height": c.secondary_height,
                "ratio": round(c.ratio, 3),
                "iupac": c.iupac,
            }
            for c in calls
            if c.flagged
        ],
    }

    if args.reference:
        reference_sequence = _load_reference(args.reference)
        report["reference_mapping"] = map_flagged_to_reference(
            calls, reference_sequence, trace.called_bases, context=args.context
        )

    return report


def _print_text_report(report: dict) -> None:
    summary = report["summary"]
    print(f"Sample: {report['sample'] or '(unnamed)'}")
    print(
        f"Bases called: {summary['num_bases']}  |  "
        f"Flagged positions: {summary['num_flagged']}  |  "
        f"SNR: {summary['signal_to_noise_ratio']}"
    )
    print(
        f"Mean primary peak height: {summary['mean_primary_peak_height']}  |  "
        f"Estimated noise floor: {summary['estimated_noise_floor']}"
    )
    print()

    if not report["flagged_positions"]:
        print("No positions exceeded the heterozygosity peak-ratio threshold.")
    else:
        print(f"{'Pos':>6} {'Base':>4} {'IUPAC':>5} {'Primary':>9} {'Secondary':>10} {'Ratio':>6}")
        for f in report["flagged_positions"]:
            print(
                f"{f['position']:>6} {f['called_base']:>4} {f['iupac']:>5} "
                f"{f['primary_base']}:{f['primary_height']:>7} "
                f"{f['secondary_base']}:{f['secondary_height']:>7} "
                f"{f['ratio']:>6.3f}"
            )

    if "reference_mapping" in report:
        print()
        if not report["reference_mapping"]:
            print("Reference alignment: no flagged positions mapped to the reference.")
        else:
            print("Reference-mapped variants:")
            for m in report["reference_mapping"]:
                pos = m["reference_position"] if m["reference_position"] is not None else "unaligned"
                print(f"  ref pos {pos}: {m['iupac']} ({m['primary_base']}/{m['secondary_base']}) "
                      f"context={m['reference_context']}")


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    try:
        record = SeqIO.read(args.ab1_file, "abi")
    except Exception as exc:  # noqa: BLE001 - surface parse errors as a clean CLI failure
        print(f"Error reading '{args.ab1_file}': {exc}", file=sys.stderr)
        return 1

    trace = parse_ab1(record)
    calls = call_peaks(trace, threshold=args.threshold, window=args.window)

    if args.min_quality:
        calls = [c for c in calls if trace.qualities[c.position] >= args.min_quality]

    summary = trace_quality_summary(trace, calls)
    report = _build_report(trace, calls, summary, args)

    if args.output:
        with open(args.output, "w") as f:
            json.dump(report, f, indent=2)
        print(f"Report written to {args.output}")
    else:
        _print_text_report(report)

    return 0


if __name__ == "__main__":
    sys.exit(main())
