"""Core algorithm for Sanger heterozygote calling.

Implements the peak-ratio heuristic used by chromatogram QC tools such as
Mutation Surveyor and Poly Peak Parser (Hill et al., BioTechniques 2014):
at each base-called position in a Sanger trace, the four fluorescence
channels (A, C, G, T) are read at the called peak location. When a second
channel's intensity exceeds a configurable fraction of the tallest
channel's intensity at that position, the site is flagged as a candidate
heterozygous call (or, at a lower threshold, low-level mosaicism) and
reported as the corresponding IUPAC ambiguity code.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np

# IUPAC nucleotide ambiguity codes keyed by the set of bases they represent.
IUPAC_CODES = {
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("GC"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
    frozenset("CGT"): "B",
    frozenset("AGT"): "D",
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}


def iupac_code(bases: Sequence[str]) -> str:
    """Return the IUPAC ambiguity code for a set of one or more unique bases."""
    unique = frozenset(b.upper() for b in bases)
    if len(unique) == 1:
        return next(iter(unique))
    code = IUPAC_CODES.get(unique)
    if code is None:
        raise ValueError(f"No IUPAC code for base set: {sorted(unique)}")
    return code


@dataclass
class TraceData:
    """Parsed four-channel trace plus per-base peak metadata from an .ab1 file."""

    channels: Dict[str, np.ndarray]  # base -> raw intensity array across the whole trace
    peak_locations: List[int]        # trace-sample index of the called peak, one per base
    called_bases: str                # primary base call string, one character per peak
    qualities: List[int]             # per-base quality score, one per peak
    sample_name: str = ""


@dataclass
class PeakCall:
    position: int          # 0-based index into called_bases / peak_locations
    trace_index: int       # sample index in the raw trace this call was made at
    called_base: str       # base originally called by the instrument software
    primary_base: str      # tallest channel at this position (peak-ratio primary)
    primary_height: int
    secondary_base: str    # second-tallest channel at this position
    secondary_height: int
    ratio: float           # secondary_height / primary_height
    iupac: str             # IUPAC code: primary_base alone, or primary+secondary if flagged
    flagged: bool          # True if ratio >= threshold


def parse_ab1(record) -> TraceData:
    """Extract four-channel trace data from a Biopython SeqIO 'abi' record.

    `record` is the object returned by ``Bio.SeqIO.read(path, "abi")`` (or any
    object exposing the same ``annotations["abif_raw"]`` mapping of raw ABIF
    tags, which is convenient for constructing test fixtures).
    """
    raw = record.annotations["abif_raw"]

    base_order = raw["FWO_1"]
    if isinstance(base_order, bytes):
        base_order = base_order.decode("ascii")

    # DATA9-12 are the "analyzed" (baseline-subtracted) traces, in the
    # channel order given by FWO_1. DATA1-4 hold the unprocessed raw traces.
    data_tags = ["DATA9", "DATA10", "DATA11", "DATA12"]
    channels = {}
    for base, tag in zip(base_order, data_tags):
        channels[base] = np.asarray(raw[tag], dtype=np.float64)

    # PLOC2 holds the (base-caller-finalized) peak locations; fall back to PLOC1.
    ploc_tag = "PLOC2" if "PLOC2" in raw else "PLOC1"
    peak_locations = [int(p) for p in raw[ploc_tag]]

    pbas_tag = "PBAS2" if "PBAS2" in raw else "PBAS1"
    called_bases = raw[pbas_tag]
    if isinstance(called_bases, bytes):
        called_bases = called_bases.decode("ascii")

    pcon_tag = "PCON2" if "PCON2" in raw else "PCON1"
    qual_raw = raw.get(pcon_tag)
    if qual_raw is None:
        qualities = [0] * len(peak_locations)
    elif isinstance(qual_raw, bytes):
        qualities = [int(b) for b in qual_raw]
    else:
        qualities = [int(q) for q in qual_raw]

    sample_name = getattr(record, "id", "") or ""

    return TraceData(
        channels=channels,
        peak_locations=peak_locations,
        called_bases=called_bases,
        qualities=qualities,
        sample_name=sample_name,
    )


def call_peaks(trace: TraceData, threshold: float = 0.25, window: int = 2) -> List[PeakCall]:
    """Apply the secondary/primary peak-height ratio heuristic at every base position.

    For each called base, the four channel traces are searched in a small
    window around the recorded peak location (to tolerate the true local
    maximum sitting a sample or two off the base caller's pick). The two
    tallest channels in that window become the primary and secondary peak;
    a position is flagged when secondary_height / primary_height >= threshold.
    """
    if not 0.0 < threshold <= 1.0:
        raise ValueError("threshold must be in (0, 1]")
    if window < 0:
        raise ValueError("window must be >= 0")

    bases = sorted(trace.channels.keys())
    n_points = len(trace.channels[bases[0]]) if bases else 0

    calls: List[PeakCall] = []
    for i, (loc, called_base) in enumerate(zip(trace.peak_locations, trace.called_bases)):
        lo = max(0, loc - window)
        hi = min(n_points, loc + window + 1)

        heights = {}
        for base in bases:
            segment = trace.channels[base][lo:hi]
            heights[base] = float(segment.max()) if segment.size else 0.0

        ranked = sorted(heights.items(), key=lambda kv: kv[1], reverse=True)
        primary_base, primary_height = ranked[0]
        secondary_base, secondary_height = ranked[1]

        ratio = (secondary_height / primary_height) if primary_height > 0 else 0.0
        flagged = ratio >= threshold and secondary_height > 0

        code = iupac_code([primary_base, secondary_base]) if flagged else primary_base

        calls.append(
            PeakCall(
                position=i,
                trace_index=loc,
                called_base=called_base.upper(),
                primary_base=primary_base,
                primary_height=int(round(primary_height)),
                secondary_base=secondary_base,
                secondary_height=int(round(secondary_height)),
                ratio=ratio,
                iupac=code,
                flagged=flagged,
            )
        )

    return calls


def trace_quality_summary(trace: TraceData, calls: List[PeakCall]) -> dict:
    """Summarize overall trace quality: estimated signal-to-noise ratio and flag count."""
    if trace.channels:
        all_signal = np.concatenate([arr for arr in trace.channels.values()])
        # The bottom quartile of all channel intensities approximates the
        # baseline noise floor, since most of a trace sits at baseline
        # between peaks across all four channels.
        noise_floor = max(float(np.percentile(all_signal, 25)), 1e-6)
    else:
        noise_floor = 1e-6

    primary_heights = [c.primary_height for c in calls]
    mean_signal = float(statistics.mean(primary_heights)) if primary_heights else 0.0
    snr = mean_signal / noise_floor

    return {
        "num_bases": len(calls),
        "num_flagged": sum(1 for c in calls if c.flagged),
        "mean_primary_peak_height": round(mean_signal, 1),
        "estimated_noise_floor": round(noise_floor, 1),
        "signal_to_noise_ratio": round(snr, 2),
    }


def align_to_reference(query_sequence: str, reference_sequence: str):
    """Global-align a called base sequence against a reference. Returns a Bio.Align.Alignment."""
    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(reference_sequence.upper(), query_sequence.upper())
    return alignments[0]


def map_flagged_to_reference(
    calls: List[PeakCall],
    reference_sequence: str,
    query_sequence: str,
    context: int = 5,
) -> List[dict]:
    """Align the called sequence to a reference and report flagged variants in reference coordinates."""
    alignment = align_to_reference(query_sequence, reference_sequence)
    indices = alignment.indices  # shape (2, ncols): row 0 = reference index or -1, row 1 = query index or -1

    query_to_ref: Dict[int, int] = {}
    for ref_idx, q_idx in zip(indices[0], indices[1]):
        if ref_idx != -1 and q_idx != -1:
            query_to_ref[int(q_idx)] = int(ref_idx)

    results = []
    for call in calls:
        if not call.flagged:
            continue

        ref_pos = query_to_ref.get(call.position)
        entry = {
            "query_position": call.position,
            "reference_position": (ref_pos + 1) if ref_pos is not None else None,  # 1-based
            "iupac": call.iupac,
            "primary_base": call.primary_base,
            "secondary_base": call.secondary_base,
            "ratio": round(call.ratio, 3),
            "reference_context": None,
        }
        if ref_pos is not None:
            lo = max(0, ref_pos - context)
            hi = min(len(reference_sequence), ref_pos + context + 1)
            entry["reference_context"] = (
                reference_sequence[lo:ref_pos]
                + "[" + reference_sequence[ref_pos] + "]"
                + reference_sequence[ref_pos + 1:hi]
            )
        results.append(entry)

    return results
