"""Sanger Heterozygote Caller.

Peak-ratio heterozygous base-calling for Sanger .ab1 chromatograms.
"""

from .core import (
    IUPAC_CODES,
    PeakCall,
    TraceData,
    align_to_reference,
    call_peaks,
    iupac_code,
    map_flagged_to_reference,
    parse_ab1,
    trace_quality_summary,
)

__version__ = "1.0.0"

__all__ = [
    "IUPAC_CODES",
    "PeakCall",
    "TraceData",
    "align_to_reference",
    "call_peaks",
    "iupac_code",
    "map_flagged_to_reference",
    "parse_ab1",
    "trace_quality_summary",
]
