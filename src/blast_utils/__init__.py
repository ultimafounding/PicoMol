"""
BLAST utilities package for PicoMol.

This package contains all BLAST-related functionality including:
- Online BLAST search capabilities
- NCBI-style user interface layouts
- BLAST results parsing and formatting
"""

# Import main functions for easy access
from .blast_utils_ncbi_layout import (
    create_ncbi_style_blastp_tab,
    create_ncbi_style_blastn_tab,
    create_ncbi_style_blastx_tab,
    create_ncbi_style_tblastn_tab,
    create_ncbi_style_tblastx_tab
)

from .blast_utils import (
    OnlineBlastWorker,
    validate_sequence,
    clean_fasta_sequence,
    format_blast_output
)

from .blast_results_parser import (
    BlastHit,
    BlastResultsParser,
    enhanced_format_blast_output
)

__all__ = [
    'create_ncbi_style_blastp_tab',
    'create_ncbi_style_blastn_tab',
    'create_ncbi_style_blastx_tab',
    'create_ncbi_style_tblastn_tab',
    'create_ncbi_style_tblastx_tab',
    'OnlineBlastWorker',
    'validate_sequence',
    'clean_fasta_sequence',
    'format_blast_output',
    'BlastHit',
    'BlastResultsParser',
    'enhanced_format_blast_output'
]