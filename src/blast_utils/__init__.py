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
    create_blastp_tab,
    create_blastn_tab,
    create_blastx_tab,
    create_tblastn_tab,
    create_tblastx_tab,
    run_online_blast_search,
    cancel_blast_search,
    on_blast_finished,
    on_blast_error,
    on_blast_progress,
    validate_sequence,
    clean_fasta_sequence,
    format_blast_output,
    save_blast_results,
    open_ncbi_blast,
    open_blast_file,
    show_blast_help,
    BLAST_CONFIG
)

from .accession_validator import (
    is_accession_number,
    validate_sequence_or_accession,
    clean_fasta_sequence_or_accession
)

from .blast_results_parser import (
    BlastHit,
    BlastResultsParser,
    enhanced_format_blast_output
)

__all__ = [
    'OnlineBlastWorker',
    'create_blastp_tab',
    'create_blastn_tab', 
    'create_blastx_tab',
    'create_tblastn_tab',
    'create_tblastx_tab',
    'run_online_blast_search',
    'cancel_blast_search',
    'on_blast_finished',
    'on_blast_error',
    'on_blast_progress',
    'validate_sequence',
    'clean_fasta_sequence',
    'format_blast_output',
    'save_blast_results',
    'open_ncbi_blast',
    'open_blast_file',
    'show_blast_help',
    'BLAST_CONFIG',
    'enhanced_format_blast_output',
    'is_accession_number',
    'validate_sequence_or_accession',
    'clean_fasta_sequence_or_accession'
]