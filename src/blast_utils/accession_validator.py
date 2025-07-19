#!/usr/bin/env python3
"""
Accession number validation module for PicoMol BLAST utilities.

This module provides functions to detect and validate various types of 
biological database accession numbers and identifiers.
"""

import re


def is_accession_number(sequence):
    """Check if the input is likely an accession number or GI number."""
    
    # Remove any whitespace
    seq = sequence.strip()
    
    # Common accession number patterns
    accession_patterns = [
        # GenBank accession numbers
        r'^[A-Z]{1,2}[0-9]{5,8}(\.[0-9]+)?$',  # e.g., NM_000001.1, AF123456
        r'^[A-Z]{4}[0-9]{8,10}(\.[0-9]+)?$',   # e.g., AAAA01000001.1
        
        # RefSeq accession numbers
        r'^(NM|NR|XM|XR|NP|XP|YP|AP|NC|NG|NT|NW|NZ)_[0-9]{6,9}(\.[0-9]+)?$',
        
        # UniProt accession numbers
        r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$',       # e.g., O43175
        r'^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$',  # e.g., P12345
        
        # EMBL/ENA accession numbers
        r'^[A-Z]{1,2}[0-9]{5,8}$',             # e.g., X12345, AB123456
        
        # PDB IDs
        r'^[0-9][A-Z0-9]{3}$',                 # e.g., 1ABC, 2XYZ
        
        # GI numbers
        r'^[0-9]{6,12}$',                      # e.g., 123456789
        
        # Ensembl IDs
        r'^ENS[A-Z]*[GPRT][0-9]{11}(\.[0-9]+)?$',  # e.g., ENSG00000139618
        
        # WormBase IDs
        r'^WBGene[0-9]{8}$',                   # e.g., WBGene00000001
        
        # FlyBase IDs
        r'^FBgn[0-9]{7}$',                     # e.g., FBgn0000001
        
        # MGI IDs
        r'^MGI:[0-9]{7}$',                     # e.g., MGI:1234567
        
        # NCBI Gene IDs
        r'^[0-9]{4,8}$',                       # e.g., 12345
        
        # NCBI Protein IDs
        r'^[A-Z]{3}[0-9]{5}(\.[0-9]+)?$',     # e.g., AAA12345.1
        
        # Swiss-Prot/TrEMBL IDs
        r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$',      # e.g., P12345
        r'^[A-Z][0-9][A-Z0-9]{3}[0-9]_[A-Z0-9]+$',  # e.g., P12345_HUMAN
        
        # PIR IDs
        r'^[A-Z][0-9]{5}$',                   # e.g., A12345
        
        # PRF IDs
        r'^[0-9]{7}[A-Z]$',                   # e.g., 1234567A
    ]
    
    # Check against patterns
    for pattern in accession_patterns:
        if re.match(pattern, seq, re.IGNORECASE):
            return True
    
    # Additional heuristics
    # If it's short (< 50 chars), contains no typical sequence characters, 
    # and has numbers/letters in typical accession format
    if len(seq) < 50:
        # Check if it looks like an identifier rather than a sequence
        has_numbers = any(c.isdigit() for c in seq)
        has_letters = any(c.isalpha() for c in seq)
        has_special_chars = any(c in '._-|:' for c in seq)
        
        # If it has both letters and numbers, and is relatively short, likely an accession
        if has_numbers and has_letters and len(seq) <= 30:
            # Make sure it doesn't look like a short sequence
            sequence_chars = set('ACDEFGHIKLMNPQRSTVWYATCGRYSWKMBDHVNXBZJUO*-')
            non_sequence_chars = set(seq.upper()) - sequence_chars
            
            # If it has characters that aren't typical sequence characters, likely an accession
            if non_sequence_chars or (has_special_chars and len(seq) <= 20):
                return True
    
    return False


def validate_sequence_or_accession(sequence, sequence_type="protein"):
    """Validate a biological sequence or accession number."""
    
    if not sequence or not sequence.strip():
        return False, "Sequence is empty"
    
    # Remove FASTA header if present
    lines = sequence.strip().split('\n')
    seq_lines = [line.strip() for line in lines if not line.startswith('>')]
    clean_sequence = ''.join(seq_lines).upper().replace(' ', '').replace('\t', '')
    
    if not clean_sequence:
        return False, "No sequence found after removing headers"
    
    # Check if this might be an accession number or GI number
    if is_accession_number(clean_sequence):
        return True, f"Valid accession number: {clean_sequence}"
    
    # Check minimum length for sequences
    if len(clean_sequence) < 3:
        return False, f"Sequence too short ({len(clean_sequence)} characters). Minimum length is 3."
    
    # Check maximum length for online BLAST
    if len(clean_sequence) > 20000:
        return False, f"Sequence too long ({len(clean_sequence)} characters). Maximum length for online BLAST is 20,000."
    
    if sequence_type == "protein":
        # Standard amino acid codes (including ambiguous codes)
        valid_chars = set("ACDEFGHIKLMNPQRSTVWYXBZJUO*-")
        invalid_chars = set(clean_sequence) - valid_chars
        if invalid_chars:
            return False, f"Invalid amino acid characters found: {', '.join(sorted(invalid_chars))}"
        
        # Check for reasonable protein composition
        if len(clean_sequence) > 10:
            nucleotide_chars = set("ATCG")
            nucleotide_count = sum(1 for char in clean_sequence if char in nucleotide_chars)
            if nucleotide_count / len(clean_sequence) > 0.8:
                return False, "Sequence appears to be nucleotide, not protein"
    elif sequence_type == "nucleotide":
        # Standard nucleotide codes (including ambiguous codes)
        valid_chars = set("ATCGRYSWKMBDHVN-")
        invalid_chars = set(clean_sequence) - valid_chars
        if invalid_chars:
            return False, f"Invalid nucleotide characters found: {', '.join(sorted(invalid_chars))}"
        
        # Check for reasonable nucleotide composition
        if len(clean_sequence) > 10:
            protein_chars = set("DEFHIKLMNPQRSVWY")
            protein_count = sum(1 for char in clean_sequence if char in protein_chars)
            if protein_count / len(clean_sequence) > 0.3:
                return False, "Sequence appears to be protein, not nucleotide"
    
    return True, f"Valid {sequence_type} sequence ({len(clean_sequence)} characters)"


def clean_fasta_sequence_or_accession(sequence):
    """Clean FASTA sequence by removing headers and formatting."""
    
    lines = sequence.strip().split('\n')
    seq_lines = [line.strip() for line in lines if not line.startswith('>')]
    clean_seq = ''.join(seq_lines)
    
    # If it looks like an accession number, return it as-is
    if is_accession_number(clean_seq):
        return clean_seq
    
    return clean_seq