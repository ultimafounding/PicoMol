#!/usr/bin/env python3
"""
Enhanced BLAST results parser for PicoMol.

This module provides comprehensive parsing of BLAST results to extract all standard
columns and format them as proper tables with all the information users expect.
"""

import re
import xml.etree.ElementTree as ET
from typing import List, Dict, Any, Optional, Tuple


class BlastHit:
    """Represents a single BLAST hit with all relevant information."""
    
    def __init__(self):
        self.query_id = ""
        self.subject_id = ""
        self.subject_title = ""
        self.subject_length = 0
        self.percent_identity = 0.0
        self.alignment_length = 0
        self.mismatches = 0
        self.gap_opens = 0
        self.query_start = 0
        self.query_end = 0
        self.subject_start = 0
        self.subject_end = 0
        self.evalue = 0.0
        self.bit_score = 0.0
        self.score = 0.0
        self.query_coverage = 0.0
        self.subject_coverage = 0.0
        self.query_frame = None
        self.subject_frame = None
        self.query_sequence = ""
        self.subject_sequence = ""
        self.alignment_midline = ""


class BlastResultsParser:
    """Parser for BLAST results in various formats."""
    
    def __init__(self):
        self.hits = []
        self.query_length = 0
        self.query_id = ""
        self.database = ""
        self.program = ""
        self.version = ""
        self.search_stats = {}
    
    def parse_xml_results(self, xml_content: str) -> List[BlastHit]:
        """Parse BLAST results from XML format."""
        try:
            root = ET.fromstring(xml_content)
            hits = []
            
            # Extract query information
            query_def = root.find('.//BlastOutput_query-def')
            if query_def is not None:
                self.query_id = query_def.text
            
            query_len = root.find('.//BlastOutput_query-len')
            if query_len is not None:
                self.query_length = int(query_len.text)
            
            # Extract program and database info
            program = root.find('.//BlastOutput_program')
            if program is not None:
                self.program = program.text
            
            database = root.find('.//BlastOutput_db')
            if database is not None:
                self.database = database.text
            
            # Parse hits
            for hit_elem in root.findall('.//Hit'):
                hit = BlastHit()
                
                # Basic hit information
                hit_id = hit_elem.find('Hit_id')
                if hit_id is not None:
                    hit.subject_id = hit_id.text
                
                hit_def = hit_elem.find('Hit_def')
                if hit_def is not None:
                    hit.subject_title = hit_def.text
                
                hit_len = hit_elem.find('Hit_len')
                if hit_len is not None:
                    hit.subject_length = int(hit_len.text)
                
                # Parse HSPs (High-scoring Segment Pairs)
                for hsp_elem in hit_elem.findall('.//Hsp'):
                    hsp_hit = BlastHit()
                    hsp_hit.subject_id = hit.subject_id
                    hsp_hit.subject_title = hit.subject_title
                    hsp_hit.subject_length = hit.subject_length
                    hsp_hit.query_id = self.query_id
                    
                    # HSP-specific data
                    score = hsp_elem.find('Hsp_score')
                    if score is not None:
                        hsp_hit.score = float(score.text)
                    
                    bit_score = hsp_elem.find('Hsp_bit-score')
                    if bit_score is not None:
                        hsp_hit.bit_score = float(bit_score.text)
                    
                    evalue = hsp_elem.find('Hsp_evalue')
                    if evalue is not None:
                        hsp_hit.evalue = float(evalue.text)
                    
                    identity = hsp_elem.find('Hsp_identity')
                    align_len = hsp_elem.find('Hsp_align-len')
                    if identity is not None and align_len is not None:
                        hsp_hit.alignment_length = int(align_len.text)
                        identities = int(identity.text)
                        hsp_hit.percent_identity = (identities / hsp_hit.alignment_length) * 100
                    
                    gaps = hsp_elem.find('Hsp_gaps')
                    if gaps is not None:
                        hsp_hit.gap_opens = int(gaps.text)
                    
                    # Query coordinates
                    query_from = hsp_elem.find('Hsp_query-from')
                    query_to = hsp_elem.find('Hsp_query-to')
                    if query_from is not None and query_to is not None:
                        hsp_hit.query_start = int(query_from.text)
                        hsp_hit.query_end = int(query_to.text)
                    
                    # Subject coordinates
                    hit_from = hsp_elem.find('Hsp_hit-from')
                    hit_to = hsp_elem.find('Hsp_hit-to')
                    if hit_from is not None and hit_to is not None:
                        hsp_hit.subject_start = int(hit_from.text)
                        hsp_hit.subject_end = int(hit_to.text)
                    
                    # Frames (for translated searches)
                    query_frame = hsp_elem.find('Hsp_query-frame')
                    if query_frame is not None:
                        hsp_hit.query_frame = int(query_frame.text)
                    
                    hit_frame = hsp_elem.find('Hsp_hit-frame')
                    if hit_frame is not None:
                        hsp_hit.subject_frame = int(hit_frame.text)
                    
                    # Alignment sequences
                    qseq = hsp_elem.find('Hsp_qseq')
                    if qseq is not None:
                        hsp_hit.query_sequence = qseq.text
                    
                    hseq = hsp_elem.find('Hsp_hseq')
                    if hseq is not None:
                        hsp_hit.subject_sequence = hseq.text
                    
                    midline = hsp_elem.find('Hsp_midline')
                    if midline is not None:
                        hsp_hit.alignment_midline = midline.text
                    
                    # Calculate coverage
                    if self.query_length > 0:
                        query_span = abs(hsp_hit.query_end - hsp_hit.query_start) + 1
                        hsp_hit.query_coverage = (query_span / self.query_length) * 100
                    
                    if hsp_hit.subject_length > 0:
                        subject_span = abs(hsp_hit.subject_end - hsp_hit.subject_start) + 1
                        hsp_hit.subject_coverage = (subject_span / hsp_hit.subject_length) * 100
                    
                    # Calculate mismatches
                    if hsp_hit.alignment_length > 0:
                        identities = int((hsp_hit.percent_identity / 100) * hsp_hit.alignment_length)
                        hsp_hit.mismatches = hsp_hit.alignment_length - identities - hsp_hit.gap_opens
                    
                    hits.append(hsp_hit)
            
            return hits
            
        except ET.ParseError as e:
            print(f"XML parsing error: {e}")
            return []
        except Exception as e:
            print(f"Error parsing XML results: {e}")
            return []

    def parse_text_results(self, text_content: str) -> List[BlastHit]:
        """Parse BLAST results from text format (fallback method)."""
        hits = []
        lines = text_content.split('\n')
        
        current_hit = None
        in_alignment = False
        
        for line in lines:
            line = line.strip()
            
            # Extract query information
            if line.startswith('Query='):
                self.query_id = line[6:].strip()
            elif 'Length=' in line and not current_hit:
                match = re.search(r'Length=(\d+)', line)
                if match:
                    self.query_length = int(match.group(1))
            
            # Detect new hit
            elif line.startswith('>'):
                if current_hit:
                    hits.append(current_hit)
                current_hit = BlastHit()
                current_hit.query_id = self.query_id
                current_hit.subject_title = line[1:].strip()
                # Extract subject ID (first word)
                parts = current_hit.subject_title.split()
                if parts:
                    current_hit.subject_id = parts[0]
            
            # Extract subject length
            elif current_hit and 'Length=' in line:
                match = re.search(r'Length=(\d+)', line)
                if match:
                    current_hit.subject_length = int(match.group(1))
            
            # Extract scores
            elif current_hit and 'Score =' in line and 'Expect =' in line:
                # Parse score line: "Score = 123 bits (456), Expect = 1e-10"
                score_match = re.search(r'Score\s*=\s*([\d.]+)\s*bits\s*\(([\d.]+)\)', line)
                if score_match:
                    current_hit.bit_score = float(score_match.group(1))
                    current_hit.score = float(score_match.group(2))
                
                expect_match = re.search(r'Expect\s*=\s*([\de.-]+)', line)
                if expect_match:
                    current_hit.evalue = float(expect_match.group(1))
            
            # Extract identities and alignment info
            elif current_hit and 'Identities =' in line:
                # Parse: "Identities = 123/456 (78%), Positives = 234/456 (89%), Gaps = 12/456 (3%)"
                identity_match = re.search(r'Identities\s*=\s*(\d+)/(\d+)\s*\(([\d.]+)%\)', line)
                if identity_match:
                    identities = int(identity_match.group(1))
                    current_hit.alignment_length = int(identity_match.group(2))
                    current_hit.percent_identity = float(identity_match.group(3))
                    current_hit.mismatches = current_hit.alignment_length - identities
                
                gap_match = re.search(r'Gaps\s*=\s*(\d+)/(\d+)', line)
                if gap_match:
                    current_hit.gap_opens = int(gap_match.group(1))
            
            # Extract coordinates from alignment lines
            elif current_hit and line.startswith('Query'):
                # Parse: "Query  123  ATCGATCG  130"
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        if not current_hit.query_start:
                            current_hit.query_start = int(parts[1])
                        current_hit.query_end = int(parts[-1])
                    except ValueError:
                        pass
            
            elif current_hit and line.startswith('Sbjct'):
                # Parse: "Sbjct  456  ATCGATCG  463"
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        if not current_hit.subject_start:
                            current_hit.subject_start = int(parts[1])
                        current_hit.subject_end = int(parts[-1])
                    except ValueError:
                        pass
        
        # Add the last hit
        if current_hit:
            hits.append(current_hit)
        
        # Calculate coverage for all hits
        for hit in hits:
            if self.query_length > 0 and hit.query_start and hit.query_end:
                query_span = abs(hit.query_end - hit.query_start) + 1
                hit.query_coverage = (query_span / self.query_length) * 100
            
            if hit.subject_length > 0 and hit.subject_start and hit.subject_end:
                subject_span = abs(hit.subject_end - hit.subject_start) + 1
                hit.subject_coverage = (subject_span / hit.subject_length) * 100
        
        return hits
    
    def format_results_table(self, hits: List[BlastHit], program_type: str = "blast") -> str:
        """Format BLAST hits as a comprehensive table with all standard columns."""
        
        if not hits:
            return self._format_no_hits_message()
        
        # Sort hits by E-value (best first)
        sorted_hits = sorted(hits, key=lambda x: x.evalue)
        
        # Create header
        result = " BLAST Search Results\n"
        result += "=" * 80 + "\n\n"
        
        # Add search summary
        result += " Search Summary:\n"
        result += f"   Program: {self.program or program_type.upper()}\n"
        result += f"   Database: {self.database or 'Unknown'}\n"
        result += f"   Query: {self.query_id or 'Unknown'} (Length: {self.query_length})\n"
        result += f"   Total Hits: {len(sorted_hits)}\n\n"
        
        # Create the main results table
        result += " Hits Table (sorted by E-value):\n"
        result += "-" * 80 + "\n"
        
        # Table headers
        headers = [
            "#", "Subject ID", "% Identity", "Align Length", 
            "Mismatches", "Gap Opens", "Query Start", "Query End",
            "Subject Start", "Subject End", "E-value", "Bit Score"
        ]
        
        # Add frame columns for translated searches
        if program_type.lower() in ['blastx', 'tblastn', 'tblastx']:
            if program_type.lower() in ['blastx', 'tblastx']:
                headers.insert(-2, "Query Frame")
            if program_type.lower() in ['tblastn', 'tblastx']:
                headers.insert(-2, "Subject Frame")
        
        # Format table with proper spacing
        col_widths = [3, 15, 10, 12, 10, 9, 11, 9, 13, 11, 12, 10]
        if "Query Frame" in headers:
            col_widths.insert(-2, 11)
        if "Subject Frame" in headers:
            col_widths.insert(-2, 13)
        
        # Print headers
        header_line = ""
        for i, header in enumerate(headers):
            header_line += f"{header:<{col_widths[i]}} "
        result += header_line + "\n"
        result += "-" * len(header_line) + "\n"
        
        # Print data rows
        for i, hit in enumerate(sorted_hits[:50], 1):  # Limit to top 50 hits
            row_data = [
                str(i),
                hit.subject_id[:14] if hit.subject_id else "N/A",
                f"{hit.percent_identity:.1f}%",
                str(hit.alignment_length),
                str(hit.mismatches),
                str(hit.gap_opens),
                str(hit.query_start),
                str(hit.query_end),
                str(hit.subject_start),
                str(hit.subject_end),
                f"{hit.evalue:.2e}" if hit.evalue < 0.01 else f"{hit.evalue:.3f}",
                f"{hit.bit_score:.1f}"
            ]
            
            # Add frame information for translated searches
            if program_type.lower() in ['blastx', 'tblastx'] and "Query Frame" in headers:
                frame_pos = headers.index("Query Frame")
                row_data.insert(frame_pos, str(hit.query_frame) if hit.query_frame else "N/A")
            
            if program_type.lower() in ['tblastn', 'tblastx'] and "Subject Frame" in headers:
                frame_pos = headers.index("Subject Frame")
                row_data.insert(frame_pos, str(hit.subject_frame) if hit.subject_frame else "N/A")
            
            # Format row
            row_line = ""
            for j, data in enumerate(row_data):
                row_line += f"{data:<{col_widths[j]}} "
            result += row_line + "\n"
        
        # Add detailed alignments for top hits
        result += "\n" + "=" * 80 + "\n"
        result += " Detailed Alignments (Top 10 Hits):\n"
        result += "=" * 80 + "\n\n"
        
        for i, hit in enumerate(sorted_hits[:10], 1):
            result += f"Hit #{i}: {hit.subject_id}\n"
            result += f"Description: {hit.subject_title[:70]}...\n" if len(hit.subject_title) > 70 else f"Description: {hit.subject_title}\n"
            result += f"Length: {hit.subject_length}\n"
            result += f"Score: {hit.bit_score:.1f} bits ({hit.score:.0f}), E-value: {hit.evalue:.2e}\n"
            result += f"Identities: {hit.percent_identity:.1f}% ({int(hit.percent_identity * hit.alignment_length / 100)}/{hit.alignment_length})\n"
            
            if hit.query_coverage > 0:
                result += f"Query Coverage: {hit.query_coverage:.1f}%\n"
            if hit.subject_coverage > 0:
                result += f"Subject Coverage: {hit.subject_coverage:.1f}%\n"
            
            # Add frame information for translated searches
            if program_type.lower() in ['blastx', 'tblastx'] and hit.query_frame:
                result += f"Query Frame: {hit.query_frame}\n"
            if program_type.lower() in ['tblastn', 'tblastx'] and hit.subject_frame:
                result += f"Subject Frame: {hit.subject_frame}\n"
            
            result += f"\nQuery: {hit.query_start:>6} {hit.query_sequence} {hit.query_end}\n"
            if hit.alignment_midline:
                result += f"       {' ' * 6} {hit.alignment_midline}\n"
            result += f"Sbjct: {hit.subject_start:>6} {hit.subject_sequence} {hit.subject_end}\n\n"
        
        # Add summary statistics
        if len(sorted_hits) > 10:
            result += f"... and {len(sorted_hits) - 10} more hits\n\n"
        
        result += " Search Statistics:\n"
        result += f"   Best E-value: {sorted_hits[0].evalue:.2e}\n"
        result += f"   Best Bit Score: {sorted_hits[0].bit_score:.1f}\n"
        result += f"   Average Identity: {sum(h.percent_identity for h in sorted_hits) / len(sorted_hits):.1f}%\n"
        
        return result
    
    def _format_no_hits_message(self) -> str:
        """Format a message when no hits are found."""
        result = " BLAST Search Results\n"
        result += "=" * 50 + "\n\n"
        result += " No significant hits found.\n\n"
        result += " Suggestions to improve your search:\n"
        result += "   • Increase the E-value threshold (try 1.0 or 10.0)\n"
        result += "   • Use a different scoring matrix\n"
        result += "   • Check your query sequence for errors\n"
        result += "   • Try a different database\n"
        result += "   • Reduce word size for more sensitive search\n"
        result += "   • Remove low complexity filtering\n\n"
        return result

def enhanced_format_blast_output(raw_output: str, program_type: str = "blast") -> str:
    """Enhanced BLAST output formatter with comprehensive table parsing."""
    
    if not raw_output or not raw_output.strip():
        return "❌ No results found or empty output."
    
    parser = BlastResultsParser()
    
    # Try to parse as XML first (if available), then fall back to text
    if raw_output.strip().startswith('<?xml') or '<BlastOutput>' in raw_output:
        hits = parser.parse_xml_results(raw_output)
    else:
        hits = parser.parse_text_results(raw_output)
    
    return parser.format_results_table(hits, program_type)


if __name__ == "__main__":
    # Test the parser with sample data
    sample_text = """
BLASTP 2.12.0+

Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.
Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.
Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of
protein database search programs", Nucleic Acids Res. 25:3389-3402.

Database: Non-redundant protein sequences (nr)
           540,732,278 sequences; 196,078,718,434 total letters

Query= test_protein

Length=100
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens ...   85.5    2e-18
sp|P02340|P53_MOUSE Cellular tumor antigen p53 OS=Mus musculus O...   82.0    2e-17

>sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 
GN=TP53 PE=1 SV=4
Length=393

Score = 85.5 bits (210),  Expect = 2e-18
Identities = 45/50 (90%), Positives = 47/50 (94%), Gaps = 0/50 (0%)

Query  1   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP  60
           MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
Sbjct  1   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP  60

Query  61  DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK  120
           DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
Sbjct  61  DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK  120
"""
    
    parser = BlastResultsParser()
    hits = parser.parse_text_results(sample_text)
    formatted = parser.format_results_table(hits, "blastp")
    print(formatted)