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
                    
                    hits.append(hsp_hit)\n            \n            return hits\n            \n        except ET.ParseError as e:\n            print(f\"XML parsing error: {e}\")\n            return []\n        except Exception as e:\n            print(f\"Error parsing XML results: {e}\")\n            return []\n    \n    def parse_text_results(self, text_content: str) -> List[BlastHit]:\n        \"\"\"Parse BLAST results from text format (fallback method).\"\"\"\n        hits = []\n        lines = text_content.split('\\n')\n        \n        current_hit = None\n        in_alignment = False\n        \n        for line in lines:\n            line = line.strip()\n            \n            # Extract query information\n            if line.startswith('Query='):\n                self.query_id = line[6:].strip()\n            elif 'Length=' in line and not current_hit:\n                match = re.search(r'Length=(\\d+)', line)\n                if match:\n                    self.query_length = int(match.group(1))\n            \n            # Detect new hit\n            elif line.startswith('>'):\n                if current_hit:\n                    hits.append(current_hit)\n                current_hit = BlastHit()\n                current_hit.query_id = self.query_id\n                current_hit.subject_title = line[1:].strip()\n                # Extract subject ID (first word)\n                parts = current_hit.subject_title.split()\n                if parts:\n                    current_hit.subject_id = parts[0]\n            \n            # Extract subject length\n            elif current_hit and 'Length=' in line:\n                match = re.search(r'Length=(\\d+)', line)\n                if match:\n                    current_hit.subject_length = int(match.group(1))\n            \n            # Extract scores\n            elif current_hit and 'Score =' in line and 'Expect =' in line:\n                # Parse score line: \"Score = 123 bits (456), Expect = 1e-10\"\n                score_match = re.search(r'Score\\s*=\\s*([\\d.]+)\\s*bits\\s*\\(([\\d.]+)\\)', line)\n                if score_match:\n                    current_hit.bit_score = float(score_match.group(1))\n                    current_hit.score = float(score_match.group(2))\n                \n                expect_match = re.search(r'Expect\\s*=\\s*([\\de.-]+)', line)\n                if expect_match:\n                    current_hit.evalue = float(expect_match.group(1))\n            \n            # Extract identities and alignment info\n            elif current_hit and 'Identities =' in line:\n                # Parse: \"Identities = 123/456 (78%), Positives = 234/456 (89%), Gaps = 12/456 (3%)\"\n                identity_match = re.search(r'Identities\\s*=\\s*(\\d+)/(\\d+)\\s*\\(([\\d.]+)%\\)', line)\n                if identity_match:\n                    identities = int(identity_match.group(1))\n                    current_hit.alignment_length = int(identity_match.group(2))\n                    current_hit.percent_identity = float(identity_match.group(3))\n                    current_hit.mismatches = current_hit.alignment_length - identities\n                \n                gap_match = re.search(r'Gaps\\s*=\\s*(\\d+)/(\\d+)', line)\n                if gap_match:\n                    current_hit.gap_opens = int(gap_match.group(1))\n            \n            # Extract coordinates from alignment lines\n            elif current_hit and line.startswith('Query'):\n                # Parse: \"Query  123  ATCGATCG  130\"\n                parts = line.split()\n                if len(parts) >= 4:\n                    try:\n                        if not current_hit.query_start:\n                            current_hit.query_start = int(parts[1])\n                        current_hit.query_end = int(parts[-1])\n                    except ValueError:\n                        pass\n            \n            elif current_hit and line.startswith('Sbjct'):\n                # Parse: \"Sbjct  456  ATCGATCG  463\"\n                parts = line.split()\n                if len(parts) >= 4:\n                    try:\n                        if not current_hit.subject_start:\n                            current_hit.subject_start = int(parts[1])\n                        current_hit.subject_end = int(parts[-1])\n                    except ValueError:\n                        pass\n        \n        # Add the last hit\n        if current_hit:\n            hits.append(current_hit)\n        \n        # Calculate coverage for all hits\n        for hit in hits:\n            if self.query_length > 0 and hit.query_start and hit.query_end:\n                query_span = abs(hit.query_end - hit.query_start) + 1\n                hit.query_coverage = (query_span / self.query_length) * 100\n            \n            if hit.subject_length > 0 and hit.subject_start and hit.subject_end:\n                subject_span = abs(hit.subject_end - hit.subject_start) + 1\n                hit.subject_coverage = (subject_span / hit.subject_length) * 100\n        \n        return hits\n    \n    def format_results_table(self, hits: List[BlastHit], program_type: str = \"blast\") -> str:\n        \"\"\"Format BLAST hits as a comprehensive table with all standard columns.\"\"\"\n        \n        if not hits:\n            return self._format_no_hits_message()\n        \n        # Sort hits by E-value (best first)\n        sorted_hits = sorted(hits, key=lambda x: x.evalue)\n        \n        # Create header\n        result = \"🧬 BLAST Search Results\\n\"\n        result += \"=\" * 80 + \"\\n\\n\"\n        \n        # Add search summary\n        result += f\"📊 Search Summary:\\n\"\n        result += f\"   Program: {self.program or program_type.upper()}\\n\"\n        result += f\"   Database: {self.database or 'Unknown'}\\n\"\n        result += f\"   Query: {self.query_id or 'Unknown'} (Length: {self.query_length})\\n\"\n        result += f\"   Total Hits: {len(sorted_hits)}\\n\\n\"\n        \n        # Create the main results table\n        result += \"📋 Hits Table (sorted by E-value):\\n\"\n        result += \"-\" * 80 + \"\\n\"\n        \n        # Table headers\n        headers = [\n            \"#\", \"Subject ID\", \"% Identity\", \"Align Length\", \n            \"Mismatches\", \"Gap Opens\", \"Query Start\", \"Query End\",\n            \"Subject Start\", \"Subject End\", \"E-value\", \"Bit Score\"\n        ]\n        \n        # Add frame columns for translated searches\n        if program_type.lower() in ['blastx', 'tblastn', 'tblastx']:\n            if program_type.lower() in ['blastx', 'tblastx']:\n                headers.insert(-2, \"Query Frame\")\n            if program_type.lower() in ['tblastn', 'tblastx']:\n                headers.insert(-2, \"Subject Frame\")\n        \n        # Format table with proper spacing\n        col_widths = [3, 15, 10, 12, 10, 9, 11, 9, 13, 11, 12, 10]\n        if \"Query Frame\" in headers:\n            col_widths.insert(-2, 11)\n        if \"Subject Frame\" in headers:\n            col_widths.insert(-2, 13)\n        \n        # Print headers\n        header_line = \"\"\n        for i, header in enumerate(headers):\n            header_line += f\"{header:<{col_widths[i]}} \"\n        result += header_line + \"\\n\"\n        result += \"-\" * len(header_line) + \"\\n\"\n        \n        # Print data rows\n        for i, hit in enumerate(sorted_hits[:50], 1):  # Limit to top 50 hits\n            row_data = [\n                str(i),\n                hit.subject_id[:14] if hit.subject_id else \"N/A\",\n                f\"{hit.percent_identity:.1f}%\",\n                str(hit.alignment_length),\n                str(hit.mismatches),\n                str(hit.gap_opens),\n                str(hit.query_start),\n                str(hit.query_end),\n                str(hit.subject_start),\n                str(hit.subject_end),\n                f\"{hit.evalue:.2e}\" if hit.evalue < 0.01 else f\"{hit.evalue:.3f}\",\n                f\"{hit.bit_score:.1f}\"\n            ]\n            \n            # Add frame information for translated searches\n            if program_type.lower() in ['blastx', 'tblastx'] and \"Query Frame\" in headers:\n                frame_pos = headers.index(\"Query Frame\")\n                row_data.insert(frame_pos, str(hit.query_frame) if hit.query_frame else \"N/A\")\n            \n            if program_type.lower() in ['tblastn', 'tblastx'] and \"Subject Frame\" in headers:\n                frame_pos = headers.index(\"Subject Frame\")\n                row_data.insert(frame_pos, str(hit.subject_frame) if hit.subject_frame else \"N/A\")\n            \n            # Format row\n            row_line = \"\"\n            for j, data in enumerate(row_data):\n                row_line += f\"{data:<{col_widths[j]}} \"\n            result += row_line + \"\\n\"\n        \n        # Add detailed alignments for top hits\n        result += \"\\n\" + \"=\" * 80 + \"\\n\"\n        result += \"🎯 Detailed Alignments (Top 10 Hits):\\n\"\n        result += \"=\" * 80 + \"\\n\\n\"\n        \n        for i, hit in enumerate(sorted_hits[:10], 1):\n            result += f\"Hit #{i}: {hit.subject_id}\\n\"\n            result += f\"Description: {hit.subject_title[:70]}...\\n\" if len(hit.subject_title) > 70 else f\"Description: {hit.subject_title}\\n\"\n            result += f\"Length: {hit.subject_length}\\n\"\n            result += f\"Score: {hit.bit_score:.1f} bits ({hit.score:.0f}), E-value: {hit.evalue:.2e}\\n\"\n            result += f\"Identities: {hit.percent_identity:.1f}% ({int(hit.percent_identity * hit.alignment_length / 100)}/{hit.alignment_length})\\n\"\n            \n            if hit.query_coverage > 0:\n                result += f\"Query Coverage: {hit.query_coverage:.1f}%\\n\"\n            if hit.subject_coverage > 0:\n                result += f\"Subject Coverage: {hit.subject_coverage:.1f}%\\n\"\n            \n            # Add frame information for translated searches\n            if program_type.lower() in ['blastx', 'tblastx'] and hit.query_frame:\n                result += f\"Query Frame: {hit.query_frame}\\n\"\n            if program_type.lower() in ['tblastn', 'tblastx'] and hit.subject_frame:\n                result += f\"Subject Frame: {hit.subject_frame}\\n\"\n            \n            result += f\"\\nQuery: {hit.query_start:>6} {hit.query_sequence} {hit.query_end}\\n\"\n            if hit.alignment_midline:\n                result += f\"       {' ' * 6} {hit.alignment_midline}\\n\"\n            result += f\"Sbjct: {hit.subject_start:>6} {hit.subject_sequence} {hit.subject_end}\\n\\n\"\n        \n        # Add summary statistics\n        if len(sorted_hits) > 10:\n            result += f\"... and {len(sorted_hits) - 10} more hits\\n\\n\"\n        \n        result += \"📈 Search Statistics:\\n\"\n        result += f\"   Best E-value: {sorted_hits[0].evalue:.2e}\\n\"\n        result += f\"   Best Bit Score: {sorted_hits[0].bit_score:.1f}\\n\"\n        result += f\"   Average Identity: {sum(h.percent_identity for h in sorted_hits) / len(sorted_hits):.1f}%\\n\"\n        \n        return result\n    \n    def _format_no_hits_message(self) -> str:\n        \"\"\"Format a message when no hits are found.\"\"\"\n        result = \"🧬 BLAST Search Results\\n\"\n        result += \"=\" * 50 + \"\\n\\n\"\n        result += \"❌ No significant hits found.\\n\\n\"\n        result += \"💡 Suggestions to improve your search:\\n\"\n        result += \"   • Increase the E-value threshold (try 1.0 or 10.0)\\n\"\n        result += \"   • Use a different scoring matrix\\n\"\n        result += \"   • Check your query sequence for errors\\n\"\n        result += \"   • Try a different database\\n\"\n        result += \"   • Reduce word size for more sensitive search\\n\"\n        result += \"   • Remove low complexity filtering\\n\\n\"\n        return result\n\n\ndef enhanced_format_blast_output(raw_output: str, program_type: str = \"blast\") -> str:\n    \"\"\"Enhanced BLAST output formatter with comprehensive table parsing.\"\"\"\n    \n    if not raw_output or not raw_output.strip():\n        return \"❌ No results found or empty output.\"\n    \n    parser = BlastResultsParser()\n    \n    # Try to parse as XML first (if available), then fall back to text\n    if raw_output.strip().startswith('<?xml') or '<BlastOutput>' in raw_output:\n        hits = parser.parse_xml_results(raw_output)\n    else:\n        hits = parser.parse_text_results(raw_output)\n    \n    return parser.format_results_table(hits, program_type)\n\n\nif __name__ == \"__main__\":\n    # Test the parser with sample data\n    sample_text = \"\"\"\nBLASTP 2.12.0+\n\n\nReference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.\nSchäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.\nLipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of\nprotein database search programs\", Nucleic Acids Res. 25:3389-3402.\n\n\nDatabase: Non-redundant protein sequences (nr)\n           540,732,278 sequences; 196,078,718,434 total letters\n\n\n\nQuery= test_protein\n\nLength=100\n                                                                      Score     E\nSequences producing significant alignments:                          (Bits)  Value\n\nsp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens ...   85.5    2e-18\nsp|P02340|P53_MOUSE Cellular tumor antigen p53 OS=Mus musculus O...   82.0    2e-17\n\n\n>sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 \nGN=TP53 PE=1 SV=4\nLength=393\n\n Score = 85.5 bits (210),  Expect = 2e-18\n Identities = 45/50 (90%), Positives = 47/50 (94%), Gaps = 0/50 (0%)\n\nQuery  1   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP  60\n           MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\nSbjct  1   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP  60\n\nQuery  61  DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK  120\n           DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK\nSbjct  61  DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK  120\n\n\"\"\"\n    \n    parser = BlastResultsParser()\n    hits = parser.parse_text_results(sample_text)\n    formatted = parser.format_results_table(hits, \"blastp\")\n    print(formatted)\n