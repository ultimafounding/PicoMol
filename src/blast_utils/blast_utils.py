"""
BLAST utilities for PicoMol - Comprehensive online BLAST functionality using NCBI's API.

This module contains functions and classes to create BLAST interface tabs and handle
online BLAST operations using NCBI's web services with complete parameter support
extracted from the official NCBI BLAST interface, eliminating the need for local
BLAST+ installation.

Supported BLAST programs:
- BLASTN: nucleotide vs nucleotide
- BLASTP: protein vs protein  
- BLASTX: translated nucleotide vs protein
- TBLASTN: protein vs translated nucleotide
- TBLASTX: translated nucleotide vs translated nucleotide

Features:
- Complete parameter sets for all BLAST programs
- All NCBI databases and scoring matrices
- Genetic code support for translated searches
- Advanced filtering and compositional adjustments
- Real-time progress tracking
- Results saving and NCBI integration
"""

import os
import sys
import json
import time
import webbrowser
import xml.etree.ElementTree as ET
from urllib.parse import quote, urlencode
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QTextEdit, QPushButton,
    QComboBox, QLineEdit, QCheckBox, QRadioButton, QButtonGroup,
    QGroupBox, QFormLayout, QFileDialog, QMessageBox, QApplication,
    QScrollArea, QSplitter, QTabWidget, QProgressBar, QDialogButtonBox, QDialog
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont


class OnlineBlastWorker(QThread):
    """Worker thread for running online BLAST searches via NCBI API."""
    
    finished = pyqtSignal(str)  # Emits the result
    error = pyqtSignal(str)     # Emits error message
    progress = pyqtSignal(str)  # Emits progress updates
    
    def __init__(self, sequence, database, program, parameters):
        super().__init__()
        self.sequence = sequence
        self.database = database
        self.program = program
        self.parameters = parameters
        self.request_id = None
        self.cancelled = False
        
    def run(self):
        try:
            self.progress.emit("Submitting BLAST search to NCBI...")
            
            # Submit BLAST search
            self.request_id = self.submit_blast_search()
            if not self.request_id:
                self.error.emit("Failed to submit BLAST search to NCBI.")
                return
            
            self.progress.emit(f"Search submitted (ID: {self.request_id}). Waiting for results...")
            
            # Poll for results
            max_wait_time = 300  # 5 minutes
            wait_time = 0
            poll_interval = 10  # seconds
            
            while wait_time < max_wait_time and not self.cancelled:
                if self.check_search_status():
                    # Search is complete, get results
                    results = self.get_blast_results()
                    if results:
                        # Handle tuple format (format_type, content) or string format
                        if isinstance(results, tuple) and len(results) == 2:
                            format_type, content = results
                            self.finished.emit(content)  # Emit only the content string
                        else:
                            self.finished.emit(results)  # Backward compatibility
                    else:
                        self.error.emit("Failed to retrieve BLAST results.")
                    return
                
                # Wait before next poll
                self.msleep(poll_interval * 1000)
                wait_time += poll_interval
                
                if not self.cancelled:
                    self.progress.emit(f"Still waiting for results... ({wait_time}s elapsed)")
            
            if wait_time >= max_wait_time:
                self.error.emit("BLAST search timed out. The search may still be running on NCBI servers.")
            
        except Exception as e:
            self.error.emit(f"An error occurred during BLAST search: {str(e)}")
    
    def cancel_search(self):
        """Cancel the search."""
        self.cancelled = True
    
    def submit_blast_search(self):
        """Submit BLAST search to NCBI and return request ID."""
        try:
            # Prepare the basic search parameters
            params = {
                'CMD': 'Put',
                'PROGRAM': self.program,
                'DATABASE': self.database,
                'QUERY': self.sequence,
                'FORMAT_TYPE': 'XML',
                'EXPECT': self.parameters.get('evalue', '0.05'),
                'HITLIST_SIZE': self.parameters.get('max_target_seqs', '100'),
                'WORD_SIZE': self.parameters.get('word_size', '3'),
                'FILTER': 'L' if self.parameters.get('low_complexity', True) else 'F',
                'FORMAT_OBJECT': 'SearchInfo'
            }
            
            # Add program-specific parameters
            if self.program in ['blastp', 'blastx', 'tblastn']:
                # Protein-based searches (with compositional adjustments)
                params['MATRIX_NAME'] = self.parameters.get('matrix', 'BLOSUM62')
                params['COMPOSITION_BASED_STATISTICS'] = self.parameters.get('comp_adjust', '2')
                if self.parameters.get('word_threshold'):
                    params['WORD_SCORE_THRESHOLD'] = self.parameters.get('word_threshold', '11')
                if self.parameters.get('ungapped_alignment'):
                    params['UNGAPPED_ALIGNMENT'] = 'T' if self.parameters.get('ungapped_alignment', False) else 'F'
            
            elif self.program == 'tblastx':
                # TBLASTX: matrix only, no compositional adjustments
                params['MATRIX_NAME'] = self.parameters.get('matrix', 'BLOSUM62')
                # TBLASTX also uses nucleotide-style match/mismatch scores
                if self.parameters.get('match_scores'):
                    match_scores = self.parameters.get('match_scores', '2,-3')
                    if ',' in match_scores:
                        match, mismatch = match_scores.split(',', 1)
                        params['MATCH_SCORES'] = f"{match},{mismatch}"
            
            elif self.program == 'blastn':
                # Nucleotide-based searches
                if self.parameters.get('match_scores'):
                    match_scores = self.parameters.get('match_scores', '2,-3')
                    if ',' in match_scores:
                        match, mismatch = match_scores.split(',', 1)
                        params['MATCH_SCORES'] = f"{match},{mismatch}"
                
                # Algorithm-specific parameters for BLASTN
                if self.program == 'blastn':
                    algorithm = self.parameters.get('algorithm', 'megaBlast')
                    if algorithm == 'megaBlast':
                        params['MEGABLAST'] = 'on'
                    elif algorithm == 'discoMegablast':
                        params['TEMPLATE_TYPE'] = 'coding'
                        params['TEMPLATE_LENGTH'] = '18'
            
            # Gap costs (for all programs that support them)
            if self.parameters.get('gapopen') and self.parameters.get('gapextend'):
                gapopen = self.parameters.get('gapopen', '11')
                gapextend = self.parameters.get('gapextend', '1')
                
                # For nucleotide searches, validate gap costs against match/mismatch scores
                if self.program in ['blastn', 'tblastx']:
                    match_scores = self.parameters.get('match_scores', '2,-3')
                    # Use NCBI-compatible gap cost combinations
                    if match_scores == '2,-3':
                        # Valid combinations for 2,-3: (5,2), (4,4), (2,4), (0,4)
                        if f"{gapopen} {gapextend}" not in ["5 2", "4 4", "2 4", "0 4"]:
                            gapopen, gapextend = "5", "2"  # Use default compatible values
                    elif match_scores == '1,-2':
                        # Valid combinations for 1,-2: (2,1), (1,1), (0,1)
                        if f"{gapopen} {gapextend}" not in ["2 1", "1 1", "0 1"]:
                            gapopen, gapextend = "2", "1"
                    elif match_scores == '1,-3':
                        # Valid combinations for 1,-3: (2,1), (1,1), (0,1)
                        if f"{gapopen} {gapextend}" not in ["2 1", "1 1", "0 1"]:
                            gapopen, gapextend = "2", "1"
                    elif match_scores == '1,-4':
                        # Valid combinations for 1,-4: (1,1), (0,1)
                        if f"{gapopen} {gapextend}" not in ["1 1", "0 1"]:
                            gapopen, gapextend = "1", "1"
                    elif match_scores == '4,-5':
                        # Valid combinations for 4,-5: (12,8), (8,6), (6,5), (4,4)
                        if f"{gapopen} {gapextend}" not in ["12 8", "8 6", "6 5", "4 4"]:
                            gapopen, gapextend = "12", "8"
                    elif match_scores == '1,-1':
                        # Valid combinations for 1,-1: (3,1), (2,1), (1,1)
                        if f"{gapopen} {gapextend}" not in ["3 1", "2 1", "1 1"]:
                            gapopen, gapextend = "3", "1"
                
                params['GAPCOSTS'] = f"{gapopen} {gapextend}"
            
            # Genetic code (for translated searches)
            if self.program in ['blastx', 'tblastn', 'tblastx']:
                params['QUERY_GENETIC_CODE'] = self.parameters.get('genetic_code', '1')
                if self.program in ['tblastn', 'tblastx']:
                    params['DB_GENETIC_CODE'] = self.parameters.get('genetic_code', '1')
            
            # Algorithm-specific parameters
            if self.parameters.get('algorithm'):
                algorithm = self.parameters.get('algorithm')
                if algorithm == 'psiBlast':
                    params['PSI_BLAST'] = 'on'
                elif algorithm == 'phiBlast':
                    params['PHI_BLAST'] = 'on'
                    if self.parameters.get('phi_pattern'):
                        params['PHI_PATTERN'] = self.parameters.get('phi_pattern')
                elif algorithm == 'deltaBlast':
                    params['DELTA_BLAST'] = 'on'
            
            # Add organism if specified
            if self.parameters.get('organism'):
                params['ENTREZ_QUERY'] = f'"{self.parameters["organism"]}"[Organism]'
            
            # Encode parameters
            data = urlencode(params).encode('utf-8')
            
            # Submit request
            url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
            request = Request(url, data=data)
            request.add_header('User-Agent', 'PicoMol/1.0 (https://github.com/ultimafounding/PicoMol)')
            
            with urlopen(request, timeout=30) as response:
                result = response.read().decode('utf-8')
            
            # Extract request ID from response - try multiple patterns
            rid = None
            
            # Pattern 1: RID = value
            for line in result.split('\n'):
                if 'RID =' in line:
                    rid = line.split('=')[1].strip()
                    break
            
            # Pattern 2: RID=value (no spaces)
            if not rid:
                for line in result.split('\n'):
                    if 'RID=' in line and 'RID =' not in line:
                        rid = line.split('=')[1].strip()
                        break
            
            # Pattern 3: Look for RID in HTML attributes
            if not rid:
                import re
                rid_match = re.search(r'RID[\s]*=[\s]*["\']?([A-Z0-9]+)["\']?', result, re.IGNORECASE)
                if rid_match:
                    rid = rid_match.group(1)
            
            return rid
            
        except (URLError, HTTPError, Exception) as e:
            raise Exception(f"Failed to submit BLAST search: {str(e)}")
    
    def check_search_status(self):
        """Check if BLAST search is complete."""
        try:
            params = {
                'CMD': 'Get',
                'FORMAT_OBJECT': 'SearchInfo',
                'RID': self.request_id
            }
            
            url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?' + urlencode(params)
            request = Request(url)
            request.add_header('User-Agent', 'PicoMol/1.0 (https://github.com/ultimafounding/PicoMol)')
            
            with urlopen(request, timeout=30) as response:
                result = response.read().decode('utf-8')
            
            # Check status
            if 'Status=WAITING' in result:
                return False
            elif 'Status=READY' in result:
                return True
            elif 'Status=UNKNOWN' in result:
                raise Exception("BLAST search failed or expired")
            
            return False
            
        except (URLError, HTTPError, Exception) as e:
            raise Exception(f"Failed to check search status: {str(e)}")
    
    def get_blast_results(self):
        """Retrieve BLAST search results."""
        try:
            # First try to get XML results for better parsing
            xml_params = {
                'CMD': 'Get',
                'FORMAT_TYPE': 'XML',
                'RID': self.request_id
            }
            
            xml_url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?' + urlencode(xml_params)
            xml_request = Request(xml_url)
            xml_request.add_header('User-Agent', 'PicoMol/1.0 (https://github.com/ultimafounding/PicoMol)')
            
            try:
                with urlopen(xml_request, timeout=60) as response:
                    xml_result = response.read().decode('utf-8')
                
                # Check if XML is valid and contains results
                if xml_result and '<BlastOutput>' in xml_result and 'No hits found' not in xml_result:
                    return ('XML', xml_result)
            except:
                pass  # Fall back to text format
            
            # Fall back to text format
            params = {
                'CMD': 'Get',
                'FORMAT_TYPE': 'Text',
                'RID': self.request_id,
                'SHOW_OVERVIEW': 'on',
                'DESCRIPTIONS': self.parameters.get('max_target_seqs', '100'),
                'ALIGNMENTS': self.parameters.get('max_target_seqs', '100')
            }
            
            url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?' + urlencode(params)
            request = Request(url)
            request.add_header('User-Agent', 'PicoMol/1.0 (https://github.com/ultimafounding/PicoMol)')
            
            with urlopen(request, timeout=60) as response:
                result = response.read().decode('utf-8')
            
            return ('Text', result)
            
        except (URLError, HTTPError, Exception) as e:
            raise Exception(f"Failed to retrieve results: {str(e)}")


# Load comprehensive BLAST configuration
BLAST_CONFIG = {
    "genetic_codes": {
        "1": "Standard (1)",
        "2": "Vertebrate Mitochondrial (2)",
        "3": "Yeast Mitochondrial (3)",
        "4": "Mold Mitochondrial; ... (4)",
        "5": "Invertebrate Mitochondrial (5)",
        "6": "Ciliate Nuclear; ... (6)",
        "9": "Echinoderm Mitochondrial (9)",
        "10": "Euplotid Nuclear (10)",
        "11": "Bacteria and Archaea (11)",
        "12": "Alternative Yeast Nuclear (12)",
        "13": "Ascidian Mitochondrial (13)",
        "14": "Flatworm Mitochondrial (14)",
        "16": "Chlorophycean Mitochondrial (16)",
        "21": "Trematode Mitochondrial (21)",
        "22": "Scenedesmus obliquus Mitochondrial (22)",
        "23": "Thraustochytrium Mitochondrial (23)",
        "24": "Pterobranchia Mitochondrial (24)",
        "25": "Candidate Division SR1 and Gracilibacteria (25)",
        "26": "Pachysolen tannophilus Nuclear (26)",
        "27": "Karyorelict Nuclear (27)",
        "28": "Condylostoma Nuclear (28)",
        "29": "Mesodinium Nuclear (29)",
        "30": "Peritrich Nuclear (30)",
        "31": "Blastocrithidia Nuclear (31)"
    },
    "matrices": ["PAM30", "PAM70", "BLOSUM90", "BLOSUM80", "BLOSUM62", "BLOSUM50", "BLOSUM45", "PAM250"],
    "gap_costs": {
        "protein": [
            {"value": "11 2", "description": "Existence: 11 Extension: 2"},
            {"value": "10 2", "description": "Existence: 10 Extension: 2"},
            {"value": "9 2", "description": "Existence: 9 Extension: 2"},
            {"value": "8 2", "description": "Existence: 8 Extension: 2"},
            {"value": "7 2", "description": "Existence: 7 Extension: 2"},
            {"value": "6 2", "description": "Existence: 6 Extension: 2"},
            {"value": "13 1", "description": "Existence: 13 Extension: 1"},
            {"value": "12 1", "description": "Existence: 12 Extension: 1"},
            {"value": "11 1", "description": "Existence: 11 Extension: 1"},
            {"value": "10 1", "description": "Existence: 10 Extension: 1"},
            {"value": "9 1", "description": "Existence: 9 Extension: 1"}
        ],
        "nucleotide": {
            "2,-3": [
                {"value": "5 2", "description": "Existence: 5 Extension: 2"},
                {"value": "4 4", "description": "Existence: 4 Extension: 4"},
                {"value": "2 4", "description": "Existence: 2 Extension: 4"},
                {"value": "0 4", "description": "Existence: 0 Extension: 4"}
            ],
            "1,-2": [
                {"value": "2 1", "description": "Existence: 2 Extension: 1"},
                {"value": "1 1", "description": "Existence: 1 Extension: 1"},
                {"value": "0 1", "description": "Existence: 0 Extension: 1"}
            ],
            "1,-3": [
                {"value": "2 1", "description": "Existence: 2 Extension: 1"},
                {"value": "1 1", "description": "Existence: 1 Extension: 1"},
                {"value": "0 1", "description": "Existence: 0 Extension: 1"}
            ],
            "1,-4": [
                {"value": "1 1", "description": "Existence: 1 Extension: 1"},
                {"value": "0 1", "description": "Existence: 0 Extension: 1"}
            ],
            "4,-5": [
                {"value": "12 8", "description": "Existence: 12 Extension: 8"},
                {"value": "8 6", "description": "Existence: 8 Extension: 6"},
                {"value": "6 5", "description": "Existence: 6 Extension: 5"},
                {"value": "4 4", "description": "Existence: 4 Extension: 4"}
            ],
            "1,-1": [
                {"value": "3 1", "description": "Existence: 3 Extension: 1"},
                {"value": "2 1", "description": "Existence: 2 Extension: 1"},
                {"value": "1 1", "description": "Existence: 1 Extension: 1"}
            ]
        }
    },
    "compositional_adjustments": [
        {"value": "0", "description": "No adjustment"},
        {"value": "1", "description": "Composition-based statistics"},
        {"value": "2", "description": "Conditional compositional score matrix adjustment"},
        {"value": "3", "description": "Universal compositional score matrix adjustment"}
    ],
    "match_scores": ["1,-2", "1,-3", "1,-4", "2,-3", "4,-5", "1,-1"],
    "databases": {
        "protein": [
            ("nr", "Non-redundant protein sequences (nr)"),
            ("refseq_select_prot", "RefSeq Select proteins (refseq_select)"),
            ("refseq_protein", "Reference proteins (refseq_protein)"),
            ("SMARTBLAST/landmark", "Model Organisms (landmark)"),
            ("swissprot", "UniProtKB/Swiss-Prot(swissprot)"),
            ("pataa", "Patented protein sequences(pataa)"),
            ("pdb", "Protein Data Bank proteins(pdb)"),
            ("env_nr", "Metagenomic proteins(env_nr)"),
            ("tsa_nr", "Transcriptome Shotgun Assembly proteins (tsa_nr)"),
            ("nr_cluster_seq", "ClusteredNR (nr_cluster_seq)")
        ],
        "nucleotide": [
            ("nt", "Nucleotide collection (nt)"),
            ("refseq_rna", "Reference RNA sequences (refseq_rna)"),
            ("refseq_genomic", "Reference genomic sequences (refseq_genomic)"),
            ("chromosome", "Chromosome sequences (chromosome)"),
            ("wgs", "Whole-genome shotgun reads (wgs)"),
            ("est", "Expressed sequence tags (est)"),
            ("gss", "Genome survey sequences (gss)"),
            ("htgs", "High throughput genomic sequences (htgs)"),
            ("pat", "Patent sequences (pat)"),
            ("pdb", "Protein Data Bank sequences (pdb)"),
            ("month", "Month sequences (month)"),
            ("env_nt", "Environmental samples (env_nt)")
        ]
    }
}


def create_blastp_tab(parent):
    """Create and return the BLASTP tab widget with all controls."""
    
    # Main widget for the tab
    blastp_widget = QWidget()
    main_layout = QVBoxLayout(blastp_widget)
    main_layout.setContentsMargins(10, 10, 10, 10)
    main_layout.setSpacing(10)
    
    # Initialize BLAST worker thread
    parent.blast_worker = None
    
    # Create splitter for input and results
    splitter = QSplitter(Qt.Vertical)
    main_layout.addWidget(splitter)
    
    # Top section - Input controls
    input_widget = QWidget()
    input_layout = QVBoxLayout(input_widget)
    
    # Online BLAST notice
    notice_label = QLabel(
        "🌐 <b>Online BLAST</b> - Uses NCBI's web services (no local installation required)"
    )
    notice_label.setStyleSheet(
        "background-color: #e8f4fd; border: 1px solid #bee5eb; "
        "border-radius: 4px; padding: 8px; margin: 5px;"
    )
    notice_label.setWordWrap(True)
    input_layout.addWidget(notice_label)
    
    # Query sequence section
    query_group = QGroupBox("Query Sequence")
    query_layout = QVBoxLayout(query_group)
    
    # File input
    file_layout = QHBoxLayout()
    file_button = QPushButton("Choose File")
    file_button.setToolTip("Select a FASTA file containing your protein sequence")
    parent.blastp_file_label = QLabel("No file selected")
    file_button.clicked.connect(lambda: open_blast_file(parent))
    file_layout.addWidget(file_button)
    file_layout.addWidget(parent.blastp_file_label)
    file_layout.addStretch()
    query_layout.addLayout(file_layout)
    
    # Query text area
    parent.blastp_query_input = QTextEdit()
    parent.blastp_query_input.setPlaceholderText(
        "Enter protein sequence in FASTA format or paste sequence directly...\n\n"
        "Example:\n"
        ">My_Protein\n"
        "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
    )
    parent.blastp_query_input.setMaximumHeight(120)
    parent.blastp_query_input.setFont(QFont("Courier", 10))
    query_layout.addWidget(parent.blastp_query_input)
    
    input_layout.addWidget(query_group)
    
    # Database selection
    db_group = QGroupBox("Database")
    db_layout = QVBoxLayout(db_group)
    
    parent.blastp_database_combo = QComboBox()
    for db_id, db_name in BLAST_CONFIG["databases"]["protein"]:
        parent.blastp_database_combo.addItem(f"{db_id} - {db_name}")
    parent.blastp_database_combo.setToolTip("Select the protein database to search")
    db_layout.addWidget(parent.blastp_database_combo)
    
    input_layout.addWidget(db_group)
    
    # Algorithm parameters
    algo_group = QGroupBox("Search Parameters")
    algo_layout = QFormLayout(algo_group)
    
    # Max target sequences
    parent.blastp_max_targets = QComboBox()
    parent.blastp_max_targets.addItems(["10", "50", "100", "250", "500"])
    parent.blastp_max_targets.setCurrentText("100")
    parent.blastp_max_targets.setToolTip("Maximum number of aligned sequences to display")
    algo_layout.addRow("Max target sequences:", parent.blastp_max_targets)
    
    # Expect threshold
    parent.blastp_expect_threshold = QLineEdit("0.05")
    parent.blastp_expect_threshold.setToolTip("Expected number of chance matches (lower = more stringent)")
    algo_layout.addRow("E-value threshold:", parent.blastp_expect_threshold)
    
    # Word size
    parent.blastp_word_size = QComboBox()
    parent.blastp_word_size.addItems(["2", "3", "6"])
    parent.blastp_word_size.setCurrentText("3")
    parent.blastp_word_size.setToolTip("Size of word for initial matches")
    algo_layout.addRow("Word size:", parent.blastp_word_size)
    
    # Scoring matrix
    parent.blastp_matrix = QComboBox()
    parent.blastp_matrix.addItems(BLAST_CONFIG["matrices"])
    parent.blastp_matrix.setCurrentText("BLOSUM62")
    parent.blastp_matrix.setToolTip("Scoring matrix for amino acid substitutions")
    algo_layout.addRow("Matrix:", parent.blastp_matrix)
    
    # Gap costs
    parent.blastp_gap_costs = QComboBox()
    for gap_cost in BLAST_CONFIG["gap_costs"]["protein"]:
        parent.blastp_gap_costs.addItem(gap_cost["description"])
    parent.blastp_gap_costs.setCurrentText("Existence: 11 Extension: 1")
    parent.blastp_gap_costs.setToolTip("Gap opening and extension penalties")
    algo_layout.addRow("Gap costs:", parent.blastp_gap_costs)

    # Compositional Adjustments
    parent.blastp_comp_adjust = QComboBox()
    for adj in BLAST_CONFIG["compositional_adjustments"]:
        parent.blastp_comp_adjust.addItem(f"{adj['value']} - {adj['description']}")
    parent.blastp_comp_adjust.setCurrentText("2 - Conditional compositional score matrix adjustment")
    parent.blastp_comp_adjust.setToolTip("Method for compositional adjustment of scoring")
    algo_layout.addRow("Compositional Adjustments:", parent.blastp_comp_adjust)

    # Word Score Threshold
    parent.blastp_word_threshold = QLineEdit("11")
    parent.blastp_word_threshold.setToolTip("Word score threshold for initial matches")
    algo_layout.addRow("Word score threshold:", parent.blastp_word_threshold)

    # Query Genetic Code
    parent.blastp_genetic_code = QComboBox()
    for code, desc in BLAST_CONFIG["genetic_codes"].items():
        parent.blastp_genetic_code.addItem(f"{code} - {desc}")
    parent.blastp_genetic_code.setCurrentText("1 - Standard (1)")
    parent.blastp_genetic_code.setToolTip("Genetic code to use for query sequence (if translated)")
    algo_layout.addRow("Query Genetic Code:", parent.blastp_genetic_code)

    # Ungapped Alignment
    parent.blastp_ungapped_alignment = QCheckBox("Perform ungapped alignment only")
    parent.blastp_ungapped_alignment.setChecked(False)
    parent.blastp_ungapped_alignment.setToolTip("Only perform ungapped alignments")
    algo_layout.addRow("", parent.blastp_ungapped_alignment)
    
    input_layout.addWidget(algo_group)
    
    # Filters and organism
    filters_group = QGroupBox("Filters and Organism")
    filters_layout = QFormLayout(filters_group)
    
    parent.blastp_low_complexity = QCheckBox("Filter low complexity regions")
    parent.blastp_low_complexity.setChecked(True)
    parent.blastp_low_complexity.setToolTip("Filter out low complexity regions")
    filters_layout.addRow("", parent.blastp_low_complexity)
    
    parent.blastp_organism_input = QLineEdit()
    parent.blastp_organism_input.setPlaceholderText("e.g., Homo sapiens, Escherichia coli")
    parent.blastp_organism_input.setToolTip("Restrict search to specific organism (optional)")
    filters_layout.addRow("Organism:", parent.blastp_organism_input)
    
    input_layout.addWidget(filters_group)
    
    # Run button and progress
    run_layout = QVBoxLayout()
    
    button_layout = QHBoxLayout()
    parent.blastp_run_button = QPushButton("🚀 Run Online BLAST Search")
    parent.blastp_run_button.setStyleSheet("""
        QPushButton {
            background-color: #4CAF50;
            color: white;
            border: none;
            padding: 12px 24px;
            font-size: 14px;
            font-weight: bold;
            border-radius: 6px;
        }
        QPushButton:hover {
            background-color: #45a049;
        }
        QPushButton:pressed {
            background-color: #3d8b40;
        }
        QPushButton:disabled {
            background-color: #cccccc;
            color: #666666;
        }
    """)
    parent.blastp_run_button.clicked.connect(lambda: run_online_blast_search(parent))
    
    # Cancel button (initially hidden)
    parent.blastp_cancel_button = QPushButton("Cancel Search")
    parent.blastp_cancel_button.setStyleSheet("""
        QPushButton {
            background-color: #f44336;
            color: white;
            border: none;
            padding: 12px 24px;
            font-size: 14px;
            font-weight: bold;
            border-radius: 6px;
        }
        QPushButton:hover {
            background-color: #da190b;
        }
    """)
    parent.blastp_cancel_button.clicked.connect(lambda: cancel_blast_search(parent))
    parent.blastp_cancel_button.hide()
    
    button_layout.addStretch()
    button_layout.addWidget(parent.blastp_run_button)
    button_layout.addWidget(parent.blastp_cancel_button)
    button_layout.addStretch()
    run_layout.addLayout(button_layout)
    
    # Progress bar
    parent.blastp_progress_bar = QProgressBar()
    parent.blastp_progress_bar.setVisible(False)
    parent.blastp_progress_bar.setRange(0, 0)  # Indeterminate progress
    run_layout.addWidget(parent.blastp_progress_bar)
    
    input_layout.addLayout(run_layout)
    
    # Add scroll area for input controls
    scroll_area = QScrollArea()
    scroll_area.setWidget(input_widget)
    scroll_area.setWidgetResizable(True)
    scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
    scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
    
    splitter.addWidget(scroll_area)
    
    # Results section
    results_group = QGroupBox("BLAST Results")
    results_layout = QVBoxLayout(results_group)
    
    # Results toolbar
    results_toolbar = QHBoxLayout()
    
    save_results_button = QPushButton("💾 Save Results")
    save_results_button.setToolTip("Save BLAST results to file")
    save_results_button.clicked.connect(lambda: save_blast_results(parent))
    results_toolbar.addWidget(save_results_button)
    
    clear_results_button = QPushButton("🗑️ Clear Results")
    clear_results_button.setToolTip("Clear the results display")
    clear_results_button.clicked.connect(lambda: parent.blastp_results_display.clear())
    results_toolbar.addWidget(clear_results_button)
    
    view_ncbi_button = QPushButton("🌐 View on NCBI")
    view_ncbi_button.setToolTip("Open NCBI BLAST with current query")
    view_ncbi_button.clicked.connect(lambda: open_ncbi_blast(parent))
    results_toolbar.addWidget(view_ncbi_button)
    
    help_button = QPushButton("❓ Help")
    help_button.setToolTip("Show BLAST parameter help")
    help_button.clicked.connect(lambda: show_blast_help(parent))
    results_toolbar.addWidget(help_button)
    
    results_toolbar.addStretch()
    results_layout.addLayout(results_toolbar)
    
    # Results display
    parent.blastp_results_display = QTextEdit()
    parent.blastp_results_display.setReadOnly(True)
    parent.blastp_results_display.setFont(QFont("Courier", 10))
    parent.blastp_results_display.setPlaceholderText(
        "BLAST results will appear here...\n\n"
        "• Enter a protein sequence above\n"
        "• Click 'Run Online BLAST Search' to begin\n"
        "• Results typically take 30 seconds to 2 minutes\n"
        "• No local BLAST+ installation required!"
    )
    results_layout.addWidget(parent.blastp_results_display)
    
    splitter.addWidget(results_group)
    
    # Set splitter proportions (60% input, 40% results)
    splitter.setSizes([600, 400])
    
    return blastp_widget


def run_online_blast_search(parent):
    """Run online BLAST search using NCBI's web services."""
    
    query = parent.blastp_query_input.toPlainText().strip()
    if not query:
        parent.blastp_results_display.setText("❌ Error: Query sequence is empty.")
        return
    
    # Validate sequence
    is_valid, message = validate_sequence(query, "protein")
    if not is_valid:
        parent.blastp_results_display.setText(f"❌ Error: {message}")
        return
    
    # Clean sequence (remove FASTA headers)
    clean_sequence = clean_fasta_sequence(query)
    
    # Get database
    db_text = parent.blastp_database_combo.currentText()
    database = db_text.split(" - ")[0]  # Extract database name
    
    # Get parameters
    parameters = {
        'evalue': parent.blastp_expect_threshold.text().strip() or "0.05",
        'max_target_seqs': parent.blastp_max_targets.currentText(),
        'word_size': parent.blastp_word_size.currentText(),
        'matrix': parent.blastp_matrix.currentText(),
        'low_complexity': parent.blastp_low_complexity.isChecked(),
        'organism': parent.blastp_organism_input.text().strip(),
        'comp_adjust': parent.blastp_comp_adjust.currentText().split(" - ")[0],
        'word_threshold': parent.blastp_word_threshold.text().strip() or "11",
        'genetic_code': parent.blastp_genetic_code.currentText().split(" - ")[0],
        'ungapped_alignment': parent.blastp_ungapped_alignment.isChecked()
    }
    
    # Parse gap costs
    gap_text = parent.blastp_gap_costs.currentText()
    for gap_cost in BLAST_CONFIG["gap_costs"]["protein"]:
        if gap_cost["description"] == gap_text:
            gap_values = gap_cost["value"].split(' ')
            if len(gap_values) >= 2:
                parameters['gapopen'] = gap_values[0]
                parameters['gapextend'] = gap_values[1]
            break
    
    # Update UI
    parent.blastp_run_button.setEnabled(False)
    parent.blastp_cancel_button.show()
    parent.blastp_progress_bar.setVisible(True)
    parent.blastp_results_display.setText("🚀 Starting online BLAST search...")
    parent.statusBar().showMessage("Running online BLAST search...")
    
    # Create and start worker thread
    parent.blast_worker = OnlineBlastWorker(clean_sequence, database, "blastp", parameters)
    parent.blast_worker.finished.connect(lambda result: on_blast_finished(parent, result))
    parent.blast_worker.error.connect(lambda error: on_blast_error(parent, error))
    parent.blast_worker.progress.connect(lambda msg: on_blast_progress(parent, msg))
    parent.blast_worker.start()


def cancel_blast_search(parent):
    """Cancel the running BLAST search."""
    
    if parent.blast_worker and parent.blast_worker.isRunning():
        parent.blast_worker.cancel_search()
        parent.blast_worker.wait()
        parent.blastp_results_display.setText("⏹️ BLAST search cancelled.")
        parent.statusBar().showMessage("BLAST search cancelled.")
    
    # Reset UI
    parent.blastp_run_button.setEnabled(True)
    parent.blastp_cancel_button.hide()
    parent.blastp_progress_bar.setVisible(False)


def on_blast_finished(parent, result):
    """Handle successful BLAST completion."""
    
    # Get job title and request ID (for blastp)
    job_title = None  # blastp doesn't have job title in the current interface
    request_id = parent.blast_worker.request_id if parent.blast_worker else None
    
    formatted_result = format_blast_output(result, "blastp", job_title, request_id)
    parent.blastp_results_display.setText(formatted_result)
    parent.statusBar().showMessage("✅ Online BLAST search completed successfully.")
    
    # Reset UI
    parent.blastp_run_button.setEnabled(True)
    parent.blastp_cancel_button.hide()
    parent.blastp_progress_bar.setVisible(False)


def on_blast_error(parent, error):
    """Handle BLAST error."""
    
    parent.blastp_results_display.setText(f"❌ BLAST search failed:\n\n{error}")
    parent.statusBar().showMessage("❌ Online BLAST search failed.")
    
    # Reset UI
    parent.blastp_run_button.setEnabled(True)
    parent.blastp_cancel_button.hide()
    parent.blastp_progress_bar.setVisible(False)


def on_blast_progress(parent, message):
    """Handle BLAST progress updates."""
    
    parent.blastp_results_display.setText(f"⏳ {message}")
    parent.statusBar().showMessage(message)


def validate_sequence(sequence, sequence_type="protein"):
    """Validate a biological sequence or accession number."""
    from .accession_validator import validate_sequence_or_accession
    return validate_sequence_or_accession(sequence, sequence_type)


def clean_fasta_sequence(sequence):
    """Clean FASTA sequence by removing headers and formatting."""
    from .accession_validator import clean_fasta_sequence_or_accession, is_accession_number
    return clean_fasta_sequence_or_accession(sequence)


def format_blast_output(result_data, program_type="blast", job_title=None, request_id=None):
    """Format BLAST output using enhanced parser for comprehensive tables."""
    
    # Import the enhanced parser
    try:
        from .blast_results_parser import enhanced_format_blast_output
        
        # Handle the new tuple format (format_type, content) or just content string
        if isinstance(result_data, tuple) and len(result_data) == 2:
            format_type, raw_output = result_data
        else:
            # Backward compatibility - assume it's just the content string
            format_type, raw_output = 'Text', result_data
        
        if not raw_output or not raw_output.strip():
            return "❌ No results found or empty output."
        
        # Use enhanced parser for better formatting
        return enhanced_format_blast_output(raw_output, program_type, job_title, request_id)
        
    except ImportError:
        # Fall back to basic formatting if enhanced parser is not available
        return format_blast_output_basic(result_data)


def format_blast_output_basic(result_data):
    """Basic BLAST output formatter (fallback)."""
    
    # Handle the new tuple format or just content string
    if isinstance(result_data, tuple) and len(result_data) == 2:
        format_type, raw_output = result_data
    else:
        raw_output = result_data
    
    if not raw_output or not raw_output.strip():
        return "❌ No results found or empty output."
    
    # Add header
    formatted = "🧬 BLAST Search Results\n"
    formatted += "=" * 50 + "\n\n"
    
    # Process the output
    lines = raw_output.split('\n')
    in_alignments = False
    hit_count = 0
    
    for line in lines:
        # Skip empty lines at start
        if not line.strip() and not formatted.endswith('\n\n'):
            continue
        
        # Detect sections
        if "Sequences producing significant alignments:" in line:
            formatted += "📊 " + line + "\n"
            formatted += "-" * 50 + "\n"
            in_alignments = True
            continue
        elif line.startswith(">"):
            hit_count += 1
            formatted += f"\n🎯 Hit #{hit_count}: {line[1:]}\n"
            continue
        elif "Length=" in line:
            formatted += f"   📏 {line.strip()}\n"
            continue
        elif "Score =" in line and "Expect =" in line:
            formatted += f"   📈 {line.strip()}\n"
            continue
        elif "Identities =" in line:
            formatted += f"   🔍 {line.strip()}\n"
            continue
        elif line.strip().startswith("Query") or line.strip().startswith("Sbjct"):
            formatted += f"   {line}\n"
            continue
        
        # Add other lines
        formatted += line + "\n"
    
    # Add summary if no hits
    if hit_count == 0 and "No hits found" not in raw_output:
        if "No significant similarity found" in raw_output or "No hits found" in raw_output:
            formatted += "\n❌ No significant hits found.\n\n"
            formatted += "💡 Try adjusting your search parameters:\n"
            formatted += "   • Increase the E-value threshold\n"
            formatted += "   • Use a different scoring matrix\n"
            formatted += "   • Check your query sequence\n"
            formatted += "   • Try a different database\n"
    
    return formatted


def save_blast_results(parent):
    """Save BLAST results to a file."""
    
    if not hasattr(parent, 'blastp_results_display') or not parent.blastp_results_display.toPlainText().strip():
        QMessageBox.warning(parent, "No Results", "No BLAST results to save.")
        return
    
    # Generate default filename with timestamp
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    default_filename = f"online_blast_results_{timestamp}.txt"
    
    file_path, _ = QFileDialog.getSaveFileName(
        parent, 
        "Save BLAST Results", 
        default_filename, 
        "Text Files (*.txt);;All Files (*)"
    )
    
    if file_path:
        try:
            content = parent.blastp_results_display.toPlainText()
            
            # Add header information
            header = f"# Online BLAST Results saved on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            header += f"# Query: {parent.blastp_query_input.toPlainText()[:100]}...\n"
            header += f"# Database: {parent.blastp_database_combo.currentText()}\n"
            if hasattr(parent, 'blast_worker') and parent.blast_worker and parent.blast_worker.request_id:
                header += f"# Request ID: {parent.blast_worker.request_id}\n"
            header += "# Generated by PicoMol using NCBI's online BLAST service\n"
            header += "# " + "="*60 + "\n\n"
            
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(header + content)
            
            QMessageBox.information(parent, "Results Saved", f"BLAST results saved to:\n{file_path}")
            parent.statusBar().showMessage(f"Results saved to {os.path.basename(file_path)}")
        except Exception as e:
            QMessageBox.critical(parent, "Save Error", f"Failed to save results:\n{str(e)}")


def open_ncbi_blast(parent):
    """Open NCBI BLAST website with current query."""
    
    if not hasattr(parent, 'blastp_query_input'):
        return
        
    query = parent.blastp_query_input.toPlainText().strip()
    if not query:
        QMessageBox.warning(parent, "No Query", "Please enter a sequence to search.")
        return
    
    # Clean sequence for URL
    clean_query = clean_fasta_sequence(query)
    
    if len(clean_query) > 8000:  # NCBI has URL length limits
        QMessageBox.warning(
            parent, 
            "Sequence Too Long", 
            "The sequence is too long to pass via URL. Please copy and paste it manually into NCBI BLAST."
        )
        # Open NCBI BLAST without query
        ncbi_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch"
    else:
        # Encode the query for URL
        encoded_query = quote(clean_query)
        ncbi_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&QUERY={encoded_query}"
    
    try:
        webbrowser.open_new_tab(ncbi_url)
        parent.statusBar().showMessage("Opening NCBI BLAST in web browser...")
    except Exception as e:
        QMessageBox.warning(parent, "Browser Error", f"Could not open web browser: {str(e)}")


def open_blast_file(parent):
    """Open a FASTA file and load it into the query input."""
    
    file_path, _ = QFileDialog.getOpenFileName(
        parent, "Open Sequence File", "", "FASTA Files (*.fasta *.fa *.fas);;All Files (*)"
    )
    
    if file_path:
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            parent.blastp_file_label.setText(os.path.basename(file_path))
            parent.blastp_query_input.setText(content)
            
            # Validate the loaded sequence
            is_valid, message = validate_sequence(content, "protein")
            if not is_valid:
                QMessageBox.warning(
                    parent, 
                    "Invalid Sequence", 
                    f"The loaded file contains invalid sequence data:\n\n{message}"
                )
        except Exception as e:
            QMessageBox.critical(
                parent, 
                "File Error", 
                f"Could not read the selected file:\n\n{str(e)}"
            )


def show_blast_help(parent):
    """Show BLAST help dialog."""
    
    help_dialog = QDialog(parent)
    help_dialog.setWindowTitle("Online BLAST Help")
    help_dialog.setMinimumSize(600, 500)
    
    layout = QVBoxLayout(help_dialog)
    
    help_text = QTextEdit()
    help_text.setReadOnly(True)
    help_text.setHtml("""
    <h2>🧬 Online BLAST Help</h2>
    
    <h3>📖 Overview</h3>
    <p>This tool uses NCBI's online BLAST service to search protein databases. 
    No local BLAST+ installation is required!</p>
    
    <h3>⚙️ Parameters Guide</h3>
    
    <h4>E-value Threshold</h4>
    <ul>
        <li><b>10</b> - Default, finds most matches</li>
        <li><b>0.001</b> - More stringent, fewer false positives</li>
        <li><b>1e-10</b> - Very stringent, only close matches</li>
    </ul>
    
    <h4>Word Size</h4>
    <ul>
        <li><b>3</b> - Default for proteins, good balance</li>
        <li><b>2</b> - More sensitive, slower</li>
        <li><b>6</b> - Less sensitive, faster</li>
    </ul>
    
    <h4>Scoring Matrix</h4>
    <ul>
        <li><b>BLOSUM62</b> - Default, good for most searches</li>
        <li><b>BLOSUM45</b> - More sensitive for distant relationships</li>
        <li><b>BLOSUM80</b> - Less sensitive, for close relationships</li>
    </ul>
    
    <h4>Gap Costs</h4>
    <ul>
        <li><b>Existence 11, Extension 1</b> - Default, common for protein alignments</li>
        <li>Adjust to control penalties for opening and extending gaps.</li>
    </ul>

    <h4>Compositional Adjustments</h4>
    <ul>
        <li><b>No adjustment</b> - For sequences with typical amino acid composition.</li>
        <li><b>Composition-based statistics</b> - Adjusts scores for sequences with unusual compositions.</li>
        <li><b>Conditional compositional score matrix adjustment</b> - Default, suitable for most protein searches.</li>
        <li><b>Universal compositional score matrix adjustment</b> - More aggressive adjustment.</li>
    </ul>

    <h4>Word Score Threshold</h4>
    <ul>
        <li><b>11</b> - Default, determines how significant a match must be to be considered a "word".</li>
        <li>Lower values increase sensitivity but also search time.</li>
    </ul>

    <h4>Query Genetic Code</h4>
    <ul>
        <li>Specifies the genetic code to use for translating the query sequence.</li>
        <li>Important for accurate searches if the query is a nucleotide sequence that needs translation (though this is blastp, it's a general BLAST option).</li>
    </ul>

    <h4>Ungapped Alignment Only</h4>
    <ul>
        <li>If checked, only ungapped alignments will be reported.</li>
        <li>Useful for very fast searches or specific types of analyses where gaps are not expected.</li>
    </ul>
    
    <h4>Databases</h4>
    <ul>
        <li><b>nr</b> - Non-redundant protein sequences (largest)</li>
        <li><b>refseq_protein</b> - High-quality reference proteins</li>
        <li><b>swissprot</b> - Manually annotated, high quality</li>
        <li><b>pdb</b> - Protein structures from PDB</li>
    </ul>
    
    <h3>⏱️ Search Times</h3>
    <p>Online searches typically take 30 seconds to 2 minutes depending on:</p>
    <ul>
        <li>Sequence length</li>
        <li>Database size</li>
        <li>NCBI server load</li>
        <li>Search parameters</li>
    </ul>
    
    <h3>📏 Limitations</h3>
    <ul>
        <li>Maximum sequence length: 20,000 amino acids</li>
        <li>Requires internet connection</li>
        <li>Subject to NCBI usage policies</li>
    </ul>
    
    <h3>💡 Tips</h3>
    <ul>
        <li>Start with default parameters</li>
        <li>Use organism filter to narrow results</li>
        <li>Save results before starting new searches</li>
        <li>Check sequence format if search fails</li>
    </ul>
    """)
    layout.addWidget(help_text)
    
    button_box = QDialogButtonBox(QDialogButtonBox.Ok)
    button_box.accepted.connect(help_dialog.accept)
    layout.addWidget(button_box)
    
    help_dialog.exec_()


# Placeholder functions for other BLAST types
def create_blastn_tab(parent):
    """Create and return the BLASTN tab widget with comprehensive options."""
    
    # Main widget for the tab
    blastn_widget = QWidget()
    main_layout = QVBoxLayout(blastn_widget)
    main_layout.setContentsMargins(10, 10, 10, 10)
    main_layout.setSpacing(10)
    
    # Initialize BLAST worker thread
    parent.blastn_worker = None
    
    # Create splitter for input and results
    splitter = QSplitter(Qt.Vertical)
    main_layout.addWidget(splitter)
    
    # Top section - Input controls
    input_widget = QWidget()
    input_layout = QVBoxLayout(input_widget)
    
    # Online BLAST notice
    notice_label = QLabel(
        "🌐 <b>Online BLASTN</b> - Search nucleotide databases using a nucleotide query<br>"
        "Uses NCBI's web services (no local installation required)"
    )
    notice_label.setStyleSheet(
        "background-color: #e8f4fd; border: 1px solid #bee5eb; "
        "border-radius: 4px; padding: 8px; margin: 5px;"
    )
    notice_label.setWordWrap(True)
    input_layout.addWidget(notice_label)
    
    # Query sequence section
    query_group = QGroupBox("Query Sequence")
    query_layout = QVBoxLayout(query_group)
    
    # File input
    file_layout = QHBoxLayout()
    file_button = QPushButton("Choose File")
    file_button.setToolTip("Select a FASTA file containing your nucleotide sequence")
    parent.blastn_file_label = QLabel("No file selected")
    file_button.clicked.connect(lambda: open_blast_file_generic(parent, 'blastn'))
    file_layout.addWidget(file_button)
    file_layout.addWidget(parent.blastn_file_label)
    file_layout.addStretch()
    query_layout.addLayout(file_layout)
    
    # Query text area
    parent.blastn_query_input = QTextEdit()
    parent.blastn_query_input.setPlaceholderText(
        "Enter nucleotide sequence in FASTA format or paste sequence directly...\n\n"
        "Example:\n"
        ">My_Nucleotide\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    )
    parent.blastn_query_input.setMaximumHeight(120)
    parent.blastn_query_input.setFont(QFont("Courier", 10))
    query_layout.addWidget(parent.blastn_query_input)
    
    input_layout.addWidget(query_group)
    
    # Algorithm selection
    algo_group = QGroupBox("Algorithm")
    algo_layout = QVBoxLayout(algo_group)
    
    parent.blastn_algorithm_combo = QComboBox()
    parent.blastn_algorithm_combo.addItems(["megaBlast", "discoMegablast", "blastn"])
    parent.blastn_algorithm_combo.setToolTip("Select the BLASTN algorithm variant")
    algo_layout.addWidget(parent.blastn_algorithm_combo)
    
    input_layout.addWidget(algo_group)
    
    # Database selection
    db_group = QGroupBox("Database")
    db_layout = QVBoxLayout(db_group)
    
    parent.blastn_database_combo = QComboBox()
    for db_id, db_name in BLAST_CONFIG["databases"]["nucleotide"]:
        parent.blastn_database_combo.addItem(f"{db_id} - {db_name}")
    parent.blastn_database_combo.setToolTip("Select the nucleotide database to search")
    db_layout.addWidget(parent.blastn_database_combo)
    
    input_layout.addWidget(db_group)
    
    # Parameters section
    params_group = QGroupBox("Search Parameters")
    params_layout = QFormLayout(params_group)
    
    # Max target sequences
    parent.blastn_max_targets = QComboBox()
    parent.blastn_max_targets.addItems(["10", "50", "100", "250", "500", "1000", "5000"])
    parent.blastn_max_targets.setCurrentText("100")
    parent.blastn_max_targets.setToolTip("Maximum number of aligned sequences to display")
    params_layout.addRow("Max target sequences:", parent.blastn_max_targets)
    
    # Expect threshold
    parent.blastn_expect_threshold = QLineEdit("0.05")
    parent.blastn_expect_threshold.setToolTip("Expected number of chance matches (lower = more stringent)")
    params_layout.addRow("E-value threshold:", parent.blastn_expect_threshold)
    
    # Word size
    parent.blastn_word_size = QComboBox()
    parent.blastn_word_size.addItems(["7", "11", "15", "20", "24", "28"])
    parent.blastn_word_size.setCurrentText("11")
    parent.blastn_word_size.setToolTip("Size of word for initial matches")
    params_layout.addRow("Word size:", parent.blastn_word_size)
    
    # Match/Mismatch scores
    parent.blastn_match_scores = QComboBox()
    for score in BLAST_CONFIG["match_scores"]:
        parent.blastn_match_scores.addItem(f"{score} (Match,Mismatch)")
    parent.blastn_match_scores.setToolTip("Reward and penalty for matching and mismatching bases")
    params_layout.addRow("Match/Mismatch scores:", parent.blastn_match_scores)
    
    # Gap costs (will be populated based on match/mismatch scores)
    parent.blastn_gap_costs = QComboBox()
    # Default to 2,-3 match/mismatch scores gap costs
    for gap_cost in BLAST_CONFIG["gap_costs"]["nucleotide"]["2,-3"]:
        parent.blastn_gap_costs.addItem(gap_cost["description"])
    parent.blastn_gap_costs.setCurrentText("Existence: 5 Extension: 2")
    parent.blastn_gap_costs.setToolTip("Gap opening and extension penalties")
    params_layout.addRow("Gap costs:", parent.blastn_gap_costs)
    
    # Function to update gap costs when match/mismatch scores change
    def update_gap_costs():
        """Update gap costs dropdown based on selected match/mismatch scores."""
        match_scores = parent.blastn_match_scores.currentText().split(" ")[0]
        parent.blastn_gap_costs.clear()
        
        if match_scores in BLAST_CONFIG["gap_costs"]["nucleotide"]:
            for gap_cost in BLAST_CONFIG["gap_costs"]["nucleotide"][match_scores]:
                parent.blastn_gap_costs.addItem(gap_cost["description"])
            # Set default to first option
            if parent.blastn_gap_costs.count() > 0:
                parent.blastn_gap_costs.setCurrentIndex(0)
    
    # Connect match/mismatch scores change to gap costs update
    parent.blastn_match_scores.currentTextChanged.connect(update_gap_costs)
    setattr(parent, 'blastn_update_gap_costs', update_gap_costs)
    
    input_layout.addWidget(params_group)
    
    # Filters and organism
    filters_group = QGroupBox("Filters and Organism")
    filters_layout = QFormLayout(filters_group)
    
    parent.blastn_low_complexity = QCheckBox("Filter low complexity regions")
    parent.blastn_low_complexity.setChecked(True)
    parent.blastn_low_complexity.setToolTip("Filter out low complexity regions")
    filters_layout.addRow("", parent.blastn_low_complexity)
    
    parent.blastn_organism_input = QLineEdit()
    parent.blastn_organism_input.setPlaceholderText("e.g., Homo sapiens, Escherichia coli")
    parent.blastn_organism_input.setToolTip("Restrict search to specific organism (optional)")
    filters_layout.addRow("Organism:", parent.blastn_organism_input)
    
    input_layout.addWidget(filters_group)
    
    # Run button and progress
    run_layout = QVBoxLayout()
    
    button_layout = QHBoxLayout()
    parent.blastn_run_button = QPushButton("🚀 Run Online BLASTN Search")
    parent.blastn_run_button.setStyleSheet("""
        QPushButton {
            background-color: #4CAF50;
            color: white;
            border: none;
            padding: 12px 24px;
            font-size: 14px;
            font-weight: bold;
            border-radius: 6px;
        }
        QPushButton:hover {
            background-color: #45a049;
        }
        QPushButton:pressed {
            background-color: #3d8b40;
        }
        QPushButton:disabled {
            background-color: #cccccc;
            color: #666666;
        }
    """)
    parent.blastn_run_button.clicked.connect(lambda: run_online_blast_search_generic(parent, 'blastn'))
    
    # Cancel button (initially hidden)
    parent.blastn_cancel_button = QPushButton("Cancel Search")
    parent.blastn_cancel_button.setStyleSheet("""
        QPushButton {
            background-color: #f44336;
            color: white;
            border: none;
            padding: 12px 24px;
            font-size: 14px;
            font-weight: bold;
            border-radius: 6px;
        }
        QPushButton:hover {
            background-color: #da190b;
        }
    """)
    parent.blastn_cancel_button.clicked.connect(lambda: cancel_blast_search_generic(parent, 'blastn'))
    parent.blastn_cancel_button.hide()
    
    button_layout.addStretch()
    button_layout.addWidget(parent.blastn_run_button)
    button_layout.addWidget(parent.blastn_cancel_button)
    button_layout.addStretch()
    run_layout.addLayout(button_layout)
    
    # Progress bar
    parent.blastn_progress_bar = QProgressBar()
    parent.blastn_progress_bar.setVisible(False)
    parent.blastn_progress_bar.setRange(0, 0)  # Indeterminate progress
    run_layout.addWidget(parent.blastn_progress_bar)
    
    input_layout.addLayout(run_layout)
    
    # Add scroll area for input controls
    scroll_area = QScrollArea()
    scroll_area.setWidget(input_widget)
    scroll_area.setWidgetResizable(True)
    scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
    scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
    
    splitter.addWidget(scroll_area)
    
    # Results section
    results_group = QGroupBox("BLASTN Results")
    results_layout = QVBoxLayout(results_group)
    
    # Results toolbar
    results_toolbar = QHBoxLayout()
    
    save_results_button = QPushButton("💾 Save Results")
    save_results_button.setToolTip("Save BLAST results to file")
    save_results_button.clicked.connect(lambda: save_blast_results_generic(parent, 'blastn'))
    results_toolbar.addWidget(save_results_button)
    
    clear_results_button = QPushButton("🗑️ Clear Results")
    clear_results_button.setToolTip("Clear the results display")
    clear_results_button.clicked.connect(lambda: parent.blastn_results_display.clear())
    results_toolbar.addWidget(clear_results_button)
    
    view_ncbi_button = QPushButton("🌐 View on NCBI")
    view_ncbi_button.setToolTip("Open NCBI BLAST with current query")
    view_ncbi_button.clicked.connect(lambda: open_ncbi_blast_generic(parent, 'blastn'))
    results_toolbar.addWidget(view_ncbi_button)
    
    help_button = QPushButton("❓ Help")
    help_button.setToolTip("Show BLAST parameter help")
    help_button.clicked.connect(lambda: show_blast_help_generic(parent, 'blastn'))
    results_toolbar.addWidget(help_button)
    
    results_toolbar.addStretch()
    results_layout.addLayout(results_toolbar)
    
    # Results display
    parent.blastn_results_display = QTextEdit()
    parent.blastn_results_display.setReadOnly(True)
    parent.blastn_results_display.setFont(QFont("Courier", 10))
    parent.blastn_results_display.setPlaceholderText(
        "BLASTN results will appear here...\n\n"
        "• Enter a nucleotide sequence above\n"
        "• Click 'Run Online BLASTN Search' to begin\n"
        "• Results typically take 30 seconds to 2 minutes\n"
        "• No local BLAST+ installation required!"
    )
    results_layout.addWidget(parent.blastn_results_display)
    
    splitter.addWidget(results_group)
    
    # Set splitter proportions (60% input, 40% results)
    splitter.setSizes([600, 400])
    
    return blastn_widget

def create_blastx_tab(parent):
    """Create and return the BLASTX tab widget."""
    return create_placeholder_tab("BLASTX", "translated nucleotide vs protein search")

def create_tblastn_tab(parent):
    """Create and return the TBLASTN tab widget."""
    return create_placeholder_tab("TBLASTN", "protein vs translated nucleotide search")

def create_tblastx_tab(parent):
    """Create and return the TBLASTX tab widget."""
    return create_placeholder_tab("TBLASTX", "translated nucleotide vs translated nucleotide search")

# Generic functions for all BLAST types
def run_online_blast_search_generic(parent, program_type):
    """Generic function to run online BLAST search for any program type."""
    
    query_input = getattr(parent, f'{program_type}_query_input')
    query = query_input.toPlainText().strip()
    if not query:
        results_display = getattr(parent, f'{program_type}_results_display')
        results_display.setText("❌ Error: Query sequence is empty.")
        return
    
    # Validate sequence
    sequence_type = "protein" if program_type in ["blastp", "tblastn"] else "nucleotide"
    is_valid, message = validate_sequence(query, sequence_type)
    if not is_valid:
        results_display = getattr(parent, f'{program_type}_results_display')
        results_display.setText(f"❌ Error: {message}")
        return
    
    # Clean sequence (remove FASTA headers)
    clean_sequence = clean_fasta_sequence(query)
    
    # Get database
    database_combo = getattr(parent, f'{program_type}_database_combo')
    db_text = database_combo.currentText()
    database = db_text.split(" - ")[0]  # Extract database name
    
    # Get parameters
    expect_threshold = getattr(parent, f'{program_type}_expect_threshold')
    parameters = {
        'evalue': expect_threshold.text().strip() or "0.05",
    }
    
    # Add common parameters
    if hasattr(parent, f'{program_type}_max_targets'):
        max_targets = getattr(parent, f'{program_type}_max_targets')
        parameters['max_target_seqs'] = max_targets.currentText()
    
    if hasattr(parent, f'{program_type}_word_size'):
        word_size = getattr(parent, f'{program_type}_word_size')
        parameters['word_size'] = word_size.currentText()
    
    if hasattr(parent, f'{program_type}_matrix'):
        matrix = getattr(parent, f'{program_type}_matrix')
        parameters['matrix'] = matrix.currentText()
    
    if hasattr(parent, f'{program_type}_match_scores'):
        match_scores = getattr(parent, f'{program_type}_match_scores')
        parameters['match_scores'] = match_scores.currentText().split(" ")[0]
    
    if hasattr(parent, f'{program_type}_gap_costs'):
        gap_costs = getattr(parent, f'{program_type}_gap_costs')
        gap_text = gap_costs.currentText()
        
        # Determine which gap costs list to use
        if program_type in ["blastp", "blastx", "tblastn"]:
            gap_costs_list = BLAST_CONFIG["gap_costs"]["protein"]
        else:
            # For nucleotide searches (including tblastx), get the match/mismatch scores first
            match_scores = parameters.get('match_scores', '2,-3')
            if match_scores in BLAST_CONFIG["gap_costs"]["nucleotide"]:
                gap_costs_list = BLAST_CONFIG["gap_costs"]["nucleotide"][match_scores]
            else:
                gap_costs_list = BLAST_CONFIG["gap_costs"]["nucleotide"]["2,-3"]  # Default
        
        for gap_cost in gap_costs_list:
            if gap_cost["description"] == gap_text:
                gap_values = gap_cost["value"].split(' ')
                if len(gap_values) >= 2:
                    parameters['gapopen'] = gap_values[0]
                    parameters['gapextend'] = gap_values[1]
                break
    
    if hasattr(parent, f'{program_type}_comp_adjust'):
        comp_adjust = getattr(parent, f'{program_type}_comp_adjust')
        parameters['comp_adjust'] = comp_adjust.currentText().split(" - ")[0]
    
    if hasattr(parent, f'{program_type}_genetic_code'):
        genetic_code = getattr(parent, f'{program_type}_genetic_code')
        parameters['genetic_code'] = genetic_code.currentText().split(" - ")[0]
    
    if hasattr(parent, f'{program_type}_algorithm_combo'):
        algorithm_combo = getattr(parent, f'{program_type}_algorithm_combo')
        parameters['algorithm'] = algorithm_combo.currentText()
    
    if hasattr(parent, f'{program_type}_low_complexity'):
        low_complexity = getattr(parent, f'{program_type}_low_complexity')
        parameters['low_complexity'] = low_complexity.isChecked()
    
    if hasattr(parent, f'{program_type}_organism_input'):
        organism_input = getattr(parent, f'{program_type}_organism_input')
        parameters['organism'] = organism_input.text().strip()
    
    # Update UI
    run_button = getattr(parent, f'{program_type}_run_button')
    cancel_button = getattr(parent, f'{program_type}_cancel_button')
    progress_bar = getattr(parent, f'{program_type}_progress_bar')
    results_display = getattr(parent, f'{program_type}_results_display')
    
    run_button.setEnabled(False)
    cancel_button.show()
    progress_bar.setVisible(True)
    results_display.setText("🚀 Starting online BLAST search...")
    parent.statusBar().showMessage(f"Running online {program_type.upper()} search...")
    
    # Create and start worker thread
    worker = OnlineBlastWorker(clean_sequence, database, program_type, parameters)
    setattr(parent, f'{program_type}_worker', worker)
    worker.finished.connect(lambda result: on_blast_finished_generic(parent, program_type, result))
    worker.error.connect(lambda error: on_blast_error_generic(parent, program_type, error))
    worker.progress.connect(lambda msg: on_blast_progress_generic(parent, program_type, msg))
    worker.start()


def cancel_blast_search_generic(parent, program_type):
    """Generic function to cancel BLAST search."""
    
    worker = getattr(parent, f'{program_type}_worker', None)
    if worker and worker.isRunning():
        worker.cancel_search()
        worker.wait()
        results_display = getattr(parent, f'{program_type}_results_display')
        results_display.setText("⏹️ BLAST search cancelled.")
        parent.statusBar().showMessage("BLAST search cancelled.")
    
    # Reset UI
    run_button = getattr(parent, f'{program_type}_run_button')
    cancel_button = getattr(parent, f'{program_type}_cancel_button')
    progress_bar = getattr(parent, f'{program_type}_progress_bar')
    
    run_button.setEnabled(True)
    cancel_button.hide()
    progress_bar.setVisible(False)


def on_blast_finished_generic(parent, program_type, result):
    """Generic function to handle successful BLAST completion."""
    
    # Get job title and request ID
    job_title = None  # Generic interface doesn't have job title
    worker = getattr(parent, f'{program_type}_worker', None)
    request_id = worker.request_id if worker else None
    
    formatted_result = format_blast_output(result, program_type, job_title, request_id)
    results_display = getattr(parent, f'{program_type}_results_display')
    results_display.setText(formatted_result)
    parent.statusBar().showMessage(f"✅ Online {program_type.upper()} search completed successfully.")
    
    # Reset UI
    run_button = getattr(parent, f'{program_type}_run_button')
    cancel_button = getattr(parent, f'{program_type}_cancel_button')
    progress_bar = getattr(parent, f'{program_type}_progress_bar')
    
    run_button.setEnabled(True)
    cancel_button.hide()
    progress_bar.setVisible(False)


def on_blast_error_generic(parent, program_type, error):
    """Generic function to handle BLAST error."""
    
    results_display = getattr(parent, f'{program_type}_results_display')
    results_display.setText(f"❌ BLAST search failed:\n\n{error}")
    parent.statusBar().showMessage(f"❌ Online {program_type.upper()} search failed.")
    
    # Reset UI
    run_button = getattr(parent, f'{program_type}_run_button')
    cancel_button = getattr(parent, f'{program_type}_cancel_button')
    progress_bar = getattr(parent, f'{program_type}_progress_bar')
    
    run_button.setEnabled(True)
    cancel_button.hide()
    progress_bar.setVisible(False)


def on_blast_progress_generic(parent, program_type, message):
    """Generic function to handle BLAST progress updates."""
    
    results_display = getattr(parent, f'{program_type}_results_display')
    results_display.setText(f"⏳ {message}")
    parent.statusBar().showMessage(message)


def save_blast_results_generic(parent, program_type):
    """Generic function to save BLAST results."""
    
    results_display = getattr(parent, f'{program_type}_results_display')
    if not results_display.toPlainText().strip():
        QMessageBox.warning(parent, "No Results", f"No {program_type.upper()} results to save.")
        return
    
    # Generate default filename with timestamp
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    default_filename = f"online_{program_type}_results_{timestamp}.txt"
    
    file_path, _ = QFileDialog.getSaveFileName(
        parent, 
        f"Save {program_type.upper()} Results", 
        default_filename, 
        "Text Files (*.txt);;All Files (*)"
    )
    
    if file_path:
        try:
            content = results_display.toPlainText()
            
            # Add header information
            query_input = getattr(parent, f'{program_type}_query_input')
            database_combo = getattr(parent, f'{program_type}_database_combo')
            
            header = f"# Online {program_type.upper()} Results saved on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            header += f"# Query: {query_input.toPlainText()[:100]}...\n"
            header += f"# Database: {database_combo.currentText()}\n"
            
            # Add job title if available
            job_title_input = getattr(parent, f'{program_type}_job_title', None)
            if job_title_input and job_title_input.text().strip():
                header += f"# Job Title: {job_title_input.text().strip()}\n"
            
            # Add request ID if available
            worker = getattr(parent, f'{program_type}_worker', None)
            if worker and worker.request_id:
                header += f"# Request ID: {worker.request_id}\n"
            
            header += "# Generated by PicoMol using NCBI's online BLAST service\n"
            header += "# " + "="*60 + "\n\n"
            
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(header + content)
            
            QMessageBox.information(parent, "Results Saved", f"{program_type.upper()} results saved to:\n{file_path}")
            parent.statusBar().showMessage(f"Results saved to {os.path.basename(file_path)}")
        except Exception as e:
            QMessageBox.critical(parent, "Save Error", f"Failed to save results:\n{str(e)}")


def open_ncbi_blast_generic(parent, program_type):
    """Generic function to open NCBI BLAST website."""
    
    query_input = getattr(parent, f'{program_type}_query_input')
    query = query_input.toPlainText().strip()
    if not query:
        QMessageBox.warning(parent, "No Query", "Please enter a sequence to search.")
        return
    
    # Clean sequence for URL
    clean_query = clean_fasta_sequence(query)
    
    if len(clean_query) > 8000:  # NCBI has URL length limits
        QMessageBox.warning(
            parent, 
            "Sequence Too Long", 
            "The sequence is too long to pass via URL. Please copy and paste it manually into NCBI BLAST."
        )
        # Open NCBI BLAST without query
        ncbi_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM={program_type}&PAGE_TYPE=BlastSearch"
    else:
        # Encode the query for URL
        encoded_query = quote(clean_query)
        ncbi_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM={program_type}&PAGE_TYPE=BlastSearch&QUERY={encoded_query}"
    
    try:
        webbrowser.open_new_tab(ncbi_url)
        parent.statusBar().showMessage("Opening NCBI BLAST in web browser...")
    except Exception as e:
        QMessageBox.warning(parent, "Browser Error", f"Could not open web browser: {str(e)}")


def open_blast_file_generic(parent, program_type):
    """Generic function to open BLAST file."""
    
    file_path, _ = QFileDialog.getOpenFileName(
        parent, "Open Sequence File", "", "FASTA Files (*.fasta *.fa *.fas);;All Files (*)"
    )
    
    if file_path:
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            file_label = getattr(parent, f'{program_type}_file_label')
            query_input = getattr(parent, f'{program_type}_query_input')
            
            file_label.setText(os.path.basename(file_path))
            query_input.setText(content)
            
            # Validate the loaded sequence
            sequence_type = "protein" if program_type in ["blastp", "tblastn"] else "nucleotide"
            is_valid, message = validate_sequence(content, sequence_type)
            if not is_valid:
                QMessageBox.warning(
                    parent, 
                    "Invalid Sequence", 
                    f"The loaded file contains invalid sequence data:\n\n{message}"
                )
        except Exception as e:
            QMessageBox.critical(
                parent, 
                "File Error", 
                f"Could not read the selected file:\n\n{str(e)}"
            )


def show_blast_help_generic(parent, program_type):
    """Generic function to show BLAST help."""
    
    help_dialog = QDialog(parent)
    help_dialog.setWindowTitle(f"Online {program_type.upper()} Help")
    help_dialog.setMinimumSize(600, 500)
    
    layout = QVBoxLayout(help_dialog)
    
    help_text = QTextEdit()
    help_text.setReadOnly(True)
    
    sequence_type = "protein" if program_type in ["blastp", "tblastn"] else "nucleotide"
    
    help_text.setHtml(f"""
    <h2>🧬 Online {program_type.upper()} Help</h2>
    
    <h3>📖 Overview</h3>
    <p>This tool uses NCBI's online BLAST service to search {sequence_type} databases. 
    No local BLAST+ installation is required!</p>
    
    <h3>⚙️ Parameters Guide</h3>
    
    <h4>E-value Threshold</h4>
    <ul>
        <li><b>10</b> - Default, finds most matches</li>
        <li><b>0.001</b> - More stringent, fewer false positives</li>
        <li><b>1e-10</b> - Very stringent, only close matches</li>
    </ul>
    
    <h4>Word Size</h4>
    <ul>
        <li>Smaller values are more sensitive but slower</li>
        <li>Larger values are faster but less sensitive</li>
        <li>Default varies by program type</li>
    </ul>
    
    {"<h4>Scoring Matrix</h4><ul><li><b>BLOSUM62</b> - Default, good for most searches</li><li><b>BLOSUM45</b> - More sensitive for distant relationships</li><li><b>BLOSUM80</b> - Less sensitive, for close relationships</li></ul>" if program_type in ["blastp", "blastx", "tblastn", "tblastx"] else ""}
    
    {"<h4>Match/Mismatch Scores</h4><ul><li>First number is reward for match</li><li>Second number is penalty for mismatch</li><li>Higher ratios are more stringent</li></ul>" if program_type in ["blastn", "tblastx"] else ""}
    
    <h4>Gap Costs</h4>
    <ul>
        <li><b>Existence</b> - Penalty for opening a gap</li>
        <li><b>Extension</b> - Penalty for extending a gap</li>
        <li>Higher values discourage gaps</li>
    </ul>
    
    {"<h4>Compositional Adjustments</h4><ul><li><b>Conditional</b> - Default, adjusts for unusual amino acid compositions</li><li><b>No adjustment</b> - For typical protein sequences</li></ul>" if program_type in ["blastp", "blastx", "tblastn"] else ""}
    
    {"<h4>Genetic Code</h4><ul><li>Used for translating nucleotide sequences</li><li><b>Standard (1)</b> - Most common genetic code</li><li>Choose organism-specific codes when appropriate</li></ul>" if program_type in ["blastx", "tblastn", "tblastx"] else ""}
    
    <h3>⏱️ Search Times</h3>
    <p>Online searches typically take 30 seconds to 2 minutes depending on:</p>
    <ul>
        <li>Sequence length</li>
        <li>Database size</li>
        <li>NCBI server load</li>
        <li>Search parameters</li>
    </ul>
    
    <h3>📏 Limitations</h3>
    <ul>
        <li>Maximum sequence length: 20,000 characters</li>
        <li>Requires internet connection</li>
        <li>Subject to NCBI usage policies</li>
    </ul>
    
    <h3>💡 Tips</h3>
    <ul>
        <li>Start with default parameters</li>
        <li>Use organism filter to narrow results</li>
        <li>Save results before starting new searches</li>
        <li>Check sequence format if search fails</li>
    </ul>
    """)
    layout.addWidget(help_text)
    
    button_box = QDialogButtonBox(QDialogButtonBox.Ok)
    button_box.accepted.connect(help_dialog.accept)
    layout.addWidget(button_box)
    
    help_dialog.exec_()


def create_placeholder_tab(name, description):
    """Create a placeholder tab for future BLAST types."""
    widget = QWidget()
    layout = QVBoxLayout(widget)
    
    placeholder = QLabel(f"🚧 {name} functionality coming soon!\n\n"
                        f"This will include:\n"
                        f"• {description}\n"
                        f"• Online NCBI integration\n"
                        f"• Parameter customization\n"
                        f"• Results visualization")
    placeholder.setAlignment(Qt.AlignCenter)
    placeholder.setStyleSheet("font-size: 14px; color: #666; margin: 50px;")
    layout.addWidget(placeholder)
    
    return widget


def save_blast_results(parent):
    """Save BLAST results to a file."""
    
    if not hasattr(parent, 'blastp_results_display') or not parent.blastp_results_display.toPlainText().strip():
        QMessageBox.warning(parent, "No Results", "No BLAST results to save.")
        return
    
    # Generate default filename with timestamp
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    default_filename = f"online_blast_results_{timestamp}.txt"
    
    file_path, _ = QFileDialog.getSaveFileName(
        parent, 
        "Save BLAST Results", 
        default_filename, 
        "Text Files (*.txt);;All Files (*)"
    )
    
    if file_path:
        try:
            content = parent.blastp_results_display.toPlainText()
            
            # Add header information
            header = f"# Online BLAST Results saved on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            header += f"# Query: {parent.blastp_query_input.toPlainText()[:100]}...\n"
            header += f"# Database: {parent.blastp_database_combo.currentText()}\n"
            header += "# Generated by PicoMol using NCBI's online BLAST service\n"
            header += "# " + "="*60 + "\n\n"
            
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(header + content)
            
            QMessageBox.information(parent, "Results Saved", f"BLAST results saved to:\n{file_path}")
            parent.statusBar().showMessage(f"Results saved to {os.path.basename(file_path)}")
        except Exception as e:
            QMessageBox.critical(parent, "Save Error", f"Failed to save results:\n{str(e)}")


def open_ncbi_blast(parent):
    """Open NCBI BLAST website with current query."""
    
    if not hasattr(parent, 'blastp_query_input'):
        return
        
    query = parent.blastp_query_input.toPlainText().strip()
    if not query:
        QMessageBox.warning(parent, "No Query", "Please enter a sequence to search.")
        return
    
    # Clean sequence for URL
    clean_query = clean_fasta_sequence(query)
    
    if len(clean_query) > 8000:  # NCBI has URL length limits
        QMessageBox.warning(
            parent, 
            "Sequence Too Long", 
            "The sequence is too long to pass via URL. Please copy and paste it manually into NCBI BLAST."
        )
        # Open NCBI BLAST without query
        ncbi_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch"
    else:
        # Encode the query for URL
        encoded_query = quote(clean_query)
        ncbi_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&QUERY={encoded_query}"
    
    try:
        webbrowser.open_new_tab(ncbi_url)
        parent.statusBar().showMessage("Opening NCBI BLAST in web browser...")
    except Exception as e:
        QMessageBox.warning(parent, "Browser Error", f"Could not open web browser: {str(e)}")


def open_blast_file(parent):
    """Open a FASTA file and load it into the query input."""
    
    file_path, _ = QFileDialog.getOpenFileName(
        parent, "Open Sequence File", "", "FASTA Files (*.fasta *.fa *.fas);;All Files (*)"
    )
    
    if file_path:
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            parent.blastp_file_label.setText(os.path.basename(file_path))
            parent.blastp_query_input.setText(content)
            
            # Validate the loaded sequence
            is_valid, message = validate_sequence(content, "protein")
            if not is_valid:
                QMessageBox.warning(
                    parent, 
                    "Invalid Sequence", 
                    f"The loaded file contains invalid sequence data:\n\n{message}"
                )
        except Exception as e:
            QMessageBox.critical(
                parent, 
                "File Error", 
                f"Could not read the selected file:\n\n{str(e)}"
            )


def show_blast_help(parent):
    """Show BLAST help dialog."""
    
    help_dialog = QDialog(parent)
    help_dialog.setWindowTitle("Online BLAST Help")
    help_dialog.setMinimumSize(600, 500)
    
    layout = QVBoxLayout(help_dialog)
    
    help_text = QTextEdit()
    help_text.setReadOnly(True)
    help_text.setHtml("""
    <h2>🧬 Online BLAST Help</h2>
    
    <h3>📖 Overview</h3>
    <p>This tool uses NCBI's online BLAST service to search protein databases. 
    No local BLAST+ installation is required!</p>
    
    <h3>⚙️ Parameters Guide</h3>
    
    <h4>E-value Threshold</h4>
    <ul>
        <li><b>10</b> - Default, finds most matches</li>
        <li><b>0.001</b> - More stringent, fewer false positives</li>
        <li><b>1e-10</b> - Very stringent, only close matches</li>
    </ul>
    
    <h4>Word Size</h4>
    <ul>
        <li><b>3</b> - Default for proteins, good balance</li>
        <li><b>2</b> - More sensitive, slower</li>
        <li><b>6</b> - Less sensitive, faster</li>
    </ul>
    
    <h4>Scoring Matrix</h4>
    <ul>
        <li><b>BLOSUM62</b> - Default, good for most searches</li>
        <li><b>BLOSUM45</b> - More sensitive for distant relationships</li>
        <li><b>BLOSUM80</b> - Less sensitive, for close relationships</li>
    </ul>
    
    <h4>Gap Costs</h4>
    <ul>
        <li><b>Existence 11, Extension 1</b> - Default, common for protein alignments</li>
        <li>Adjust to control penalties for opening and extending gaps.</li>
    </ul>

    <h4>Compositional Adjustments</h4>
    <ul>
        <li><b>No adjustment</b> - For sequences with typical amino acid composition.</li>
        <li><b>Composition-based statistics</b> - Adjusts scores for sequences with unusual compositions.</li>
        <li><b>Conditional compositional score matrix adjustment</b> - Default, suitable for most protein searches.</li>
        <li><b>Universal compositional score matrix adjustment</b> - More aggressive adjustment.</li>
    </ul>

    <h4>Word Score Threshold</h4>
    <ul>
        <li><b>11</b> - Default, determines how significant a match must be to be considered a "word".</li>
        <li>Lower values increase sensitivity but also search time.</li>
    </ul>

    <h4>Query Genetic Code</h4>
    <ul>
        <li>Specifies the genetic code to use for translating the query sequence.</li>
        <li>Important for accurate searches if the query is a nucleotide sequence that needs translation (though this is blastp, it's a general BLAST option).</li>
    </ul>

    <h4>Ungapped Alignment Only</h4>
    <ul>
        <li>If checked, only ungapped alignments will be reported.</li>
        <li>Useful for very fast searches or specific types of analyses where gaps are not expected.</li>
    </ul>
    
    <h4>Databases</h4>
    <ul>
        <li><b>nr</b> - Non-redundant protein sequences (largest)</li>
        <li><b>refseq_protein</b> - High-quality reference proteins</li>
        <li><b>swissprot</b> - Manually annotated, high quality</li>
        <li><b>pdb</b> - Protein structures from PDB</li>
    </ul>
    
    <h3>⏱️ Search Times</h3>
    <p>Online searches typically take 30 seconds to 2 minutes depending on:</p>
    <ul>
        <li>Sequence length</li>
        <li>Database size</li>
        <li>NCBI server load</li>
        <li>Search parameters</li>
    </ul>
    
    <h3>📏 Limitations</h3>
    <ul>
        <li>Maximum sequence length: 20,000 amino acids</li>
        <li>Requires internet connection</li>
        <li>Subject to NCBI usage policies</li>
    </ul>
    
    <h3>💡 Tips</h3>
    <ul>
        <li>Start with default parameters</li>
        <li>Use organism filter to narrow results</li>
        <li>Save results before starting new searches</li>
        <li>Check sequence format if search fails</li>
    </ul>
    """)
    layout.addWidget(help_text)
    
    button_box = QDialogButtonBox(QDialogButtonBox.Ok)
    button_box.accepted.connect(help_dialog.accept)
    layout.addWidget(button_box)
    
    help_dialog.exec_()


