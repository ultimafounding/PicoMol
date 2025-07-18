"""
Fixed BLASTN utilities for PicoMol - Exact NCBI BLASTN layout matching blastn.html.

This module recreates the exact BLASTN layout and functionality from NCBI's blastn.html,
providing an authentic user experience that matches the official NCBI BLASTN interface.
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
    QScrollArea, QSplitter, QTabWidget, QProgressBar, QDialogButtonBox, QDialog,
    QFrame, QGridLayout, QSpacerItem, QSizePolicy
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont, QPalette

# Import the worker class from the original blast_utils
from blast_utils import OnlineBlastWorker, validate_sequence, clean_fasta_sequence, format_blast_output


def create_ncbi_style_blastn_tab_fixed(parent):
    """Create BLASTN tab with exact NCBI layout from blastn.html."""
    
    # Main widget
    main_widget = QWidget()
    main_layout = QVBoxLayout(main_widget)
    main_layout.setContentsMargins(10, 10, 10, 10)
    main_layout.setSpacing(0)
    
    # Initialize worker
    parent.blastn_worker = None
    
    # Page title
    title_label = QLabel("Standard Nucleotide BLAST")
    title_label.setStyleSheet("""
        QLabel {
            font-size: 18px;
            font-weight: bold;
            color: #333;
            padding: 10px 0;
            border-bottom: 1px solid #ddd;
            margin-bottom: 15px;
        }
    """)
    main_layout.addWidget(title_label)
    
    # Program description
    desc_label = QLabel("BLASTN searches nucleotide databases using a nucleotide query.")
    desc_label.setStyleSheet("""
        QLabel {
            font-size: 14px;
            color: #666;
            padding: 5px 0;
            margin-bottom: 10px;
        }
    """)
    main_layout.addWidget(desc_label)
    
    # Old description removed - using consistent format above
    
    # Create scroll area for the form
    scroll_area = QScrollArea()
    scroll_area.setWidgetResizable(True)
    scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
    scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
    
    # Form widget
    form_widget = QWidget()
    form_layout = QVBoxLayout(form_widget)
    form_layout.setContentsMargins(0, 0, 0, 0)
    form_layout.setSpacing(20)
    
    # 1. Enter Query Sequence Section
    query_section = create_blastn_query_section(parent)
    form_layout.addWidget(query_section)
    
    # 2. Choose Search Set Section  
    search_set_section = create_blastn_search_set_section(parent)
    form_layout.addWidget(search_set_section)
    
    # 3. Program Selection Section
    program_section = create_blastn_program_selection_section(parent)
    form_layout.addWidget(program_section)
    
    # 4. Algorithm Parameters Section (collapsible)
    params_section = create_blastn_algorithm_parameters_section(parent)
    form_layout.addWidget(params_section)
    
    # Add stretch
    form_layout.addStretch()
    
    # Set form widget to scroll area
    scroll_area.setWidget(form_widget)
    main_layout.addWidget(scroll_area)
    
    # Results section
    results_section = create_blastn_results_section(parent)
    main_layout.addWidget(results_section)
    
    return main_widget


def create_blastn_query_section(parent):
    """Create the 'Enter Query Sequence' section exactly like NCBI BLASTN."""
    
    section = QFrame()
    section.setFrameStyle(QFrame.Box)
    section.setStyleSheet("""
        QFrame {
            border: 1px solid #ccc;
            border-radius: 4px;
            background-color: #fafafa;
        }
    """)
    
    layout = QVBoxLayout(section)
    layout.setContentsMargins(15, 10, 15, 15)
    
    # Legend/Title
    legend = QLabel("Enter Query Sequence")
    legend.setStyleSheet("""
        QLabel {
            font-weight: bold;
            font-size: 14px;
            color: #333;
            background-color: #fafafa;
            padding: 0 5px;
            margin-bottom: 10px;
        }
    """)
    layout.addWidget(legend)
    
    # Query input area
    query_layout = QVBoxLayout()
    
    # Label with help link
    label_layout = QHBoxLayout()
    query_label = QLabel("Enter accession number(s), gi(s), or FASTA sequence(s)")
    query_label.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    label_layout.addWidget(query_label)
    
    help_button = QPushButton("?")
    help_button.setFixedSize(20, 20)
    help_button.setStyleSheet("""
        QPushButton {
            background-color: #007cba;
            color: white;
            border: none;
            border-radius: 10px;
            font-size: 12px;
            font-weight: bold;
        }
        QPushButton:hover {
            background-color: #005a87;
        }
    """)
    help_button.setToolTip("How to enter queries")
    label_layout.addWidget(help_button)
    
    clear_button = QPushButton("Clear")
    clear_button.setStyleSheet("""
        QPushButton {
            background: none;
            border: none;
            color: #007cba;
            text-decoration: underline;
            font-size: 12px;
        }
        QPushButton:hover {
            color: #005a87;
        }
    """)
    clear_button.clicked.connect(lambda: parent.blastn_query_input.clear())
    label_layout.addWidget(clear_button)
    label_layout.addStretch()
    
    query_layout.addLayout(label_layout)
    
    # Text area
    query_input = QTextEdit()
    query_input.setFixedHeight(120)
    query_input.setStyleSheet("""
        QTextEdit {
            border: 1px solid #ccc;
            font-family: monospace;
            font-size: 12px;
            background-color: white;
        }
    """)
    query_input.setPlaceholderText("Enter nucleotide sequence here...\n\nExample:\n>My_sequence\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
    parent.blastn_query_input = query_input
    query_layout.addWidget(query_input)
    
    layout.addLayout(query_layout)
    
    # Query subrange
    subrange_layout = QHBoxLayout()
    subrange_label = QLabel("Query subrange")
    subrange_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
    subrange_layout.addWidget(subrange_label)
    
    from_label = QLabel("From")
    from_input = QLineEdit()
    from_input.setFixedWidth(80)
    from_input.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    parent.blastn_query_from = from_input
    
    to_label = QLabel("To")
    to_input = QLineEdit()
    to_input.setFixedWidth(80)
    to_input.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    parent.blastn_query_to = to_input
    
    subrange_layout.addWidget(from_label)
    subrange_layout.addWidget(from_input)
    subrange_layout.addWidget(to_label)
    subrange_layout.addWidget(to_input)
    subrange_layout.addStretch()
    
    layout.addLayout(subrange_layout)
    
    # File upload
    upload_layout = QHBoxLayout()
    upload_label = QLabel("Or, upload file")
    upload_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
    upload_layout.addWidget(upload_label)
    
    upload_button = QPushButton("Choose File")
    upload_button.setStyleSheet("""
        QPushButton {
            background-color: #f8f9fa;
            border: 1px solid #ccc;
            padding: 4px 8px;
            border-radius: 3px;
        }
        QPushButton:hover {
            background-color: #e9ecef;
        }
    """)
    upload_button.clicked.connect(lambda: open_blastn_file(parent))
    upload_layout.addWidget(upload_button)
    
    file_label = QLabel("No file selected")
    file_label.setStyleSheet("color: #666; margin-left: 10px;")
    parent.blastn_file_label = file_label
    upload_layout.addWidget(file_label)
    upload_layout.addStretch()
    
    layout.addLayout(upload_layout)
    
    # Job Title
    job_layout = QHBoxLayout()
    job_label = QLabel("Job Title")
    job_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
    job_layout.addWidget(job_label)
    
    job_input = QLineEdit()
    job_input.setFixedWidth(400)
    job_input.setStyleSheet("border: 1px solid #ccc; padding: 4px;")
    job_input.setPlaceholderText("Enter a descriptive title for your BLAST search")
    parent.blastn_job_title = job_input
    job_layout.addWidget(job_input)
    job_layout.addStretch()
    
    layout.addLayout(job_layout)
    
    return section


def create_blastn_search_set_section(parent):
    """Create the 'Choose Search Set' section exactly like NCBI BLASTN."""
    
    section = QFrame()
    section.setFrameStyle(QFrame.Box)
    section.setStyleSheet("""
        QFrame {
            border: 1px solid #ccc;
            border-radius: 4px;
            background-color: #fafafa;
        }
    """)
    
    layout = QVBoxLayout(section)
    layout.setContentsMargins(15, 10, 15, 15)
    
    # Legend
    legend = QLabel("Choose Search Set")
    legend.setStyleSheet("""
        QLabel {
            font-weight: bold;
            font-size: 14px;
            color: #333;
            background-color: #fafafa;
            padding: 0 5px;
            margin-bottom: 10px;
        }
    """)
    layout.addWidget(legend)
    
    # Database type selection (radio buttons like NCBI)
    db_type_layout = QVBoxLayout()
    
    # Create button group for database type selection
    db_type_group = QButtonGroup()
    
    # Standard databases radio button
    standard_radio = QRadioButton("Standard databases (nr etc.):")
    standard_radio.setChecked(True)
    standard_radio.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    parent.blastn_standard_radio = standard_radio
    db_type_group.addButton(standard_radio, 0)
    db_type_layout.addWidget(standard_radio)
    
    # rRNA/ITS databases radio button
    rrna_radio = QRadioButton("rRNA/ITS databases")
    rrna_radio.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    parent.blastn_rrna_radio = rrna_radio
    db_type_group.addButton(rrna_radio, 1)
    db_type_layout.addWidget(rrna_radio)
    
    # Genomic + transcript databases radio button
    genomic_radio = QRadioButton("Genomic + transcript databases")
    genomic_radio.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    parent.blastn_genomic_radio = genomic_radio
    db_type_group.addButton(genomic_radio, 2)
    db_type_layout.addWidget(genomic_radio)
    
    # Betacoronavirus radio button
    virus_radio = QRadioButton("Betacoronavirus")
    virus_radio.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    parent.blastn_virus_radio = virus_radio
    db_type_group.addButton(virus_radio, 3)
    db_type_layout.addWidget(virus_radio)
    
    parent.blastn_db_type_group = db_type_group
    layout.addLayout(db_type_layout)
    
    # Database dropdown
    db_layout = QHBoxLayout()
    db_label = QLabel("Database")
    db_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
    db_layout.addWidget(db_label)
    
    database_combo = QComboBox()
    database_combo.setFixedWidth(500)
    database_combo.setStyleSheet("""
        QComboBox {
            border: 1px solid #ccc;
            padding: 4px;
            background-color: white;
        }
        QComboBox::drop-down {
            border: none;
        }
        QComboBox::down-arrow {
            width: 12px;
            height: 12px;
        }
    """)
    
    # Define all database categories from blastn.html
    all_databases = {
        'standard': [
            ("core_nt", "Core nucleotide database (core_nt)"),
            ("refseq_select_rna", "RefSeq Select RNA sequences (refseq_select)"),
            ("refseq_rna", "Reference RNA sequences (refseq_rna)"),
            ("refseq_representative_genomes", "RefSeq Reference genomes (refseq_reference_genomes)"),
            ("refseq_genomes", "RefSeq Genome Database (refseq_genomes)"),
            ("nt", "Nucleotide collection (nr/nt)"),
            ("Whole_Genome_Shotgun_contigs", "Whole-genome shotgun contigs (wgs)"),
            ("est", "Expressed sequence tags (est)"),
            ("sra", "Sequence Read Archive (SRA)"),
            ("tsa_nt", "Transcriptome Shotgun Assembly (TSA)"),
            ("tls_nt", "Targeted Loci(TLS)"),
            ("htgs", "High throughput genomic sequences (HTGS)"),
            ("patnt", "Patent sequences(pat)"),
            ("pdbnt", "PDB nucleotide database (pdb)"),
            ("genomic/9606/RefSeqGene", "Human RefSeqGene sequences(RefSeq_Gene)"),
            ("gss", "Genomic survey sequences (gss)"),
            ("sts", "Sequence tagged sites (dbsts)")
        ],
        'rrna': [
            ("rRNA_typestrains/16S_ribosomal_RNA", "16S ribosomal RNA sequences (Bacteria and Archaea)"),
            ("TL/18S_fungal_sequences", "18S ribosomal RNA sequences (SSU) from Fungi type and reference material"),
            ("TL/28S_fungal_sequences", "28S ribosomal RNA sequences (LSU) from Fungi type and reference material"),
            ("rRNA_typestrains/ITS_RefSeq_Fungi", "Internal transcribed spacer region (ITS) from Fungi type and reference material")
        ],
        'genomic': [
            ("GPIPE/9606/current/ref_top_level GPIPE/9606/current/rna", "Human genomic plus transcript (Human G+T)"),
            ("GPIPE/10090/current/ref_top_level GPIPE/10090/current/rna", "Mouse genomic plus transcript (Mouse G+T)")
        ],
        'virus': [
            ("genomic/Viruses/Betacoronavirus", "Betacoronavirus Genbank")
        ]
    }
    
    # Function to update database list based on radio selection
    def update_blastn_database_list():
        """Update database dropdown based on selected database type."""
        database_combo.clear()
        
        if parent.blastn_standard_radio.isChecked():
            databases = all_databases['standard']
        elif parent.blastn_rrna_radio.isChecked():
            databases = all_databases['rrna']
        elif parent.blastn_genomic_radio.isChecked():
            databases = all_databases['genomic']
        elif parent.blastn_virus_radio.isChecked():
            databases = all_databases['virus']
        else:
            databases = all_databases['standard']
        
        for db_id, db_name in databases:
            database_combo.addItem(db_name, db_id)
    
    # Initial population with standard databases
    for db_id, db_name in all_databases['standard']:
        database_combo.addItem(db_name, db_id)
    
    # Connect radio buttons to update function
    for button in db_type_group.buttons():
        button.toggled.connect(update_blastn_database_list)
    
    parent.blastn_update_database_list = update_blastn_database_list
    
    parent.blastn_database_combo = database_combo
    db_layout.addWidget(database_combo)
    
    help_db = QPushButton("?")
    help_db.setFixedSize(20, 20)
    help_db.setStyleSheet("""
        QPushButton {
            background-color: #007cba;
            color: white;
            border: none;
            border-radius: 10px;
            font-size: 12px;
            font-weight: bold;
        }
        QPushButton:hover {
            background-color: #005a87;
        }
    """)
    help_db.setToolTip("Database help")
    db_layout.addWidget(help_db)
    db_layout.addStretch()
    
    layout.addLayout(db_layout)
    
    # Organism
    org_layout = QHBoxLayout()
    org_label = QLabel("Organism\nOptional")
    org_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
    org_layout.addWidget(org_label)
    
    organism_input = QLineEdit()
    organism_input.setFixedWidth(300)
    organism_input.setStyleSheet("border: 1px solid #ccc; padding: 4px;")
    organism_input.setPlaceholderText("Enter organism common name, binomial, or tax id")
    parent.blastn_organism_input = organism_input
    org_layout.addWidget(organism_input)
    
    exclude_checkbox = QCheckBox("exclude")
    exclude_checkbox.setStyleSheet("margin-left: 10px;")
    parent.blastn_exclude_organism = exclude_checkbox
    org_layout.addWidget(exclude_checkbox)
    org_layout.addStretch()
    
    layout.addLayout(org_layout)
    
    # Exclude options
    exclude_layout = QVBoxLayout()
    exclude_label = QLabel("Exclude\nOptional")
    exclude_label.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    exclude_layout.addWidget(exclude_label)
    
    exclude_models = QCheckBox("Models (XM/XP)")
    exclude_models.setStyleSheet("margin-left: 20px;")
    parent.blastn_exclude_models = exclude_models
    exclude_layout.addWidget(exclude_models)
    
    exclude_uncult = QCheckBox("Uncultured/environmental sample sequences")
    exclude_uncult.setStyleSheet("margin-left: 20px;")
    parent.blastn_exclude_uncult = exclude_uncult
    exclude_layout.addWidget(exclude_uncult)
    
    layout.addLayout(exclude_layout)
    
    return section


def create_blastn_program_selection_section(parent):
    """Create the 'Program Selection' section exactly like NCBI BLASTN."""
    
    section = QFrame()
    section.setFrameStyle(QFrame.Box)
    section.setStyleSheet("""
        QFrame {
            border: 1px solid #ccc;
            border-radius: 4px;
            background-color: #fafafa;
        }
    """)
    
    layout = QVBoxLayout(section)
    layout.setContentsMargins(15, 10, 15, 15)
    
    # Legend
    legend = QLabel("Program Selection")
    legend.setStyleSheet("""
        QLabel {
            font-weight: bold;
            font-size: 14px;
            color: #333;
            background-color: #fafafa;
            padding: 0 5px;
            margin-bottom: 10px;
        }
    """)
    layout.addWidget(legend)
    
    # Optimize for label
    optimize_label = QLabel("Optimize for")
    optimize_label.setStyleSheet("font-weight: normal; margin-bottom: 10px;")
    layout.addWidget(optimize_label)
    
    # Radio buttons for algorithms (exactly from blastn.html)
    algo_group = QButtonGroup()
    algo_layout = QVBoxLayout()
    
    algorithms = [
        ("megaBlast", "Highly similar sequences (megablast)"),
        ("discoMegablast", "More dissimilar sequences (discontiguous megablast)"),
        ("blastn", "Somewhat similar sequences (blastn)")
    ]
    
    for i, (algo_id, algo_desc) in enumerate(algorithms):
        radio = QRadioButton(algo_desc)
        radio.setStyleSheet("margin-bottom: 5px;")
        if algo_id == "megaBlast":
            radio.setChecked(True)
        algo_group.addButton(radio, i)
        setattr(parent, f'blastn_{algo_id}_radio', radio)
        algo_layout.addWidget(radio)
    
    parent.blastn_algorithm_group = algo_group
    layout.addLayout(algo_layout)
    
    # Help text
    help_text = QLabel("Choose a BLAST algorithm")
    help_text.setStyleSheet("color: #666; font-size: 12px; margin-top: 10px;")
    layout.addWidget(help_text)
    
    # Connect algorithm radio buttons to update summary
    def update_blastn_algorithm_ui():
        """Update UI elements based on selected algorithm."""
        checked_button = algo_group.checkedButton()
        if checked_button:
            button_text = checked_button.text()
            
            # Update summary text
            summary_label = parent.blastn_summary_label
            if "megablast" in button_text.lower():
                if "discontiguous" in button_text.lower():
                    summary_text = "Search nt using Discontiguous Megablast (More dissimilar sequences)"
                else:
                    summary_text = "Search nt using Megablast (Optimize for highly similar sequences)"
            elif "blastn" in button_text.lower():
                summary_text = "Search nt using Blastn (Somewhat similar sequences)"
            else:
                summary_text = "Search nt using Blastn (nucleotide-nucleotide BLAST)"
            
            summary_label.setText(summary_text)
    
    # Connect all radio buttons to the update function
    for button in algo_group.buttons():
        button.toggled.connect(update_blastn_algorithm_ui)
    
    # Store the update function for external access
    parent.blastn_update_algorithm_ui = update_blastn_algorithm_ui
    
    # Search button and summary
    search_layout = QVBoxLayout()
    
    blast_button = QPushButton("BLAST")
    blast_button.setFixedSize(100, 40)
    blast_button.setStyleSheet("""
        QPushButton {
            background-color: #007cba;
            color: white;
            border: none;
            border-radius: 4px;
            font-size: 16px;
            font-weight: bold;
        }
        QPushButton:hover {
            background-color: #005a87;
        }
        QPushButton:pressed {
            background-color: #004a73;
        }
    """)
    blast_button.clicked.connect(lambda: run_ncbi_blastn_search(parent))
    parent.blastn_blast_button = blast_button
    
    # Cancel button (initially hidden)
    cancel_button = QPushButton("Cancel")
    cancel_button.setFixedSize(100, 40)
    cancel_button.setStyleSheet("""
        QPushButton {
            background-color: #dc3545;
            color: white;
            border: none;
            border-radius: 4px;
            font-size: 14px;
            font-weight: bold;
        }
        QPushButton:hover {
            background-color: #c82333;
        }
    """)
    cancel_button.setVisible(False)
    cancel_button.clicked.connect(lambda: cancel_ncbi_blastn_search(parent))
    parent.blastn_cancel_button = cancel_button
    
    button_layout = QHBoxLayout()
    button_layout.addWidget(blast_button)
    button_layout.addWidget(cancel_button)
    button_layout.addStretch()
    
    search_layout.addLayout(button_layout)
    
    # Progress bar
    progress_bar = QProgressBar()
    progress_bar.setVisible(False)
    progress_bar.setRange(0, 0)  # Indeterminate
    parent.blastn_progress_bar = progress_bar
    search_layout.addWidget(progress_bar)
    
    # Search summary
    summary_layout = QVBoxLayout()
    summary_label = QLabel("Search nt using Megablast (Optimize for highly similar sequences)")
    summary_label.setStyleSheet("color: #666; font-size: 12px; margin-top: 10px;")
    parent.blastn_summary_label = summary_label
    summary_layout.addWidget(summary_label)
    
    # Removed "Show results in a new window" checkbox
    
    search_layout.addLayout(summary_layout)
    layout.addLayout(search_layout)
    
    return section


def create_blastn_algorithm_parameters_section(parent):
    """Create the collapsible 'Algorithm parameters' section exactly like NCBI BLASTN."""
    
    section = QFrame()
    section.setFrameStyle(QFrame.Box)
    section.setStyleSheet("""
        QFrame {
            border: 1px solid #ccc;
            border-radius: 4px;
            background-color: #fafafa;
        }
    """)
    
    layout = QVBoxLayout(section)
    layout.setContentsMargins(15, 10, 15, 15)
    
    # Collapsible header
    header_layout = QHBoxLayout()
    
    toggle_button = QPushButton("▼ Algorithm parameters")
    toggle_button.setStyleSheet("""
        QPushButton {
            background: none;
            border: none;
            font-weight: bold;
            font-size: 14px;
            color: #333;
            text-align: left;
            padding: 5px 0;
        }
        QPushButton:hover {
            color: #007cba;
        }
    """)
    
    # Parameters content (initially hidden)
    params_content = QWidget()
    params_content.setVisible(False)
    params_layout = QVBoxLayout(params_content)
    
    def toggle_params():
        """Toggle parameters visibility."""
        visible = params_content.isVisible()
        params_content.setVisible(not visible)
        toggle_button.setText("▲ Algorithm parameters" if not visible else "▼ Algorithm parameters")
    
    toggle_button.clicked.connect(toggle_params)
    parent.blastn_toggle_params = toggle_button
    parent.blastn_params_content = params_content
    
    header_layout.addWidget(toggle_button)
    
    reset_button = QPushButton("Restore default search parameters")
    reset_button.setStyleSheet("""
        QPushButton {
            background-color: #f8f9fa;
            border: 1px solid #ccc;
            padding: 4px 8px;
            border-radius: 3px;
            font-size: 12px;
        }
        QPushButton:hover {
            background-color: #e9ecef;
        }
    """)
    header_layout.addWidget(reset_button)
    header_layout.addStretch()
    
    layout.addLayout(header_layout)
    
    # Note about non-default parameters
    note_label = QLabel("Note: Parameter values that differ from the default are highlighted in yellow and marked with ♦ sign")
    note_label.setStyleSheet("color: #666; font-size: 11px; margin: 5px 0;")
    layout.addWidget(note_label)
    
    # General Parameters
    general_group = create_blastn_general_parameters_group(parent)
    params_layout.addWidget(general_group)
    
    # Scoring Parameters
    scoring_group = create_blastn_scoring_parameters_group(parent)
    params_layout.addWidget(scoring_group)
    
    # Filters and Masking
    filters_group = create_blastn_filters_masking_group(parent)
    params_layout.addWidget(filters_group)
    
    layout.addWidget(params_content)
    
    return section


def create_blastn_general_parameters_group(parent):
    """Create General Parameters group for BLASTN."""
    
    group = QFrame()
    group.setFrameStyle(QFrame.StyledPanel)
    group.setStyleSheet("""
        QFrame {
            border: 1px solid #ddd;
            border-radius: 4px;
            background-color: white;
            margin: 5px 0;
        }
    """)
    
    layout = QVBoxLayout(group)
    layout.setContentsMargins(10, 10, 10, 10)
    
    # Title
    title = QLabel("General Parameters")
    title.setStyleSheet("font-weight: bold; font-size: 13px; margin-bottom: 10px;")
    layout.addWidget(title)
    
    # Max target sequences
    max_layout = QHBoxLayout()
    max_label = QLabel("Max target sequences")
    max_label.setFixedWidth(150)
    max_layout.addWidget(max_label)
    
    max_combo = QComboBox()
    max_combo.setFixedWidth(100)
    max_combo.addItems(["10", "50", "100", "250", "500", "1000", "5000"])
    max_combo.setCurrentText("100")
    max_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    parent.blastn_max_targets = max_combo
    max_layout.addWidget(max_combo)
    max_layout.addStretch()
    
    layout.addLayout(max_layout)
    
    # Short queries
    short_layout = QHBoxLayout()
    short_checkbox = QCheckBox("Automatically adjust parameters for short input sequences")
    short_checkbox.setChecked(True)
    short_checkbox.setStyleSheet("margin: 5px 0;")
    parent.blastn_short_adjust = short_checkbox
    short_layout.addWidget(short_checkbox)
    short_layout.addStretch()
    
    layout.addLayout(short_layout)
    
    # Expect threshold
    expect_layout = QHBoxLayout()
    expect_label = QLabel("Expect threshold")
    expect_label.setFixedWidth(150)
    expect_layout.addWidget(expect_label)
    
    expect_input = QLineEdit("0.05")
    expect_input.setFixedWidth(100)
    expect_input.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    parent.blastn_expect_threshold = expect_input
    expect_layout.addWidget(expect_input)
    expect_layout.addStretch()
    
    layout.addLayout(expect_layout)
    
    # Word size (BLASTN specific values from blastn.html)
    word_layout = QHBoxLayout()
    word_label = QLabel("Word size")
    word_label.setFixedWidth(150)
    word_layout.addWidget(word_label)
    
    word_combo = QComboBox()
    word_combo.setFixedWidth(100)
    # From blastn.html: 16, 20, 24, 28 (default), 32, 48, 64
    word_combo.addItems(["16", "20", "24", "28", "32", "48", "64"])
    word_combo.setCurrentText("28")
    word_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    parent.blastn_word_size = word_combo
    word_layout.addWidget(word_combo)
    word_layout.addStretch()
    
    layout.addLayout(word_layout)
    
    # Max matches in query range
    range_layout = QHBoxLayout()
    range_label = QLabel("Max matches in a query range")
    range_label.setFixedWidth(150)
    range_layout.addWidget(range_label)
    
    range_input = QLineEdit("0")
    range_input.setFixedWidth(100)
    range_input.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    parent.blastn_hsp_range_max = range_input
    range_layout.addWidget(range_input)
    range_layout.addStretch()
    
    layout.addLayout(range_layout)
    
    return group


def create_blastn_scoring_parameters_group(parent):
    """Create Scoring Parameters group for BLASTN."""
    
    group = QFrame()
    group.setFrameStyle(QFrame.StyledPanel)
    group.setStyleSheet("""
        QFrame {
            border: 1px solid #ddd;
            border-radius: 4px;
            background-color: white;
            margin: 5px 0;
        }
    """)
    
    layout = QVBoxLayout(group)
    layout.setContentsMargins(10, 10, 10, 10)
    
    # Title
    title = QLabel("Scoring Parameters")
    title.setStyleSheet("font-weight: bold; font-size: 13px; margin-bottom: 10px;")
    layout.addWidget(title)
    
    # Match/Mismatch scores (from blastn.html)
    match_layout = QHBoxLayout()
    match_label = QLabel("Match/Mismatch Scores")
    match_label.setFixedWidth(150)
    match_layout.addWidget(match_label)
    
    match_combo = QComboBox()
    match_combo.setFixedWidth(120)
    # From blastn.html: 1,-2 (default), 1,-3, 1,-4, 2,-3, 4,-5, 1,-1
    match_scores = ["1,-2", "1,-3", "1,-4", "2,-3", "4,-5", "1,-1"]
    match_combo.addItems(match_scores)
    match_combo.setCurrentText("1,-2")
    match_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    parent.blastn_match_scores = match_combo
    match_layout.addWidget(match_combo)
    match_layout.addStretch()
    
    layout.addLayout(match_layout)
    
    # Gap costs (from blastn.html)
    gap_layout = QHBoxLayout()
    gap_label = QLabel("Gap Costs")
    gap_label.setFixedWidth(150)
    gap_layout.addWidget(gap_label)
    
    gap_combo = QComboBox()
    gap_combo.setFixedWidth(200)
    # From blastn.html: Linear, Existence: 5 Extension: 2 (default), etc.
    gap_costs = [
        "Linear",
        "Existence: 5 Extension: 2",
        "Existence: 2 Extension: 2",
        "Existence: 1 Extension: 2",
        "Existence: 0 Extension: 2",
        "Existence: 3 Extension: 1",
        "Existence: 2 Extension: 1",
        "Existence: 1 Extension: 1"
    ]
    gap_combo.addItems(gap_costs)
    gap_combo.setCurrentText("Existence: 5 Extension: 2")
    gap_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    parent.blastn_gap_costs = gap_combo
    gap_layout.addWidget(gap_combo)
    gap_layout.addStretch()
    
    layout.addLayout(gap_layout)
    
    return group


def create_blastn_filters_masking_group(parent):
    """Create Filters and Masking group for BLASTN."""
    
    group = QFrame()
    group.setFrameStyle(QFrame.StyledPanel)
    group.setStyleSheet("""
        QFrame {
            border: 1px solid #ddd;
            border-radius: 4px;
            background-color: white;
            margin: 5px 0;
        }
    """)
    
    layout = QVBoxLayout(group)
    layout.setContentsMargins(10, 10, 10, 10)
    
    # Title
    title = QLabel("Filters and Masking")
    title.setStyleSheet("font-weight: bold; font-size: 13px; margin-bottom: 10px;")
    layout.addWidget(title)
    
    # Filter label
    filter_label = QLabel("Filter")
    filter_label.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    layout.addWidget(filter_label)
    
    # Low complexity regions
    low_comp_layout = QHBoxLayout()
    low_comp_checkbox = QCheckBox("Low complexity regions")
    low_comp_checkbox.setStyleSheet("margin-left: 20px;")
    parent.blastn_low_complexity = low_comp_checkbox
    low_comp_layout.addWidget(low_comp_checkbox)
    low_comp_layout.addStretch()
    
    layout.addLayout(low_comp_layout)
    
    # Species-specific repeats
    repeats_layout = QHBoxLayout()
    repeats_checkbox = QCheckBox("Species-specific repeats filter for:")
    repeats_checkbox.setStyleSheet("margin-left: 20px;")
    parent.blastn_repeats_filter = repeats_checkbox
    repeats_layout.addWidget(repeats_checkbox)
    
    # Species selection (simplified)
    species_combo = QComboBox()
    species_combo.setFixedWidth(200)
    species_combo.addItems(["Human", "Mouse", "Rat", "Drosophila", "C. elegans"])
    species_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    parent.blastn_species_repeats = species_combo
    repeats_layout.addWidget(species_combo)
    repeats_layout.addStretch()
    
    layout.addLayout(repeats_layout)
    
    return group


def create_blastn_results_section(parent):
    """Create the results display section for BLASTN."""
    
    section = QFrame()
    section.setFrameStyle(QFrame.Box)
    section.setStyleSheet("""
        QFrame {
            border: 1px solid #ccc;
            border-radius: 4px;
            background-color: white;
        }
    """)
    
    layout = QVBoxLayout(section)
    layout.setContentsMargins(10, 10, 10, 10)
    
    # Title
    title = QLabel("BLASTN Results")
    title.setStyleSheet("font-weight: bold; font-size: 14px; margin-bottom: 10px;")
    layout.addWidget(title)
    
    # Toolbar
    toolbar_layout = QHBoxLayout()
    
    save_button = QPushButton("💾 Save Results")
    save_button.setStyleSheet("""
        QPushButton {
            background-color: #f8f9fa;
            border: 1px solid #ccc;
            padding: 6px 12px;
            border-radius: 3px;
            font-size: 12px;
        }
        QPushButton:hover {
            background-color: #e9ecef;
        }
    """)
    save_button.clicked.connect(lambda: save_blastn_results(parent))
    toolbar_layout.addWidget(save_button)
    
    clear_button = QPushButton("🗑️ Clear")
    clear_button.setStyleSheet("""
        QPushButton {
            background-color: #f8f9fa;
            border: 1px solid #ccc;
            padding: 6px 12px;
            border-radius: 3px;
            font-size: 12px;
        }
        QPushButton:hover {
            background-color: #e9ecef;
        }
    """)
    clear_button.clicked.connect(lambda: parent.blastn_results_display.clear())
    toolbar_layout.addWidget(clear_button)
    
    ncbi_button = QPushButton("🌐 View on NCBI")
    ncbi_button.setStyleSheet("""
        QPushButton {
            background-color: #f8f9fa;
            border: 1px solid #ccc;
            padding: 6px 12px;
            border-radius: 3px;
            font-size: 12px;
        }
        QPushButton:hover {
            background-color: #e9ecef;
        }
    """)
    ncbi_button.clicked.connect(lambda: open_ncbi_blastn(parent))
    toolbar_layout.addWidget(ncbi_button)
    
    toolbar_layout.addStretch()
    layout.addLayout(toolbar_layout)
    
    # Results display
    results_display = QTextEdit()
    results_display.setReadOnly(True)
    results_display.setFont(QFont("Courier", 10))
    results_display.setStyleSheet("""
        QTextEdit {
            border: 1px solid #ccc;
            background-color: #fafafa;
            font-family: monospace;
        }
    """)
    results_display.setPlaceholderText(
        "BLASTN results will appear here...\n\n"
        "• Enter a nucleotide sequence above\n"
        "• Configure your search parameters\n"
        "• Click the BLAST button to begin\n"
        "• Results typically take 30 seconds to 2 minutes"
    )
    parent.blastn_results_display = results_display
    layout.addWidget(results_display)
    
    return section


# Helper functions for BLASTN
def run_ncbi_blastn_search(parent):
    """Run BLASTN search with NCBI-style interface."""
    
    query = parent.blastn_query_input.toPlainText().strip()
    
    if not query:
        QMessageBox.warning(parent, "No Query", "Please enter a sequence to search.")
        return
    
    # Validate sequence
    is_valid, message = validate_sequence(query, "nucleotide")
    if not is_valid:
        QMessageBox.warning(parent, "Invalid Sequence", message)
        return
    
    # Get parameters
    database_combo = parent.blastn_database_combo
    database = database_combo.currentData()
    
    expect_input = parent.blastn_expect_threshold
    max_targets = parent.blastn_max_targets
    word_size = parent.blastn_word_size
    
    parameters = {
        'evalue': expect_input.text().strip() or "0.05",
        'max_target_seqs': max_targets.currentText(),
        'word_size': word_size.currentText(),
    }
    
    # Add BLASTN-specific parameters
    match_scores = parent.blastn_match_scores
    gap_costs = parent.blastn_gap_costs
    
    parameters['match_scores'] = match_scores.currentText()
    
    # Parse gap costs
    gap_text = gap_costs.currentText()
    if gap_text == "Linear":
        parameters['gapcosts'] = "0 0"
    elif "Existence:" in gap_text:
        parts = gap_text.replace("Existence:", "").replace("Extension:", "").split()
        if len(parts) >= 2:
            parameters['gapopen'] = parts[0]
            parameters['gapextend'] = parts[1]
    
    # Get algorithm
    algorithm_group = parent.blastn_algorithm_group
    checked_button = algorithm_group.checkedButton()
    if checked_button:
        button_text = checked_button.text()
        if "megablast" in button_text.lower():
            if "discontiguous" in button_text.lower():
                parameters['algorithm'] = 'discoMegablast'
            else:
                parameters['algorithm'] = 'megaBlast'
        else:
            parameters['algorithm'] = 'blastn'
    
    # Get organism filter
    organism_input = parent.blastn_organism_input
    if organism_input.text().strip():
        parameters['organism'] = organism_input.text().strip()
    
    # Update UI
    blast_button = parent.blastn_blast_button
    cancel_button = parent.blastn_cancel_button
    progress_bar = parent.blastn_progress_bar
    results_display = parent.blastn_results_display
    
    blast_button.setVisible(False)
    cancel_button.setVisible(True)
    progress_bar.setVisible(True)
    results_display.setText("🚀 Submitting BLASTN search to NCBI...")
    
    # Clean sequence
    clean_sequence = clean_fasta_sequence(query)
    
    # Create and start worker
    worker = OnlineBlastWorker(clean_sequence, database, "blastn", parameters)
    parent.blastn_worker = worker
    
    worker.finished.connect(lambda result: on_ncbi_blastn_finished(parent, result))
    worker.error.connect(lambda error: on_ncbi_blastn_error(parent, error))
    worker.progress.connect(lambda msg: on_ncbi_blastn_progress(parent, msg))
    worker.start()


def cancel_ncbi_blastn_search(parent):
    """Cancel BLASTN search."""
    
    worker = getattr(parent, 'blastn_worker', None)
    if worker and worker.isRunning():
        worker.cancel_search()
        worker.wait()
    
    # Reset UI
    blast_button = parent.blastn_blast_button
    cancel_button = parent.blastn_cancel_button
    progress_bar = parent.blastn_progress_bar
    results_display = parent.blastn_results_display
    
    blast_button.setVisible(True)
    cancel_button.setVisible(False)
    progress_bar.setVisible(False)
    results_display.setText("⏹️ BLASTN search cancelled.")


def on_ncbi_blastn_finished(parent, result):
    """Handle successful BLASTN completion."""
    
    formatted_result = format_blast_output(result)
    results_display = parent.blastn_results_display
    results_display.setText(formatted_result)
    
    # Reset UI
    blast_button = parent.blastn_blast_button
    cancel_button = parent.blastn_cancel_button
    progress_bar = parent.blastn_progress_bar
    
    blast_button.setVisible(True)
    cancel_button.setVisible(False)
    progress_bar.setVisible(False)


def on_ncbi_blastn_error(parent, error):
    """Handle BLASTN error."""
    
    results_display = parent.blastn_results_display
    results_display.setText(f"❌ BLASTN search failed:\n\n{error}")
    
    # Reset UI
    blast_button = parent.blastn_blast_button
    cancel_button = parent.blastn_cancel_button
    progress_bar = parent.blastn_progress_bar
    
    blast_button.setVisible(True)
    cancel_button.setVisible(False)
    progress_bar.setVisible(False)


def on_ncbi_blastn_progress(parent, message):
    """Handle BLASTN progress updates."""
    
    results_display = parent.blastn_results_display
    results_display.setText(f"⏳ {message}")


# Utility functions (reuse from original blast_utils)
def open_blastn_file(parent):
    """Open a FASTA file for BLASTN."""
    
    file_path, _ = QFileDialog.getOpenFileName(
        parent, "Open Sequence File", "", "FASTA Files (*.fasta *.fa *.fas);;All Files (*)"
    )
    
    if file_path:
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            parent.blastn_file_label.setText(os.path.basename(file_path))
            parent.blastn_query_input.setText(content)
            
        except Exception as e:
            QMessageBox.critical(parent, "File Error", f"Could not read file: {str(e)}")


def save_blastn_results(parent):
    """Save BLASTN results to file."""
    
    if not parent.blastn_results_display.toPlainText().strip():
        QMessageBox.warning(parent, "No Results", "No results to save.")
        return
    
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    default_filename = f"blastn_results_{timestamp}.txt"
    
    file_path, _ = QFileDialog.getSaveFileName(
        parent, "Save BLASTN Results", default_filename, "Text Files (*.txt);;All Files (*)"
    )
    
    if file_path:
        try:
            with open(file_path, 'w') as f:
                f.write(parent.blastn_results_display.toPlainText())
            QMessageBox.information(parent, "Results Saved", f"Results saved to: {file_path}")
        except Exception as e:
            QMessageBox.critical(parent, "Save Error", f"Failed to save: {str(e)}")


def open_ncbi_blastn(parent):
    """Open NCBI BLASTN website."""
    
    query = parent.blastn_query_input.toPlainText().strip()
    
    if not query:
        QMessageBox.warning(parent, "No Query", "Please enter a sequence first.")
        return
    
    clean_query = clean_fasta_sequence(query)
    
    if len(clean_query) > 8000:
        ncbi_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch"
    else:
        encoded_query = quote(clean_query)
        ncbi_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&QUERY={encoded_query}"
    
    try:
        webbrowser.open_new_tab(ncbi_url)
    except Exception as e:
        QMessageBox.warning(parent, "Browser Error", f"Could not open browser: {str(e)}")