"""
BLAST utilities for PicoMol - NCBI-style layout matching the official BLAST interface.

This module recreates the exact layout and functionality of NCBI's BLAST interface
as seen in blast.html, providing an authentic user experience.
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


def create_ncbi_style_blastp_tab(parent):
    """Create BLASTP tab with exact NCBI layout."""
    
    # Main widget
    main_widget = QWidget()
    main_layout = QVBoxLayout(main_widget)
    main_layout.setContentsMargins(10, 10, 10, 10)
    main_layout.setSpacing(0)
    
    # Initialize worker
    parent.blastp_worker = None
    
    # Page title
    title_label = QLabel("Standard Protein BLAST")
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
    query_section = create_query_section(parent, "blastp")
    form_layout.addWidget(query_section)
    
    # 2. Choose Search Set Section  
    search_set_section = create_search_set_section(parent, "blastp")
    form_layout.addWidget(search_set_section)
    
    # 3. Program Selection Section
    program_section = create_program_selection_section(parent, "blastp")
    form_layout.addWidget(program_section)
    
    # 4. Algorithm Parameters Section (collapsible)
    params_section = create_algorithm_parameters_section(parent, "blastp")
    form_layout.addWidget(params_section)
    
    # Add stretch
    form_layout.addStretch()
    
    # Set form widget to scroll area
    scroll_area.setWidget(form_widget)
    main_layout.addWidget(scroll_area)
    
    # Results section
    results_section = create_results_section(parent, "blastp")
    main_layout.addWidget(results_section)
    
    return main_widget


def create_query_section(parent, program_type):
    """Create the 'Enter Query Sequence' section exactly like NCBI."""
    
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
    clear_button.clicked.connect(lambda: getattr(parent, f'{program_type}_query_input').clear())
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
    query_input.setPlaceholderText("Enter protein sequence here...")
    setattr(parent, f'{program_type}_query_input', query_input)
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
    setattr(parent, f'{program_type}_query_from', from_input)
    
    to_label = QLabel("To")
    to_input = QLineEdit()
    to_input.setFixedWidth(80)
    to_input.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    setattr(parent, f'{program_type}_query_to', to_input)
    
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
    upload_button.clicked.connect(lambda: open_blast_file(parent, program_type))
    upload_layout.addWidget(upload_button)
    
    file_label = QLabel("No file selected")
    file_label.setStyleSheet("color: #666; margin-left: 10px;")
    setattr(parent, f'{program_type}_file_label', file_label)
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
    setattr(parent, f'{program_type}_job_title', job_input)
    job_layout.addWidget(job_input)
    job_layout.addStretch()
    
    layout.addLayout(job_layout)
    
    return section


def create_search_set_section(parent, program_type):
    """Create the 'Choose Search Set' section exactly like NCBI."""
    
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
    
    # Database selection with radio buttons (like NCBI)
    db_group_layout = QVBoxLayout()
    
    # Standard databases radio button
    standard_radio = QRadioButton("Standard databases (nr etc.):")
    standard_radio.setChecked(True)
    standard_radio.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    setattr(parent, f'{program_type}_standard_radio', standard_radio)
    db_group_layout.addWidget(standard_radio)
    
    # ClusteredNR radio button with recommendation badge
    clustered_layout = QHBoxLayout()
    clustered_radio = QRadioButton("ClusteredNR database")
    clustered_radio.setStyleSheet("font-weight: normal;")
    setattr(parent, f'{program_type}_clustered_radio', clustered_radio)
    clustered_layout.addWidget(clustered_radio)
    
    recommended_badge = QLabel("Recommended")
    recommended_badge.setStyleSheet("""
        QLabel {
            background-color: #28a745;
            color: white;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 11px;
            font-weight: bold;
            margin-left: 5px;
        }
    """)
    clustered_layout.addWidget(recommended_badge)
    
    learn_more = QLabel('<a href="#" style="color: #007cba;">Learn more...</a>')
    learn_more.setOpenExternalLinks(False)
    clustered_layout.addWidget(learn_more)
    clustered_layout.addStretch()
    
    db_group_layout.addLayout(clustered_layout)
    
    layout.addLayout(db_group_layout)
    
    # Database dropdown
    db_layout = QHBoxLayout()
    db_label = QLabel("Database")
    db_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
    db_layout.addWidget(db_label)
    
    database_combo = QComboBox()
    database_combo.setFixedWidth(400)
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
    
    # Add protein databases
    protein_databases = [
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
    ]
    
    for db_id, db_name in protein_databases:
        database_combo.addItem(db_name, db_id)
    
    setattr(parent, f'{program_type}_database_combo', database_combo)
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
    org_label = QLabel("Organism")
    org_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
    org_layout.addWidget(org_label)
    
    organism_input = QLineEdit()
    organism_input.setFixedWidth(300)
    organism_input.setStyleSheet("border: 1px solid #ccc; padding: 4px;")
    organism_input.setPlaceholderText("Enter organism common name, binomial, or tax id")
    setattr(parent, f'{program_type}_organism_input', organism_input)
    org_layout.addWidget(organism_input)
    
    exclude_checkbox = QCheckBox("exclude")
    exclude_checkbox.setStyleSheet("margin-left: 10px;")
    setattr(parent, f'{program_type}_exclude_organism', exclude_checkbox)
    org_layout.addWidget(exclude_checkbox)
    org_layout.addStretch()
    
    layout.addLayout(org_layout)
    
    # Exclude options
    exclude_layout = QVBoxLayout()
    exclude_label = QLabel("Exclude")
    exclude_label.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    exclude_layout.addWidget(exclude_label)
    
    exclude_models = QCheckBox("Models (XM/XP)")
    exclude_models.setStyleSheet("margin-left: 20px;")
    setattr(parent, f'{program_type}_exclude_models', exclude_models)
    exclude_layout.addWidget(exclude_models)
    
    wp_proteins = QCheckBox("Non-redundant RefSeq proteins (WP)")
    wp_proteins.setStyleSheet("margin-left: 20px;")
    setattr(parent, f'{program_type}_wp_proteins', wp_proteins)
    exclude_layout.addWidget(wp_proteins)
    
    exclude_uncult = QCheckBox("Uncultured/environmental sample sequences")
    exclude_uncult.setStyleSheet("margin-left: 20px;")
    setattr(parent, f'{program_type}_exclude_uncult', exclude_uncult)
    exclude_layout.addWidget(exclude_uncult)
    
    layout.addLayout(exclude_layout)
    
    return section


def create_program_selection_section(parent, program_type):
    """Create the 'Program Selection' section exactly like NCBI."""
    
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
    
    # Algorithm selection
    algo_label = QLabel("Algorithm")
    algo_label.setStyleSheet("font-weight: normal; margin-bottom: 10px;")
    layout.addWidget(algo_label)
    
    # Radio buttons for algorithms
    algo_group = QButtonGroup()
    algo_layout = QVBoxLayout()
    
    algorithms = [
        ("kmerBlastp", "Quick BLASTP (Accelerated protein-protein BLAST)"),
        ("blastp", "blastp (protein-protein BLAST)"),
        ("psiBlast", "PSI-BLAST (Position-Specific Iterated BLAST)"),
        ("phiBlast", "PHI-BLAST (Pattern Hit Initiated BLAST)"),
        ("deltaBlast", "DELTA-BLAST (Domain Enhanced Lookup Time Accelerated BLAST)")
    ]
    
    for i, (algo_id, algo_desc) in enumerate(algorithms):
        radio = QRadioButton(algo_desc)
        radio.setStyleSheet("margin-bottom: 5px;")
        if algo_id == "blastp":
            radio.setChecked(True)
        algo_group.addButton(radio, i)
        setattr(parent, f'{program_type}_{algo_id}_radio', radio)
        algo_layout.addWidget(radio)
    
    setattr(parent, f'{program_type}_algorithm_group', algo_group)
    layout.addLayout(algo_layout)
    
    # PHI pattern input (initially hidden)
    phi_layout = QHBoxLayout()
    phi_label = QLabel("Enter a PHI pattern")
    phi_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
    phi_layout.addWidget(phi_label)
    
    phi_input = QLineEdit()
    phi_input.setFixedWidth(300)
    phi_input.setStyleSheet("border: 1px solid #ccc; padding: 4px;")
    phi_input.setVisible(False)
    setattr(parent, f'{program_type}_phi_pattern', phi_input)
    phi_layout.addWidget(phi_input)
    phi_layout.addStretch()
    
    layout.addLayout(phi_layout)
    
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
    blast_button.clicked.connect(lambda: run_ncbi_blast_search(parent, program_type))
    setattr(parent, f'{program_type}_blast_button', blast_button)
    
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
    cancel_button.clicked.connect(lambda: cancel_ncbi_blast_search(parent, program_type))
    setattr(parent, f'{program_type}_cancel_button', cancel_button)
    
    button_layout = QHBoxLayout()
    button_layout.addWidget(blast_button)
    button_layout.addWidget(cancel_button)
    button_layout.addStretch()
    
    search_layout.addLayout(button_layout)
    
    # Progress bar
    progress_bar = QProgressBar()
    progress_bar.setVisible(False)
    progress_bar.setRange(0, 0)  # Indeterminate
    setattr(parent, f'{program_type}_progress_bar', progress_bar)
    search_layout.addWidget(progress_bar)
    
    # Search summary
    summary_layout = QVBoxLayout()
    summary_label = QLabel("Search nr using Blastp (protein-protein BLAST)")
    summary_label.setStyleSheet("color: #666; font-size: 12px; margin-top: 10px;")
    setattr(parent, f'{program_type}_summary_label', summary_label)
    summary_layout.addWidget(summary_label)
    
    new_window_checkbox = QCheckBox("Show results in a new window")
    new_window_checkbox.setStyleSheet("font-size: 12px; margin-top: 5px;")
    setattr(parent, f'{program_type}_new_window', new_window_checkbox)
    summary_layout.addWidget(new_window_checkbox)
    
    search_layout.addLayout(summary_layout)
    layout.addLayout(search_layout)
    
    return section


def create_algorithm_parameters_section(parent, program_type):
    """Create the collapsible 'Algorithm parameters' section exactly like NCBI."""
    
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
        visible = params_content.isVisible()
        params_content.setVisible(not visible)
        toggle_button.setText("▲ Algorithm parameters" if not visible else "▼ Algorithm parameters")
    
    toggle_button.clicked.connect(toggle_params)
    setattr(parent, f'{program_type}_toggle_params', toggle_button)
    setattr(parent, f'{program_type}_params_content', params_content)
    
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
    general_group = create_general_parameters_group(parent, program_type)
    params_layout.addWidget(general_group)
    
    # Scoring Parameters
    scoring_group = create_scoring_parameters_group(parent, program_type)
    params_layout.addWidget(scoring_group)
    
    # Filters and Masking
    filters_group = create_filters_masking_group(parent, program_type)
    params_layout.addWidget(filters_group)
    
    layout.addWidget(params_content)
    
    return section


def create_general_parameters_group(parent, program_type):
    """Create General Parameters group."""
    
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
    setattr(parent, f'{program_type}_max_targets', max_combo)
    max_layout.addWidget(max_combo)
    max_layout.addStretch()
    
    layout.addLayout(max_layout)
    
    # Short queries
    short_layout = QHBoxLayout()
    short_checkbox = QCheckBox("Automatically adjust parameters for short input sequences")
    short_checkbox.setChecked(True)
    short_checkbox.setStyleSheet("margin: 5px 0;")
    setattr(parent, f'{program_type}_short_adjust', short_checkbox)
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
    setattr(parent, f'{program_type}_expect_threshold', expect_input)
    expect_layout.addWidget(expect_input)
    expect_layout.addStretch()
    
    layout.addLayout(expect_layout)
    
    # Word size
    word_layout = QHBoxLayout()
    word_label = QLabel("Word size")
    word_label.setFixedWidth(150)
    word_layout.addWidget(word_label)
    
    word_combo = QComboBox()
    word_combo.setFixedWidth(100)
    word_combo.addItems(["2", "3", "5", "6"])
    word_combo.setCurrentText("3")
    word_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    setattr(parent, f'{program_type}_word_size', word_combo)
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
    setattr(parent, f'{program_type}_hsp_range_max', range_input)
    range_layout.addWidget(range_input)
    range_layout.addStretch()
    
    layout.addLayout(range_layout)
    
    return group


def create_scoring_parameters_group(parent, program_type):
    """Create Scoring Parameters group."""
    
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
    
    # Matrix
    matrix_layout = QHBoxLayout()
    matrix_label = QLabel("Matrix")
    matrix_label.setFixedWidth(150)
    matrix_layout.addWidget(matrix_label)
    
    matrix_combo = QComboBox()
    matrix_combo.setFixedWidth(120)
    matrices = ["PAM30", "PAM70", "BLOSUM90", "BLOSUM80", "BLOSUM62", "BLOSUM50", "BLOSUM45", "PAM250"]
    matrix_combo.addItems(matrices)
    matrix_combo.setCurrentText("BLOSUM62")
    matrix_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    setattr(parent, f'{program_type}_matrix', matrix_combo)
    matrix_layout.addWidget(matrix_combo)
    matrix_layout.addStretch()
    
    layout.addLayout(matrix_layout)
    
    # Gap costs
    gap_layout = QHBoxLayout()
    gap_label = QLabel("Gap Costs")
    gap_label.setFixedWidth(150)
    gap_layout.addWidget(gap_label)
    
    gap_combo = QComboBox()
    gap_combo.setFixedWidth(200)
    gap_costs = [
        "Existence: 11 Extension: 2",
        "Existence: 10 Extension: 2", 
        "Existence: 9 Extension: 2",
        "Existence: 8 Extension: 2",
        "Existence: 7 Extension: 2",
        "Existence: 6 Extension: 2",
        "Existence: 13 Extension: 1",
        "Existence: 12 Extension: 1",
        "Existence: 11 Extension: 1",
        "Existence: 10 Extension: 1",
        "Existence: 9 Extension: 1"
    ]
    gap_combo.addItems(gap_costs)
    gap_combo.setCurrentText("Existence: 11 Extension: 1")
    gap_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    setattr(parent, f'{program_type}_gap_costs', gap_combo)
    gap_layout.addWidget(gap_combo)
    gap_layout.addStretch()
    
    layout.addLayout(gap_layout)
    
    # Compositional adjustments
    comp_layout = QHBoxLayout()
    comp_label = QLabel("Compositional adjustments")
    comp_label.setFixedWidth(150)
    comp_layout.addWidget(comp_label)
    
    comp_combo = QComboBox()
    comp_combo.setFixedWidth(300)
    comp_adjustments = [
        "No adjustment",
        "Composition-based statistics",
        "Conditional compositional score matrix adjustment",
        "Universal compositional score matrix adjustment"
    ]
    comp_combo.addItems(comp_adjustments)
    comp_combo.setCurrentText("Conditional compositional score matrix adjustment")
    comp_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    setattr(parent, f'{program_type}_comp_adjust', comp_combo)
    comp_layout.addWidget(comp_combo)
    comp_layout.addStretch()
    
    layout.addLayout(comp_layout)
    
    return group


def create_filters_masking_group(parent, program_type):
    """Create Filters and Masking group."""
    
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
    setattr(parent, f'{program_type}_low_complexity', low_comp_checkbox)
    low_comp_layout.addWidget(low_comp_checkbox)
    low_comp_layout.addStretch()
    
    layout.addLayout(low_comp_layout)
    
    return group


def create_results_section(parent, program_type):
    """Create the results display section."""
    
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
    title = QLabel("BLAST Results")
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
    save_button.clicked.connect(lambda: save_blast_results(parent, program_type))
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
    clear_button.clicked.connect(lambda: getattr(parent, f'{program_type}_results_display').clear())
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
    ncbi_button.clicked.connect(lambda: open_ncbi_blast(parent, program_type))
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
        "BLAST results will appear here...\n\n"
        "• Enter a protein sequence above\n"
        "• Configure your search parameters\n"
        "• Click the BLAST button to begin\n"
        "• Results typically take 30 seconds to 2 minutes"
    )
    setattr(parent, f'{program_type}_results_display', results_display)
    layout.addWidget(results_display)
    
    return section


# Helper functions
def run_ncbi_blast_search(parent, program_type):
    """Run BLAST search with NCBI-style interface."""
    
    query_input = getattr(parent, f'{program_type}_query_input')
    query = query_input.toPlainText().strip()
    
    if not query:
        QMessageBox.warning(parent, "No Query", "Please enter a sequence to search.")
        return
    
    # Validate sequence
    sequence_type = "protein" if program_type in ["blastp", "tblastn"] else "nucleotide"
    is_valid, message = validate_sequence(query, sequence_type)
    if not is_valid:
        QMessageBox.warning(parent, "Invalid Sequence", message)
        return
    
    # Get parameters
    database_combo = getattr(parent, f'{program_type}_database_combo')
    database = database_combo.currentData()
    
    expect_input = getattr(parent, f'{program_type}_expect_threshold')
    max_targets = getattr(parent, f'{program_type}_max_targets')
    word_size = getattr(parent, f'{program_type}_word_size')
    
    parameters = {
        'evalue': expect_input.text().strip() or "0.05",
        'max_target_seqs': max_targets.currentText(),
        'word_size': word_size.currentText(),
    }
    
    # Add protein-specific parameters
    if program_type == "blastp":
        matrix = getattr(parent, f'{program_type}_matrix')
        gap_costs = getattr(parent, f'{program_type}_gap_costs')
        comp_adjust = getattr(parent, f'{program_type}_comp_adjust')
        
        parameters['matrix'] = matrix.currentText()
        
        # Parse gap costs
        gap_text = gap_costs.currentText()
        if "Existence:" in gap_text:
            parts = gap_text.replace("Existence:", "").replace("Extension:", "").split()
            if len(parts) >= 2:
                parameters['gapopen'] = parts[0]
                parameters['gapextend'] = parts[1]
        
        # Map compositional adjustment
        comp_map = {
            "No adjustment": "0",
            "Composition-based statistics": "1", 
            "Conditional compositional score matrix adjustment": "2",
            "Universal compositional score matrix adjustment": "3"
        }
        parameters['comp_adjust'] = comp_map.get(comp_adjust.currentText(), "2")
    
    # Get organism filter
    organism_input = getattr(parent, f'{program_type}_organism_input')
    if organism_input.text().strip():
        parameters['organism'] = organism_input.text().strip()
    
    # Get algorithm
    algorithm_group = getattr(parent, f'{program_type}_algorithm_group')
    checked_button = algorithm_group.checkedButton()
    if checked_button:
        button_text = checked_button.text()
        if "Quick BLASTP" in button_text:
            parameters['algorithm'] = 'kmerBlastp'
        elif "PSI-BLAST" in button_text:
            parameters['algorithm'] = 'psiBlast'
        elif "PHI-BLAST" in button_text:
            parameters['algorithm'] = 'phiBlast'
        elif "DELTA-BLAST" in button_text:
            parameters['algorithm'] = 'deltaBlast'
        else:
            parameters['algorithm'] = 'blastp'
    
    # Update UI
    blast_button = getattr(parent, f'{program_type}_blast_button')
    cancel_button = getattr(parent, f'{program_type}_cancel_button')
    progress_bar = getattr(parent, f'{program_type}_progress_bar')
    results_display = getattr(parent, f'{program_type}_results_display')
    
    blast_button.setVisible(False)
    cancel_button.setVisible(True)
    progress_bar.setVisible(True)
    results_display.setText("🚀 Submitting BLAST search to NCBI...")
    
    # Clean sequence
    clean_sequence = clean_fasta_sequence(query)
    
    # Create and start worker
    worker = OnlineBlastWorker(clean_sequence, database, program_type, parameters)
    setattr(parent, f'{program_type}_worker', worker)
    
    worker.finished.connect(lambda result: on_ncbi_blast_finished(parent, program_type, result))
    worker.error.connect(lambda error: on_ncbi_blast_error(parent, program_type, error))
    worker.progress.connect(lambda msg: on_ncbi_blast_progress(parent, program_type, msg))
    worker.start()


def cancel_ncbi_blast_search(parent, program_type):
    """Cancel BLAST search."""
    
    worker = getattr(parent, f'{program_type}_worker', None)
    if worker and worker.isRunning():
        worker.cancel_search()
        worker.wait()
    
    # Reset UI
    blast_button = getattr(parent, f'{program_type}_blast_button')
    cancel_button = getattr(parent, f'{program_type}_cancel_button')
    progress_bar = getattr(parent, f'{program_type}_progress_bar')
    results_display = getattr(parent, f'{program_type}_results_display')
    
    blast_button.setVisible(True)
    cancel_button.setVisible(False)
    progress_bar.setVisible(False)
    results_display.setText("⏹️ BLAST search cancelled.")


def on_ncbi_blast_finished(parent, program_type, result):
    """Handle successful BLAST completion."""
    
    formatted_result = format_blast_output(result)
    results_display = getattr(parent, f'{program_type}_results_display')
    results_display.setText(formatted_result)
    
    # Reset UI
    blast_button = getattr(parent, f'{program_type}_blast_button')
    cancel_button = getattr(parent, f'{program_type}_cancel_button')
    progress_bar = getattr(parent, f'{program_type}_progress_bar')
    
    blast_button.setVisible(True)
    cancel_button.setVisible(False)
    progress_bar.setVisible(False)


def on_ncbi_blast_error(parent, program_type, error):
    """Handle BLAST error."""
    
    results_display = getattr(parent, f'{program_type}_results_display')
    results_display.setText(f"❌ BLAST search failed:\n\n{error}")
    
    # Reset UI
    blast_button = getattr(parent, f'{program_type}_blast_button')
    cancel_button = getattr(parent, f'{program_type}_cancel_button')
    progress_bar = getattr(parent, f'{program_type}_progress_bar')
    
    blast_button.setVisible(True)
    cancel_button.setVisible(False)
    progress_bar.setVisible(False)


def on_ncbi_blast_progress(parent, program_type, message):
    """Handle BLAST progress updates."""
    
    results_display = getattr(parent, f'{program_type}_results_display')
    results_display.setText(f"⏳ {message}")


# Utility functions (reuse from original blast_utils)
def open_blast_file(parent, program_type):
    """Open a FASTA file."""
    
    file_path, _ = QFileDialog.getOpenFileName(
        parent, "Open Sequence File", "", "FASTA Files (*.fasta *.fa *.fas);;All Files (*)"
    )
    
    if file_path:
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            file_label = getattr(parent, f'{program_type}_file_label')
            query_input = getattr(parent, f'{program_type}_query_input')
            
            file_label.setText(os.path.basename(file_path))
            query_input.setText(content)
            
        except Exception as e:
            QMessageBox.critical(parent, "File Error", f"Could not read file: {str(e)}")


def save_blast_results(parent, program_type):
    """Save BLAST results to file."""
    
    results_display = getattr(parent, f'{program_type}_results_display')
    if not results_display.toPlainText().strip():
        QMessageBox.warning(parent, "No Results", "No results to save.")
        return
    
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    default_filename = f"blast_results_{timestamp}.txt"
    
    file_path, _ = QFileDialog.getSaveFileName(
        parent, "Save BLAST Results", default_filename, "Text Files (*.txt);;All Files (*)"
    )
    
    if file_path:
        try:
            with open(file_path, 'w') as f:
                f.write(results_display.toPlainText())
            QMessageBox.information(parent, "Results Saved", f"Results saved to: {file_path}")
        except Exception as e:
            QMessageBox.critical(parent, "Save Error", f"Failed to save: {str(e)}")


def open_ncbi_blast(parent, program_type):
    """Open NCBI BLAST website."""
    
    query_input = getattr(parent, f'{program_type}_query_input')
    query = query_input.toPlainText().strip()
    
    if not query:
        QMessageBox.warning(parent, "No Query", "Please enter a sequence first.")
        return
    
    clean_query = clean_fasta_sequence(query)
    
    if len(clean_query) > 8000:
        ncbi_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM={program_type}&PAGE_TYPE=BlastSearch"
    else:
        encoded_query = quote(clean_query)
        ncbi_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM={program_type}&PAGE_TYPE=BlastSearch&QUERY={encoded_query}"
    
    try:
        webbrowser.open_new_tab(ncbi_url)
    except Exception as e:
        QMessageBox.warning(parent, "Browser Error", f"Could not open browser: {str(e)}")