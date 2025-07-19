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
from .blast_utils import OnlineBlastWorker, validate_sequence, clean_fasta_sequence, format_blast_output
from .accession_validator import is_accession_number


def create_ncbi_style_blastp_tab(parent):
    """Create BLASTP tab with exact NCBI layout."""
    return create_ncbi_style_blast_tab(parent, "blastp", "Standard Protein BLAST", "protein")


def create_ncbi_style_blastn_tab(parent):
    """Create BLASTN tab with exact NCBI layout."""
    return create_ncbi_style_blastn_tab_unified(parent)


def create_ncbi_style_blastx_tab(parent):
    """Create BLASTX tab with exact NCBI layout."""
    return create_ncbi_style_blast_tab(parent, "blastx", "Translated BLAST: blastx", "nucleotide")


def create_ncbi_style_tblastn_tab(parent):
    """Create TBLASTN tab with exact NCBI layout."""
    return create_ncbi_style_blast_tab(parent, "tblastn", "Translated BLAST: tblastn", "protein")


def create_ncbi_style_tblastx_tab(parent):
    """Create TBLASTX tab with exact NCBI layout."""
    return create_ncbi_style_blast_tab(parent, "tblastx", "Translated BLAST: tblastx", "nucleotide")


def create_ncbi_style_blast_tab(parent, program_type, title, sequence_type):
    """Create a BLAST tab with exact NCBI layout for any program type."""
    
    # Main widget
    main_widget = QWidget()
    main_layout = QVBoxLayout(main_widget)
    main_layout.setContentsMargins(10, 10, 10, 10)
    main_layout.setSpacing(0)
    
    # Initialize worker
    setattr(parent, f'{program_type}_worker', None)
    
    # Page title
    title_label = QLabel(title)
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
    
    # Program description for all programs
    if program_type == "blastp":
        desc_label = QLabel("BLASTP searches protein databases using a protein query.")
    elif program_type == "blastn":
        desc_label = QLabel("BLASTN searches nucleotide databases using a nucleotide query.")
    elif program_type == "blastx":
        desc_label = QLabel("BLASTX searches protein databases using a translated nucleotide query.")
    elif program_type == "tblastn":
        desc_label = QLabel("TBLASTN searches translated nucleotide databases using a protein query.")
    elif program_type == "tblastx":
        desc_label = QLabel("TBLASTX searches translated nucleotide databases using a translated nucleotide query.")
    else:
        desc_label = QLabel(f"{program_type.upper()} sequence search.")
    
    desc_label.setStyleSheet("""
        QLabel {
            font-size: 14px;
            color: #666;
            padding: 5px 0;
            margin-bottom: 10px;
        }
    """)
    main_layout.addWidget(desc_label)
    
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
    query_section = create_query_section(parent, program_type, sequence_type)
    form_layout.addWidget(query_section)
    
    # 2. Subject Sequence Section (for BL2SEQ mode - initially hidden)
    subject_section = create_subject_sequence_section(parent, program_type, sequence_type)
    subject_section.setVisible(False)  # Initially hidden
    setattr(parent, f'{program_type}_subject_section', subject_section)
    form_layout.addWidget(subject_section)
    
    # 3. Choose Search Set Section  
    search_set_section = create_search_set_section(parent, program_type, sequence_type)
    form_layout.addWidget(search_set_section)
    
    # 4. Program Selection Section
    program_section = create_program_selection_section(parent, program_type)
    form_layout.addWidget(program_section)
    
    # 5. Algorithm Parameters Section (collapsible)
    params_section = create_algorithm_parameters_section(parent, program_type)
    form_layout.addWidget(params_section)
    
    # Add stretch
    form_layout.addStretch()
    
    # Set form widget to scroll area
    scroll_area.setWidget(form_widget)
    main_layout.addWidget(scroll_area)
    
    # Results section
    results_section = create_results_section(parent, program_type)
    main_layout.addWidget(results_section)
    
    return main_widget


def create_query_section(parent, program_type, sequence_type="protein"):
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
    placeholder_text = f"Enter {sequence_type} sequence here..."
    if sequence_type == "nucleotide":
        placeholder_text += "\n\nExample:\n>My_sequence\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    else:
        placeholder_text += "\n\nExample:\n>My_protein\nMKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
    
    query_input.setPlaceholderText(placeholder_text)
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
    
    # BL2SEQ (Align two sequences) option
    bl2seq_layout = QHBoxLayout()
    bl2seq_checkbox = QCheckBox("Align two or more sequences")
    bl2seq_checkbox.setStyleSheet("font-weight: normal; margin: 10px 0;")
    bl2seq_checkbox.setToolTip("Enable pairwise sequence alignment mode")
    setattr(parent, f'{program_type}_bl2seq', bl2seq_checkbox)
    
    # Connect checkbox to show/hide subject sequence section
    def toggle_subject_section():
        """Show/hide subject sequence section based on BL2SEQ checkbox."""
        subject_section = getattr(parent, f'{program_type}_subject_section', None)
        if subject_section:
            subject_section.setVisible(bl2seq_checkbox.isChecked())
    
    bl2seq_checkbox.toggled.connect(toggle_subject_section)
    setattr(parent, f'{program_type}_toggle_subject_section', toggle_subject_section)
    
    bl2seq_layout.addWidget(bl2seq_checkbox)
    bl2seq_layout.addStretch()
    
    layout.addLayout(bl2seq_layout)
    
    # Genetic Code (for translated searches - in query section like NCBI)
    # Note: tblastn does NOT show genetic code in query section according to HTML
    if program_type in ["blastx", "tblastx"]:
        genetic_layout = QHBoxLayout()
        genetic_label = QLabel("Genetic code")
        genetic_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
        genetic_layout.addWidget(genetic_label)
        
        genetic_combo = QComboBox()
        genetic_combo.setFixedWidth(250)
        genetic_codes = [
            ("1", "Standard (1)"),
            ("2", "Vertebrate Mitochondrial (2)"),
            ("3", "Yeast Mitochondrial (3)"),
            ("4", "Mold Mitochondrial; ... (4)"),
            ("5", "Invertebrate Mitochondrial (5)"),
            ("6", "Ciliate Nuclear; ... (6)"),
            ("9", "Echinoderm Mitochondrial (9)"),
            ("10", "Euplotid Nuclear (10)"),
            ("11", "Bacteria and Archaea (11)"),
            ("12", "Alternative Yeast Nuclear (12)"),
            ("13", "Ascidian Mitochondrial (13)"),
            ("14", "Flatworm Mitochondrial (14)"),
            ("16", "Chlorophycean Mitochondrial (16)"),
            ("21", "Trematode Mitochondrial (21)"),
            ("22", "Scenedesmus obliquus Mitochondrial (22)"),
            ("23", "Thraustochytrium Mitochondrial (23)"),
            ("24", "Pterobranchia Mitochondrial (24)"),
            ("25", "Candidate Division SR1 and Gracilibacteria (25)"),
            ("26", "Pachysolen tannophilus Nuclear (26)"),
            ("27", "Karyorelict Nuclear (27)"),
            ("28", "Condylostoma Nuclear (28)"),
            ("29", "Mesodinium Nuclear (29)"),
            ("30", "Peritrich Nuclear (30)"),
            ("31", "Blastocrithidia Nuclear (31)")
        ]
        for code, desc in genetic_codes:
            genetic_combo.addItem(desc, code)
        genetic_combo.setCurrentText("Standard (1)")
        genetic_combo.setStyleSheet("border: 1px solid #ccc; padding: 4px;")
        setattr(parent, f'{program_type}_genetic_code', genetic_combo)
        genetic_layout.addWidget(genetic_combo)
        genetic_layout.addStretch()
        
        layout.addLayout(genetic_layout)
    
    return section


def create_subject_sequence_section(parent, program_type, sequence_type="protein"):
    """Create the 'Enter Subject Sequence' section for BL2SEQ mode."""
    
    section = QFrame()
    section.setFrameStyle(QFrame.Box)
    section.setStyleSheet("""
        QFrame {
            border: 1px solid #ccc;
            border-radius: 4px;
            background-color: #fff3cd;
        }
    """)
    
    layout = QVBoxLayout(section)
    layout.setContentsMargins(15, 10, 15, 15)
    
    # Legend/Title
    legend = QLabel("Enter Subject Sequence (BL2SEQ Mode)")
    legend.setStyleSheet("""
        QLabel {
            font-weight: bold;
            font-size: 14px;
            color: #333;
            background-color: #fff3cd;
            padding: 0 5px;
            margin-bottom: 10px;
        }
    """)
    layout.addWidget(legend)
    
    # Subject input area
    subject_layout = QVBoxLayout()
    
    # Label with help link
    label_layout = QHBoxLayout()
    subject_label = QLabel("Enter accession number(s), gi(s), or FASTA sequence(s)")
    subject_label.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    label_layout.addWidget(subject_label)
    
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
    help_button.setToolTip("How to enter subject sequences")
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
    clear_button.clicked.connect(lambda: getattr(parent, f'{program_type}_subject_input').clear())
    label_layout.addWidget(clear_button)
    label_layout.addStretch()
    
    subject_layout.addLayout(label_layout)
    
    # Text area
    subject_input = QTextEdit()
    subject_input.setFixedHeight(120)
    subject_input.setStyleSheet("""
        QTextEdit {
            border: 1px solid #ccc;
            font-family: monospace;
            font-size: 12px;
            background-color: white;
        }
    """)
    
    # Set appropriate placeholder based on sequence type
    if sequence_type == "nucleotide" or program_type in ["blastn", "blastx", "tblastx"]:
        placeholder_text = "Enter nucleotide subject sequence here...\n\nExample:\n>Subject_sequence\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    else:
        placeholder_text = "Enter protein subject sequence here...\n\nExample:\n>Subject_protein\nMKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
    
    subject_input.setPlaceholderText(placeholder_text)
    setattr(parent, f'{program_type}_subject_input', subject_input)
    subject_layout.addWidget(subject_input)
    
    layout.addLayout(subject_layout)
    
    # Subject subrange
    subrange_layout = QHBoxLayout()
    subrange_label = QLabel("Subject subrange")
    subrange_label.setStyleSheet("font-weight: normal; margin-right: 10px;")
    subrange_layout.addWidget(subrange_label)
    
    from_label = QLabel("From")
    from_input = QLineEdit()
    from_input.setFixedWidth(80)
    from_input.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    setattr(parent, f'{program_type}_subject_from', from_input)
    
    to_label = QLabel("To")
    to_input = QLineEdit()
    to_input.setFixedWidth(80)
    to_input.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    setattr(parent, f'{program_type}_subject_to', to_input)
    
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
    upload_button.clicked.connect(lambda: open_subject_file(parent, program_type))
    upload_layout.addWidget(upload_button)
    
    file_label = QLabel("No file selected")
    file_label.setStyleSheet("color: #666; margin-left: 10px;")
    setattr(parent, f'{program_type}_subject_file_label', file_label)
    upload_layout.addWidget(file_label)
    upload_layout.addStretch()
    
    layout.addLayout(upload_layout)
    
    return section


def create_search_set_section(parent, program_type, sequence_type="protein"):
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
    
    # Database selection with radio buttons exactly like NCBI HTML
    # CRITICAL: genomeDbs div has class="blastn megaBlast discoMegablast" - NOT for tblastn!
    # Only show for blastn (megaBlast, discoMegablast), NOT for tblastn or tblastx
    if program_type == "blastn":
        # Create genomeDbs div equivalent
        genomic_group = QFrame()
        genomic_group.setStyleSheet("border: none; margin: 10px 0;")
        genomic_layout = QVBoxLayout(genomic_group)
        genomic_layout.setContentsMargins(0, 0, 0, 0)
        
        # Database type radio button group
        db_type_group = QButtonGroup()
        
        # Standard databases radio button (Rgen)
        standard_radio = QRadioButton("Standard databases (nr etc.):")
        standard_radio.setChecked(True)
        standard_radio.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
        setattr(parent, f'{program_type}_Rgen', standard_radio)
        db_type_group.addButton(standard_radio, 0)
        genomic_layout.addWidget(standard_radio)
        
        # rRNA/ITS databases radio button (Rrna)
        rrna_radio = QRadioButton("rRNA/ITS databases")
        rrna_radio.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
        setattr(parent, f'{program_type}_Rrna', rrna_radio)
        db_type_group.addButton(rrna_radio, 1)
        genomic_layout.addWidget(rrna_radio)
        
        # Genomic + transcript databases radio button (Rhc)
        genomic_radio = QRadioButton("Genomic + transcript databases")
        genomic_radio.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
        setattr(parent, f'{program_type}_Rhc', genomic_radio)
        db_type_group.addButton(genomic_radio, 2)
        genomic_layout.addWidget(genomic_radio)
        
        # Betacoronavirus radio button (Rvir)
        beta_radio = QRadioButton("Betacoronavirus")
        beta_radio.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
        setattr(parent, f'{program_type}_Rvir', beta_radio)
        db_type_group.addButton(beta_radio, 3)
        genomic_layout.addWidget(beta_radio)
        
        setattr(parent, f'{program_type}_db_type_group', db_type_group)
        layout.addWidget(genomic_group)
        
        # rRNA Extra info (hidden by default)
        rrna_extra = QLabel('<a href="https://www.ncbi.nlm.nih.gov/refseq/targetedloci/" style="color: #007cba;">Targeted Loci Project Information</a>')
        rrna_extra.setOpenExternalLinks(True)
        rrna_extra.setVisible(False)
        setattr(parent, f'{program_type}_rrna_extra', rrna_extra)
        layout.addWidget(rrna_extra)
        
        # Connect radio buttons to show/hide extra info
        def update_rrna_extra():
            rrna_extra.setVisible(rrna_radio.isChecked())
        
        rrna_radio.toggled.connect(update_rrna_extra)
        
    elif program_type in ["blastp", "blastx"]:
        # For protein searches only, show ClusteredNR option
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
    # For tblastn and tblastx: NO radio buttons at all - just direct database dropdown
    
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
    
    # Add databases based on program type - exact match from HTML
    if program_type == "tblastn":
        # TBLASTN searches nucleotide databases using protein query - from HTML
        databases = [
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
            ("sts", "Sequence tagged sites (dbsts)"),
            ("genomic/Viruses/Betacoronavirus", "Betacoronavirus Genbank")
        ]
    elif program_type in ["blastp", "blastx"]:
        # Protein databases for protein searches only (blastp, blastx)
        databases = [
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
    elif program_type == "tblastx":
        # TBLASTX searches nucleotide databases using translated nucleotide query - from HTML
        databases = [
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
            ("sts", "Sequence tagged sites (dbsts)"),
            ("genomic/Viruses/Betacoronavirus", "Betacoronavirus Genbank")
        ]
    else:
        # Nucleotide databases for other nucleotide searches
        databases = [
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
    
    # Store all database options
    all_databases = {
        'standard': databases,
        'clustered': [("nr_cluster_seq", "ClusteredNR (nr_cluster_seq)")],
        'nucleotide': [
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
    
    # Function to update database list based on radio selection
    def update_database_list():
        """Update database dropdown based on selected radio button."""
        database_combo.clear()
        
        if program_type == "blastn":
            # For blastn, check which database type is selected
            db_type_group = getattr(parent, f'{program_type}_db_type_group', None)
            if db_type_group:
                checked_id = db_type_group.checkedId()
                if checked_id == 0:  # Standard (Rgen)
                    db_group = "Rgen"
                elif checked_id == 1:  # rRNA (Rrna)
                    db_group = "Rrna"
                elif checked_id == 2:  # Genomic + transcript (Rhc)
                    db_group = "Rhc"
                elif checked_id == 3:  # Betacoronavirus (Rvir)
                    db_group = "Rvir"
                else:
                    db_group = "Rgen"  # Default
            else:
                db_group = "Rgen"  # Default
            
            # Filter databases by group
            for db_id, db_name in all_databases['standard']:
                # For now, show all standard databases (would need dbgroup mapping for full filtering)
                database_combo.addItem(db_name, db_id)
        elif program_type in ["tblastn", "tblastx"]:
            # For tblastn/tblastx, no radio buttons - just show all nucleotide databases
            for db_id, db_name in all_databases['standard']:
                database_combo.addItem(db_name, db_id)
        else:
            # For protein searches
            if hasattr(parent, f'{program_type}_clustered_radio') and getattr(parent, f'{program_type}_clustered_radio').isChecked():
                # ClusteredNR selected
                for db_id, db_name in all_databases['clustered']:
                    database_combo.addItem(db_name, db_id)
            else:
                # Standard databases selected
                for db_id, db_name in all_databases['standard']:
                    database_combo.addItem(db_name, db_id)
    
    # Initial population
    for db_id, db_name in databases:
        database_combo.addItem(db_name, db_id)
    
    # Set default database based on program type (from HTML defVal)
    if program_type in ["tblastn", "tblastx"]:
        # Find and set "nt" as default for tblastn and tblastx
        for i in range(database_combo.count()):
            if database_combo.itemData(i) == "nt":
                database_combo.setCurrentIndex(i)
                break
    
    # Connect radio buttons to update function
    if program_type == "blastn":
        # Connect nucleotide database type radio buttons (only for blastn)
        db_type_group = getattr(parent, f'{program_type}_db_type_group', None)
        if db_type_group:
            for button in db_type_group.buttons():
                button.toggled.connect(update_database_list)
    elif program_type in ["blastp", "blastx"]:
        # Connect protein database radio buttons
        standard_radio.toggled.connect(update_database_list)
        if hasattr(parent, f'{program_type}_clustered_radio'):
            getattr(parent, f'{program_type}_clustered_radio').toggled.connect(update_database_list)
    
    setattr(parent, f'{program_type}_database_combo', database_combo)
    setattr(parent, f'{program_type}_update_database_list', update_database_list)
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
    
    # WP proteins checkbox (only for protein searches - blastp, blastx)
    if program_type in ["blastp", "blastx"]:
        wp_proteins = QCheckBox("Non-redundant RefSeq proteins (WP)")
        wp_proteins.setStyleSheet("margin-left: 20px;")
        setattr(parent, f'{program_type}_wp_proteins', wp_proteins)
        exclude_layout.addWidget(wp_proteins)
    
    exclude_uncult = QCheckBox("Uncultured/environmental sample sequences")
    exclude_uncult.setStyleSheet("margin-left: 20px;")
    setattr(parent, f'{program_type}_exclude_uncult', exclude_uncult)
    exclude_layout.addWidget(exclude_uncult)
    
    layout.addLayout(exclude_layout)
    
    # Limit to sequences from type material (seqFTLimitBy section)
    type_material_layout = QVBoxLayout()
    type_material_label = QLabel("Limit to")
    type_material_label.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    type_material_hint = QLabel("Optional")
    type_material_hint.setStyleSheet("color: #666; font-size: 11px;")
    type_material_layout.addWidget(type_material_label)
    type_material_layout.addWidget(type_material_hint)
    
    type_material_checkbox = QCheckBox("Sequences from type material")
    type_material_checkbox.setStyleSheet("margin-left: 20px;")
    setattr(parent, f'{program_type}_type_material', type_material_checkbox)
    type_material_layout.addWidget(type_material_checkbox)
    
    layout.addLayout(type_material_layout)
    
    # Note: Genetic Code is now in the query section for NCBI consistency
    
    # Entrez Query
    entrez_layout = QVBoxLayout()
    entrez_label = QLabel("Entrez Query")
    entrez_label.setStyleSheet("font-weight: normal; margin-bottom: 5px;")
    entrez_layout.addWidget(entrez_label)
    
    entrez_input = QLineEdit()
    entrez_input.setFixedWidth(400)
    entrez_input.setStyleSheet("border: 1px solid #ccc; padding: 4px;")
    entrez_input.setPlaceholderText("Enter an Entrez query to limit search (optional)")
    entrez_input.setToolTip("Use Entrez query syntax to search a subset of the selected database")
    setattr(parent, f'{program_type}_entrez_query', entrez_input)
    
    entrez_help_layout = QHBoxLayout()
    entrez_help_layout.addWidget(entrez_input)
    
    entrez_help = QPushButton("?")
    entrez_help.setFixedSize(20, 20)
    entrez_help.setStyleSheet("""
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
    entrez_help.setToolTip("Entrez query help")
    entrez_help_layout.addWidget(entrez_help)
    entrez_help_layout.addStretch()
    
    entrez_layout.addLayout(entrez_help_layout)
    
    entrez_help_text = QLabel("Enter an Entrez query to limit search to a subset of the database")
    entrez_help_text.setStyleSheet("color: #666; font-size: 11px; margin-top: 5px;")
    entrez_layout.addWidget(entrez_help_text)
    
    layout.addLayout(entrez_layout)
    
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
    
    # Algorithm options based on program type
    if program_type == "blastp":
        algorithms = [
            ("kmerBlastp", "Quick BLASTP (Accelerated protein-protein BLAST)"),
            ("blastp", "blastp (protein-protein BLAST)"),
            ("psiBlast", "PSI-BLAST (Position-Specific Iterated BLAST)"),
            ("phiBlast", "PHI-BLAST (Pattern Hit Initiated BLAST)"),
            ("deltaBlast", "DELTA-BLAST (Domain Enhanced Lookup Time Accelerated BLAST)")
        ]
        default_algo = "blastp"
    elif program_type == "blastn":
        algorithms = [
            ("megaBlast", "Highly similar sequences (megablast)"),
            ("discoMegablast", "More dissimilar sequences (discontiguous megablast)"),
            ("blastn", "Somewhat similar sequences (blastn)")
        ]
        default_algo = "megaBlast"
    elif program_type == "blastx":
        algorithms = [
            ("blastx", "blastx (translated nucleotide vs protein database)")
        ]
        default_algo = "blastx"
    elif program_type == "tblastn":
        algorithms = [
            ("tblastn", "tblastn (protein vs translated nucleotide database)")
        ]
        default_algo = "tblastn"
    elif program_type == "tblastx":
        algorithms = [
            ("tblastx", "tblastx (translated nucleotide vs translated nucleotide)")
        ]
        default_algo = "tblastx"
    else:
        algorithms = [(program_type, f"{program_type} search")]
        default_algo = program_type
    
    for i, (algo_id, algo_desc) in enumerate(algorithms):
        radio = QRadioButton(algo_desc)
        radio.setStyleSheet("margin-bottom: 5px;")
        if algo_id == default_algo:
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
    phi_label.setVisible(False)
    phi_layout.addWidget(phi_label)
    
    phi_input = QLineEdit()
    phi_input.setFixedWidth(300)
    phi_input.setStyleSheet("border: 1px solid #ccc; padding: 4px;")
    phi_input.setVisible(False)
    setattr(parent, f'{program_type}_phi_pattern', phi_input)
    setattr(parent, f'{program_type}_phi_label', phi_label)
    phi_layout.addWidget(phi_input)
    phi_layout.addStretch()
    
    layout.addLayout(phi_layout)
    
    # Connect algorithm radio buttons to update UI
    def update_algorithm_ui():
        """Update UI elements based on selected algorithm."""
        checked_button = algo_group.checkedButton()
        if checked_button:
            button_text = checked_button.text()
            
            # Show/hide PHI pattern input
            phi_visible = "PHI-BLAST" in button_text
            phi_input.setVisible(phi_visible)
            phi_label.setVisible(phi_visible)
            
            # Update summary text
            summary_label = getattr(parent, f'{program_type}_summary_label')
            if "Quick BLASTP" in button_text:
                summary_text = f"Search nr using Quick BLASTP (Accelerated protein-protein BLAST)"
            elif "PSI-BLAST" in button_text:
                summary_text = f"Search nr using PSI-BLAST (Position-Specific Iterated BLAST)"
            elif "PHI-BLAST" in button_text:
                summary_text = f"Search nr using PHI-BLAST (Pattern Hit Initiated BLAST)"
            elif "DELTA-BLAST" in button_text:
                summary_text = f"Search nr using DELTA-BLAST (Domain Enhanced Lookup Time Accelerated BLAST)"
            elif "megablast" in button_text.lower():
                if "discontiguous" in button_text.lower():
                    summary_text = f"Search nt using Discontiguous Megablast (More dissimilar sequences)"
                else:
                    summary_text = f"Search nt using Megablast (Optimize for highly similar sequences)"
            elif "blastn" in button_text.lower() and program_type == "blastn":
                summary_text = f"Search nt using Blastn (Somewhat similar sequences)"
            elif program_type == "blastp":
                summary_text = f"Search nr using Blastp (protein-protein BLAST)"
            elif program_type == "blastn":
                summary_text = f"Search nt using Blastn (nucleotide-nucleotide BLAST)"
            elif program_type == "blastx":
                summary_text = f"Search nr using Blastx (translated nucleotide vs protein)"
            elif program_type == "tblastn":
                summary_text = f"Search nt using Tblastn (protein vs translated nucleotide)"
            elif program_type == "tblastx":
                summary_text = f"Search nt using Tblastx (translated nucleotide vs translated nucleotide)"
            else:
                summary_text = f"Search using {program_type.upper()}"
            
            summary_label.setText(summary_text)
    
    # Connect all radio buttons to the update function
    for button in algo_group.buttons():
        button.toggled.connect(update_algorithm_ui)
    
    # Store the update function for external access
    setattr(parent, f'{program_type}_update_algorithm_ui', update_algorithm_ui)
    
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
    # Create summary based on program type
    if program_type == "blastp":
        summary_text = "Search nr using Blastp (protein-protein BLAST)"
    elif program_type == "blastn":
        summary_text = "Search nt using Blastn (nucleotide-nucleotide BLAST)"
    elif program_type == "blastx":
        summary_text = "Search nr using Blastx (translated nucleotide vs protein)"
    elif program_type == "tblastn":
        summary_text = "Search nt using Tblastn (protein vs translated nucleotide)"
    elif program_type == "tblastx":
        summary_text = "Search nt using Tblastx (translated nucleotide vs translated nucleotide)"
    else:
        summary_text = f"Search using {program_type.upper()}"
    
    summary_label = QLabel(summary_text)
    summary_label.setStyleSheet("color: #666; font-size: 12px; margin-top: 10px;")
    setattr(parent, f'{program_type}_summary_label', summary_label)
    summary_layout.addWidget(summary_label)
    
    # Removed "Show results in a new window" checkbox
    
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
        new_visible = not visible
        params_content.setVisible(new_visible)
        toggle_button.setText("▲ Algorithm parameters" if new_visible else "▼ Algorithm parameters")
    
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
    
    # Set default expect threshold based on program type
    if program_type == "tblastn":
        default_expect = "0.05"
    else:
        default_expect = "0.05"
    
    expect_input = QLineEdit(default_expect)
    expect_input.setFixedWidth(100)
    expect_input.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    setattr(parent, f'{program_type}_expect_threshold', expect_input)
    expect_layout.addWidget(expect_input)
    expect_layout.addStretch()
    
    layout.addLayout(expect_layout)
    
    # Word size (program-specific options)
    word_layout = QHBoxLayout()
    word_label = QLabel("Word size")
    word_label.setFixedWidth(150)
    word_layout.addWidget(word_label)
    
    word_combo = QComboBox()
    word_combo.setFixedWidth(100)
    
    # Set word size options based on program type
    if program_type == "blastp":
        word_options = ["2", "3", "5", "6"]
        default_word = "3"
    elif program_type == "blastn":
        word_options = ["7", "11", "15", "20", "24", "28"]
        default_word = "11"
    elif program_type == "tblastn":
        word_options = ["2", "3", "5", "6"]
        default_word = "5"
    elif program_type == "blastx":
        word_options = ["2", "3", "6"]
        default_word = "3"
    elif program_type == "tblastx":
        word_options = ["2", "3"]
        default_word = "3"
    else:
        word_options = ["3"]
        default_word = "3"
    
    word_combo.addItems(word_options)
    word_combo.setCurrentText(default_word)
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
    
    # Matrix (only for protein-based searches)
    if program_type in ["blastp", "blastx", "tblastn"]:
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
    
    # Match/Mismatch scores (for nucleotide searches)
    if program_type in ["blastn", "tblastx"]:
        match_layout = QHBoxLayout()
        match_label = QLabel("Match/Mismatch Scores")
        match_label.setFixedWidth(150)
        match_layout.addWidget(match_label)
        
        match_combo = QComboBox()
        match_combo.setFixedWidth(120)
        match_scores = ["1,-2", "1,-3", "1,-4", "2,-3", "4,-5", "1,-1"]
        match_combo.addItems([f"{score} (Match,Mismatch)" for score in match_scores])
        match_combo.setCurrentText("2,-3 (Match,Mismatch)")
        match_combo.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
        setattr(parent, f'{program_type}_match_scores', match_combo)
        match_layout.addWidget(match_combo)
        match_layout.addStretch()
        
        layout.addLayout(match_layout)
    
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
    
    # Compositional adjustments (only for protein-based searches)
    if program_type in ["blastp", "blastx", "tblastn"]:
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
    
    # Genetic code is now in the search set section, not here
    
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
    
    # View on NCBI button removed
    
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
    # Set appropriate placeholder text based on program type
    if program_type == "blastp":
        placeholder_text = (
            "BLAST results will appear here...\n\n"
            "• Enter a protein sequence above\n"
            "• Configure your search parameters\n"
            "• Click the BLAST button to begin\n"
            "• Results typically take 30 seconds to 2 minutes"
        )
    elif program_type == "blastn":
        placeholder_text = (
            "BLAST results will appear here...\n\n"
            "• Enter a nucleotide sequence above\n"
            "• Configure your search parameters\n"
            "• Click the BLAST button to begin\n"
            "• Results typically take 30 seconds to 2 minutes"
        )
    elif program_type == "blastx":
        placeholder_text = (
            "BLASTX results will appear here...\n\n"
            "• Enter a nucleotide sequence above\n"
            "• Your sequence will be translated in all 6 reading frames\n"
            "• Configure your search parameters\n"
            "• Click the BLAST button to begin\n"
            "• Results typically take 30 seconds to 2 minutes"
        )
    elif program_type == "tblastn":
        placeholder_text = (
            "TBLASTN results will appear here...\n\n"
            "• Enter a protein sequence above\n"
            "• Database sequences will be translated in all 6 reading frames\n"
            "• Configure your search parameters\n"
            "• Click the BLAST button to begin\n"
            "• Results typically take 30 seconds to 2 minutes"
        )
    elif program_type == "tblastx":
        placeholder_text = (
            "TBLASTX results will appear here...\n\n"
            "• Enter a nucleotide sequence above\n"
            "• Both query and database will be translated in all 6 reading frames\n"
            "• Configure your search parameters\n"
            "• Click the BLAST button to begin\n"
            "• Results typically take 30 seconds to 2 minutes"
        )
    else:
        placeholder_text = (
            "BLAST results will appear here...\n\n"
            "• Enter a sequence above\n"
            "• Configure your search parameters\n"
            "• Click the BLAST button to begin\n"
            "• Results typically take 30 seconds to 2 minutes"
        )
    
    results_display.setPlaceholderText(placeholder_text)
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
    
    # Validate sequence based on query type
    if program_type in ["blastp", "tblastn"]:
        sequence_type = "protein"
    else:
        sequence_type = "nucleotide"
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
    
    # Add program-specific parameters
    
    # Matrix (for protein-based searches)
    if hasattr(parent, f'{program_type}_matrix'):
        matrix = getattr(parent, f'{program_type}_matrix')
        parameters['matrix'] = matrix.currentText()
    
    # Match/Mismatch scores (for nucleotide searches)
    if hasattr(parent, f'{program_type}_match_scores'):
        match_scores = getattr(parent, f'{program_type}_match_scores')
        parameters['match_scores'] = match_scores.currentText().split(" ")[0]
    
    # Gap costs
    if hasattr(parent, f'{program_type}_gap_costs'):
        gap_costs = getattr(parent, f'{program_type}_gap_costs')
        gap_text = gap_costs.currentText()
        if "Existence:" in gap_text:
            parts = gap_text.replace("Existence:", "").replace("Extension:", "").split()
            if len(parts) >= 2:
                parameters['gapopen'] = parts[0]
                parameters['gapextend'] = parts[1]
    
    # Compositional adjustment (for protein-based searches)
    if hasattr(parent, f'{program_type}_comp_adjust'):
        comp_adjust = getattr(parent, f'{program_type}_comp_adjust')
        comp_map = {
            "No adjustment": "0",
            "Composition-based statistics": "1", 
            "Conditional compositional score matrix adjustment": "2",
            "Universal compositional score matrix adjustment": "3"
        }
        parameters['comp_adjust'] = comp_map.get(comp_adjust.currentText(), "2")
    
    # Genetic code (for translated searches)
    if hasattr(parent, f'{program_type}_genetic_code'):
        genetic_code = getattr(parent, f'{program_type}_genetic_code')
        parameters['genetic_code'] = genetic_code.currentData()
    
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
        elif "megablast" in button_text.lower():
            parameters['algorithm'] = 'megaBlast'
        elif "discontiguous" in button_text.lower():
            parameters['algorithm'] = 'discoMegablast'
        else:
            parameters['algorithm'] = program_type
    
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
    
    # Get job title and request ID
    job_title_input = getattr(parent, f'{program_type}_job_title', None)
    job_title = job_title_input.text().strip() if job_title_input else None
    
    worker = getattr(parent, f'{program_type}_worker', None)
    request_id = worker.request_id if worker else None
    
    formatted_result = format_blast_output(result, program_type, job_title, request_id)
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
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            file_label = getattr(parent, f'{program_type}_file_label')
            query_input = getattr(parent, f'{program_type}_query_input')
            
            file_label.setText(os.path.basename(file_path))
            query_input.setText(content)
            
        except Exception as e:
            QMessageBox.critical(parent, "File Error", f"Could not read file: {str(e)}")


def open_subject_file(parent, program_type):
    """Open a FASTA file for subject sequence (BL2SEQ mode)."""
    
    file_path, _ = QFileDialog.getOpenFileName(
        parent, "Open Subject Sequence File", "", "FASTA Files (*.fasta *.fa *.fas);;All Files (*)"
    )
    
    if file_path:
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            file_label = getattr(parent, f'{program_type}_subject_file_label')
            subject_input = getattr(parent, f'{program_type}_subject_input')
            
            file_label.setText(os.path.basename(file_path))
            subject_input.setText(content)
            
        except Exception as e:
            QMessageBox.critical(parent, "Subject File Error", f"Could not read subject file: {str(e)}")


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
            with open(file_path, 'w', encoding='utf-8') as f:
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


# ============================================================================
# BLASTN UNIFIED IMPLEMENTATION
# ============================================================================

def create_ncbi_style_blastn_tab_unified(parent):
    """Create BLASTN tab with exact NCBI layout from blastn.html."""
    
    # Main widget
    main_widget = QWidget()
    main_layout = QVBoxLayout(main_widget)
    main_layout.setContentsMargins(10, 10, 10, 10)
    main_layout.setSpacing(0)
    
    # Initialize worker
    setattr(parent, 'blastn_worker', None)
    
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
    query_section = create_blastn_query_section_unified(parent)
    form_layout.addWidget(query_section)
    
    # 2. Choose Search Set Section  
    search_set_section = create_search_set_section(parent, "blastn", "nucleotide")
    form_layout.addWidget(search_set_section)
    
    # 3. Program Selection Section
    program_section = create_program_selection_section(parent, "blastn")
    form_layout.addWidget(program_section)
    
    # 4. Algorithm Parameters Section (collapsible)
    params_section = create_algorithm_parameters_section(parent, "blastn")
    form_layout.addWidget(params_section)
    
    # Add stretch
    form_layout.addStretch()
    
    # Set form widget to scroll area
    scroll_area.setWidget(form_widget)
    main_layout.addWidget(scroll_area)
    
    # Results section
    results_section = create_results_section(parent, "blastn")
    main_layout.addWidget(results_section)
    
    return main_widget


def create_blastn_query_section_unified(parent):
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
    clear_button.clicked.connect(lambda: getattr(parent, 'blastn_query_input').clear())
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
    setattr(parent, 'blastn_query_input', query_input)
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
    setattr(parent, 'blastn_query_from', from_input)
    
    to_label = QLabel("To")
    to_input = QLineEdit()
    to_input.setFixedWidth(80)
    to_input.setStyleSheet("border: 1px solid #ccc; padding: 2px;")
    setattr(parent, 'blastn_query_to', to_input)
    
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
    upload_button.clicked.connect(lambda: open_blast_file(parent, "blastn"))
    upload_layout.addWidget(upload_button)
    
    file_label = QLabel("No file selected")
    file_label.setStyleSheet("color: #666; margin-left: 10px;")
    setattr(parent, 'blastn_file_label', file_label)
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
    setattr(parent, 'blastn_job_title', job_input)
    job_layout.addWidget(job_input)
    job_layout.addStretch()
    
    layout.addLayout(job_layout)
    
    return section