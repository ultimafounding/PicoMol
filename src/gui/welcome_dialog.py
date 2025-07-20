#!/usr/bin/env python3
"""
Enhanced Welcome Dialog for PicoMol.

This module provides a comprehensive welcome screen with feature overview,
quick start guide, and modern design.
"""

import webbrowser
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QCheckBox, QPushButton,
    QScrollArea, QWidget, QGroupBox
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont


class WelcomeDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Welcome to PicoMol!")
        self.setMinimumSize(700, 600)
        self.setMaximumSize(800, 700)
        
        # Create scroll area for content
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        
        # Main content widget
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)
        content_layout.setContentsMargins(30, 20, 30, 20)
        content_layout.setSpacing(15)
        
        # Header section with logo and title
        header_widget = QWidget()
        header_layout = QVBoxLayout(header_widget)
        header_layout.setAlignment(Qt.AlignCenter)
        
        # Title
        title_label = QLabel("🧬 PicoMol")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("""
            QLabel {
                font-size: 32px;
                font-weight: bold;
                color: #2c3e50;
                margin: 10px 0;
            }
        """)
        header_layout.addWidget(title_label)
        
        # Subtitle
        subtitle_label = QLabel("Molecular Visualization & Bioinformatics Suite")
        subtitle_label.setAlignment(Qt.AlignCenter)
        subtitle_label.setStyleSheet("""
            QLabel {
                font-size: 16px;
                color: #7f8c8d;
                font-style: italic;
                margin-bottom: 20px;
            }
        """)
        header_layout.addWidget(subtitle_label)
        
        content_layout.addWidget(header_widget)
        
        # Description
        desc_label = QLabel(
            "<p style='font-size: 14px; line-height: 1.6; color: #34495e; text-align: center;'>"
            "A comprehensive desktop application for molecular visualization, structural analysis, "
            "and bioinformatics research. Built for researchers, students, and developers working "
            "with protein structures and biomolecular data."
            "</p>"
        )
        desc_label.setWordWrap(True)
        desc_label.setAlignment(Qt.AlignCenter)
        content_layout.addWidget(desc_label)
        
        # Features section
        features_group = QGroupBox("🔬 Key Features")
        features_group.setStyleSheet("""
            QGroupBox {
                font-size: 16px;
                font-weight: bold;
                color: #2c3e50;
                border: 2px solid #bdc3c7;
                border-radius: 8px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 10px 0 10px;
                background-color: white;
            }
        """)
        features_layout = QVBoxLayout(features_group)
        
        # Create feature items in a grid-like layout
        features_data = [
            ("🧪", "3D Molecular Visualization", "Interactive protein structure viewing with multiple representations, themes, and export options"),
            ("🔍", "BLAST Integration", "Complete BLAST suite (BLASTP, BLASTN, BLASTX, TBLASTN, TBLASTX) with direct NCBI connectivity"),
            ("🧬", "Bioinformatics Analysis", "Comprehensive sequence analysis for proteins, DNA, and RNA with detailed statistics"),
            ("🎯", "Motif & Domain Analysis", "InterPro and PROSITE integration for protein domain annotation and motif identification"),
            ("📄", "FASTA Support", "Smart parsing of single and multi-sequence FASTA files with automatic type detection"),
            ("🎨", "Customizable Interface", "Multiple built-in themes, preferences system, and responsive design"),
            ("💾", "Data Management", "PDB fetching, local file support, drag-and-drop, and export capabilities"),
            ("🌐", "API Integration", "Direct connections to NCBI, EBI InterPro, and ExPASy PROSITE web services")
        ]
        
        for i in range(0, len(features_data), 2):
            row_widget = QWidget()
            row_layout = QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            
            # Add first feature
            feature1 = self.create_feature_widget(*features_data[i])
            row_layout.addWidget(feature1)
            
            # Add second feature if it exists
            if i + 1 < len(features_data):
                feature2 = self.create_feature_widget(*features_data[i + 1])
                row_layout.addWidget(feature2)
            else:
                row_layout.addStretch()
            
            features_layout.addWidget(row_widget)
        
        content_layout.addWidget(features_group)
        
        # Quick start section
        quickstart_group = QGroupBox("🚀 Quick Start Guide")
        quickstart_group.setStyleSheet(features_group.styleSheet())
        quickstart_layout = QVBoxLayout(quickstart_group)
        
        quickstart_steps = [
            "1. **Load a Structure:** Enter a PDB ID (e.g., '1CRN') or drag & drop a local PDB file",
            "2. **Explore Visualization:** Use the 3D Viewer tab to adjust representations and colors",
            "3. **Analyze Sequences:** Switch to the Bioinformatics tab for comprehensive sequence analysis",
            "4. **Find Motifs & Domains:** Use the Motifs & Domains tab for protein functional annotation",
            "5. **Run BLAST Searches:** Use the BLAST tab to search against NCBI databases",
            "6. **Customize Experience:** Access Preferences (Tools → Preferences) to personalize your workflow"
        ]
        
        for step in quickstart_steps:
            step_label = QLabel(step)
            step_label.setWordWrap(True)
            step_label.setStyleSheet("""
                QLabel {
                    font-size: 13px;
                    color: #2c3e50;
                    padding: 5px 10px;
                    margin: 2px 0;
                    background-color: #ecf0f1;
                    border-radius: 4px;
                }
            """)
            quickstart_layout.addWidget(step_label)
        
        content_layout.addWidget(quickstart_group)
        
        # Tips section
        tips_label = QLabel(
            "<div style='background-color: #e8f5e8; border: 1px solid #27ae60; border-radius: 6px; padding: 15px; margin: 10px 0;'>"
            "<p style='margin: 0; font-size: 13px; color: #27ae60;'>"
            "<b>💡 Pro Tips:</b><br>"
            "• Hover over any control for helpful tooltips<br>"
            "• Use Ctrl+Z/Ctrl+Y for undo/redo in the 3D viewer<br>"
            "• Access recent files through File → Recent Files<br>"
            "• Try 'Search All Databases' for comprehensive motif analysis<br>"
            "• Export analysis results in multiple formats<br>"
            "• Try different themes in Tools → Preferences"
            "</p>"
            "</div>"
        )
        tips_label.setWordWrap(True)
        content_layout.addWidget(tips_label)
        
        # Set scroll widget
        scroll.setWidget(content_widget)
        
        # Main dialog layout
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(scroll)
        
        # Bottom section with checkbox and buttons
        bottom_widget = QWidget()
        bottom_layout = QVBoxLayout(bottom_widget)
        bottom_layout.setContentsMargins(20, 10, 20, 20)
        
        # Checkbox
        self.checkbox = QCheckBox("Show this welcome screen on startup")
        self.checkbox.setChecked(True)
        self.checkbox.setStyleSheet("""
            QCheckBox {
                font-size: 13px;
                color: #2c3e50;
                spacing: 8px;
            }
            QCheckBox::indicator {
                width: 16px;
                height: 16px;
            }
        """)
        bottom_layout.addWidget(self.checkbox)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        # Help button
        help_button = QPushButton("📖 Documentation")
        help_button.setToolTip("Open online documentation and help resources")
        help_button.clicked.connect(self.open_documentation)
        help_button.setStyleSheet("""
            QPushButton {
                background-color: #3498db;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                font-size: 13px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #2980b9;
            }
        """)
        
        # Get started button
        get_started_button = QPushButton("🚀 Get Started")
        get_started_button.setToolTip("Close this dialog and start using PicoMol")
        get_started_button.clicked.connect(self.accept)
        get_started_button.setDefault(True)
        get_started_button.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                border: none;
                padding: 8px 20px;
                border-radius: 4px;
                font-size: 13px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #229954;
            }
        """)
        
        button_layout.addWidget(help_button)
        button_layout.addStretch()
        button_layout.addWidget(get_started_button)
        
        bottom_layout.addLayout(button_layout)
        main_layout.addWidget(bottom_widget)
    
    def create_feature_widget(self, icon, title, description):
        """Create a feature widget with icon, title, and description."""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        
        # Icon and title
        header_layout = QHBoxLayout()
        header_layout.setContentsMargins(0, 0, 0, 0)
        
        icon_label = QLabel(icon)
        icon_label.setStyleSheet("font-size: 24px;")
        icon_label.setFixedSize(32, 32)
        icon_label.setAlignment(Qt.AlignCenter)
        header_layout.addWidget(icon_label)
        
        title_label = QLabel(title)
        title_label.setStyleSheet("""
            QLabel {
                font-size: 14px;
                font-weight: bold;
                color: #2c3e50;
            }
        """)
        header_layout.addWidget(title_label)
        header_layout.addStretch()
        
        layout.addLayout(header_layout)
        
        # Description
        desc_label = QLabel(description)
        desc_label.setWordWrap(True)
        desc_label.setStyleSheet("""
            QLabel {
                font-size: 12px;
                color: #7f8c8d;
                line-height: 1.4;
            }
        """)
        layout.addWidget(desc_label)
        
        # Style the widget
        widget.setStyleSheet("""
            QWidget {
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 6px;
                margin: 2px;
            }
            QWidget:hover {
                background-color: #e9ecef;
                border-color: #3498db;
            }
        """)
        
        widget.setMinimumHeight(100)
        widget.setMaximumHeight(120)
        
        return widget
    
    def open_documentation(self):
        """Open documentation in web browser."""
        webbrowser.open_new_tab("https://github.com/ultimafounding/PicoMol/blob/main/README.md")
    
    def should_show_next_time(self):
        return self.checkbox.isChecked()