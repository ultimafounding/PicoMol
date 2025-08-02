from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QLabel, 
                             QScrollArea, QFrame, QSplitter, QToolBar, QAction,
                             QComboBox, QSpinBox, QCheckBox, QPushButton, QMenu,
                             QTableWidget, QTableWidgetItem, QHeaderView, QTabWidget,
                             QGroupBox, QFormLayout)
from PyQt5.QtGui import (QFont, QColor, QTextCharFormat, QTextCursor, QPainter, 
                         QPen, QBrush, QFontMetrics, QTextDocument, QTextBlock)
from PyQt5.QtCore import Qt, pyqtSignal, QRect, QPoint
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import math

# Import custom sequence widget
from src.gui.custom_sequence_widget import CustomSequenceWidget

# Try to import Bio.Restriction for enzyme analysis
try:
    from Bio.Restriction import RestrictionBatch, Analysis
    from Bio.Restriction import EcoRI, BamHI, HindIII, XhoI, SalI, XbaI, SpeI, NotI
    from Bio.Restriction import PstI, SacI, KpnI, SmaI, BglII, NcoI, NdeI, ApaI
    RESTRICTION_AVAILABLE = True
except ImportError:
    RESTRICTION_AVAILABLE = False

# Old SnapGeneSequenceWidget class removed - now using CustomSequenceWidget

class SnapGeneSequenceView(QWidget):
    """Main SnapGene-style sequence view widget"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.record = None
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the SnapGene-style interface"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)
        
        # Create toolbar
        toolbar_layout = QHBoxLayout()
        
        # Show options
        toolbar_layout.addWidget(QLabel("Show:"))
        
        self.features_check = QCheckBox("Features")
        self.features_check.setChecked(True)
        self.features_check.toggled.connect(self.on_features_toggled)
        toolbar_layout.addWidget(self.features_check)
        
        self.translation_check = QCheckBox("Translation")
        self.translation_check.setChecked(True)
        self.translation_check.toggled.connect(self.on_translation_toggled)
        toolbar_layout.addWidget(self.translation_check)
        
        self.complement_check = QCheckBox("Complement")
        self.complement_check.setChecked(False)
        self.complement_check.toggled.connect(self.on_complement_toggled)
        toolbar_layout.addWidget(self.complement_check)
        
        self.enzymes_check = QCheckBox("Enzymes")
        self.enzymes_check.setChecked(False)
        self.enzymes_check.toggled.connect(self.on_enzymes_toggled)
        toolbar_layout.addWidget(self.enzymes_check)
        
        toolbar_layout.addWidget(QLabel("  |  Bases per line:"))
        
        self.bases_spin = QSpinBox()
        self.bases_spin.setRange(20, 100)
        self.bases_spin.setValue(60)
        self.bases_spin.valueChanged.connect(self.on_bases_changed)
        toolbar_layout.addWidget(self.bases_spin)
        
        toolbar_layout.addStretch()
        
        # Add toolbar to layout
        toolbar_widget = QWidget()
        toolbar_widget.setLayout(toolbar_layout)
        toolbar_widget.setStyleSheet("""
            QWidget {
                background-color: #f0f0f0;
                border-bottom: 1px solid #cccccc;
                padding: 5px;
            }
        """)
        layout.addWidget(toolbar_widget)
        
        # Create main content area
        content_layout = QHBoxLayout()
        
        # Left side - sequence display using custom painted widget
        self.sequence_widget = CustomSequenceWidget()
        self.sequence_widget.positionChanged.connect(self.on_position_changed)
        
        # Add scroll area for sequence
        scroll_area = QScrollArea()
        scroll_area.setWidget(self.sequence_widget)
        scroll_area.setWidgetResizable(True)
        scroll_area.setMinimumWidth(600)
        
        content_layout.addWidget(scroll_area, 2)
        
        # Right side - info panel
        info_panel = QWidget()
        info_panel.setMaximumWidth(250)
        info_panel.setStyleSheet("""
            QWidget {
                background-color: #f8f8f8;
                border-left: 1px solid #cccccc;
            }
        """)
        
        info_layout = QVBoxLayout(info_panel)
        
        # Sequence info
        info_group = QGroupBox("Sequence Information")
        info_form = QFormLayout(info_group)
        
        self.name_label = QLabel("-")
        self.length_label = QLabel("-")
        self.type_label = QLabel("-")
        self.gc_label = QLabel("-")
        
        info_form.addRow("Name:", self.name_label)
        info_form.addRow("Length:", self.length_label)
        info_form.addRow("Type:", self.type_label)
        info_form.addRow("GC Content:", self.gc_label)
        
        info_layout.addWidget(info_group)
        
        # Position info
        pos_group = QGroupBox("Current Position")
        pos_form = QFormLayout(pos_group)
        
        self.position_label = QLabel("0")
        self.base_label = QLabel("-")
        
        pos_form.addRow("Position:", self.position_label)
        pos_form.addRow("Base:", self.base_label)
        
        info_layout.addWidget(pos_group)
        
        # Enzyme info (only show if enzymes are enabled)
        self.enzyme_group = QGroupBox("Restriction Enzymes")
        self.enzyme_list = QLabel("Enable enzymes to see cut sites")
        self.enzyme_list.setWordWrap(True)
        self.enzyme_list.setStyleSheet("font-size: 10px;")
        
        enzyme_layout = QVBoxLayout(self.enzyme_group)
        enzyme_layout.addWidget(self.enzyme_list)
        
        info_layout.addWidget(self.enzyme_group)
        info_layout.addStretch()
        
        content_layout.addWidget(info_panel)
        
        # Add content to main layout
        content_widget = QWidget()
        content_widget.setLayout(content_layout)
        layout.addWidget(content_widget)
        
        # Status bar
        self.status_label = QLabel("Ready")
        self.status_label.setStyleSheet("""
            QLabel {
                background-color: #e0e0e0;
                border-top: 1px solid #cccccc;
                padding: 3px 8px;
                font-size: 11px;
            }
        """)
        layout.addWidget(self.status_label)
    
    def display_sequence(self, record):
        """Display a sequence record"""
        self.record = record
        
        # Update sequence display
        self.sequence_widget.display_sequence(record)
        
        # Update enzyme analysis if enzymes are enabled
        if self.enzymes_check.isChecked():
            self.analyze_enzymes()
            self.sequence_widget.enzyme_sites = self.enzyme_sites
            self.update_enzyme_info()
        
        # Update info panel
        if record:
            self.name_label.setText(record.id or "Unnamed")
            self.length_label.setText(f"{len(record.seq)} bp")
            
            # Determine sequence type
            seq_type = "DNA"
            if 'U' in str(record.seq).upper():
                seq_type = "RNA"
            self.type_label.setText(seq_type)
            
            # Calculate GC content
            try:
                gc_content = gc_fraction(record.seq) * 100
                self.gc_label.setText(f"{gc_content:.1f}%")
            except:
                self.gc_label.setText("N/A")
            
            self.status_label.setText(f"Loaded: {record.id} ({len(record.seq)} bp)")
        else:
            self.name_label.setText("-")
            self.length_label.setText("-")
            self.type_label.setText("-")
            self.gc_label.setText("-")
            self.status_label.setText("Ready")
    
    def on_features_toggled(self, checked):
        """Handle features toggle"""
        self.sequence_widget.set_show_features(checked)
    
    def on_translation_toggled(self, checked):
        """Handle translation toggle"""
        self.sequence_widget.set_show_translation(checked)
    
    def on_complement_toggled(self, checked):
        """Handle complement toggle"""
        self.sequence_widget.set_show_complement(checked)
    
    def on_enzymes_toggled(self, checked):
        """Handle enzymes toggle"""
        if checked and self.record:
            self.analyze_enzymes()
            self.sequence_widget.enzyme_sites = self.enzyme_sites
        self.sequence_widget.set_show_enzymes(checked)
        self.update_enzyme_info()
    
    def on_bases_changed(self, bases):
        """Handle bases per line change"""
        self.sequence_widget.set_bases_per_line(bases)
    
    def on_position_changed(self, position):
        """Handle position changes"""
        self.position_label.setText(str(position + 1))
        
        if self.record and position < len(self.record.seq):
            base = str(self.record.seq)[position]
            self.base_label.setText(base)
            self.status_label.setText(f"Position {position + 1} | Base: {base}")
        else:
            self.base_label.setText("-")
    
    def analyze_enzymes(self):
        """Analyze restriction enzyme cut sites"""
        if not RESTRICTION_AVAILABLE or not self.record:
            self.enzyme_sites = {}
            return
        
        try:
            # Create a common set of restriction enzymes
            common_enzymes = RestrictionBatch([
                EcoRI, BamHI, HindIII, XhoI, SalI, XbaI, SpeI, NotI,
                PstI, SacI, KpnI, SmaI, BglII, NcoI, NdeI, ApaI
            ])
            
            # Store cut sites - get individual enzyme results
            self.enzyme_sites = {}
            for enzyme in common_enzymes:
                sites = enzyme.search(self.record.seq)
                if sites:
                    self.enzyme_sites[str(enzyme)] = sites
                    
        except Exception as e:
            print(f"Error analyzing enzymes: {e}")
            self.enzyme_sites = {}
    
    def update_enzyme_info(self):
        """Update enzyme information display"""
        if not hasattr(self, 'enzyme_sites') or not self.enzyme_sites:
            self.enzyme_list.setText("No enzyme cut sites found" if self.enzymes_check.isChecked() else "Enable enzymes to see cut sites")
            return
        
        enzyme_info = []
        for enzyme_name, sites in self.enzyme_sites.items():
            if sites:
                site_list = ", ".join(str(site) for site in sorted(sites))
                enzyme_info.append(f"{enzyme_name}: {site_list}")
        
        if enzyme_info:
            self.enzyme_list.setText("\n".join(enzyme_info[:10]))  # Show first 10 enzymes
        else:
            self.enzyme_list.setText("No enzyme cut sites found")