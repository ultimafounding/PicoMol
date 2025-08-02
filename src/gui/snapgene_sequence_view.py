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

class SnapGeneSequenceWidget(QTextEdit):
    """Authentic SnapGene-style sequence display widget"""
    
    positionChanged = pyqtSignal(int)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.record = None
        self.show_features = True
        self.show_translation = True
        self.show_complement = False
        self.show_enzymes = False
        self.bases_per_line = 60  # SnapGene default
        self.current_position = 0
        
        # Setup authentic SnapGene appearance
        self.setFont(QFont("Consolas", 11))  # SnapGene uses Consolas
        self.setReadOnly(True)
        self.setLineWrapMode(QTextEdit.NoWrap)
        self.setStyleSheet("""
            QTextEdit {
                background-color: white;
                border: none;
                selection-background-color: #316AC5;
                color: black;
                padding: 10px;
                line-height: 1.2;
            }
        """)
        
        # Authentic SnapGene feature colors
        self.feature_colors = {
            "gene": QColor(255, 255, 0, 100),        # Yellow
            "CDS": QColor(255, 165, 0, 100),         # Orange  
            "promoter": QColor(255, 0, 0, 100),      # Red
            "terminator": QColor(128, 0, 128, 100),  # Purple
            "rep_origin": QColor(0, 255, 0, 100),    # Green
            "misc_feature": QColor(192, 192, 192, 100), # Gray
            "regulatory": QColor(255, 192, 203, 100), # Pink
            "enhancer": QColor(255, 255, 224, 100),   # Light yellow
            "primer_bind": QColor(0, 255, 255, 100),  # Cyan
            "protein_bind": QColor(144, 238, 144, 100) # Light green
        }
        
        self.cursorPositionChanged.connect(self.on_cursor_changed)
    
    def display_sequence(self, record):
        """Display sequence in authentic SnapGene style"""
        self.record = record
        self.update_display()
    
    def update_display(self):
        """Update the sequence display with authentic SnapGene formatting"""
        if not self.record:
            self.clear()
            return
        
        self.clear()
        
        sequence = str(self.record.seq).upper()
        display_lines = []
        
        # Add header with sequence info
        display_lines.append(f"Sequence: {self.record.id} ({len(sequence)} bp)")
        display_lines.append("")
        
        # Process sequence in chunks
        for i in range(0, len(sequence), self.bases_per_line):
            line_start = i + 1
            line_end = min(i + self.bases_per_line, len(sequence))
            line_seq = sequence[i:line_end]
            
            # Add feature annotation line if features are shown
            if self.show_features and hasattr(self.record, 'features'):
                feature_line = self.create_feature_line(i, line_seq)
                if feature_line.strip():
                    display_lines.append(f"       {feature_line}")
            
            # Add position ruler
            ruler = self.create_ruler(len(line_seq))
            display_lines.append(f"       {ruler}")
            
            # Add main sequence line with position number
            formatted_seq = self.format_sequence(line_seq)
            display_lines.append(f"{line_start:>6} {formatted_seq}")
            
            # Add translation if enabled
            if self.show_translation:
                translation = self.create_translation(line_seq, line_start)
                if translation.strip():
                    display_lines.append(f"       {translation}")
            
            # Add complement if enabled
            if self.show_complement:
                complement = self.create_complement(line_seq)
                display_lines.append(f"       {complement}")
            
            # Add spacing between blocks
            display_lines.append("")
        
        # Set the formatted text
        self.setPlainText("\n".join(display_lines))
        
        # Apply feature highlighting
        if self.show_features:
            self.apply_feature_highlighting()
    
    def create_ruler(self, length):
        """Create position ruler like SnapGene"""
        ruler = ""
        for i in range(length):
            if i % 10 == 0:
                ruler += "|"
            elif i % 5 == 0:
                ruler += "."
            else:
                ruler += " "
        return ruler
    
    def create_feature_line(self, start_pos, sequence):
        """Create feature annotation line"""
        feature_chars = [" "] * len(sequence)
        
        if hasattr(self.record, 'features'):
            for feature in self.record.features:
                if feature.type == 'source':
                    continue
                    
                feat_start = int(feature.location.start)
                feat_end = int(feature.location.end)
                
                # Check overlap with current line
                line_end = start_pos + len(sequence)
                if feat_start < line_end and feat_end > start_pos:
                    # Calculate overlap positions
                    overlap_start = max(feat_start - start_pos, 0)
                    overlap_end = min(feat_end - start_pos, len(sequence))
                    
                    # Get feature character
                    char = self.get_feature_char(feature.type)
                    
                    # Mark the feature
                    for pos in range(overlap_start, overlap_end):
                        if 0 <= pos < len(feature_chars):
                            feature_chars[pos] = char
        
        return self.format_sequence("".join(feature_chars))
    
    def get_feature_char(self, feature_type):
        """Get character for feature type"""
        chars = {
            "gene": "=",
            "CDS": "-", 
            "promoter": ">",
            "terminator": "|",
            "rep_origin": "o",
            "misc_feature": ".",
            "regulatory": "~",
            "enhancer": "+",
            "primer_bind": "^",
            "protein_bind": "#"
        }
        return chars.get(feature_type, ".")
    
    def format_sequence(self, sequence):
        """Format sequence with SnapGene-style spacing"""
        formatted = ""
        for i, char in enumerate(sequence):
            if i > 0 and i % 10 == 0:
                formatted += " "
            formatted += char
        return formatted
    
    def create_translation(self, sequence, start_pos):
        """Create translation line"""
        # Calculate reading frame
        frame = (start_pos - 1) % 3
        
        # Pad sequence to align with reading frame
        padded_seq = " " * frame + sequence
        
        translation = ""
        for i in range(0, len(padded_seq), 3):
            codon = padded_seq[i:i+3]
            if len(codon) == 3 and all(c in 'ATCG ' for c in codon):
                if ' ' not in codon:
                    try:
                        aa = str(Seq(codon).translate())
                        translation += aa + "  "
                    except:
                        translation += "   "
                else:
                    translation += "   "
            else:
                translation += " " * (3 - len(codon)) if len(codon) < 3 else ""
        
        # Format with spacing
        return self.format_sequence(translation.rstrip())
    
    def create_complement(self, sequence):
        """Create complement line"""
        comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        complement = "".join(comp_map.get(base, base) for base in sequence[::-1])
        return self.format_sequence(complement)
    
    def apply_feature_highlighting(self):
        """Apply feature highlighting"""
        if not self.record or not hasattr(self.record, 'features'):
            return
        
        # This is a simplified version - full implementation would be more complex
        pass
    
    def on_cursor_changed(self):
        """Handle cursor position changes"""
        cursor = self.textCursor()
        position = cursor.position()
        
        # Convert to sequence position (simplified)
        seq_pos = max(0, position // 10)  # Rough estimate
        if seq_pos != self.current_position:
            self.current_position = seq_pos
            self.positionChanged.emit(seq_pos)
    
    def set_show_features(self, show):
        """Toggle feature display"""
        self.show_features = show
        self.update_display()
    
    def set_show_translation(self, show):
        """Toggle translation display"""
        self.show_translation = show
        self.update_display()
    
    def set_show_complement(self, show):
        """Toggle complement display"""
        self.show_complement = show
        self.update_display()
    
    def set_show_enzymes(self, show):
        """Toggle enzyme display"""
        self.show_enzymes = show
        self.update_display()
    
    def set_bases_per_line(self, bases):
        """Set bases per line"""
        self.bases_per_line = bases
        self.update_display()

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
        
        # Left side - sequence display
        self.sequence_widget = SnapGeneSequenceWidget()
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
        self.sequence_widget.set_show_enzymes(checked)
    
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