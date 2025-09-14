from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QTextEdit, QLabel, 
                             QPushButton, QComboBox, QSpinBox, QCheckBox, QGroupBox, 
                             QFormLayout, QTabWidget, QWidget, QTableWidget, QTableWidgetItem,
                             QHeaderView, QMessageBox, QLineEdit, QColorDialog, QListWidget,
                             QListWidgetItem, QSplitter, QToolBar, QAction, QMenu, QInputDialog)
from PyQt5.QtGui import QFont, QColor, QTextCharFormat, QTextCursor, QIcon
from PyQt5.QtCore import Qt, pyqtSignal
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import copy

class FeatureEditor(QDialog):
    """Dialog for editing sequence features"""
    
    def __init__(self, feature=None, sequence_length=0, parent=None):
        super().__init__(parent)
        self.feature = feature
        self.sequence_length = sequence_length
        self.setWindowTitle("Edit Feature" if feature else "Add Feature")
        self.setMinimumSize(500, 400)
        self.setup_ui()
        
        if feature:
            self.load_feature_data()
    
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Basic information
        basic_group = QGroupBox("Basic Information")
        basic_layout = QFormLayout(basic_group)
        
        self.type_combo = QComboBox()
        self.type_combo.addItems([
            "gene", "CDS", "promoter", "terminator", "rep_origin", 
            "misc_feature", "regulatory", "enhancer", "primer_bind", "protein_bind"
        ])
        basic_layout.addRow("Type:", self.type_combo)
        
        self.label_edit = QLineEdit()
        basic_layout.addRow("Label:", self.label_edit)
        
        self.gene_edit = QLineEdit()
        basic_layout.addRow("Gene:", self.gene_edit)
        
        self.product_edit = QLineEdit()
        basic_layout.addRow("Product:", self.product_edit)
        
        layout.addWidget(basic_group)
        
        # Location
        location_group = QGroupBox("Location")
        location_layout = QFormLayout(location_group)
        
        self.start_spin = QSpinBox()
        self.start_spin.setRange(1, self.sequence_length or 10000)
        location_layout.addRow("Start:", self.start_spin)
        
        self.end_spin = QSpinBox()
        self.end_spin.setRange(1, self.sequence_length or 10000)
        location_layout.addRow("End:", self.end_spin)
        
        self.strand_combo = QComboBox()
        self.strand_combo.addItems(["Forward (+1)", "Reverse (-1)", "Both (0)"])
        location_layout.addRow("Strand:", self.strand_combo)
        
        layout.addWidget(location_group)
        
        # Additional qualifiers
        qualifiers_group = QGroupBox("Additional Qualifiers")
        qualifiers_layout = QVBoxLayout(qualifiers_group)
        
        self.qualifiers_table = QTableWidget(0, 2)
        self.qualifiers_table.setHorizontalHeaderLabels(["Qualifier", "Value"])
        self.qualifiers_table.horizontalHeader().setStretchLastSection(True)
        qualifiers_layout.addWidget(self.qualifiers_table)
        
        # Buttons for qualifiers
        qual_buttons = QHBoxLayout()
        add_qual_btn = QPushButton("Add Qualifier")
        add_qual_btn.clicked.connect(self.add_qualifier)
        qual_buttons.addWidget(add_qual_btn)
        
        remove_qual_btn = QPushButton("Remove Selected")
        remove_qual_btn.clicked.connect(self.remove_qualifier)
        qual_buttons.addWidget(remove_qual_btn)
        
        qual_buttons.addStretch()
        qualifiers_layout.addLayout(qual_buttons)
        
        layout.addWidget(qualifiers_group)
        
        # Dialog buttons
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        button_layout.addWidget(ok_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
    
    def load_feature_data(self):
        """Load existing feature data into the form"""
        if not self.feature:
            return
        
        # Basic info
        self.type_combo.setCurrentText(self.feature.type)
        
        qualifiers = self.feature.qualifiers
        if 'label' in qualifiers:
            self.label_edit.setText(qualifiers['label'][0])
        if 'gene' in qualifiers:
            self.gene_edit.setText(qualifiers['gene'][0])
        if 'product' in qualifiers:
            self.product_edit.setText(qualifiers['product'][0])
        
        # Location
        self.start_spin.setValue(int(self.feature.location.start) + 1)  # Convert to 1-based
        self.end_spin.setValue(int(self.feature.location.end))
        
        strand = self.feature.location.strand
        if strand == 1:
            self.strand_combo.setCurrentIndex(0)
        elif strand == -1:
            self.strand_combo.setCurrentIndex(1)
        else:
            self.strand_combo.setCurrentIndex(2)
        
        # Additional qualifiers
        for key, values in qualifiers.items():
            if key not in ['label', 'gene', 'product']:
                self.add_qualifier(key, values[0] if values else "")
    
    def add_qualifier(self, key="", value=""):
        """Add a qualifier row"""
        row = self.qualifiers_table.rowCount()
        self.qualifiers_table.insertRow(row)
        self.qualifiers_table.setItem(row, 0, QTableWidgetItem(key))
        self.qualifiers_table.setItem(row, 1, QTableWidgetItem(value))
    
    def remove_qualifier(self):
        """Remove selected qualifier"""
        current_row = self.qualifiers_table.currentRow()
        if current_row >= 0:
            self.qualifiers_table.removeRow(current_row)
    
    def get_feature(self):
        """Get the feature from the form data"""
        # Create location
        start = self.start_spin.value() - 1  # Convert to 0-based
        end = self.end_spin.value()
        
        strand_text = self.strand_combo.currentText()
        if "Forward" in strand_text:
            strand = 1
        elif "Reverse" in strand_text:
            strand = -1
        else:
            strand = 0
        
        location = FeatureLocation(start, end, strand)
        
        # Create qualifiers
        qualifiers = {}
        
        if self.label_edit.text():
            qualifiers['label'] = [self.label_edit.text()]
        if self.gene_edit.text():
            qualifiers['gene'] = [self.gene_edit.text()]
        if self.product_edit.text():
            qualifiers['product'] = [self.product_edit.text()]
        
        # Add additional qualifiers
        for row in range(self.qualifiers_table.rowCount()):
            key_item = self.qualifiers_table.item(row, 0)
            value_item = self.qualifiers_table.item(row, 1)
            if key_item and value_item and key_item.text():
                qualifiers[key_item.text()] = [value_item.text()]
        
        return SeqFeature(location, type=self.type_combo.currentText(), qualifiers=qualifiers)

class SequenceEditor(QDialog):
    """Advanced sequence editor with feature editing capabilities"""
    
    sequenceModified = pyqtSignal(object)  # Emits modified SeqRecord
    
    def __init__(self, record, parent=None):
        super().__init__(parent)
        self.original_record = record
        self.record = copy.deepcopy(record) if record else None
        self.setWindowTitle("Sequence Editor")
        self.setMinimumSize(900, 700)
        self.setup_ui()
        
        if self.record:
            self.load_sequence()
    
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Toolbar
        toolbar = QToolBar()
        layout.addWidget(toolbar)
        
        # File operations
        save_action = QAction("💾 Save", self)
        save_action.setToolTip("Save changes")
        save_action.triggered.connect(self.save_changes)
        toolbar.addAction(save_action)
        
        revert_action = QAction("↶ Revert", self)
        revert_action.setToolTip("Revert to original")
        revert_action.triggered.connect(self.revert_changes)
        toolbar.addAction(revert_action)
        
        toolbar.addSeparator()
        
        # Sequence operations
        find_action = QAction("🔍 Find", self)
        find_action.setToolTip("Find sequence")
        find_action.triggered.connect(self.find_sequence)
        toolbar.addAction(find_action)
        
        replace_action = QAction("🔄 Replace", self)
        replace_action.setToolTip("Find and replace")
        replace_action.triggered.connect(self.find_replace)
        toolbar.addAction(replace_action)
        
        toolbar.addSeparator()
        
        # Feature operations
        add_feature_action = QAction("➕ Add Feature", self)
        add_feature_action.setToolTip("Add new feature")
        add_feature_action.triggered.connect(self.add_feature)
        toolbar.addAction(add_feature_action)
        
        edit_feature_action = QAction("✏️ Edit Feature", self)
        edit_feature_action.setToolTip("Edit selected feature")
        edit_feature_action.triggered.connect(self.edit_selected_feature)
        toolbar.addAction(edit_feature_action)
        
        delete_feature_action = QAction("🗑️ Delete Feature", self)
        delete_feature_action.setToolTip("Delete selected feature")
        delete_feature_action.triggered.connect(self.delete_selected_feature)
        toolbar.addAction(delete_feature_action)
        
        # Main content
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)
        
        # Left panel - sequence editing
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # Sequence info
        info_group = QGroupBox("Sequence Information")
        info_layout = QFormLayout(info_group)
        
        self.name_edit = QLineEdit()
        info_layout.addRow("Name:", self.name_edit)
        
        self.description_edit = QLineEdit()
        info_layout.addRow("Description:", self.description_edit)
        
        self.length_label = QLabel("0")
        info_layout.addRow("Length:", self.length_label)
        
        left_layout.addWidget(info_group)
        
        # Sequence editor
        seq_group = QGroupBox("Sequence")
        seq_layout = QVBoxLayout(seq_group)
        
        self.sequence_edit = QTextEdit()
        self.sequence_edit.setFont(QFont("Consolas", 11))
        self.sequence_edit.textChanged.connect(self.on_sequence_changed)
        seq_layout.addWidget(self.sequence_edit)
        
        # Sequence tools
        seq_tools = QHBoxLayout()
        
        self.case_combo = QComboBox()
        self.case_combo.addItems(["Keep Case", "Upper Case", "Lower Case"])
        self.case_combo.currentTextChanged.connect(self.change_case)
        seq_tools.addWidget(QLabel("Case:"))
        seq_tools.addWidget(self.case_combo)
        
        complement_btn = QPushButton("Reverse Complement")
        complement_btn.clicked.connect(self.reverse_complement)
        seq_tools.addWidget(complement_btn)
        
        seq_tools.addStretch()
        seq_layout.addLayout(seq_tools)
        
        left_layout.addWidget(seq_group)
        
        splitter.addWidget(left_panel)
        
        # Right panel - features
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        features_group = QGroupBox("Features")
        features_layout = QVBoxLayout(features_group)
        
        self.features_list = QListWidget()
        self.features_list.itemDoubleClicked.connect(self.edit_selected_feature)
        self.features_list.itemSelectionChanged.connect(self.on_feature_selection_changed)
        features_layout.addWidget(self.features_list)
        
        # Feature buttons
        feature_buttons = QHBoxLayout()
        
        add_btn = QPushButton("Add")
        add_btn.clicked.connect(self.add_feature)
        feature_buttons.addWidget(add_btn)
        
        edit_btn = QPushButton("Edit")
        edit_btn.clicked.connect(self.edit_selected_feature)
        feature_buttons.addWidget(edit_btn)
        
        delete_btn = QPushButton("Delete")
        delete_btn.clicked.connect(self.delete_selected_feature)
        feature_buttons.addWidget(delete_btn)
        
        features_layout.addLayout(feature_buttons)
        
        right_layout.addWidget(features_group)
        
        # Feature details
        details_group = QGroupBox("Feature Details")
        details_layout = QVBoxLayout(details_group)
        
        self.feature_details = QTextEdit()
        self.feature_details.setReadOnly(True)
        self.feature_details.setMaximumHeight(150)
        details_layout.addWidget(self.feature_details)
        
        right_layout.addWidget(details_group)
        
        splitter.addWidget(right_panel)
        
        # Set splitter proportions
        splitter.setSizes([600, 300])
        
        # Dialog buttons
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        apply_button = QPushButton("Apply Changes")
        apply_button.clicked.connect(self.apply_changes)
        button_layout.addWidget(apply_button)
        
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept_changes)
        button_layout.addWidget(ok_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
    
    def load_sequence(self):
        """Load sequence data into the editor"""
        if not self.record:
            return
        
        self.name_edit.setText(self.record.id or "")
        self.description_edit.setText(self.record.description or "")
        self.sequence_edit.setPlainText(str(self.record.seq))
        self.update_length()
        self.load_features()
    
    def load_features(self):
        """Load features into the list"""
        self.features_list.clear()
        
        if not hasattr(self.record, 'features'):
            return
        
        for i, feature in enumerate(self.record.features):
            if feature.type == 'source':
                continue
            
            # Create display text
            label = feature.qualifiers.get('label', [feature.type])[0]
            location_str = f"{feature.location.start+1}..{feature.location.end}"
            strand_str = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "?"
            
            display_text = f"{feature.type}: {label} [{location_str}] ({strand_str})"
            
            item = QListWidgetItem(display_text)
            item.setData(Qt.UserRole, i)  # Store feature index
            self.features_list.addItem(item)
    
    def on_sequence_changed(self):
        """Handle sequence text changes"""
        self.update_length()
        
        # Update sequence in record
        if self.record:
            sequence_text = self.sequence_edit.toPlainText().replace(" ", "").replace("\n", "")
            self.record.seq = Seq(sequence_text.upper())
    
    def update_length(self):
        """Update sequence length display"""
        text = self.sequence_edit.toPlainText().replace(" ", "").replace("\n", "")
        self.length_label.setText(f"{len(text)} bp")
    
    def change_case(self, case_option):
        """Change sequence case"""
        text = self.sequence_edit.toPlainText()
        
        if case_option == "Upper Case":
            text = text.upper()
        elif case_option == "Lower Case":
            text = text.lower()
        
        self.sequence_edit.setPlainText(text)
    
    def reverse_complement(self):
        """Create reverse complement of sequence"""
        try:
            current_seq = self.sequence_edit.toPlainText().replace(" ", "").replace("\n", "")
            seq_obj = Seq(current_seq)
            rev_comp = seq_obj.reverse_complement()
            self.sequence_edit.setPlainText(str(rev_comp))
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not create reverse complement: {str(e)}")
    
    def find_sequence(self):
        """Find sequence dialog"""
        search_text, ok = QInputDialog.getText(self, "Find Sequence", "Enter sequence to find:")
        if ok and search_text:
            self.highlight_sequence(search_text)
    
    def find_replace(self):
        """Find and replace dialog"""
        dialog = QDialog(self)
        dialog.setWindowTitle("Find and Replace")
        layout = QFormLayout(dialog)
        
        find_edit = QLineEdit()
        replace_edit = QLineEdit()
        
        layout.addRow("Find:", find_edit)
        layout.addRow("Replace:", replace_edit)
        
        buttons = QHBoxLayout()
        replace_btn = QPushButton("Replace All")
        cancel_btn = QPushButton("Cancel")
        
        buttons.addWidget(replace_btn)
        buttons.addWidget(cancel_btn)
        layout.addRow(buttons)
        
        def replace_all():
            find_text = find_edit.text()
            replace_text = replace_edit.text()
            if find_text:
                current_text = self.sequence_edit.toPlainText()
                new_text = current_text.replace(find_text, replace_text)
                self.sequence_edit.setPlainText(new_text)
                dialog.accept()
        
        replace_btn.clicked.connect(replace_all)
        cancel_btn.clicked.connect(dialog.reject)
        
        dialog.exec_()
    
    def highlight_sequence(self, search_text):
        """Highlight found sequence"""
        cursor = self.sequence_edit.textCursor()
        document = self.sequence_edit.document()
        
        # Clear previous highlights
        cursor.select(QTextCursor.Document)
        format = QTextCharFormat()
        cursor.setCharFormat(format)
        
        # Find and highlight
        cursor.movePosition(QTextCursor.Start)
        color = QColor(255, 255, 0)  # Yellow highlight
        
        while True:
            cursor = document.find(search_text, cursor)
            if cursor.isNull():
                break
            
            format = QTextCharFormat()
            format.setBackground(color)
            cursor.setCharFormat(format)
    
    def add_feature(self):
        """Add new feature"""
        if not self.record:
            return
        
        seq_length = len(str(self.record.seq))
        dialog = FeatureEditor(sequence_length=seq_length, parent=self)
        
        if dialog.exec_() == QDialog.Accepted:
            new_feature = dialog.get_feature()
            if not hasattr(self.record, 'features'):
                self.record.features = []
            self.record.features.append(new_feature)
            self.load_features()
    
    def edit_selected_feature(self):
        """Edit selected feature"""
        current_item = self.features_list.currentItem()
        if not current_item or not self.record:
            return
        
        feature_index = current_item.data(Qt.UserRole)
        if feature_index is None or feature_index >= len(self.record.features):
            return
        
        feature = self.record.features[feature_index]
        seq_length = len(str(self.record.seq))
        
        dialog = FeatureEditor(feature, seq_length, parent=self)
        
        if dialog.exec_() == QDialog.Accepted:
            updated_feature = dialog.get_feature()
            self.record.features[feature_index] = updated_feature
            self.load_features()
    
    def delete_selected_feature(self):
        """Delete selected feature"""
        current_item = self.features_list.currentItem()
        if not current_item or not self.record:
            return
        
        feature_index = current_item.data(Qt.UserRole)
        if feature_index is None or feature_index >= len(self.record.features):
            return
        
        reply = QMessageBox.question(
            self, "Delete Feature", 
            "Are you sure you want to delete this feature?",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            del self.record.features[feature_index]
            self.load_features()
    
    def on_feature_selection_changed(self):
        """Handle feature selection changes"""
        current_item = self.features_list.currentItem()
        if not current_item or not self.record:
            self.feature_details.clear()
            return
        
        feature_index = current_item.data(Qt.UserRole)
        if feature_index is None or feature_index >= len(self.record.features):
            return
        
        feature = self.record.features[feature_index]
        
        # Display feature details
        details = f"Type: {feature.type}\n"
        details += f"Location: {feature.location.start+1}..{feature.location.end}\n"
        details += f"Strand: {'+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '?'}\n"
        details += f"Length: {len(feature.location)} bp\n\n"
        
        details += "Qualifiers:\n"
        for key, values in feature.qualifiers.items():
            details += f"  {key}: {', '.join(values)}\n"
        
        self.feature_details.setPlainText(details)
    
    def save_changes(self):
        """Save changes to record"""
        if not self.record:
            return
        
        # Update basic info
        self.record.id = self.name_edit.text()
        self.record.description = self.description_edit.text()
        
        # Sequence is already updated in on_sequence_changed
        
        QMessageBox.information(self, "Saved", "Changes saved successfully.")
    
    def revert_changes(self):
        """Revert to original sequence"""
        reply = QMessageBox.question(
            self, "Revert Changes", 
            "Are you sure you want to revert all changes?",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            self.record = copy.deepcopy(self.original_record)
            self.load_sequence()
    
    def apply_changes(self):
        """Apply changes and emit signal"""
        self.save_changes()
        self.sequenceModified.emit(self.record)
    
    def accept_changes(self):
        """Accept and close dialog"""
        self.apply_changes()
        self.accept()
    
    def get_modified_record(self):
        """Get the modified sequence record"""
        return self.record