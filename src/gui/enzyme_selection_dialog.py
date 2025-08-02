from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QListWidget, QDialogButtonBox, 
                             QAbstractItemView, QLabel, QLineEdit, QHBoxLayout,
                             QPushButton, QMessageBox, QCheckBox, QGroupBox)
from PyQt5.QtCore import Qt
from Bio.Restriction import Restriction

class EnzymeSelectionDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Restriction Enzymes")
        self.setMinimumSize(500, 600)

        layout = QVBoxLayout(self)
        
        # Instructions
        instructions = QLabel("Select restriction enzymes to analyze. Use Ctrl+Click to select multiple enzymes.")
        instructions.setWordWrap(True)
        instructions.setStyleSheet("QLabel { background-color: #f0f0f0; padding: 8px; border-radius: 4px; }")
        layout.addWidget(instructions)
        
        # Search box
        search_layout = QHBoxLayout()
        search_layout.addWidget(QLabel("Search:"))
        self.search_box = QLineEdit()
        self.search_box.setPlaceholderText("Type enzyme name to filter...")
        self.search_box.textChanged.connect(self.filter_enzymes)
        search_layout.addWidget(self.search_box)
        layout.addLayout(search_layout)
        
        # Quick selection buttons
        quick_select_group = QGroupBox("Quick Selection")
        quick_layout = QHBoxLayout(quick_select_group)
        
        common_button = QPushButton("Common Enzymes")
        common_button.setToolTip("Select commonly used restriction enzymes")
        common_button.clicked.connect(self.select_common_enzymes)
        quick_layout.addWidget(common_button)
        
        clear_button = QPushButton("Clear All")
        clear_button.clicked.connect(self.clear_selection)
        quick_layout.addWidget(clear_button)
        
        select_all_button = QPushButton("Select All Visible")
        select_all_button.clicked.connect(self.select_all_visible)
        quick_layout.addWidget(select_all_button)
        
        layout.addWidget(quick_select_group)

        # Enzyme list
        self.enzyme_list = QListWidget()
        self.enzyme_list.setSelectionMode(QAbstractItemView.MultiSelection)
        layout.addWidget(self.enzyme_list)
        
        # Store all enzymes for filtering
        self.all_enzymes = []
        
        try:
            # Populate with all commercially available enzymes from Biopython
            for enzyme in Restriction.CommOnly:
                enzyme_name = str(enzyme)
                self.all_enzymes.append(enzyme_name)
                self.enzyme_list.addItem(enzyme_name)
        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Could not load all restriction enzymes: {str(e)}")
            # Fallback to a basic set
            basic_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI', 'NotI', 'XbaI', 'SpeI', 'PstI']
            for enzyme in basic_enzymes:
                self.all_enzymes.append(enzyme)
                self.enzyme_list.addItem(enzyme)
        
        # Selection info
        self.selection_label = QLabel("0 enzymes selected")
        self.selection_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.selection_label)
        
        # Update selection count when selection changes
        self.enzyme_list.itemSelectionChanged.connect(self.update_selection_count)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        # Set initial selection count
        self.update_selection_count()
    
    def filter_enzymes(self):
        """Filter enzyme list based on search text"""
        search_text = self.search_box.text().lower()
        
        # Clear and repopulate list
        self.enzyme_list.clear()
        
        for enzyme in self.all_enzymes:
            if search_text in enzyme.lower():
                self.enzyme_list.addItem(enzyme)
    
    def select_common_enzymes(self):
        """Select a set of commonly used restriction enzymes"""
        common_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI', 'NotI', 'XbaI']
        
        # Clear current selection
        self.enzyme_list.clearSelection()
        
        # Select common enzymes if they exist in the list
        for i in range(self.enzyme_list.count()):
            item = self.enzyme_list.item(i)
            if item.text() in common_enzymes:
                item.setSelected(True)
    
    def clear_selection(self):
        """Clear all selections"""
        self.enzyme_list.clearSelection()
    
    def select_all_visible(self):
        """Select all visible enzymes in the list"""
        for i in range(self.enzyme_list.count()):
            self.enzyme_list.item(i).setSelected(True)
    
    def update_selection_count(self):
        """Update the selection count label"""
        count = len(self.enzyme_list.selectedItems())
        self.selection_label.setText(f"{count} enzyme{'s' if count != 1 else ''} selected")

    def get_selected_enzymes(self):
        """Get list of selected enzyme names"""
        return [item.text() for item in self.enzyme_list.selectedItems()]
