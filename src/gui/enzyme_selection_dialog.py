from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QListWidget, QDialogButtonBox, 
                             QAbstractItemView, QLabel, QLineEdit, QHBoxLayout,
                             QPushButton, QMessageBox, QCheckBox, QGroupBox)
from PyQt6.QtCore import Qt
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
        quick_select_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 13px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        quick_layout = QHBoxLayout(quick_select_group)
        quick_layout.setContentsMargins(10, 15, 10, 10)
        quick_layout.setSpacing(8)
        
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
        self.enzyme_list.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        layout.addWidget(self.enzyme_list)
        
        # Store all enzymes for filtering
        self.all_enzymes = []
        
        try:
            # Try to load enzymes more comprehensively
            all_enzyme_names = []
            
            # First try to get all commercially available enzymes
            try:
                if hasattr(Restriction, 'CommOnly'):
                    for enzyme in Restriction.CommOnly:
                        all_enzyme_names.append(str(enzyme))
                else:
                    # Fallback: try to get enzymes from Restriction.__dict__
                    for name in dir(Restriction):
                        obj = getattr(Restriction, name)
                        if hasattr(obj, 'site') and hasattr(obj, 'fst5'):
                            all_enzyme_names.append(name)
            except:
                pass
            
            # If we still don't have enough enzymes, add common ones manually
            if len(all_enzyme_names) < 20:
                common_enzymes = [
                    'EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI', 'NotI', 
                    'XbaI', 'SpeI', 'PstI', 'SalI', 'SmaI', 'BglII', 'NcoI', 
                    'NdeI', 'ApaI', 'EcoRV', 'NruI', 'SphI', 'ClaI'
                ]
                for enzyme in common_enzymes:
                    if enzyme not in all_enzyme_names:
                        all_enzyme_names.append(enzyme)
            
            # Sort the enzymes alphabetically
            all_enzyme_names.sort()
            
            # Add to list
            for enzyme_name in all_enzyme_names:
                self.all_enzymes.append(enzyme_name)
                self.enzyme_list.addItem(enzyme_name)
                
        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Could not load restriction enzymes: {str(e)}")
            # Final fallback to a basic set
            basic_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI', 'NotI', 'XbaI', 'SpeI', 'PstI']
            for enzyme in basic_enzymes:
                self.all_enzymes.append(enzyme)
                self.enzyme_list.addItem(enzyme)
        
        # Selection info
        self.selection_label = QLabel("0 enzymes selected")
        self.selection_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.selection_label)
        
        # Update selection count when selection changes
        self.enzyme_list.itemSelectionChanged.connect(self.update_selection_count)

        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
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
