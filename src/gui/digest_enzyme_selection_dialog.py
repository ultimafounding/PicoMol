from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QListWidget, 
                             QDialogButtonBox, QLabel, QPushButton, QGroupBox,
                             QLineEdit, QTextEdit, QCheckBox, QSplitter, QWidget,
                             QAbstractItemView)
from PyQt6.QtCore import Qt
from Bio.Restriction import Restriction, RestrictionBatch
import os

class DigestEnzymeSelectionDialog(QDialog):
    def __init__(self, current_enzymes, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Enzymes for Virtual Digest")
        self.setMinimumSize(600, 500)
        self.current_enzymes = current_enzymes or []
        
        layout = QVBoxLayout(self)
        
        # Instructions
        instructions = QLabel(
            "Select restriction enzymes for virtual digest. "
            "Choose from currently selected enzymes or add custom enzymes."
        )
        instructions.setWordWrap(True)
        layout.addWidget(instructions)
        
        # Create splitter for enzyme lists
        splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(splitter)
        
        # Left panel - Current enzymes
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        current_group = QGroupBox("Current Enzyme Selection")
        current_group.setStyleSheet("""
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
        current_layout = QVBoxLayout(current_group)
        current_layout.setContentsMargins(10, 15, 10, 10)
        current_layout.setSpacing(8)
        
        self.current_list = QListWidget()
        self.current_list.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        current_layout.addWidget(self.current_list)
        
        # Quick selection buttons for current enzymes
        quick_layout = QHBoxLayout()
        select_all_current = QPushButton("Select All")
        select_all_current.clicked.connect(self.select_all_current)
        clear_current = QPushButton("Clear")
        clear_current.clicked.connect(self.clear_current)
        quick_layout.addWidget(select_all_current)
        quick_layout.addWidget(clear_current)
        current_layout.addLayout(quick_layout)
        
        left_layout.addWidget(current_group)
        splitter.addWidget(left_panel)
        
        # Right panel - Additional enzymes
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        additional_group = QGroupBox("Additional Enzymes")
        additional_group.setStyleSheet(current_group.styleSheet())
        additional_layout = QVBoxLayout(additional_group)
        additional_layout.setContentsMargins(10, 15, 10, 10)
        additional_layout.setSpacing(8)
        
        # Search box
        self.search_box = QLineEdit()
        self.search_box.setPlaceholderText("Search enzymes...")
        self.search_box.textChanged.connect(self.filter_enzymes)
        additional_layout.addWidget(QLabel("Search:"))
        additional_layout.addWidget(self.search_box)
        
        self.additional_list = QListWidget()
        self.additional_list.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        additional_layout.addWidget(self.additional_list)
        
        # Quick selection buttons for additional enzymes
        additional_quick_layout = QHBoxLayout()
        select_common = QPushButton("Common Enzymes")
        select_common.clicked.connect(self.select_common_enzymes)
        select_all_additional = QPushButton("Select All")
        select_all_additional.clicked.connect(self.select_all_additional)
        clear_additional = QPushButton("Clear")
        clear_additional.clicked.connect(self.clear_additional)
        additional_quick_layout.addWidget(select_common)
        additional_quick_layout.addWidget(select_all_additional)
        additional_quick_layout.addWidget(clear_additional)
        additional_layout.addLayout(additional_quick_layout)
        
        right_layout.addWidget(additional_group)
        splitter.addWidget(right_panel)
        
        # Set splitter sizes
        splitter.setSizes([300, 300])
        
        # Selection info
        self.selection_label = QLabel("0 enzymes selected")
        self.selection_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.selection_label)
        
        # Button box
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        # Load enzymes
        self.load_enzymes()
        
        # Connect selection changes
        self.current_list.itemSelectionChanged.connect(self.update_selection_count)
        self.additional_list.itemSelectionChanged.connect(self.update_selection_count)
        
        # Set initial selection count
        self.update_selection_count()
    
    def load_enzymes(self):
        """Load enzymes into the lists"""
        # Load current enzymes
        for enzyme in self.current_enzymes:
            self.current_list.addItem(str(enzyme))
        
        # Load additional enzymes
        try:
            all_enzyme_names = []
            
            # Try to get all commercially available enzymes
            try:
                if hasattr(Restriction, 'CommOnly'):
                    for enzyme in Restriction.CommOnly:
                        enzyme_name = str(enzyme)
                        if enzyme_name not in self.current_enzymes:
                            all_enzyme_names.append(enzyme_name)
                else:
                    # Fallback: try to get enzymes from Restriction.__dict__
                    for name in dir(Restriction):
                        obj = getattr(Restriction, name)
                        if hasattr(obj, 'site') and hasattr(obj, 'fst5'):
                            if name not in self.current_enzymes:
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
                    if enzyme not in self.current_enzymes and enzyme not in all_enzyme_names:
                        all_enzyme_names.append(enzyme)
            
            # Sort the enzymes alphabetically
            all_enzyme_names.sort()
            
            # Add to additional list
            for enzyme_name in all_enzyme_names:
                self.additional_list.addItem(enzyme_name)
                
        except Exception as e:
            # Fallback to basic set
            basic_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI', 'NotI', 'XbaI', 'SpeI', 'PstI']
            for enzyme in basic_enzymes:
                if enzyme not in self.current_enzymes:
                    self.additional_list.addItem(enzyme)
    
    def filter_enzymes(self):
        """Filter the additional enzymes list based on search"""
        search_text = self.search_box.text().lower()
        self.additional_list.clear()
        
        try:
            all_enzyme_names = []
            
            # Try to get all commercially available enzymes
            try:
                if hasattr(Restriction, 'CommOnly'):
                    for enzyme in Restriction.CommOnly:
                        enzyme_name = str(enzyme)
                        if enzyme_name not in self.current_enzymes and search_text in enzyme_name.lower():
                            all_enzyme_names.append(enzyme_name)
                else:
                    # Fallback: try to get enzymes from Restriction.__dict__
                    for name in dir(Restriction):
                        obj = getattr(Restriction, name)
                        if hasattr(obj, 'site') and hasattr(obj, 'fst5'):
                            if name not in self.current_enzymes and search_text in name.lower():
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
                    if enzyme not in self.current_enzymes and search_text in enzyme.lower():
                        all_enzyme_names.append(enzyme)
            
            # Sort the enzymes alphabetically
            all_enzyme_names.sort()
            
            # Add to additional list
            for enzyme_name in all_enzyme_names:
                self.additional_list.addItem(enzyme_name)
                
        except Exception as e:
            # Fallback to basic set
            basic_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI', 'NotI', 'XbaI', 'SpeI', 'PstI']
            for enzyme in basic_enzymes:
                if enzyme not in self.current_enzymes and search_text in enzyme.lower():
                    self.additional_list.addItem(enzyme)
    
    def select_all_current(self):
        """Select all current enzymes"""
        for i in range(self.current_list.count()):
            item = self.current_list.item(i)
            item.setSelected(True)
    
    def clear_current(self):
        """Clear current enzyme selection"""
        self.current_list.clearSelection()
    
    def select_common_enzymes(self):
        """Select common enzymes from additional list"""
        common_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI', 'NotI', 'XbaI', 'SpeI', 'PstI']
        for i in range(self.additional_list.count()):
            item = self.additional_list.item(i)
            enzyme_name = item.text()
            if enzyme_name in common_enzymes:
                item.setSelected(True)
    
    def select_all_additional(self):
        """Select all additional enzymes"""
        for i in range(self.additional_list.count()):
            item = self.additional_list.item(i)
            item.setSelected(True)
    
    def clear_additional(self):
        """Clear additional enzyme selection"""
        self.additional_list.clearSelection()
    
    def update_selection_count(self):
        """Update the selection count label"""
        current_selected = len(self.current_list.selectedItems())
        additional_selected = len(self.additional_list.selectedItems())
        total_selected = current_selected + additional_selected
        self.selection_label.setText(f"{total_selected} enzymes selected ({current_selected} current, {additional_selected} additional)")
    
    def get_selected_enzymes(self):
        """Get the list of selected enzymes"""
        selected_enzymes = []
        
        # Get selected from current list
        for item in self.current_list.selectedItems():
            selected_enzymes.append(item.text())
        
        # Get selected from additional list
        for item in self.additional_list.selectedItems():
            selected_enzymes.append(item.text())
        
        return selected_enzymes
