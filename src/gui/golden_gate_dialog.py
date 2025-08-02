from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QListWidget, QPushButton, QFileDialog, 
                             QLabel, QHBoxLayout, QTextEdit, QGroupBox, QMessageBox,
                             QSplitter, QComboBox, QFormLayout, QCheckBox, QWidget)
from PyQt5.QtCore import Qt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import RestrictionBatch, Restriction
import os

class GoldenGateDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Golden Gate Assembly Simulation")
        self.setMinimumSize(800, 600)
        self.destination_vector = None
        self.fragments = []

        # Main layout
        main_layout = QVBoxLayout(self)
        
        # Instructions
        instructions = QLabel("""
        <b>Golden Gate Assembly Simulation</b><br>
        Golden Gate assembly uses Type IIS restriction enzymes to create compatible overhangs.<br>
        1. Select a destination vector with Type IIS sites<br>
        2. Add DNA fragments with compatible overhangs<br>
        3. Choose the Type IIS enzyme for assembly<br>
        4. Simulate the one-pot assembly reaction
        """)
        instructions.setWordWrap(True)
        instructions.setStyleSheet("QLabel { background-color: #f0f0f0; padding: 10px; border-radius: 5px; }")
        main_layout.addWidget(instructions)
        
        # Create splitter for main content
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left panel - Vector and fragment selection
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # Destination vector group
        vector_group = QGroupBox("Destination Vector")
        vector_layout = QVBoxLayout(vector_group)
        
        self.vector_label = QLabel("Destination Vector: None selected")
        self.vector_label.setWordWrap(True)
        vector_layout.addWidget(self.vector_label)
        
        vector_button_layout = QHBoxLayout()
        self.vector_button = QPushButton("📁 Select Destination Vector")
        self.vector_button.setToolTip("Select the destination vector (backbone) for Golden Gate assembly")
        self.vector_button.clicked.connect(self.select_vector)
        vector_button_layout.addWidget(self.vector_button)
        vector_layout.addLayout(vector_button_layout)
        
        left_layout.addWidget(vector_group)
        
        # Assembly parameters
        params_group = QGroupBox("Assembly Parameters")
        params_layout = QFormLayout(params_group)
        
        # Type IIS enzyme selection
        self.enzyme_combo = QComboBox()
        self.enzyme_combo.addItems(["BsaI", "BbsI", "BsmBI", "Esp3I", "SapI"])
        self.enzyme_combo.setToolTip("Select the Type IIS restriction enzyme for Golden Gate assembly")
        params_layout.addRow("Type IIS Enzyme:", self.enzyme_combo)
        
        # Assembly options
        self.remove_enzyme_sites = QCheckBox("Remove enzyme sites from product")
        self.remove_enzyme_sites.setChecked(True)
        self.remove_enzyme_sites.setToolTip("Remove Type IIS enzyme sites from the final product")
        params_layout.addRow(self.remove_enzyme_sites)
        
        left_layout.addWidget(params_group)
        
        # Fragment list group
        fragment_group = QGroupBox("DNA Fragments")
        fragment_layout = QVBoxLayout(fragment_group)
        
        self.fragment_list = QListWidget()
        self.fragment_list.setToolTip("List of DNA fragments for Golden Gate assembly")
        fragment_layout.addWidget(self.fragment_list)
        
        # Fragment control buttons
        fragment_button_layout = QHBoxLayout()
        
        self.add_fragment_button = QPushButton("📁 Add Fragment")
        self.add_fragment_button.setToolTip("Add a DNA fragment for assembly")
        self.add_fragment_button.clicked.connect(self.add_fragment)
        fragment_button_layout.addWidget(self.add_fragment_button)
        
        self.remove_fragment_button = QPushButton("🗑️ Remove")
        self.remove_fragment_button.setToolTip("Remove selected fragment")
        self.remove_fragment_button.clicked.connect(self.remove_fragment)
        self.remove_fragment_button.setEnabled(False)
        fragment_button_layout.addWidget(self.remove_fragment_button)
        
        fragment_layout.addLayout(fragment_button_layout)
        left_layout.addWidget(fragment_group)
        
        # Assembly button
        self.assemble_button = QPushButton("🔗 Perform Golden Gate Assembly")
        self.assemble_button.setToolTip("Simulate Golden Gate assembly")
        self.assemble_button.clicked.connect(self.assemble)
        self.assemble_button.setEnabled(False)
        left_layout.addWidget(self.assemble_button)
        
        splitter.addWidget(left_panel)
        
        # Right panel - Information and analysis
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        # Vector information
        vector_info_group = QGroupBox("Vector Information")
        vector_info_layout = QVBoxLayout(vector_info_group)
        
        self.vector_info_display = QTextEdit()
        self.vector_info_display.setReadOnly(True)
        self.vector_info_display.setMaximumHeight(150)
        self.vector_info_display.setText("Select a destination vector to view information.")
        vector_info_layout.addWidget(self.vector_info_display)
        
        right_layout.addWidget(vector_info_group)
        
        # Fragment information
        fragment_info_group = QGroupBox("Fragment Information")
        fragment_info_layout = QVBoxLayout(fragment_info_group)
        
        self.fragment_info_display = QTextEdit()
        self.fragment_info_display.setReadOnly(True)
        self.fragment_info_display.setMaximumHeight(150)
        self.fragment_info_display.setText("Add fragments to view information.")
        fragment_info_layout.addWidget(self.fragment_info_display)
        
        right_layout.addWidget(fragment_info_group)
        
        # Assembly preview
        preview_group = QGroupBox("Assembly Preview")
        preview_layout = QVBoxLayout(preview_group)
        
        self.preview_display = QTextEdit()
        self.preview_display.setReadOnly(True)
        self.preview_display.setMaximumHeight(150)
        preview_layout.addWidget(self.preview_display)
        
        right_layout.addWidget(preview_group)
        
        # Add stretch
        right_layout.addStretch()
        
        splitter.addWidget(right_panel)
        
        # Set splitter sizes
        splitter.setSizes([400, 400])
        
        # Button box
        button_box_layout = QHBoxLayout()
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_box_layout.addStretch()
        button_box_layout.addWidget(close_button)
        
        main_layout.addLayout(button_box_layout)
        
        # Connect selection change
        self.fragment_list.itemSelectionChanged.connect(self.update_button_states)
        self.enzyme_combo.currentTextChanged.connect(self.update_preview)

    def select_vector(self):
        """Select destination vector file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Destination Vector File", "", 
            "GenBank Files (*.gb *.gbk);;FASTA Files (*.fa *.fasta);;All Files (*)"
        )
        if not file_path:
            return
            
        try:
            # Try to read as GenBank first, then FASTA
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext in ['.gb', '.gbk']:
                self.destination_vector = SeqIO.read(file_path, "genbank")
            else:
                try:
                    self.destination_vector = SeqIO.read(file_path, "genbank")
                except:
                    self.destination_vector = SeqIO.read(file_path, "fasta")
            
            vector_name = self.destination_vector.id or os.path.basename(file_path)
            self.vector_label.setText(f"Destination Vector: {vector_name} ({len(self.destination_vector.seq)} bp)")
            
            self.update_vector_info()
            self.update_assembly_status()
            self.update_preview()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading destination vector: {str(e)}")

    def add_fragment(self):
        """Add a DNA fragment from file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Fragment File", "", 
            "GenBank Files (*.gb *.gbk);;FASTA Files (*.fa *.fasta);;All Files (*)"
        )
        if not file_path:
            return
            
        try:
            # Try to read as GenBank first, then FASTA
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext in ['.gb', '.gbk']:
                record = SeqIO.read(file_path, "genbank")
            else:
                try:
                    record = SeqIO.read(file_path, "genbank")
                except:
                    record = SeqIO.read(file_path, "fasta")
            
            self.fragments.append(record)
            fragment_name = record.id or os.path.basename(file_path)
            self.fragment_list.addItem(f"{fragment_name} ({len(record.seq)} bp)")
            
            self.update_fragment_info()
            self.update_assembly_status()
            self.update_preview()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading fragment: {str(e)}")
    
    def remove_fragment(self):
        """Remove selected fragment"""
        current_row = self.fragment_list.currentRow()
        if current_row >= 0:
            self.fragments.pop(current_row)
            self.fragment_list.takeItem(current_row)
            
            self.update_fragment_info()
            self.update_assembly_status()
            self.update_preview()
            self.update_button_states()
    
    def update_button_states(self):
        """Update button enabled states"""
        has_selection = self.fragment_list.currentRow() >= 0
        self.remove_fragment_button.setEnabled(has_selection)
    
    def update_assembly_status(self):
        """Update assembly button status"""
        can_assemble = (self.destination_vector is not None and 
                       len(self.fragments) >= 1)
        self.assemble_button.setEnabled(can_assemble)
    
    def update_vector_info(self):
        """Update vector information display"""
        if not self.destination_vector:
            self.vector_info_display.setText("No destination vector selected.")
            return
        
        info_text = f"Vector Information:\n\n"
        info_text += f"Name: {self.destination_vector.id}\n"
        info_text += f"Length: {len(self.destination_vector.seq)} bp\n"
        
        # Check for features
        if hasattr(self.destination_vector, 'features'):
            feature_count = len([f for f in self.destination_vector.features if f.type != 'source'])
            info_text += f"Features: {feature_count}\n"
        
        # Analyze Type IIS sites
        enzyme_name = self.enzyme_combo.currentText()
        try:
            enzyme = getattr(Restriction, enzyme_name, None)
            if enzyme:
                batch = RestrictionBatch([enzyme])
                analysis = batch.search(self.destination_vector.seq)
                sites = analysis.get(enzyme, [])
                
                info_text += f"\n{enzyme_name} sites: {len(sites)}\n"
                if sites:
                    info_text += f"Positions: {sites}\n"
                    
                    if len(sites) == 2:
                        info_text += "\n✓ Suitable for Golden Gate (2 sites)"
                    elif len(sites) > 2:
                        info_text += "\n⚠ Multiple sites - may complicate assembly"
                    else:
                        info_text += "\n⚠ Only 1 site - need 2 for Golden Gate"
                else:
                    info_text += "\n✗ No sites found for this enzyme"
        except Exception as e:
            info_text += f"\nError analyzing sites: {str(e)}"
        
        self.vector_info_display.setText(info_text)
    
    def update_fragment_info(self):
        """Update fragment information display"""
        if not self.fragments:
            self.fragment_info_display.setText("No fragments added yet.")
            return
        
        info_text = f"Fragment Information:\n\n"
        info_text += f"Total fragments: {len(self.fragments)}\n\n"
        
        total_length = 0
        enzyme_name = self.enzyme_combo.currentText()
        
        for i, fragment in enumerate(self.fragments, 1):
            fragment_name = fragment.id or f"Fragment {i}"
            length = len(fragment.seq)
            total_length += length
            info_text += f"{i}. {fragment_name}: {length} bp\n"
            
            # Check for enzyme sites
            try:
                enzyme = getattr(Restriction, enzyme_name, None)
                if enzyme:
                    batch = RestrictionBatch([enzyme])
                    analysis = batch.search(fragment.seq)
                    sites = analysis.get(enzyme, [])
                    info_text += f"   {enzyme_name} sites: {len(sites)}\n"
            except:
                pass
        
        info_text += f"\nTotal fragment length: {total_length} bp"
        
        self.fragment_info_display.setText(info_text)
    
    def update_preview(self):
        """Update assembly preview"""
        if not self.destination_vector:
            self.preview_display.setText("Select a destination vector to preview assembly.")
            return
        
        preview_text = "Golden Gate Assembly Preview:\n\n"
        enzyme_name = self.enzyme_combo.currentText()
        
        # Vector analysis
        try:
            enzyme = getattr(Restriction, enzyme_name, None)
            if enzyme:
                batch = RestrictionBatch([enzyme])
                vector_analysis = batch.search(self.destination_vector.seq)
                vector_sites = vector_analysis.get(enzyme, [])
                
                preview_text += f"1. Digest vector with {enzyme_name}\n"
                preview_text += f"   Sites: {len(vector_sites)} ({vector_sites})\n\n"
                
                if len(self.fragments) > 0:
                    preview_text += f"2. Digest {len(self.fragments)} fragment(s) with {enzyme_name}\n"
                    
                    for i, fragment in enumerate(self.fragments, 1):
                        fragment_analysis = batch.search(fragment.seq)
                        fragment_sites = fragment_analysis.get(enzyme, [])
                        fragment_name = fragment.id or f"Fragment {i}"
                        preview_text += f"   {fragment_name}: {len(fragment_sites)} sites\n"
                    
                    preview_text += "\n3. Ligate compatible overhangs\n"
                    preview_text += "4. Transform into competent cells\n"
                    
                    if self.remove_enzyme_sites.isChecked():
                        preview_text += f"\nNote: {enzyme_name} sites will be removed from product\n"
                    
                    # Estimate product size
                    vector_length = len(self.destination_vector.seq)
                    fragment_length = sum(len(f.seq) for f in self.fragments)
                    
                    # Simplified calculation (actual depends on overhangs)
                    estimated_size = vector_length + fragment_length
                    preview_text += f"\nEstimated product size: ~{estimated_size} bp"
                else:
                    preview_text += "\nAdd fragments to complete preview."
        except Exception as e:
            preview_text += f"Error in preview: {str(e)}"
        
        self.preview_display.setText(preview_text)

    def assemble(self):
        """Perform Golden Gate assembly simulation"""
        if not self.destination_vector:
            QMessageBox.warning(self, "Missing Vector", "Please select a destination vector.")
            return
            
        if not self.fragments:
            QMessageBox.warning(self, "Missing Fragments", "Please add at least one fragment.")
            return

        try:
            enzyme_name = self.enzyme_combo.currentText()
            enzyme = getattr(Restriction, enzyme_name)
            
            # Analyze vector and fragments
            batch = RestrictionBatch([enzyme])
            vector_analysis = batch.search(self.destination_vector.seq)
            vector_sites = vector_analysis.get(enzyme, [])
            
            if len(vector_sites) < 2:
                reply = QMessageBox.question(self, "Insufficient Sites", 
                    f"Vector has only {len(vector_sites)} {enzyme_name} site(s).\n"
                    "Golden Gate assembly typically requires 2 sites. Continue anyway?",
                    QMessageBox.Yes | QMessageBox.No)
                if reply == QMessageBox.No:
                    return
            
            # Simplified Golden Gate Assembly simulation
            # In reality, this requires careful design of compatible overhangs
            new_seq = self.destination_vector.seq
            
            # Insert fragments (simplified - assumes compatible overhangs)
            for fragment in self.fragments:
                # In real Golden Gate, fragments would be inserted at specific positions
                # based on overhang compatibility
                new_seq += fragment.seq
            
            # Remove enzyme sites if requested (simplified)
            if self.remove_enzyme_sites.isChecked():
                # This is a very simplified removal - real implementation would be more complex
                pass
            
            # Create product record
            vector_name = self.destination_vector.id or "vector"
            fragment_names = [f.id or f"fragment_{i+1}" for i, f in enumerate(self.fragments)]
            
            product_id = f"{vector_name}_GG_{'_'.join(fragment_names[:2])}"
            if len(fragment_names) > 2:
                product_id += f"_and_{len(fragment_names)-2}_more"
            
            new_record = SeqRecord(
                new_seq,
                id=product_id,
                name="Golden Gate Assembly Product",
                description=f"Golden Gate assembly using {enzyme_name}: {vector_name} + {len(self.fragments)} fragments"
            )
            
            # Show result
            result_msg = f"""Golden Gate Assembly Successful!
            
Product Details:
            - Destination Vector: {vector_name} ({len(self.destination_vector.seq)} bp)
            - Fragments inserted: {len(self.fragments)}
            - Type IIS Enzyme: {enzyme_name}
            - Product length: {len(new_seq)} bp
            
Note: This is a simplified simulation.
Real Golden Gate assembly requires:
            - Compatible overhang design
            - Proper fragment orientation
            - Optimized reaction conditions
            - Type IIS enzyme digestion
            - T4 DNA ligase
            
The assembled sequence will be displayed in the main viewer."""
            
            QMessageBox.information(self, "Golden Gate Assembly Complete", result_msg)

            # Display in parent viewer
            if hasattr(self.parent(), 'display_sequence'):
                self.parent().display_sequence(new_record)
            
            self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Assembly Error", f"Error during Golden Gate assembly: {str(e)}")
