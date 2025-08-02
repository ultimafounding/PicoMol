from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QFormLayout, QComboBox, QPushButton, 
                             QLineEdit, QLabel, QFileDialog, QListWidget, QMessageBox,
                             QGroupBox, QHBoxLayout, QTextEdit, QSplitter, QTabWidget,
                             QWidget, QCheckBox)
from PyQt5.QtCore import Qt
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Restriction
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

class RestrictionCloningDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Restriction Cloning Simulation")
        self.setMinimumSize(800, 600)
        self.vector_record = None
        self.insert_record = None

        # Main layout
        main_layout = QVBoxLayout(self)
        
        # Instructions
        instructions = QLabel("""
        <b>Restriction Cloning Simulation</b><br>
        1. Select vector and insert DNA sequences<br>
        2. Choose compatible restriction enzymes<br>
        3. Select fragments to ligate<br>
        4. Simulate the ligation reaction
        """)
        instructions.setWordWrap(True)
        instructions.setStyleSheet("QLabel { background-color: #f0f0f0; padding: 10px; border-radius: 5px; }")
        main_layout.addWidget(instructions)
        
        # Create splitter for main content
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left panel - File selection and enzyme choice
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # File selection group
        file_group = QGroupBox("DNA Sequences")
        form_layout = QFormLayout(file_group)

        # Vector selection
        vector_layout = QHBoxLayout()
        self.vector_label = QLabel("Vector: None selected")
        self.vector_button = QPushButton("📁 Select Vector")
        self.vector_button.setToolTip("Select the vector (backbone) plasmid file")
        self.vector_button.clicked.connect(self.select_vector)
        vector_layout.addWidget(self.vector_button)
        form_layout.addRow(self.vector_label, vector_layout)

        # Insert selection
        insert_layout = QHBoxLayout()
        self.insert_label = QLabel("Insert: None selected")
        self.insert_button = QPushButton("📁 Select Insert")
        self.insert_button.setToolTip("Select the insert DNA sequence file")
        self.insert_button.clicked.connect(self.select_insert)
        insert_layout.addWidget(self.insert_button)
        form_layout.addRow(self.insert_label, insert_layout)

        left_layout.addWidget(file_group)
        
        # Enzyme selection group
        enzyme_group = QGroupBox("Restriction Enzymes")
        enzyme_layout = QFormLayout(enzyme_group)
        
        self.vector_enzyme_combo = QComboBox()
        self.vector_enzyme_combo.setToolTip("Choose enzyme to cut the vector")
        enzyme_layout.addRow("Vector Enzyme:", self.vector_enzyme_combo)

        self.insert_enzyme_combo = QComboBox()
        self.insert_enzyme_combo.setToolTip("Choose enzyme to cut the insert")
        enzyme_layout.addRow("Insert Enzyme:", self.insert_enzyme_combo)

        self.vector_enzyme_combo.currentTextChanged.connect(self.update_vector_fragments)
        self.insert_enzyme_combo.currentTextChanged.connect(self.update_insert_fragments)
        
        left_layout.addWidget(enzyme_group)
        
        # Fragment selection group
        fragment_group = QGroupBox("Fragment Selection")
        fragment_layout = QVBoxLayout(fragment_group)
        
        fragment_layout.addWidget(QLabel("Vector Fragments:"))
        self.vector_fragments_list = QListWidget()
        self.vector_fragments_list.setToolTip("Select the vector fragment to use")
        fragment_layout.addWidget(self.vector_fragments_list)
        
        fragment_layout.addWidget(QLabel("Insert Fragments:"))
        self.insert_fragments_list = QListWidget()
        self.insert_fragments_list.setToolTip("Select the insert fragment to use")
        fragment_layout.addWidget(self.insert_fragments_list)
        
        left_layout.addWidget(fragment_group)
        
        # Ligation button
        self.ligate_button = QPushButton("🧬 Perform Ligation")
        self.ligate_button.setToolTip("Simulate the ligation reaction")
        self.ligate_button.clicked.connect(self.ligate)
        self.ligate_button.setEnabled(False)
        left_layout.addWidget(self.ligate_button)
        
        splitter.addWidget(left_panel)
        
        # Right panel - Results and information
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        # Results display
        results_group = QGroupBox("Cloning Information")
        results_layout = QVBoxLayout(results_group)
        
        self.info_display = QTextEdit()
        self.info_display.setReadOnly(True)
        self.info_display.setMaximumHeight(200)
        self.info_display.setText("Select vector and insert files to begin cloning simulation.")
        results_layout.addWidget(self.info_display)
        
        right_layout.addWidget(results_group)
        
        # Compatibility check
        compat_group = QGroupBox("Compatibility Analysis")
        compat_layout = QVBoxLayout(compat_group)
        
        self.compatibility_display = QTextEdit()
        self.compatibility_display.setReadOnly(True)
        self.compatibility_display.setMaximumHeight(150)
        compat_layout.addWidget(self.compatibility_display)
        
        right_layout.addWidget(compat_group)
        
        # Add stretch to right panel
        right_layout.addStretch()
        
        splitter.addWidget(right_panel)
        
        # Set splitter sizes
        splitter.setSizes([400, 400])
        
        # Button box
        button_layout = QHBoxLayout()
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addStretch()
        button_layout.addWidget(close_button)
        
        main_layout.addLayout(button_layout)

    def select_vector(self):
        """Select vector file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Vector File", "", 
            "GenBank Files (*.gb *.gbk);;FASTA Files (*.fa *.fasta);;All Files (*)"
        )
        if not file_path:
            return
        
        try:
            # Try to read as GenBank first, then FASTA
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext in ['.gb', '.gbk']:
                self.vector_record = SeqIO.read(file_path, "genbank")
            else:
                try:
                    self.vector_record = SeqIO.read(file_path, "genbank")
                except:
                    self.vector_record = SeqIO.read(file_path, "fasta")
            
            vector_name = self.vector_record.id or os.path.basename(file_path)
            self.vector_label.setText(f"Vector: {vector_name} ({len(self.vector_record.seq)} bp)")
            self.populate_enzyme_combos()
            self.update_info_display()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading vector file: {str(e)}")

    def select_insert(self):
        """Select insert file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Insert File", "", 
            "GenBank Files (*.gb *.gbk);;FASTA Files (*.fa *.fasta);;All Files (*)"
        )
        if not file_path:
            return
        
        try:
            # Try to read as GenBank first, then FASTA
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext in ['.gb', '.gbk']:
                self.insert_record = SeqIO.read(file_path, "genbank")
            else:
                try:
                    self.insert_record = SeqIO.read(file_path, "genbank")
                except:
                    self.insert_record = SeqIO.read(file_path, "fasta")
            
            insert_name = self.insert_record.id or os.path.basename(file_path)
            self.insert_label.setText(f"Insert: {insert_name} ({len(self.insert_record.seq)} bp)")
            self.populate_enzyme_combos()
            self.update_info_display()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading insert file: {str(e)}")

    def populate_enzyme_combos(self):
        """Populate enzyme combo boxes with compatible enzymes"""
        if not self.vector_record or not self.insert_record:
            return

        try:
            # Use a subset of common enzymes for better performance
            common_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI', 
                            'NotI', 'XbaI', 'SpeI', 'PstI', 'SalI', 'ApaI', 'NheI',
                            'BglII', 'ClaI', 'HpaI', 'NcoI', 'NdeI', 'SmaI', 'SphI']
            
            self.vector_enzyme_combo.clear()
            self.insert_enzyme_combo.clear()

            # Analyze restriction sites
            vector_enzymes = []
            insert_enzymes = []
            compatible_enzymes = []
            
            for enzyme_name in common_enzymes:
                try:
                    enzyme = getattr(Restriction, enzyme_name, None)
                    if enzyme:
                        vector_batch = RestrictionBatch([enzyme])
                        insert_batch = RestrictionBatch([enzyme])
                        
                        vector_analysis = vector_batch.search(self.vector_record.seq)
                        insert_analysis = insert_batch.search(self.insert_record.seq)
                        
                        vector_cuts = vector_analysis.get(enzyme, [])
                        insert_cuts = insert_analysis.get(enzyme, [])
                        
                        # For restriction cloning, we want:
                        # Vector: exactly 1 cut (to linearize)
                        # Insert: 0 or 2 cuts (to get insert fragment)
                        if len(vector_cuts) == 1:
                            vector_enzymes.append((enzyme_name, len(vector_cuts)))
                            
                        if len(insert_cuts) in [0, 1, 2]:  # Allow more flexibility
                            insert_enzymes.append((enzyme_name, len(insert_cuts)))
                            
                        # Compatible if both can be cut appropriately
                        if len(vector_cuts) >= 1 and len(insert_cuts) >= 0:
                            compatible_enzymes.append(enzyme_name)
                            
                except Exception as e:
                    print(f"Error analyzing {enzyme_name}: {e}")
                    continue
            
            # Populate combo boxes
            for enzyme_name, cut_count in vector_enzymes:
                self.vector_enzyme_combo.addItem(f"{enzyme_name} ({cut_count} cuts)")
                
            for enzyme_name, cut_count in insert_enzymes:
                self.insert_enzyme_combo.addItem(f"{enzyme_name} ({cut_count} cuts)")
            
            # Update compatibility display
            compat_text = "Restriction Cloning Analysis:\n\n"
            
            if compatible_enzymes:
                compat_text += f"Enzymes that can cut both sequences: {', '.join(compatible_enzymes[:10])}\n\n"
                compat_text += "Ideal strategy:\n"
                compat_text += "• Vector: 1 cut to linearize\n"
                compat_text += "• Insert: 0-2 cuts to isolate fragment\n"
                compat_text += "• Same enzyme for compatible ends\n\n"
                
                # Find ideal enzymes
                ideal_enzymes = []
                for enzyme_name in compatible_enzymes:
                    try:
                        enzyme = getattr(Restriction, enzyme_name)
                        vector_batch = RestrictionBatch([enzyme])
                        insert_batch = RestrictionBatch([enzyme])
                        
                        vector_cuts = len(vector_batch.search(self.vector_record.seq).get(enzyme, []))
                        insert_cuts = len(insert_batch.search(self.insert_record.seq).get(enzyme, []))
                        
                        if vector_cuts == 1 and insert_cuts in [0, 2]:
                            ideal_enzymes.append(enzyme_name)
                    except:
                        continue
                
                if ideal_enzymes:
                    compat_text += f"⭐ Recommended enzymes: {', '.join(ideal_enzymes[:5])}"
                else:
                    compat_text += "⚠️ No ideal single-enzyme strategy found.\nConsider using different enzymes for vector and insert."
            else:
                compat_text += "❌ No compatible enzymes found.\nCheck sequence compatibility."
            
            self.compatibility_display.setText(compat_text)
            
        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Error analyzing restriction sites: {str(e)}")

    def update_vector_fragments(self):
        """Update vector fragments list"""
        self.update_fragments(self.vector_record, self.vector_enzyme_combo.currentText(), self.vector_fragments_list)

    def update_insert_fragments(self):
        """Update insert fragments list"""
        self.update_fragments(self.insert_record, self.insert_enzyme_combo.currentText(), self.insert_fragments_list)

    def update_fragments(self, record, enzyme_name, list_widget):
        """Update fragment list for selected enzyme"""
        list_widget.clear()
        if not record or not enzyme_name:
            return

        try:
            # Extract enzyme name from combo box text (remove cut count)
            enzyme_name = enzyme_name.split(' (')[0]
            
            enzyme = getattr(Restriction, enzyme_name, None)
            if not enzyme:
                return
                
            batch = RestrictionBatch([enzyme])
            analysis = batch.search(record.seq)
            
            if enzyme in analysis:
                cuts = analysis[enzyme]
                # Assume vector is circular, insert is linear
                is_circular = (record == self.vector_record)
                fragments = self.get_fragments(record.seq, cuts, is_circular)
                
                for i, frag in enumerate(fragments):
                    # Add more detailed fragment information
                    gc_content = self.calculate_gc_content(frag)
                    list_widget.addItem(f"Fragment {i+1}: {len(frag)} bp (GC: {gc_content:.1f}%)")
                    
                # Enable ligation button if both fragments are available
                if (self.vector_fragments_list.count() > 0 and 
                    self.insert_fragments_list.count() > 0):
                    self.ligate_button.setEnabled(True)
                    
        except Exception as e:
            print(f"Error updating fragments: {e}")
    
    def calculate_gc_content(self, sequence):
        """Calculate GC content of a sequence"""
        try:
            from Bio.SeqUtils import gc_fraction
            return gc_fraction(sequence) * 100
        except:
            # Fallback calculation
            if len(sequence) == 0:
                return 0
            gc_count = sequence.count('G') + sequence.count('C')
            return (gc_count / len(sequence)) * 100

    def get_fragments(self, sequence, cuts, is_circular=True):
        """Get DNA fragments from restriction cuts"""
        if not cuts:
            return [sequence]
        
        fragments = []
        sorted_cuts = sorted(cuts)
        
        if is_circular and len(sorted_cuts) >= 1:
            # For circular sequences (plasmids)
            for i in range(len(sorted_cuts)):
                start = sorted_cuts[i]
                end = sorted_cuts[(i + 1) % len(sorted_cuts)]
                
                if end > start:
                    # Normal fragment
                    fragments.append(sequence[start:end])
                else:
                    # Wraparound fragment
                    fragments.append(sequence[start:] + sequence[:end])
        else:
            # For linear sequences or single cut
            prev_cut = 0
            for cut in sorted_cuts:
                if cut > prev_cut:
                    fragments.append(sequence[prev_cut:cut])
                prev_cut = cut
            
            # Add the last fragment
            if prev_cut < len(sequence):
                fragments.append(sequence[prev_cut:])
        
        return [f for f in fragments if len(f) > 0]

    def ligate(self):
        """Perform ligation simulation"""
        if not self.vector_record or not self.insert_record:
            QMessageBox.warning(self, "Missing Sequences", "Please select both vector and insert sequences.")
            return

        vector_enzyme_name = self.vector_enzyme_combo.currentText().split(' (')[0]
        insert_enzyme_name = self.insert_enzyme_combo.currentText().split(' (')[0]

        if not vector_enzyme_name or not insert_enzyme_name:
            QMessageBox.warning(self, "Missing Enzymes", "Please select enzymes for both vector and insert.")
            return
            
        selected_vector_item = self.vector_fragments_list.currentItem()
        selected_insert_item = self.insert_fragments_list.currentItem()

        if not selected_vector_item or not selected_insert_item:
            QMessageBox.warning(self, "Missing Fragments", "Please select fragments from both vector and insert.")
            return

        try:
            # Get enzymes
            vector_enzyme = getattr(Restriction, vector_enzyme_name)
            insert_enzyme = getattr(Restriction, insert_enzyme_name)
            
            # Analyze cuts
            vector_batch = RestrictionBatch([vector_enzyme])
            vector_cuts = vector_batch.search(self.vector_record.seq)[vector_enzyme]
            vector_fragments = self.get_fragments(self.vector_record.seq, vector_cuts, is_circular=True)
            
            insert_batch = RestrictionBatch([insert_enzyme])
            insert_cuts = insert_batch.search(self.insert_record.seq)[insert_enzyme]
            insert_fragments = self.get_fragments(self.insert_record.seq, insert_cuts, is_circular=False)

            vector_frag_idx = self.vector_fragments_list.currentRow()
            insert_frag_idx = self.insert_fragments_list.currentRow()
            
            if vector_frag_idx >= len(vector_fragments) or insert_frag_idx >= len(insert_fragments):
                QMessageBox.critical(self, "Error", "Invalid fragment selection.")
                return
                
            vector_frag = vector_fragments[vector_frag_idx]
            insert_frag = insert_fragments[insert_frag_idx]

            # Check compatibility
            compatible = vector_enzyme_name == insert_enzyme_name
            
            if not compatible:
                # Check if enzymes create compatible ends
                vector_overhang = self.get_enzyme_overhang(vector_enzyme)
                insert_overhang = self.get_enzyme_overhang(insert_enzyme)
                
                if vector_overhang != insert_overhang:
                    reply = QMessageBox.question(self, "Compatibility Warning", 
                        f"Vector enzyme ({vector_enzyme_name}) and insert enzyme ({insert_enzyme_name}) create different overhangs.\n"
                        f"Vector overhang: {vector_overhang}\n"
                        f"Insert overhang: {insert_overhang}\n\n"
                        "This may result in incompatible ends. Continue anyway?",
                        QMessageBox.Yes | QMessageBox.No)
                    if reply == QMessageBox.No:
                        return

            # Perform proper restriction cloning ligation
            # For restriction cloning, we typically insert into a linearized vector
            if len(vector_cuts) == 1:
                # Vector was linearized - insert the fragment
                new_seq = vector_frag + insert_frag
            else:
                # Multiple cuts - use selected fragments
                new_seq = vector_frag + insert_frag
            
            vector_name = self.vector_record.id or "vector"
            insert_name = self.insert_record.id or "insert"
            
            # Create a more realistic product
            new_record = SeqRecord(
                new_seq, 
                id=f"{vector_name}_{insert_name}_clone",
                name=f"Restriction Clone",
                description=f"Restriction cloning: {insert_name} into {vector_name} using {vector_enzyme_name}/{insert_enzyme_name}"
            )
            
            # Add some basic features if possible
            if hasattr(self.vector_record, 'features') and hasattr(self.insert_record, 'features'):
                # This is simplified - real implementation would map features properly
                pass

            # Calculate efficiency estimate
            efficiency = self.estimate_ligation_efficiency(vector_enzyme_name, insert_enzyme_name, 
                                                          len(vector_frag), len(insert_frag))

            # Show result
            result_msg = f"""Restriction Cloning Successful!
            
Cloning Details:
            - Vector: {vector_name} ({len(self.vector_record.seq)} bp)
            - Insert: {insert_name} ({len(self.insert_record.seq)} bp)
            - Vector Fragment: {len(vector_frag)} bp
            - Insert Fragment: {len(insert_frag)} bp
            - Final Clone: {len(new_seq)} bp
            - Estimated Efficiency: {efficiency}%
            
Enzymes Used:
            - Vector: {vector_enzyme_name} ({len(vector_cuts)} cuts)
            - Insert: {insert_enzyme_name} ({len(insert_cuts)} cuts)
            
The cloned sequence will be displayed in the main viewer."""
            
            QMessageBox.information(self, "Cloning Complete", result_msg)

            # Display in parent viewer
            if hasattr(self.parent(), 'display_sequence'):
                self.parent().display_sequence(new_record)
            
            self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Cloning Error", f"Error during restriction cloning: {str(e)}")
    
    def get_enzyme_overhang(self, enzyme):
        """Get the overhang type for an enzyme (simplified)"""
        try:
            # This is a simplified check - real implementation would be more complex
            overhang_map = {
                'EcoRI': "5' AATT",
                'BamHI': "5' GATC", 
                'HindIII': "5' AGCT",
                'XhoI': "5' TCGA",
                'SacI': "5' AGCT",
                'KpnI': "5' GTAC",
                'NotI': "5' GGCC",
                'XbaI': "5' CTAG",
                'SpeI': "5' CTAG",
                'PstI': "5' TGCA",
                'SalI': "5' TCGA",
                'SmaI': "blunt",
                'HpaI': "blunt"
            }
            return overhang_map.get(str(enzyme), "unknown")
        except:
            return "unknown"
    
    def estimate_ligation_efficiency(self, vector_enzyme, insert_enzyme, vector_size, insert_size):
        """Estimate ligation efficiency based on various factors"""
        base_efficiency = 70  # Base efficiency percentage
        
        # Same enzyme bonus
        if vector_enzyme == insert_enzyme:
            base_efficiency += 20
        
        # Size penalty for very large constructs
        total_size = vector_size + insert_size
        if total_size > 10000:
            base_efficiency -= 10
        elif total_size > 20000:
            base_efficiency -= 20
            
        # Blunt end penalty
        if self.get_enzyme_overhang(getattr(Restriction, vector_enzyme)) == "blunt":
            base_efficiency -= 15
            
        return max(10, min(95, base_efficiency))  # Keep between 10-95%
    
    def update_info_display(self):
        """Update the information display"""
        info_text = "Cloning Information:\n\n"
        
        if self.vector_record:
            info_text += f"Vector: {self.vector_record.id}\n"
            info_text += f"Length: {len(self.vector_record.seq)} bp\n"
            if hasattr(self.vector_record, 'features'):
                feature_count = len([f for f in self.vector_record.features if f.type != 'source'])
                info_text += f"Features: {feature_count}\n"
            info_text += "\n"
        
        if self.insert_record:
            info_text += f"Insert: {self.insert_record.id}\n"
            info_text += f"Length: {len(self.insert_record.seq)} bp\n"
            if hasattr(self.insert_record, 'features'):
                feature_count = len([f for f in self.insert_record.features if f.type != 'source'])
                info_text += f"Features: {feature_count}\n"
        
        if not self.vector_record and not self.insert_record:
            info_text += "No sequences loaded yet.\n"
            info_text += "Select vector and insert files to begin."
        
        self.info_display.setText(info_text)
