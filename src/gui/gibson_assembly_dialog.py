from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QListWidget, QPushButton, QFileDialog,
                             QHBoxLayout, QLabel, QTextEdit, QGroupBox, QMessageBox,
                             QSplitter, QSpinBox, QFormLayout, QCheckBox, QWidget)
from PyQt5.QtCore import Qt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

class GibsonAssemblyDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Gibson Assembly Simulation")
        self.setMinimumSize(700, 500)
        self.fragments = []

        # Main layout
        main_layout = QVBoxLayout(self)
        
        # Instructions
        instructions = QLabel("""
        <b>Gibson Assembly Simulation</b><br>
        Gibson Assembly allows seamless joining of DNA fragments with overlapping ends.<br>
        1. Add DNA fragments in the correct order<br>
        2. Specify overlap length (default: 20 bp)<br>
        3. Simulate the assembly reaction
        """)
        instructions.setWordWrap(True)
        instructions.setStyleSheet("QLabel { background-color: #f0f0f0; padding: 10px; border-radius: 5px; }")
        main_layout.addWidget(instructions)
        
        # Create splitter for main content
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left panel - Fragment management
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # Fragment list group
        fragment_group = QGroupBox("DNA Fragments")
        fragment_layout = QVBoxLayout(fragment_group)
        
        self.fragment_list = QListWidget()
        self.fragment_list.setToolTip("List of DNA fragments for assembly")
        fragment_layout.addWidget(self.fragment_list)
        
        # Fragment control buttons
        button_layout = QHBoxLayout()
        
        self.add_fragment_button = QPushButton("📁 Add Fragment")
        self.add_fragment_button.setToolTip("Add a DNA fragment from file")
        self.add_fragment_button.clicked.connect(self.add_fragment)
        button_layout.addWidget(self.add_fragment_button)
        
        self.remove_fragment_button = QPushButton("🗑️ Remove")
        self.remove_fragment_button.setToolTip("Remove selected fragment")
        self.remove_fragment_button.clicked.connect(self.remove_fragment)
        self.remove_fragment_button.setEnabled(False)
        button_layout.addWidget(self.remove_fragment_button)
        
        self.move_up_button = QPushButton("⬆️ Up")
        self.move_up_button.setToolTip("Move fragment up in order")
        self.move_up_button.clicked.connect(self.move_fragment_up)
        self.move_up_button.setEnabled(False)
        button_layout.addWidget(self.move_up_button)
        
        self.move_down_button = QPushButton("⬇️ Down")
        self.move_down_button.setToolTip("Move fragment down in order")
        self.move_down_button.clicked.connect(self.move_fragment_down)
        self.move_down_button.setEnabled(False)
        button_layout.addWidget(self.move_down_button)
        
        fragment_layout.addLayout(button_layout)
        left_layout.addWidget(fragment_group)
        
        # Assembly parameters
        params_group = QGroupBox("Assembly Parameters")
        params_layout = QFormLayout(params_group)
        
        self.overlap_length = QSpinBox()
        self.overlap_length.setRange(10, 100)
        self.overlap_length.setValue(20)
        self.overlap_length.setSuffix(" bp")
        self.overlap_length.setToolTip("Length of overlapping sequences between fragments")
        params_layout.addRow("Overlap Length:", self.overlap_length)
        
        self.circular_assembly = QCheckBox("Circular Assembly")
        self.circular_assembly.setToolTip("Create a circular product (plasmid)")
        self.circular_assembly.setChecked(True)
        params_layout.addRow(self.circular_assembly)
        
        left_layout.addWidget(params_group)
        
        # Assembly button
        self.assemble_button = QPushButton("🧬 Perform Gibson Assembly")
        self.assemble_button.setToolTip("Simulate Gibson assembly of fragments")
        self.assemble_button.clicked.connect(self.assemble)
        self.assemble_button.setEnabled(False)
        left_layout.addWidget(self.assemble_button)
        
        splitter.addWidget(left_panel)
        
        # Right panel - Information and results
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        # Fragment information
        info_group = QGroupBox("Fragment Information")
        info_layout = QVBoxLayout(info_group)
        
        self.info_display = QTextEdit()
        self.info_display.setReadOnly(True)
        self.info_display.setMaximumHeight(200)
        self.info_display.setText("Add DNA fragments to begin Gibson assembly.")
        info_layout.addWidget(self.info_display)
        
        right_layout.addWidget(info_group)
        
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
        splitter.setSizes([350, 350])
        
        # Button box
        button_box_layout = QHBoxLayout()
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_box_layout.addStretch()
        button_box_layout.addWidget(close_button)
        
        main_layout.addLayout(button_box_layout)
        
        # Connect selection change
        self.fragment_list.itemSelectionChanged.connect(self.update_button_states)

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
            
            # Enable assembly button if we have at least 2 fragments
            if len(self.fragments) >= 2:
                self.assemble_button.setEnabled(True)
            
            self.update_info_display()
            self.update_assembly_preview()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading fragment: {str(e)}")
    
    def remove_fragment(self):
        """Remove selected fragment"""
        current_row = self.fragment_list.currentRow()
        if current_row >= 0:
            self.fragments.pop(current_row)
            self.fragment_list.takeItem(current_row)
            
            # Disable assembly button if less than 2 fragments
            if len(self.fragments) < 2:
                self.assemble_button.setEnabled(False)
            
            self.update_info_display()
            self.update_assembly_preview()
            self.update_button_states()
    
    def move_fragment_up(self):
        """Move selected fragment up in order"""
        current_row = self.fragment_list.currentRow()
        if current_row > 0:
            # Swap fragments
            self.fragments[current_row], self.fragments[current_row - 1] = \
                self.fragments[current_row - 1], self.fragments[current_row]
            
            # Update list
            self.update_fragment_list()
            self.fragment_list.setCurrentRow(current_row - 1)
            
            self.update_assembly_preview()
    
    def move_fragment_down(self):
        """Move selected fragment down in order"""
        current_row = self.fragment_list.currentRow()
        if current_row < len(self.fragments) - 1:
            # Swap fragments
            self.fragments[current_row], self.fragments[current_row + 1] = \
                self.fragments[current_row + 1], self.fragments[current_row]
            
            # Update list
            self.update_fragment_list()
            self.fragment_list.setCurrentRow(current_row + 1)
            
            self.update_assembly_preview()
    
    def update_fragment_list(self):
        """Update the fragment list display"""
        self.fragment_list.clear()
        for fragment in self.fragments:
            fragment_name = fragment.id or "Unknown"
            self.fragment_list.addItem(f"{fragment_name} ({len(fragment.seq)} bp)")
    
    def update_button_states(self):
        """Update button enabled states based on selection"""
        current_row = self.fragment_list.currentRow()
        has_selection = current_row >= 0
        
        self.remove_fragment_button.setEnabled(has_selection)
        self.move_up_button.setEnabled(has_selection and current_row > 0)
        self.move_down_button.setEnabled(has_selection and current_row < len(self.fragments) - 1)
    
    def update_info_display(self):
        """Update fragment information display"""
        if not self.fragments:
            self.info_display.setText("No fragments added yet.\nAdd at least 2 fragments for Gibson assembly.")
            return
        
        info_text = f"Gibson Assembly Information:\n\n"
        info_text += f"Total fragments: {len(self.fragments)}\n\n"
        
        total_length = 0
        for i, fragment in enumerate(self.fragments, 1):
            fragment_name = fragment.id or f"Fragment {i}"
            length = len(fragment.seq)
            total_length += length
            info_text += f"{i}. {fragment_name}: {length} bp\n"
        
        info_text += f"\nTotal length: {total_length} bp\n"
        
        if len(self.fragments) >= 2:
            overlap = self.overlap_length.value()
            expected_length = total_length - (len(self.fragments) - 1) * overlap
            if self.circular_assembly.isChecked():
                expected_length -= overlap  # Additional overlap for circularization
            info_text += f"Expected product: {expected_length} bp"
        
        self.info_display.setText(info_text)
    
    def update_assembly_preview(self):
        """Update assembly preview"""
        if len(self.fragments) < 2:
            self.preview_display.setText("Add at least 2 fragments to preview assembly.")
            return
        
        preview_text = "Gibson Assembly Preview:\n\n"
        overlap = self.overlap_length.value()
        
        # Check for potential issues
        issues = []
        min_length = min(len(f.seq) for f in self.fragments)
        if overlap >= min_length / 2:
            issues.append(f"⚠️ Large overlap vs. smallest fragment ({min_length} bp)")
        
        if len(self.fragments) > 6:
            issues.append("⚠️ Many fragments may reduce efficiency")
        
        total_original = sum(len(f.seq) for f in self.fragments)
        total_overlaps = (len(self.fragments) - 1) * overlap
        if self.circular_assembly.isChecked():
            total_overlaps += overlap
        
        expected_length = total_original - total_overlaps
        
        preview_text += f"Expected product: {expected_length} bp\n"
        preview_text += f"Efficiency estimate: {self.estimate_gibson_efficiency(len(self.fragments), overlap, expected_length)}%\n\n"
        
        if issues:
            preview_text += "Issues detected:\n"
            for issue in issues:
                preview_text += f"{issue}\n"
            preview_text += "\n"
        
        preview_text += "Assembly order:\n"
        for i, fragment in enumerate(self.fragments):
            fragment_name = fragment.id or f"Fragment {i+1}"
            preview_text += f"{i+1}. {fragment_name} ({len(fragment.seq)} bp)\n"
            
            if i < len(self.fragments) - 1:
                preview_text += f"   ↓ {overlap} bp overlap removed\n"
        
        if self.circular_assembly.isChecked():
            preview_text += f"   ↓ {overlap} bp overlap (circularization)\n"
            preview_text += f"→ Circular product\n"
        
        preview_text += "\nGibson Assembly Requirements:\n"
        preview_text += f"• Primers with {overlap} bp overlaps\n"
        preview_text += "• 5' exonuclease, DNA polymerase, ligase\n"
        preview_text += "• Optimized temperature and timing"
        
        self.preview_display.setText(preview_text)

    def assemble(self):
        """Perform Gibson assembly simulation"""
        if len(self.fragments) < 2:
            QMessageBox.warning(self, "Insufficient Fragments", "At least 2 fragments are required for Gibson assembly.")
            return

        try:
            overlap = self.overlap_length.value()
            
            # Validate overlap length
            min_fragment_length = min(len(f.seq) for f in self.fragments)
            if overlap >= min_fragment_length / 2:
                reply = QMessageBox.question(self, "Large Overlap Warning", 
                    f"Overlap length ({overlap} bp) is large compared to smallest fragment ({min_fragment_length} bp).\n"
                    "This may result in significant sequence loss. Continue?",
                    QMessageBox.Yes | QMessageBox.No)
                if reply == QMessageBox.No:
                    return
            
            # More realistic Gibson Assembly simulation
            assembled_fragments = []
            total_overlap_removed = 0
            
            for i, fragment in enumerate(self.fragments):
                fragment_seq = fragment.seq
                
                if i > 0:
                    # Remove overlap from beginning of fragment (except first)
                    if len(fragment_seq) > overlap:
                        fragment_seq = fragment_seq[overlap:]
                        total_overlap_removed += overlap
                    else:
                        # Fragment too short - use half
                        overlap_to_remove = len(fragment_seq) // 2
                        fragment_seq = fragment_seq[overlap_to_remove:]
                        total_overlap_removed += overlap_to_remove
                
                assembled_fragments.append(fragment_seq)
            
            # Join all fragments
            new_seq = Seq("")
            for frag_seq in assembled_fragments:
                new_seq += frag_seq
            
            # For circular assembly, check if ends can be joined
            if self.circular_assembly.isChecked():
                # Check if first and last fragments have compatible ends
                first_fragment = self.fragments[0].seq
                last_fragment = self.fragments[-1].seq
                
                # Remove final overlap for circularization
                if len(new_seq) > overlap:
                    new_seq = new_seq[:-overlap]
                    total_overlap_removed += overlap
            
            # Validate result
            if len(new_seq) < 100:
                QMessageBox.warning(self, "Assembly Warning", 
                    f"Assembled sequence is very short ({len(new_seq)} bp).\n"
                    "Check overlap settings and fragment compatibility.")
            
            # Create product record
            fragment_names = [f.id or f"Fragment_{i+1}" for i, f in enumerate(self.fragments)]
            product_id = "_".join(fragment_names[:3])  # Use first 3 names
            if len(fragment_names) > 3:
                product_id += f"_and_{len(fragment_names)-3}_more"
            
            topology = "circular" if self.circular_assembly.isChecked() else "linear"
            
            new_record = SeqRecord(
                new_seq,
                id=f"{product_id}_gibson",
                name="Gibson Assembly Product",
                description=f"Gibson assembly of {len(self.fragments)} fragments ({topology}), {overlap}bp overlaps"
            )
            
            # Calculate efficiency estimate
            efficiency = self.estimate_gibson_efficiency(len(self.fragments), overlap, len(new_seq))
            
            # Show result
            original_length = sum(len(f.seq) for f in self.fragments)
            result_msg = f"""Gibson Assembly Successful!
            
Assembly Details:
            - Fragments assembled: {len(self.fragments)}
            - Topology: {topology.title()}
            - Overlap length: {overlap} bp
            - Original total length: {original_length} bp
            - Overlap removed: {total_overlap_removed} bp
            - Final product length: {len(new_seq)} bp
            - Estimated efficiency: {efficiency}%
            
Gibson Assembly Process:
            ✓ 5' exonuclease creates overhangs
            ✓ DNA polymerase fills gaps
            ✓ DNA ligase seals nicks
            
Note: Real Gibson assembly requires:
            - Primers designed with {overlap}bp overlaps
            - Optimized reaction conditions
            - Quality control of fragments
            
The assembled sequence will be displayed in the main viewer."""
            
            QMessageBox.information(self, "Gibson Assembly Complete", result_msg)

            # Display in parent viewer
            if hasattr(self.parent(), 'display_sequence'):
                self.parent().display_sequence(new_record)
            
            self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Assembly Error", f"Error during Gibson assembly: {str(e)}")
    
    def estimate_gibson_efficiency(self, fragment_count, overlap_length, product_length):
        """Estimate Gibson assembly efficiency"""
        base_efficiency = 80  # Base efficiency for Gibson assembly
        
        # Fragment count penalty
        if fragment_count > 4:
            base_efficiency -= (fragment_count - 4) * 5
        
        # Overlap length optimization
        if 15 <= overlap_length <= 25:
            base_efficiency += 10  # Optimal overlap range
        elif overlap_length < 15:
            base_efficiency -= 15  # Too short
        elif overlap_length > 40:
            base_efficiency -= 10  # Too long
        
        # Product size penalty
        if product_length > 15000:
            base_efficiency -= 10
        elif product_length > 25000:
            base_efficiency -= 20
        
        return max(20, min(95, base_efficiency))
