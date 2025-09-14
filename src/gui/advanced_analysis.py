from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
                             QComboBox, QSpinBox, QCheckBox, QGroupBox, QFormLayout,
                             QTabWidget, QWidget, QTextEdit, QTableWidget, QTableWidgetItem,
                             QHeaderView, QProgressDialog, QMessageBox, QSlider,
                             QSplitter, QListWidget, QListWidgetItem, QLineEdit)
from PyQt5.QtGui import QFont, QColor, QPainter, QPen, QBrush
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
try:
    from Bio.Restriction import RestrictionBatch, Analysis
    RESTRICTION_AVAILABLE = True
except ImportError:
    RESTRICTION_AVAILABLE = False

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    import numpy as np
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

import re

class AnalysisWorker(QThread):
    """Worker thread for analysis operations"""
    progress = pyqtSignal(int)
    result = pyqtSignal(object)
    error = pyqtSignal(str)
    
    def __init__(self, analysis_func, *args, **kwargs):
        super().__init__()
        self.analysis_func = analysis_func
        self.args = args
        self.kwargs = kwargs
    
    def run(self):
        try:
            result = self.analysis_func(*self.args, **self.kwargs)
            self.result.emit(result)
        except Exception as e:
            self.error.emit(str(e))

class ORFAnalysisWidget(QWidget):
    """Widget for Open Reading Frame analysis"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.record = None
        self.setup_ui()
    
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Parameters
        params_group = QGroupBox("ORF Parameters")
        params_layout = QFormLayout(params_group)
        
        self.min_length_spin = QSpinBox()
        self.min_length_spin.setRange(30, 3000)
        self.min_length_spin.setValue(100)
        self.min_length_spin.setSuffix(" bp")
        params_layout.addRow("Minimum length:", self.min_length_spin)
        
        self.start_codons_edit = QLineEdit("ATG")
        params_layout.addRow("Start codons:", self.start_codons_edit)
        
        self.stop_codons_edit = QLineEdit("TAA,TAG,TGA")
        params_layout.addRow("Stop codons:", self.stop_codons_edit)
        
        self.all_frames_check = QCheckBox("Analyze all 6 reading frames")
        self.all_frames_check.setChecked(True)
        params_layout.addRow("", self.all_frames_check)
        
        analyze_btn = QPushButton("Find ORFs")
        analyze_btn.clicked.connect(self.find_orfs)
        params_layout.addRow("", analyze_btn)
        
        layout.addWidget(params_group)
        
        # Results
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels([
            "Frame", "Start", "End", "Length (bp)", "Length (aa)", "Sequence"
        ])
        self.results_table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.results_table)
        
        # Statistics
        self.stats_text = QTextEdit()
        self.stats_text.setMaximumHeight(100)
        self.stats_text.setReadOnly(True)
        layout.addWidget(self.stats_text)
    
    def set_sequence(self, record):
        """Set the sequence for analysis"""
        self.record = record
        self.results_table.setRowCount(0)
        self.stats_text.clear()
    
    def find_orfs(self):
        """Find Open Reading Frames"""
        if not self.record:
            QMessageBox.warning(self, "No Sequence", "No sequence loaded for analysis.")
            return
        
        sequence = str(self.record.seq).upper()
        min_length = self.min_length_spin.value()
        start_codons = [c.strip().upper() for c in self.start_codons_edit.text().split(',')]
        stop_codons = [c.strip().upper() for c in self.stop_codons_edit.text().split(',')]
        
        orfs = []
        
        # Analyze reading frames
        frames_to_analyze = 6 if self.all_frames_check.isChecked() else 3
        
        for frame in range(frames_to_analyze):
            if frame < 3:
                # Forward frames
                seq = sequence[frame:]
                strand = 1
                frame_name = f"+{frame + 1}"
            else:
                # Reverse frames
                rev_seq = str(Seq(sequence).reverse_complement())
                seq = rev_seq[frame - 3:]
                strand = -1
                frame_name = f"-{frame - 2}"
            
            # Find ORFs in this frame
            frame_orfs = self.find_orfs_in_frame(seq, frame, strand, min_length, start_codons, stop_codons)
            orfs.extend([(frame_name, orf) for orf in frame_orfs])
        
        # Display results
        self.display_orfs(orfs)
        self.display_statistics(orfs)
    
    def find_orfs_in_frame(self, sequence, frame, strand, min_length, start_codons, stop_codons):
        """Find ORFs in a specific reading frame"""
        orfs = []
        
        i = 0
        while i < len(sequence) - 2:
            codon = sequence[i:i+3]
            
            if len(codon) == 3 and codon in start_codons:
                # Found start codon, look for stop
                start_pos = i
                j = i + 3
                
                while j < len(sequence) - 2:
                    stop_codon = sequence[j:j+3]
                    if len(stop_codon) == 3 and stop_codon in stop_codons:
                        # Found stop codon
                        end_pos = j + 3
                        orf_length = end_pos - start_pos
                        
                        if orf_length >= min_length:
                            # Calculate actual positions in original sequence
                            if strand == 1:
                                actual_start = frame + start_pos
                                actual_end = frame + end_pos
                            else:
                                # For reverse strand, calculate from end
                                seq_len = len(str(self.record.seq))
                                actual_start = seq_len - (frame + end_pos)
                                actual_end = seq_len - (frame + start_pos)
                            
                            orf_seq = sequence[start_pos:end_pos]
                            orfs.append({
                                'start': actual_start + 1,  # 1-based
                                'end': actual_end,
                                'length_bp': orf_length,
                                'length_aa': orf_length // 3,
                                'sequence': orf_seq,
                                'strand': strand
                            })
                        
                        i = j + 3
                        break
                    j += 3
                else:
                    # No stop codon found
                    i += 3
            else:
                i += 3
        
        return orfs
    
    def display_orfs(self, orfs):
        """Display ORFs in the table"""
        self.results_table.setRowCount(len(orfs))
        
        for row, (frame, orf) in enumerate(orfs):
            self.results_table.setItem(row, 0, QTableWidgetItem(frame))
            self.results_table.setItem(row, 1, QTableWidgetItem(str(orf['start'])))
            self.results_table.setItem(row, 2, QTableWidgetItem(str(orf['end'])))
            self.results_table.setItem(row, 3, QTableWidgetItem(str(orf['length_bp'])))
            self.results_table.setItem(row, 4, QTableWidgetItem(str(orf['length_aa'])))
            self.results_table.setItem(row, 5, QTableWidgetItem(orf['sequence'][:50] + "..." if len(orf['sequence']) > 50 else orf['sequence']))
    
    def display_statistics(self, orfs):
        """Display ORF statistics"""
        if not orfs:
            self.stats_text.setText("No ORFs found with the specified criteria.")
            return
        
        total_orfs = len(orfs)
        lengths = [orf[1]['length_bp'] for orf in orfs]
        avg_length = sum(lengths) / len(lengths)
        max_length = max(lengths)
        min_length = min(lengths)
        
        forward_orfs = len([orf for orf in orfs if orf[1]['strand'] == 1])
        reverse_orfs = len([orf for orf in orfs if orf[1]['strand'] == -1])
        
        stats = f"""ORF Analysis Results:
        Total ORFs found: {total_orfs}
        Forward strand: {forward_orfs}
        Reverse strand: {reverse_orfs}
        Average length: {avg_length:.1f} bp
        Longest ORF: {max_length} bp
        Shortest ORF: {min_length} bp"""
        
        self.stats_text.setText(stats)

class PrimerDesignWidget(QWidget):
    """Widget for primer design"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.record = None
        self.setup_ui()
    
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Target region
        target_group = QGroupBox("Target Region")
        target_layout = QFormLayout(target_group)
        
        self.start_spin = QSpinBox()
        self.start_spin.setRange(1, 10000)
        self.start_spin.setValue(100)
        target_layout.addRow("Start position:", self.start_spin)
        
        self.end_spin = QSpinBox()
        self.end_spin.setRange(1, 10000)
        self.end_spin.setValue(500)
        target_layout.addRow("End position:", self.end_spin)
        
        # Add amplicon length display
        self.amplicon_label = QLabel("Amplicon length: 400 bp")
        target_layout.addRow("", self.amplicon_label)
        
        # Connect signals to update amplicon length
        self.start_spin.valueChanged.connect(self.update_amplicon_length)
        self.end_spin.valueChanged.connect(self.update_amplicon_length)
        
        layout.addWidget(target_group)
        
        # Primer parameters
        params_group = QGroupBox("Primer Parameters")
        params_layout = QFormLayout(params_group)
        
        self.primer_length_spin = QSpinBox()
        self.primer_length_spin.setRange(15, 35)
        self.primer_length_spin.setValue(20)
        params_layout.addRow("Target primer length:", self.primer_length_spin)
        
        self.length_tolerance_spin = QSpinBox()
        self.length_tolerance_spin.setRange(1, 10)
        self.length_tolerance_spin.setValue(3)
        params_layout.addRow("Length tolerance (±):", self.length_tolerance_spin)
        
        self.tm_min_spin = QSpinBox()
        self.tm_min_spin.setRange(45, 75)
        self.tm_min_spin.setValue(55)
        self.tm_min_spin.setSuffix("°C")
        params_layout.addRow("Min Tm:", self.tm_min_spin)
        
        self.tm_max_spin = QSpinBox()
        self.tm_max_spin.setRange(50, 80)
        self.tm_max_spin.setValue(65)
        self.tm_max_spin.setSuffix("°C")
        params_layout.addRow("Max Tm:", self.tm_max_spin)
        
        self.gc_min_spin = QSpinBox()
        self.gc_min_spin.setRange(20, 80)
        self.gc_min_spin.setValue(40)
        self.gc_min_spin.setSuffix("%")
        params_layout.addRow("Min GC%:", self.gc_min_spin)
        
        self.gc_max_spin = QSpinBox()
        self.gc_max_spin.setRange(20, 80)
        self.gc_max_spin.setValue(60)
        self.gc_max_spin.setSuffix("%")
        params_layout.addRow("Max GC%:", self.gc_max_spin)
        
        # Advanced options
        self.avoid_runs_check = QCheckBox("Avoid nucleotide runs (4+ same bases)")
        self.avoid_runs_check.setChecked(True)
        params_layout.addRow("", self.avoid_runs_check)
        
        self.gc_clamp_check = QCheckBox("Prefer GC clamp at 3' end")
        self.gc_clamp_check.setChecked(True)
        params_layout.addRow("", self.gc_clamp_check)
        
        design_btn = QPushButton("Design Primers")
        design_btn.clicked.connect(self.design_primers)
        design_btn.setStyleSheet("QPushButton { background-color: #4CAF50; color: white; font-weight: bold; padding: 8px; }")
        params_layout.addRow("", design_btn)
        
        layout.addWidget(params_group)
        
        # Results
        results_group = QGroupBox("Primer Results")
        results_layout = QVBoxLayout(results_group)
        
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(7)
        self.results_table.setHorizontalHeaderLabels([
            "Type", "Sequence (5' → 3')", "Length", "Tm (°C)", "GC%", "Position", "Score"
        ])
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.setAlternatingRowColors(True)
        self.results_table.setSelectionBehavior(QTableWidget.SelectRows)
        results_layout.addWidget(self.results_table)
        
        # Add primer pair analysis
        self.pair_analysis = QTextEdit()
        self.pair_analysis.setMaximumHeight(120)
        self.pair_analysis.setReadOnly(True)
        self.pair_analysis.setPlaceholderText("Primer pair analysis will appear here...")
        results_layout.addWidget(QLabel("Best Primer Pair Analysis:"))
        results_layout.addWidget(self.pair_analysis)
        
        # Add export buttons
        export_layout = QHBoxLayout()
        
        export_primers_btn = QPushButton("📋 Copy Primers")
        export_primers_btn.setToolTip("Copy primer sequences to clipboard")
        export_primers_btn.clicked.connect(self.copy_primers_to_clipboard)
        export_layout.addWidget(export_primers_btn)
        
        export_csv_btn = QPushButton("💾 Export CSV")
        export_csv_btn.setToolTip("Export primer data to CSV file")
        export_csv_btn.clicked.connect(self.export_primers_csv)
        export_layout.addWidget(export_csv_btn)
        
        export_layout.addStretch()
        
        # Add primer validation button
        validate_btn = QPushButton("🔍 Validate Selected")
        validate_btn.setToolTip("Detailed validation of selected primer")
        validate_btn.clicked.connect(self.validate_selected_primer)
        export_layout.addWidget(validate_btn)
        
        results_layout.addLayout(export_layout)
        
        layout.addWidget(results_group)
    
    def set_sequence(self, record):
        """Set the sequence for primer design"""
        self.record = record
        if record:
            seq_len = len(record.seq)
            self.start_spin.setRange(1, seq_len)
            self.end_spin.setRange(1, seq_len)
            
            # Set reasonable default values
            default_start = min(100, seq_len // 4)
            default_end = min(seq_len - 100, seq_len * 3 // 4)
            
            self.start_spin.setValue(default_start)
            self.end_spin.setValue(default_end)
            
            self.update_amplicon_length()
        self.results_table.setRowCount(0)
        self.pair_analysis.clear()
    
    def update_amplicon_length(self):
        """Update the amplicon length display"""
        if self.record:
            start = self.start_spin.value()
            end = self.end_spin.value()
            length = max(0, end - start + 1)
            self.amplicon_label.setText(f"Amplicon length: {length} bp")
    
    def design_primers(self):
        """Design primers for the specified region"""
        if not self.record:
            QMessageBox.warning(self, "No Sequence", "No sequence loaded for primer design.")
            return
        
        sequence = str(self.record.seq).upper()
        start = self.start_spin.value() - 1  # Convert to 0-based
        end = self.end_spin.value() - 1  # Convert to 0-based
        
        if start >= end:
            QMessageBox.warning(self, "Invalid Region", "Start position must be less than end position.")
            return
        
        if end - start < 50:
            QMessageBox.warning(self, "Region Too Small", "Target region should be at least 50 bp for meaningful primer design.")
            return
        
        # Get parameters
        target_length = self.primer_length_spin.value()
        length_tolerance = self.length_tolerance_spin.value()
        tm_min = self.tm_min_spin.value()
        tm_max = self.tm_max_spin.value()
        gc_min = self.gc_min_spin.value()
        gc_max = self.gc_max_spin.value()
        avoid_runs = self.avoid_runs_check.isChecked()
        prefer_gc_clamp = self.gc_clamp_check.isChecked()
        
        try:
            # Design forward primers (around start position)
            forward_primers = self.find_primers(
                sequence, start, target_length, length_tolerance, 
                tm_min, tm_max, gc_min, gc_max, avoid_runs, prefer_gc_clamp, forward=True
            )
            
            # Design reverse primers (around end position)
            reverse_primers = self.find_primers(
                sequence, end, target_length, length_tolerance,
                tm_min, tm_max, gc_min, gc_max, avoid_runs, prefer_gc_clamp, forward=False
            )
            
            # Combine and display results
            all_primers = [("Forward", p) for p in forward_primers] + [("Reverse", p) for p in reverse_primers]
            
            if not all_primers:
                QMessageBox.information(self, "No Primers Found", 
                    "No primers found matching the specified criteria. Try relaxing the parameters.")
                return
            
            self.display_primers(all_primers)
            
            # Analyze primer pairs
            if forward_primers and reverse_primers:
                self.analyze_primer_pairs(forward_primers, reverse_primers, start, end)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error during primer design: {str(e)}")
    
    def find_primers(self, sequence, position, target_length, length_tolerance, 
                    tm_min, tm_max, gc_min, gc_max, avoid_runs, prefer_gc_clamp, forward=True):
        """Find suitable primers at a given position with enhanced scoring"""
        primers = []
        
        # Define search range around the target position
        search_range = 50  # Look 50 bp in each direction
        
        if forward:
            # For forward primers, search around the start position
            search_start = max(0, position - search_range)
            search_end = min(len(sequence), position + search_range)
        else:
            # For reverse primers, search around the end position
            search_start = max(0, position - search_range)
            search_end = min(len(sequence), position + search_range)
        
        # Try different lengths
        min_length = max(15, target_length - length_tolerance)
        max_length = min(35, target_length + length_tolerance)
        
        for start_pos in range(search_start, search_end):
            for length in range(min_length, max_length + 1):
                if forward:
                    if start_pos + length <= len(sequence):
                        primer_seq = sequence[start_pos:start_pos + length]
                        primer_pos = start_pos + 1  # 1-based
                    else:
                        continue
                else:
                    # For reverse primers, we want the reverse complement
                    if start_pos + length <= len(sequence):
                        forward_seq = sequence[start_pos:start_pos + length]
                        primer_seq = str(Seq(forward_seq).reverse_complement())
                        primer_pos = start_pos + 1  # 1-based position of the forward sequence
                    else:
                        continue
                
                # Calculate properties
                tm = self.calculate_tm_improved(primer_seq)
                gc_content = gc_fraction(primer_seq) * 100
                
                # Check basic criteria
                if not (tm_min <= tm <= tm_max and gc_min <= gc_content <= gc_max):
                    continue
                
                # Check for problematic sequences
                if avoid_runs and self.has_problematic_sequences(primer_seq):
                    continue
                
                # Calculate primer score
                score = self.calculate_primer_score(
                    primer_seq, tm, gc_content, target_length, prefer_gc_clamp
                )
                
                primers.append({
                    'sequence': primer_seq,
                    'length': len(primer_seq),
                    'tm': tm,
                    'gc': gc_content,
                    'position': primer_pos,
                    'score': score,
                    'start_pos': start_pos  # Store for pair analysis
                })
        
        # Sort by score (higher is better)
        primers.sort(key=lambda p: p['score'], reverse=True)
        return primers[:10]  # Return top 10
    
    def calculate_tm_improved(self, sequence):
        """Calculate melting temperature using improved method"""
        # More accurate Tm calculation using salt-adjusted formula
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')
        length = len(sequence)
        
        if length < 14:
            # Wallace rule for short primers
            tm = (at_count * 2) + (gc_count * 4)
        else:
            # Salt-adjusted formula (assuming 50mM Na+, 1.5mM Mg2+)
            tm = 81.5 + 0.41 * (gc_count / length * 100) - (675 / length) - 0.65
        
        return round(tm, 1)
    
    def calculate_primer_score(self, sequence, tm, gc_content, target_length, prefer_gc_clamp):
        """Calculate a comprehensive primer score (0-100, higher is better)"""
        score = 50  # Base score
        
        # Tm score (optimal around 60°C)
        tm_diff = abs(tm - 60)
        if tm_diff <= 2:
            score += 20
        elif tm_diff <= 5:
            score += 15 - tm_diff
        else:
            score -= tm_diff
        
        # GC content score (optimal 40-60%)
        if 45 <= gc_content <= 55:
            score += 15
        elif 40 <= gc_content <= 60:
            score += 10
        else:
            score -= abs(gc_content - 50) / 2
        
        # Length score (prefer target length)
        length_diff = abs(len(sequence) - target_length)
        if length_diff == 0:
            score += 10
        elif length_diff <= 2:
            score += 5
        else:
            score -= length_diff
        
        # GC clamp bonus
        if prefer_gc_clamp:
            last_base = sequence[-1]
            if last_base in 'GC':
                score += 5
            # Prefer 1-2 GC bases in last 5 positions
            last_five = sequence[-5:]
            gc_in_last_five = last_five.count('G') + last_five.count('C')
            if 1 <= gc_in_last_five <= 3:
                score += 3
        
        # Penalize runs of same nucleotide
        for base in 'ATCG':
            if base * 4 in sequence:
                score -= 10
            elif base * 3 in sequence:
                score -= 5
        
        # Penalize strong secondary structures (simplified)
        if 'GGGG' in sequence or 'CCCC' in sequence:
            score -= 15
        
        # Penalize palindromes (simplified check)
        if len(sequence) >= 6:
            for i in range(len(sequence) - 5):
                substr = sequence[i:i+6]
                rev_comp = str(Seq(substr).reverse_complement())
                if substr == rev_comp:
                    score -= 10
                    break
        
        return max(0, min(100, score))
    
    def copy_primers_to_clipboard(self):
        """Copy primer sequences to clipboard"""
        try:
            from PyQt5.QtWidgets import QApplication
            
            clipboard_text = "Primer Design Results\n" + "="*50 + "\n\n"
            
            # Add primer table data
            for row in range(self.results_table.rowCount()):
                primer_type = self.results_table.item(row, 0).text()
                sequence = self.results_table.item(row, 1).text()
                length = self.results_table.item(row, 2).text()
                tm = self.results_table.item(row, 3).text()
                gc = self.results_table.item(row, 4).text()
                position = self.results_table.item(row, 5).text()
                score = self.results_table.item(row, 6).text()
                
                clipboard_text += f"{primer_type} Primer:\n"
                clipboard_text += f"Sequence: {sequence}\n"
                clipboard_text += f"Length: {length} bp, Tm: {tm}°C, GC: {gc}%, Score: {score}\n"
                clipboard_text += f"Position: {position}\n\n"
            
            # Add pair analysis if available
            pair_text = self.pair_analysis.toPlainText()
            if pair_text:
                clipboard_text += "\n" + pair_text
            
            QApplication.clipboard().setText(clipboard_text)
            QMessageBox.information(self, "Copied", "Primer data copied to clipboard!")
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not copy to clipboard: {str(e)}")
    
    def export_primers_csv(self):
        """Export primer data to CSV file"""
        try:
            from PyQt5.QtWidgets import QFileDialog
            import csv
            
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Export Primers", "primers.csv", "CSV Files (*.csv)"
            )
            
            if not file_path:
                return
            
            with open(file_path, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write header
                headers = ["Type", "Sequence", "Length", "Tm_C", "GC_percent", "Position", "Score"]
                writer.writerow(headers)
                
                # Write primer data
                for row in range(self.results_table.rowCount()):
                    row_data = []
                    for col in range(self.results_table.columnCount()):
                        item = self.results_table.item(row, col)
                        row_data.append(item.text() if item else "")
                    writer.writerow(row_data)
            
            QMessageBox.information(self, "Exported", f"Primer data exported to {file_path}")
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not export CSV: {str(e)}")
    
    def validate_selected_primer(self):
        """Show detailed validation for selected primer"""
        current_row = self.results_table.currentRow()
        if current_row < 0:
            QMessageBox.information(self, "No Selection", "Please select a primer to validate.")
            return
        
        try:
            # Get primer data
            primer_type = self.results_table.item(current_row, 0).text()
            sequence = self.results_table.item(current_row, 1).text()
            
            # Perform detailed validation
            validation_results = self.detailed_primer_validation(sequence)
            
            # Create validation dialog
            dialog = QDialog(self)
            dialog.setWindowTitle(f"Primer Validation - {primer_type}")
            dialog.setMinimumSize(500, 400)
            
            layout = QVBoxLayout(dialog)
            
            # Primer info
            info_text = QTextEdit()
            info_text.setReadOnly(True)
            info_text.setPlainText(validation_results)
            layout.addWidget(info_text)
            
            # Close button
            close_btn = QPushButton("Close")
            close_btn.clicked.connect(dialog.accept)
            layout.addWidget(close_btn)
            
            dialog.exec_()
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not validate primer: {str(e)}")
    
    def detailed_primer_validation(self, sequence):
        """Perform detailed primer validation"""
        results = f"Detailed Primer Validation\n{'='*40}\n\n"
        results += f"Sequence: {sequence}\n"
        results += f"Length: {len(sequence)} bp\n\n"
        
        # Basic properties
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')
        gc_percent = (gc_count / len(sequence)) * 100
        tm = self.calculate_tm_improved(sequence)
        
        results += f"Composition Analysis:\n"
        results += f"  A: {sequence.count('A')} ({sequence.count('A')/len(sequence)*100:.1f}%)\n"
        results += f"  T: {sequence.count('T')} ({sequence.count('T')/len(sequence)*100:.1f}%)\n"
        results += f"  G: {sequence.count('G')} ({sequence.count('G')/len(sequence)*100:.1f}%)\n"
        results += f"  C: {sequence.count('C')} ({sequence.count('C')/len(sequence)*100:.1f}%)\n"
        results += f"  GC Content: {gc_percent:.1f}%\n"
        results += f"  Melting Temperature: {tm:.1f}°C\n\n"
        
        # Secondary structure analysis
        results += f"Secondary Structure Analysis:\n"
        
        # Check for hairpins (simplified)
        hairpin_risk = "Low"
        for i in range(len(sequence) - 7):
            substr = sequence[i:i+8]
            rev_comp = str(Seq(substr).reverse_complement())
            if substr in sequence[i+8:]:
                hairpin_risk = "High"
                break
            elif substr[:6] in rev_comp:
                hairpin_risk = "Medium"
        
        results += f"  Hairpin Risk: {hairpin_risk}\n"
        
        # Check for runs
        max_run = 1
        current_run = 1
        for i in range(1, len(sequence)):
            if sequence[i] == sequence[i-1]:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
        
        results += f"  Longest Nucleotide Run: {max_run} {'(Good)' if max_run <= 3 else '(Warning)' if max_run <= 4 else '(Poor)'}\n"
        
        # 3' end analysis
        three_prime_end = sequence[-5:]
        gc_in_3prime = three_prime_end.count('G') + three_prime_end.count('C')
        results += f"  3' End (last 5 bp): {three_prime_end}\n"
        results += f"  GC in 3' end: {gc_in_3prime}/5 {'(Good)' if 1 <= gc_in_3prime <= 3 else '(Suboptimal)'}\n\n"
        
        # Primer-dimer analysis
        results += f"Self-Complementarity Analysis:\n"
        rev_comp_full = str(Seq(sequence).reverse_complement())
        
        # Check for self-complementarity
        max_complement = 0
        for i in range(len(sequence)):
            for j in range(len(rev_comp_full)):
                matches = 0
                k = 0
                while (i + k < len(sequence) and j + k < len(rev_comp_full) and 
                       sequence[i + k] == rev_comp_full[j + k]):
                    matches += 1
                    k += 1
                max_complement = max(max_complement, matches)
        
        results += f"  Max Self-Complementarity: {max_complement} bp "
        if max_complement <= 3:
            results += "(Low risk)\n"
        elif max_complement <= 5:
            results += "(Medium risk)\n"
        else:
            results += "(High risk)\n"
        
        # Overall assessment
        results += f"\nOverall Assessment:\n"
        score = self.calculate_primer_score(sequence, tm, gc_percent, 20, True)
        
        if score >= 80:
            assessment = "Excellent - This primer should work very well"
        elif score >= 60:
            assessment = "Good - This primer should work well with minor considerations"
        elif score >= 40:
            assessment = "Fair - This primer may work but consider alternatives"
        else:
            assessment = "Poor - Consider redesigning this primer"
        
        results += f"  Score: {score:.0f}/100\n"
        results += f"  Assessment: {assessment}\n"
        
        return results
    
    def has_problematic_sequences(self, sequence):
        """Check for problematic sequences in primer"""
        # Check for runs of same nucleotide (4 or more)
        if re.search(r'([ATCG])\1{3,}', sequence):
            return True
        
        # Check for strong secondary structures
        if 'GGGG' in sequence or 'CCCC' in sequence:
            return True
        
        # Check for low complexity regions
        if len(set(sequence)) < 3:  # Less than 3 different bases
            return True
        
        # Check for extreme AT or GC content in any 10-base window
        if len(sequence) >= 10:
            for i in range(len(sequence) - 9):
                window = sequence[i:i+10]
                gc_window = (window.count('G') + window.count('C')) / 10 * 100
                if gc_window < 20 or gc_window > 80:
                    return True
        
        return False
    
    def display_primers(self, primers):
        """Display primers in the table with enhanced formatting"""
        self.results_table.setRowCount(len(primers))
        
        for row, (primer_type, primer) in enumerate(primers):
            # Type
            type_item = QTableWidgetItem(primer_type)
            if primer_type == "Forward":
                type_item.setBackground(QColor(200, 255, 200))  # Light green
            else:
                type_item.setBackground(QColor(255, 200, 200))  # Light red
            self.results_table.setItem(row, 0, type_item)
            
            # Sequence (formatted)
            seq_item = QTableWidgetItem(primer['sequence'])
            seq_item.setFont(QFont("Consolas", 10))
            self.results_table.setItem(row, 1, seq_item)
            
            # Length
            self.results_table.setItem(row, 2, QTableWidgetItem(str(primer['length'])))
            
            # Tm with color coding
            tm_item = QTableWidgetItem(f"{primer['tm']:.1f}")
            tm = primer['tm']
            if 58 <= tm <= 62:
                tm_item.setBackground(QColor(200, 255, 200))  # Good Tm
            elif 55 <= tm <= 65:
                tm_item.setBackground(QColor(255, 255, 200))  # OK Tm
            else:
                tm_item.setBackground(QColor(255, 200, 200))  # Poor Tm
            self.results_table.setItem(row, 3, tm_item)
            
            # GC% with color coding
            gc_item = QTableWidgetItem(f"{primer['gc']:.1f}")
            gc = primer['gc']
            if 45 <= gc <= 55:
                gc_item.setBackground(QColor(200, 255, 200))  # Good GC
            elif 40 <= gc <= 60:
                gc_item.setBackground(QColor(255, 255, 200))  # OK GC
            else:
                gc_item.setBackground(QColor(255, 200, 200))  # Poor GC
            self.results_table.setItem(row, 4, gc_item)
            
            # Position
            self.results_table.setItem(row, 5, QTableWidgetItem(str(primer['position'])))
            
            # Score
            score_item = QTableWidgetItem(f"{primer['score']:.0f}")
            score = primer['score']
            if score >= 80:
                score_item.setBackground(QColor(200, 255, 200))  # Excellent
            elif score >= 60:
                score_item.setBackground(QColor(255, 255, 200))  # Good
            else:
                score_item.setBackground(QColor(255, 200, 200))  # Poor
            self.results_table.setItem(row, 6, score_item)
        
        # Resize columns to content
        self.results_table.resizeColumnsToContents()
    
    def analyze_primer_pairs(self, forward_primers, reverse_primers, start_pos, end_pos):
        """Analyze primer pairs and find the best combination"""
        if not forward_primers or not reverse_primers:
            return
        
        best_pair = None
        best_score = 0
        
        # Analyze top 3 forward and reverse primers
        for fwd in forward_primers[:3]:
            for rev in reverse_primers[:3]:
                pair_score = self.calculate_pair_score(fwd, rev, start_pos, end_pos)
                if pair_score > best_score:
                    best_score = pair_score
                    best_pair = (fwd, rev)
        
        if best_pair:
            fwd, rev = best_pair
            
            # Calculate amplicon properties
            amplicon_start = fwd['start_pos']
            amplicon_end = rev['start_pos'] + rev['length']
            amplicon_length = amplicon_end - amplicon_start
            
            # Tm difference
            tm_diff = abs(fwd['tm'] - rev['tm'])
            
            analysis = f"""Best Primer Pair Analysis (Score: {best_score:.0f}/100):

Forward:  {fwd['sequence']} (Tm: {fwd['tm']:.1f}°C, GC: {fwd['gc']:.1f}%)
Reverse:  {rev['sequence']} (Tm: {rev['tm']:.1f}°C, GC: {rev['gc']:.1f}%)

Amplicon length: {amplicon_length} bp
Tm difference: {tm_diff:.1f}°C {'✓' if tm_diff <= 3 else '⚠' if tm_diff <= 5 else '✗'}
Pair compatibility: {'Excellent' if best_score >= 80 else 'Good' if best_score >= 60 else 'Fair'}"""
            
            self.pair_analysis.setText(analysis)
    
    def calculate_pair_score(self, fwd_primer, rev_primer, start_pos, end_pos):
        """Calculate compatibility score for a primer pair"""
        score = (fwd_primer['score'] + rev_primer['score']) / 2
        
        # Tm difference penalty
        tm_diff = abs(fwd_primer['tm'] - rev_primer['tm'])
        if tm_diff <= 2:
            score += 10
        elif tm_diff <= 5:
            score += 5 - tm_diff
        else:
            score -= tm_diff * 2
        
        # Amplicon length consideration
        amplicon_start = fwd_primer['start_pos']
        amplicon_end = rev_primer['start_pos'] + rev_primer['length']
        amplicon_length = amplicon_end - amplicon_start
        
        target_length = end_pos - start_pos
        length_diff = abs(amplicon_length - target_length)
        
        if length_diff <= 50:
            score += 5
        elif length_diff <= 100:
            score += 2
        else:
            score -= length_diff / 50
        
        # Check for primer-dimer potential (simplified)
        fwd_seq = fwd_primer['sequence']
        rev_seq = rev_primer['sequence']
        
        # Check 3' end complementarity
        fwd_3prime = fwd_seq[-6:]
        rev_3prime = rev_seq[-6:]
        rev_3prime_comp = str(Seq(rev_3prime).reverse_complement())
        
        # Simple check for complementarity
        matches = sum(1 for i, base in enumerate(fwd_3prime) 
                     if i < len(rev_3prime_comp) and base == rev_3prime_comp[i])
        
        if matches >= 4:
            score -= 20  # High primer-dimer risk
        elif matches >= 3:
            score -= 10  # Moderate risk
        
        return max(0, min(100, score))

class SequenceComparisonWidget(QWidget):
    """Widget for sequence comparison and alignment"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.record1 = None
        self.record2 = None
        self.setup_ui()
    
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Sequence inputs
        input_group = QGroupBox("Sequences to Compare")
        input_layout = QVBoxLayout(input_group)
        
        # Sequence 1
        seq1_layout = QHBoxLayout()
        seq1_layout.addWidget(QLabel("Sequence 1:"))
        self.seq1_button = QPushButton("Load from current")
        self.seq1_button.clicked.connect(lambda: self.load_sequence(1))
        seq1_layout.addWidget(self.seq1_button)
        self.seq1_label = QLabel("No sequence loaded")
        seq1_layout.addWidget(self.seq1_label)
        seq1_layout.addStretch()
        input_layout.addLayout(seq1_layout)
        
        # Sequence 2
        seq2_layout = QHBoxLayout()
        seq2_layout.addWidget(QLabel("Sequence 2:"))
        self.seq2_button = QPushButton("Load from file...")
        self.seq2_button.clicked.connect(lambda: self.load_sequence(2))
        seq2_layout.addWidget(self.seq2_button)
        self.seq2_label = QLabel("No sequence loaded")
        seq2_layout.addWidget(self.seq2_label)
        seq2_layout.addStretch()
        input_layout.addLayout(seq2_layout)
        
        layout.addWidget(input_group)
        
        # Comparison options
        options_group = QGroupBox("Comparison Options")
        options_layout = QFormLayout(options_group)
        
        self.comparison_type = QComboBox()
        self.comparison_type.addItems(["Global Alignment", "Local Alignment", "Dot Plot", "Identity Analysis"])
        options_layout.addRow("Analysis type:", self.comparison_type)
        
        compare_btn = QPushButton("Compare Sequences")
        compare_btn.clicked.connect(self.compare_sequences)
        options_layout.addRow("", compare_btn)
        
        layout.addWidget(options_group)
        
        # Results
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        layout.addWidget(self.results_text)
    
    def load_sequence(self, seq_num):
        """Load sequence for comparison"""
        if seq_num == 1:
            # Load from current sequence viewer
            parent = self.parent()
            while parent and not hasattr(parent, 'record'):
                parent = parent.parent()
            
            if parent and parent.record:
                self.record1 = parent.record
                self.seq1_label.setText(f"{parent.record.id} ({len(parent.record.seq)} bp)")
            else:
                QMessageBox.warning(self, "No Sequence", "No sequence loaded in the main viewer.")
        
        else:
            # Load from file
            from PyQt5.QtWidgets import QFileDialog
            from Bio import SeqIO
            
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Load Sequence File", "",
                "Sequence Files (*.gb *.gbk *.fa *.fasta);;All Files (*)"
            )
            
            if file_path:
                try:
                    # Try different formats
                    for fmt in ['genbank', 'fasta']:
                        try:
                            self.record2 = SeqIO.read(file_path, fmt)
                            self.seq2_label.setText(f"{self.record2.id} ({len(self.record2.seq)} bp)")
                            break
                        except:
                            continue
                    else:
                        QMessageBox.warning(self, "Error", "Could not read sequence file.")
                except Exception as e:
                    QMessageBox.warning(self, "Error", f"Error loading file: {str(e)}")
    
    def compare_sequences(self):
        """Compare the loaded sequences"""
        if not self.record1 or not self.record2:
            QMessageBox.warning(self, "Missing Sequences", "Please load both sequences for comparison.")
            return
        
        comparison_type = self.comparison_type.currentText()
        
        if comparison_type == "Identity Analysis":
            self.analyze_identity()
        elif comparison_type == "Dot Plot":
            self.create_dot_plot()
        else:
            self.results_text.setText(f"{comparison_type} not yet implemented.\nThis feature will be available in a future update.")
    
    def analyze_identity(self):
        """Analyze sequence identity"""
        seq1 = str(self.record1.seq).upper()
        seq2 = str(self.record2.seq).upper()
        
        # Simple identity analysis
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))
        
        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
        identity = (matches / min_len) * 100 if min_len > 0 else 0
        
        # Find longest common substring
        longest_match = self.find_longest_common_substring(seq1, seq2)
        
        results = f"""Sequence Identity Analysis
        
        Sequence 1: {self.record1.id} ({len(seq1)} bp)
        Sequence 2: {self.record2.id} ({len(seq2)} bp)
        
        Length difference: {abs(len(seq1) - len(seq2))} bp
        
        Identity (first {min_len} bp): {identity:.2f}%
        Matches: {matches}/{min_len}
        
        Longest common substring: {len(longest_match)} bp
        
        Note: This is a simple position-by-position comparison.
        For more sophisticated alignment, use specialized tools.
        """
        
        self.results_text.setText(results)
    
    def find_longest_common_substring(self, seq1, seq2):
        """Find longest common substring"""
        longest = ""
        
        for i in range(len(seq1)):
            for j in range(len(seq2)):
                k = 0
                while (i + k < len(seq1) and j + k < len(seq2) and 
                       seq1[i + k] == seq2[j + k]):
                    k += 1
                
                if k > len(longest):
                    longest = seq1[i:i + k]
        
        return longest
    
    def create_dot_plot(self):
        """Create a simple dot plot"""
        # This would create a matplotlib dot plot
        # For now, just show a placeholder
        self.results_text.setText("Dot plot visualization will be implemented with matplotlib integration.")

class AdvancedAnalysisDialog(QDialog):
    """Main dialog for advanced sequence analysis"""
    
    def __init__(self, sequence_viewer, parent=None):
        super().__init__(parent)
        self.sequence_viewer = sequence_viewer
        self.setWindowTitle("Advanced Sequence Analysis")
        self.setMinimumSize(800, 600)
        self.setup_ui()
    
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Create tab widget
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)
        
        # ORF Analysis tab
        self.orf_widget = ORFAnalysisWidget()
        self.tab_widget.addTab(self.orf_widget, "ORF Analysis")
        
        # Primer Design tab
        self.primer_widget = PrimerDesignWidget()
        self.tab_widget.addTab(self.primer_widget, "Primer Design")
        
        # Sequence Comparison tab
        self.comparison_widget = SequenceComparisonWidget()
        self.tab_widget.addTab(self.comparison_widget, "Sequence Comparison")
        
        # Set sequence for all widgets
        if hasattr(self.sequence_viewer, 'record') and self.sequence_viewer.record:
            self.orf_widget.set_sequence(self.sequence_viewer.record)
            self.primer_widget.set_sequence(self.sequence_viewer.record)
        
        # Close button
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)