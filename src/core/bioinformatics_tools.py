#!/usr/bin/env python3
"""
Bioinformatics Tools for PicoMol.

This module provides various bioinformatics analysis tools including:
- Sequence analysis and statistics
- Protein property calculations
- Secondary structure prediction
- Sequence alignment tools
- Motif finding
- And more...
"""

import os
import re
from io import StringIO
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTabWidget, QGroupBox, QFormLayout,
    QTextEdit, QLineEdit, QPushButton, QLabel, QComboBox, QSpinBox,
    QCheckBox, QTableWidget, QTableWidgetItem, QHeaderView, QScrollArea,
    QMessageBox, QFileDialog, QProgressBar, QSplitter, QDialog, QDialogButtonBox,
    QListWidget, QListWidgetItem, QSplitter as QSplitterWidget
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont

try:
    from Bio.Seq import Seq
    from Bio.SeqUtils import molecular_weight
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    
    # Handle different Biopython versions for GC content
    try:
        from Bio.SeqUtils import GC
        gc_content_func = GC
    except ImportError:
        from Bio.SeqUtils import gc_fraction
        gc_content_func = lambda seq: gc_fraction(seq) * 100  # Convert to percentage
    
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    gc_content_func = None


class FastaSequence:
    """Class to represent a FASTA sequence with metadata."""
    def __init__(self, header, sequence, sequence_type=None):
        self.header = header
        self.sequence = sequence.upper().replace(' ', '').replace('\n', '')
        self.sequence_type = sequence_type or self._guess_sequence_type()
        self.length = len(self.sequence)
    
    def _guess_sequence_type(self):
        """Guess the sequence type based on composition."""
        sequence_set = set(self.sequence)
        
        # Check for DNA
        if sequence_set <= set('ATGCNRYSWKMBDHV'):
            return 'dna'
        # Check for RNA
        elif sequence_set <= set('AUGCNRYSWKMBDHV'):
            return 'rna'
        # Assume protein
        else:
            return 'protein'
    
    def get_display_name(self):
        """Get a display name for the sequence."""
        # Extract meaningful part of header
        header_parts = self.header.split('|')
        if len(header_parts) > 1:
            # Handle common formats like >gi|123456|ref|NP_123456.1|
            return header_parts[-2] if header_parts[-1] == '' else header_parts[-1]
        else:
            # Simple header
            return self.header[:50] + '...' if len(self.header) > 50 else self.header
    
    def __str__(self):
        return f"{self.get_display_name()} ({self.sequence_type}, {self.length} residues)"


class SequenceSelectionDialog(QDialog):
    """Dialog for selecting sequences from a FASTA file."""
    
    def __init__(self, sequences, parent=None):
        super().__init__(parent)
        self.sequences = sequences
        self.selected_sequence = None
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Select Sequence")
        self.setMinimumSize(600, 400)
        
        layout = QVBoxLayout(self)
        
        # Info label
        info_label = QLabel(f"Found {len(self.sequences)} sequences. Select one to analyze:")
        layout.addWidget(info_label)
        
        # Create splitter for list and preview
        splitter = QSplitterWidget(Qt.Horizontal)
        
        # Sequence list
        self.sequence_list = QListWidget()
        self.sequence_list.setMinimumWidth(300)
        
        for i, seq in enumerate(self.sequences):
            item = QListWidgetItem(str(seq))
            item.setData(Qt.UserRole, i)  # Store index
            self.sequence_list.addItem(item)
        
        self.sequence_list.currentItemChanged.connect(self.on_selection_changed)
        splitter.addWidget(self.sequence_list)
        
        # Preview area
        preview_widget = QWidget()
        preview_layout = QVBoxLayout(preview_widget)
        
        preview_layout.addWidget(QLabel("Preview:"))
        
        self.preview_text = QTextEdit()
        self.preview_text.setReadOnly(True)
        self.preview_text.setFont(QFont("Courier", 9))
        self.preview_text.setMaximumHeight(150)
        preview_layout.addWidget(self.preview_text)
        
        # Sequence info
        self.info_label = QLabel("Select a sequence to see details")
        self.info_label.setWordWrap(True)
        preview_layout.addWidget(self.info_label)
        
        splitter.addWidget(preview_widget)
        splitter.setSizes([300, 300])
        
        layout.addWidget(splitter)
        
        # Buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        # Select first item by default
        if self.sequences:
            self.sequence_list.setCurrentRow(0)
    
    def on_selection_changed(self, current, previous):
        """Handle selection change."""
        if current is None:
            return
        
        index = current.data(Qt.UserRole)
        sequence = self.sequences[index]
        
        # Update preview
        preview_seq = sequence.sequence[:200] + ('...' if len(sequence.sequence) > 200 else '')
        self.preview_text.setPlainText(f">{sequence.header}\n{preview_seq}")
        
        # Update info
        info_text = f"<b>Header:</b> {sequence.header}<br>"
        info_text += f"<b>Type:</b> {sequence.sequence_type.upper()}<br>"
        info_text += f"<b>Length:</b> {sequence.length} residues<br>"
        
        # Add composition info
        if sequence.sequence_type == 'protein':
            aa_count = len([c for c in sequence.sequence if c in 'ACDEFGHIKLMNPQRSTVWY'])
            info_text += f"<b>Valid amino acids:</b> {aa_count}/{sequence.length}"
        elif sequence.sequence_type in ['dna', 'rna']:
            if sequence.sequence_type == 'dna':
                valid_chars = 'ATGCN'
            else:
                valid_chars = 'AUGCN'
            valid_count = len([c for c in sequence.sequence if c in valid_chars])
            info_text += f"<b>Valid nucleotides:</b> {valid_count}/{sequence.length}"
            
            # Calculate GC content for preview
            if sequence.sequence_type == 'dna':
                gc_count = sequence.sequence.count('G') + sequence.sequence.count('C')
            else:
                gc_count = sequence.sequence.count('G') + sequence.sequence.count('C')
            gc_percent = (gc_count / sequence.length * 100) if sequence.length > 0 else 0
            info_text += f"<br><b>GC Content:</b> {gc_percent:.1f}%"
        
        self.info_label.setText(info_text)
        self.selected_sequence = sequence
    
    def get_selected_sequence(self):
        """Get the selected sequence."""
        return self.selected_sequence


def parse_fasta_file(file_path):
    """Parse a FASTA file and return a list of FastaSequence objects."""
    sequences = []
    
    try:
        if BIOPYTHON_AVAILABLE:
            # Use Biopython for robust parsing
            for record in SeqIO.parse(file_path, "fasta"):
                header = record.description
                sequence = str(record.seq)
                if sequence:  # Only add non-empty sequences
                    fasta_seq = FastaSequence(header, sequence)
                    sequences.append(fasta_seq)
        else:
            # Fallback manual parsing
            with open(file_path, 'r', encoding='utf-8') as f:
                current_header = None
                current_sequence = []
                
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_header and current_sequence:
                            sequence = ''.join(current_sequence)
                            if sequence:  # Only add non-empty sequences
                                fasta_seq = FastaSequence(current_header, sequence)
                                sequences.append(fasta_seq)
                        
                        # Start new sequence
                        current_header = line[1:]  # Remove >
                        current_sequence = []
                    else:
                        # Add to current sequence
                        current_sequence.append(line)
                
                # Don't forget the last sequence
                if current_header and current_sequence:
                    sequence = ''.join(current_sequence)
                    if sequence:  # Only add non-empty sequences
                        fasta_seq = FastaSequence(current_header, sequence)
                        sequences.append(fasta_seq)
    
    except Exception as e:
        raise Exception(f"Error parsing FASTA file: {str(e)}")
    
    return sequences


class SequenceAnalysisWorker(QThread):
    """Worker thread for sequence analysis to prevent UI freezing."""
    
    analysis_complete = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)
    
    def __init__(self, sequence, sequence_type):
        super().__init__()
        self.sequence = sequence
        self.sequence_type = sequence_type
    
    def run(self):
        try:
            results = {}
            
            if BIOPYTHON_AVAILABLE:
                # Use Biopython for full analysis
                if self.sequence_type == "protein":
                    results = self.analyze_protein_sequence(self.sequence)
                elif self.sequence_type == "dna":
                    results = self.analyze_dna_sequence(self.sequence)
                elif self.sequence_type == "rna":
                    results = self.analyze_rna_sequence(self.sequence)
            else:
                # Use basic analysis without Biopython
                if self.sequence_type == "protein":
                    results = self.analyze_protein_sequence_basic(self.sequence)
                elif self.sequence_type == "dna":
                    results = self.analyze_dna_sequence_basic(self.sequence)
                elif self.sequence_type == "rna":
                    results = self.analyze_rna_sequence_basic(self.sequence)
            
            self.analysis_complete.emit(results)
            
        except Exception as e:
            self.error_occurred.emit(f"Error during analysis: {str(e)}")
    
    def analyze_protein_sequence(self, sequence):
        """Analyze protein sequence properties."""
        seq_obj = Seq(sequence)
        protein_analysis = ProteinAnalysis(str(seq_obj))
        
        results = {
            'length': len(sequence),
            'molecular_weight': protein_analysis.molecular_weight(),
            'aromaticity': protein_analysis.aromaticity(),
            'instability_index': protein_analysis.instability_index(),
            'isoelectric_point': protein_analysis.isoelectric_point(),
            'gravy': protein_analysis.gravy(),  # Grand average of hydropathy
            'amino_acid_percent': protein_analysis.amino_acids_percent,
            'secondary_structure': protein_analysis.secondary_structure_fraction(),
            'flexibility': protein_analysis.flexibility(),
            'charge_at_pH': {
                'pH_7': protein_analysis.charge_at_pH(7.0),
                'pH_5': protein_analysis.charge_at_pH(5.0),
                'pH_9': protein_analysis.charge_at_pH(9.0)
            }
        }
        
        # Add note about partial codons if sequence length is not divisible by 3
        if len(sequence) % 3 != 0:
            results['translation_note'] = f"Note: Sequence length ({len(sequence)} bp) is not divisible by 3. Translation may include partial codons."
        
        return results
    
    def analyze_protein_sequence_basic(self, sequence):
        """Basic protein sequence analysis without Biopython."""
        # Count amino acids
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        aa_counts = {aa: sequence.count(aa) for aa in amino_acids}
        aa_percent = {aa: (count / len(sequence)) * 100 for aa, count in aa_counts.items()}
        
        # Basic molecular weight calculation (approximate)
        aa_weights = {
            'A': 89.1, 'C': 121.0, 'D': 133.1, 'E': 147.1, 'F': 165.2,
            'G': 75.1, 'H': 155.2, 'I': 131.2, 'K': 146.2, 'L': 131.2,
            'M': 149.2, 'N': 132.1, 'P': 115.1, 'Q': 146.2, 'R': 174.2,
            'S': 105.1, 'T': 119.1, 'V': 117.1, 'W': 204.2, 'Y': 181.2
        }
        
        mol_weight = sum(aa_weights.get(aa, 0) * count for aa, count in aa_counts.items())
        mol_weight -= (len(sequence) - 1) * 18.015  # Subtract water for peptide bonds
        
        results = {
            'length': len(sequence),
            'molecular_weight': mol_weight,
            'amino_acid_percent': aa_percent,
            'note': 'Basic analysis - install Biopython for advanced features'
        }
        
        return results
    
    def analyze_dna_sequence_basic(self, sequence):
        """Basic DNA sequence analysis without Biopython."""
        # Count nucleotides
        nucleotide_counts = {
            'A': sequence.count('A'),
            'T': sequence.count('T'),
            'G': sequence.count('G'),
            'C': sequence.count('C'),
            'N': sequence.count('N')
        }
        
        # Calculate GC content
        gc_count = nucleotide_counts['G'] + nucleotide_counts['C']
        total_count = sum(nucleotide_counts.values())
        gc_content = (gc_count / total_count) * 100 if total_count > 0 else 0
        
        # Basic complement
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        complement = ''.join(complement_map.get(base, base) for base in sequence)
        reverse_complement = complement[::-1]
        
        # Basic translation (first frame only)
        codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        translation = ''
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            translation += codon_table.get(codon, 'X')
        
        results = {
            'length': len(sequence),
            'gc_content': gc_content,
            'nucleotide_counts': nucleotide_counts,
            'nucleotide_percent': {
                base: (count / len(sequence)) * 100 
                for base, count in nucleotide_counts.items()
            },
            'complement': complement,
            'reverse_complement': reverse_complement,
            'translation_frames': {
                'frame_1': translation
            },
            'note': 'Basic analysis - install Biopython for advanced features'
        }
        
        return results
    
    def analyze_rna_sequence_basic(self, sequence):
        """Basic RNA sequence analysis without Biopython."""
        # Count nucleotides
        nucleotide_counts = {
            'A': sequence.count('A'),
            'U': sequence.count('U'),
            'G': sequence.count('G'),
            'C': sequence.count('C'),
            'N': sequence.count('N')
        }
        
        # Calculate GC content
        gc_count = nucleotide_counts['G'] + nucleotide_counts['C']
        total_count = sum(nucleotide_counts.values())
        gc_content = (gc_count / total_count) * 100 if total_count > 0 else 0
        
        # Basic complement
        complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        complement = ''.join(complement_map.get(base, base) for base in sequence)
        reverse_complement = complement[::-1]
        
        # Basic translation
        codon_table = {
            'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
            'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
            'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        translation = ''
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            translation += codon_table.get(codon, 'X')
        
        results = {
            'length': len(sequence),
            'gc_content': gc_content,
            'nucleotide_counts': nucleotide_counts,
            'nucleotide_percent': {
                base: (count / len(sequence)) * 100 
                for base, count in nucleotide_counts.items()
            },
            'complement': complement,
            'reverse_complement': reverse_complement,
            'translation': translation,
            'note': 'Basic analysis - install Biopython for advanced features'
        }
        
        return results
    
    def analyze_dna_sequence(self, sequence):
        """Analyze DNA sequence properties."""
        seq_obj = Seq(sequence)
        
        # Count nucleotides
        nucleotide_counts = {
            'A': sequence.count('A'),
            'T': sequence.count('T'),
            'G': sequence.count('G'),
            'C': sequence.count('C'),
            'N': sequence.count('N')  # Unknown nucleotides
        }
        
        # Helper function to safely translate sequences
        def safe_translate(seq):
            try:
                # Suppress the warning by using to_stop=False and handling partial codons
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    return str(seq.translate(to_stop=False))
            except Exception:
                return "Translation failed"
        
        results = {
            'length': len(sequence),
            'gc_content': gc_content_func(seq_obj),
            'molecular_weight': molecular_weight(seq_obj, seq_type='DNA'),
            'nucleotide_counts': nucleotide_counts,
            'nucleotide_percent': {
                base: (count / len(sequence)) * 100 
                for base, count in nucleotide_counts.items()
            },
            'complement': str(seq_obj.complement()),
            'reverse_complement': str(seq_obj.reverse_complement()),
            'translation_frames': {
                'frame_1': safe_translate(seq_obj),
                'frame_2': safe_translate(seq_obj[1:]),
                'frame_3': safe_translate(seq_obj[2:]),
                'frame_-1': safe_translate(seq_obj.reverse_complement()),
                'frame_-2': safe_translate(seq_obj.reverse_complement()[1:]),
                'frame_-3': safe_translate(seq_obj.reverse_complement()[2:])
            }
        }
        
        # Add note about partial codons if sequence length is not divisible by 3
        if len(sequence) % 3 != 0:
            results['translation_note'] = f"Note: Sequence length ({len(sequence)} bp) is not divisible by 3. Translation may include partial codons."
        
        return results
    
    def analyze_rna_sequence(self, sequence):
        """Analyze RNA sequence properties."""
        seq_obj = Seq(sequence)
        
        # Count nucleotides
        nucleotide_counts = {
            'A': sequence.count('A'),
            'U': sequence.count('U'),
            'G': sequence.count('G'),
            'C': sequence.count('C'),
            'N': sequence.count('N')  # Unknown nucleotides
        }
        
        # Helper function to safely translate sequences
        def safe_translate(seq):
            try:
                # Suppress the warning by using to_stop=False and handling partial codons
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    return str(seq.translate(to_stop=False))
            except Exception:
                return "Translation failed"
        
        results = {
            'length': len(sequence),
            'gc_content': gc_content_func(seq_obj),
            'molecular_weight': molecular_weight(seq_obj, seq_type='RNA'),
            'nucleotide_counts': nucleotide_counts,
            'nucleotide_percent': {
                base: (count / len(sequence)) * 100 
                for base, count in nucleotide_counts.items()
            },
            'complement': str(seq_obj.complement()),
            'reverse_complement': str(seq_obj.reverse_complement()),
            'translation': safe_translate(seq_obj)
        }
        
        # Add note about partial codons if sequence length is not divisible by 3
        if len(sequence) % 3 != 0:
            results['translation_note'] = f"Note: Sequence length ({len(sequence)} bp) is not divisible by 3. Translation may include partial codons."
        
        return results


class SequenceAnalysisTab(QWidget):
    """Tab for sequence analysis tools."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_app = parent
        self.analysis_worker = None
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)
        
        # Input section
        input_group = QGroupBox("Sequence Input")
        input_layout = QVBoxLayout(input_group)
        
        # Sequence type selection
        type_layout = QHBoxLayout()
        type_layout.addWidget(QLabel("Sequence Type:"))
        self.sequence_type_combo = QComboBox()
        self.sequence_type_combo.addItems(["protein", "dna", "rna"])
        self.sequence_type_combo.setToolTip("Select the type of sequence you want to analyze")
        type_layout.addWidget(self.sequence_type_combo)
        
        # Load from structure button
        load_structure_btn = QPushButton("Load from Current Structure")
        load_structure_btn.setToolTip("Load sequence from the currently displayed protein structure")
        load_structure_btn.clicked.connect(self.load_from_structure)
        type_layout.addWidget(load_structure_btn)
        
        # Load from file button
        load_file_btn = QPushButton("Load from File...")
        load_file_btn.setToolTip("Load sequence from a FASTA file")
        load_file_btn.clicked.connect(self.load_from_file)
        type_layout.addWidget(load_file_btn)
        
        type_layout.addStretch()
        input_layout.addLayout(type_layout)
        
        # Sequence input area
        self.sequence_input = QTextEdit()
        self.sequence_input.setPlaceholderText(
            "Enter your sequence here or use the buttons above to load from structure/file...\n\n"
            "Examples:\n"
            "Protein: MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL\n"
            "DNA: ATGAAATGGGTAACATTTTCGATCCTTCTGCTCTTCCTGTTCTCAAGCGCCTACAGCCGTGGAGTATTCCGCCGTGATGCCCACAAATCAGAAGTGGCCCACCGCTTCAAAGACCTGGGAGAAGAAAACTTCAAAGCCCTGGTGCTGATCGCCTTCGCCCAGTACCTGCAGCAGTGCCCCTTTGAAGACCACGTGAAACTGGTGAACGAAGTGACCGAGTTCGCCAAGACCTGCGTGGCCGACGAATCCGCCGAGAACTGCGACAAATCCCTGCACACCCTGTTCGGAGACAAACTGTGCACCGTGGCCACCCTGCGTGAGACCTACGGAGAAATGGCCGACTGCTGCGCCAAGCAGGAACCAGAGCGTAACGAGTGCTTCCTGCAGCACAAAGACGACAACCCAAACCTGCCCCGTCTGGTGCGTCCAGAAGTGGACGTGATGTGCACCGCCTTCCACGACAACGAAGAGACCTTCCTGAAGAAATACCTGTACGAGATCGCCCGTCGTCACCCATACTTCTACGCCCCAGAACTGCTGTTCTTCGCCAAGCGCTACAAAGCCGCCTTCACCGAGTGCTGCCAGGCCGCCGACAAAGCCGCCTGCCTGCTGCCAAAGCTGGACGAACTGCGTGACGAAGGAAAAGCCTCCTCCGCCAAGCAGCGTCTGAAGTGCGCCTCCCTGCAGAAGTTCGGAGAGCGTGCCTTCAAAGCCTGGGCCGTGGCCCGTCTGTCCCAGCGTTTCCCAAAAGCCGAGTTCGCCGAAGTGTCCAAGCTGGTGACCGACCTGACCAAGGTGCACACCGAGTGCTGCCACGGAGACCTGCTGGAGTGCGCCGACGACCGTGCCGACCTGGCCAAGTACATCTGCGAGAACCAGGACTCCATCTCCTCCAAGCTGAAAGAGTGCTGCGAGAAGCCACTGCTGGAGAAATCCCACTGCATCGCCGAAGTGGAGAACGACGAGATGCCAGCCGACCTGCCCTCCCTGGCCGCCGACTTCGTGGAATCCAAAGACGTGTGCAAGAACTACGCCGAAGCCAAAGACGTGTTCCTGGGAATGTTCCTGTACGAGTACGCCCGTCGTCACCCAGACTACTCCGTGGTGCTGCTGCTGCGTCTGGCCAAGACCTACGAGACCACCCTGGAGAAGTGCTGCGCCGCCGCCGACCCACACGAGTGCTACGCCAAGGTGTTCGACGAGTTCAAGCCACTGGTGGAAGAACCACAGAACCTGATCAAGCAGAACTGCGAACTGTTCGAGCAGCTGGGAGAGTACAAGTTCCAGAACGCCCTGCTGGTGCGCTACACCAAGAAGGTGCCACAGGTGTCCACCCCAACCCTGGTGGAAGTGTCCCGTAACCTGGGAAAGGTGGGCTCCAAGTGCTGCAAGCACCCAGAAGCCAAGCGTATGCCATGCGCCGAGGACTACCTGTCCGTGGTGCTGAACCAGCTGTGCGTGCTGCACGAGAAGACCCCAGTGTCCGACCGTGTGACCAAGTGCTGCACCGAATCCCTGGTGAACCGTCGTCCATGCTTCTCCGCCCTGGAAGTGGACGAGACCTACGTGCCAAAAGAGTTCAACGCCGAGACCTTCACCTTCCACGCCGACATCTGCACCCTGTCCGAGAAAGAGCGTCAGATCAAGAAGCAGACCGCCCTGGTGGAACTGGTGAAGCACAAGCCAAAAGCCACCAAAGAGCAGCTGAAAGCCGTGATGGACGACTTCGCCGCCTTCGTGGAGAAGTGCTGCAAAGCCGACGACAAAGAGACCTGCTTCGCCGAAGAAGGAAAGAAGCTGGTGGCCGCCTCCCAGGCCGCCCTGGGACTGTAG\n"
            "RNA: AUGAAAUGGGUAACAUUUUCGAUCCUUCUGCUCUUCCUGUUCUCAAGCGCCUACAGCCGUGGAGUAUUCCGCCGUGAUGCCCACAAAUCAGAAGUGGGCCCACCGCUUCAAAGACCUGGGAGAAGAAAACUUCAAAGCCCUGGUGGUGAUCGCCUUCGCCCAGUACCUGCAGCAGUGCCCUUUGAAGACCACGUGAAACUGGUUGAACGAAGUGACCGAGUUCGCCAAGACCUGCGUGGCCGACGAAUCCGCCGAGAACUGCGACAAAUCCCUGCACACCCUGUUCGGAGACAAACUGUGCACCGUGGCCACCCUGCGUGAGACCUACGGAGAAAUGGCCGACUGCUGCGCCAAGCAGGAACCAGAGCGUAACGAGUGCUUCCUGCAGCACAAAGACGACAACCCAAACCUGCCCCGUCUGGUUGCGUCCAGAAGUGGGACGUGAUGUGCACCGCCUUCCACGACAACGAAGAGACCUUCCUGAAGAAAUACCUGUACGAGAUCGCCCGUCGUCACCCAUACUUCUACGCCCCAGAACUGCUGUUCUUCGCCAAGCGCUACAAAGCCGCCUUCACCGAGUGCUGCCAGGCCGCCGACAAAGCCGCCUGCCUGCUGCCAAAGCUGGACGAACUGCGUGACGAAGGAAAAGCCUCCUCCGCCAAGCAGCGUCUGAAGUGCGCCUCCCCUGCAGAAGUUCGGAGAGCGUGCCUUCAAAGCCUGGGCCGUGGCCCGUCUGUCCCAGCGUUUCCCAAAAGCCGAGUUCGCCGAAGUGUCCAAGCUGGUUGACCGACCUGACCAAGUGCACACCGAGUGCUGCCACGGAGACCUGCUGGAGUGCGCCGACGACCGUGCCGACCUGGCCAAGUACAUCUGCGAGAACCAGGACUCCAUCUCCCUCCAAGCUGAAAGAGUGCUGCGAGAAGCCACUGCUGGAGAAAUCCCACUGCAUCGCCGAAGUGGGAGAACGACGAGAUGCCAGCCGACCUGCCCUCCCCUGGCCGCCGACUUCGUGGAAUCCCAAAGACGUGUGCAAGAACUACGCCGAAGCCAAAGACGUGUUCCUGGGAAUGUUCCUGUACGAGUACGCCCGUCGUCACCCAGACUACUCCGUGGUUGCUGCUGCUGCGUCUGGCCAAGACCUACGAGACCACCCUGGAGAAGUGCUGCGCCGCCGCCGACCCACACGAGUGCUACGCCAAGUGGUUCGACGAGUUCAAGCCACUGGUGGAAGAACCACAGAACCUGAUCAAGCAGAACUGCGAACUGUUCGAGCAGCUGGGAGAGUACAAGUUCCAGAACGCCCUGCUGGUUGCGUACACCAAGAAGUGUCCACAGGUGUCCACCCCAACCCUGGUGGAAGUGUCCCGUAACCUGGGAAAGGUGGGGCUCCAAGUGCUGCAAGCACCCAGAAGCCAAGCGUAUGCCAUGCGCCGAGGACUACCUGUCCGUGGUUGCUGAACCAGCUGUGCGUGCUGCACGAGAAGACCCCAGUGUCCGACCGUGUGACCAAGUGCUGCACCGAAUCCCCUGGUUGAACCGUCGUCCAUGCUUCUCCGCCCUGGAAGUGGGACGAGACCUACGUGCCAAAAGAGUUCAACGCCGAGACCUUCACCUUCCACGCCGACAUCUGCACCCUGUCCGAGAAAGAGCGUCAGAUCAAGAAGCAGACCGCCCUGGUGGAACUGGUUGAAGCACAAGCCAAAAGCCACCAAAGAGCAGCUGAAAGCCGUGAUGGACGACUUCGCCGCCUUCGUGGAGAAGUGCUGCAAAGCCGACGACAAAGAGACCUGCUUCGCCGAAGAAGGAAAGAAGCUGGUGGCCGCCUCCCAGGCCGCCCUGGGACUGUAG"
        )
        self.sequence_input.setMaximumHeight(150)
        self.sequence_input.setFont(QFont("Courier", 10))
        input_layout.addWidget(self.sequence_input)
        
        # Analysis button
        analyze_btn = QPushButton("Analyze Sequence")
        analyze_btn.setToolTip("Perform comprehensive analysis of the input sequence")
        analyze_btn.clicked.connect(self.analyze_sequence)
        input_layout.addWidget(analyze_btn)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        input_layout.addWidget(self.progress_bar)
        
        layout.addWidget(input_group)
        
        # Results section
        self.results_area = QScrollArea()
        self.results_widget = QWidget()
        self.results_layout = QVBoxLayout(self.results_widget)
        self.results_area.setWidget(self.results_widget)
        self.results_area.setWidgetResizable(True)
        
        layout.addWidget(self.results_area, 1)
        
        # Initially show a placeholder
        self.show_placeholder()
    
    def show_placeholder(self):
        """Show placeholder text when no analysis has been performed."""
        placeholder = QLabel(
            "<h3>Sequence Analysis Tools</h3>"
            "<p>Enter a sequence above and click 'Analyze Sequence' to get:</p>"
            "<ul>"
            "<li><b>Basic Properties:</b> Length, molecular weight, composition</li>"
            "<li><b>Protein Analysis:</b> Isoelectric point, hydropathy, secondary structure</li>"
            "<li><b>DNA/RNA Analysis:</b> GC content, translation frames, complements</li>"
            "<li><b>Detailed Statistics:</b> Amino acid/nucleotide frequencies</li>"
            "</ul>"
            "<p><i>Tip: You can load sequences from the current structure or from FASTA files.</i></p>"
        )
        placeholder.setAlignment(Qt.AlignTop)
        placeholder.setWordWrap(True)
        placeholder.setStyleSheet("color: #666; padding: 20px;")
        
        # Clear existing layout
        for i in reversed(range(self.results_layout.count())):
            item = self.results_layout.itemAt(i)
            if item is not None:
                widget = item.widget()
                if widget is not None:
                    widget.setParent(None)
                else:
                    # Handle spacer items
                    self.results_layout.removeItem(item)
        
        self.results_layout.addWidget(placeholder)
        self.results_layout.addStretch()
    
    def load_from_structure(self):
        """Load sequence from the currently displayed structure."""
        if hasattr(self.parent_app, 'sequence_display'):
            sequence_text = self.parent_app.sequence_display.toPlainText()
            if sequence_text and not sequence_text.startswith("Sequence will appear"):
                # Extract just the sequence part (remove FASTA headers)
                lines = sequence_text.split('\n')
                sequence_lines = [line for line in lines if not line.startswith('>')]
                sequence = ''.join(sequence_lines).replace(' ', '').upper()
                
                if sequence:
                    self.sequence_input.setPlainText(sequence)
                    self.sequence_type_combo.setCurrentText("protein")  # Assume protein from PDB
                    self.parent_app.statusBar().showMessage("Sequence loaded from current structure")
                else:
                    QMessageBox.warning(self, "No Sequence", "No valid sequence found in the current structure.")
            else:
                QMessageBox.warning(self, "No Structure", "No structure is currently loaded.")
        else:
            QMessageBox.warning(self, "Error", "Cannot access sequence display.")
    
    def load_from_file(self):
        """Load sequence from a FASTA file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Sequence File", "", 
            "FASTA Files (*.fasta *.fa *.fas);;GenBank Files (*.gb *.gbk);;Text Files (*.txt);;All Files (*)"
        )
        
        if file_path:
            try:
                # Parse the file
                sequences = parse_fasta_file(file_path)
                
                if not sequences:
                    QMessageBox.warning(self, "No Sequences", "No valid sequences found in the file.")
                    return
                
                # If only one sequence, load it directly
                if len(sequences) == 1:
                    selected_sequence = sequences[0]
                else:
                    # Multiple sequences - show selection dialog
                    dialog = SequenceSelectionDialog(sequences, self)
                    if dialog.exec_() != QDialog.Accepted:
                        return
                    
                    selected_sequence = dialog.get_selected_sequence()
                    if not selected_sequence:
                        return
                
                # Load the selected sequence
                self.sequence_input.setPlainText(selected_sequence.sequence)
                self.sequence_type_combo.setCurrentText(selected_sequence.sequence_type)
                
                # Update status
                filename = os.path.basename(file_path)
                if len(sequences) == 1:
                    status_msg = f"Loaded sequence from {filename}"
                else:
                    status_msg = f"Loaded sequence '{selected_sequence.get_display_name()}' from {filename}"
                
                if hasattr(self.parent_app, 'statusBar'):
                    self.parent_app.statusBar().showMessage(status_msg)
                
                # Show info about the loaded sequence
                info_msg = f"Loaded {selected_sequence.sequence_type.upper()} sequence\n"
                info_msg += f"Length: {selected_sequence.length} residues\n"
                info_msg += f"Header: {selected_sequence.header[:100]}{'...' if len(selected_sequence.header) > 100 else ''}"
                
                QMessageBox.information(self, "Sequence Loaded", info_msg)
                        
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load sequence file:\n{str(e)}")
    
    def parse_input_text(self, input_text):
        """Parse input text and return (sequence, sequence_type).
        
        Handles both plain sequences and FASTA format.
        Returns (sequence_string, detected_type) or (sequence_string, None) if no type detected.
        """
        lines = input_text.strip().split('\n')
        
        # Check if it looks like FASTA format (starts with >)
        if lines[0].startswith('>'):
            # Parse as FASTA
            try:
                # Create a temporary file-like object
                from io import StringIO
                fasta_buffer = StringIO(input_text)
                
                if BIOPYTHON_AVAILABLE:
                    # Use Biopython to parse
                    sequences = []
                    for record in SeqIO.parse(fasta_buffer, "fasta"):
                        if record.seq:
                            fasta_seq = FastaSequence(record.description, str(record.seq))
                            sequences.append(fasta_seq)
                else:
                    # Manual parsing
                    sequences = []
                    current_header = None
                    current_sequence = []
                    
                    for line in lines:
                        line = line.strip()
                        if not line:
                            continue
                        
                        if line.startswith('>'):
                            # Save previous sequence if exists
                            if current_header and current_sequence:
                                sequence = ''.join(current_sequence)
                                if sequence:
                                    fasta_seq = FastaSequence(current_header, sequence)
                                    sequences.append(fasta_seq)
                            
                            # Start new sequence
                            current_header = line[1:]  # Remove >
                            current_sequence = []
                        else:
                            # Add to current sequence
                            current_sequence.append(line)
                    
                    # Don't forget the last sequence
                    if current_header and current_sequence:
                        sequence = ''.join(current_sequence)
                        if sequence:
                            fasta_seq = FastaSequence(current_header, sequence)
                            sequences.append(fasta_seq)
                
                if sequences:
                    if len(sequences) == 1:
                        # Single sequence - use it
                        selected_seq = sequences[0]
                        return selected_seq.sequence, selected_seq.sequence_type
                    else:
                        # Multiple sequences - show selection dialog
                        dialog = SequenceSelectionDialog(sequences, self)
                        if dialog.exec_() == QDialog.Accepted:
                            selected_seq = dialog.get_selected_sequence()
                            if selected_seq:
                                return selected_seq.sequence, selected_seq.sequence_type
                        return None, None
                else:
                    return None, None
                    
            except Exception as e:
                # If FASTA parsing fails, treat as plain sequence
                print(f"FASTA parsing failed: {e}")
                pass
        
        # Treat as plain sequence
        # Remove any remaining FASTA headers and clean up
        sequence_lines = []
        for line in lines:
            line = line.strip()
            if line and not line.startswith('>'):
                sequence_lines.append(line)
        
        sequence = ''.join(sequence_lines).replace(' ', '').upper()
        
        # Try to auto-detect sequence type
        if sequence:
            sequence_set = set(sequence)
            if sequence_set <= set('ATGCNRYSWKMBDHV'):
                detected_type = 'dna'
            elif sequence_set <= set('AUGCNRYSWKMBDHV'):
                detected_type = 'rna'
            else:
                detected_type = 'protein'
            
            return sequence, detected_type
        
        return None, None
    
    def analyze_sequence(self):
        """Analyze the input sequence."""
        input_text = self.sequence_input.toPlainText().strip()
        
        if not input_text:
            QMessageBox.warning(self, "No Sequence", "Please enter a sequence to analyze.")
            return
        
        # Check if input looks like FASTA format
        sequence, sequence_type = self.parse_input_text(input_text)
        
        if not sequence:
            QMessageBox.warning(self, "No Sequence", "No valid sequence found in the input.")
            return
        
        # Update sequence type if auto-detected from FASTA
        if sequence_type:
            self.sequence_type_combo.setCurrentText(sequence_type)
        
        # Get final sequence type
        final_sequence_type = self.sequence_type_combo.currentText()
        
        # Validate sequence
        if not self.validate_sequence(sequence, final_sequence_type):
            return
        
        # Show progress
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        
        # Start analysis in worker thread
        self.analysis_worker = SequenceAnalysisWorker(sequence, final_sequence_type)
        self.analysis_worker.analysis_complete.connect(self.display_results)
        self.analysis_worker.error_occurred.connect(self.handle_analysis_error)
        self.analysis_worker.start()
    
    def validate_sequence(self, sequence, sequence_type):
        """Validate sequence based on type."""
        if sequence_type == "protein":
            valid_chars = set('ACDEFGHIKLMNPQRSTVWYXBZJUO*')
            invalid_chars = set(sequence) - valid_chars
            if invalid_chars:
                QMessageBox.warning(
                    self, "Invalid Sequence", 
                    f"Invalid characters for protein sequence: {', '.join(invalid_chars)}"
                )
                return False
        elif sequence_type == "dna":
            valid_chars = set('ATGCNRYSWKMBDHV')
            invalid_chars = set(sequence) - valid_chars
            if invalid_chars:
                QMessageBox.warning(
                    self, "Invalid Sequence", 
                    f"Invalid characters for DNA sequence: {', '.join(invalid_chars)}"
                )
                return False
        elif sequence_type == "rna":
            valid_chars = set('AUGCNRYSWKMBDHV')
            invalid_chars = set(sequence) - valid_chars
            if invalid_chars:
                QMessageBox.warning(
                    self, "Invalid Sequence", 
                    f"Invalid characters for RNA sequence: {', '.join(invalid_chars)}"
                )
                return False
        
        return True
    
    def handle_analysis_error(self, error_message):
        """Handle analysis errors."""
        self.progress_bar.setVisible(False)
        QMessageBox.critical(self, "Analysis Error", error_message)
    
    def display_results(self, results):
        """Display analysis results."""
        self.progress_bar.setVisible(False)
        
        # Clear existing results
        for i in reversed(range(self.results_layout.count())):
            item = self.results_layout.itemAt(i)
            if item is not None:
                widget = item.widget()
                if widget is not None:
                    widget.setParent(None)
                else:
                    # Handle spacer items
                    self.results_layout.removeItem(item)
        
        sequence_type = self.sequence_type_combo.currentText()
        
        if sequence_type == "protein":
            self.display_protein_results(results)
        elif sequence_type in ["dna", "rna"]:
            self.display_nucleic_acid_results(results, sequence_type)
    
    def display_protein_results(self, results):
        """Display protein analysis results."""
        # Show note if using basic analysis
        if 'note' in results:
            note_label = QLabel(f"<i>{results['note']}</i>")
            note_label.setStyleSheet("color: #666; font-style: italic; padding: 5px;")
            self.results_layout.addWidget(note_label)
        
        # Basic properties
        basic_group = QGroupBox("Basic Properties")
        basic_layout = QFormLayout(basic_group)
        
        basic_layout.addRow("Length:", QLabel(f"{results['length']} amino acids"))
        basic_layout.addRow("Molecular Weight:", QLabel(f"{results['molecular_weight']:.2f} Da"))
        
        # Advanced properties (only if available)
        if 'isoelectric_point' in results:
            basic_layout.addRow("Isoelectric Point:", QLabel(f"{results['isoelectric_point']:.2f}"))
        if 'aromaticity' in results:
            basic_layout.addRow("Aromaticity:", QLabel(f"{results['aromaticity']:.3f}"))
        if 'instability_index' in results:
            basic_layout.addRow("Instability Index:", QLabel(f"{results['instability_index']:.2f}"))
        if 'gravy' in results:
            basic_layout.addRow("GRAVY (Hydropathy):", QLabel(f"{results['gravy']:.3f}"))
        
        self.results_layout.addWidget(basic_group)
        
        # Charge at different pH (only if available)
        if 'charge_at_pH' in results:
            charge_group = QGroupBox("Charge at Different pH")
            charge_layout = QFormLayout(charge_group)
            
            for ph, charge in results['charge_at_pH'].items():
                charge_layout.addRow(f"{ph.replace('_', ' ')}:", QLabel(f"{charge:.2f}"))
            
            self.results_layout.addWidget(charge_group)
        
        # Secondary structure (only if available)
        if 'secondary_structure' in results:
            ss_group = QGroupBox("Secondary Structure Prediction")
            ss_layout = QFormLayout(ss_group)
            
            ss_fractions = results['secondary_structure']
            ss_layout.addRow("Helix:", QLabel(f"{ss_fractions[0]:.1%}"))
            ss_layout.addRow("Turn:", QLabel(f"{ss_fractions[1]:.1%}"))
            ss_layout.addRow("Sheet:", QLabel(f"{ss_fractions[2]:.1%}"))
            
            self.results_layout.addWidget(ss_group)
        
        # Amino acid composition
        aa_group = QGroupBox("Amino Acid Composition")
        aa_layout = QVBoxLayout(aa_group)
        
        aa_table = QTableWidget()
        aa_table.setColumnCount(3)
        aa_table.setHorizontalHeaderLabels(["Amino Acid", "Count", "Percentage"])
        
        aa_percent = results['amino_acid_percent']
        aa_table.setRowCount(len(aa_percent))
        
        for i, (aa, percent) in enumerate(sorted(aa_percent.items())):
            count = int(percent * results['length'] / 100)
            aa_table.setItem(i, 0, QTableWidgetItem(aa))
            aa_table.setItem(i, 1, QTableWidgetItem(str(count)))
            aa_table.setItem(i, 2, QTableWidgetItem(f"{percent:.1f}%"))
        
        aa_table.resizeColumnsToContents()
        aa_table.setMaximumHeight(300)
        aa_layout.addWidget(aa_table)
        
        self.results_layout.addWidget(aa_group)
        
        self.results_layout.addStretch()
    
    def display_nucleic_acid_results(self, results, sequence_type):
        """Display DNA/RNA analysis results."""
        # Show note if using basic analysis
        if 'note' in results:
            note_label = QLabel(f"<i>{results['note']}</i>")
            note_label.setStyleSheet("color: #666; font-style: italic; padding: 5px;")
            self.results_layout.addWidget(note_label)
        
        # Basic properties
        basic_group = QGroupBox("Basic Properties")
        basic_layout = QFormLayout(basic_group)
        
        basic_layout.addRow("Length:", QLabel(f"{results['length']} nucleotides"))
        if 'molecular_weight' in results:
            basic_layout.addRow("Molecular Weight:", QLabel(f"{results['molecular_weight']:.2f} Da"))
        basic_layout.addRow("GC Content:", QLabel(f"{results['gc_content']:.1f}%"))
        
        self.results_layout.addWidget(basic_group)
        
        # Nucleotide composition
        comp_group = QGroupBox("Nucleotide Composition")
        comp_layout = QVBoxLayout(comp_group)
        
        comp_table = QTableWidget()
        comp_table.setColumnCount(3)
        comp_table.setHorizontalHeaderLabels(["Nucleotide", "Count", "Percentage"])
        
        nucleotide_counts = results['nucleotide_counts']
        nucleotide_percent = results['nucleotide_percent']
        comp_table.setRowCount(len(nucleotide_counts))
        
        # Define conventional nucleotide order
        if sequence_type == "dna":
            nucleotide_order = ['A', 'T', 'G', 'C', 'N']
        else:  # RNA
            nucleotide_order = ['A', 'U', 'G', 'C', 'N']
        
        # Display nucleotides in conventional order
        for i, nt in enumerate(nucleotide_order):
            if nt in nucleotide_counts:
                count = nucleotide_counts[nt]
                percent = nucleotide_percent[nt]
                comp_table.setItem(i, 0, QTableWidgetItem(nt))
                comp_table.setItem(i, 1, QTableWidgetItem(str(count)))
                comp_table.setItem(i, 2, QTableWidgetItem(f"{percent:.1f}%"))
        
        comp_table.resizeColumnsToContents()
        comp_table.setMaximumHeight(200)
        comp_layout.addWidget(comp_table)
        
        self.results_layout.addWidget(comp_group)
        
        # Complement sequences
        comp_seq_group = QGroupBox("Complement Sequences")
        comp_seq_layout = QVBoxLayout(comp_seq_group)
        
        complement_text = QTextEdit()
        complement_text.setMaximumHeight(100)
        complement_text.setFont(QFont("Courier", 9))
        complement_content = f"Complement: {results['complement'][:100]}{'...' if len(results['complement']) > 100 else ''}\n"
        complement_content += f"Reverse Complement: {results['reverse_complement'][:100]}{'...' if len(results['reverse_complement']) > 100 else ''}"
        complement_text.setPlainText(complement_content)
        complement_text.setReadOnly(True)
        comp_seq_layout.addWidget(complement_text)
        
        self.results_layout.addWidget(comp_seq_group)
        
        # Translation (for DNA/RNA)
        if sequence_type == "dna":
            trans_group = QGroupBox("Translation (All 6 Reading Frames)")
            trans_layout = QVBoxLayout(trans_group)
            
            # Add translation note if present
            if 'translation_note' in results:
                note_label = QLabel(f"<i>{results['translation_note']}</i>")
                note_label.setStyleSheet("color: #888; font-style: italic; padding: 5px;")
                note_label.setWordWrap(True)
                trans_layout.addWidget(note_label)
            
            trans_text = QTextEdit()
            trans_text.setMaximumHeight(150)
            trans_text.setFont(QFont("Courier", 9))
            
            trans_content = ""
            for frame, translation in results['translation_frames'].items():
                trans_content += f"{frame}: {translation[:50]}{'...' if len(translation) > 50 else ''}\n"
            
            trans_text.setPlainText(trans_content)
            trans_text.setReadOnly(True)
            trans_layout.addWidget(trans_text)
            
            self.results_layout.addWidget(trans_group)
        
        elif sequence_type == "rna":
            trans_group = QGroupBox("Translation")
            trans_layout = QVBoxLayout(trans_group)
            
            # Add translation note if present
            if 'translation_note' in results:
                note_label = QLabel(f"<i>{results['translation_note']}</i>")
                note_label.setStyleSheet("color: #888; font-style: italic; padding: 5px;")
                note_label.setWordWrap(True)
                trans_layout.addWidget(note_label)
            
            trans_text = QTextEdit()
            trans_text.setMaximumHeight(80)
            trans_text.setFont(QFont("Courier", 9))
            translation = results['translation']
            trans_text.setPlainText(f"Translation: {translation[:100]}{'...' if len(translation) > 100 else ''}")
            trans_text.setReadOnly(True)
            trans_layout.addWidget(trans_text)
            
            self.results_layout.addWidget(trans_group)
        
        self.results_layout.addStretch()


def create_bioinformatics_tab(parent=None):
    """Create the main bioinformatics tab with sub-tabs for different tools."""
    main_widget = QWidget()
    main_layout = QVBoxLayout(main_widget)
    main_layout.setContentsMargins(5, 5, 5, 5)
    
    # Create tab widget for different bioinformatics tools
    bio_tabs = QTabWidget()
    
    # Sequence Analysis tab
    seq_analysis_tab = SequenceAnalysisTab(parent)
    bio_tabs.addTab(seq_analysis_tab, "Sequence Analysis")
    
    # Placeholder tabs for future tools
    # Alignment tab
    alignment_tab = QWidget()
    alignment_layout = QVBoxLayout(alignment_tab)
    alignment_placeholder = QLabel(
        "<h3>Sequence Alignment Tools</h3>"
        "<p>Coming soon:</p>"
        "<ul>"
        "<li>Pairwise sequence alignment</li>"
        "<li>Multiple sequence alignment</li>"
        "<li>Alignment visualization</li>"
        "<li>Phylogenetic analysis</li>"
        "</ul>"
    )
    alignment_placeholder.setAlignment(Qt.AlignTop)
    alignment_placeholder.setWordWrap(True)
    alignment_placeholder.setStyleSheet("color: #666; padding: 20px;")
    alignment_layout.addWidget(alignment_placeholder)
    alignment_layout.addStretch()
    bio_tabs.addTab(alignment_tab, "Alignment")
    
    # Motif Finding tab
    motif_tab = QWidget()
    motif_layout = QVBoxLayout(motif_tab)
    motif_placeholder = QLabel(
        "<h3>Motif and Pattern Analysis</h3>"
        "<p>Coming soon:</p>"
        "<ul>"
        "<li>Motif discovery</li>"
        "<li>Pattern matching</li>"
        "<li>Regulatory element prediction</li>"
        "<li>Domain analysis</li>"
        "</ul>"
    )
    motif_placeholder.setAlignment(Qt.AlignTop)
    motif_placeholder.setWordWrap(True)
    motif_placeholder.setStyleSheet("color: #666; padding: 20px;")
    motif_layout.addWidget(motif_placeholder)
    motif_layout.addStretch()
    bio_tabs.addTab(motif_tab, "Motifs")
    
    # Structure Analysis tab
    structure_tab = QWidget()
    structure_layout = QVBoxLayout(structure_tab)
    structure_placeholder = QLabel(
        "<h3>Structure Analysis Tools</h3>"
        "<p>Coming soon:</p>"
        "<ul>"
        "<li>Secondary structure prediction</li>"
        "<li>Protein fold classification</li>"
        "<li>Active site prediction</li>"
        "<li>Structural comparison</li>"
        "</ul>"
    )
    structure_placeholder.setAlignment(Qt.AlignTop)
    structure_placeholder.setWordWrap(True)
    structure_placeholder.setStyleSheet("color: #666; padding: 20px;")
    structure_layout.addWidget(structure_placeholder)
    structure_layout.addStretch()
    bio_tabs.addTab(structure_tab, "Structure")
    
    main_layout.addWidget(bio_tabs)
    
    return main_widget


if __name__ == "__main__":
    # Test the bioinformatics tools
    import sys
    from PyQt5.QtWidgets import QApplication
    
    app = QApplication(sys.argv)
    
    widget = create_bioinformatics_tab()
    widget.setWindowTitle("PicoMol Bioinformatics Tools")
    widget.resize(800, 600)
    widget.show()
    
    sys.exit(app.exec_())