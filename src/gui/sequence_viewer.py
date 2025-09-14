from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QToolBar, QAction, 
                             QGraphicsView, QTextEdit, QListWidget, QSplitter, 
                             QFileDialog, QActionGroup, QGraphicsScene, 
                             QGraphicsRectItem, QGraphicsPathItem, QGraphicsTextItem,
                             QDialog, QMessageBox, QLabel, QPushButton, QComboBox,
                             QFormLayout, QGroupBox, QCheckBox, QSpinBox, QTabWidget,
                             QGraphicsPolygonItem)
from PyQt5.QtGui import QPen, QColor, QFont, QBrush, QPainterPath, QPainter, QPolygonF
from PyQt5.QtCore import Qt, QPointF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Restriction
from src.gui.enzyme_selection_dialog import EnzymeSelectionDialog
from src.gui.sequence_text_view import SequenceTextView
from src.gui.snapgene_sequence_view import SnapGeneSequenceView
from src.gui.restriction_cloning_dialog import RestrictionCloningDialog
from src.gui.gibson_assembly_dialog import GibsonAssemblyDialog
from src.gui.virtual_gel_dialog import VirtualGelDialog
from src.gui.golden_gate_dialog import GoldenGateDialog
from itertools import combinations
# Import new enhanced features
try:
    from src.gui.export_dialog import ExportDialog
    from src.gui.sequence_editor import SequenceEditor
    from src.gui.advanced_analysis import AdvancedAnalysisDialog
    ENHANCED_FEATURES_AVAILABLE = True
except ImportError as e:
    print(f"Enhanced features not available: {e}")
    ENHANCED_FEATURES_AVAILABLE = False
import math
import os

class SequenceViewer(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("SnapGene-style Sequence Viewer")
        self.feature_items = {}
        self.record = None
        self.parent_app = parent

        # Define a default set of restriction enzymes
        try:
            self.restriction_batch = RestrictionBatch(['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI'])
        except Exception as e:
            print(f"Warning: Could not initialize restriction enzymes: {e}")
            self.restriction_batch = None

        # Main layout
        main_layout = QVBoxLayout(self)
        
        # Create tab widget for different views
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Create main sequence viewer tab
        self.create_sequence_tab()
        
        # Create SnapGene-style sequence view tab
        self.create_snapgene_sequence_tab()
        
        # Create analysis tab
        self.create_analysis_tab()
        
        # Create tools tab
        self.create_tools_tab()




    def create_sequence_tab(self):
        """Create the main sequence viewer tab"""
        sequence_tab = QWidget()
        layout = QVBoxLayout(sequence_tab)
        
        # Toolbar
        toolbar = QToolBar()
        layout.addWidget(toolbar)

        self.load_action = QAction("📁 Load Sequence", self)
        self.load_action.setToolTip("Load a sequence file (GenBank, FASTA, etc.)")
        self.load_action.triggered.connect(self.load_sequence_file)
        toolbar.addAction(self.load_action)

        self.enzyme_action = QAction("🧬 Select Enzymes", self)
        self.enzyme_action.setToolTip("Choose restriction enzymes to display")
        self.enzyme_action.triggered.connect(self.select_enzymes)
        toolbar.addAction(self.enzyme_action)

        toolbar.addSeparator()
        
        # Enhanced features if available
        if ENHANCED_FEATURES_AVAILABLE:
            self.export_action = QAction("📤 Export", self)
            self.export_action.setToolTip("Export sequence maps and data")
            self.export_action.triggered.connect(self.open_export_dialog)
            toolbar.addAction(self.export_action)
            
            self.edit_action = QAction("✏️ Edit Sequence", self)
            self.edit_action.setToolTip("Edit sequence and features")
            self.edit_action.triggered.connect(self.open_sequence_editor)
            toolbar.addAction(self.edit_action)
            
            self.advanced_action = QAction("🔬 Advanced Analysis", self)
            self.advanced_action.setToolTip("Advanced sequence analysis tools")
            self.advanced_action.triggered.connect(self.open_advanced_analysis)
            toolbar.addAction(self.advanced_action)

        toolbar.addSeparator()

        self.circular_action = QAction("⭕ Circular View", self)
        self.circular_action.setCheckable(True)
        self.circular_action.setChecked(True)
        self.circular_action.setToolTip("Show plasmid in circular view")
        toolbar.addAction(self.circular_action)

        self.linear_action = QAction("📏 Linear View", self)
        self.linear_action.setCheckable(True)
        self.linear_action.setToolTip("Show sequence in linear view")
        toolbar.addAction(self.linear_action)

        self.circular_action.triggered.connect(self.update_view)
        self.linear_action.triggered.connect(self.update_view)

        # Action group for view switching
        self.view_action_group = QActionGroup(self)
        self.view_action_group.addAction(self.circular_action)
        self.view_action_group.addAction(self.linear_action)
        self.view_action_group.setExclusive(True)

        # Main content area
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        # Left panel: Feature list and info
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # Sequence info
        info_group = QGroupBox("Sequence Information")
        info_layout = QVBoxLayout(info_group)
        self.info_label = QLabel("No sequence loaded")
        self.info_label.setWordWrap(True)
        info_layout.addWidget(self.info_label)
        left_layout.addWidget(info_group)
        
        # Feature list
        feature_group = QGroupBox("Features")
        feature_layout = QVBoxLayout(feature_group)
        self.feature_list = QListWidget()
        self.feature_list.itemSelectionChanged.connect(self.highlight_feature)
        feature_layout.addWidget(self.feature_list)
        left_layout.addWidget(feature_group)
        
        splitter.addWidget(left_panel)

        # Center panel: Map and sequence views
        center_panel = QSplitter(Qt.Vertical)
        splitter.addWidget(center_panel)

        self.map_view = QGraphicsView()
        self.scene = QGraphicsScene()
        self.map_view.setScene(self.scene)
        self.map_view.setRenderHint(QPainter.Antialiasing)
        self.sequence_view = SequenceTextView()
        
        center_panel.addWidget(self.map_view)
        center_panel.addWidget(self.sequence_view)

        # Translation view
        translation_group = QGroupBox("Translation")
        translation_layout = QVBoxLayout(translation_group)
        self.translation_view = QTextEdit()
        self.translation_view.setReadOnly(True)
        self.translation_view.setMaximumHeight(120)
        self.translation_view.setFont(QFont("monospace", 10))
        translation_layout.addWidget(self.translation_view)
        center_panel.addWidget(translation_group)

        # Set initial sizes
        splitter.setSizes([300, 700])
        center_panel.setSizes([400, 300, 120])

        self.tab_widget.addTab(sequence_tab, "Sequence Map")
    
    def create_snapgene_sequence_tab(self):
        """Create the SnapGene-style sequence view tab"""
        self.snapgene_view = SnapGeneSequenceView()
        self.tab_widget.addTab(self.snapgene_view, "Sequence View")
        # Enable/disable tools based on whether sequence is loaded
        has_sequence = self.record is not None
        if hasattr(self, 'enzyme_action'):
            self.enzyme_action.setEnabled(has_sequence)
        if hasattr(self, 'cloning_button'):
            self.cloning_button.setEnabled(has_sequence)
        if hasattr(self, 'gibson_button'):
            self.gibson_button.setEnabled(has_sequence)
        if hasattr(self, 'golden_gate_button'):
            self.golden_gate_button.setEnabled(has_sequence)
        if hasattr(self, 'digest_button'):
            self.digest_button.setEnabled(has_sequence)
        if hasattr(self, 'analyze_button'):
            self.analyze_button.setEnabled(has_sequence)
        
        self.update_view()
    
    def create_analysis_tab(self):
        """Create the analysis tab"""
        analysis_tab = QWidget()
        layout = QVBoxLayout(analysis_tab)
        
        # Analysis tools
        tools_group = QGroupBox("Analysis Tools")
        tools_layout = QVBoxLayout(tools_group)
        
        # Restriction analysis
        restriction_group = QGroupBox("Restriction Analysis")
        restriction_layout = QFormLayout(restriction_group)
        
        self.min_cuts_spin = QSpinBox()
        self.min_cuts_spin.setMinimum(0)
        self.min_cuts_spin.setMaximum(10)
        self.min_cuts_spin.setValue(1)
        restriction_layout.addRow("Minimum cuts:", self.min_cuts_spin)
        
        self.max_cuts_spin = QSpinBox()
        self.max_cuts_spin.setMinimum(1)
        self.max_cuts_spin.setMaximum(20)
        self.max_cuts_spin.setValue(3)
        restriction_layout.addRow("Maximum cuts:", self.max_cuts_spin)
        
        self.analyze_button = QPushButton("🔍 Analyze Restriction Sites")
        self.analyze_button.clicked.connect(self.analyze_restriction_sites)
        restriction_layout.addRow(self.analyze_button)
        
        tools_layout.addWidget(restriction_group)
        
        # Results area
        self.analysis_results = QTextEdit()
        self.analysis_results.setReadOnly(True)
        tools_layout.addWidget(self.analysis_results)
        
        layout.addWidget(tools_group)
        
        self.tab_widget.addTab(analysis_tab, "Analysis")
    
    def create_tools_tab(self):
        """Create the molecular cloning tools tab"""
        tools_tab = QWidget()
        layout = QVBoxLayout(tools_tab)
        
        # Cloning tools
        cloning_group = QGroupBox("Molecular Cloning Tools")
        cloning_layout = QVBoxLayout(cloning_group)
        
        # Restriction cloning
        self.cloning_button = QPushButton("✂️ Restriction Cloning")
        self.cloning_button.setToolTip("Perform traditional restriction enzyme cloning")
        self.cloning_button.clicked.connect(self.open_cloning_dialog)
        cloning_layout.addWidget(self.cloning_button)
        
        # Gibson assembly
        self.gibson_button = QPushButton("🧬 Gibson Assembly")
        self.gibson_button.setToolTip("Perform Gibson assembly of DNA fragments")
        self.gibson_button.clicked.connect(self.open_gibson_dialog)
        cloning_layout.addWidget(self.gibson_button)
        
        # Golden Gate
        self.golden_gate_button = QPushButton("🔗 Golden Gate Assembly")
        self.golden_gate_button.setToolTip("Perform Golden Gate assembly using Type IIS enzymes")
        self.golden_gate_button.clicked.connect(self.open_golden_gate_dialog)
        cloning_layout.addWidget(self.golden_gate_button)
        
        # Virtual digest
        self.digest_button = QPushButton("⚡ Virtual Digest")
        self.digest_button.setToolTip("Simulate restriction enzyme digestion")
        self.digest_button.clicked.connect(self.virtual_digest)
        cloning_layout.addWidget(self.digest_button)
        
        layout.addWidget(cloning_group)
        
        # Enhanced tools if available
        if ENHANCED_FEATURES_AVAILABLE:
            enhanced_group = QGroupBox("Enhanced Tools")
            enhanced_layout = QVBoxLayout(enhanced_group)
            
            # Export tools
            export_button = QPushButton("📤 Export Maps & Sequences")
            export_button.setToolTip("Export high-quality maps and sequence data")
            export_button.clicked.connect(self.open_export_dialog)
            enhanced_layout.addWidget(export_button)
            
            # Sequence editor
            edit_button = QPushButton("✏️ Sequence Editor")
            edit_button.setToolTip("Edit sequence and manage features")
            edit_button.clicked.connect(self.open_sequence_editor)
            enhanced_layout.addWidget(edit_button)
            
            # Advanced analysis
            analysis_button = QPushButton("🔬 Advanced Analysis")
            analysis_button.setToolTip("ORF finding, primer design, sequence comparison")
            analysis_button.clicked.connect(self.open_advanced_analysis)
            enhanced_layout.addWidget(analysis_button)
            
            layout.addWidget(enhanced_group)
        
        # Add spacer
        layout.addStretch()
        
        # Instructions
        instructions = QLabel("""
        <h3>Molecular Cloning Tools</h3>
        <p>This tab provides tools for molecular cloning simulation:</p>
        <ul>
        <li><b>Restriction Cloning:</b> Traditional cut-and-paste cloning</li>
        <li><b>Gibson Assembly:</b> Seamless assembly of DNA fragments</li>
        <li><b>Golden Gate:</b> Type IIS enzyme-based assembly</li>
        <li><b>Virtual Digest:</b> Simulate restriction digestion</li>
        </ul>
        <p>Load a sequence file first to enable these tools.</p>
        """)
        instructions.setWordWrap(True)
        instructions.setStyleSheet("QLabel { background-color: #f0f0f0; padding: 10px; border-radius: 5px; }")
        layout.addWidget(instructions)
        
        self.tab_widget.addTab(tools_tab, "Cloning Tools")
    
    def analyze_restriction_sites(self):
        """Analyze restriction sites in the current sequence"""
        if not self.record or not self.restriction_batch:
            self.analysis_results.setText("No sequence loaded or restriction enzymes selected.")
            return
        
        try:
            analysis = self.restriction_batch.search(self.record.seq)
            
            min_cuts = self.min_cuts_spin.value()
            max_cuts = self.max_cuts_spin.value()
            
            results = []
            results.append(f"Restriction Analysis for {self.record.id}")
            results.append(f"Sequence length: {len(self.record.seq)} bp")
            results.append(f"Showing enzymes with {min_cuts}-{max_cuts} cut sites\n")
            
            suitable_enzymes = []
            
            for enzyme, sites in analysis.items():
                num_cuts = len(sites)
                if min_cuts <= num_cuts <= max_cuts:
                    suitable_enzymes.append((enzyme, sites, num_cuts))
            
            if suitable_enzymes:
                results.append("Suitable enzymes:")
                for enzyme, sites, num_cuts in sorted(suitable_enzymes, key=lambda x: x[2]):
                    results.append(f"  {enzyme}: {num_cuts} cuts at positions {sites}")
                    
                    # Calculate fragment sizes
                    if sites:
                        fragments = self.get_fragments(self.record.seq, sites)
                        fragment_sizes = [len(f) for f in fragments]
                        results.append(f"    Fragment sizes: {fragment_sizes} bp")
            else:
                results.append("No enzymes found with the specified number of cuts.")
            
            self.analysis_results.setText("\n".join(results))
            
        except Exception as e:
            self.analysis_results.setText(f"Error during analysis: {str(e)}")
    
    def open_golden_gate_dialog(self):
        """Open Golden Gate assembly dialog"""
        try:
            dialog = GoldenGateDialog(self)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening Golden Gate dialog: {str(e)}")

    def open_gibson_dialog(self):
        """Open Gibson assembly dialog"""
        try:
            dialog = GibsonAssemblyDialog(self)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening Gibson assembly dialog: {str(e)}")

    def virtual_digest(self):
        """Perform virtual restriction digest with enzyme combinations"""
        if not self.record:
            QMessageBox.warning(self, "No Sequence", "Please load a sequence first.")
            return
        
        if not self.restriction_batch:
            QMessageBox.warning(self, "No Enzymes", "Please select restriction enzymes first.")
            return

        try:
            analysis = self.restriction_batch.search(self.record.seq)
            all_fragments = []
            enzyme_fragments = {}
            
            # Get enzymes that actually cut the sequence
            cutting_enzymes = [(enzyme, sites) for enzyme, sites in analysis.items() if sites]
            
            if not cutting_enzymes:
                QMessageBox.information(self, "No Cuts", "Selected enzymes do not cut this sequence.")
                return
            
            # Add uncut sequence
            enzyme_fragments["Uncut"] = [self.record.seq]
            all_fragments.extend([self.record.seq])
            
            # Single enzyme digests
            for enzyme, sites in cutting_enzymes:
                fragments = self.get_fragments(self.record.seq, sites)
                all_fragments.extend(fragments)
                enzyme_fragments[str(enzyme)] = fragments
            
            # Double digests (combinations of 2 enzymes)
            if len(cutting_enzymes) >= 2:
                from itertools import combinations
                for enzyme_pair in combinations(cutting_enzymes, 2):
                    enzyme1, sites1 = enzyme_pair[0]
                    enzyme2, sites2 = enzyme_pair[1]
                    
                    # Combine cut sites from both enzymes
                    combined_sites = sorted(set(sites1 + sites2))
                    if combined_sites:
                        fragments = self.get_fragments(self.record.seq, combined_sites)
                        combo_name = f"{enzyme1} + {enzyme2}"
                        enzyme_fragments[combo_name] = fragments
                        all_fragments.extend(fragments)
            
            # Triple digests (combinations of 3 enzymes)
            if len(cutting_enzymes) >= 3:
                for enzyme_triple in combinations(cutting_enzymes, 3):
                    enzyme1, sites1 = enzyme_triple[0]
                    enzyme2, sites2 = enzyme_triple[1]
                    enzyme3, sites3 = enzyme_triple[2]
                    
                    # Combine cut sites from all three enzymes
                    combined_sites = sorted(set(sites1 + sites2 + sites3))
                    if combined_sites:
                        fragments = self.get_fragments(self.record.seq, combined_sites)
                        combo_name = f"{enzyme1} + {enzyme2} + {enzyme3}"
                        enzyme_fragments[combo_name] = fragments
                        all_fragments.extend(fragments)
            
            # Quadruple digests (combinations of 4 enzymes) - only if we have 4+ enzymes
            if len(cutting_enzymes) >= 4:
                for enzyme_quad in combinations(cutting_enzymes, 4):
                    enzymes = [e[0] for e in enzyme_quad]
                    all_sites = []
                    for _, sites in enzyme_quad:
                        all_sites.extend(sites)
                    
                    combined_sites = sorted(set(all_sites))
                    if combined_sites:
                        fragments = self.get_fragments(self.record.seq, combined_sites)
                        combo_name = f"{enzymes[0]} + {enzymes[1]} + {enzymes[2]} + {enzymes[3]}"
                        enzyme_fragments[combo_name] = fragments
                        all_fragments.extend(fragments)
            
            # Complete digest (all enzymes together)
            if len(cutting_enzymes) > 1:
                all_sites = []
                for _, sites in cutting_enzymes:
                    all_sites.extend(sites)
                
                combined_sites = sorted(set(all_sites))
                if combined_sites:
                    fragments = self.get_fragments(self.record.seq, combined_sites)
                    enzyme_names = [str(e[0]) for e in cutting_enzymes]
                    combo_name = "Complete digest (" + " + ".join(enzyme_names) + ")"
                    enzyme_fragments[combo_name] = fragments
                    all_fragments.extend(fragments)
            
            if all_fragments:
                dialog = VirtualGelDialog(all_fragments, self, enzyme_fragments=enzyme_fragments)
                dialog.exec_()
            else:
                QMessageBox.information(self, "No Cuts", "Selected enzymes do not cut this sequence.")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error during virtual digest: {str(e)}")

    def get_fragments(self, sequence, cuts):
        """Get DNA fragments from restriction cuts"""
        if not cuts:
            return [sequence]
        
        fragments = []
        sorted_cuts = sorted(cuts)
        
        # For circular sequences (plasmids), create fragments between consecutive cuts
        if len(sorted_cuts) == 1:
            # Single cut site: linearizes the plasmid into one fragment
            linearized = sequence[sorted_cuts[0]:] + sequence[:sorted_cuts[0]]
            return [linearized]
        else:
            # Multiple cut sites: create fragments between consecutive cuts
            for i in range(len(sorted_cuts)):
                start = sorted_cuts[i]
                end = sorted_cuts[(i + 1) % len(sorted_cuts)]  # Wrap around for circular
                
                if end > start:
                    # Normal fragment
                    fragments.append(sequence[start:end])
                else:
                    # Wraparound fragment (from last cut to first cut)
                    fragments.append(sequence[start:] + sequence[:end])
        
        return [f for f in fragments if len(f) > 0]

    def open_cloning_dialog(self):
        """Open restriction cloning dialog"""
        try:
            dialog = RestrictionCloningDialog(self)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening cloning dialog: {str(e)}")

    def select_enzymes(self):
        """Open enzyme selection dialog"""
        try:
            dialog = EnzymeSelectionDialog(self)
            if dialog.exec_() == QDialog.Accepted:
                selected_enzymes = dialog.get_selected_enzymes()
                if selected_enzymes:
                    self.restriction_batch = RestrictionBatch(selected_enzymes)
                    self.update_view()
                    
                    # Update sequence view highlighting
                    if hasattr(self, 'sequence_view'):
                        self.sequence_view.update_restriction_highlighting(self.restriction_batch)
                    
                    # Update SnapGene view highlighting
                    if hasattr(self, 'snapgene_view'):
                        self.snapgene_view.update_restriction_highlighting(self.restriction_batch)
                    
                    if hasattr(self, 'parent_app') and self.parent_app:
                        self.parent_app.statusBar().showMessage(f"Selected {len(selected_enzymes)} restriction enzymes", 3000)
                else:
                    QMessageBox.information(self, "No Selection", "No enzymes were selected.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error selecting enzymes: {str(e)}")

    def load_sequence_file(self):
        """Load a sequence file with improved error handling"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Sequence File", "", 
            "GenBank Files (*.gb *.gbk);;FASTA Files (*.fa *.fasta);;SnapGene Files (*.dna);;All Files (*)"
        )
        if not file_path:
            return

        try:
            # Determine file format and load accordingly
            file_ext = os.path.splitext(file_path)[1].lower()
            
            if file_ext in ['.gb', '.gbk', '.genbank']:
                record = SeqIO.read(file_path, "genbank")
            elif file_ext in ['.fa', '.fasta', '.fas']:
                record = SeqIO.read(file_path, "fasta")
            elif file_ext == '.dna':
                # SnapGene files - try to read as GenBank first
                try:
                    record = SeqIO.read(file_path, "genbank")
                except:
                    QMessageBox.warning(self, "Unsupported Format", 
                                       "SnapGene .dna files are not fully supported. Please export as GenBank or FASTA format.")
                    return
            else:
                # Try to auto-detect format
                try:
                    record = SeqIO.read(file_path, "genbank")
                except:
                    try:
                        record = SeqIO.read(file_path, "fasta")
                    except:
                        QMessageBox.critical(self, "Error", "Could not determine file format. Please use GenBank (.gb) or FASTA (.fa) format.")
                        return
            
            self.display_sequence(record)
            
            if hasattr(self, 'parent_app') and self.parent_app:
                self.parent_app.statusBar().showMessage(f"Loaded sequence: {record.id} ({len(record.seq)} bp)", 5000)
                
        except Exception as e:
            QMessageBox.critical(self, "Error Loading File", f"Error parsing file: {str(e)}")
            if hasattr(self, 'sequence_view'):
                self.sequence_view.setPlainText(f"Error parsing file: {e}")

    def display_sequence(self, record):
        """Display a sequence record with enhanced information"""
        self.record = record
        
        # Update sequence view
        if hasattr(self, 'sequence_view'):
            self.sequence_view.display_sequence(record, self.restriction_batch)
        
        # Update SnapGene sequence view
        if hasattr(self, 'snapgene_view'):
            self.snapgene_view.display_sequence(record, self.restriction_batch)
        
        # Update sequence information
        if hasattr(self, 'info_label'):
            info_text = f"<b>ID:</b> {record.id}<br>"
            info_text += f"<b>Length:</b> {len(record.seq)} bp<br>"
            if hasattr(record, 'description') and record.description:
                info_text += f"<b>Description:</b> {record.description}<br>"
            
            # Calculate GC content
            try:
                from Bio.SeqUtils import gc_fraction
                gc_content = gc_fraction(record.seq) * 100
                info_text += f"<b>GC Content:</b> {gc_content:.1f}%<br>"
            except:
                pass
            
            # Count features
            if hasattr(record, 'features'):
                feature_count = len([f for f in record.features if f.type != 'source'])
                info_text += f"<b>Features:</b> {feature_count}"
            
            self.info_label.setText(info_text)

        # Update feature list
        if hasattr(self, 'feature_list'):
            self.feature_list.clear()
            self.feature_items.clear()
            if hasattr(self.record, 'features'):
                for i, feature in enumerate(self.record.features):
                    if feature.type != 'source':  # Skip source features
                        # Get feature label
                        label = feature.qualifiers.get('label', [feature.type])[0] if 'label' in feature.qualifiers else feature.type
                        gene = feature.qualifiers.get('gene', [''])[0]
                        product = feature.qualifiers.get('product', [''])[0]
                        
                        display_name = label
                        if gene and gene != label:
                            display_name += f" ({gene})"
                        elif product and product != label:
                            display_name += f" ({product})"
                        
                        item_text = f"{feature.type}: {display_name} [{feature.location.start}:{feature.location.end}]"
                        self.feature_list.addItem(item_text)

        self.update_view()

    def update_view(self):
        """Update the sequence view based on current settings"""
        if hasattr(self, 'circular_action') and self.circular_action.isChecked():
            self.show_circular_view()
        else:
            self.show_linear_view()
        
        # Fit the view to show all content
        if hasattr(self, 'scene') and hasattr(self, 'map_view'):
            self.map_view.fitInView(self.scene.itemsBoundingRect(), Qt.KeepAspectRatio)

    def show_circular_view(self):
        """Display sequence in circular plasmid view"""
        self.scene.clear()
        
        if not self.record:
            radius = 150
            self.scene.addEllipse(-radius, -radius, radius * 2, radius * 2, QPen(Qt.gray, 2))
            # Add placeholder text
            text = self.scene.addText("Load a sequence\nto view plasmid map", QFont("Arial", 12))
            text.setPos(-60, -10)
            return

        plasmid_len = len(self.record.seq)
        radius = 150
        center_x, center_y = 0, 0

        # Draw plasmid backbone
        pen = QPen(Qt.black, 4)
        self.scene.addEllipse(center_x - radius, center_y - radius, radius * 2, radius * 2, pen)

        # Draw features
        self.feature_items.clear()
        if hasattr(self.record, 'features'):
            feature_index = 0
            for i, feature in enumerate(self.record.features):
                if feature.type in ["gene", "promoter", "CDS", "terminator", "rep_origin", "misc_feature"]:
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    
                    # Handle wraparound features
                    if start > end:
                        end += plasmid_len
                    
                    start_angle = (start / plasmid_len) * 360
                    span_angle = ((end - start) / plasmid_len) * 360

                    pen_color = self.get_feature_color(feature.type)
                    
                    pen = QPen(pen_color, 12)
                    path = QPainterPath()
                    path.arcMoveTo(center_x - radius, center_y - radius, radius * 2, radius * 2, -start_angle)
                    path.arcTo(center_x - radius, center_y - radius, radius * 2, radius * 2, -start_angle, -span_angle)
                    arc = self.scene.addPath(path, pen)
                    self.feature_items[feature_index] = arc
                    
                    # Add directional arrow for the feature
                    self.add_circular_arrow(feature, center_x, center_y, radius, start_angle, span_angle, pen_color)
                    
                    feature_index += 1

                    # Add feature label
                    if span_angle > 10:  # Only label features that are large enough
                        label_angle = start_angle + (span_angle / 2)
                        rad_angle = math.radians(label_angle)
                        label_radius = radius + 25
                        label_x = center_x + label_radius * math.cos(rad_angle)
                        label_y = center_y + label_radius * math.sin(rad_angle)
                        
                        # Get feature name
                        feature_name = feature.qualifiers.get('label', [feature.type])[0]
                        if 'gene' in feature.qualifiers:
                            feature_name = feature.qualifiers['gene'][0]
                        
                        label = self.scene.addText(feature_name, QFont("Arial", 8))
                        label.setPos(label_x - 20, label_y - 5)
                        if abs(label_angle) > 90 and abs(label_angle) < 270:
                            label.setRotation(label_angle + 180)
                        else:
                            label.setRotation(label_angle)

        # Draw restriction sites with arrows
        if self.restriction_batch:
            try:
                analysis = self.restriction_batch.search(self.record.seq)
                for enzyme, sites in analysis.items():
                    for site in sites:
                        angle = (site / plasmid_len) * 360
                        rad_angle = math.radians(angle)
                        
                        x1 = center_x + (radius - 8) * math.cos(rad_angle)
                        y1 = center_y + (radius - 8) * math.sin(rad_angle)
                        x2 = center_x + (radius + 8) * math.cos(rad_angle)
                        y2 = center_y + (radius + 8) * math.sin(rad_angle)
                        
                        self.scene.addLine(x1, y1, x2, y2, QPen(Qt.darkRed, 2))
                        
                        # Add restriction cut arrow
                        self.add_circular_restriction_arrow(site, enzyme)
                        
                        # Add enzyme label with cut position indicator
                        if len(sites) <= 5:  # Show labels for enzymes with few cut sites
                            # Get cut position information
                            cut_info = ""
                            try:
                                recognition_site = str(enzyme.site)
                                # Use fst5 for 5' cut position (1-based, convert to 0-based)
                                cut_pos = getattr(enzyme, 'fst5', None)
                                if cut_pos is not None and cut_pos > 0:
                                    cut_pos_0based = cut_pos - 1
                                    if 0 <= cut_pos_0based < len(recognition_site):
                                        # Show cut position within recognition site
                                        cut_site_display = recognition_site[:cut_pos_0based] + "↓" + recognition_site[cut_pos_0based:]
                                        cut_info = f" ({cut_site_display})"
                            except:
                                pass
                            
                            label_text = f"{str(enzyme)}{cut_info}"
                            label = self.scene.addText(label_text, QFont("Arial", 6))
                            label.setDefaultTextColor(Qt.darkRed)
                            label_x = center_x + (radius + 20) * math.cos(rad_angle)
                            label_y = center_y + (radius + 20) * math.sin(rad_angle)
                            label.setPos(label_x - 15, label_y - 8)
            except Exception as e:
                print(f"Error drawing restriction sites: {e}")
        
        # Add sequence info in center
        info_text = f"{self.record.id}\n{len(self.record.seq)} bp"
        center_label = self.scene.addText(info_text, QFont("Arial", 10))
        center_label.setPos(-30, -10)


    def show_linear_view(self):
        """Display sequence in linear view"""
        self.scene.clear()

        if not self.record:
            self.scene.addLine(0, 0, 800, 0, QPen(Qt.gray, 2))
            text = self.scene.addText("Load a sequence to view linear map", QFont("Arial", 12))
            text.setPos(200, -30)
            return

        plasmid_len = len(self.record.seq)
        scale_factor = min(1.0, 1000.0 / plasmid_len) if plasmid_len > 1000 else 1.0
        scaled_len = plasmid_len * scale_factor
        
        # Draw main sequence line
        self.scene.addLine(0, 0, scaled_len, 0, QPen(Qt.black, 4))

        # Draw features
        self.feature_items.clear()
        if hasattr(self.record, 'features'):
            feature_index = 0
            y_offset = 0
            for i, feature in enumerate(self.record.features):
                if feature.type in ["gene", "promoter", "CDS", "terminator", "rep_origin", "misc_feature"]:
                    start = int(feature.location.start) * scale_factor
                    end = int(feature.location.end) * scale_factor
                    width = end - start
                    
                    if width < 2:  # Minimum width for visibility
                        width = 2
                        end = start + width

                    color = self.get_feature_color(feature.type)
                    
                    # Alternate feature positions to avoid overlap
                    y_pos = -20 - (feature_index % 3) * 25
                    
                    rect = self.scene.addRect(start, y_pos, width, 15, QPen(color, 1), QBrush(color))
                    self.feature_items[feature_index] = rect
                    
                    # Add directional arrow for the feature
                    self.add_linear_arrow(feature, start, end, y_pos, color)
                    
                    feature_index += 1

                    # Add feature label
                    if width > 20:  # Only label if feature is wide enough
                        feature_name = feature.qualifiers.get('label', [feature.type])[0]
                        if 'gene' in feature.qualifiers:
                            feature_name = feature.qualifiers['gene'][0]
                        
                        label = self.scene.addText(feature_name, QFont("Arial", 8))
                        label.setPos(start, y_pos - 20)

        # Draw restriction sites with arrows
        if self.restriction_batch:
            try:
                analysis = self.restriction_batch.search(self.record.seq)
                for enzyme, sites in analysis.items():
                    for site in sites:
                        scaled_site = site * scale_factor
                        self.scene.addLine(scaled_site, -30, scaled_site, 30, QPen(Qt.darkRed, 2))
                        
                        # Add restriction cut arrow
                        self.add_linear_restriction_arrow(site, enzyme)
                        
                        # Add enzyme label with cut position indicator
                        try:
                            recognition_site = str(enzyme.site)
                            # Use fst5 for 5' cut position (1-based, convert to 0-based)
                            cut_pos = getattr(enzyme, 'fst5', None)
                            if cut_pos is not None and cut_pos > 0:
                                cut_pos_0based = cut_pos - 1
                                if 0 <= cut_pos_0based < len(recognition_site):
                                    # Show cut position within recognition site
                                    cut_site_display = recognition_site[:cut_pos_0based] + "↓" + recognition_site[cut_pos_0based:]
                                    label_text = f"{str(enzyme)} ({cut_site_display})"
                                else:
                                    label_text = str(enzyme)
                            else:
                                label_text = str(enzyme)
                        except:
                            label_text = str(enzyme)
                        
                        label = self.scene.addText(label_text, QFont("Arial", 6))
                        label.setDefaultTextColor(Qt.darkRed)
                        label.setPos(scaled_site - 20, 35)
            except Exception as e:
                print(f"Error drawing restriction sites: {e}")

        # Add scale information
        scale_info = f"Scale: 1 pixel = {1/scale_factor:.0f} bp" if scale_factor < 1 else "Scale: 1:1"
        scale_label = self.scene.addText(scale_info, QFont("Arial", 8))
        scale_label.setPos(0, 60)
        
        # Add sequence info
        info_text = f"{self.record.id} - {len(self.record.seq)} bp"
        info_label = self.scene.addText(info_text, QFont("Arial", 10))
        info_label.setPos(0, -60)
        
        # Add directional indicator for linear sequence
        self.add_linear_direction_indicator(scaled_len)


    def highlight_feature(self):
        """Highlight selected features and show translation"""
        if not self.record or not hasattr(self.record, 'features'):
            return
            
        # Reset all feature colors first
        for i, item in self.feature_items.items():
            if isinstance(item, QGraphicsRectItem):
                # Find the corresponding feature
                feature_index = 0
                for j, feature in enumerate(self.record.features):
                    if feature.type in ["gene", "promoter", "CDS", "terminator", "rep_origin", "misc_feature"]:
                        if feature_index == i:
                            color = self.get_feature_color(feature.type)
                            item.setBrush(QBrush(color))
                            item.setPen(QPen(color))
                            break
                        feature_index += 1
            elif isinstance(item, QGraphicsPathItem):
                # Find the corresponding feature for arcs
                feature_index = 0
                for j, feature in enumerate(self.record.features):
                    if feature.type in ["gene", "promoter", "CDS", "terminator", "rep_origin", "misc_feature"]:
                        if feature_index == i:
                            pen = item.pen()
                            pen.setColor(self.get_feature_color(feature.type))
                            pen.setWidth(12)
                            item.setPen(pen)
                            break
                        feature_index += 1
        
        # Highlight selected features
        selected_items = self.feature_list.selectedItems()
        if selected_items:
            for selected_item in selected_items:
                row = self.feature_list.row(selected_item)
                if row in self.feature_items:
                    item = self.feature_items[row]
                    if isinstance(item, QGraphicsRectItem):
                        item.setBrush(QBrush(Qt.yellow))
                        item.setPen(QPen(Qt.black, 3))
                    elif isinstance(item, QGraphicsPathItem):
                        pen = item.pen()
                        pen.setColor(Qt.yellow)
                        pen.setWidth(15)
                        item.setPen(pen)

        # Display translation for the first selected CDS feature
        self.translation_view.clear()
        if not selected_items:
            return
            
        # Find the first CDS feature in selection
        for selected_item in selected_items:
            row = self.feature_list.row(selected_item)
            # Map row to actual feature
            feature_index = 0
            for feature in self.record.features:
                if feature.type != 'source':
                    if feature_index == row:
                        if feature.type == "CDS":
                            self.display_translation(feature)
                            return
                    feature_index += 1

    def display_translation(self, feature):
        """Display protein translation for a CDS feature"""
        try:
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand

            seq = self.record.seq[start:end]
            if strand == -1:
                seq = seq.reverse_complement()
            
            # Try to get the translation from the feature first
            if 'translation' in feature.qualifiers:
                protein = feature.qualifiers['translation'][0]
                translation_text = f"Translation (from annotation):\n{protein}"
            else:
                # Translate the sequence
                protein = seq.translate()
                translation_text = f"Translation (computed):\n{protein}"
            
            # Add additional information
            translation_text += f"\n\nFeature: {feature.type}"
            translation_text += f"\nLocation: {start}..{end}"
            translation_text += f"\nStrand: {'+' if strand == 1 else '-'}"
            translation_text += f"\nLength: {len(protein)} amino acids"
            
            self.translation_view.setText(translation_text)
            
        except Exception as e:
            self.translation_view.setText(f"Error translating sequence: {str(e)}")

    def get_feature_color(self, feature_type):
        """Get color for different feature types"""
        colors = {
            "gene": QColor("#3498db"),        # Blue
            "CDS": QColor("#2ecc71"),         # Green  
            "promoter": QColor("#e74c3c"),     # Red
            "terminator": QColor("#f39c12"),   # Orange
            "rep_origin": QColor("#9b59b6"),   # Purple
            "misc_feature": QColor("#95a5a6"), # Gray
            "regulatory": QColor("#e67e22"),   # Dark orange
            "enhancer": QColor("#f1c40f"),     # Yellow
            "primer_bind": QColor("#1abc9c"),  # Turquoise
            "protein_bind": QColor("#34495e") # Dark gray
        }
        return colors.get(feature_type, QColor("#bdc3c7"))  # Light gray for others
    
    def add_circular_arrow(self, feature, center_x, center_y, radius, start_angle, span_angle, color):
        """Add directional indicator pointing around the plasmid circle"""
        try:
            # Get feature strand information
            strand = getattr(feature.location, 'strand', 1)
            if strand is None:
                strand = 1
            
            # Only add arrows for features large enough and with strand info
            if span_angle < 15 or strand == 0:  # Only for features spanning more than 15 degrees
                return
            
            # Calculate arrow position (towards the end of the feature)
            if strand >= 0:  # Forward strand (clockwise)
                arrow_angle = start_angle + (span_angle * 0.8)  # Near the end
            else:  # Reverse strand (counterclockwise)
                arrow_angle = start_angle + (span_angle * 0.2)  # Near the beginning
                
            rad_angle = math.radians(arrow_angle)
            
            # Arrow positioned on the feature arc
            arrow_radius = radius
            arrow_size = 6  # Slightly larger for better visibility
            
            # Arrow tip position (on the circle)
            tip_x = center_x + arrow_radius * math.cos(rad_angle)
            tip_y = center_y + arrow_radius * math.sin(rad_angle)
            
            # Create arrow pointing tangentially around the circle
            if strand >= 0:  # Forward strand (clockwise direction)
                # Arrow pointing in clockwise direction (tangent to circle)
                # Tangent direction is perpendicular to radius
                tangent_angle = rad_angle + math.radians(90)  # 90 degrees ahead of radius
                
                # Arrow base points (behind the tip in tangent direction)
                base_angle1 = tangent_angle + math.radians(150)  # 150 degrees back
                base_angle2 = tangent_angle + math.radians(210)  # 210 degrees back
                
                base1_x = tip_x + arrow_size * math.cos(base_angle1)
                base1_y = tip_y + arrow_size * math.sin(base_angle1)
                
                base2_x = tip_x + arrow_size * math.cos(base_angle2)
                base2_y = tip_y + arrow_size * math.sin(base_angle2)
                
            else:  # Reverse strand (counterclockwise direction)
                # Arrow pointing in counterclockwise direction
                tangent_angle = rad_angle - math.radians(90)  # 90 degrees behind radius
                
                # Arrow base points (behind the tip in tangent direction)
                base_angle1 = tangent_angle - math.radians(150)  # 150 degrees back
                base_angle2 = tangent_angle - math.radians(210)  # 210 degrees back
                
                base1_x = tip_x + arrow_size * math.cos(base_angle1)
                base1_y = tip_y + arrow_size * math.sin(base_angle1)
                
                base2_x = tip_x + arrow_size * math.cos(base_angle2)
                base2_y = tip_y + arrow_size * math.sin(base_angle2)
            
            arrow_points = [QPointF(tip_x, tip_y), QPointF(base1_x, base1_y), QPointF(base2_x, base2_y)]
            
            # Create arrow with contrasting color for visibility
            arrow_color = color.darker(150)  # Darker for better contrast
            arrow_polygon = QPolygonF(arrow_points)
            arrow_item = self.scene.addPolygon(arrow_polygon, QPen(arrow_color, 1), QBrush(arrow_color))
                
        except Exception as e:
            print(f"Error adding circular arrow: {e}")
    
    def add_linear_arrow(self, feature, start, end, y_pos, color):
        """Add elegant directional indicator to linear feature"""
        try:
            # Get feature strand information
            strand = getattr(feature.location, 'strand', 1)
            if strand is None:
                strand = 1
            
            width = end - start
            
            # Only add arrow if feature is wide enough and has strand info
            if width < 20 or strand == 0:
                return
            
            # Create a subtle chevron-style arrow
            arrow_height = 3
            arrow_width = min(6, width / 6)
            
            # Arrow position (slightly above center)
            arrow_y = y_pos + 5
            
            # Create chevron-style arrow
            if strand >= 0:  # Forward strand (left to right)
                # Chevron pointing right
                arrow_x = end - arrow_width - 2
                
                # Two lines forming a chevron
                line1_start = QPointF(arrow_x, arrow_y - arrow_height)
                line1_end = QPointF(arrow_x + arrow_width, arrow_y)
                
                line2_start = QPointF(arrow_x, arrow_y + arrow_height)
                line2_end = QPointF(arrow_x + arrow_width, arrow_y)
                
                # Draw chevron lines
                arrow_color = color.darker(130)
                self.scene.addLine(line1_start.x(), line1_start.y(), 
                                 line1_end.x(), line1_end.y(), 
                                 QPen(arrow_color, 2))
                self.scene.addLine(line2_start.x(), line2_start.y(), 
                                 line2_end.x(), line2_end.y(), 
                                 QPen(arrow_color, 2))
                
            else:  # Reverse strand (right to left)
                # Chevron pointing left
                arrow_x = start + arrow_width + 2
                
                # Two lines forming a chevron
                line1_start = QPointF(arrow_x, arrow_y - arrow_height)
                line1_end = QPointF(arrow_x - arrow_width, arrow_y)
                
                line2_start = QPointF(arrow_x, arrow_y + arrow_height)
                line2_end = QPointF(arrow_x - arrow_width, arrow_y)
                
                # Draw chevron lines
                arrow_color = color.darker(130)
                self.scene.addLine(line1_start.x(), line1_start.y(), 
                                 line1_end.x(), line1_end.y(), 
                                 QPen(arrow_color, 2))
                self.scene.addLine(line2_start.x(), line2_start.y(), 
                                 line2_end.x(), line2_end.y(), 
                                 QPen(arrow_color, 2))
            
        except Exception as e:
            print(f"Error adding linear arrow: {e}")
    
    def add_restriction_site_arrows(self):
        """Add arrows to show restriction site cut directions"""
        if not self.restriction_batch or not self.record:
            return
            
        try:
            analysis = self.restriction_batch.search(self.record.seq)
            
            for enzyme, sites in analysis.items():
                for site in sites:
                    if hasattr(self, 'circular_action') and self.circular_action.isChecked():
                        self.add_circular_restriction_arrow(site, enzyme)
                    else:
                        self.add_linear_restriction_arrow(site, enzyme)
                        
        except Exception as e:
            print(f"Error adding restriction arrows: {e}")
    
    def add_circular_restriction_arrow(self, site, enzyme):
        """Add restriction cut indicator in circular view"""
        try:
            plasmid_len = len(self.record.seq)
            radius = 150
            center_x, center_y = 0, 0
            
            angle = (site / plasmid_len) * 360
            rad_angle = math.radians(angle)
            
            # Smaller cut indicator - tiny triangle pointing radially inward
            arrow_size = 1.5  # Reduced from 3
            cut_radius = radius + 4  # Reduced from 6
            
            # Triangle pointing toward center (correct for restriction cuts)
            tip_x = center_x + (radius - 1) * math.cos(rad_angle)  # Reduced from -2
            tip_y = center_y + (radius - 1) * math.sin(rad_angle)
            
            # Base points with smaller angle spread
            base_angle1 = rad_angle + math.radians(5)  # Reduced from 8
            base_angle2 = rad_angle - math.radians(5)  # Reduced from 8
            
            base1_x = center_x + cut_radius * math.cos(base_angle1)
            base1_y = center_y + cut_radius * math.sin(base_angle1)
            
            base2_x = center_x + cut_radius * math.cos(base_angle2)
            base2_y = center_y + cut_radius * math.sin(base_angle2)
            
            # Create small triangle
            arrow_points = [
                QPointF(tip_x, tip_y),
                QPointF(base1_x, base1_y),
                QPointF(base2_x, base2_y)
            ]
            
            arrow_polygon = QPolygonF(arrow_points)
            arrow_item = self.scene.addPolygon(arrow_polygon, 
                                             QPen(Qt.darkRed, 1), 
                                             QBrush(QColor(200, 50, 50, 120)))  # More transparent
            
        except Exception as e:
            print(f"Error adding circular restriction arrow: {e}")
    
    def add_linear_restriction_arrow(self, site, enzyme):
        """Add subtle restriction cut indicator in linear view"""
        try:
            plasmid_len = len(self.record.seq)
            scale_factor = min(1.0, 1000.0 / plasmid_len) if plasmid_len > 1000 else 1.0
            scaled_site = site * scale_factor
            
            # Smaller cut indicator - tiny downward triangle
            arrow_size = 1.5  # Reduced from 2.5
            arrow_y_top = -30  # Moved closer to sequence line
            arrow_y_bottom = -27  # Moved closer to sequence line
            
            # Small triangle pointing down
            arrow_points = [
                QPointF(scaled_site, arrow_y_bottom),           # Tip (bottom)
                QPointF(scaled_site - arrow_size, arrow_y_top), # Left wing
                QPointF(scaled_site + arrow_size, arrow_y_top)  # Right wing
            ]
            
            arrow_polygon = QPolygonF(arrow_points)
            arrow_item = self.scene.addPolygon(arrow_polygon, 
                                             QPen(Qt.darkRed, 1), 
                                             QBrush(QColor(200, 50, 50, 120)))  # More transparent
            
        except Exception as e:
            print(f"Error adding linear restriction arrow: {e}")
    
    def add_linear_direction_indicator(self, sequence_length):
        """Add elegant 5' to 3' direction indicator for linear view"""
        try:
            # Add 5' and 3' labels with better styling
            label_5prime = self.scene.addText("5'", QFont("Arial", 10, QFont.Bold))
            label_5prime.setDefaultTextColor(QColor(70, 130, 180))  # Steel blue
            label_5prime.setPos(-25, -8)
            
            label_3prime = self.scene.addText("3'", QFont("Arial", 10, QFont.Bold))
            label_3prime.setDefaultTextColor(QColor(70, 130, 180))
            label_3prime.setPos(sequence_length + 8, -8)
            
            # Add subtle directional indicator
            if sequence_length > 100:  # Only for longer sequences
                arrow_y = 25
                arrow_start = sequence_length * 0.2
                arrow_end = sequence_length * 0.8
                
                # Subtle dashed line
                pen = QPen(QColor(70, 130, 180, 100), 1)  # Semi-transparent
                pen.setStyle(Qt.DashLine)
                self.scene.addLine(arrow_start, arrow_y, arrow_end, arrow_y, pen)
                
                # Small arrow head
                arrow_head_size = 4
                arrow_points = [
                    QPointF(arrow_end, arrow_y),
                    QPointF(arrow_end - arrow_head_size, arrow_y - 2),
                    QPointF(arrow_end - arrow_head_size, arrow_y + 2)
                ]
                
                arrow_polygon = QPolygonF(arrow_points)
                arrow_item = self.scene.addPolygon(arrow_polygon, 
                                                 QPen(QColor(70, 130, 180), 1), 
                                                 QBrush(QColor(70, 130, 180, 150)))
            
        except Exception as e:
            print(f"Error adding direction indicator: {e}")
    
    # Enhanced feature methods
    def open_export_dialog(self):
        """Open the enhanced export dialog"""
        if not ENHANCED_FEATURES_AVAILABLE:
            QMessageBox.warning(self, "Feature Unavailable", "Enhanced export features are not available.")
            return
        
        if not self.record:
            QMessageBox.warning(self, "No Sequence", "Please load a sequence first.")
            return
        
        try:
            dialog = ExportDialog(self)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening export dialog: {str(e)}")
    
    def open_sequence_editor(self):
        """Open the sequence editor"""
        if not ENHANCED_FEATURES_AVAILABLE:
            QMessageBox.warning(self, "Feature Unavailable", "Sequence editing features are not available.")
            return
        
        if not self.record:
            QMessageBox.warning(self, "No Sequence", "Please load a sequence first.")
            return
        
        try:
            dialog = SequenceEditor(self.record, self)
            dialog.sequenceModified.connect(self.on_sequence_modified)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening sequence editor: {str(e)}")
    
    def open_advanced_analysis(self):
        """Open advanced analysis dialog"""
        if not ENHANCED_FEATURES_AVAILABLE:
            QMessageBox.warning(self, "Feature Unavailable", "Advanced analysis features are not available.")
            return
        
        if not self.record:
            QMessageBox.warning(self, "No Sequence", "Please load a sequence first.")
            return
        
        try:
            dialog = AdvancedAnalysisDialog(self)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening advanced analysis: {str(e)}")
    
    def on_sequence_modified(self, modified_record):
        """Handle sequence modifications from the editor"""
        if modified_record:
            self.display_sequence(modified_record)
            if hasattr(self, 'parent_app') and self.parent_app:
                self.parent_app.statusBar().showMessage("Sequence modified and updated", 3000)
