from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QToolBar, 
                             QGraphicsView, QTextEdit, QListWidget, QSplitter, 
                             QFileDialog, QGraphicsScene, 
                             QGraphicsRectItem, QGraphicsPathItem, QGraphicsTextItem,
                             QGraphicsItem, QGraphicsPolygonItem, QListWidgetItem,
                             QDialog, QMessageBox, QLabel, QPushButton, QComboBox,
                             QFormLayout, QGroupBox, QCheckBox, QSpinBox, QTabWidget)
# QAction may be located in QtWidgets or QtGui depending on PyQt6 build
try:
    from PyQt6.QtWidgets import QAction
except ImportError:
    from PyQt6.QtGui import QAction
from PyQt6.QtGui import QPen, QColor, QFont, QBrush, QPainterPath, QPainter, QPolygonF, QActionGroup
from PyQt6.QtCore import Qt, QPointF
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
from itertools import combinations, product
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
            default_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SacI', 'KpnI', 'NotI', 
                              'XbaI', 'SpeI', 'PstI', 'SalI', 'SmaI', 'BglII', 'NcoI', 
                              'NdeI', 'ApaI', 'EcoRV', 'NruI', 'SphI', 'ClaI']
            self.restriction_batch = RestrictionBatch(default_enzymes)
            self.current_enzyme_names = default_enzymes  # Store the names for easy access
        except Exception as e:
            print(f"Warning: Could not initialize restriction enzymes: {e}")
            self.restriction_batch = None
            self.current_enzyme_names = []

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

        toolbar.addSeparator()
        
        # Zoom controls
        self.zoom_in_action = QAction("🔍+", self)
        self.zoom_in_action.setToolTip("Zoom in (Ctrl++)")
        self.zoom_in_action.triggered.connect(self.zoom_in)
        toolbar.addAction(self.zoom_in_action)
        
        self.zoom_out_action = QAction("🔍-", self)
        self.zoom_out_action.setToolTip("Zoom out (Ctrl+-)")
        self.zoom_out_action.triggered.connect(self.zoom_out)
        toolbar.addAction(self.zoom_out_action)
        
        self.reset_zoom_action = QAction("🔄 Reset Zoom", self)
        self.reset_zoom_action.setToolTip("Reset zoom to fit view (Ctrl+0)")
        self.reset_zoom_action.triggered.connect(self.reset_zoom)
        toolbar.addAction(self.reset_zoom_action)
        
        toolbar.addSeparator()

        # View controls
        self.circular_action = QAction("⭕ Circular View", self)
        self.circular_action.setToolTip("Show circular plasmid map")
        self.circular_action.setCheckable(True)
        self.circular_action.setChecked(True)
        self.circular_action.triggered.connect(self.toggle_circular_view)
        toolbar.addAction(self.circular_action)
        
        self.linear_action = QAction("📏 Linear View", self)
        self.linear_action.setToolTip("Show linear sequence map")
        self.linear_action.setCheckable(True)
        self.linear_action.setChecked(False)
        self.linear_action.triggered.connect(self.toggle_linear_view)
        toolbar.addAction(self.linear_action)
        
        # Create button group for mutual exclusion
        self.view_group = QActionGroup(self)
        self.view_group.addAction(self.circular_action)
        self.view_group.addAction(self.linear_action)
        self.view_group.setExclusive(True)
        
        toolbar.addSeparator()

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
        splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(splitter)

        # Left panel: Feature list and info
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # Sequence info
        info_group = QGroupBox("Sequence Information")
        info_group.setStyleSheet("""
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
        info_layout = QVBoxLayout(info_group)
        info_layout.setContentsMargins(10, 15, 10, 10)
        info_layout.setSpacing(8)
        self.info_label = QLabel("No sequence loaded")
        self.info_label.setWordWrap(True)
        info_layout.addWidget(self.info_label)
        left_layout.addWidget(info_group)
        
        # Feature list
        feature_group = QGroupBox("Features")
        feature_group.setStyleSheet("""
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
        feature_layout = QVBoxLayout(feature_group)
        feature_layout.setContentsMargins(10, 15, 10, 10)
        feature_layout.setSpacing(8)
        self.feature_list = QListWidget()
        self.feature_list.itemSelectionChanged.connect(self.on_feature_selection_changed)
        feature_layout.addWidget(self.feature_list)
        left_layout.addWidget(feature_group)
        
        splitter.addWidget(left_panel)

        # Center panel: Map and sequence views
        center_panel = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(center_panel)

        self.map_view = QGraphicsView()
        self.scene = QGraphicsScene()
        self.map_view.setScene(self.scene)
        self.map_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.map_view.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        
        # Enhanced map view features
        self.map_view.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
        self.map_view.setRubberBandSelectionMode(Qt.ItemSelectionMode.ContainsItemShape)
        self.map_view.setOptimizationFlag(QGraphicsView.OptimizationFlag.DontAdjustForAntialiasing, False)
        self.map_view.setOptimizationFlag(QGraphicsView.OptimizationFlag.DontSavePainterState, False)
        self.map_view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
        
        # Enable mouse tracking for better interaction
        self.map_view.setMouseTracking(True)
        self.map_view.viewport().setMouseTracking(True)
        
        # Custom zoom controls
        self.map_view.wheelEvent = self.map_view_wheel_event
        self.map_view.mousePressEvent = self.map_view_mouse_press
        self.map_view.mouseMoveEvent = self.map_view_mouse_move
        self.map_view.mouseReleaseEvent = self.map_view_mouse_release
        
        # Zoom state
        self.zoom_factor = 1.0
        self.min_zoom = 0.1
        self.max_zoom = 10.0
        self.last_mouse_pos = None
        self.is_panning = False
        
        self.sequence_view = SequenceTextView()
        
        center_panel.addWidget(self.map_view)
        center_panel.addWidget(self.sequence_view)

        # Translation view
        translation_group = QGroupBox("Translation")
        translation_group.setStyleSheet("""
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
        translation_layout = QVBoxLayout(translation_group)
        translation_layout.setContentsMargins(10, 15, 10, 10)
        translation_layout.setSpacing(8)
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
        tools_group.setStyleSheet("""
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
        tools_layout = QVBoxLayout(tools_group)
        tools_layout.setContentsMargins(10, 15, 10, 10)
        tools_layout.setSpacing(8)
        
        # Restriction analysis
        restriction_group = QGroupBox("Restriction Analysis")
        restriction_group.setStyleSheet("""
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
        restriction_layout = QFormLayout(restriction_group)
        restriction_layout.setContentsMargins(10, 15, 10, 10)
        restriction_layout.setSpacing(8)
        
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
        cloning_group.setStyleSheet("""
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
        cloning_layout = QVBoxLayout(cloning_group)
        cloning_layout.setContentsMargins(10, 15, 10, 10)
        cloning_layout.setSpacing(8)
        
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
        self.digest_button = QPushButton("⚡ Virtual Digest (Select Enzymes)")
        self.digest_button.setToolTip("Select specific enzymes for virtual restriction digest")
        self.digest_button.clicked.connect(self.virtual_digest)
        cloning_layout.addWidget(self.digest_button)
        
        layout.addWidget(cloning_group)
        
            # Enhanced tools if available
        if ENHANCED_FEATURES_AVAILABLE:
            enhanced_group = QGroupBox("Enhanced Tools")
            enhanced_group.setStyleSheet("""
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
            enhanced_layout = QVBoxLayout(enhanced_group)
            enhanced_layout.setContentsMargins(10, 15, 10, 10)
            enhanced_layout.setSpacing(8)
            
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
            dialog.exec()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening Golden Gate dialog: {str(e)}")

    def open_gibson_dialog(self):
        """Open Gibson assembly dialog"""
        try:
            dialog = GibsonAssemblyDialog(self)
            dialog.exec()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening Gibson assembly dialog: {str(e)}")

    def virtual_digest(self):
        """Perform virtual restriction digest with enzyme selection"""
        if not self.record:
            QMessageBox.warning(self, "No Sequence", "Please load a sequence first.")
            return
        
        # Get current restriction enzymes
        current_enzymes = []
        if self.restriction_batch:
            try:
                current_enzymes = [str(enzyme) for enzyme in self.restriction_batch]
            except:
                current_enzymes = getattr(self, 'current_enzyme_names', [])
        
        # Show enzyme selection dialog
        from src.gui.digest_enzyme_selection_dialog import DigestEnzymeSelectionDialog
        dialog = DigestEnzymeSelectionDialog(current_enzymes, self)
        
        if dialog.exec() == QDialog.DialogCode.Accepted:
            selected_enzymes = dialog.get_selected_enzymes()
            
            if not selected_enzymes:
                QMessageBox.information(self, "No Enzymes", "No enzymes selected for digest.")
                return
            
            try:
                # Create restriction batch with selected enzymes
                digest_batch = RestrictionBatch(selected_enzymes)
                analysis = digest_batch.search(self.record.seq)
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
                
                # Double enzyme digests (combinations of 2)
                if len(cutting_enzymes) >= 2:
                    for i, (enzyme1, sites1) in enumerate(cutting_enzymes):
                        for enzyme2, sites2 in cutting_enzymes[i+1:]:
                            combined_sites = sorted(set(sites1 + sites2))
                            fragments = self.get_fragments(self.record.seq, combined_sites)
                            all_fragments.extend(fragments)
                            enzyme_fragments[f"{enzyme1} + {enzyme2}"] = fragments
                
                # Show virtual gel
                from src.gui.virtual_gel_dialog import VirtualGelDialog
                gel_dialog = VirtualGelDialog(all_fragments, self, enzyme_fragments)
                gel_dialog.exec()
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error during virtual digest: {str(e)}")
    
    def map_view_wheel_event(self, event):
        """Handle mouse wheel for zooming"""
        try:
            # Get zoom delta
            zoom_in_factor = 1.25
            zoom_out_factor = 1 / zoom_in_factor
            
            # Save the scene pos
            old_pos = self.map_view.mapToScene(event.position().toPoint())
            
            # Zoom
            if event.angleDelta().y() > 0:
                zoom_factor = zoom_in_factor
            else:
                zoom_factor = zoom_out_factor
            
            # Calculate new zoom factor
            new_zoom = self.zoom_factor * zoom_factor
            if self.min_zoom <= new_zoom <= self.max_zoom:
                self.zoom_factor = new_zoom
                self.map_view.scale(zoom_factor, zoom_factor)
                
                # Get the new position
                new_pos = self.map_view.mapToScene(event.position().toPoint())
                
                # Move scene to old position
                delta = new_pos - old_pos
                self.map_view.translate(delta.x(), delta.y())
        except Exception as e:
            # Fallback to default behavior
            QGraphicsView.wheelEvent(self.map_view, event)
    
    def map_view_mouse_press(self, event):
        """Handle mouse press for panning and selection"""
        if event.button() == Qt.MouseButton.MiddleButton:
            self.is_panning = True
            self.last_mouse_pos = event.position().toPoint()
            self.map_view.setCursor(Qt.CursorShape.ClosedHandCursor)
        elif event.button() == Qt.MouseButton.LeftButton:
            # Check if clicking on a feature
            scene_pos = self.map_view.mapToScene(event.position().toPoint())
            item = self.scene.itemAt(scene_pos, self.map_view.transform())
            if item and item in self.feature_items.values():
                # Find the feature index
                for idx, feature_item in self.feature_items.items():
                    if feature_item == item:
                        self.highlight_feature(idx)
                        break
        else:
            QGraphicsView.mousePressEvent(self.map_view, event)
    
    def map_view_mouse_move(self, event):
        """Handle mouse move for panning"""
        if self.is_panning and self.last_mouse_pos:
            current_pos = event.position().toPoint()
            delta = current_pos - self.last_mouse_pos
            self.last_mouse_pos = current_pos
            
            # Pan the view
            self.map_view.translate(delta.x() / self.zoom_factor, delta.y() / self.zoom_factor)
        else:
            QGraphicsView.mouseMoveEvent(self.map_view, event)
    
    def map_view_mouse_release(self, event):
        """Handle mouse release"""
        if event.button() == Qt.MouseButton.MiddleButton:
            self.is_panning = False
            self.last_mouse_pos = None
            self.map_view.setCursor(Qt.CursorShape.ArrowCursor)
        else:
            QGraphicsView.mouseReleaseEvent(self.map_view, event)
    
    def on_feature_selection_changed(self):
        """Handle feature selection change from the feature list"""
        selected_items = self.feature_list.selectedItems()
        if selected_items:
            selected_item = selected_items[0]
            feature_index = selected_item.data(Qt.ItemDataRole.UserRole)
            
            # Check if this is a restriction site ( UserRole is None )
            if feature_index is None:
                # This is a restriction site - don't highlight anything on the map
                # but we could potentially highlight the restriction site marker
                pass
            else:
                # This is a regular feature - highlight it
                self.highlight_feature(feature_index)

    def highlight_feature(self, feature_index):
        """Highlight a specific feature"""
        if not self.record or not hasattr(self.record, 'features'):
            return
            
        # Reset all features to normal appearance
        for idx, item in self.feature_items.items():
            if isinstance(item, QGraphicsPathItem):  # Circular view
                pen = item.pen()
                color = self.get_feature_color(self.record.features[idx].type)
                pen.setColor(color)
                pen.setWidth(12)
                item.setPen(pen)
            elif isinstance(item, QGraphicsRectItem):  # Linear view
                color = self.get_feature_color(self.record.features[idx].type)
                item.setBrush(QBrush(color))
                item.setPen(QPen(color))
        
        # Highlight selected feature
        if feature_index < len(self.record.features):
            selected_item = self.feature_items.get(feature_index)
            if selected_item:
                if isinstance(selected_item, QGraphicsPathItem):
                    pen = selected_item.pen()
                    pen.setColor(QColor(255, 255, 0))  # Yellow highlight
                    pen.setWidth(16)
                    selected_item.setPen(pen)
                elif isinstance(selected_item, QGraphicsRectItem):
                    selected_item.setBrush(QBrush(QColor(255, 255, 0, 180)))  # Yellow highlight
                    selected_item.setPen(QPen(QColor(255, 200, 0), 3))
                
                # Update feature list selection
                if hasattr(self, 'feature_list') and feature_index < self.feature_list.count():
                    self.feature_list.setCurrentRow(feature_index)
    
    def zoom_in(self):
        """Zoom in the map view"""
        if self.zoom_factor < self.max_zoom:
            self.map_view.scale(1.25, 1.25)
            self.zoom_factor *= 1.25
    
    def zoom_out(self):
        """Zoom out the map view"""
        if self.zoom_factor > self.min_zoom:
            self.map_view.scale(0.8, 0.8)
            self.zoom_factor *= 0.8
    
    def reset_zoom(self):
        """Reset zoom to fit view"""
        self.zoom_factor = 1.0
        self.map_view.resetTransform()
        if hasattr(self, 'scene') and self.scene.items():
            self.map_view.fitInView(self.scene.itemsBoundingRect(), Qt.AspectRatioMode.KeepAspectRatio)
    
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
            dialog.exec()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening cloning dialog: {str(e)}")

    def select_enzymes(self):
        """Open enzyme selection dialog"""
        try:
            dialog = EnzymeSelectionDialog(self)
            if dialog.exec() == QDialog.DialogCode.Accepted:
                selected_enzymes = dialog.get_selected_enzymes()
                if selected_enzymes:
                    self.restriction_batch = RestrictionBatch(selected_enzymes)
                    self.current_enzyme_names = selected_enzymes  # Store the names
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
            "GenBank Files (*.gb *.gbk);;FASTA Files (*.fa *.fasta);;All Files (*)"
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
            feature_index = 0
            
            # Add regular features
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
                        item = QListWidgetItem(item_text)
                        item.setData(Qt.ItemDataRole.UserRole, feature_index)
                        self.feature_list.addItem(item)
                        feature_index += 1
            
            # Add restriction sites
            if self.restriction_batch:
                try:
                    analysis = self.restriction_batch.search(self.record.seq)
                    for enzyme, sites in analysis.items():
                        if sites:  # Only add enzymes that actually cut
                            enzyme_name = str(enzyme)
                            for site in sites:
                                item_text = f"Restriction Site: {enzyme_name} [{site}]"
                                item = QListWidgetItem(item_text)
                                item.setData(Qt.ItemDataRole.UserRole, None)  # Restriction sites don't correspond to feature items
                                self.feature_list.addItem(item)
                except Exception as e:
                    print(f"Error adding restriction sites to list: {e}")

        self.update_view()

    def toggle_circular_view(self):
        """Switch to circular view"""
        self.circular_action.setChecked(True)
        self.linear_action.setChecked(False)
        self.show_circular_view()
    
    def toggle_linear_view(self):
        """Switch to linear view"""
        self.linear_action.setChecked(True)
        self.circular_action.setChecked(False)
        self.show_linear_view()

    def update_view(self):
        """Update the sequence view based on current settings"""
        if hasattr(self, 'circular_action') and self.circular_action.isChecked():
            self.show_circular_view()
        else:
            self.show_linear_view()
        
        # Fit the view to show all content
        if hasattr(self, 'scene') and hasattr(self, 'map_view'):
            self.map_view.fitInView(self.scene.itemsBoundingRect(), Qt.AspectRatioMode.KeepAspectRatio)

    def show_circular_view(self):
        """Display sequence in circular plasmid view - SnapGene/Benchling style"""
        self.scene.clear()
        
        if not self.record:
            radius = 150
            self.scene.addEllipse(-radius, -radius, radius * 2, radius * 2, QPen(Qt.GlobalColor.gray, 2))
            text = self.scene.addText("Load a sequence\nto view plasmid map", QFont("Arial", 12))
            text.setPos(-60, -10)
            return

        plasmid_len = len(self.record.seq)
        
        # Define radii for different elements
        main_radius = 150
        outer_radius = main_radius + 40  # For restriction sites
        label_radius = main_radius + 60  # For restriction site labels
        center_x, center_y = 0, 0

        # Draw main plasmid backbone
        pen = QPen(Qt.GlobalColor.black, 3)
        self.scene.addEllipse(center_x - main_radius, center_y - main_radius, main_radius * 2, main_radius * 2, pen)

        # Draw outer circle for restriction sites
        outer_pen = QPen(Qt.GlobalColor.gray, 1, Qt.PenStyle.DashLine)
        self.scene.addEllipse(center_x - outer_radius, center_y - outer_radius, outer_radius * 2, outer_radius * 2, outer_pen)

        # Draw features with SnapGene-style layering system
        self.feature_items.clear()
        all_labels = []  # Collect all labels for smart positioning
        
        if hasattr(self.record, 'features'):
            feature_index = 0
            
            # Collect and sort features by position
            feature_data = []
            for i, feature in enumerate(self.record.features):
                if feature.type in ["gene", "promoter", "CDS", "terminator", "rep_origin", "misc_feature", 
                           "regulatory", "enhancer", "primer_bind", "protein_bind", "RBS", 
                           "5'UTR", "3'UTR", "signal_peptide", "stem_loop", "polyA_signal"]:
                    
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    
                    # Handle wraparound features
                    if start > end:
                        end += plasmid_len
                    
                    start_angle = (start / plasmid_len) * 360
                    span_angle = ((end - start) / plasmid_len) * 360
                    
                    # Only include features large enough to be visible
                    if span_angle < 2:
                        continue
                    
                    feature_data.append({
                        'start_angle': start_angle,
                        'span_angle': span_angle,
                        'feature': feature,
                        'feature_type': feature.type,
                        'start': start,
                        'end': end,
                        'index': i
                    })
            
            # Sort features by start angle
            feature_data.sort(key=lambda x: x['start_angle'])
            
            # Assign features to layers (tracks) to avoid overlaps
            layers = []  # Each layer will have features that don't overlap
            layer_spacing = 12  # Pixels between layers
            
            for data in feature_data:
                assigned_layer = None
                
                # Try to fit in existing layers
                for layer_idx, layer in enumerate(layers):
                    can_fit = True
                    
                    for existing_feature in layer:
                        # Check angular overlap
                        angle_diff = abs(data['start_angle'] - existing_feature['start_angle'])
                        if angle_diff > 180:
                            angle_diff = 360 - angle_diff
                        
                        # Features overlap if they're too close
                        min_separation = (data['span_angle'] + existing_feature['span_angle']) / 2 + 5
                        if angle_diff < min_separation:
                            can_fit = False
                            break
                    
                    if can_fit:
                        assigned_layer = layer_idx
                        layer.append(data)
                        break
                
                # If no existing layer fits, create new layer
                if assigned_layer is None:
                    layers.append([data])
                    assigned_layer = len(layers) - 1
            
            # Draw features in their assigned layers
            for layer_idx, layer in enumerate(layers):
                # Calculate radius for this layer
                layer_radius = main_radius - 10 - (layer_idx * layer_spacing)
                
                # Only draw if radius is reasonable
                if layer_radius < 30:
                    break
                
                for data in layer:
                    start_angle = data['start_angle']
                    span_angle = data['span_angle']
                    feature = data['feature']
                    pen_color = self.get_feature_color(feature.type)
                    
                    # Draw feature as arc on this layer
                    feature_width = max(6, min(15, span_angle / 4))
                    pen = QPen(pen_color, feature_width)
                    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                    
                    path = QPainterPath()
                    path.arcMoveTo(center_x - layer_radius, center_y - layer_radius, layer_radius * 2, layer_radius * 2, -start_angle)
                    path.arcTo(center_x - layer_radius, center_y - layer_radius, layer_radius * 2, layer_radius * 2, -start_angle, -span_angle)
                    arc = self.scene.addPath(path, pen)
                    
                    # Make feature clickable
                    arc.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, True)
                    self.feature_items[feature_index] = arc
                    
                    # Add arrow for direction
                    self.add_circular_arrow(feature, center_x, center_y, layer_radius, start_angle, span_angle, pen_color)
                    
                    # Collect labels for smart positioning
                    if span_angle > 15:  # Only label larger features
                        label_angle = start_angle + (span_angle / 2)
                        
                        # Get feature name
                        feature_name = feature.qualifiers.get('label', [feature.type])[0]
                        if 'gene' in feature.qualifiers:
                            feature_name = feature.qualifiers['gene'][0]
                        
                        # Truncate long names
                        if len(feature_name) > 12:
                            feature_name = feature_name[:10] + "..."
                        
                        # Create label item
                        label = self.scene.addText(feature_name, QFont("Arial", 7, QFont.Weight.Bold))
                        label.setDefaultTextColor(Qt.GlobalColor.black)
                        
                        # Store label data for smart positioning
                        all_labels.append({
                            'label': label,
                            'angle': label_angle,
                            'feature_name': feature_name,
                            'layer_idx': layer_idx,
                            'feature_angle': start_angle,
                            'feature_span': span_angle
                        })
                    
                    feature_index += 1
        
        # Smart label positioning system
        if all_labels:
            # Sort labels by angle
            all_labels.sort(key=lambda x: x['angle'])
            
            # Define label radii for different layers
            label_radii = [main_radius + 35, main_radius + 55, main_radius + 75, main_radius + 95]
            
            # Assign labels to concentric label layers to avoid overlaps
            label_layers = [[] for _ in label_radii]
            
            for label_data in all_labels:
                angle = label_data['angle']
                assigned_label_layer = 0
                
                # Try to fit in existing label layers
                for layer_idx, layer in enumerate(label_layers):
                    can_fit = True
                    
                    for existing_label in layer:
                        angle_diff = abs(angle - existing_label['angle'])
                        if angle_diff > 180:
                            angle_diff = 360 - angle_diff
                        
                        # Labels overlap if they're too close (consider label width)
                        if angle_diff < 20:  # Minimum 20 degrees separation
                            can_fit = False
                            break
                    
                    if can_fit:
                        assigned_label_layer = layer_idx
                        layer.append(label_data)
                        break
                
                # If no existing layer fits, try next available layer
                if assigned_label_layer >= len(label_layers):
                    continue  # Skip this label if no space available
                
                # Position the label
                radius = label_radii[assigned_label_layer]
                rad_angle = math.radians(angle)
                label_x = center_x + radius * math.cos(rad_angle)
                label_y = center_y + radius * math.sin(rad_angle)
                
                label = label_data['label']
                
                # Adjust position for readability
                if 90 <= angle <= 270:
                    # Left side - align right
                    label.setPos(label_x - 25, label_y - 5)
                else:
                    # Right side - align left
                    label.setPos(label_x + 5, label_y - 5)
                
                # Rotate label to be readable
                if 90 <= angle <= 270:
                    label.setRotation(angle + 180)
                else:
                    label.setRotation(angle)

        # Draw restriction sites on outer circle with labels
        if self.restriction_batch:
            try:
                analysis = self.restriction_batch.search(self.record.seq)
                site_positions = []
                
                # Collect all restriction sites
                for enzyme, sites in analysis.items():
                    for site in sites:
                        site_positions.append((site, str(enzyme)))  # Convert to string
                
                # Sort by position
                site_positions.sort(key=lambda x: x[0])
                
                # Draw restriction sites with labels
                for site, enzyme_name in site_positions:  # Use enzyme_name
                    angle = (site / plasmid_len) * 360
                    rad_angle = math.radians(angle)
                    
                    # Draw tick mark on outer circle
                    tick_x = center_x + outer_radius * math.cos(rad_angle)
                    tick_y = center_y + outer_radius * math.sin(rad_angle)  # Fixed: was cos instead of sin
                    
                    inner_x = center_x + (main_radius + 5) * math.cos(rad_angle)
                    inner_y = center_y + (main_radius + 5) * math.sin(rad_angle)
                    
                    self.scene.addLine(inner_x, inner_y, tick_x, tick_y, QPen(Qt.GlobalColor.darkRed, 2))
                    
                    # Add enzyme label
                    label_x = center_x + label_radius * math.cos(rad_angle)
                    label_y = center_y + label_radius * math.sin(rad_angle)
                    
                    # Position label to be readable
                    label = self.scene.addText(enzyme_name, QFont("Arial", 7, QFont.Weight.Bold))  # Use enzyme_name
                    label.setDefaultTextColor(Qt.GlobalColor.darkRed)
                    
                    # Adjust label position for readability
                    if 90 <= angle <= 270:
                        # Left side - align right
                        label.setPos(label_x - 25, label_y - 5)
                    else:
                        # Right side - align left
                        label.setPos(label_x + 5, label_y - 5)
                        
            except Exception as e:
                print(f"Error drawing restriction sites: {e}")

        # Add center info
        info_text = f"{self.record.id}\n{plasmid_len:,} bp"
        center_label = self.scene.addText(info_text, QFont("Arial", 10, QFont.Weight.Bold))
        center_label.setDefaultTextColor(Qt.GlobalColor.black)
        center_label.setPos(-30, -15)

    def add_scale_bar(self, sequence_length, radius):
        """Add a scale bar to the plasmid map"""
        # Calculate appropriate scale
        if sequence_length < 1000:
            scale_size = 100
            scale_text = f"{scale_size} bp"
        elif sequence_length < 10000:
            scale_size = 1000
            scale_text = f"{scale_size/1000:.0f} kb"
        else:
            scale_size = 10000
            scale_text = f"{scale_size/1000:.0f} kb"
        
        # Draw scale bar at bottom
        scale_length = (scale_size / sequence_length) * 2 * math.pi * radius
        scale_y = radius + 50
        
        # Scale bar line
        scale_start_x = -scale_length / 2
        scale_end_x = scale_length / 2
        
        self.scene.addLine(scale_start_x, scale_y, scale_end_x, scale_y, QPen(Qt.GlobalColor.black, 2))
        
        # Scale bar ticks
        self.scene.addLine(scale_start_x, scale_y - 5, scale_start_x, scale_y + 5, QPen(Qt.GlobalColor.black, 2))
        self.scene.addLine(scale_end_x, scale_y - 5, scale_end_x, scale_y + 5, QPen(Qt.GlobalColor.black, 2))
        
        # Scale label
        scale_label = self.scene.addText(scale_text, QFont("Arial", 8))
        scale_label.setDefaultTextColor(Qt.GlobalColor.black)
        scale_label.setPos(-20, scale_y + 8)

    def show_linear_view(self):
        """Display sequence in linear view - professional style with reduced overlaps"""
        self.scene.clear()

        if not self.record:
            self.scene.addLine(0, 0, 800, 0, QPen(Qt.GlobalColor.gray, 2))
            text = self.scene.addText("Load a sequence to view linear map", QFont("Arial", 12))
            text.setPos(200, -30)
            return

        plasmid_len = len(self.record.seq)
        
        # Calculate scale for display
        max_width = 1400
        scale_factor = max(0.05, min(1.0, max_width / plasmid_len))
        scaled_len = plasmid_len * scale_factor
        
        # Center the map
        start_x = -scaled_len / 2
        center_y = 0
        
        # Define track positions to avoid overlaps
        tracks = {
            'gene': center_y - 40,
            'CDS': center_y - 40, 
            'promoter': center_y - 60,
            'enhancer': center_y - 60,
            'regulatory': center_y - 60,
            'terminator': center_y + 25,
            'rep_origin': center_y + 25,
            'misc_feature': center_y + 45,
            'primer_bind': center_y + 65,
            'protein_bind': center_y + 65,
            'RBS': center_y + 45,
            '5\'UTR': center_y + 45,
            '3\'UTR': center_y + 45,
            'signal_peptide': center_y + 65,
            'stem_loop': center_y + 65,
            'polyA_signal': center_y + 65
        }
        
        # Draw main sequence line - thick black line
        self.scene.addLine(start_x, center_y, start_x + scaled_len, center_y, QPen(Qt.GlobalColor.black, 4))

        # Draw features with professional layout
        self.feature_items.clear()
        if hasattr(self.record, 'features'):
            feature_index = 0
            
            for i, feature in enumerate(self.record.features):
                if feature.type in tracks:
                    
                    start = int(feature.location.start) * scale_factor + start_x
                    end = int(feature.location.end) * scale_factor + start_x
                    
                    # Ensure minimum width for visibility
                    if end - start < 3:
                        end = start + 3
                    
                    pen_color = self.get_feature_color(feature.type)
                    y_pos = tracks[feature.type]
                    
                    # Draw feature as rectangle
                    feature_height = 6
                    rect = self.scene.addRect(start, y_pos - feature_height/2, end - start, feature_height, 
                                           QPen(pen_color, 1), QBrush(pen_color))
                    
                    # Make feature clickable
                    rect.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, True)
                    self.feature_items[feature_index] = rect
                    
                    # Add thin connecting line to main sequence
                    line_start_x = start + (end - start) / 2
                    self.scene.addLine(line_start_x, center_y, line_start_x, y_pos, 
                                     QPen(pen_color, 1, Qt.PenStyle.DotLine))
                    
                    # Add feature label only for larger features
                    feature_width = end - start
                    if feature_width > 25:
                        # Get feature name
                        feature_name = feature.qualifiers.get('label', [feature.type])[0]
                        if 'gene' in feature.qualifiers:
                            feature_name = feature.qualifiers['gene'][0]
                        
                        # Truncate long names
                        if len(feature_name) > 12:
                            feature_name = feature_name[:10] + "..."
                        
                        # Add label
                        label = self.scene.addText(feature_name, QFont("Arial", 6, QFont.Weight.Bold))
                        label.setDefaultTextColor(Qt.GlobalColor.black)
                        
                        # Position label above or below based on track position
                        label_width = label.boundingRect().width()
                        label_x = start + (feature_width - label_width) / 2
                        
                        if y_pos < center_y:
                            label_y = y_pos - 15  # Above for upper tracks
                        else:
                            label_y = y_pos + 10  # Below for lower tracks
                        
                        # Only show label if it fits
                        if label_width <= feature_width * 2:
                            label.setPos(label_x, label_y)
                    
                    feature_index += 1

        # Draw scale markers with better spacing
        scale_intervals = [100, 500, 1000, 5000, 10000]
        scale_interval = 100
        for interval in scale_intervals:
            if plasmid_len / interval <= 20:  # Not too many markers
                scale_interval = interval
                break
        
        # Draw scale ticks and labels
        for pos in range(0, plasmid_len + 1, scale_interval):
            scaled_pos = pos * scale_factor + start_x
            
            # Draw tick mark
            self.scene.addLine(scaled_pos, center_y - 10, scaled_pos, center_y + 10, QPen(Qt.GlobalColor.black, 2))
            
            # Add label
            if plasmid_len <= 1000:
                label_text = f"{pos}"
            elif plasmid_len <= 10000:
                label_text = f"{pos//1000}.{(pos%1000)//100}kb"
            else:
                label_text = f"{pos//1000}kb"
            
            label = self.scene.addText(label_text, QFont("Arial", 7, QFont.Weight.Bold))
            label.setDefaultTextColor(Qt.GlobalColor.black)
            label.setPos(scaled_pos - 15, center_y + 20)

        # Draw restriction sites on separate track below
        if self.restriction_batch:
            try:
                analysis = self.restriction_batch.search(self.record.seq)
                restriction_y = center_y + 85
                site_count = 0
                
                # Draw restriction site line
                self.scene.addLine(start_x, restriction_y, start_x + scaled_len, restriction_y, 
                                 QPen(Qt.GlobalColor.darkRed, 1, Qt.PenStyle.DashLine))
                
                for enzyme, sites in analysis.items():
                    for site in sites:
                        if site_count >= 50:  # Limit to prevent overcrowding
                            break
                        
                        scaled_site = site * scale_factor + start_x
                        
                        # Draw restriction site as tick mark
                        self.scene.addLine(scaled_site, restriction_y - 5, scaled_site, restriction_y + 5, 
                                         QPen(Qt.GlobalColor.darkRed, 2))
                        
                        # Add enzyme label for major sites
                        if len(sites) <= 10:
                            label = self.scene.addText(str(enzyme), QFont("Arial", 6))  # Convert to string
                            label.setDefaultTextColor(Qt.GlobalColor.darkRed)
                            label.setPos(scaled_site - 10, restriction_y + 8)
                        
                        site_count += 1
            except Exception as e:
                print(f"Error drawing restriction sites: {e}")

        # Add sequence info at bottom
        info_text = f"{self.record.id} - {plasmid_len:,} bp"
        info_label = self.scene.addText(info_text, QFont("Arial", 11, QFont.Weight.Bold))
        info_label.setDefaultTextColor(Qt.GlobalColor.black)
        info_label.setPos(start_x, center_y + 110)

    def get_fragments(self, sequence, cuts):
        """Get DNA fragments from restriction cuts"""
        if not cuts:
            return [sequence]
        
        fragments = []
        sorted_cuts = sorted(cuts)
        
        # For circular sequences (plasmids), create fragments between consecutive cuts
        if len(sorted_cuts) == 1:
            # Single cut - linearize the plasmid
            cut_pos = sorted_cuts[0]
            fragments.append(sequence[cut_pos:] + sequence[:cut_pos])
        else:
            # Multiple cuts - create fragments between cuts
            for i in range(len(sorted_cuts)):
                start = sorted_cuts[i]
                end = sorted_cuts[(i + 1) % len(sorted_cuts)]
                
                if i < len(sorted_cuts) - 1:
                    # Normal fragment
                    fragment = sequence[start:end]
                else:
                    # Wraparound fragment
                    fragment = sequence[start:] + sequence[:end]
                
                if fragment:  # Only add non-empty fragments
                    fragments.append(fragment)
        
        return fragments

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
            "protein_bind": QColor("#34495e"), # Dark gray
            "RBS": QColor("#16a085"),          # Dark turquoise
            "5'UTR": QColor("#27ae60"),        # Dark green
            "3'UTR": QColor("#2980b9"),        # Dark blue
            "signal_peptide": QColor("#8e44ad"), # Dark purple
            "stem_loop": QColor("#d35400"),     # Dark orange
            "polyA_signal": QColor("#c0392b")   # Dark red
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
                                             QPen(Qt.GlobalColor.darkRed, 1), 
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
                                             QPen(Qt.GlobalColor.darkRed, 1), 
                                             QBrush(QColor(200, 50, 50, 120)))  # More transparent
            
        except Exception as e:
            print(f"Error adding linear restriction arrow: {e}")
    
    def add_linear_direction_indicator(self, sequence_length):
        """Add elegant 5' to 3' direction indicator for linear view"""
        try:
            # Add 5' and 3' labels with better styling
            label_5prime = self.scene.addText("5'", QFont("Arial", 10, QFont.Weight.Bold))
            label_5prime.setDefaultTextColor(QColor(70, 130, 180))  # Steel blue
            label_5prime.setPos(-25, -8)
            
            label_3prime = self.scene.addText("3'", QFont("Arial", 10, QFont.Weight.Bold))
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
            dialog.exec()
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
            dialog.exec()
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
            dialog.exec()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening advanced analysis: {str(e)}")
    
    def on_sequence_modified(self, modified_record):
        """Handle sequence modifications from the editor"""
        if modified_record:
            self.display_sequence(modified_record)
            if hasattr(self, 'parent_app') and self.parent_app:
                self.parent_app.statusBar().showMessage("Sequence modified and updated", 3000)
