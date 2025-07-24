#!/usr/bin/env python3
"""
Structural Analysis Tools for PicoMol.

This module provides protein structural analysis including:
- Secondary structure analysis
- Geometric property calculations
- Basic structural properties
"""

import os
import math
import numpy as np
from collections import defaultdict, Counter
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTabWidget, QGroupBox, QFormLayout,
    QTextEdit, QLineEdit, QPushButton, QLabel, QComboBox, QSpinBox,
    QCheckBox, QTableWidget, QTableWidgetItem, QHeaderView, QScrollArea,
    QMessageBox, QFileDialog, QProgressBar, QSplitter, QDialog, QDialogButtonBox,
    QListWidget, QListWidgetItem, QSplitter as QSplitterWidget, QFrame,
    QGridLayout, QTextBrowser, QApplication, QDoubleSpinBox
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont, QPixmap, QPainter, QPen, QBrush, QColor, QPainterPath

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    from Bio.PDB import PDBParser, NeighborSearch, PDBIO, PDBList
    from Bio.PDB.vectors import Vector, calc_angle, calc_dihedral
    from Bio.PDB.Polypeptide import PPBuilder, is_aa
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    import ramachandraw
    RAMACHANDRAW_AVAILABLE = True
except ImportError:
    RAMACHANDRAW_AVAILABLE = False

try:
    from .pdb_fetch_worker import PDBFetchManager
    OPTIMIZED_FETCH_AVAILABLE = True
except ImportError:
    try:
        # Add project root to path for absolute imports
        import sys
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        if project_root not in sys.path:
            sys.path.insert(0, project_root)
        from src.core.pdb_fetch_worker import PDBFetchManager
        OPTIMIZED_FETCH_AVAILABLE = True
    except ImportError:
        OPTIMIZED_FETCH_AVAILABLE = False

try:
    from .enhanced_pdb_puller_fixed import EnhancedPDBPuller
    ENHANCED_PDB_AVAILABLE = True
except ImportError:
    try:
        # Add project root to path for absolute imports
        import sys
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        if project_root not in sys.path:
            sys.path.insert(0, project_root)
        from src.core.enhanced_pdb_puller_fixed import EnhancedPDBPuller
        ENHANCED_PDB_AVAILABLE = True
    except ImportError:
        try:
            from .enhanced_pdb_puller import EnhancedPDBPuller
            ENHANCED_PDB_AVAILABLE = True
        except ImportError:
            try:
                from src.core.enhanced_pdb_puller import EnhancedPDBPuller
                ENHANCED_PDB_AVAILABLE = True
            except ImportError:
                ENHANCED_PDB_AVAILABLE = False

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    # Fallback implementations
    class np:
        @staticmethod
        def array(data):
            return data
        @staticmethod
        def mean(data):
            return sum(data) / len(data) if data else 0
        @staticmethod
        def std(data):
            if not data:
                return 0
            mean_val = sum(data) / len(data)
            variance = sum((x - mean_val) ** 2 for x in data) / len(data)
            return math.sqrt(variance)




class StructuralAnalysisWorker(QThread):
    """Worker thread for structural analysis calculations."""
    
    analysis_complete = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)
    progress_update = pyqtSignal(str)
    
    def __init__(self, structure_path, analysis_types=None, comprehensive_data=None):
        super().__init__()
        self.structure_path = structure_path
        self.analysis_types = analysis_types or ['basic', 'secondary', 'geometry']
        self.structure = None
        self.comprehensive_data = comprehensive_data
    
    def run(self):
        try:
            if not BIOPYTHON_AVAILABLE:
                self.error_occurred.emit("Biopython is required for structural analysis. Please install it with: pip install biopython")
                return
            
            self.progress_update.emit("Loading structure...")
            
            # Load structure
            parser = PDBParser(QUIET=True)
            structure_id = os.path.basename(self.structure_path).split('.')[0]
            self.structure = parser.get_structure(structure_id, self.structure_path)
            
            results = {
                'structure_id': structure_id,
                'structure_path': self.structure_path
            }
            
            # Perform requested analyses
            if 'basic' in self.analysis_types:
                self.progress_update.emit("Analyzing basic properties...")
                results['basic'] = self.analyze_basic_properties()
            
            if 'secondary' in self.analysis_types:
                self.progress_update.emit("Analyzing secondary structure...")
                results['secondary'] = self.analyze_secondary_structure()
            
            if 'geometry' in self.analysis_types:
                self.progress_update.emit("Analyzing geometric properties...")
                results['geometry'] = self.analyze_geometry()
            
            self.progress_update.emit("Analysis complete!")
            self.analysis_complete.emit(results)
            
        except Exception as e:
            self.error_occurred.emit(f"Structural analysis error: {str(e)}")
    
    def analyze_basic_properties(self):
        """Analyze basic structural properties."""
        results = {
            'chains': [],
            'total_residues': 0,
            'total_atoms': 0,
            'molecular_weight': 0.0,
            'composition': {}
        }
        
        # Amino acid molecular weights (average)
        aa_weights = {
            'ALA': 89.1, 'ARG': 174.2, 'ASN': 132.1, 'ASP': 133.1, 'CYS': 121.0,
            'GLN': 146.2, 'GLU': 147.1, 'GLY': 75.1, 'HIS': 155.2, 'ILE': 131.2,
            'LEU': 131.2, 'LYS': 146.2, 'MET': 149.2, 'PHE': 165.2, 'PRO': 115.1,
            'SER': 105.1, 'THR': 119.1, 'TRP': 204.2, 'TYR': 181.2, 'VAL': 117.1
        }
        
        composition = Counter()
        
        for model in self.structure:
            for chain in model:
                chain_info = {
                    'id': chain.id,
                    'residues': 0,
                    'atoms': 0,
                    'sequence': ''
                }
                
                for residue in chain:
                    if is_aa(residue):
                        chain_info['residues'] += 1
                        results['total_residues'] += 1
                        
                        # Add to composition
                        res_name = residue.get_resname()
                        composition[res_name] += 1
                        
                        # Add to molecular weight
                        results['molecular_weight'] += aa_weights.get(res_name, 110.0)
                        
                        # Build sequence
                        try:
                            from Bio.PDB.Polypeptide import protein_letters_3to1
                            chain_info['sequence'] += protein_letters_3to1.get(res_name, 'X')
                        except:
                            chain_info['sequence'] += 'X'
                    
                    # Count atoms
                    for atom in residue:
                        chain_info['atoms'] += 1
                        results['total_atoms'] += 1
                
                results['chains'].append(chain_info)
        
        results['composition'] = dict(composition)
        
        # Use more accurate molecular weight from comprehensive metadata if available
        if self.comprehensive_data and self.comprehensive_data.get('metadata'):
            metadata = self.comprehensive_data['metadata']
            if 'entry' in metadata:
                entry = metadata['entry']
                rcsb_info = entry.get('rcsb_entry_info', {})
                comp_mol_weight = rcsb_info.get('molecular_weight')
                if comp_mol_weight:
                    results['molecular_weight'] = comp_mol_weight * 1000  # Convert kDa to Da for consistency
        
        return results
    
    def analyze_secondary_structure(self):
        """Analyze secondary structure using simple geometric rules and collect Ramachandran data."""
        results = {
            'helix_content': 0.0,
            'sheet_content': 0.0,
            'loop_content': 0.0,
            'secondary_elements': [],
            'ramachandran_data': [],
            'pdb_id': None  # Add PDB ID to results
        }
        
        total_residues = 0
        helix_residues = 0
        sheet_residues = 0
        ramachandran_data = []
        
        for model in self.structure:
            for chain in model:
                residues = [res for res in chain if is_aa(res)]
                total_residues += len(residues)
                
                # Calculate phi/psi angles and classify secondary structure
                for i in range(1, len(residues) - 1):
                    try:
                        phi, psi = self.calculate_phi_psi(residues, i)
                        
                        # Collect Ramachandran data
                        residue = residues[i]
                        res_name = residue.get_resname()
                        res_id = residue.id[1]
                        chain_id = chain.id
                        
                        # Classify Ramachandran region
                        rama_region = self.classify_ramachandran_region(phi, psi, res_name)
                        
                        ramachandran_data.append({
                            'residue': res_name,
                            'chain': chain_id,
                            'position': res_id,
                            'phi': phi,
                            'psi': psi,
                            'region': rama_region
                        })
                        
                        # Classify secondary structure based on phi/psi
                        ss_type = self.classify_secondary_structure(phi, psi)
                        
                        if ss_type == 'H':  # Helix
                            helix_residues += 1
                        elif ss_type == 'E':  # Sheet
                            sheet_residues += 1
                        
                    except Exception as e:
                        continue
        
        if total_residues > 0:
            results['helix_content'] = (helix_residues / total_residues) * 100
            results['sheet_content'] = (sheet_residues / total_residues) * 100
            results['loop_content'] = 100 - results['helix_content'] - results['sheet_content']
        
        results['ramachandran_data'] = ramachandran_data
        
        # Add PDB ID from comprehensive data if available
        if self.comprehensive_data and self.comprehensive_data.get('pdb_id'):
            results['pdb_id'] = self.comprehensive_data['pdb_id']
        
        return results
    
    def calculate_phi_psi(self, residues, i):
        """Calculate phi and psi dihedral angles for residue i."""
        try:
            # Get atoms for dihedral calculation
            # Phi: C(i-1) - N(i) - CA(i) - C(i)
            # Psi: N(i) - CA(i) - C(i) - N(i+1)
            
            prev_c = residues[i-1]['C'].get_vector()
            curr_n = residues[i]['N'].get_vector()
            curr_ca = residues[i]['CA'].get_vector()
            curr_c = residues[i]['C'].get_vector()
            next_n = residues[i+1]['N'].get_vector()
            
            phi = calc_dihedral(prev_c, curr_n, curr_ca, curr_c)
            psi = calc_dihedral(curr_n, curr_ca, curr_c, next_n)
            
            # Convert to degrees
            phi = math.degrees(phi)
            psi = math.degrees(psi)
            
            return phi, psi
            
        except Exception as e:
            raise Exception(f"Error calculating phi/psi angles: {e}")
    
    def classify_secondary_structure(self, phi, psi):
        """Classify secondary structure based on phi/psi angles."""
        # Simple classification rules
        if -180 <= phi <= -30 and -70 <= psi <= 50:
            return 'H'  # Alpha helix
        elif -180 <= phi <= -30 and 90 <= psi <= 180:
            return 'E'  # Beta sheet
        elif -90 <= phi <= 30 and -180 <= psi <= -90:
            return 'E'  # Beta sheet
        else:
            return 'C'  # Coil/loop
    
    def classify_ramachandran_region(self, phi, psi, residue_name=None):
        """Classify phi/psi angles into Ramachandran regions (simplified version)."""
        # Simplified classification for secondary structure analysis
        if -180 <= phi <= -30 and -70 <= psi <= 50:
            return 'favored'  # Alpha helix region
        elif -180 <= phi <= -30 and 90 <= psi <= 180:
            return 'favored'  # Beta sheet region
        elif -90 <= phi <= 30 and -180 <= psi <= -90:
            return 'favored'  # Beta sheet region
        else:
            return 'outlier'  # Everything else
    
    def analyze_geometry(self):
        """Analyze geometric properties of the structure."""
        results = {
            'bond_lengths': [],
            'bond_angles': [],
            'b_factors': [],
            'center_of_mass': [0, 0, 0],
            'radius_of_gyration': 0.0,
            'geometric_center': [0, 0, 0]
        }
        
        all_atoms = []
        ca_atoms = []
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        for atom in residue:
                            all_atoms.append(atom)
                            results['b_factors'].append(atom.get_bfactor())
                            
                            # Collect CA atoms for radius of gyration
                            if atom.get_name() == 'CA':
                                ca_atoms.append(atom)
        
        # Calculate center of mass and geometric center
        if all_atoms:
            coords = np.array([atom.get_coord() for atom in all_atoms])
            results['geometric_center'] = np.mean(coords, axis=0).tolist()
            
            # Simple center of mass (assuming equal masses)
            results['center_of_mass'] = results['geometric_center']
        
        # Calculate radius of gyration using CA atoms
        if ca_atoms:
            ca_coords = np.array([atom.get_coord() for atom in ca_atoms])
            center = np.mean(ca_coords, axis=0)
            distances_sq = np.sum((ca_coords - center) ** 2, axis=1)
            results['radius_of_gyration'] = math.sqrt(np.mean(distances_sq))
        
        # Sample bond lengths and angles
        self.calculate_sample_geometry(results)
        
        return results
    
    def calculate_sample_geometry(self, results):
        """Calculate sample bond lengths and angles."""
        bond_lengths = []
        bond_angles = []
        
        for model in self.structure:
            for chain in model:
                residues = [res for res in chain if is_aa(res)]
                
                for residue in residues:
                    try:
                        # Sample backbone bond lengths
                        if 'N' in residue and 'CA' in residue:
                            n_ca = residue['N'] - residue['CA']
                            bond_lengths.append(n_ca)
                        
                        if 'CA' in residue and 'C' in residue:
                            ca_c = residue['CA'] - residue['C']
                            bond_lengths.append(ca_c)
                        
                        # Sample backbone bond angles
                        if 'N' in residue and 'CA' in residue and 'C' in residue:
                            n_vec = residue['N'].get_vector()
                            ca_vec = residue['CA'].get_vector()
                            c_vec = residue['C'].get_vector()
                            
                            angle = calc_angle(n_vec, ca_vec, c_vec)
                            bond_angles.append(math.degrees(angle))
                            
                    except Exception:
                        continue
        
        results['bond_lengths'] = bond_lengths[:100]  # Sample first 100
        results['bond_angles'] = bond_angles[:100]   # Sample first 100
    

    









    


    








class StructuralAnalysisTab(QWidget):
    """Tab for structural analysis tools."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_app = parent
        self.analysis_worker = None
        self.current_structure_path = None
        self.pdb_list = PDBList() if BIOPYTHON_AVAILABLE else None
        self.pdb_parser = PDBParser(QUIET=True) if BIOPYTHON_AVAILABLE else None
        self.current_results = None  # Store analysis results for export
        self.current_comprehensive_data = None  # Store comprehensive data
        
        # Create directories if not exist
        if hasattr(parent, 'pulled_structures_dir'):
            self.pulled_structures_dir = parent.pulled_structures_dir
        else:
            # Fallback directory
            project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            self.pulled_structures_dir = os.path.join(project_root, "data", "pulled_structures")
            os.makedirs(self.pulled_structures_dir, exist_ok=True)
        
        # Initialize PDB fetch manager (optimized or fallback)
        if OPTIMIZED_FETCH_AVAILABLE:
            self.pdb_fetch_manager = PDBFetchManager(self.pulled_structures_dir)
            self.enhanced_pdb_puller = None  # For compatibility
        elif ENHANCED_PDB_AVAILABLE:
            self.enhanced_pdb_puller = EnhancedPDBPuller(self.pulled_structures_dir)
            self.pdb_fetch_manager = None
        else:
            self.enhanced_pdb_puller = None
            self.pdb_fetch_manager = None
        
        self.init_ui()
    
    def update_comprehensive_data(self):
        """Update comprehensive data from parent app if available."""
        if hasattr(self.parent_app, 'current_pdb_data') and self.parent_app.current_pdb_data:
            self.current_comprehensive_data = self.parent_app.current_pdb_data
            print(f"[DEBUG] Updated comprehensive data in structural analysis tab: {self.current_comprehensive_data.get('pdb_id', 'Unknown')}")
        else:
            self.current_comprehensive_data = None
            print(f"[DEBUG] No comprehensive data available from parent app")
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)
        
        # Control section
        control_group = QGroupBox("Structure Analysis")
        control_layout = QVBoxLayout(control_group)
        
        # PDB ID Input Section
        pdb_id_layout = QHBoxLayout()
        pdb_id_layout.addWidget(QLabel("PDB ID:"))
        self.pdb_id_entry = QLineEdit()
        self.pdb_id_entry.setPlaceholderText("e.g., 1CRN, 4HHB")
        self.pdb_id_entry.setToolTip("Enter a valid PDB ID to fetch a protein structure from the PDB database.")
        pdb_id_layout.addWidget(self.pdb_id_entry)
        
        fetch_button = QPushButton("Fetch PDB")
        fetch_button.setToolTip("Download and load a protein structure using the entered PDB ID.")
        fetch_button.clicked.connect(self.fetch_pdb_id)
        pdb_id_layout.addWidget(fetch_button)
        
        control_layout.addLayout(pdb_id_layout)
        
        # Structure selection
        structure_layout = QHBoxLayout()
        structure_layout.addWidget(QLabel("Current Structure:"))
        
        self.structure_label = QLabel("No structure loaded")
        self.structure_label.setStyleSheet("color: #666; font-style: italic;")
        structure_layout.addWidget(self.structure_label)
        
        load_current_btn = QPushButton("Load from 3D Viewer")
        load_current_btn.setToolTip("Load the currently displayed structure from the 3D viewer for analysis")
        load_current_btn.clicked.connect(self.load_current_structure)
        structure_layout.addWidget(load_current_btn)
        
        load_file_btn = QPushButton("Load PDB File...")
        load_file_btn.setToolTip("Load a PDB file for structural analysis")
        load_file_btn.clicked.connect(self.load_structure_file)
        structure_layout.addWidget(load_file_btn)
        
        structure_layout.addStretch()
        control_layout.addLayout(structure_layout)
        
        # Analysis options
        options_layout = QHBoxLayout()
        options_layout.addWidget(QLabel("Analysis Types:"))
        
        self.basic_checkbox = QCheckBox("Basic Properties")
        self.basic_checkbox.setChecked(True)
        self.basic_checkbox.setToolTip("Analyze basic structural properties")
        options_layout.addWidget(self.basic_checkbox)
        
        self.secondary_checkbox = QCheckBox("Secondary Structure")
        self.secondary_checkbox.setChecked(True)
        self.secondary_checkbox.setToolTip("Analyze secondary structure content")
        options_layout.addWidget(self.secondary_checkbox)
        
        self.geometry_checkbox = QCheckBox("Geometry")
        self.geometry_checkbox.setChecked(True)
        self.geometry_checkbox.setToolTip("Analyze geometric properties")
        options_layout.addWidget(self.geometry_checkbox)
        

        

        
        options_layout.addStretch()
        control_layout.addLayout(options_layout)
        
        # Analysis and Export buttons
        button_layout = QHBoxLayout()
        
        analyze_btn = QPushButton("Analyze Structure")
        analyze_btn.setToolTip("Perform structural analysis with selected options")
        analyze_btn.clicked.connect(self.analyze_structure)
        button_layout.addWidget(analyze_btn)
        
        self.export_btn = QPushButton("Export Results")
        self.export_btn.setToolTip("Export analysis results to various formats")
        self.export_btn.setEnabled(False)  # Initially disabled
        self.export_btn.clicked.connect(self.export_results)
        button_layout.addWidget(self.export_btn)
        
        control_layout.addLayout(button_layout)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        control_layout.addWidget(self.progress_bar)
        
        layout.addWidget(control_group)
        
        # Results section
        self.results_area = QScrollArea()
        self.results_widget = QWidget()
        self.results_layout = QVBoxLayout(self.results_widget)
        self.results_area.setWidget(self.results_widget)
        self.results_area.setWidgetResizable(True)
        
        # Set responsive minimum size for results area
        self.results_area.setMinimumHeight(400)  # Reduced for smaller screens
        
        layout.addWidget(self.results_area, 1)
        
        # Initially show placeholder
        self.show_placeholder()
    
    def show_placeholder(self):
        """Show placeholder text when no analysis has been performed."""
        matplotlib_note = ""
        if not MATPLOTLIB_AVAILABLE:
            matplotlib_note = "<p><b>Note:</b> Install matplotlib for enhanced visualizations: <code>pip install matplotlib</code></p>"
        
        placeholder = QLabel(
            "<h3>Structural Analysis Tools</h3>"
            "<p>Load a protein structure and select analysis options to get:</p>"
            "<ul>"
            "<li><b>Basic Properties:</b> Chain composition, molecular weight</li>"
            "<li><b>Secondary Structure:</b> Helix/sheet/loop content, Ramachandran plots</li>"
            "<li><b>Geometric Analysis:</b> Bond lengths, angles, radius of gyration</li>"


            "</ul>"
            "<p><i>Tip: Enter a PDB ID to fetch from the database, load from the 3D viewer, or select a local PDB file.</i></p>"
            + matplotlib_note +
            "<p><i>Features include interactive charts, plots, and comprehensive data tables.</i></p>"
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
                    self.results_layout.removeItem(item)
        
        self.results_layout.addWidget(placeholder)
        self.results_layout.addStretch()
    
    def load_current_structure(self):
        """Load the currently displayed structure."""
        if not hasattr(self.parent_app, 'current_structure_id') or not self.parent_app.current_structure_id:
            QMessageBox.warning(self, "No Structure", "No structure is currently loaded in the 3D viewer.")
            return
        
        structure_id = self.parent_app.current_structure_id
        structure_path = os.path.join(self.parent_app.pulled_structures_dir, f"{structure_id}.pdb")
        
        if not os.path.exists(structure_path):
            QMessageBox.warning(self, "File Not Found", f"Structure file not found: {structure_path}")
            return
        
        self.current_structure_path = structure_path
        self.structure_label.setText(f"Current: {structure_id}")
        self.structure_label.setStyleSheet("color: #000; font-weight: bold;")
        
        # Try to fetch comprehensive metadata for this structure
        self._fetch_comprehensive_metadata_for_structure(structure_id)
        
        if hasattr(self.parent_app, 'statusBar'):
            self.parent_app.statusBar().showMessage(f"Loaded structure {structure_id} for analysis")
    
    def fetch_pdb_id(self):
        """Fetch a PDB structure by ID using enhanced puller."""
        pdb_id = self.pdb_id_entry.text().strip().upper()
        if not pdb_id:
            QMessageBox.warning(self, "No PDB ID", "Please enter a PDB ID.")
            return
        
        if not BIOPYTHON_AVAILABLE:
            QMessageBox.critical(
                self, "Missing Dependency", 
                "Biopython is required for PDB fetching.\n\n"
                "Please install it with: pip install biopython"
            )
            return
        
        try:
            if hasattr(self.parent_app, 'statusBar'):
                self.parent_app.statusBar().showMessage(f"🔄 Fetching comprehensive metadata for PDB ID {pdb_id}...")
            
            # Always fetch comprehensive metadata for structural analysis
            # (No popup - using status bar for progress updates)
            
            # Use the comprehensive GraphQL-based puller for advanced metadata
            try:
                from .enhanced_pdb_puller_fixed import EnhancedPDBPuller
            except ImportError:
                # Fallback to absolute import if relative import fails
                import sys
                # Add the project root to Python path
                project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
                if project_root not in sys.path:
                    sys.path.insert(0, project_root)
                from src.core.enhanced_pdb_puller_fixed import EnhancedPDBPuller
            
            temp_puller = EnhancedPDBPuller(self.pulled_structures_dir)
            comprehensive_data = temp_puller.fetch_comprehensive_pdb_data(
                pdb_id,
                include_validation=True,   # Include validation for comprehensive analysis
                include_sequences=True,
                include_mmcif=False
            )
            
            if comprehensive_data['errors']:
                error_details = "\n".join(comprehensive_data['errors'])
                raise Exception(f"Errors occurred during data fetching:\n{error_details}")
            
            pdb_path = comprehensive_data['files'].get('pdb')
            if not pdb_path or not os.path.exists(pdb_path):
                raise FileNotFoundError(f"PDB file not found for {pdb_id}")
            
            self.current_comprehensive_data = comprehensive_data
            
            # Also set in parent app for other tabs to access
            if hasattr(self.parent_app, '__dict__'):
                self.parent_app.current_pdb_data = comprehensive_data
            
            # Set as current structure
            self.current_structure_path = pdb_path
            self.structure_label.setText(f"PDB: {pdb_id}")
            self.structure_label.setStyleSheet("color: #000; font-weight: bold;")
            
            if hasattr(self.parent_app, 'statusBar'):
                self.parent_app.statusBar().showMessage(f"✅ Successfully fetched comprehensive metadata for {pdb_id}")
            
            # Success - comprehensive data fetched (no popup, just status bar update)
            
        except Exception as e:
            if hasattr(self.parent_app, 'statusBar'):
                self.parent_app.statusBar().showMessage(f"Error fetching PDB ID {pdb_id}")
            
            QMessageBox.critical(
                self, "Error Fetching PDB", 
                f"Could not fetch data for PDB ID {pdb_id}.\n\n"
                f"Please check your internet connection and verify the PDB ID is correct.\n\n"
                f"Error details: {str(e)}"
            )
    
    def load_structure_file(self):
        """Load a PDB file for analysis."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open PDB File", "", 
            "PDB Files (*.pdb *.ent);;All Files (*)"
        )
        
        if file_path:
            self.current_structure_path = file_path
            filename = os.path.basename(file_path)
            self.structure_label.setText(f"File: {filename}")
            self.structure_label.setStyleSheet("color: #000; font-weight: bold;")
            
            # Try to extract PDB ID from filename and fetch comprehensive metadata
            pdb_id = self._extract_pdb_id_from_filename(filename)
            if pdb_id:
                self._fetch_comprehensive_metadata_for_structure(pdb_id)
            else:
                # No PDB ID found, clear comprehensive data
                self.current_comprehensive_data = None
            
            if hasattr(self.parent_app, 'statusBar'):
                self.parent_app.statusBar().showMessage(f"Loaded {filename} for analysis")
    
    def _extract_pdb_id_from_filename(self, filename):
        """Extract PDB ID from filename if possible."""
        import re
        
        # Remove file extension
        base_name = os.path.splitext(filename)[0]
        
        # Common PDB filename patterns
        patterns = [
            r'^([0-9][a-zA-Z0-9]{3})$',  # Standard PDB ID (e.g., 1CRN)
            r'pdb([0-9][a-zA-Z0-9]{3})\.',  # pdb1xxx.ent or similar
            r'([0-9][a-zA-Z0-9]{3})_model'  # 1xxx_model format
        ]
        
        for pattern in patterns:
            match = re.search(pattern, base_name, re.IGNORECASE)
            if match:
                return match.group(1).upper()
                
        return None
    
    def _fetch_comprehensive_metadata_for_structure(self, structure_id):
        """Fetch comprehensive metadata for a structure ID."""
        try:

            
            # Initialize comprehensive data if not exists
            if not hasattr(self, 'current_comprehensive_data'):
                self.current_comprehensive_data = None
            
            # Use optimized fetch manager if available
            if OPTIMIZED_FETCH_AVAILABLE and self.pdb_fetch_manager:

                # Check if data is already cached
                cached_info = self.pdb_fetch_manager.get_structure_info(structure_id)
                if cached_info['available'] and cached_info.get('metadata'):

                    self.current_comprehensive_data = {
                        'pdb_id': structure_id,
                        'files': cached_info['files'],
                        'metadata': cached_info.get('metadata', {}),
                        'sequences': {},
                        'validation': {},
                        'errors': []
                    }
                    
                    # Also set in parent app for other tabs to access
                    if hasattr(self.parent_app, '__dict__'):
                        self.parent_app.current_pdb_data = self.current_comprehensive_data

                    return
            
            # Use enhanced PDB puller if available
            if ENHANCED_PDB_AVAILABLE and self.enhanced_pdb_puller:

                try:
                    comprehensive_data = self.enhanced_pdb_puller.fetch_comprehensive_pdb_data(
                        structure_id,
                        include_validation=True,
                        include_sequences=True,
                        include_mmcif=False  # Skip mmCIF for speed
                    )
                    if comprehensive_data and not comprehensive_data.get('errors'):
                        self.current_comprehensive_data = comprehensive_data
                        
                        # Also set in parent app for other tabs to access
                        if hasattr(self.parent_app, '__dict__'):
                            self.parent_app.current_pdb_data = comprehensive_data

                        return
                    else:
                        pass
                except Exception as e:
                    pass
            # Fallback: create minimal metadata structure

            self.current_comprehensive_data = {
                'pdb_id': structure_id,
                'files': {},
                'metadata': {
                    'entry': {
                        'struct': {
                            'title': f"PDB Structure {structure_id}"
                        },
                        'rcsb_entry_info': {
                            'experimental_method': 'Unknown'
                        }
                    }
                },
                'sequences': {},
                'validation': {},
                'errors': ['Comprehensive metadata not available']
            }
            
            # Also set in parent app for other tabs to access
            if hasattr(self.parent_app, '__dict__'):
                self.parent_app.current_pdb_data = self.current_comprehensive_data

            
        except Exception as e:
            # Set minimal fallback data
            self.current_comprehensive_data = {
                'pdb_id': structure_id,
                'files': {},
                'metadata': {},
                'sequences': {},
                'validation': {},
                'errors': [f'Error fetching metadata: {str(e)}']
            }
            
            # Also set in parent app for other tabs to access
            if hasattr(self.parent_app, '__dict__'):
                self.parent_app.current_pdb_data = self.current_comprehensive_data
        
        if hasattr(self.parent_app, 'statusBar'):
            if self.current_comprehensive_data and self.current_comprehensive_data.get('metadata'):
                self.parent_app.statusBar().showMessage(f"Comprehensive metadata loaded for {structure_id}")
            else:
                self.parent_app.statusBar().showMessage(f"Basic metadata only for {structure_id}")
    

    
    def on_analysis_complete(self, results):
        """Handle completion of structural analysis."""
        self.progress_bar.setVisible(False)
        self.current_results = results
        self.export_btn.setEnabled(True)
        
        # Display results
        self.display_results(results)
        
        if hasattr(self.parent_app, 'statusBar'):
            self.parent_app.statusBar().showMessage("Structural analysis complete!")
    
    def on_analysis_error(self, error_message):
        """Handle analysis error."""
        self.progress_bar.setVisible(False)
        QMessageBox.critical(self, "Analysis Error", error_message)
        
        if hasattr(self.parent_app, 'statusBar'):
            self.parent_app.statusBar().showMessage("Analysis failed")
    
    def on_progress_update(self, message):
        """Handle progress updates from analysis worker."""
        if hasattr(self.parent_app, 'statusBar'):
            self.parent_app.statusBar().showMessage(message)
        
    def analyze_structure(self):
        """Perform structural analysis on the current structure."""
        if not self.current_structure_path:
            QMessageBox.warning(self, "No Structure", "Please load a structure first.")
            return
        
        if not BIOPYTHON_AVAILABLE:
            QMessageBox.critical(
                self, "Missing Dependency", 
                "Biopython is required for structural analysis.\n\n"
                "Please install it with: pip install biopython"
            )
            return
        
        # Get selected analysis types
        analysis_types = []
        if self.basic_checkbox.isChecked():
            analysis_types.append('basic')
        if self.secondary_checkbox.isChecked():
            analysis_types.append('secondary')
        if self.geometry_checkbox.isChecked():
            analysis_types.append('geometry')


        
        if not analysis_types:
            QMessageBox.warning(self, "No Analysis Selected", "Please select at least one analysis type.")
            return
        
        # Show progress
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        
        # Update comprehensive data from parent app before analysis
        self.update_comprehensive_data()
        
        # Start analysis in worker thread with comprehensive metadata if available
        comprehensive_data = getattr(self, 'current_comprehensive_data', None)
        self.analysis_worker = StructuralAnalysisWorker(
            self.current_structure_path, 
            analysis_types, 
            comprehensive_data=comprehensive_data
        )
        

        self.analysis_worker.analysis_complete.connect(self.display_results)
        self.analysis_worker.error_occurred.connect(self.handle_analysis_error)
        self.analysis_worker.progress_update.connect(self.update_progress)
        self.analysis_worker.start()
    

    
    def update_progress(self, message):
        """Update progress message."""
        if hasattr(self.parent_app, 'statusBar'):
            self.parent_app.statusBar().showMessage(message)
    
    def handle_analysis_error(self, error_message):
        """Handle analysis errors."""
        self.progress_bar.setVisible(False)
        QMessageBox.critical(self, "Analysis Error", error_message)
        if hasattr(self.parent_app, 'statusBar'):
            self.parent_app.statusBar().showMessage("Analysis failed")
    
    def create_pie_chart(self, data, title, colors=None):
        """Create a pie chart widget with robust error handling."""
        if not MATPLOTLIB_AVAILABLE or not data:
            return None
        
        try:
            # Filter out zero, negative, and NaN values
            filtered_data = {}
            for k, v in data.items():
                if isinstance(v, (int, float)) and v > 0 and not (math.isnan(v) or math.isinf(v)):
                    filtered_data[k] = v
            
            if not filtered_data:
                # Return a simple label if no valid data
                label = QLabel(f"<i>No data available for {title}</i>")
                label.setAlignment(Qt.AlignCenter)
                label.setStyleSheet("color: #666; padding: 20px;")
                return label
            
            # Ensure we have at least some meaningful data
            total_value = sum(filtered_data.values())
            if total_value <= 0:
                label = QLabel(f"<i>No valid data for {title}</i>")
                label.setAlignment(Qt.AlignCenter)
                label.setStyleSheet("color: #666; padding: 20px;")
                return label
            
            # Use larger figure for amino acid composition
            if "Amino Acid" in title:
                fig = Figure(figsize=(12, 10), dpi=100)
            else:
                fig = Figure(figsize=(8, 6), dpi=100)
            
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)
            
            # Prepare data
            labels = list(filtered_data.keys())
            values = list(filtered_data.values())
            
            # Double-check values are valid
            values = [max(0.001, v) for v in values if not (math.isnan(v) or math.isinf(v))]
            if len(values) != len(labels):
                # Mismatch, rebuild both lists
                valid_items = [(k, v) for k, v in filtered_data.items() 
                              if isinstance(v, (int, float)) and v > 0 and not (math.isnan(v) or math.isinf(v))]
                if not valid_items:
                    label = QLabel(f"<i>No valid data for {title}</i>")
                    label.setAlignment(Qt.AlignCenter)
                    label.setStyleSheet("color: #666; padding: 20px;")
                    return label
                labels, values = zip(*valid_items)
                labels, values = list(labels), list(values)
            
            # Create distinct colors for amino acids
            if colors is None:
                if "Amino Acid" in title and len(labels) <= 20:
                    # Use a colormap that provides distinct colors for all 20 amino acids
                    colors = plt.cm.tab20(range(len(labels)))
                else:
                    colors = plt.cm.Set3(range(len(labels)))
            
            # Special handling for amino acid composition
            if "Amino Acid" in title:
                # Only show percentages for slices > 2% to avoid clutter
                def autopct_func(pct):
                    return f'{pct:.1f}%' if pct > 2 else ''
                
                wedges, texts, autotexts = ax.pie(values, labels=labels, autopct=autopct_func, 
                                                  colors=colors, startangle=90, 
                                                  textprops={'fontsize': 9})
                
                # Position labels outside the pie
                for text in texts:
                    text.set_fontsize(10)
                    text.set_fontweight('bold')
                
                # Make percentage text more readable
                for autotext in autotexts:
                    autotext.set_color('black')
                    autotext.set_fontweight('bold')
                    autotext.set_fontsize(9)
                
                # Add a legend for better readability
                legend_labels = [f'{label}: {value} ({value/sum(values)*100:.1f}%)' 
                               for label, value in zip(labels, values)]
                ax.legend(wedges, legend_labels, title="Amino Acids", 
                         loc="center left", bbox_to_anchor=(1, 0, 0.5, 1),
                         fontsize=9)
            else:
                # Standard pie chart for other data
                wedges, texts, autotexts = ax.pie(values, labels=labels, autopct='%1.1f%%', 
                                                  colors=colors, startangle=90)
                
                # Make percentage text more readable
                for autotext in autotexts:
                    autotext.set_color('white')
                    autotext.set_fontweight('bold')
                    autotext.set_fontsize(10)
                
                # Make labels more readable
                for text in texts:
                    text.set_fontsize(10)
            
            ax.set_title(title, fontsize=14, fontweight='bold')
            
            fig.tight_layout()
            
            if "Amino Acid" in title:
                canvas.setMinimumSize(1200, 900)
            else:
                canvas.setMinimumSize(800, 600)
            
            return canvas
            
        except Exception as e:
            print(f"[DEBUG] Error creating pie chart '{title}': {e}")
            # Return error message widget
            label = QLabel(f"<i>Error creating chart: {title}</i>")
            label.setAlignment(Qt.AlignCenter)
            label.setStyleSheet("color: #cc0000; padding: 20px;")
            return label
    
    def create_bar_chart(self, data, title, xlabel, ylabel):
        """Create a bar chart widget."""
        if not MATPLOTLIB_AVAILABLE:
            return None
        
        fig = Figure(figsize=(14, 8), dpi=100)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        
        labels = list(data.keys())
        values = list(data.values())
        
        bars = ax.bar(labels, values, color='skyblue', alpha=0.7)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.1f}',
                   ha='center', va='bottom', fontsize=10)
        
        # Rotate x-axis labels if needed
        if len(labels) > 10:
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        
        fig.tight_layout()
        canvas.setMinimumSize(1000, 600)
        return canvas
    
    def create_ramachandran_plot(self, ramachandran_data, pdb_id=None):
        """Create a Ramachandran plot using ramachandraw."""
        if not RAMACHANDRAW_AVAILABLE or not MATPLOTLIB_AVAILABLE or not ramachandran_data:
            return None, {}
        
        try:
            import ramachandraw.utils
            import tempfile
            
            # We need to use the current structure file with ramachandraw
            if not hasattr(self, 'current_structure_path') or not self.current_structure_path:
                # Fallback to manual plot if no structure file available
                return self._create_manual_ramachandran_plot(ramachandran_data, pdb_id)
            
            # Create a temporary directory for ramachandraw output
            with tempfile.TemporaryDirectory() as temp_dir:
                # Use ramachandraw to create the plot
                ax = ramachandraw.utils.plot(
                    pdb_filepath=self.current_structure_path,
                    cmap='viridis',
                    alpha=0.7,
                    dpi=100,
                    save=False,  # Don't save to file
                    show=False   # Don't show interactively
                )
                
                # Override the title set by ramachandraw with PDB ID
                if pdb_id:
                    plot_title = f'Ramachandran Plot - {pdb_id}'
                else:
                    plot_title = 'Ramachandran Plot'
                ax.set_title(plot_title, fontsize=16, fontweight='bold', pad=20)
                
                # Get the figure from the axes
                fig = ax.get_figure()
                
                # Create a Qt canvas from the matplotlib figure
                canvas = FigureCanvas(fig)
                canvas.setMinimumSize(800, 600)
                
                # Count residues by region using our classification
                region_counts = {'favored': 0, 'allowed': 0, 'outlier': 0}
                for data in ramachandran_data:
                    region = data.get('region', 'outlier')
                    region_counts[region] = region_counts.get(region, 0) + 1
                
                return canvas, region_counts
            
        except Exception as e:
            print(f"[DEBUG] Error creating Ramachandran plot with ramachandraw: {e}")
            print(f"[DEBUG] Falling back to manual plot")
            # Fallback to manual implementation
            return self._create_manual_ramachandran_plot(ramachandran_data, pdb_id)
    
    def _create_manual_ramachandran_plot(self, ramachandran_data, pdb_id=None):
        """Fallback manual Ramachandran plot creation."""
        try:
            # Create figure
            fig = Figure(figsize=(10, 8), dpi=100, facecolor='white')
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)
            
            # Extract phi and psi angles
            phi_angles = [data['phi'] for data in ramachandran_data]
            psi_angles = [data['psi'] for data in ramachandran_data]
            
            # Plot background regions (simplified)
            # Alpha helix region
            alpha_region = plt.Rectangle((-180, -70), 130, 120, 
                                       alpha=0.3, color='blue', label='Alpha helix')
            ax.add_patch(alpha_region)
            
            # Beta sheet region  
            beta_region = plt.Rectangle((-180, 90), 130, 90, 
                                      alpha=0.3, color='green', label='Beta sheet')
            ax.add_patch(beta_region)
            
            # Plot the data points
            colors = []
            for data in ramachandran_data:
                region = data.get('region', 'outlier')
                if region == 'favored':
                    colors.append('darkblue')
                elif region == 'allowed':
                    colors.append('green')
                else:
                    colors.append('red')
            
            scatter = ax.scatter(phi_angles, psi_angles, c=colors, alpha=0.7, s=30, 
                               edgecolors='black', linewidth=0.5)
            
            # Customize the plot
            # Set title based on PDB ID if available
            if pdb_id:
                plot_title = f'Ramachandran Plot - {pdb_id}'
            else:
                plot_title = 'Ramachandran Plot'
            ax.set_title(plot_title, fontsize=16, fontweight='bold', pad=20)
            ax.set_xlabel('Phi (φ) angle (degrees)', fontsize=14)
            ax.set_ylabel('Psi (ψ) angle (degrees)', fontsize=14)
            ax.grid(True, alpha=0.3)
            
            # Set axis limits
            ax.set_xlim(-180, 180)
            ax.set_ylim(-180, 180)
            
            # Add ticks
            ax.set_xticks(range(-180, 181, 60))
            ax.set_yticks(range(-180, 181, 60))
            
            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='darkblue', alpha=0.7, label='Favored'),
                Patch(facecolor='green', alpha=0.7, label='Allowed'),
                Patch(facecolor='red', alpha=0.7, label='Outlier')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
            
            fig.tight_layout()
            canvas.setMinimumSize(800, 600)
            
            # Count residues by region
            region_counts = {'favored': 0, 'allowed': 0, 'outlier': 0}
            for data in ramachandran_data:
                region = data.get('region', 'outlier')
                region_counts[region] = region_counts.get(region, 0) + 1
            
            return canvas, region_counts
            
        except Exception as e:
            print(f"[DEBUG] Error creating manual Ramachandran plot: {e}")
            # Return error message widget
            label = QLabel(f"<i>Error creating Ramachandran plot: {str(e)}</i>")
            label.setAlignment(Qt.AlignCenter)
            label.setStyleSheet("color: #cc0000; padding: 20px;")
            return label, {}
    

    

    

    
    def is_in_favored_region(self, phi, psi, residue_name=None):
        """Check if phi/psi angles are in favored regions according to MolProbity standard.
        
        Args:
            phi: Phi angle in degrees
            psi: Psi angle in degrees
            residue_name: Optional three-letter residue code for special cases
            
        Returns:
            bool: True if in favored region, False otherwise
        """
        # Handle special cases first
        if residue_name == 'GLY':
            return self._is_glycine_favored(phi, psi)
        elif residue_name == 'PRO':
            return self._is_proline_favored(phi, psi)
            
        # General case for standard amino acids
        # Core alpha-helix region (αR)
        if (-100 <= phi <= -30 and -60 <= psi <= 10):
            return True
            
        # Core beta-sheet region (β)
        if (-160 <= phi <= -60 and 90 <= psi <= 180):  # Upper beta
            return True
        if (-160 <= phi <= -60 and -180 <= psi <= -80):  # Lower beta
            return True
            
        # Core left-handed alpha-helix (αL)
        if (30 <= phi <= 80 and 30 <= psi <= 80):
            return True
            
        return False
    
    def is_in_allowed_region(self, phi, psi, residue_name=None):
        """Check if phi/psi angles are in allowed regions according to MolProbity standard.
        
        Args:
            phi: Phi angle in degrees
            psi: Psi angle in degrees
            residue_name: Optional three-letter residue code for special cases
            
        Returns:
            bool: True if in allowed region, False otherwise
        """
        # Handle special cases first
        if residue_name == 'GLY':
            return self._is_glycine_allowed(phi, psi)
        elif residue_name == 'PRO':
            return self._is_proline_allowed(phi, psi)
            
        # General case for standard amino acids
        # If already in favored region, it's also in allowed region
        if self.is_in_favored_region(phi, psi, residue_name):
            return True
            
        # Allowed alpha-helix region
        if (-180 <= phi <= -25 and 20 <= psi <= 80):
            return True
        if (-25 <= phi <= 0 and -80 <= psi <= 80):
            return True
            
        # Allowed beta-sheet region
        if (-100 <= phi <= -25 and 80 <= psi <= 180):
            return True
        if (-100 <= phi <= -25 and -180 <= psi <= -100):
            return True
            
        # Left-handed alpha-helix allowed region
        if (0 <= phi <= 90 and 0 <= psi <= 90):
            return True
            
        # Bridge region (connects alpha and beta)
        if (-100 <= phi <= -25 and -100 <= psi <= 80):
            return True
            
        return False
    
    def _is_glycine_favored(self, phi, psi):
        """Check if glycine phi/psi angles are in favored regions."""
        # Glycine has much larger allowed regions due to its flexibility
        # Core alpha-helix region (αR)
        if (-180 <= phi <= -25 and -120 <= psi <= 80):
            return True
            
        # Core beta-sheet region (β)
        if (-180 <= phi <= -25 and 80 <= psi <= 180):
            return True
            
        # Core left-handed alpha-helix (αL)
        if (25 <= phi <= 180 and 0 <= psi <= 180):
            return True
            
        return False
    
    def _is_glycine_allowed(self, phi, psi):
        """Check if glycine phi/psi angles are in allowed regions."""
        # For glycine, almost all regions are allowed
        # Only exclude the most sterically disfavored regions
        if (-180 <= phi <= 180 and -180 <= psi <= 180):
            # Only exclude the most disfavored region
            if not (-30 <= phi <= 30 and -60 <= psi <= -30):
                return True
        return False
    
    def _is_proline_favored(self, phi, psi):
        """Check if proline phi/psi angles are in favored regions."""
        # Proline is highly constrained due to its ring structure
        # Core alpha-helix region (αR)
        if (-75 <= phi <= -45 and -60 <= psi <= -10):
            return True
            
        # Polyproline II helix (common for proline)
        if (-75 <= phi <= -45 and 120 <= psi <= 180):
            return True
            
        return False
    
    def _is_proline_allowed(self, phi, psi):
        """Check if proline phi/psi angles are in allowed regions."""
        # Proline has very restricted allowed regions
        if self._is_proline_favored(phi, psi):
            return True
            
        # Slightly extended alpha-helix region
        if (-90 <= phi <= -30 and -90 <= psi <= 30):
            return True
            
        # Extended polyproline II region
        if (-90 <= phi <= -30 and 90 <= psi <= 180):
            return True
            
        return False
    
    def create_histogram(self, data, title, xlabel, ylabel, bins=20):
        """Create a histogram widget."""
        if not MATPLOTLIB_AVAILABLE or not data:
            return None
        
        fig = Figure(figsize=(10, 6), dpi=100)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        
        n, bins_edges, patches_hist = ax.hist(data, bins=bins, alpha=0.7, color='skyblue', edgecolor='black')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        mean_val = np.mean(data)
        std_val = np.std(data)
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2f}')
        ax.axvline(mean_val + std_val, color='orange', linestyle=':', linewidth=1, label=f'+1σ: {mean_val + std_val:.2f}')
        ax.axvline(mean_val - std_val, color='orange', linestyle=':', linewidth=1, label=f'-1σ: {mean_val - std_val:.2f}')
        ax.legend(fontsize=10)
        
        # Improve tick labels
        ax.tick_params(axis='both', which='major', labelsize=10)
        
        fig.tight_layout()
        canvas.setMinimumSize(1000, 600)
        return canvas
    
    def display_enhanced_metadata(self, basic_results):
        """Display simplified structure information with essential details only."""
        try:
            comprehensive_data = self.current_comprehensive_data
            if not comprehensive_data:

                self.display_basic_results_fallback(basic_results)
                return
            
            metadata = comprehensive_data.get('metadata', {})
            
            # If no metadata, try to show what we have
            if not metadata:

                metadata = {
                    'entry': {
                        'struct': {
                            'title': f"PDB Structure {basic_results.get('structure_id', 'Unknown')}"
                        },
                        'rcsb_entry_info': {
                            'experimental_method': 'Unknown'
                        }
                    }
                }
        except Exception as e:

            self.display_basic_results_fallback(basic_results)
            return
        
        # Simplified Structure Information with metadata source indicator
        metadata_source = "📊 Comprehensive Metadata" if metadata.get('entry') else "📄 Basic PDB Header"
        enhanced_basic_group = QGroupBox(f"Structure Information ({metadata_source})")
        enhanced_basic_layout = QVBoxLayout(enhanced_basic_group)
        
        # Create simplified tabs with only essential information
        metadata_tabs = QTabWidget()
        
        # 1. Basic Information Tab (Simplified)
        self._create_simple_basic_tab(metadata_tabs, basic_results, metadata)
        
        # 2. Experimental Details Tab - Removed
        
        # 3. Sequence Information Tab (If available)
        if basic_results.get('chains'):
            self._create_simple_sequence_tab(metadata_tabs, basic_results)
        
        enhanced_basic_layout.addWidget(metadata_tabs)
        self.results_layout.addWidget(enhanced_basic_group)
        
        # Add chain comparison charts if multiple chains
        if basic_results.get('chains') and len(basic_results['chains']) > 1:
            chain_charts_group = QGroupBox("Chain Size Comparison")
            chain_charts_layout = QHBoxLayout(chain_charts_group)
            
            # Create bar charts for chain comparison
            chain_residues = {f"Chain {chain['id']}": chain['residues'] for chain in basic_results['chains']}
            chain_atoms = {f"Chain {chain['id']}": chain['atoms'] for chain in basic_results['chains']}
            
            residue_chart = self.create_bar_chart(chain_residues, "Residues per Chain", "Chain", "Number of Residues")
            if residue_chart:
                chain_charts_layout.addWidget(residue_chart)
            
            atom_chart = self.create_bar_chart(chain_atoms, "Atoms per Chain", "Chain", "Number of Atoms")
            if atom_chart:
                chain_charts_layout.addWidget(atom_chart)
            
            self.results_layout.addWidget(chain_charts_group)
        
        # Also display standard composition analysis
        if basic_results.get('composition'):
            self.display_composition_analysis(basic_results)
    
    def _create_simple_basic_tab(self, metadata_tabs, basic_results, metadata):
        """Create simplified basic information tab."""
        basic_tab = QWidget()
        basic_tab_layout = QVBoxLayout(basic_tab)
        
        # Essential Structure Information
        info_group = QGroupBox("Structure Information")
        info_layout = QFormLayout(info_group)
        
        # Basic structure data
        info_layout.addRow("Chains:", QLabel(str(len(basic_results.get('chains', [])))))
        info_layout.addRow("Residues:", QLabel(str(basic_results.get('total_residues', 0))))
        info_layout.addRow("Atoms:", QLabel(str(basic_results.get('total_atoms', 0))))
        info_layout.addRow("Molecular Weight:", QLabel(f"{basic_results.get('molecular_weight', 0):.1f} Da"))
        
        # Title from metadata if available
        if 'entry' in metadata:
            entry = metadata['entry']
            if entry.get('struct') and entry['struct'].get('title'):
                title = entry['struct']['title']
                title_label = QLabel(title)
                title_label.setWordWrap(True)
                # No width limit - allow full title display
                info_layout.addRow("Title:", title_label)
            
            # Authors (first 3 only)
            if entry.get('audit_author'):
                authors = [author.get('name', '') for author in entry['audit_author'][:3] if author.get('name')]
                if authors:
                    authors_text = ', '.join(authors)
                    if len(entry.get('audit_author', [])) > 3:
                        authors_text += " et al."
                    info_layout.addRow("Authors:", QLabel(authors_text))
            
            # Key dates
            if entry.get('rcsb_accession_info'):
                accession_info = entry['rcsb_accession_info']
                deposit_date = accession_info.get('deposit_date', 'N/A')
                release_date = accession_info.get('initial_release_date', 'N/A')
                info_layout.addRow("Deposited:", QLabel(deposit_date))
                info_layout.addRow("Released:", QLabel(release_date))
            
            # Experimental method and resolution
            if entry.get('rcsb_entry_info'):
                rcsb_info = entry['rcsb_entry_info']
                
                # Experimental method
                experimental_method = rcsb_info.get('experimental_method', 'N/A')
                info_layout.addRow("Experimental Method:", QLabel(experimental_method))
                
                # Resolution
                resolution = rcsb_info.get('resolution_combined', 'N/A')
                if isinstance(resolution, list) and resolution:
                    resolution = resolution[0]
                if resolution != 'N/A' and isinstance(resolution, (int, float)):
                    resolution_text = f"{resolution:.2f} Å"
                else:
                    resolution_text = str(resolution)
                info_layout.addRow("Resolution:", QLabel(resolution_text))
        
        basic_tab_layout.addWidget(info_group)
        basic_tab_layout.addStretch()
        metadata_tabs.addTab(basic_tab, "Basic Info")
    
    # Experimental tab creation method removed
    
    def _create_simple_sequence_tab(self, metadata_tabs, basic_results):
        """Create simplified sequence information tab."""
        sequence_tab = QWidget()
        sequence_tab_layout = QVBoxLayout(sequence_tab)
        
        # Chain Sequences
        sequences_group = QGroupBox("Chain Sequences")
        sequences_layout = QVBoxLayout(sequences_group)
        
        chains = basic_results.get('chains', [])
        for chain in chains:
            chain_info = QLabel(f"Chain {chain['id']}: {chain['residues']} residues")
            sequences_layout.addWidget(chain_info)
            
            # Display sequence if available
            if chain.get('sequence'):
                sequence = chain['sequence']  # Show full sequence
                
                # Format sequence in blocks of 50
                formatted_sequence = '\n'.join([sequence[i:i+50] for i in range(0, len(sequence), 50)])
                sequence_text = QTextEdit()
                sequence_text.setPlainText(formatted_sequence)
                sequence_text.setMaximumHeight(80)
                sequence_text.setReadOnly(True)
                sequence_text.setFont(QFont("Courier", 9))
                sequences_layout.addWidget(sequence_text)
        
        sequence_tab_layout.addWidget(sequences_group)
        sequence_tab_layout.addStretch()
        metadata_tabs.addTab(sequence_tab, "Sequences")
    
    def _create_basic_info_tab(self, metadata_tabs, basic_results, metadata):
        """Create enhanced basic information tab."""
        basic_tab = QWidget()
        basic_tab_layout = QVBoxLayout(basic_tab)
        
        # Structure Overview Section
        overview_group = QGroupBox("Structure Overview")
        overview_layout = QFormLayout(overview_group)
        
        # Standard basic info
        overview_layout.addRow("PDB ID:", QLabel(basic_results.get('structure_id', 'N/A')))
        overview_layout.addRow("Total Chains:", QLabel(str(len(basic_results.get('chains', [])))))
        overview_layout.addRow("Total Residues:", QLabel(str(basic_results.get('total_residues', 0))))
        overview_layout.addRow("Total Atoms:", QLabel(str(basic_results.get('total_atoms', 0))))
        overview_layout.addRow("Molecular Weight:", QLabel(f"{basic_results.get('molecular_weight', 0):.1f} Da"))
        
        # Enhanced info from metadata
        if 'entry' in metadata:
            entry = metadata['entry']
            
            # Title and description
            title = 'N/A'
            if entry.get('struct') and entry['struct'].get('title'):
                title = entry['struct']['title']
            title_label = QLabel(title)
            title_label.setWordWrap(True)
            # No width limit - allow full title display
            overview_layout.addRow("Title:", title_label)
            
            # Classification
            if entry.get('struct_keywords') and entry['struct_keywords'].get('pdbx_keywords'):
                keywords = entry['struct_keywords']['pdbx_keywords']
                keywords_label = QLabel(keywords)
                keywords_label.setWordWrap(True)
                keywords_label.setMaximumWidth(400)
                overview_layout.addRow("Keywords:", keywords_label)
            
            # Organism information
            if entry.get('rcsb_entry_info') and entry['rcsb_entry_info'].get('selected_polymer_entity_types'):
                entity_types = ', '.join(entry['rcsb_entry_info']['selected_polymer_entity_types'])
                overview_layout.addRow("Entity Types:", QLabel(entity_types))
        
        basic_tab_layout.addWidget(overview_group)
        
        # Authors and Dates Section
        authors_group = QGroupBox("Authors and Dates")
        authors_layout = QFormLayout(authors_group)
        
        if 'entry' in metadata:
            entry = metadata['entry']
            
            # Authors
            authors = []
            if entry.get('audit_author'):
                authors = [author.get('name', '') for author in entry['audit_author'] if author.get('name')]
            if authors:
                authors_text = ', '.join(authors[:5])  # Show first 5 authors
                if len(authors) > 5:
                    authors_text += f" and {len(authors) - 5} others"
                authors_label = QLabel(authors_text)
                authors_label.setWordWrap(True)
                authors_label.setMaximumWidth(400)
                authors_layout.addRow("Authors:", authors_label)
            
            # Dates
            if entry.get('rcsb_accession_info'):
                accession_info = entry['rcsb_accession_info']
                deposit_date = accession_info.get('deposit_date', 'N/A')
                release_date = accession_info.get('initial_release_date', 'N/A')
                revision_date = accession_info.get('revision_date', 'N/A')
                
                authors_layout.addRow("Deposition Date:", QLabel(deposit_date))
                authors_layout.addRow("Release Date:", QLabel(release_date))
                if revision_date != 'N/A':
                    authors_layout.addRow("Last Revision:", QLabel(revision_date))
        
        basic_tab_layout.addWidget(authors_group)
        basic_tab_layout.addStretch()
        metadata_tabs.addTab(basic_tab, "Basic Info")
    

    
    def _create_quality_metrics_tab(self, metadata_tabs, metadata):
        """Create quality metrics and validation tab."""
        quality_tab = QWidget()
        quality_tab_layout = QVBoxLayout(quality_tab)
        
        if 'entry' not in metadata:
            quality_tab_layout.addWidget(QLabel("No quality data available"))
            metadata_tabs.addTab(quality_tab, "Quality")
            return
        
        entry = metadata['entry']
        
        # Overall Quality Section
        overall_group = QGroupBox("Overall Quality Assessment")
        overall_layout = QFormLayout(overall_group)
        
        # Resolution-based quality assessment
        if entry.get('rcsb_entry_info') and entry['rcsb_entry_info'].get('resolution_combined'):
            res_data = entry['rcsb_entry_info']['resolution_combined']
            if isinstance(res_data, list) and res_data:
                resolution = res_data[0]
            elif isinstance(res_data, (int, float)):
                resolution = res_data
            else:
                resolution = None
            
            if resolution:
                if resolution <= 1.5:
                    quality = "Excellent (≤1.5Å)"
                elif resolution <= 2.0:
                    quality = "Very Good (1.5-2.0Å)"
                elif resolution <= 2.5:
                    quality = "Good (2.0-2.5Å)"
                elif resolution <= 3.0:
                    quality = "Moderate (2.5-3.0Å)"
                else:
                    quality = "Low (>3.0Å)"
                overall_layout.addRow("Resolution Quality:", QLabel(quality))
        
        # R-factor based assessment
        if entry.get('refine') and entry['refine']:
            refine_data = entry['refine'][0]
            r_work = refine_data.get('ls_R_factor_R_work')
            r_free = refine_data.get('ls_R_factor_R_free')
            
            if r_work and r_free:
                r_gap = r_free - r_work
                if r_gap < 0.05:
                    r_quality = "Excellent (gap < 0.05)"
                elif r_gap < 0.07:
                    r_quality = "Good (gap < 0.07)"
                elif r_gap < 0.10:
                    r_quality = "Acceptable (gap < 0.10)"
                else:
                    r_quality = "Poor (gap ≥ 0.10)"
                overall_layout.addRow("R-factor Quality:", QLabel(r_quality))
                overall_layout.addRow("R-free - R-work:", QLabel(f"{r_gap:.3f}"))
        
        quality_tab_layout.addWidget(overall_group)
        
        # Validation Information
        validation_group = QGroupBox("Validation Information")
        validation_layout = QVBoxLayout(validation_group)
        validation_layout.addWidget(QLabel("Detailed validation data would be available from wwPDB validation reports."))
        validation_layout.addWidget(QLabel("This includes Ramachandran plot analysis, bond geometry, and clash analysis."))
        
        quality_tab_layout.addWidget(validation_group)
        quality_tab_layout.addStretch()
        metadata_tabs.addTab(quality_tab, "Quality")
    
    def _create_biological_assembly_tab(self, metadata_tabs, metadata):
        """Create biological assembly information tab."""
        assembly_tab = QWidget()
        assembly_tab_layout = QVBoxLayout(assembly_tab)
        
        if 'entry' not in metadata:
            assembly_tab_layout.addWidget(QLabel("No biological assembly data available"))
            metadata_tabs.addTab(assembly_tab, "Assembly")
            return
        
        entry = metadata['entry']
        
        # Assembly Information
        assembly_group = QGroupBox("Biological Assembly Information")
        assembly_layout = QFormLayout(assembly_group)
        
        # Basic assembly info
        if entry.get('rcsb_entry_info'):
            entry_info = entry['rcsb_entry_info']
            
            # Assembly count
            if entry_info.get('assembly_count'):
                assembly_layout.addRow("Number of Assemblies:", QLabel(str(entry_info['assembly_count'])))
            
            # Polymer entity instances
            polymer_count = entry_info.get('deposited_polymer_entity_instance_count', 0)
            if polymer_count > 0:
                assembly_layout.addRow("Polymer Chains:", QLabel(str(polymer_count)))
            
            # Non-polymer entities
            nonpolymer_count = entry_info.get('deposited_nonpolymer_entity_instance_count', 0)
            if nonpolymer_count > 0:
                assembly_layout.addRow("Non-polymer Entities:", QLabel(str(nonpolymer_count)))
        
        # Symmetry information
        if entry.get('symmetry'):
            symmetry = entry['symmetry']
            if symmetry.get('space_group_name_H_M'):
                assembly_layout.addRow("Space Group:", QLabel(symmetry['space_group_name_H_M']))
            if symmetry.get('Int_Tables_number'):
                assembly_layout.addRow("Space Group Number:", QLabel(str(symmetry['Int_Tables_number'])))
        
        assembly_tab_layout.addWidget(assembly_group)
        
        # Quaternary Structure
        quaternary_group = QGroupBox("Quaternary Structure")
        quaternary_layout = QVBoxLayout(quaternary_group)
        quaternary_layout.addWidget(QLabel("Detailed quaternary structure analysis would include:"))
        quaternary_layout.addWidget(QLabel("• Oligomeric state and stoichiometry"))
        quaternary_layout.addWidget(QLabel("• Interface analysis between chains"))
        quaternary_layout.addWidget(QLabel("• Symmetry operations and transformations"))
        
        assembly_tab_layout.addWidget(quaternary_group)
        assembly_tab_layout.addStretch()
        metadata_tabs.addTab(assembly_tab, "Assembly")
    
    def _create_ligands_tab(self, metadata_tabs, metadata):
        """Create ligands and binding sites tab."""
        ligands_tab = QWidget()
        ligands_tab_layout = QVBoxLayout(ligands_tab)
        
        if 'entry' not in metadata:
            ligands_tab_layout.addWidget(QLabel("No ligand data available"))
            metadata_tabs.addTab(ligands_tab, "Ligands")
            return
        
        entry = metadata['entry']
        
        # Ligand Information
        ligands_group = QGroupBox("Ligand and Small Molecule Information")
        ligands_layout = QVBoxLayout(ligands_group)
        
        # Non-polymer entity count
        if entry.get('rcsb_entry_info'):
            entry_info = entry['rcsb_entry_info']
            nonpolymer_count = entry_info.get('deposited_nonpolymer_entity_instance_count', 0)
            
            if nonpolymer_count > 0:
                count_label = QLabel(f"This structure contains {nonpolymer_count} non-polymer entity instances.")
                ligands_layout.addWidget(count_label)
                
                info_label = QLabel("Non-polymer entities may include:")
                ligands_layout.addWidget(info_label)
                
                types_label = QLabel("• Small molecule ligands and inhibitors\n"
                                   "• Cofactors and prosthetic groups\n"
                                   "• Metal ions and coordination complexes\n"
                                   "• Solvent molecules and crystallization additives")
                ligands_layout.addWidget(types_label)
            else:
                no_ligands_label = QLabel("No non-polymer entities detected in this structure.")
                ligands_layout.addWidget(no_ligands_label)
        
        ligands_tab_layout.addWidget(ligands_group)
        
        # Binding Sites Information
        binding_group = QGroupBox("Binding Sites Analysis")
        binding_layout = QVBoxLayout(binding_group)
        binding_layout.addWidget(QLabel("Detailed binding site analysis would include:"))
        binding_layout.addWidget(QLabel("• Identification of ligand binding pockets"))
        binding_layout.addWidget(QLabel("• Analysis of protein-ligand interactions"))
        binding_layout.addWidget(QLabel("• Binding affinity and thermodynamic data"))
        binding_layout.addWidget(QLabel("• Comparison with known binding sites"))
        
        ligands_tab_layout.addWidget(binding_group)
        ligands_tab_layout.addStretch()
        metadata_tabs.addTab(ligands_tab, "Ligands")
    
    def _create_sequence_tab(self, metadata_tabs, metadata, basic_results):
        """Create sequence information tab."""
        sequence_tab = QWidget()
        sequence_tab_layout = QVBoxLayout(sequence_tab)
        
        # Chain Sequences
        sequences_group = QGroupBox("Chain Sequences")
        sequences_layout = QVBoxLayout(sequences_group)
        
        chains = basic_results.get('chains', [])
        if chains:
            for i, chain in enumerate(chains):
                chain_group = QGroupBox(f"Chain {chain['id']}")
                chain_layout = QFormLayout(chain_group)
                
                chain_layout.addRow("Length:", QLabel(f"{chain['residues']} residues"))
                chain_layout.addRow("Atoms:", QLabel(str(chain['atoms'])))
                
                # Display sequence if available
                if chain.get('sequence'):
                    sequence = chain['sequence']
                    # Format sequence in blocks of 50
                    formatted_sequence = '\n'.join([sequence[i:i+50] for i in range(0, len(sequence), 50)])
                    sequence_text = QTextEdit()
                    sequence_text.setPlainText(formatted_sequence)
                    sequence_text.setMaximumHeight(100)
                    sequence_text.setReadOnly(True)
                    sequence_text.setFont(QFont("Courier", 9))
                    chain_layout.addRow("Sequence:", sequence_text)
                
                sequences_layout.addWidget(chain_group)
        else:
            sequences_layout.addWidget(QLabel("No chain sequence information available"))
        
        sequence_tab_layout.addWidget(sequences_group)
        
        # Sequence Features
        features_group = QGroupBox("Sequence Features")
        features_layout = QVBoxLayout(features_group)
        features_layout.addWidget(QLabel("Advanced sequence analysis would include:"))
        features_layout.addWidget(QLabel("• Domain and motif identification"))
        features_layout.addWidget(QLabel("• Post-translational modifications"))
        features_layout.addWidget(QLabel("• Signal peptides and transmembrane regions"))
        features_layout.addWidget(QLabel("• Homology and evolutionary relationships"))
        
        sequence_tab_layout.addWidget(features_group)
        sequence_tab_layout.addStretch()
        metadata_tabs.addTab(sequence_tab, "Sequences")
    
    def _create_cross_references_tab(self, metadata_tabs, metadata, basic_results):
        """Create cross-references tab."""
        xref_tab = QWidget()
        xref_tab_layout = QVBoxLayout(xref_tab)
        
        # Database Cross-References
        xref_group = QGroupBox("Database Cross-References")
        xref_layout = QVBoxLayout(xref_group)
        
        if 'entry' in metadata:
            entry = metadata['entry']
            pdb_id = basic_results.get('structure_id', 'N/A')
            
            # Create clickable links (conceptual - would need proper implementation)
            links_text = f"""Related database entries for {pdb_id}:
            
• RCSB PDB: https://www.rcsb.org/structure/{pdb_id}
• PDBe: https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}
• PDBj: https://pdbj.org/mine/summary/{pdb_id}
• UniProt: Protein sequence database entries
• Pfam: Protein family classifications
• SCOP: Structural classification
• CATH: Class, Architecture, Topology, Homology
• ChEMBL: Bioactivity database
• DrugBank: Drug and drug target database"""
            
            links_label = QLabel(links_text)
            links_label.setWordWrap(True)
            xref_layout.addWidget(links_label)
        else:
            xref_layout.addWidget(QLabel("No cross-reference data available"))
        
        xref_tab_layout.addWidget(xref_group)
        xref_tab_layout.addStretch()
        metadata_tabs.addTab(xref_tab, "Cross-Refs")
    
    def _create_citations_tab(self, metadata_tabs, metadata):
        """Create citations and publications tab."""
        citations_tab = QWidget()
        citations_tab_layout = QVBoxLayout(citations_tab)
        
        # Publications
        pubs_group = QGroupBox("Publications and Citations")
        pubs_layout = QVBoxLayout(pubs_group)
        
        if 'entry' in metadata:
            entry = metadata['entry']
            
            # Primary citation information
            if entry.get('citation'):
                citations = entry['citation']
                for i, citation in enumerate(citations[:3]):  # Show first 3 citations
                    citation_group = QGroupBox(f"Citation {i+1}")
                    citation_layout = QFormLayout(citation_group)
                    
                    if citation.get('title'):
                        title_label = QLabel(citation['title'])
                        title_label.setWordWrap(True)
                        # No width limit - allow full citation title display
                        citation_layout.addRow("Title:", title_label)
                    
                    if citation.get('journal_abbrev'):
                        citation_layout.addRow("Journal:", QLabel(citation['journal_abbrev']))
                    
                    if citation.get('year'):
                        citation_layout.addRow("Year:", QLabel(str(citation['year'])))
                    
                    if citation.get('pdbx_database_id_PubMed'):
                        pubmed_id = citation['pdbx_database_id_PubMed']
                        pubmed_label = QLabel(f"PubMed ID: {pubmed_id}")
                        citation_layout.addRow("PubMed:", pubmed_label)
                    
                    if citation.get('pdbx_database_id_DOI'):
                        doi = citation['pdbx_database_id_DOI']
                        doi_label = QLabel(f"DOI: {doi}")
                        doi_label.setWordWrap(True)
                        citation_layout.addRow("DOI:", doi_label)
                    
                    pubs_layout.addWidget(citation_group)
            else:
                pubs_layout.addWidget(QLabel("No citation information available"))
        else:
            pubs_layout.addWidget(QLabel("No publication data available"))
        
        citations_tab_layout.addWidget(pubs_group)
        citations_tab_layout.addStretch()
        metadata_tabs.addTab(citations_tab, "Citations")
    
    def _create_structural_features_tab(self, metadata_tabs, metadata):
        """Create structural features tab."""
        features_tab = QWidget()
        features_tab_layout = QVBoxLayout(features_tab)
        
        # Structural Features
        features_group = QGroupBox("Structural Features and Annotations")
        features_layout = QVBoxLayout(features_group)
        
        features_layout.addWidget(QLabel("Comprehensive structural feature analysis includes:"))
        
        # Secondary Structure
        ss_group = QGroupBox("Secondary Structure")
        ss_layout = QVBoxLayout(ss_group)
        ss_layout.addWidget(QLabel("• Alpha helices, beta sheets, and loop regions"))
        ss_layout.addWidget(QLabel("• Secondary structure assignment (DSSP/STRIDE)"))
        ss_layout.addWidget(QLabel("• Ramachandran plot analysis"))
        features_layout.addWidget(ss_group)
        
        # Domains and Motifs
        domains_group = QGroupBox("Domains and Motifs")
        domains_layout = QVBoxLayout(domains_group)
        domains_layout.addWidget(QLabel("• Protein domain identification and boundaries"))
        domains_layout.addWidget(QLabel("• Functional motifs and active sites"))
        domains_layout.addWidget(QLabel("• Structural motifs and supersecondary structure"))
        features_layout.addWidget(domains_group)
        
        # Geometric Features
        geometry_group = QGroupBox("Geometric Features")
        geometry_layout = QVBoxLayout(geometry_group)
        geometry_layout.addWidget(QLabel("• Disulfide bonds and cross-links"))
        geometry_layout.addWidget(QLabel("• Hydrogen bonding patterns"))
        geometry_layout.addWidget(QLabel("• Surface accessibility and cavities"))
        features_layout.addWidget(geometry_group)
        
        features_tab_layout.addWidget(features_group)
        features_tab_layout.addStretch()
        metadata_tabs.addTab(features_tab, "Features")
    
    # Crystallographic tab creation method removed
    
    def display_basic_results_fallback(self, results):
        """Display basic results when enhanced metadata is not available."""
        basic_group = QGroupBox("Basic Properties")
        basic_layout = QFormLayout(basic_group)
        
        basic_layout.addRow("Structure ID:", QLabel(results.get('structure_id', 'N/A')))
        basic_layout.addRow("Total Chains:", QLabel(str(len(results.get('chains', [])))))
        basic_layout.addRow("Total Residues:", QLabel(str(results.get('total_residues', 0))))
        basic_layout.addRow("Total Atoms:", QLabel(str(results.get('total_atoms', 0))))
        basic_layout.addRow("Molecular Weight:", QLabel(f"{results.get('molecular_weight', 0):.1f} Da"))
        
        self.results_layout.addWidget(basic_group)
        
        # Also display composition analysis
        if results.get('composition'):
            self.display_composition_analysis(results)
    
    def display_composition_analysis(self, results):
        """Display amino acid composition analysis."""
        comp_group = QGroupBox("Amino Acid Composition")
        comp_layout = QVBoxLayout(comp_group)
        
        # Create horizontal layout for table and chart
        content_layout = QHBoxLayout()
        
        # Table
        table_widget = QWidget()
        table_layout = QVBoxLayout(table_widget)
        
        comp_table = QTableWidget()
        comp_table.setColumnCount(3)
        comp_table.setHorizontalHeaderLabels(["Residue", "Count", "Percentage"])
        
        composition = results['composition']
        total_residues = sum(composition.values())
        comp_table.setRowCount(len(composition))
        
        for i, (residue, count) in enumerate(sorted(composition.items())):
            percentage = (count / total_residues) * 100 if total_residues > 0 else 0
            comp_table.setItem(i, 0, QTableWidgetItem(residue))
            comp_table.setItem(i, 1, QTableWidgetItem(str(count)))
            comp_table.setItem(i, 2, QTableWidgetItem(f"{percentage:.1f}%"))
        
        comp_table.resizeColumnsToContents()
        comp_table.setMinimumHeight(200)
        comp_table.setMaximumWidth(400)
        comp_table.setMaximumHeight(600)
        table_layout.addWidget(comp_table)
        table_widget.setMaximumWidth(400)
        content_layout.addWidget(table_widget)
        
        # Pie chart with legend
        pie_chart = self.create_pie_chart(composition, "Amino Acid Distribution")
        if pie_chart:
            content_layout.addWidget(pie_chart, 2)  # Give chart much more space for legend
        
        comp_layout.addLayout(content_layout)
        self.results_layout.addWidget(comp_group)
    
    def display_results(self, results):
        """Display analysis results."""
        self.progress_bar.setVisible(False)
        
        # Store results for export
        self.current_results = results
        self.export_btn.setEnabled(True)
        
        # Clear existing results
        for i in reversed(range(self.results_layout.count())):
            item = self.results_layout.itemAt(i)
            if item is not None:
                widget = item.widget()
                if widget is not None:
                    widget.setParent(None)
                else:
                    self.results_layout.removeItem(item)
        
        # Display results for each analysis type
        if 'basic' in results:
            self.display_basic_results(results['basic'])
        
        if 'secondary' in results:
            self.display_secondary_results(results['secondary'])
        
        if 'geometry' in results:
            self.display_geometry_results(results['geometry'])
        

        

        
        self.results_layout.addStretch()
        
        if hasattr(self.parent_app, 'statusBar'):
            self.parent_app.statusBar().showMessage("Structural analysis complete")
    
    def display_basic_results(self, results):
        """Display basic properties results."""
        basic_group = QGroupBox("Basic Properties")
        basic_layout = QFormLayout(basic_group)
        
        # Add enhanced metadata if available
        if hasattr(self, 'current_comprehensive_data') and self.current_comprehensive_data:
            self.display_enhanced_metadata(results)
            return
        
        basic_layout.addRow("Structure ID:", QLabel(results.get('structure_id', 'N/A')))
        basic_layout.addRow("Total Chains:", QLabel(str(len(results.get('chains', [])))))
        basic_layout.addRow("Total Residues:", QLabel(str(results.get('total_residues', 0))))
        basic_layout.addRow("Total Atoms:", QLabel(str(results.get('total_atoms', 0))))
        basic_layout.addRow("Molecular Weight:", QLabel(f"{results.get('molecular_weight', 0):.1f} Da"))
        # Resolution and space group removed
        
        self.results_layout.addWidget(basic_group)
        
        # Chain details with comparison chart
        if results.get('chains'):
            chain_group = QGroupBox("Chain Analysis")
            chain_layout = QVBoxLayout(chain_group)
            
            # Chain comparison chart
            if len(results['chains']) > 1:
                chain_comparison_group = QGroupBox("Chain Size Comparison")
                chain_comparison_layout = QHBoxLayout(chain_comparison_group)
                
                # Create bar charts for chain comparison
                chain_residues = {f"Chain {chain['id']}": chain['residues'] for chain in results['chains']}
                chain_atoms = {f"Chain {chain['id']}": chain['atoms'] for chain in results['chains']}
                
                residue_chart = self.create_bar_chart(chain_residues, "Residues per Chain", "Chain", "Number of Residues")
                if residue_chart:
                    chain_comparison_layout.addWidget(residue_chart)
                
                atom_chart = self.create_bar_chart(chain_atoms, "Atoms per Chain", "Chain", "Number of Atoms")
                if atom_chart:
                    chain_comparison_layout.addWidget(atom_chart)
                
                chain_layout.addWidget(chain_comparison_group)
            
            # Chain details table
            table_group = QGroupBox("Chain Details")
            table_layout = QVBoxLayout(table_group)
            
            chain_table = QTableWidget()
            chain_table.setColumnCount(4)
            chain_table.setHorizontalHeaderLabels(["Chain ID", "Residues", "Atoms", "Sequence"])
            chain_table.setRowCount(len(results['chains']))
            
            for i, chain in enumerate(results['chains']):
                chain_table.setItem(i, 0, QTableWidgetItem(chain['id']))
                chain_table.setItem(i, 1, QTableWidgetItem(str(chain['residues'])))
                chain_table.setItem(i, 2, QTableWidgetItem(str(chain['atoms'])))
                
                # Show full sequence
                sequence = chain['sequence']
                chain_table.setItem(i, 3, QTableWidgetItem(sequence))
            
            chain_table.resizeColumnsToContents()
            chain_table.setMinimumHeight(100)
            chain_table.setMaximumHeight(400)
            
            # Set word wrap for sequence column and add tooltips
            for row in range(chain_table.rowCount()):
                item = chain_table.item(row, 3)  # Sequence column
                if item:
                    item.setToolTip(f"Full sequence ({len(item.text())} residues):\n{item.text()}")
            
            # Make sequence column wider and allow text wrapping
            chain_table.setColumnWidth(3, 400)
            chain_table.horizontalHeader().setStretchLastSection(True)
            
            # Enable horizontal scrolling for long sequences
            chain_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            chain_table.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            table_layout.addWidget(chain_table)
            
            chain_layout.addWidget(table_group)
            self.results_layout.addWidget(chain_group)
        
        # Composition
        if results.get('composition'):
            comp_group = QGroupBox("Amino Acid Composition")
            comp_layout = QVBoxLayout(comp_group)
            
            # Create horizontal layout for table and chart
            content_layout = QHBoxLayout()
            
            # Table
            table_widget = QWidget()
            table_layout = QVBoxLayout(table_widget)
            
            comp_table = QTableWidget()
            comp_table.setColumnCount(3)
            comp_table.setHorizontalHeaderLabels(["Residue", "Count", "Percentage"])
            
            composition = results['composition']
            total_residues = sum(composition.values())
            comp_table.setRowCount(len(composition))
            
            for i, (residue, count) in enumerate(sorted(composition.items())):
                percentage = (count / total_residues) * 100 if total_residues > 0 else 0
                comp_table.setItem(i, 0, QTableWidgetItem(residue))
                comp_table.setItem(i, 1, QTableWidgetItem(str(count)))
                comp_table.setItem(i, 2, QTableWidgetItem(f"{percentage:.1f}%"))
            
            comp_table.resizeColumnsToContents()
            comp_table.setMinimumHeight(200)
            comp_table.setMaximumWidth(400)
            comp_table.setMaximumHeight(600)
            table_layout.addWidget(comp_table)
            table_widget.setMaximumWidth(400)
            content_layout.addWidget(table_widget)
            
            # Pie chart with legend
            pie_chart = self.create_pie_chart(composition, "Amino Acid Distribution")
            if pie_chart:
                content_layout.addWidget(pie_chart, 2)  # Give chart much more space for legend
            
            comp_layout.addLayout(content_layout)
            self.results_layout.addWidget(comp_group)
    
    def display_secondary_results(self, results):
        """Display secondary structure results."""
        ss_group = QGroupBox("Secondary Structure Analysis")
        ss_layout = QVBoxLayout(ss_group)
        
        # Summary statistics and pie chart
        summary_layout = QHBoxLayout()
        
        # Statistics
        stats_widget = QWidget()
        stats_layout = QFormLayout(stats_widget)
        stats_layout.addRow("Helix Content:", QLabel(f"{results.get('helix_content', 0):.1f}%"))
        stats_layout.addRow("Sheet Content:", QLabel(f"{results.get('sheet_content', 0):.1f}%"))
        stats_layout.addRow("Loop Content:", QLabel(f"{results.get('loop_content', 0):.1f}%"))
        stats_widget.setMaximumWidth(300)
        stats_widget.setMinimumHeight(150)
        summary_layout.addWidget(stats_widget)
        
        # Pie chart for secondary structure distribution
        ss_data = {
            'Helix': results.get('helix_content', 0),
            'Sheet': results.get('sheet_content', 0),
            'Loop': results.get('loop_content', 0)
        }
        ss_colors = ['#ff9999', '#66b3ff', '#99ff99']  # Red, Blue, Green
        pie_chart = self.create_pie_chart(ss_data, "Secondary Structure Distribution", ss_colors)
        if pie_chart:
            summary_layout.addWidget(pie_chart, 1)  # Give chart more space
        
        ss_layout.addLayout(summary_layout)
        self.results_layout.addWidget(ss_group)
        
        # Ramachandran data
        if results.get('ramachandran_data'):
            rama_group = QGroupBox("Ramachandran Analysis")
            rama_layout = QVBoxLayout(rama_group)
            
            # Create horizontal layout for plot and table
            rama_content_layout = QHBoxLayout()
            
            # Ramachandran plot
            # Get PDB ID from results (passed from analysis worker)
            pdb_id = results.get('pdb_id')
            
            rama_plot, rama_counts = self.create_ramachandran_plot(results['ramachandran_data'], pdb_id)
            if rama_plot:
                rama_content_layout.addWidget(rama_plot, 2)  # Give plot more space
            
            # Table widget
            table_widget = QWidget()
            table_layout = QVBoxLayout(table_widget)
            
            rama_table = QTableWidget()
            rama_table.setColumnCount(5)
            rama_table.setHorizontalHeaderLabels(["Residue", "Chain", "Position", "Phi (°)", "Psi (°)"])
            
            # Show ALL Ramachandran data
            rama_data = results['ramachandran_data']
            rama_table.setRowCount(len(rama_data))
            
            for i, data in enumerate(rama_data):
                rama_table.setItem(i, 0, QTableWidgetItem(data['residue']))
                rama_table.setItem(i, 1, QTableWidgetItem(data['chain']))
                rama_table.setItem(i, 2, QTableWidgetItem(str(data['position'])))
                rama_table.setItem(i, 3, QTableWidgetItem(f"{data['phi']:.1f}"))
                rama_table.setItem(i, 4, QTableWidgetItem(f"{data['psi']:.1f}"))
            
            rama_table.resizeColumnsToContents()
            rama_table.setMinimumHeight(300)
            rama_table.setMaximumWidth(550)
            rama_table.setMaximumHeight(600)
            table_layout.addWidget(rama_table)
            
            # Show total count
            count_label = QLabel(f"<i>Total: {len(results['ramachandran_data'])} residues</i>")
            count_label.setStyleSheet("color: #666; font-style: italic;")
            table_layout.addWidget(count_label)
            
            table_widget.setMaximumWidth(550)
            rama_content_layout.addWidget(table_widget, 1)  # Table gets less space than plot
            
            rama_layout.addLayout(rama_content_layout)
            self.results_layout.addWidget(rama_group)
    
    def display_geometry_results(self, results):
        """Display geometric analysis results."""
        geom_group = QGroupBox("Geometric Properties")
        geom_layout = QFormLayout(geom_group)
        
        # Center of mass and geometric center
        com = results.get('center_of_mass', [0, 0, 0])
        geom_layout.addRow("Center of Mass:", QLabel(f"({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})"))
        
        gc = results.get('geometric_center', [0, 0, 0])
        geom_layout.addRow("Geometric Center:", QLabel(f"({gc[0]:.2f}, {gc[1]:.2f}, {gc[2]:.2f})"))
        
        geom_layout.addRow("Radius of Gyration:", QLabel(f"{results.get('radius_of_gyration', 0):.2f} Å"))
        
        self.results_layout.addWidget(geom_group)
        
        # Bond statistics with histograms
        if results.get('bond_lengths') or results.get('bond_angles'):
            bond_group = QGroupBox("Bond Geometry Analysis")
            bond_layout = QVBoxLayout(bond_group)
            
            # Bond lengths analysis
            if results.get('bond_lengths'):
                bond_lengths = results['bond_lengths']
                if bond_lengths:
                    bl_group = QGroupBox("Bond Lengths")
                    bl_layout = QHBoxLayout(bl_group)
                    
                    # Statistics
                    stats_widget = QWidget()
                    stats_layout = QFormLayout(stats_widget)
                    mean_length = np.mean(bond_lengths)
                    std_length = np.std(bond_lengths)
                    min_length = min(bond_lengths)
                    max_length = max(bond_lengths)
                    
                    stats_layout.addRow("Mean Length:", QLabel(f"{mean_length:.3f} Å"))
                    stats_layout.addRow("Std Dev:", QLabel(f"{std_length:.3f} Å"))
                    stats_layout.addRow("Range:", QLabel(f"{min_length:.3f} - {max_length:.3f} Å"))
                    stats_layout.addRow("Sample Size:", QLabel(str(len(bond_lengths))))
                    
                    stats_widget.setMaximumWidth(300)
                    bl_layout.addWidget(stats_widget)
                    
                    # Histogram
                    bl_histogram = self.create_histogram(bond_lengths, "Bond Length Distribution", 
                                                       "Bond Length (Å)", "Frequency", bins=25)
                    if bl_histogram:
                        bl_layout.addWidget(bl_histogram, 2)
                    
                    bond_layout.addWidget(bl_group)
            
            # Bond angles analysis
            if results.get('bond_angles'):
                bond_angles = results['bond_angles']
                if bond_angles:
                    ba_group = QGroupBox("Bond Angles")
                    ba_layout = QHBoxLayout(ba_group)
                    
                    # Statistics
                    stats_widget = QWidget()
                    stats_layout = QFormLayout(stats_widget)
                    mean_angle = np.mean(bond_angles)
                    std_angle = np.std(bond_angles)
                    min_angle = min(bond_angles)
                    max_angle = max(bond_angles)
                    
                    stats_layout.addRow("Mean Angle:", QLabel(f"{mean_angle:.1f}°"))
                    stats_layout.addRow("Std Dev:", QLabel(f"{std_angle:.1f}°"))
                    stats_layout.addRow("Range:", QLabel(f"{min_angle:.1f}° - {max_angle:.1f}°"))
                    stats_layout.addRow("Sample Size:", QLabel(str(len(bond_angles))))
                    
                    stats_widget.setMaximumWidth(300)
                    ba_layout.addWidget(stats_widget)
                    
                    # Histogram
                    ba_histogram = self.create_histogram(bond_angles, "Bond Angle Distribution", 
                                                       "Bond Angle (degrees)", "Frequency", bins=25)
                    if ba_histogram:
                        ba_layout.addWidget(ba_histogram, 2)
                    
                    bond_layout.addWidget(ba_group)
            
            self.results_layout.addWidget(bond_group)
        
        # B-factor statistics
        if results.get('b_factors'):
            b_factors = results['b_factors']
            if b_factors:
                bf_group = QGroupBox("B-factor Analysis")
                bf_layout = QVBoxLayout(bf_group)
                
                # Create horizontal layout for stats and histogram
                bf_content_layout = QHBoxLayout()
                
                # Statistics
                stats_widget = QWidget()
                stats_layout = QFormLayout(stats_widget)
                
                mean_bf = np.mean(b_factors)
                std_bf = np.std(b_factors)
                min_bf = min(b_factors)
                max_bf = max(b_factors)
                
                stats_layout.addRow("Mean B-factor:", QLabel(f"{mean_bf:.2f} Ų"))
                stats_layout.addRow("Std B-factor:", QLabel(f"{std_bf:.2f} Ų"))
                stats_layout.addRow("Range:", QLabel(f"{min_bf:.2f} - {max_bf:.2f} Ų"))
                stats_layout.addRow("Total Atoms:", QLabel(str(len(b_factors))))
                
                stats_widget.setMaximumWidth(350)
                stats_widget.setMinimumHeight(200)
                bf_content_layout.addWidget(stats_widget)
                
                # Histogram
                bf_histogram = self.create_histogram(b_factors, "B-factor Distribution", 
                                                   "B-factor (Ų)", "Frequency", bins=30)
                if bf_histogram:
                    bf_content_layout.addWidget(bf_histogram, 2)  # Give histogram more space
                
                bf_layout.addLayout(bf_content_layout)
                self.results_layout.addWidget(bf_group)
    

    
    def display_surface_results(self, results):
        """Display surface analysis results."""
        # Create title with timing information
        method_used = results.get('method_used', 'unknown')
        calc_time = results.get('calculation_time', 0)
        title = f"Surface Properties Analysis ({method_used} method, {calc_time:.2f}s)"
        
        surface_group = QGroupBox(title)
        surface_layout = QVBoxLayout(surface_group)
        
        # Add performance info
        if calc_time > 0:
            perf_label = QLabel(f"<i>Calculation completed in {calc_time:.2f} seconds using {method_used} method</i>")
            perf_label.setStyleSheet("color: #666; font-style: italic; margin-bottom: 10px;")
            surface_layout.addWidget(perf_label)
        
        # Surface area summary with pie chart
        area_summary_group = QGroupBox("Surface Area Distribution")
        area_summary_layout = QHBoxLayout(area_summary_group)
        
        # Statistics
        stats_widget = QWidget()
        stats_layout = QFormLayout(stats_widget)
        
        accessible_area = results.get('accessible_surface_area', 0)
        buried_area = results.get('buried_surface_area', 0)
        total_area = accessible_area + buried_area
        
        stats_layout.addRow("Accessible Surface Area:", QLabel(f"{accessible_area:.1f} Ų"))
        stats_layout.addRow("Buried Surface Area:", QLabel(f"{buried_area:.1f} Ų"))
        stats_layout.addRow("Total Surface Area:", QLabel(f"{total_area:.1f} Ų"))
        
        if total_area > 0:
            accessible_pct = (accessible_area / total_area) * 100
            buried_pct = (buried_area / total_area) * 100
            stats_layout.addRow("Accessible %:", QLabel(f"{accessible_pct:.1f}%"))
            stats_layout.addRow("Buried %:", QLabel(f"{buried_pct:.1f}%"))
        
        stats_widget.setMaximumWidth(350)
        area_summary_layout.addWidget(stats_widget)
        
        # Pie chart for surface area distribution
        if total_area > 0:
            area_data = {
                'Accessible': accessible_area,
                'Buried': buried_area
            }
            area_colors = ['#87CEEB', '#8B4513']  # Sky blue, Brown
            area_pie = self.create_pie_chart(area_data, "Surface Area Distribution", area_colors)
            if area_pie:
                area_summary_layout.addWidget(area_pie, 1)
        
        surface_layout.addWidget(area_summary_group)
        
        # Residue distribution with pie chart
        surface_residues = results.get('surface_residues', [])
        buried_residues = results.get('buried_residues', [])
        total_residues = len(surface_residues) + len(buried_residues)
        
        if total_residues > 0:
            residue_summary_group = QGroupBox("Residue Distribution")
            residue_summary_layout = QHBoxLayout(residue_summary_group)
            
            # Statistics
            res_stats_widget = QWidget()
            res_stats_layout = QFormLayout(res_stats_widget)
            
            surface_count = len(surface_residues)
            buried_count = len(buried_residues)
            surface_pct = (surface_count / total_residues) * 100
            buried_pct = (buried_count / total_residues) * 100
            
            res_stats_layout.addRow("Surface Residues:", QLabel(f"{surface_count} ({surface_pct:.1f}%)"))
            res_stats_layout.addRow("Buried Residues:", QLabel(f"{buried_count} ({buried_pct:.1f}%)"))
            res_stats_layout.addRow("Total Residues:", QLabel(str(total_residues)))
            
            res_stats_widget.setMaximumWidth(350)
            residue_summary_layout.addWidget(res_stats_widget)
            
            # Pie chart for residue distribution
            residue_data = {
                'Surface': surface_count,
                'Buried': buried_count
            }
            residue_colors = ['#90EE90', '#CD853F']  # Light green, Peru
            residue_pie = self.create_pie_chart(residue_data, "Surface vs Buried Residues", residue_colors)
            if residue_pie:
                residue_summary_layout.addWidget(residue_pie, 1)
            
            surface_layout.addWidget(residue_summary_group)
        
        self.results_layout.addWidget(surface_group)
        
        # Surface residues
        if results.get('surface_residues'):
            surf_res_group = QGroupBox("Surface Residues (Sample)")
            surf_res_layout = QVBoxLayout(surf_res_group)
            
            surf_table = QTableWidget()
            surf_table.setColumnCount(4)
            surf_table.setHorizontalHeaderLabels(["Residue", "Chain", "Position", "RSA (%)"])
            
            # Show ALL surface residues
            surface_residues = results['surface_residues']
            surf_table.setRowCount(len(surface_residues))
            
            for i, res in enumerate(surface_residues):
                surf_table.setItem(i, 0, QTableWidgetItem(res['residue']))
                surf_table.setItem(i, 1, QTableWidgetItem(res['chain']))
                surf_table.setItem(i, 2, QTableWidgetItem(str(res['position'])))
                surf_table.setItem(i, 3, QTableWidgetItem(f"{res['rsa']:.1f}"))
            
            surf_table.resizeColumnsToContents()
            # Allow table to expand to show all data
            surf_table.setMinimumHeight(100)
            # Remove maximum height restriction to show all data
            surf_res_layout.addWidget(surf_table)
            
            # Show total count
            count_label = QLabel(f"<i>Total: {len(results['surface_residues'])} surface residues</i>")
            count_label.setStyleSheet("color: #666; font-style: italic;")
            surf_res_layout.addWidget(count_label)
            
            self.results_layout.addWidget(surf_res_group)
    
    def display_cavity_results(self, results):
        """Display cavity detection results."""
        # Create title with timing information
        calc_time = results.get('calculation_time', 0)
        cavity_points = results.get('cavity_points_found', 0)
        title = f"Binding Sites and Cavities Analysis ({calc_time:.1f}s)"
        
        cavity_group = QGroupBox(title)
        cavity_layout = QVBoxLayout(cavity_group)
        
        # Add performance and detection info
        if calc_time > 0 or cavity_points > 0:
            perf_label = QLabel(
                f"<i>Detection completed in {calc_time:.1f} seconds. "
                f"Scanned grid and found {cavity_points:,} potential cavity points.</i>"
            )
            perf_label.setStyleSheet("color: #666; font-style: italic; margin-bottom: 10px;")
            cavity_layout.addWidget(perf_label)
        
        # Cavity summary with statistics
        summary_group = QGroupBox("Cavity Summary")
        summary_layout = QHBoxLayout(summary_group)
        
        # Statistics
        stats_widget = QWidget()
        stats_layout = QFormLayout(stats_widget)
        
        cavities = results.get('potential_cavities', [])
        total_volume = results.get('cavity_volume', 0)
        
        stats_layout.addRow("Potential Cavities:", QLabel(str(len(cavities))))
        stats_layout.addRow("Total Cavity Volume:", QLabel(f"{total_volume:.1f} ų"))
        
        if results.get('largest_cavity'):
            largest = results['largest_cavity']
            stats_layout.addRow("Largest Cavity Volume:", QLabel(f"{largest.get('volume', 0):.1f} ų"))
            stats_layout.addRow("Largest Cavity Radius:", QLabel(f"{largest.get('radius', 0):.1f} Å"))
        
        if cavities:
            avg_volume = total_volume / len(cavities)
            stats_layout.addRow("Average Cavity Volume:", QLabel(f"{avg_volume:.1f} ų"))
        
        stats_widget.setMaximumWidth(350)
        summary_layout.addWidget(stats_widget)
        
        # Bar chart for cavity volumes
        if cavities and len(cavities) > 1:
            cavity_volumes = {f"Cavity {i+1}": cavity.get('volume', 0) for i, cavity in enumerate(cavities[:10])}  # Show top 10
            volume_chart = self.create_bar_chart(cavity_volumes, "Cavity Volumes", "Cavity", "Volume (ų)")
            if volume_chart:
                summary_layout.addWidget(volume_chart, 2)
        
        cavity_layout.addWidget(summary_group)
        self.results_layout.addWidget(cavity_group)
        
        # Cavity details
        if results.get('potential_cavities'):
            cav_detail_group = QGroupBox("Cavity Details")
            cav_detail_layout = QVBoxLayout(cav_detail_group)
            
            cav_table = QTableWidget()
            cav_table.setColumnCount(5)
            cav_table.setHorizontalHeaderLabels(["Cavity #", "Volume (ų)", "Radius (Å)", "Points", "Center"])
            
            cavities = results['potential_cavities']
            cav_table.setRowCount(len(cavities))
            
            for i, cavity in enumerate(cavities):
                center = cavity.get('center', [0, 0, 0])
                center_str = f"({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})"
                
                cav_table.setItem(i, 0, QTableWidgetItem(str(i + 1)))
                cav_table.setItem(i, 1, QTableWidgetItem(f"{cavity.get('volume', 0):.1f}"))
                cav_table.setItem(i, 2, QTableWidgetItem(f"{cavity.get('radius', 0):.1f}"))
                cav_table.setItem(i, 3, QTableWidgetItem(str(cavity.get('point_count', 0))))
                cav_table.setItem(i, 4, QTableWidgetItem(center_str))
            
            cav_table.resizeColumnsToContents()
            # Allow table to expand to show all cavities
            cav_table.setMinimumHeight(100)
            # Remove maximum height restriction to show all data
            cav_detail_layout.addWidget(cav_table)
            
            self.results_layout.addWidget(cav_detail_group)
    
    def export_results(self):
        """Export analysis results with all visualizations."""
        if not self.current_results:
            QMessageBox.warning(self, "Export Error", "No analysis results to export. Please run an analysis first.")
            return
        
        # Simple file dialog for HTML export
        structure_id = self.current_results.get('structure_id', 'analysis')
        suggested_name = f"{structure_id}_complete_analysis.html"
        
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Complete Analysis Report",
            suggested_name,
            "HTML Files (*.html);;PDF Files (*.pdf)"
        )
        
        if file_path:
            try:
                if file_path.endswith('.pdf'):
                    self.export_complete_pdf(file_path)
                else:
                    self.export_complete_html(file_path)
                
                QMessageBox.information(self, "Export Complete", f"Complete analysis report exported to:\n{file_path}")
                
                if hasattr(self.parent_app, 'statusBar'):
                    self.parent_app.statusBar().showMessage(f"Complete report exported successfully")
                    
            except Exception as e:
                QMessageBox.critical(self, "Export Error", f"Failed to export complete report:\n{str(e)}")
    
    def export_to_json(self, file_path, options):
        """Export results to JSON format."""
        import json
        from datetime import datetime
        
        # Prepare data for JSON serialization
        export_data = {
            'export_info': {
                'timestamp': datetime.now().isoformat(),
                'picomol_version': '1.0',
                'export_format': 'JSON'
            },
            'structure_info': {
                'structure_id': self.current_results.get('structure_id', 'Unknown'),
                'structure_path': self.current_results.get('structure_path', 'Unknown')
            },
            'analysis_results': {}
        }
        
        # Add comprehensive metadata if available
        if hasattr(self, 'current_comprehensive_data') and self.current_comprehensive_data:
            export_data['comprehensive_metadata'] = self._extract_comprehensive_metadata_for_export()
        
        # Experimental details removed
        
        # Add selected analysis types
        for analysis_type in ['basic', 'secondary', 'geometry']:
            if analysis_type in self.current_results and options.get(f'include_{analysis_type}', True):
                export_data['analysis_results'][analysis_type] = self.current_results[analysis_type]
        
        # Write JSON file
        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(export_data, f, indent=2, default=str)
    
    def export_to_csv(self, file_path, options):
        """Export results to CSV format."""
        import csv
        from datetime import datetime
        
        with open(file_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Header information
            writer.writerow(['PicoMol Structural Analysis Results'])
            writer.writerow(['Export Date:', datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
            writer.writerow(['Structure ID:', self.current_results.get('structure_id', 'Unknown')])
            writer.writerow(['Structure Path:', self.current_results.get('structure_path', 'Unknown')])
            writer.writerow([])
            
            # Add comprehensive metadata if available
            if hasattr(self, 'current_comprehensive_data') and self.current_comprehensive_data:
                comp_metadata = self._extract_comprehensive_metadata_for_export()
                if comp_metadata:
                    writer.writerow(['=== COMPREHENSIVE METADATA ==='])
                    for key, value in comp_metadata.items():
                        if isinstance(value, list):
                            value = ', '.join(str(v) for v in value)
                        writer.writerow([key.replace('_', ' ').title() + ':', str(value)])
                    writer.writerow([])
            
            # Experimental details removed
            
            # Basic Properties
            if 'basic' in self.current_results and options.get('include_basic', True):
                basic = self.current_results['basic']
                writer.writerow(['=== BASIC PROPERTIES ==='])
                writer.writerow(['Property', 'Value'])
                writer.writerow(['Total Residues', basic.get('total_residues', 0)])
                writer.writerow(['Total Atoms', basic.get('total_atoms', 0)])
                writer.writerow(['Molecular Weight (Da)', f"{basic.get('molecular_weight', 0):.2f}"])
                # Resolution and space group removed
                writer.writerow([])
                
                # Chain information
                if 'chains' in basic:
                    writer.writerow(['=== CHAIN INFORMATION ==='])
                    writer.writerow(['Chain ID', 'Residues', 'Atoms', 'Sequence Length'])
                    for chain in basic['chains']:
                        writer.writerow([chain['id'], chain['residues'], chain['atoms'], len(chain.get('sequence', ''))])
                    writer.writerow([])
                
                # Composition
                if 'composition' in basic:
                    writer.writerow(['=== AMINO ACID COMPOSITION ==='])
                    writer.writerow(['Residue', 'Count', 'Percentage'])
                    total_residues = sum(basic['composition'].values())
                    for residue, count in sorted(basic['composition'].items()):
                        percentage = (count / total_residues) * 100 if total_residues > 0 else 0
                        writer.writerow([residue, count, f"{percentage:.2f}%"])
                    writer.writerow([])
            
            # Secondary Structure
            if 'secondary' in self.current_results and options.get('include_secondary', True):
                secondary = self.current_results['secondary']
                writer.writerow(['=== SECONDARY STRUCTURE ==='])
                writer.writerow(['Structure Type', 'Percentage'])
                writer.writerow(['Alpha Helix', f"{secondary.get('helix_content', 0):.2f}%"])
                writer.writerow(['Beta Sheet', f"{secondary.get('sheet_content', 0):.2f}%"])
                writer.writerow(['Loop/Coil', f"{secondary.get('loop_content', 0):.2f}%"])
                writer.writerow([])
                
                # Ramachandran data
                if 'ramachandran_data' in secondary and options.get('include_ramachandran', True):
                    writer.writerow(['=== RAMACHANDRAN DATA ==='])
                    writer.writerow(['Residue', 'Chain', 'Position', 'Phi (degrees)', 'Psi (degrees)'])
                    for data in secondary['ramachandran_data'][:100]:  # Limit to first 100
                        writer.writerow([data['residue'], data['chain'], data['position'], 
                                       f"{data['phi']:.2f}", f"{data['psi']:.2f}"])
                    writer.writerow([])
            
            # Geometry
            if 'geometry' in self.current_results and options.get('include_geometry', True):
                geometry = self.current_results['geometry']
                writer.writerow(['=== GEOMETRIC PROPERTIES ==='])
                writer.writerow(['Property', 'Value'])
                writer.writerow(['Radius of Gyration (Å)', f"{geometry.get('radius_of_gyration', 0):.2f}"])
                
                center = geometry.get('geometric_center', [0, 0, 0])
                writer.writerow(['Geometric Center (x, y, z)', f"({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})"])
                writer.writerow([])
            

            
            # Surface Properties

            if False:  # Surface export disabled
                surface = self.current_results['surface']
                writer.writerow(['=== SURFACE PROPERTIES ==='])
                writer.writerow(['Property', 'Value'])
                writer.writerow(['Total SASA (Ų)', f"{surface.get('total_sasa', 0):.2f}"])
                writer.writerow(['Accessible Surface Area (Ų)', f"{surface.get('accessible_surface_area', 0):.2f}"])
                writer.writerow(['Buried Surface Area (Ų)', f"{surface.get('buried_surface_area', 0):.2f}"])
                writer.writerow(['Surface Percentage', f"{surface.get('surface_percentage', 0):.2f}%"])
                writer.writerow(['Buried Percentage', f"{surface.get('buried_percentage', 0):.2f}%"])
                writer.writerow(['Average RSA', f"{surface.get('average_rsa', 0):.2f}%"])
                writer.writerow([])
    
    def export_to_excel(self, file_path, options):
        """Export results to Excel format."""
        try:
            import pandas as pd
            from datetime import datetime
            
            with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
                # Summary sheet
                summary_data = {
                    'Property': ['Export Date', 'Structure ID', 'Structure Path'],
                    'Value': [
                        datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                        self.current_results.get('structure_id', 'Unknown'),
                        self.current_results.get('structure_path', 'Unknown')
                    ]
                }
                pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)
                
                # Basic Properties
                if 'basic' in self.current_results and options.get('include_basic', True):
                    basic = self.current_results['basic']
                    
                    # Basic info
                    basic_data = {
                        'Property': ['Total Residues', 'Total Atoms', 'Molecular Weight (Da)'],
                        'Value': [
                            basic.get('total_residues', 0),
                            basic.get('total_atoms', 0),
                            f"{basic.get('molecular_weight', 0):.2f}"
                        ]
                    }
                    pd.DataFrame(basic_data).to_excel(writer, sheet_name='Basic Properties', index=False)
                    
                    # Chain information
                    if 'chains' in basic:
                        chain_data = []
                        for chain in basic['chains']:
                            chain_data.append({
                                'Chain ID': chain['id'],
                                'Residues': chain['residues'],
                                'Atoms': chain['atoms'],
                                'Sequence Length': len(chain.get('sequence', ''))
                            })
                        if chain_data:
                            pd.DataFrame(chain_data).to_excel(writer, sheet_name='Chains', index=False)
                    
                    # Composition
                    if 'composition' in basic:
                        comp_data = []
                        total_residues = sum(basic['composition'].values())
                        for residue, count in sorted(basic['composition'].items()):
                            percentage = (count / total_residues) * 100 if total_residues > 0 else 0
                            comp_data.append({
                                'Residue': residue,
                                'Count': count,
                                'Percentage': f"{percentage:.2f}%"
                            })
                        if comp_data:
                            pd.DataFrame(comp_data).to_excel(writer, sheet_name='Composition', index=False)
                
                # Secondary Structure
                if 'secondary' in self.current_results and options.get('include_secondary', True):
                    secondary = self.current_results['secondary']
                    
                    ss_data = {
                        'Structure Type': ['Alpha Helix', 'Beta Sheet', 'Loop/Coil'],
                        'Percentage': [
                            f"{secondary.get('helix_content', 0):.2f}%",
                            f"{secondary.get('sheet_content', 0):.2f}%",
                            f"{secondary.get('loop_content', 0):.2f}%"
                        ]
                    }
                    pd.DataFrame(ss_data).to_excel(writer, sheet_name='Secondary Structure', index=False)
                    
                    # Ramachandran data
                    if 'ramachandran_data' in secondary and options.get('include_ramachandran', True):
                        rama_data = []
                        for data in secondary['ramachandran_data'][:500]:  # Limit to first 500
                            rama_data.append({
                                'Residue': data['residue'],
                                'Chain': data['chain'],
                                'Position': data['position'],
                                'Phi (degrees)': f"{data['phi']:.2f}",
                                'Psi (degrees)': f"{data['psi']:.2f}"
                            })
                        if rama_data:
                            pd.DataFrame(rama_data).to_excel(writer, sheet_name='Ramachandran', index=False)
                

                
                # Surface Properties

                if False:  # Surface export disabled
                    surface = self.current_results['surface']
                    
                    surface_data = {
                        'Property': ['Total SASA (Ų)', 'Accessible Surface Area (Ų)', 'Buried Surface Area (Ų)', 
                                   'Surface Percentage', 'Buried Percentage', 'Average RSA'],
                        'Value': [
                            f"{surface.get('total_sasa', 0):.2f}",
                            f"{surface.get('accessible_surface_area', 0):.2f}",
                            f"{surface.get('buried_surface_area', 0):.2f}",
                            f"{surface.get('surface_percentage', 0):.2f}%",
                            f"{surface.get('buried_percentage', 0):.2f}%",
                            f"{surface.get('average_rsa', 0):.2f}%"
                        ]
                    }
                    pd.DataFrame(surface_data).to_excel(writer, sheet_name='Surface', index=False)
                    
                    # Surface residues
                    if 'surface_residues' in surface and options.get('include_residue_details', True):
                        surf_res_data = []
                        for res in surface['surface_residues'][:1000]:  # Limit to first 1000
                            surf_res_data.append({
                                'Residue': res['residue'],
                                'Chain': res['chain'],
                                'Position': res['position'],
                                'SASA (Ų)': f"{res['sasa']:.2f}",
                                'RSA (%)': f"{res['rsa']:.2f}"
                            })
                        if surf_res_data:
                            pd.DataFrame(surf_res_data).to_excel(writer, sheet_name='Surface Residues', index=False)
        
        except ImportError:
            # Fallback to CSV if pandas/openpyxl not available
            csv_path = file_path.replace('.xlsx', '.csv')
            self.export_to_csv(csv_path, options)
            raise ImportError(f"pandas/openpyxl not available. Exported to CSV instead: {csv_path}")
    
    def export_to_text(self, file_path, options):
        """Export results to text report format."""
        from datetime import datetime
        
        with open(file_path, 'w', encoding='utf-8') as f:
            # Header
            f.write("="*80 + "\n")
            f.write("PicoMol Structural Analysis Report\n")
            f.write("="*80 + "\n")
            f.write(f"Export Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Structure ID: {self.current_results.get('structure_id', 'Unknown')}\n")
            f.write(f"Structure Path: {self.current_results.get('structure_path', 'Unknown')}\n")
            f.write("\n")
            
            # Add comprehensive metadata if available
            if hasattr(self, 'current_comprehensive_data') and self.current_comprehensive_data:
                comp_metadata = self._extract_comprehensive_metadata_for_export()
                if comp_metadata:
                    f.write("=== COMPREHENSIVE METADATA ===\n")
                    for key, value in comp_metadata.items():
                        if isinstance(value, list):
                            value = ', '.join(str(v) for v in value)
                        f.write(f"{key.replace('_', ' ').title()}: {value}\n")
                    f.write("\n")
            
            # Experimental details removed
            
            # Basic Properties
            if 'basic' in self.current_results and options.get('include_basic', True):
                basic = self.current_results['basic']
                f.write("BASIC PROPERTIES\n")
                f.write("-"*40 + "\n")
                f.write(f"Total Residues: {basic.get('total_residues', 0)}\n")
                f.write(f"Total Atoms: {basic.get('total_atoms', 0)}\n")
                f.write(f"Molecular Weight: {basic.get('molecular_weight', 0):.2f} Da\n")
                # Resolution and space group removed
                f.write("\n")
                
                # Chain information
                if 'chains' in basic:
                    f.write("CHAIN INFORMATION\n")
                    f.write("-"*40 + "\n")
                    for chain in basic['chains']:
                        f.write(f"Chain {chain['id']}: {chain['residues']} residues, {chain['atoms']} atoms\n")
                    f.write("\n")
                
                # Composition
                if 'composition' in basic:
                    f.write("AMINO ACID COMPOSITION\n")
                    f.write("-"*40 + "\n")
                    total_residues = sum(basic['composition'].values())
                    for residue, count in sorted(basic['composition'].items()):
                        percentage = (count / total_residues) * 100 if total_residues > 0 else 0
                        f.write(f"{residue}: {count} ({percentage:.2f}%)\n")
                    f.write("\n")
            
            # Secondary Structure
            if 'secondary' in self.current_results and options.get('include_secondary', True):
                secondary = self.current_results['secondary']
                f.write("SECONDARY STRUCTURE\n")
                f.write("-"*40 + "\n")
                f.write(f"Alpha Helix: {secondary.get('helix_content', 0):.2f}%\n")
                f.write(f"Beta Sheet: {secondary.get('sheet_content', 0):.2f}%\n")
                f.write(f"Loop/Coil: {secondary.get('loop_content', 0):.2f}%\n")
                f.write("\n")
            
            # Geometry
            if 'geometry' in self.current_results and options.get('include_geometry', True):
                geometry = self.current_results['geometry']
                f.write("GEOMETRIC PROPERTIES\n")
                f.write("-"*40 + "\n")
                f.write(f"Radius of Gyration: {geometry.get('radius_of_gyration', 0):.2f} Å\n")
                center = geometry.get('geometric_center', [0, 0, 0])
                f.write(f"Geometric Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})\n")
                f.write("\n")
            

            
            # Surface Properties
    
            if False:  # Surface export disabled
                surface = self.current_results['surface']
                f.write("SURFACE PROPERTIES\n")
                f.write("-"*40 + "\n")
                f.write(f"Total SASA: {surface.get('total_sasa', 0):.2f} Ų\n")
                f.write(f"Accessible Surface Area: {surface.get('accessible_surface_area', 0):.2f} Ų\n")
                f.write(f"Buried Surface Area: {surface.get('buried_surface_area', 0):.2f} Ų\n")
                f.write(f"Surface Percentage: {surface.get('surface_percentage', 0):.2f}%\n")
                f.write(f"Buried Percentage: {surface.get('buried_percentage', 0):.2f}%\n")
                f.write(f"Average RSA: {surface.get('average_rsa', 0):.2f}%\n")
                f.write("\n")
    
    def export_complete_html(self, file_path):
        """Export complete analysis with all data and visualizations as HTML."""
        import base64
        import io
        from datetime import datetime
        
        # Create a temporary directory for images
        import tempfile
        temp_dir = tempfile.mkdtemp()
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                # Simple, clean HTML
                f.write("""<!DOCTYPE html>
<html>
<head>
    <title>PicoMol Structural Analysis Report</title>
    <meta charset="UTF-8">
    <style>
        body { 
            font-family: Arial, sans-serif; 
            margin: 40px; 
            line-height: 1.4;
        }
        h1 { 
            color: #333; 
            border-bottom: 2px solid #666; 
            padding-bottom: 10px;
        }
        h2 { 
            color: #444; 
            border-bottom: 1px solid #999; 
            padding-bottom: 5px;
            margin-top: 40px;
        }
        h3 { 
            color: #555; 
            margin-top: 25px;
        }
        table { 
            border-collapse: collapse; 
            width: 100%; 
            margin: 15px 0; 
        }
        th, td { 
            border: 1px solid #ccc; 
            padding: 8px; 
            text-align: left; 
        }
        th { 
            background-color: #f0f0f0; 
            font-weight: bold;
        }
        .chart { 
            text-align: center; 
            margin: 20px 0; 
        }
        .chart img { 
            max-width: 800px; 
            height: auto; 
        }
        .summary { 
            background-color: #f5f5f5; 
            padding: 15px; 
            margin-bottom: 20px;
            border: 1px solid #ddd;
        }
    </style>
</head>
<body>
""")
                
                # Title and summary
                f.write("<h1>PicoMol Structural Analysis Report</h1>\n")
                f.write("<div class='summary'>\n")
                f.write(f"<p><strong>Export Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n")
                f.write(f"<p><strong>Structure ID:</strong> {self.current_results.get('structure_id', 'Unknown')}</p>\n")
                f.write(f"<p><strong>Structure Path:</strong> {self.current_results.get('structure_path', 'Unknown')}</p>\n")
                
                # Add analysis summary
                analysis_types = []
                if 'basic' in self.current_results:
                    analysis_types.append('Basic Properties')
                if 'secondary' in self.current_results:
                    analysis_types.append('Secondary Structure')
                if 'geometry' in self.current_results:
                    analysis_types.append('Geometric Properties')
                if 'surface' in self.current_results:
                    analysis_types.append('Surface Analysis (SASA)')
                if 'cavity' in self.current_results:
                    analysis_types.append('Cavity Analysis')
                
                if analysis_types:
                    f.write(f"<p><strong>Analysis Types:</strong> {', '.join(analysis_types)}</p>\n")
                
                # Add chart summary
                chart_count = 0
                chart_types = []
                
                # Count charts that will be included
                if 'basic' in self.current_results and self.current_results['basic'].get('composition'):
                    chart_types.append('Amino Acid Composition')
                    chart_count += 1
                
                if 'basic' in self.current_results and len(self.current_results['basic'].get('chains', [])) > 1:
                    chart_types.extend(['Chain Comparison (Residues)', 'Chain Comparison (Atoms)'])
                    chart_count += 2
                
                if 'secondary' in self.current_results:
                    chart_types.append('Secondary Structure Distribution')
                    chart_count += 1
                    if self.current_results['secondary'].get('ramachandran_data'):
                        chart_types.append('Ramachandran Plot')
                        chart_count += 1
                
                if 'geometry' in self.current_results:
                    if self.current_results['geometry'].get('b_factors'):
                        chart_types.append('B-factor Distribution')
                        chart_count += 1
                    if self.current_results['geometry'].get('bond_lengths'):
                        chart_types.append('Bond Length Distribution')
                        chart_count += 1
                    if self.current_results['geometry'].get('bond_angles'):
                        chart_types.append('Bond Angle Distribution')
                        chart_count += 1
                
                if 'surface' in self.current_results:
                    chart_types.extend(['Surface Area Distribution', 'Surface vs Buried Residues'])
                    chart_count += 2
                
                if 'cavity' in self.current_results and len(self.current_results['cavity'].get('potential_cavities', [])) > 1:
                    chart_types.append('Cavity Volumes')
                    chart_count += 1
                
                if chart_count > 0:
                    f.write(f"<p><strong>Charts Included:</strong> {chart_count} charts ({', '.join(chart_types)})</p>\n")
                
                f.write("</div>\n\n")
                
                # Add comprehensive metadata if available
                if hasattr(self, 'current_comprehensive_data') and self.current_comprehensive_data:
                    self._export_comprehensive_metadata(f)
                
                # Add experimental details with enhanced fallback logic
                if 'basic' in self.current_results:
                    self._export_experimental_details(f, self.current_results['basic'])
                
                # Enhanced Basic Properties
                if 'basic' in self.current_results:
                    basic = self.current_results['basic']
                    f.write("<h2>Structure Properties</h2>\n")
                    
                    # Detailed table with enhanced information
                    f.write("<table>\n")
                    f.write("<tr><th>Property</th><th>Value</th></tr>\n")
                    f.write(f"<tr><td>Total Chains</td><td>{len(basic.get('chains', []))}</td></tr>\n")
                    f.write(f"<tr><td>Total Residues</td><td>{basic.get('total_residues', 0)}</td></tr>\n")
                    f.write(f"<tr><td>Total Atoms</td><td>{basic.get('total_atoms', 0)}</td></tr>\n")
                    f.write(f"<tr><td>Molecular Weight</td><td>{basic.get('molecular_weight', 0):.2f} Da</td></tr>\n")
                    
                    # Resolution and space group removed
                    f.write("</table>\n")
                    
                    # Enhanced Chain information with sequences
                    if 'chains' in basic:
                        f.write("<h3>Chain Information</h3>\n")
                        
                        # Add chain comparison charts if multiple chains
                        if len(basic['chains']) > 1 and MATPLOTLIB_AVAILABLE:
                            chain_residues = {f"Chain {chain['id']}": chain['residues'] for chain in basic['chains']}
                            chain_atoms = {f"Chain {chain['id']}": chain['atoms'] for chain in basic['chains']}
                            
                            # Chain residues chart
                            residue_chart_path = self.save_chart_as_image(chain_residues, "Residues per Chain", "bar", temp_dir, "Chain", "Number of Residues")
                            if residue_chart_path:
                                f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(residue_chart_path)}' alt='Residues per Chain Chart'></div>\n")
                            
                            # Chain atoms chart
                            atom_chart_path = self.save_chart_as_image(chain_atoms, "Atoms per Chain", "bar", temp_dir, "Chain", "Number of Atoms")
                            if atom_chart_path:
                                f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(atom_chart_path)}' alt='Atoms per Chain Chart'></div>\n")
                        
                        f.write("<table>\n")
                        f.write("<tr><th>Chain ID</th><th>Residues</th><th>Atoms</th><th>Sequence</th></tr>\n")
                        for chain in basic['chains']:
                            sequence = chain.get('sequence', '')
                            # Format sequence in blocks of 50 for better readability
                            formatted_sequence = '<br>'.join([sequence[i:i+50] for i in range(0, len(sequence), 50)]) if sequence else 'N/A'
                            f.write(f"<tr><td>{chain['id']}</td><td>{chain['residues']}</td><td>{chain['atoms']}</td><td style='font-family: monospace; font-size: 10px;'>{formatted_sequence}</td></tr>\n")
                        f.write("</table>\n")
                    
                    # Amino acid composition with chart
                    if 'composition' in basic:
                        f.write("<h3>Amino Acid Composition</h3>\n")
                        
                        # Generate composition chart
                        composition_data = basic['composition']
                        if composition_data and MATPLOTLIB_AVAILABLE:
                            chart_path = self.save_chart_as_image(composition_data, "Amino Acid Composition", "pie", temp_dir)
                            if chart_path:
                                f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(chart_path)}' alt='Amino Acid Composition Chart'></div>\n")
                        
                        # Composition table
                        f.write("<table>\n")
                        f.write("<tr><th>Residue</th><th>Count</th><th>Percentage</th></tr>\n")
                        total_residues = sum(composition_data.values())
                        for residue, count in sorted(composition_data.items()):
                            percentage = (count / total_residues) * 100 if total_residues > 0 else 0
                            f.write(f"<tr><td>{residue}</td><td>{count}</td><td>{percentage:.2f}%</td></tr>\n")
                        f.write("</table>\n")
                    

                
                # Secondary Structure
                if 'secondary' in self.current_results:
                    secondary = self.current_results['secondary']
                    f.write("<h2>Secondary Structure Analysis</h2>\n")
                    
                    # Secondary structure table
                    f.write("<table>\n")
                    f.write("<tr><th>Structure Type</th><th>Percentage</th></tr>\n")
                    f.write(f"<tr><td>Alpha Helix</td><td>{secondary.get('helix_content', 0):.1f}%</td></tr>\n")
                    f.write(f"<tr><td>Beta Sheet</td><td>{secondary.get('sheet_content', 0):.1f}%</td></tr>\n")
                    f.write(f"<tr><td>Loop/Coil</td><td>{secondary.get('loop_content', 0):.1f}%</td></tr>\n")
                    f.write("</table>\n")
                    
                    # Secondary structure chart
                    ss_data = {
                        'Alpha Helix': secondary.get('helix_content', 0),
                        'Beta Sheet': secondary.get('sheet_content', 0),
                        'Loop/Coil': secondary.get('loop_content', 0)
                    }
                    if MATPLOTLIB_AVAILABLE:
                        chart_path = self.save_chart_as_image(ss_data, "Secondary Structure Distribution", "pie", temp_dir)
                        if chart_path:
                            f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(chart_path)}' alt='Secondary Structure Chart'></div>\n")
                    
                    # Ramachandran plot
                    if 'ramachandran_data' in secondary and MATPLOTLIB_AVAILABLE:
                        # Get PDB ID for Ramachandran plot
                        pdb_id = None
                        if hasattr(self, 'current_comprehensive_data') and self.current_comprehensive_data:
                            pdb_id = self.current_comprehensive_data.get('pdb_id')
                        
                        rama_chart_path = self.save_ramachandran_plot(secondary['ramachandran_data'], temp_dir, pdb_id)
                        if rama_chart_path:
                            f.write("<h3>Ramachandran Plot</h3>\n")
                            f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(rama_chart_path)}' alt='Ramachandran Plot'></div>\n")
                
                # Geometry
                if 'geometry' in self.current_results:
                    geometry = self.current_results['geometry']
                    f.write("<h2>Geometric Properties</h2>\n")
                    
                    f.write("<table>\n")
                    f.write("<tr><th>Property</th><th>Value</th></tr>\n")
                    f.write(f"<tr><td>Radius of Gyration</td><td>{geometry.get('radius_of_gyration', 0):.2f} Å</td></tr>\n")
                    center = geometry.get('geometric_center', [0, 0, 0])
                    f.write(f"<tr><td>Geometric Center</td><td>({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})</td></tr>\n")
                    f.write("</table>\n")
                    
                    # B-factor histogram
                    if 'b_factors' in geometry and MATPLOTLIB_AVAILABLE:
                        bfactor_chart_path = self.save_histogram(geometry['b_factors'], "B-factor Distribution", "B-factor", "Frequency", temp_dir)
                        if bfactor_chart_path:
                            f.write("<h3>B-factor Distribution</h3>\n")
                            f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(bfactor_chart_path)}' alt='B-factor Distribution'></div>\n")
                    
                    # Bond length histogram
                    if 'bond_lengths' in geometry and geometry['bond_lengths'] and MATPLOTLIB_AVAILABLE:
                        bond_length_chart_path = self.save_histogram(geometry['bond_lengths'], "Bond Length Distribution", "Bond Length (Å)", "Frequency", temp_dir)
                        if bond_length_chart_path:
                            f.write("<h3>Bond Length Distribution</h3>\n")
                            f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(bond_length_chart_path)}' alt='Bond Length Distribution'></div>\n")
                    
                    # Bond angle histogram
                    if 'bond_angles' in geometry and geometry['bond_angles'] and MATPLOTLIB_AVAILABLE:
                        bond_angle_chart_path = self.save_histogram(geometry['bond_angles'], "Bond Angle Distribution", "Bond Angle (°)", "Frequency", temp_dir)
                        if bond_angle_chart_path:
                            f.write("<h3>Bond Angle Distribution</h3>\n")
                            f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(bond_angle_chart_path)}' alt='Bond Angle Distribution'></div>\n")
                

                    

                
                # Surface Properties
                if 'surface' in self.current_results:
                    surface = self.current_results['surface']
                    f.write("<h2>Surface Properties (SASA Analysis)</h2>\n")
                    
                    # Surface area distribution chart
                    accessible_area = surface.get('accessible_surface_area', 0)
                    buried_area = surface.get('buried_surface_area', 0)
                    if accessible_area > 0 or buried_area > 0:
                        area_data = {
                            'Accessible': accessible_area,
                            'Buried': buried_area
                        }
                        if MATPLOTLIB_AVAILABLE:
                            area_chart_path = self.save_chart_as_image(area_data, "Surface Area Distribution", "pie", temp_dir)
                            if area_chart_path:
                                f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(area_chart_path)}' alt='Surface Area Distribution Chart'></div>\n")
                    
                    # Surface/buried residue distribution chart
                    surface_residues = len(surface.get('surface_residues', []))
                    buried_residues = len(surface.get('buried_residues', []))
                    if surface_residues > 0 or buried_residues > 0:
                        residue_data = {
                            'Surface': surface_residues,
                            'Buried': buried_residues
                        }
                        if MATPLOTLIB_AVAILABLE:
                            residue_chart_path = self.save_chart_as_image(residue_data, "Surface vs Buried Residues", "pie", temp_dir)
                            if residue_chart_path:
                                f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(residue_chart_path)}' alt='Surface vs Buried Residues Chart'></div>\n")
                    
                    # Surface properties table
                    f.write("<table>\n")
                    f.write("<tr><th>Surface Property</th><th>Value</th></tr>\n")
                    f.write(f"<tr><td>Total SASA</td><td>{surface.get('total_sasa', 0):.2f} Ų</td></tr>\n")
                    f.write(f"<tr><td>Accessible Surface Area</td><td>{surface.get('accessible_surface_area', 0):.2f} Ų</td></tr>\n")
                    f.write(f"<tr><td>Buried Surface Area</td><td>{surface.get('buried_surface_area', 0):.2f} Ų</td></tr>\n")
                    f.write(f"<tr><td>Surface Percentage</td><td>{surface.get('surface_percentage', 0):.2f}%</td></tr>\n")
                    f.write(f"<tr><td>Buried Percentage</td><td>{surface.get('buried_percentage', 0):.2f}%</td></tr>\n")
                    f.write(f"<tr><td>Average RSA</td><td>{surface.get('average_rsa', 0):.2f}%</td></tr>\n")
                    f.write(f"<tr><td>Surface Residues</td><td>{surface_residues}</td></tr>\n")
                    f.write(f"<tr><td>Buried Residues</td><td>{buried_residues}</td></tr>\n")
                    f.write("</table>\n")
                
                # Cavity Analysis
                if 'cavity' in self.current_results:
                    cavity = self.current_results['cavity']
                    f.write("<h2>Cavity Analysis</h2>\n")
                    
                    # Cavity volumes chart if multiple cavities
                    cavities = cavity.get('potential_cavities', [])
                    if len(cavities) > 1 and MATPLOTLIB_AVAILABLE:
                        cavity_volumes = {f"Cavity {i+1}": cavity_info.get('volume', 0) for i, cavity_info in enumerate(cavities[:10])}
                        if any(vol > 0 for vol in cavity_volumes.values()):
                            volume_chart_path = self.save_chart_as_image(cavity_volumes, "Cavity Volumes", "bar", temp_dir, "Cavity", "Volume (Ų)")
                            if volume_chart_path:
                                f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(volume_chart_path)}' alt='Cavity Volumes Chart'></div>\n")
                    
                    # Cavity summary table
                    f.write("<table>\n")
                    f.write("<tr><th>Cavity Property</th><th>Value</th></tr>\n")
                    f.write(f"<tr><td>Total Cavities Found</td><td>{len(cavities)}</td></tr>\n")
                    f.write(f"<tr><td>Total Cavity Volume</td><td>{cavity.get('cavity_volume', 0):.2f} Ų</td></tr>\n")
                    
                    largest_cavity = cavity.get('largest_cavity')
                    if largest_cavity:
                        f.write(f"<tr><td>Largest Cavity Volume</td><td>{largest_cavity.get('volume', 0):.2f} Ų</td></tr>\n")
                        f.write(f"<tr><td>Largest Cavity Surface Area</td><td>{largest_cavity.get('surface_area', 0):.2f} Ų</td></tr>\n")
                    
                    druggable_cavities = cavity.get('druggable_cavities', [])
                    f.write(f"<tr><td>Druggable Cavities</td><td>{len(druggable_cavities)}</td></tr>\n")
                    f.write(f"<tr><td>Calculation Time</td><td>{cavity.get('calculation_time', 0):.2f} seconds</td></tr>\n")
                    f.write("</table>\n")
                    
                    # Detailed cavity information
                    if cavities:
                        f.write("<h3>Cavity Details</h3>\n")
                        f.write("<table>\n")
                        f.write("<tr><th>Cavity ID</th><th>Volume (Ų)</th><th>Surface Area (Ų)</th><th>Max Radius (Å)</th><th>Druggability Score</th><th>Lining Residues</th></tr>\n")
                        for cavity_info in cavities[:10]:  # Show top 10 cavities
                            cavity_id = cavity_info.get('id', 'N/A')
                            volume = cavity_info.get('volume', 0)
                            surface_area = cavity_info.get('surface_area', 0)
                            max_radius = cavity_info.get('max_radius', 0)
                            druggability = cavity_info.get('druggability_score', 0)
                            lining_residues = len(cavity_info.get('lining_residues', []))
                            f.write(f"<tr><td>{cavity_id}</td><td>{volume:.2f}</td><td>{surface_area:.2f}</td><td>{max_radius:.2f}</td><td>{druggability:.3f}</td><td>{lining_residues}</td></tr>\n")
                        f.write("</table>\n")
                
                # Footer
                f.write("<hr>\n")
                f.write(f"<p><em>Generated by PicoMol v1.0 on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')}</em></p>\n")
                
                f.write("</body>\n</html>")
        
        finally:
            # Clean up temporary files
            import shutil
            try:
                shutil.rmtree(temp_dir)
            except:
                pass
    
    def _export_comprehensive_metadata(self, f):
        """Export comprehensive metadata to HTML."""
        try:
            metadata = self.current_comprehensive_data.get('metadata', {})
            if not metadata or 'entry' not in metadata:
                return
            
            entry = metadata['entry']
            f.write("<h2>Structure Information</h2>\n")
            
            # Title and basic info
            f.write("<table>\n")
            f.write("<tr><th>Property</th><th>Value</th></tr>\n")
            
            # Title
            if entry.get('struct') and entry['struct'].get('title'):
                title = entry['struct']['title']
                f.write(f"<tr><td>Title</td><td>{title}</td></tr>\n")
            
            # Authors
            if entry.get('audit_author'):
                authors = [author.get('name', '') for author in entry['audit_author'][:5] if author.get('name')]
                if authors:
                    authors_text = ', '.join(authors)
                    if len(entry.get('audit_author', [])) > 5:
                        authors_text += " et al."
                    f.write(f"<tr><td>Authors</td><td>{authors_text}</td></tr>\n")
            
            # Dates
            if entry.get('rcsb_accession_info'):
                accession_info = entry['rcsb_accession_info']
                deposit_date = accession_info.get('deposit_date', 'N/A')
                release_date = accession_info.get('initial_release_date', 'N/A')
                f.write(f"<tr><td>Deposited</td><td>{deposit_date}</td></tr>\n")
                f.write(f"<tr><td>Released</td><td>{release_date}</td></tr>\n")
            
            # Experimental method and resolution
            if entry.get('rcsb_entry_info'):
                rcsb_info = entry['rcsb_entry_info']
                
                # Experimental method
                experimental_method = rcsb_info.get('experimental_method', 'N/A')
                f.write(f"<tr><td>Experimental Method</td><td>{experimental_method}</td></tr>\n")
                
                # Resolution
                resolution = rcsb_info.get('resolution_combined', 'N/A')
                if isinstance(resolution, list) and resolution:
                    resolution = resolution[0]
                if resolution != 'N/A' and isinstance(resolution, (int, float)):
                    resolution_text = f"{resolution:.2f} Å"
                else:
                    resolution_text = str(resolution)
                f.write(f"<tr><td>Resolution</td><td>{resolution_text}</td></tr>\n")
            
            # Keywords
            if entry.get('struct_keywords') and entry['struct_keywords'].get('pdbx_keywords'):
                keywords = entry['struct_keywords']['pdbx_keywords']
                f.write(f"<tr><td>Keywords</td><td>{keywords}</td></tr>\n")
            
            # Entity types
            if entry.get('rcsb_entry_info') and entry['rcsb_entry_info'].get('selected_polymer_entity_types'):
                entity_types = ', '.join(entry['rcsb_entry_info']['selected_polymer_entity_types'])
                f.write(f"<tr><td>Entity Types</td><td>{entity_types}</td></tr>\n")
            
            f.write("</table>\n")
            
        except Exception as e:
            print(f"[DEBUG] Error exporting comprehensive metadata: {e}")
    
    def _export_experimental_details(self, f, basic_results):
        """Export experimental details with enhanced fallback logic."""
        f.write("<h2>Experimental Details</h2>\n")
        f.write("<table>\n")
        f.write("<tr><th>Property</th><th>Value</th></tr>\n")
        
        # Experimental method and resolution removed
        
        # Refinement statistics if available
        if hasattr(self, 'current_comprehensive_data') and self.current_comprehensive_data:
            metadata = self.current_comprehensive_data.get('metadata', {})
            if 'entry' in metadata:
                entry = metadata['entry']
                if entry.get('refine') and entry['refine']:
                    refine_data = entry['refine'][0]
                    r_work = refine_data.get('ls_R_factor_R_work')
                    r_free = refine_data.get('ls_R_factor_R_free')
                    
                    if r_work is not None:
                        r_work_str = f"{r_work:.3f}" if isinstance(r_work, (int, float)) else str(r_work)
                        f.write(f"<tr><td>R-work</td><td>{r_work_str}</td></tr>\n")
                    if r_free is not None:
                        r_free_str = f"{r_free:.3f}" if isinstance(r_free, (int, float)) else str(r_free)
                        f.write(f"<tr><td>R-free</td><td>{r_free_str}</td></tr>\n")
        
        f.write("</table>\n")
    
    # Enhanced experimental method helper removed
    
    # Enhanced resolution helper removed
    
    def _extract_comprehensive_metadata_for_export(self):
        """Extract comprehensive metadata for JSON/other exports."""
        try:
            metadata = self.current_comprehensive_data.get('metadata', {})
            if not metadata or 'entry' not in metadata:
                return {}
            
            entry = metadata['entry']
            result = {}
            
            # Title
            if entry.get('struct') and entry['struct'].get('title'):
                result['title'] = entry['struct']['title']
            
            # Authors
            if entry.get('audit_author'):
                authors = [author.get('name', '') for author in entry['audit_author'] if author.get('name')]
                result['authors'] = authors
            
            # Dates
            if entry.get('rcsb_accession_info'):
                accession_info = entry['rcsb_accession_info']
                result['deposit_date'] = accession_info.get('deposit_date')
                result['release_date'] = accession_info.get('initial_release_date')
            
            # Keywords
            if entry.get('struct_keywords') and entry['struct_keywords'].get('pdbx_keywords'):
                result['keywords'] = entry['struct_keywords']['pdbx_keywords']
            
            # Entity types
            if entry.get('rcsb_entry_info') and entry['rcsb_entry_info'].get('selected_polymer_entity_types'):
                result['entity_types'] = entry['rcsb_entry_info']['selected_polymer_entity_types']
            
            return result
            
        except Exception as e:
            print(f"[DEBUG] Error extracting comprehensive metadata for export: {e}")
            return {}
    
    # Experimental details extraction method removed
    
    def _export_comprehensive_metadata(self, f):
        """Export comprehensive metadata to HTML."""
        try:
            metadata = self.current_comprehensive_data.get('metadata', {})
            if not metadata or 'entry' not in metadata:
                return
            
            entry = metadata['entry']
            f.write("<h2>Structure Information</h2>\n")
            
            # Title and basic info
            f.write("<table>\n")
            f.write("<tr><th>Property</th><th>Value</th></tr>\n")
            
            # Title
            if entry.get('struct') and entry['struct'].get('title'):
                title = entry['struct']['title']
                f.write(f"<tr><td>Title</td><td>{title}</td></tr>\n")
            
            # Authors
            if entry.get('audit_author'):
                authors = [author.get('name', '') for author in entry['audit_author'][:5] if author.get('name')]
                if authors:
                    authors_text = ', '.join(authors)
                    if len(entry.get('audit_author', [])) > 5:
                        authors_text += " et al."
                    f.write(f"<tr><td>Authors</td><td>{authors_text}</td></tr>\n")
            
            # Dates
            if entry.get('rcsb_accession_info'):
                accession_info = entry['rcsb_accession_info']
                deposit_date = accession_info.get('deposit_date', 'N/A')
                release_date = accession_info.get('initial_release_date', 'N/A')
                f.write(f"<tr><td>Deposited</td><td>{deposit_date}</td></tr>\n")
                f.write(f"<tr><td>Released</td><td>{release_date}</td></tr>\n")
            
            # Experimental method and resolution
            if entry.get('rcsb_entry_info'):
                rcsb_info = entry['rcsb_entry_info']
                
                # Experimental method
                experimental_method = rcsb_info.get('experimental_method', 'N/A')
                f.write(f"<tr><td>Experimental Method</td><td>{experimental_method}</td></tr>\n")
                
                # Resolution
                resolution = rcsb_info.get('resolution_combined', 'N/A')
                if isinstance(resolution, list) and resolution:
                    resolution = resolution[0]
                if resolution != 'N/A' and isinstance(resolution, (int, float)):
                    resolution_text = f"{resolution:.2f} Å"
                else:
                    resolution_text = str(resolution)
                f.write(f"<tr><td>Resolution</td><td>{resolution_text}</td></tr>\n")
            
            # Keywords
            if entry.get('struct_keywords') and entry['struct_keywords'].get('pdbx_keywords'):
                keywords = entry['struct_keywords']['pdbx_keywords']
                f.write(f"<tr><td>Keywords</td><td>{keywords}</td></tr>\n")
            
            # Entity types
            if entry.get('rcsb_entry_info') and entry['rcsb_entry_info'].get('selected_polymer_entity_types'):
                entity_types = ', '.join(entry['rcsb_entry_info']['selected_polymer_entity_types'])
                f.write(f"<tr><td>Entity Types</td><td>{entity_types}</td></tr>\n")
            
            f.write("</table>\n")
            
        except Exception as e:
            print(f"[DEBUG] Error exporting comprehensive metadata: {e}")
    
    def _export_experimental_details(self, f, basic_results):
        """Export experimental details with enhanced fallback logic."""
        f.write("<h2>Experimental Details</h2>\n")
        f.write("<table>\n")
        f.write("<tr><th>Property</th><th>Value</th></tr>\n")
        
        # Experimental method and resolution removed
        
        # Refinement statistics if available
        if hasattr(self, 'current_comprehensive_data') and self.current_comprehensive_data:
            metadata = self.current_comprehensive_data.get('metadata', {})
            if 'entry' in metadata:
                entry = metadata['entry']
                if entry.get('refine') and entry['refine']:
                    refine_data = entry['refine'][0]
                    r_work = refine_data.get('ls_R_factor_R_work')
                    r_free = refine_data.get('ls_R_factor_R_free')
                    
                    if r_work is not None:
                        r_work_str = f"{r_work:.3f}" if isinstance(r_work, (int, float)) else str(r_work)
                        f.write(f"<tr><td>R-work</td><td>{r_work_str}</td></tr>\n")
                    if r_free is not None:
                        r_free_str = f"{r_free:.3f}" if isinstance(r_free, (int, float)) else str(r_free)
                        f.write(f"<tr><td>R-free</td><td>{r_free_str}</td></tr>\n")
        
        f.write("</table>\n")
    
    # Duplicate experimental helper methods removed
    
    def export_complete_pdf(self, file_path):
        """Export complete analysis as PDF."""
        try:
            # First create HTML
            html_path = file_path.replace('.pdf', '_temp.html')
            self.export_complete_html(html_path)
            
            # Try to convert to PDF using weasyprint
            try:
                import weasyprint
                weasyprint.HTML(filename=html_path).write_pdf(file_path)
                os.remove(html_path)  # Clean up temp HTML
            except ImportError:
                # Fallback: just rename HTML to PDF and inform user
                os.rename(html_path, file_path.replace('.pdf', '.html'))
                raise ImportError("weasyprint not available. Exported as HTML instead. Install weasyprint for PDF export: pip install weasyprint")
                
        except Exception as e:
            raise Exception(f"PDF export failed: {str(e)}")
    
    def save_chart_as_image(self, data, title, chart_type, temp_dir, xlabel=None, ylabel=None):
        """Save a chart as PNG image and return the path with robust error handling."""
        if not MATPLOTLIB_AVAILABLE or not data:
            return None
        
        try:
            import matplotlib.pyplot as plt
            
            # Filter out zero, negative, and NaN values
            filtered_data = {}
            for k, v in data.items():
                if isinstance(v, (int, float)) and v > 0 and not (math.isnan(v) or math.isinf(v)):
                    filtered_data[k] = v
            
            if not filtered_data:
                print(f"[DEBUG] No valid data for chart '{title}'")
                return None
            
            # Ensure we have at least some meaningful data
            total_value = sum(filtered_data.values())
            if total_value <= 0:
                print(f"[DEBUG] Total value is zero for chart '{title}'")
                return None
            
            # Ultra high resolution settings
            plt.rcParams['font.size'] = 14
            plt.rcParams['axes.linewidth'] = 1.5
            plt.rcParams['font.weight'] = 'bold'
            plt.rcParams['axes.labelweight'] = 'bold'
            
            fig, ax = plt.subplots(figsize=(12, 9), dpi=300)  # Much larger and higher DPI
            
            if chart_type == "pie":
                # Prepare data
                labels = list(filtered_data.keys())
                values = list(filtered_data.values())
                
                # Double-check values are valid
                values = [max(0.001, v) for v in values if not (math.isnan(v) or math.isinf(v))]
                if len(values) != len(labels):
                    # Mismatch, rebuild both lists
                    valid_items = [(k, v) for k, v in filtered_data.items() 
                                  if isinstance(v, (int, float)) and v > 0 and not (math.isnan(v) or math.isinf(v))]
                    if not valid_items:
                        print(f"[DEBUG] No valid items for chart '{title}'")
                        return None
                    labels, values = zip(*valid_items)
                    labels, values = list(labels), list(values)
                
                # Use better colors and formatting
                colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD', '#98D8C8', '#F7DC6F']
                wedges, texts, autotexts = ax.pie(values, labels=labels, autopct='%1.1f%%', 
                                                 colors=colors[:len(values)], startangle=90, 
                                                 textprops={'fontsize': 14, 'weight': 'bold'})
                for autotext in autotexts:
                    autotext.set_color('white')
                    autotext.set_fontweight('bold')
                    autotext.set_fontsize(13)
                for text in texts:
                    text.set_fontsize(14)
                    text.set_fontweight('bold')
            elif chart_type == "bar":
                bars = ax.bar(filtered_data.keys(), filtered_data.values(), color='#4ECDC4', edgecolor='black', linewidth=2)
                ax.set_ylabel(ylabel or 'Count', fontsize=16, fontweight='bold')
                if xlabel:
                    ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')
                plt.xticks(rotation=45, ha='right', fontsize=14, fontweight='bold')
                plt.yticks(fontsize=14, fontweight='bold')
                
                # Add value labels on top of bars
                for bar in bars:
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height + max(filtered_data.values())*0.01,
                           f'{height:.1f}', ha='center', va='bottom', fontweight='bold', fontsize=12)
            
            ax.set_title(title, fontsize=18, fontweight='bold', pad=20)
            plt.tight_layout()
            
            # Save with ultra high quality
            chart_path = os.path.join(temp_dir, f"{title.replace(' ', '_').replace('/', '_')}.png")
            plt.savefig(chart_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none', 
                       format='png', quality=100)
            plt.close()
            plt.rcdefaults()  # Reset matplotlib settings
            
            return chart_path
        except Exception as e:
            print(f"[DEBUG] Error saving chart '{title}': {e}")
            return None
    
    def save_ramachandran_plot(self, rama_data, temp_dir, pdb_id=None):
        """Save Ramachandran plot as PNG image."""
        if not MATPLOTLIB_AVAILABLE or not rama_data:
            return None
        
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches
            
            plt.rcParams['font.size'] = 14
            plt.rcParams['font.weight'] = 'bold'
            plt.rcParams['axes.linewidth'] = 2
            fig, ax = plt.subplots(figsize=(12, 12), dpi=300)  # Much larger and higher DPI
            
            # Extract phi and psi angles
            phi_angles = [d['phi'] for d in rama_data]
            psi_angles = [d['psi'] for d in rama_data]
            
            # Plot points with better visibility
            ax.scatter(phi_angles, psi_angles, alpha=0.8, s=25, c='#2E86AB', edgecolors='black', linewidth=0.5)
            
            # Add favored regions with better colors
            # Alpha helix region
            alpha_helix = patches.Rectangle((-180, -70), 150, 120, linewidth=2, 
                                          edgecolor='#E74C3C', facecolor='#E74C3C', alpha=0.2)
            ax.add_patch(alpha_helix)
            
            # Beta sheet regions
            beta_sheet1 = patches.Rectangle((-180, 90), 150, 90, linewidth=2, 
                                          edgecolor='#27AE60', facecolor='#27AE60', alpha=0.2)
            beta_sheet2 = patches.Rectangle((-90, -180), 120, 90, linewidth=2, 
                                          edgecolor='#27AE60', facecolor='#27AE60', alpha=0.2)
            ax.add_patch(beta_sheet1)
            ax.add_patch(beta_sheet2)
            
            ax.set_xlim(-180, 180)
            ax.set_ylim(-180, 180)
            ax.set_xlabel('Phi (degrees)', fontsize=16, fontweight='bold')
            ax.set_ylabel('Psi (degrees)', fontsize=16, fontweight='bold')
            # Set title based on PDB ID if available
            if pdb_id:
                plot_title = f'Ramachandran Plot - {pdb_id}'
            else:
                plot_title = 'Ramachandran Plot'
            ax.set_title(plot_title, fontsize=20, fontweight='bold', pad=25)
            ax.grid(True, alpha=0.4, linewidth=1)
            
            # Set tick parameters for better visibility
            ax.tick_params(axis='both', which='major', labelsize=14, width=2, length=6)
            
            # Add cleaner legend
            ax.text(-170, 160, 'Favored regions:', fontsize=14, fontweight='bold', 
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='black'))
            ax.text(-170, 135, '■ Alpha helix', fontsize=13, color='#E74C3C', fontweight='bold')
            ax.text(-170, 110, '■ Beta sheet', fontsize=13, color='#27AE60', fontweight='bold')
            
            plt.tight_layout()
            
            # Save with ultra high quality
            chart_path = os.path.join(temp_dir, "ramachandran_plot.png")
            plt.savefig(chart_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none',
                       format='png', quality=100)
            plt.close()
            plt.rcdefaults()
            
            return chart_path
        except Exception as e:
            print(f"Error saving Ramachandran plot: {e}")
            return None
    
    def save_histogram(self, data, title, xlabel, ylabel, temp_dir):
        """Save histogram as PNG image."""
        if not MATPLOTLIB_AVAILABLE or not data:
            return None
        
        try:
            import matplotlib.pyplot as plt
            
            plt.rcParams['font.size'] = 14
            plt.rcParams['font.weight'] = 'bold'
            plt.rcParams['axes.linewidth'] = 2
            fig, ax = plt.subplots(figsize=(12, 8), dpi=300)  # Much larger and higher DPI
            
            # Create histogram with better styling
            n, bins, patches = ax.hist(data, bins=30, alpha=0.8, color='#4ECDC4', edgecolor='black', linewidth=1.2)
            
            ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')
            ax.set_ylabel(ylabel, fontsize=16, fontweight='bold')
            ax.set_title(title, fontsize=20, fontweight='bold', pad=25)
            ax.grid(True, alpha=0.4, linewidth=1)
            
            # Set tick parameters for better visibility
            ax.tick_params(axis='both', which='major', labelsize=14, width=2, length=6)
            
            # Add some statistics text
            mean_val = sum(data) / len(data)
            ax.axvline(mean_val, color='red', linestyle='--', linewidth=3, alpha=0.9, 
                      label=f'Mean: {mean_val:.1f}')
            ax.legend(fontsize=14, loc='upper right')
            
            plt.tight_layout()
            
            # Save with ultra high quality
            chart_path = os.path.join(temp_dir, f"{title.replace(' ', '_')}.png")
            plt.savefig(chart_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none',
                       format='png', quality=100)
            plt.close()
            plt.rcdefaults()
            
            return chart_path
        except Exception as e:
            print(f"Error saving histogram: {e}")
            return None
    
    def image_to_base64(self, image_path):
        """Convert image to base64 string for embedding in HTML."""
        try:
            import base64
            with open(image_path, 'rb') as img_file:
                return base64.b64encode(img_file.read()).decode('utf-8')
        except Exception as e:
            print(f"Error converting image to base64: {e}")
            return ""


    def export_to_html(self, file_path, options):
        """Legacy HTML export method (simplified version)."""
        # This is kept for compatibility but the main export now uses export_complete_html
        self.export_complete_html(file_path)


class ExportDialog(QDialog):
    """Dialog for configuring export options."""
    
    def __init__(self, parent=None, results=None):
        super().__init__(parent)
        self.setWindowTitle("Export Analysis Results")
        self.setModal(True)
        self.resize(500, 600)
        self.results = results
        
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Title
        title = QLabel("<h3>Export Analysis Results</h3>")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)
        
        # Format selection
        format_group = QGroupBox("Export Format")
        format_layout = QVBoxLayout(format_group)
        
        self.format_combo = QComboBox()
        formats = ["JSON", "CSV", "Excel", "Text Report", "HTML Report"]
        self.format_combo.addItems(formats)
        self.format_combo.setCurrentText("JSON")
        self.format_combo.currentTextChanged.connect(self.on_format_changed)
        format_layout.addWidget(self.format_combo)
        
        # Format descriptions
        self.format_desc = QLabel()
        self.format_desc.setWordWrap(True)
        self.format_desc.setStyleSheet("color: #666; font-size: 10px; margin: 5px;")
        format_layout.addWidget(self.format_desc)
        
        layout.addWidget(format_group)
        
        # Content selection
        content_group = QGroupBox("Content to Export")
        content_layout = QVBoxLayout(content_group)
        
        # Analysis type checkboxes
        self.include_basic = QCheckBox("Basic Properties")
        self.include_basic.setChecked(True)
        self.include_basic.setEnabled('basic' in self.results if self.results else False)
        content_layout.addWidget(self.include_basic)
        
        self.include_secondary = QCheckBox("Secondary Structure")
        self.include_secondary.setChecked(True)
        self.include_secondary.setEnabled('secondary' in self.results if self.results else False)
        content_layout.addWidget(self.include_secondary)
        
        self.include_geometry = QCheckBox("Geometric Properties")
        self.include_geometry.setChecked(True)
        self.include_geometry.setEnabled('geometry' in self.results if self.results else False)
        content_layout.addWidget(self.include_geometry)
        

        
        # Surface Properties export removed - will be added in a future release
        
        # Additional options
        content_layout.addWidget(QLabel(""))  # Spacer
        
        self.include_ramachandran = QCheckBox("Include Ramachandran Data (detailed)")
        self.include_ramachandran.setChecked(False)
        self.include_ramachandran.setToolTip("Include detailed phi/psi angle data for each residue")
        content_layout.addWidget(self.include_ramachandran)
        
        self.include_residue_details = QCheckBox("Include Residue-level Details")
        self.include_residue_details.setChecked(False)
        self.include_residue_details.setToolTip("Include detailed per-residue data (SASA, RSA, etc.)")
        content_layout.addWidget(self.include_residue_details)
        
        layout.addWidget(content_group)
        
        # File selection
        file_group = QGroupBox("Output File")
        file_layout = QHBoxLayout(file_group)
        
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select output file...")
        file_layout.addWidget(self.file_path_edit)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self.browse_file)
        file_layout.addWidget(browse_btn)
        
        layout.addWidget(file_group)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(cancel_btn)
        
        button_layout.addStretch()
        
        export_btn = QPushButton("Export")
        export_btn.setDefault(True)
        export_btn.clicked.connect(self.accept)
        button_layout.addWidget(export_btn)
        
        layout.addLayout(button_layout)
        
        # Update format description
        self.on_format_changed("JSON")
    
    def on_format_changed(self, format_name):
        """Update format description and file extension."""
        descriptions = {
            "JSON": "JavaScript Object Notation - structured data format, good for programmatic access",
            "CSV": "Comma-Separated Values - spreadsheet compatible, good for data analysis",
            "Excel": "Microsoft Excel format with multiple sheets - comprehensive and organized",
            "Text Report": "Human-readable text report - good for documentation and sharing",
            "HTML Report": "Web page format - good for viewing and sharing with formatting"
        }
        
        extensions = {
            "JSON": ".json",
            "CSV": ".csv",
            "Excel": ".xlsx",
            "Text Report": ".txt",
            "HTML Report": ".html"
        }
        
        self.format_desc.setText(descriptions.get(format_name, ""))
        
        # Update file extension if path is set
        current_path = self.file_path_edit.text()
        if current_path:
            base_path = os.path.splitext(current_path)[0]
            new_extension = extensions.get(format_name, ".txt")
            self.file_path_edit.setText(base_path + new_extension)
    
    def browse_file(self):
        """Browse for output file."""
        format_name = self.format_combo.currentText()
        
        filters = {
            "JSON": "JSON Files (*.json)",
            "CSV": "CSV Files (*.csv)",
            "Excel": "Excel Files (*.xlsx)",
            "Text Report": "Text Files (*.txt)",
            "HTML Report": "HTML Files (*.html)"
        }
        
        extensions = {
            "JSON": ".json",
            "CSV": ".csv",
            "Excel": ".xlsx",
            "Text Report": ".txt",
            "HTML Report": ".html"
        }
        
        # Suggest filename based on structure ID
        structure_id = self.results.get('structure_id', 'analysis') if self.results else 'analysis'
        suggested_name = f"{structure_id}_structural_analysis{extensions.get(format_name, '.txt')}"
        
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            f"Export {format_name}",
            suggested_name,
            filters.get(format_name, "All Files (*)")
        )
        
        if file_path:
            self.file_path_edit.setText(file_path)
    
    def get_export_settings(self):
        """Get export settings from dialog."""
        export_format = self.format_combo.currentText()
        file_path = self.file_path_edit.text()
        
        options = {
            'include_basic': self.include_basic.isChecked(),
            'include_secondary': self.include_secondary.isChecked(),
            'include_geometry': self.include_geometry.isChecked(),


            'include_ramachandran': self.include_ramachandran.isChecked(),
            'include_residue_details': self.include_residue_details.isChecked()
        }
        
        return export_format, file_path, options
    
    def accept(self):
        """Validate settings before accepting."""
        if not self.file_path_edit.text().strip():
            QMessageBox.warning(self, "Export Error", "Please select an output file.")
            return
        
        # Check if at least one content type is selected
        if not any([
            self.include_basic.isChecked(),
            self.include_secondary.isChecked(),
            self.include_geometry.isChecked(),

            # Surface Properties export removed
        ]):
            QMessageBox.warning(self, "Export Error", "Please select at least one content type to export.")
            return
        
        super().accept()


def create_structural_analysis_tab(parent=None):
    """Create the structural analysis tab."""
    return StructuralAnalysisTab(parent)


if __name__ == "__main__":
    # Test the structural analysis tools
    import sys
    from PyQt5.QtWidgets import QApplication
    
    app = QApplication(sys.argv)
    
    widget = create_structural_analysis_tab()
    widget.setWindowTitle("PicoMol Structural Analysis Tools")
    widget.resize(800, 600)
    widget.show()
    
    sys.exit(app.exec_())
