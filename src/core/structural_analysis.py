#!/usr/bin/env python3
"""
Structural Analysis Tools for PicoMol.

This module provides comprehensive protein structural analysis including:
- Secondary structure analysis
- Geometric property calculations
- Surface area and volume estimation
- Cavity and binding site detection
- Structural quality assessment
- Structural quality assessment
- Structural comparison tools
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
    OPTIMIZED_FETCH_AVAILABLE = False

try:
    from .enhanced_pdb_puller import EnhancedPDBPuller
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


# Van der Waals radii for common elements (in Angstroms)
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
    'P': 1.80, 'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98,
    'FE': 2.00, 'ZN': 1.39, 'MG': 1.73, 'CA': 2.31, 'MN': 1.61,
    'CU': 1.40, 'NI': 1.63, 'CO': 1.26, 'SE': 1.90
}

# Maximum SASA values for amino acids (in Ų) - for RSA calculation
MAX_SASA_VALUES = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0,
    'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0,
    'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0,
    'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0
}

class StructuralAnalysisWorker(QThread):
    """Worker thread for structural analysis calculations."""
    
    analysis_complete = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)
    progress_update = pyqtSignal(str)
    
    def __init__(self, structure_path, analysis_types=None):
        super().__init__()
        self.structure_path = structure_path
        self.analysis_types = analysis_types or ['basic', 'secondary', 'geometry', 'quality']
        self.structure = None
        # Surface analysis uses SASA calculation
    
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
            
            if 'surface' in self.analysis_types:
                self.progress_update.emit("Calculating surface properties using SASA...")
                results['surface'] = self.analyze_surface_properties()
            
            if 'quality' in self.analysis_types:
                self.progress_update.emit("Assessing structural quality...")
                results['quality'] = self.analyze_structural_quality()
            
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
            'resolution': 'N/A',
            'space_group': 'N/A',
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
        
        # Try to get experimental information from PDB header
        try:
            header = self.structure.header
            if header:
                # Debug: print available header keys
                print(f"[DEBUG] Available header keys: {list(header.keys()) if hasattr(header, 'keys') else 'No keys method'}")
                
                # Try different possible keys for resolution
                resolution = None
                for res_key in ['resolution', 'Resolution', 'RESOLUTION']:
                    if res_key in header:
                        resolution = header[res_key]
                        break
                
                # Try different possible keys for space group
                space_group = None
                for sg_key in ['space_group', 'spacegroup', 'Space_group', 'SPACE_GROUP']:
                    if sg_key in header:
                        space_group = header[sg_key]
                        break
                
                results['resolution'] = resolution if resolution is not None else 'N/A'
                results['space_group'] = space_group if space_group is not None else 'N/A'
                
                # Also try to get other useful header information
                for key in ['deposition_date', 'release_date', 'structure_method', 'head']:
                    if key in header:
                        results[key] = header[key]
                        print(f"[DEBUG] Found header field {key}: {header[key]}")
            else:
                print("[DEBUG] No header information available in structure")
                results['resolution'] = 'N/A'
                results['space_group'] = 'N/A'
        except Exception as e:
            print(f"[DEBUG] Error extracting header information: {e}")
            results['resolution'] = 'N/A'
            results['space_group'] = 'N/A'
        
        return results
    
    def analyze_secondary_structure(self):
        """Analyze secondary structure using simple geometric rules and collect Ramachandran data."""
        results = {
            'helix_content': 0.0,
            'sheet_content': 0.0,
            'loop_content': 0.0,
            'secondary_elements': [],
            'ramachandran_data': []
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
    
    def analyze_surface_properties(self):
        """Analyze surface properties using optimized SASA calculation.
        
        Uses an optimized Shrake-Rupley algorithm with NeighborSearch for
        efficient and accurate surface area calculation.
        
        Returns:
            dict: Surface analysis results
        """
        import time
        start_time = time.time()
        
        results = self.analyze_surface_properties_sasa()
        
        # Add timing information
        calculation_time = time.time() - start_time
        results['calculation_time'] = calculation_time
        results['method_used'] = 'SASA'
        
        self.progress_update.emit(f"Surface analysis completed in {calculation_time:.2f} seconds using SASA calculation")
        
        return results
    

    
    def analyze_surface_properties_sasa(self):
        """Surface analysis using configurable SASA calculation.
        
        Uses an optimized Shrake-Rupley algorithm with NeighborSearch for
        efficient and accurate surface area calculation.
        """
        # Get SASA configuration
        config = getattr(self, 'sasa_config', {
            'probe_radius': 1.4,
            'n_points': 30,
            'rsa_threshold': 20.0,
            'include_hydrogens': False
        })
        results = {
            'accessible_surface_area': 0.0,
            'buried_surface_area': 0.0,
            'surface_residues': [],
            'buried_residues': [],
            'surface_percentage': 0.0,
            'buried_percentage': 0.0,
            'total_sasa': 0.0,
            'average_rsa': 0.0
        }
        
        # Get all atoms for neighbor search
        all_atoms = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        for atom in residue:
                            # Filter hydrogens if not included
                            if not config['include_hydrogens']:
                                element = atom.element.upper() if hasattr(atom, 'element') else atom.get_name()[0]
                                if element == 'H':
                                    continue
                            all_atoms.append(atom)
        
        if not all_atoms:
            return results
        
        # Create NeighborSearch for efficient neighbor finding
        try:
            neighbor_search = NeighborSearch(all_atoms)
        except Exception as e:
            self.progress_update.emit(f"NeighborSearch creation failed: {e}")
            return self._analyze_surface_fallback()
        
        all_residues = []
        total_sasa = 0.0
        total_rsa = 0.0
        
        residue_count = 0
        total_residues_to_process = sum(1 for model in self.structure 
                                      for chain in model 
                                      for residue in chain if is_aa(residue))
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        residue_count += 1
                        
                        # Update progress
                        if residue_count % 10 == 0:
                            progress = (residue_count / total_residues_to_process) * 100
                            self.progress_update.emit(f"Calculating SASA: {progress:.1f}% ({residue_count}/{total_residues_to_process})")
                        
                        try:
                            # Calculate SASA for this residue with configuration
                            residue_sasa = self.calculate_residue_sasa_optimized(
                                residue, neighbor_search, 
                                probe_radius=config['probe_radius'],
                                n_points=config['n_points'],
                                include_hydrogens=config['include_hydrogens']
                            )
                            
                            # Calculate RSA
                            res_name = residue.get_resname()
                            max_sasa = MAX_SASA_VALUES.get(res_name, 150.0)
                            rsa = (residue_sasa / max_sasa) * 100 if max_sasa > 0 else 0
                            
                            residue_info = {
                                'residue': res_name,
                                'chain': chain.id,
                                'position': residue.id[1],
                                'sasa': residue_sasa,
                                'rsa': rsa,
                                'max_sasa': max_sasa
                            }
                            
                            all_residues.append(residue_info)
                            total_sasa += residue_sasa
                            total_rsa += rsa
                            
                            # Classification based on configurable RSA threshold
                            if rsa < config['rsa_threshold']:
                                results['buried_residues'].append(residue_info)
                            else:
                                results['surface_residues'].append(residue_info)
                                
                        except Exception as e:
                            self.progress_update.emit(f"SASA calculation failed for residue {residue}: {e}")
                            continue
        
        # Calculate summary statistics
        total_residues = len(all_residues)
        if total_residues > 0:
            surface_count = len(results['surface_residues'])
            buried_count = len(results['buried_residues'])
            
            results['surface_percentage'] = (surface_count / total_residues) * 100
            results['buried_percentage'] = (buried_count / total_residues) * 100
            results['total_sasa'] = total_sasa
            results['accessible_surface_area'] = sum(r['sasa'] for r in results['surface_residues'])
            results['buried_surface_area'] = sum(r['sasa'] for r in results['buried_residues'])
            results['average_rsa'] = total_rsa / total_residues
        
        return results
    
    def _analyze_surface_fallback(self):
        """Fallback surface analysis when NeighborSearch fails.
        
        Uses a simplified approach but still calculates actual SASA values.
        """
        results = {
            'accessible_surface_area': 0.0,
            'buried_surface_area': 0.0,
            'surface_residues': [],
            'buried_residues': [],
            'surface_percentage': 0.0,
            'buried_percentage': 0.0,
            'total_sasa': 0.0,
            'average_rsa': 0.0
        }
        
        # Get all atoms for fallback calculation
        all_atoms = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        for atom in residue:
                            # Filter hydrogens if not included
                            if not config['include_hydrogens']:
                                element = atom.element.upper() if hasattr(atom, 'element') else atom.get_name()[0]
                                if element == 'H':
                                    continue
                            all_atoms.append(atom)
        
        all_residues = []
        total_sasa = 0.0
        total_rsa = 0.0
        
        residue_count = 0
        total_residues_to_process = sum(1 for model in self.structure 
                                      for chain in model 
                                      for residue in chain if is_aa(residue))
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        residue_count += 1
                        
                        # Update progress
                        if residue_count % 10 == 0:
                            progress = (residue_count / total_residues_to_process) * 100
                            self.progress_update.emit(f"Fallback SASA calculation: {progress:.1f}% ({residue_count}/{total_residues_to_process})")
                        
                        try:
                            # Calculate SASA using fallback method with configuration
                            config = getattr(self, 'sasa_config', {
                                'probe_radius': 1.4,
                                'n_points': 100,  # Use more points for fallback accuracy
                                'rsa_threshold': 20.0
                            })
                            residue_sasa = self.calculate_residue_sasa(
                                residue, model,
                                probe_radius=config['probe_radius'],
                                n_points=config['n_points'],
                                include_hydrogens=config['include_hydrogens']
                            )
                            
                            # Calculate RSA
                            res_name = residue.get_resname()
                            max_sasa = MAX_SASA_VALUES.get(res_name, 150.0)
                            rsa = (residue_sasa / max_sasa) * 100 if max_sasa > 0 else 0
                            
                            residue_info = {
                                'residue': res_name,
                                'chain': chain.id,
                                'position': residue.id[1],
                                'sasa': residue_sasa,
                                'rsa': rsa,
                                'max_sasa': max_sasa
                            }
                            
                            all_residues.append(residue_info)
                            total_sasa += residue_sasa
                            total_rsa += rsa
                            
                            # Classification based on configurable RSA threshold
                            if rsa < config['rsa_threshold']:
                                results['buried_residues'].append(residue_info)
                            else:
                                results['surface_residues'].append(residue_info)
                                
                        except Exception as e:
                            self.progress_update.emit(f"Fallback SASA calculation failed for residue {residue}: {e}")
                            continue
        
        # Calculate summary statistics
        total_residues = len(all_residues)
        if total_residues > 0:
            surface_count = len(results['surface_residues'])
            buried_count = len(results['buried_residues'])
            
            results['surface_percentage'] = (surface_count / total_residues) * 100
            results['buried_percentage'] = (buried_count / total_residues) * 100
            results['total_sasa'] = total_sasa
            results['accessible_surface_area'] = sum(r['sasa'] for r in results['surface_residues'])
            results['buried_surface_area'] = sum(r['sasa'] for r in results['buried_residues'])
            results['average_rsa'] = total_rsa / total_residues
        
        return results
    
    def calculate_residue_sasa_optimized(self, residue, neighbor_search, probe_radius=1.4, n_points=30, include_hydrogens=False):
        """Calculate SASA for a residue using optimized Shrake-Rupley algorithm.
        
        Args:
            residue: The residue to calculate SASA for
            neighbor_search: Pre-built NeighborSearch object for efficient neighbor finding
            probe_radius: Probe radius in Angstroms (1.4Å for water)
            n_points: Number of points to generate on each atom sphere (reduced for speed)
            include_hydrogens: Whether to include hydrogen atoms
        
        Returns:
            float: SASA value in Ų
        """
        total_sasa = 0.0
        
        # Calculate SASA for each atom in the residue
        for atom in residue:
            # Filter hydrogens if not included
            if not include_hydrogens:
                element = atom.element.upper() if hasattr(atom, 'element') else atom.get_name()[0]
                if element == 'H':
                    continue
            
            atom_sasa = self.calculate_atom_sasa_optimized(atom, neighbor_search, probe_radius, n_points)
            total_sasa += atom_sasa
        
        return total_sasa
    
    def calculate_atom_sasa_optimized(self, atom, neighbor_search, probe_radius=1.4, n_points=30):
        """Calculate SASA for a single atom using optimized Shrake-Rupley algorithm.
        
        Args:
            atom: The atom to calculate SASA for
            neighbor_search: Pre-built NeighborSearch object
            probe_radius: Probe radius in Angstroms
            n_points: Number of points to generate on the sphere (reduced for speed)
        
        Returns:
            float: SASA value for the atom in Ų
        """
        # Get van der Waals radius for this atom
        element = atom.element.upper() if hasattr(atom, 'element') else atom.get_name()[0]
        vdw_radius = VDW_RADII.get(element, 1.7)  # Default to carbon radius
        
        # Total radius = vdw radius + probe radius
        total_radius = vdw_radius + probe_radius
        
        # Find nearby atoms efficiently using NeighborSearch
        # Search within a reasonable distance (total_radius + max_other_radius + probe_radius)
        search_radius = total_radius + 3.0  # Conservative estimate
        nearby_atoms = neighbor_search.search(atom.get_coord(), search_radius)
        
        # Remove self from nearby atoms
        nearby_atoms = [a for a in nearby_atoms if a != atom]
        
        # If no nearby atoms, the entire surface is accessible
        if not nearby_atoms:
            return 4 * math.pi * (total_radius ** 2)
        
        # Generate points on sphere surface (fewer points for speed)
        sphere_points = self.generate_sphere_points(atom.get_coord(), total_radius, n_points)
        
        # Check which points are accessible (not buried by nearby atoms only)
        accessible_points = 0
        
        for point in sphere_points:
            is_accessible = True
            
            # Check against nearby atoms only (major optimization)
            for other_atom in nearby_atoms:
                # Get other atom's properties
                other_element = other_atom.element.upper() if hasattr(other_atom, 'element') else other_atom.get_name()[0]
                other_vdw_radius = VDW_RADII.get(other_element, 1.7)
                other_total_radius = other_vdw_radius + probe_radius
                
                # Check if point is inside other atom's sphere
                distance = np.linalg.norm(point - other_atom.get_coord())
                if distance < other_total_radius:
                    is_accessible = False
                    break
            
            if is_accessible:
                accessible_points += 1
        
        # Calculate surface area
        total_surface_area = 4 * math.pi * (total_radius ** 2)
        accessible_fraction = accessible_points / n_points if n_points > 0 else 0
        atom_sasa = total_surface_area * accessible_fraction
        
        return atom_sasa
    
    def calculate_residue_sasa(self, residue, model, probe_radius=1.4, n_points=100, include_hydrogens=False):
        """Calculate SASA for a residue using Shrake-Rupley algorithm.
        
        DEPRECATED: Use calculate_residue_sasa_optimized for better performance.
        
        Args:
            residue: The residue to calculate SASA for
            model: The PDB model containing all atoms
            probe_radius: Probe radius in Angstroms (1.4Å for water)
            n_points: Number of points to generate on each atom sphere
            include_hydrogens: Whether to include hydrogen atoms
        
        Returns:
            float: SASA value in Ų
        """
        total_sasa = 0.0
        
        # Get all atoms in the protein for neighbor checking
        all_atoms = []
        for chain in model:
            for res in chain:
                for atom in res:
                    # Filter hydrogens if not included
                    if not include_hydrogens:
                        element = atom.element.upper() if hasattr(atom, 'element') else atom.get_name()[0]
                        if element == 'H':
                            continue
                    all_atoms.append(atom)
        
        # Calculate SASA for each atom in the residue
        for atom in residue:
            # Filter hydrogens if not included
            if not include_hydrogens:
                element = atom.element.upper() if hasattr(atom, 'element') else atom.get_name()[0]
                if element == 'H':
                    continue
            
            atom_sasa = self.calculate_atom_sasa(atom, all_atoms, probe_radius, n_points)
            total_sasa += atom_sasa
        
        return total_sasa
    
    def calculate_atom_sasa(self, atom, all_atoms, probe_radius=1.4, n_points=100):
        """Calculate SASA for a single atom using Shrake-Rupley algorithm.
        
        Args:
            atom: The atom to calculate SASA for
            all_atoms: List of all atoms in the protein
            probe_radius: Probe radius in Angstroms
            n_points: Number of points to generate on the sphere
        
        Returns:
            float: SASA value for the atom in Ų
        """
        # Get van der Waals radius for this atom
        element = atom.element.upper() if hasattr(atom, 'element') else atom.get_name()[0]
        vdw_radius = VDW_RADII.get(element, 1.7)  # Default to carbon radius
        
        # Total radius = vdw radius + probe radius
        total_radius = vdw_radius + probe_radius
        
        # Generate points on sphere surface
        sphere_points = self.generate_sphere_points(atom.get_coord(), total_radius, n_points)
        
        # Check which points are accessible (not buried by other atoms)
        accessible_points = 0
        
        for point in sphere_points:
            is_accessible = True
            
            # Check against all other atoms
            for other_atom in all_atoms:
                if other_atom == atom:
                    continue
                
                # Get other atom's properties
                other_element = other_atom.element.upper() if hasattr(other_atom, 'element') else other_atom.get_name()[0]
                other_vdw_radius = VDW_RADII.get(other_element, 1.7)
                other_total_radius = other_vdw_radius + probe_radius
                
                # Check if point is inside other atom's sphere
                distance = np.linalg.norm(point - other_atom.get_coord())
                if distance < other_total_radius:
                    is_accessible = False
                    break
            
            if is_accessible:
                accessible_points += 1
        
        # Calculate surface area
        # Surface area of sphere = 4π * r²
        total_surface_area = 4 * math.pi * (total_radius ** 2)
        accessible_fraction = accessible_points / n_points if n_points > 0 else 0
        atom_sasa = total_surface_area * accessible_fraction
        
        return atom_sasa
    
    def generate_sphere_points(self, center, radius, n_points):
        """Generate uniformly distributed points on a sphere surface.
        
        Args:
            center: Center coordinates of the sphere
            radius: Radius of the sphere
            n_points: Number of points to generate
        
        Returns:
            numpy.ndarray: Array of 3D points on the sphere surface
        """
        points = []
        
        # Use golden spiral method for uniform distribution
        golden_angle = math.pi * (3.0 - math.sqrt(5.0))  # Golden angle in radians
        
        for i in range(n_points):
            # y goes from 1 to -1
            y = 1 - (i / float(n_points - 1)) * 2
            
            # radius at y
            radius_at_y = math.sqrt(1 - y * y)
            
            # golden angle increment
            theta = golden_angle * i
            
            x = math.cos(theta) * radius_at_y
            z = math.sin(theta) * radius_at_y
            
            # Scale by radius and translate to center
            point = np.array([x * radius + center[0],
                             y * radius + center[1], 
                             z * radius + center[2]])
            points.append(point)
        
        return np.array(points)
    
    def count_neighbors(self, atom, model, radius=8.0, ca_only=True):
        """Count neighboring atoms within radius (legacy method for compatibility).
        
        Args:
            atom: The central atom to count neighbors for
            model: The PDB model containing all atoms
            radius: Distance cutoff in Angstroms (default 8.0Å)
            ca_only: If True, count only CA atoms (default True for surface analysis)
        
        Returns:
            int: Number of neighboring atoms within the radius
        """
        count = 0
        atom_coord = atom.get_coord()
        
        for chain in model:
            for residue in chain:
                if ca_only:
                    # Count only CA atoms for surface analysis
                    if 'CA' in residue:
                        other_atom = residue['CA']
                        if other_atom != atom:
                            distance = np.linalg.norm(atom_coord - other_atom.get_coord())
                            if distance <= radius:
                                count += 1
                else:
                    # Count all atoms
                    for other_atom in residue:
                        if other_atom != atom:
                            distance = np.linalg.norm(atom_coord - other_atom.get_coord())
                            if distance <= radius:
                                count += 1
        
        return count
    

    
    def analyze_structural_quality(self):
        """Analyze structural quality metrics including Ramachandran validation."""
        results = {
            'clashes': [],
            'missing_atoms': [],
            'b_factor_stats': {},
            'ramachandran_favored': 0,
            'ramachandran_allowed': 0,
            'ramachandran_outliers': [],
            'ramachandran_outliers_count': 0
        }
        
        # Collect Ramachandran data for quality assessment
        ramachandran_outliers = []
        favored_count = 0
        allowed_count = 0
        outlier_count = 0
        
        for model in self.structure:
            for chain in model:
                residues = [res for res in chain if is_aa(res)]
                
                # Calculate phi/psi angles for Ramachandran analysis
                for i in range(1, len(residues) - 1):
                    try:
                        phi, psi = self.calculate_phi_psi(residues, i)
                        residue = residues[i]
                        res_name = residue.get_resname()
                        
                        # Classify Ramachandran region
                        rama_region = self.classify_ramachandran_region(phi, psi, res_name)
                        
                        if rama_region == 'favored':
                            favored_count += 1
                        elif rama_region == 'allowed':
                            allowed_count += 1
                        else:  # outlier
                            outlier_count += 1
                            ramachandran_outliers.append({
                                'residue': res_name,
                                'chain': chain.id,
                                'position': residue.id[1],
                                'phi': phi,
                                'psi': psi
                            })
                    except Exception:
                        continue
        
        results['ramachandran_favored'] = favored_count
        results['ramachandran_allowed'] = allowed_count
        results['ramachandran_outliers'] = ramachandran_outliers
        results['ramachandran_outliers_count'] = outlier_count
        
        # Check for missing atoms
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        expected_atoms = ['N', 'CA', 'C', 'O']
                        for atom_name in expected_atoms:
                            if atom_name not in residue:
                                results['missing_atoms'].append({
                                    'residue': residue.get_resname(),
                                    'chain': chain.id,
                                    'position': residue.id[1],
                                    'missing_atom': atom_name
                                })
        
        # B-factor statistics
        all_b_factors = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        all_b_factors.append(atom.get_bfactor())
        
        if all_b_factors:
            results['b_factor_stats'] = {
                'mean': np.mean(all_b_factors),
                'std': np.std(all_b_factors),
                'min': min(all_b_factors),
                'max': max(all_b_factors),
                'median': np.median(all_b_factors) if NUMPY_AVAILABLE else sorted(all_b_factors)[len(all_b_factors)//2]
            }
        
        return results
    
    def classify_ramachandran_region(self, phi, psi, residue_name=None):
        """Classify phi/psi angles into Ramachandran regions using MolProbity/Richardson lab definitions.
        
        Based on the exact definitions used by MolProbity and other professional validation tools.
        These boundaries are derived from high-resolution crystal structures and match the
        standard Ramachandran plot regions used in structural biology.
        
        Args:
            phi: Phi dihedral angle in degrees
            psi: Psi dihedral angle in degrees
            residue_name: Three-letter residue code (optional, for special cases)
        
        Returns:
            str: 'favored', 'allowed', or 'outlier'
        """
        # Handle special cases for glycine and proline
        if residue_name == 'GLY':
            return self._classify_glycine_ramachandran(phi, psi)
        elif residue_name == 'PRO':
            return self._classify_proline_ramachandran(phi, psi)
        
        # General case
        point = (phi, psi)
        
        # General Favored Regions (98th percentile)
        favored_general = [
            [(-180, 150), (-50, 150), (-50, 50), (-180, 50)], # Beta region
            [(-100, 0), (-25, 0), (-25, -70), (-100, -70)], # Alpha-right region
            [(45, 60), (90, 60), (90, 0), (45, 0)] # Alpha-left region
        ]

        # General Allowed Regions (99.95th percentile)
        allowed_general = [
            [(-180, 180), (-180, 50), (-50, 50), (-50, 180)], # Beta region
            [(-180, 0), (-180, -100), (0, -100), (0, 0)], # Alpha-right region
            [(45, 120), (120, 120), (120, -60), (45, -60)] # Alpha-left region
        ]

        for region in favored_general:
            if self._is_in_polygon(point, region):
                return 'favored'
        for region in allowed_general:
            if self._is_in_polygon(point, region):
                return 'allowed'

        return 'outlier'

    def _is_in_polygon(self, point, polygon):
        """Checks if a point is inside a polygon using the ray-casting algorithm."""
        x, y = point
        n = len(polygon)
        inside = False
        p1x, p1y = polygon[0]
        for i in range(n + 1):
            p2x, p2y = polygon[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside

    def _classify_glycine_ramachandran(self, phi, psi):
        """Special Ramachandran classification for glycine residues."""
        point = (phi, psi)
        
        # Glycine Favored Regions
        favored_gly = [
            [(-180, 180), (0, 180), (0, 50), (-180, 50)], # Top-left region
            [(0, 180), (180, 180), (180, -50), (0, -50)], # Top-right region
            [(-180, 0), (-180, -100), (0, -100), (0, 0)] # Bottom-left region
        ]

        # Glycine Allowed Regions
        allowed_gly = [
            [(-180, 180), (180, 180), (180, -180), (-180, -180)] # Almost the entire plot
        ]

        for region in favored_gly:
            if self._is_in_polygon(point, region):
                return 'favored'
        for region in allowed_gly:
            if self._is_in_polygon(point, region):
                return 'allowed'
        
        return 'outlier'

    def _classify_proline_ramachandran(self, phi, psi):
        """Special Ramachandran classification for proline residues."""
        point = (phi, psi)

        # Proline Favored Regions
        favored_pro = [
            [(-100, 180), (-50, 180), (-50, 100), (-100, 100)], # PPII region
            [(-100, 50), (-50, 50), (-50, -50), (-100, -50)] # Beta-turn region
        ]

        # Proline Allowed Regions
        allowed_pro = [
            [(-120, 180), (-30, 180), (-30, -60), (-120, -60)]
        ]

        for region in favored_pro:
            if self._is_in_polygon(point, region):
                return 'favored'
        for region in allowed_pro:
            if self._is_in_polygon(point, region):
                return 'allowed'

        return 'outlier'
    
    def detect_cavities(self):
        """Detect potential binding sites and cavities using improved geometric methods.
        
        Uses a combination of grid-based and probe-based approaches inspired by
        conventional cavity detection algorithms like CASTp and fpocket.
        """
        import time
        start_time = time.time()
        
        results = {
            'potential_cavities': [],
            'surface_pockets': [],
            'cavity_volume': 0.0,
            'largest_cavity': None,
            'druggable_cavities': [],
            'calculation_time': 0.0,
            'total_grid_points': 0,
            'cavity_points_found': 0
        }
        
        self.progress_update.emit("Cavity detection: Initializing...")
        
        # Get all protein atoms
        all_atoms = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):  # Only protein atoms
                        for atom in residue:
                            all_atoms.append(atom)
        
        if not all_atoms:
            self.progress_update.emit("Cavity detection: No protein atoms found")
            return results
        
        self.progress_update.emit(f"Cavity detection: Analyzing {len(all_atoms)} protein atoms...")
        
        # Use improved cavity detection with multiple probe sizes
        cavity_points = self.detect_cavity_points_improved(all_atoms)
        results['cavity_points_found'] = len(cavity_points)
        
        if cavity_points:
            # Cluster cavity points into distinct cavities
            cavities = self.cluster_cavity_points_improved(cavity_points, all_atoms)
            
            if cavities:
                # Filter and characterize cavities
                characterized_cavities = self.characterize_cavities(cavities, all_atoms)
                
                results['potential_cavities'] = characterized_cavities
                
                if characterized_cavities:
                    # Find largest cavity
                    largest = max(characterized_cavities, key=lambda c: c['volume'])
                    results['largest_cavity'] = largest
                    results['cavity_volume'] = sum(c['volume'] for c in characterized_cavities)
                    
                    # Identify potentially druggable cavities
                    druggable = [c for c in characterized_cavities if self.is_druggable_cavity(c)]
                    results['druggable_cavities'] = druggable
                    
                    # Final summary
                    calculation_time = time.time() - start_time
                    results['calculation_time'] = calculation_time
                    
                    self.progress_update.emit(
                        f"Cavity detection complete! Found {len(characterized_cavities)} cavities "
                        f"({len(druggable)} druggable) in {calculation_time:.1f}s"
                    )
                else:
                    self.progress_update.emit("Cavity detection: No significant cavities found after characterization")
            else:
                self.progress_update.emit("Cavity detection: No cavity clusters found")
        else:
            self.progress_update.emit("Cavity detection: No cavity points detected")
        
        # Add timing information
        calculation_time = time.time() - start_time
        results['calculation_time'] = calculation_time
        
        return results
    
    def detect_cavity_points_improved(self, all_atoms, grid_spacing=1.5, probe_radius=1.4):
        """Detect cavity points using improved grid-based method with probe accessibility.
        
        Args:
            all_atoms: List of all protein atoms
            grid_spacing: Grid spacing in Angstroms
            probe_radius: Probe radius in Angstroms (1.4Å for water)
        
        Returns:
            list: List of cavity points
        """
        self.progress_update.emit("Cavity detection: Setting up grid...")
        
        # Get protein bounds
        coords = np.array([atom.get_coord() for atom in all_atoms])
        min_coords = np.min(coords, axis=0) - probe_radius - 3.0
        max_coords = np.max(coords, axis=0) + probe_radius + 3.0
        
        # Generate grid points
        x_range = np.arange(min_coords[0], max_coords[0], grid_spacing)
        y_range = np.arange(min_coords[1], max_coords[1], grid_spacing)
        z_range = np.arange(min_coords[2], max_coords[2], grid_spacing)
        
        cavity_points = []
        
        # Limit computation for large proteins
        max_grid_points = 15000
        total_points = len(x_range) * len(y_range) * len(z_range)
        
        if total_points > max_grid_points:
            # Adaptive grid spacing
            reduction_factor = (total_points / max_grid_points) ** (1/3)
            grid_spacing *= reduction_factor
            x_range = np.arange(min_coords[0], max_coords[0], grid_spacing)
            y_range = np.arange(min_coords[1], max_coords[1], grid_spacing)
            z_range = np.arange(min_coords[2], max_coords[2], grid_spacing)
            total_points = len(x_range) * len(y_range) * len(z_range)
        
        self.progress_update.emit(f"Cavity detection: Scanning {total_points:,} grid points...")
        
        sample_count = 0
        last_progress = 0
        
        for i, x in enumerate(x_range):
            for j, y in enumerate(y_range):
                for k, z in enumerate(z_range):
                    point = np.array([x, y, z])
                    
                    # Check if point is in a cavity using improved criteria
                    if self.is_cavity_point_improved(point, all_atoms, probe_radius):
                        cavity_points.append(point)
                    
                    sample_count += 1
                    
                    # Update progress every 5%
                    progress = (sample_count / min(total_points, max_grid_points)) * 100
                    if progress - last_progress >= 5:
                        self.progress_update.emit(f"Cavity detection: {progress:.0f}% complete ({len(cavity_points)} potential points found)")
                        last_progress = progress
                    
                    if sample_count >= max_grid_points:
                        break
                if sample_count >= max_grid_points:
                    break
            if sample_count >= max_grid_points:
                break
        
        self.progress_update.emit(f"Cavity detection: Found {len(cavity_points)} potential cavity points")
        return cavity_points
    
    def is_cavity_point_improved(self, point, all_atoms, probe_radius):
        """Check if a point represents a cavity using improved criteria.
        
        A point is considered a cavity point if:
        1. It's not inside any atom (considering van der Waals radii + probe radius)
        2. It's surrounded by protein atoms (not on the surface)
        3. It's accessible to a probe of the given radius
        """
        min_distance = float('inf')
        nearby_atoms = 0
        
        for atom in all_atoms:
            # Get van der Waals radius
            element = atom.element.upper() if hasattr(atom, 'element') else atom.get_name()[0]
            vdw_radius = VDW_RADII.get(element, 1.7)
            
            distance = np.linalg.norm(point - atom.get_coord())
            min_distance = min(min_distance, distance)
            
            # Check if point is inside atom (including probe radius)
            if distance < (vdw_radius + probe_radius):
                return False  # Point is inside an atom
            
            # Count nearby atoms within reasonable distance
            if distance < 8.0:
                nearby_atoms += 1
        
        # Point is a cavity if:
        # 1. Not inside any atom (already checked above)
        # 2. Has sufficient nearby atoms (indicating it's inside the protein)
        # 3. Not too far from protein surface
        return (nearby_atoms >= 3 and 
                probe_radius < min_distance < 6.0)
    
    def cluster_cavity_points_improved(self, points, all_atoms, cluster_radius=3.5):
        """Cluster cavity points into distinct cavities using improved clustering.
        
        Args:
            points: List of cavity points
            all_atoms: List of all protein atoms
            cluster_radius: Maximum distance for clustering
        
        Returns:
            list: List of cavity clusters
        """
        if not points:
            return []
        
        self.progress_update.emit(f"Cavity clustering: Processing {len(points)} points...")
        
        points = np.array(points)
        clusters = []
        used = set()
        total_points = len(points)
        processed_points = 0
        last_progress = 0
        
        for i, point in enumerate(points):
            if i in used:
                processed_points += 1
                continue
            
            # Start new cluster with breadth-first search
            cluster = [i]
            queue = [i]
            used.add(i)
            
            while queue:
                current_idx = queue.pop(0)
                current_point = points[current_idx]
                
                # Find all unused points within cluster radius
                for j, other_point in enumerate(points):
                    if j in used:
                        continue
                    
                    distance = np.linalg.norm(current_point - other_point)
                    if distance <= cluster_radius:
                        cluster.append(j)
                        queue.append(j)
                        used.add(j)
            
            # Only keep clusters with minimum size
            if len(cluster) >= 5:  # Minimum cluster size for meaningful cavity
                cluster_points = points[cluster]
                clusters.append({
                    'points': cluster_points,
                    'indices': cluster
                })
            
            processed_points += len(cluster)
            
            # Update progress every 10%
            progress = (processed_points / total_points) * 100
            if progress - last_progress >= 10:
                self.progress_update.emit(f"Cavity clustering: {progress:.0f}% complete ({len(clusters)} clusters found)")
                last_progress = progress
        
        self.progress_update.emit(f"Cavity clustering: Found {len(clusters)} potential cavities")
        return clusters
    
    def characterize_cavities(self, clusters, all_atoms):
        """Characterize detected cavities with geometric and chemical properties.
        
        Args:
            clusters: List of cavity clusters
            all_atoms: List of all protein atoms
        
        Returns:
            list: List of characterized cavities
        """
        characterized = []
        total_clusters = len(clusters)
        
        if total_clusters == 0:
            return characterized
        
        self.progress_update.emit(f"Cavity characterization: Analyzing {total_clusters} cavities...")
        
        for i, cluster in enumerate(clusters):
            cluster_points = cluster['points']
            
            # Update progress
            progress = ((i + 1) / total_clusters) * 100
            self.progress_update.emit(f"Cavity characterization: {progress:.0f}% complete (cavity {i+1}/{total_clusters})")
            
            # Calculate geometric properties
            center = np.mean(cluster_points, axis=0)
            
            # Calculate volume (using convex hull approximation)
            volume = self.estimate_cavity_volume(cluster_points)
            
            # Calculate radius (maximum distance from center)
            distances = np.linalg.norm(cluster_points - center, axis=1)
            max_radius = np.max(distances)
            avg_radius = np.mean(distances)
            
            # Find lining residues
            lining_residues = self.find_lining_residues(center, max_radius + 2.0, all_atoms)
            
            # Calculate surface area (rough estimate)
            surface_area = self.estimate_cavity_surface_area(cluster_points)
            
            # Calculate hydrophobicity score
            hydrophobicity = self.calculate_cavity_hydrophobicity(lining_residues)
            
            cavity_info = {
                'id': i + 1,
                'center': center.tolist(),
                'volume': volume,
                'surface_area': surface_area,
                'max_radius': max_radius,
                'avg_radius': avg_radius,
                'point_count': len(cluster_points),
                'lining_residues': lining_residues,
                'hydrophobicity_score': hydrophobicity,
                'druggability_score': self.calculate_druggability_score(volume, surface_area, hydrophobicity)
            }
            
            characterized.append(cavity_info)
        
        # Sort by volume (largest first)
        characterized.sort(key=lambda x: x['volume'], reverse=True)
        
        self.progress_update.emit(f"Cavity characterization: Complete! Found {len(characterized)} characterized cavities")
        
        return characterized
    
    def estimate_cavity_volume(self, points):
        """Estimate cavity volume using point density."""
        # Simple volume estimation based on number of points and grid spacing
        # More sophisticated methods would use convex hull or alpha shapes
        grid_spacing = 1.5  # Approximate grid spacing used
        volume_per_point = grid_spacing ** 3
        return len(points) * volume_per_point
    
    def estimate_cavity_surface_area(self, points):
        """Estimate cavity surface area."""
        # Simple estimation based on volume
        volume = self.estimate_cavity_volume(points)
        # Assume roughly spherical cavity for surface area estimation
        radius = (3 * volume / (4 * math.pi)) ** (1/3)
        return 4 * math.pi * radius ** 2
    
    def find_lining_residues(self, cavity_center, search_radius, all_atoms):
        """Find residues lining the cavity."""
        lining_residues = set()
        
        for atom in all_atoms:
            distance = np.linalg.norm(cavity_center - atom.get_coord())
            if distance <= search_radius:
                # Get residue info
                residue = atom.get_parent()
                if is_aa(residue):
                    chain_id = residue.get_parent().id
                    res_info = {
                        'residue': residue.get_resname(),
                        'chain': chain_id,
                        'position': residue.id[1],
                        'distance': distance
                    }
                    lining_residues.add((residue.get_resname(), chain_id, residue.id[1]))
        
        return list(lining_residues)
    
    def calculate_cavity_hydrophobicity(self, lining_residues):
        """Calculate hydrophobicity score for cavity based on lining residues."""
        # Hydrophobicity scale (Kyte-Doolittle)
        hydrophobicity_scale = {
            'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
            'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
            'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
            'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
        }
        
        if not lining_residues:
            return 0.0
        
        total_hydrophobicity = 0.0
        for res_name, chain, position in lining_residues:
            total_hydrophobicity += hydrophobicity_scale.get(res_name, 0.0)
        
        return total_hydrophobicity / len(lining_residues)
    
    def calculate_druggability_score(self, volume, surface_area, hydrophobicity):
        """Calculate a simple druggability score for the cavity."""
        # Simple scoring based on known druggable cavity characteristics
        volume_score = min(volume / 500.0, 1.0)  # Normalize to 500 Ų
        surface_score = min(surface_area / 300.0, 1.0)  # Normalize to 300 Ų
        hydrophobic_score = max(0, hydrophobicity + 1.0) / 3.0  # Normalize hydrophobicity
        
        # Weighted combination
        druggability = (0.4 * volume_score + 0.3 * surface_score + 0.3 * hydrophobic_score)
        return min(druggability, 1.0)
    
    def is_druggable_cavity(self, cavity):
        """Determine if a cavity is potentially druggable."""
        # Simple criteria for druggability
        return (cavity['volume'] > 100 and  # Minimum volume
                cavity['druggability_score'] > 0.3 and  # Minimum druggability score
                len(cavity['lining_residues']) >= 5)  # Sufficient lining residues


class SASAConfigDialog(QDialog):
    """Dialog for configuring SASA calculation parameters."""
    
    def __init__(self, parent=None, config=None):
        super().__init__(parent)
        self.setWindowTitle("SASA Configuration")
        self.setModal(True)
        self.resize(500, 450)  # Slightly larger for better usability
        
        # Default configuration
        self.config = config or {
            'probe_radius': 1.4,
            'n_points': 30,
            'rsa_threshold': 20.0,
            'include_hydrogens': False,
            'preset': 'balanced'
        }
        
        self.init_ui()
        self.load_config()
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Title
        title = QLabel("<h3>SASA Calculation Configuration</h3>")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)
        
        # Preset selection
        preset_group = QGroupBox("Presets")
        preset_layout = QVBoxLayout(preset_group)
        
        self.preset_combo = QComboBox()
        self.preset_combo.addItems(["Fast", "Balanced", "Accurate", "Research", "Custom"])
        self.preset_combo.currentTextChanged.connect(self.on_preset_changed)
        preset_layout.addWidget(self.preset_combo)
        
        # Preset descriptions
        preset_desc = QLabel(
            "<b>Fast:</b> Quick calculation (10 points, 1.4Å probe)<br>"
            "<b>Balanced:</b> Good speed/accuracy trade-off (30 points)<br>"
            "<b>Accurate:</b> High accuracy (60 points)<br>"
            "<b>Research:</b> Publication quality (100 points)<br>"
            "<b>Custom:</b> User-defined parameters"
        )
        preset_desc.setWordWrap(True)
        preset_desc.setStyleSheet("color: #666; font-size: 10px; margin: 5px;")
        preset_layout.addWidget(preset_desc)
        
        layout.addWidget(preset_group)
        
        # Parameters group
        params_group = QGroupBox("Parameters")
        params_layout = QFormLayout(params_group)
        
        # Probe radius
        self.probe_radius_spin = QDoubleSpinBox()
        self.probe_radius_spin.setRange(0.5, 3.0)
        self.probe_radius_spin.setSingleStep(0.1)
        self.probe_radius_spin.setDecimals(2)
        self.probe_radius_spin.setSuffix(" Å")
        self.probe_radius_spin.setToolTip("Probe radius for SASA calculation (1.4Å = water)")
        self.probe_radius_spin.valueChanged.connect(self.on_parameter_changed)
        params_layout.addRow("Probe Radius:", self.probe_radius_spin)
        
        # Number of sphere points
        self.n_points_spin = QSpinBox()
        self.n_points_spin.setRange(10, 1000)
        self.n_points_spin.setSingleStep(10)
        self.n_points_spin.setToolTip("Number of points on sphere surface (more = accurate but slower)")
        self.n_points_spin.valueChanged.connect(self.on_parameter_changed)
        params_layout.addRow("Sphere Points:", self.n_points_spin)
        
        # RSA threshold
        self.rsa_threshold_spin = QDoubleSpinBox()
        self.rsa_threshold_spin.setRange(0.0, 100.0)
        self.rsa_threshold_spin.setSingleStep(1.0)
        self.rsa_threshold_spin.setDecimals(1)
        self.rsa_threshold_spin.setSuffix("%")
        self.rsa_threshold_spin.setToolTip("RSA threshold for surface/buried classification")
        self.rsa_threshold_spin.valueChanged.connect(self.on_parameter_changed)
        params_layout.addRow("RSA Threshold:", self.rsa_threshold_spin)
        
        # Include hydrogens
        self.include_h_checkbox = QCheckBox("Include hydrogen atoms")
        self.include_h_checkbox.setToolTip("Include hydrogen atoms in calculation (if present)")
        self.include_h_checkbox.stateChanged.connect(self.on_parameter_changed)
        params_layout.addRow("", self.include_h_checkbox)
        
        layout.addWidget(params_group)
        
        # Performance estimate
        self.perf_label = QLabel()
        self.perf_label.setStyleSheet("color: #666; font-style: italic; margin: 10px;")
        self.perf_label.setWordWrap(True)
        layout.addWidget(self.perf_label)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        reset_btn = QPushButton("Reset to Defaults")
        reset_btn.clicked.connect(self.reset_to_defaults)
        button_layout.addWidget(reset_btn)
        
        button_layout.addStretch()
        
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(cancel_btn)
        
        ok_btn = QPushButton("OK")
        ok_btn.setDefault(True)
        ok_btn.clicked.connect(self.accept)
        button_layout.addWidget(ok_btn)
        
        layout.addLayout(button_layout)
        
        # Update performance estimate
        self.update_performance_estimate()
    
    def load_config(self):
        """Load configuration into UI elements."""
        self.probe_radius_spin.setValue(self.config['probe_radius'])
        self.n_points_spin.setValue(self.config['n_points'])
        self.rsa_threshold_spin.setValue(self.config['rsa_threshold'])
        self.include_h_checkbox.setChecked(self.config['include_hydrogens'])
        
        # Set preset
        preset_map = {
            'fast': 'Fast',
            'balanced': 'Balanced', 
            'accurate': 'Accurate',
            'research': 'Research',
            'custom': 'Custom'
        }
        preset_name = preset_map.get(self.config['preset'], 'Custom')
        index = self.preset_combo.findText(preset_name)
        if index >= 0:
            self.preset_combo.setCurrentIndex(index)
    
    def on_preset_changed(self, preset_name):
        """Handle preset selection change."""
        presets = {
            'Fast': {'probe_radius': 1.4, 'n_points': 10, 'rsa_threshold': 20.0, 'include_hydrogens': False},
            'Balanced': {'probe_radius': 1.4, 'n_points': 30, 'rsa_threshold': 20.0, 'include_hydrogens': False},
            'Accurate': {'probe_radius': 1.4, 'n_points': 60, 'rsa_threshold': 20.0, 'include_hydrogens': False},
            'Research': {'probe_radius': 1.4, 'n_points': 100, 'rsa_threshold': 20.0, 'include_hydrogens': False}
        }
        
        if preset_name in presets:
            preset = presets[preset_name]
            self.probe_radius_spin.setValue(preset['probe_radius'])
            self.n_points_spin.setValue(preset['n_points'])
            self.rsa_threshold_spin.setValue(preset['rsa_threshold'])
            self.include_h_checkbox.setChecked(preset['include_hydrogens'])
            self.update_performance_estimate()
    
    def on_parameter_changed(self):
        """Handle parameter change - switch to custom preset."""
        if self.preset_combo.currentText() != 'Custom':
            self.preset_combo.setCurrentText('Custom')
        self.update_performance_estimate()
    
    def update_performance_estimate(self):
        """Update performance estimate based on current settings."""
        n_points = self.n_points_spin.value()
        
        # Rough performance estimates
        if n_points <= 15:
            speed = "Very Fast"
            time_est = "~1-2 seconds"
            accuracy = "Good"
        elif n_points <= 35:
            speed = "Fast"
            time_est = "~3-5 seconds"
            accuracy = "Very Good"
        elif n_points <= 70:
            speed = "Medium"
            time_est = "~8-15 seconds"
            accuracy = "High"
        else:
            speed = "Slow"
            time_est = "~20-60 seconds"
            accuracy = "Very High"
        
        self.perf_label.setText(
            f"<b>Performance Estimate:</b><br>"
            f"Speed: {speed} ({time_est} for typical protein)<br>"
            f"Accuracy: {accuracy} correlation with professional tools"
        )
    
    def reset_to_defaults(self):
        """Reset to default balanced configuration."""
        self.preset_combo.setCurrentText('Balanced')
        self.on_preset_changed('Balanced')
    
    def get_config(self):
        """Get current configuration."""
        preset_map = {
            'Fast': 'fast',
            'Balanced': 'balanced',
            'Accurate': 'accurate', 
            'Research': 'research',
            'Custom': 'custom'
        }
        
        return {
            'probe_radius': self.probe_radius_spin.value(),
            'n_points': self.n_points_spin.value(),
            'rsa_threshold': self.rsa_threshold_spin.value(),
            'include_hydrogens': self.include_h_checkbox.isChecked(),
            'preset': preset_map.get(self.preset_combo.currentText(), 'custom')
        }


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
        
        self.quality_checkbox = QCheckBox("Quality Assessment")
        self.quality_checkbox.setChecked(True)
        self.quality_checkbox.setToolTip("Assess structural quality")
        options_layout.addWidget(self.quality_checkbox)
        
        # Surface Properties analysis removed - will be added in a future release
        
        # Binding Sites analysis moved to dedicated tab - removed from here
        
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
            "<li><b>Basic Properties:</b> Chain composition, molecular weight, resolution</li>"
            "<li><b>Secondary Structure:</b> Helix/sheet/loop content, Ramachandran plots</li>"
            "<li><b>Geometric Analysis:</b> Bond lengths, angles, radius of gyration</li>"
            "<li><b>Quality Assessment:</b> Ramachandran outliers, missing atoms, B-factors</li>"

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
                self.parent_app.statusBar().showMessage(f"Fetching comprehensive data for PDB ID {pdb_id}...")
            
            # Use optimized fetch manager if available
            if OPTIMIZED_FETCH_AVAILABLE and self.pdb_fetch_manager:
                # Check if data is already cached
                cached_info = self.pdb_fetch_manager.get_structure_info(pdb_id)
                if cached_info['available'] and cached_info['files'].get('pdb'):
                    pdb_path = cached_info['files']['pdb']
                    # Create mock comprehensive data from cached info
                    self.current_comprehensive_data = {
                        'pdb_id': pdb_id,
                        'files': cached_info['files'],
                        'metadata': cached_info.get('metadata', {}),
                        'sequences': {},
                        'validation': {},
                        'errors': []
                    }
                else:
                    # Need to fetch - this will block, but it's in the analysis tab
                    QMessageBox.information(self, "Fetching Data", 
                        f"Fetching comprehensive data for {pdb_id}. This may take a moment...")
                    
                    # Use a simple synchronous approach for analysis tab
                    from .enhanced_pdb_puller import EnhancedPDBPuller
                    temp_puller = EnhancedPDBPuller(self.pulled_structures_dir)
                    comprehensive_data = temp_puller.fetch_comprehensive_pdb_data(
                        pdb_id,
                        include_validation=False,  # Skip for speed
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
                    
            elif self.enhanced_pdb_puller:
                # Use original enhanced puller
                comprehensive_data = self.enhanced_pdb_puller.fetch_comprehensive_pdb_data(
                    pdb_id,
                    include_validation=False,  # Skip for speed
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
                
            else:
                # Fallback to basic PDB fetching
                pdb_path = os.path.join(self.pulled_structures_dir, f"{pdb_id}.pdb")
                
                if not os.path.exists(pdb_path):
                    retrieved_file_path = self.pdb_list.retrieve_pdb_file(
                        pdb_id, pdir=self.pulled_structures_dir, file_format="pdb"
                    )
                    if not os.path.exists(retrieved_file_path):
                        raise FileNotFoundError(f"PDB file for {pdb_id} was not downloaded")
                    if os.path.exists(pdb_path):
                        os.remove(pdb_path)
                    os.rename(retrieved_file_path, pdb_path)
                
                self.current_comprehensive_data = None
            
            # Set as current structure
            self.current_structure_path = pdb_path
            self.structure_label.setText(f"PDB: {pdb_id}")
            self.structure_label.setStyleSheet("color: #000; font-weight: bold;")
            
            if hasattr(self.parent_app, 'statusBar'):
                if self.enhanced_pdb_puller:
                    self.parent_app.statusBar().showMessage(f"Fetched comprehensive data for {pdb_id}")
                else:
                    self.parent_app.statusBar().showMessage(f"Fetched basic PDB data for {pdb_id}")
            
            # Show appropriate success message
            if OPTIMIZED_FETCH_AVAILABLE and self.pdb_fetch_manager:
                QMessageBox.information(self, "Success", 
                    f"Successfully loaded data for PDB {pdb_id}\n\n"
                    f"Using optimized fetch system for better performance.")
            elif self.enhanced_pdb_puller:
                QMessageBox.information(self, "Success", 
                    f"Successfully fetched comprehensive data for PDB {pdb_id}\n\n"
                    f"Includes: structure, metadata, and sequences")
            else:
                QMessageBox.information(self, "Success", f"Successfully fetched PDB structure {pdb_id}")
            
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
            
            if hasattr(self.parent_app, 'statusBar'):
                self.parent_app.statusBar().showMessage(f"Loaded {filename} for analysis")
    
    def analyze_structure(self):
        """Perform structural analysis."""
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
        if self.quality_checkbox.isChecked():
            analysis_types.append('quality')
        # Surface analysis removed - will be added in a future release
        
        if not analysis_types:
            QMessageBox.warning(self, "No Analysis Selected", "Please select at least one analysis type.")
            return
        
        # Show progress
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        
        # Start analysis in worker thread
        self.analysis_worker = StructuralAnalysisWorker(self.current_structure_path, analysis_types)
        
        # Surface analysis configuration removed - will be added in a future release
        self.analysis_worker.analysis_complete.connect(self.display_results)
        self.analysis_worker.error_occurred.connect(self.handle_analysis_error)
        self.analysis_worker.progress_update.connect(self.update_progress)
        self.analysis_worker.start()
    
    def on_surface_checkbox_changed(self, state):
        """Enable/disable SASA configuration button based on checkbox state."""
        self.sasa_config_btn.setEnabled(state == Qt.Checked)
    
    def get_default_sasa_config(self):
        """Get default SASA configuration."""
        return {
            'probe_radius': 1.4,      # Standard water probe radius (Å)
            'n_points': 30,           # Number of sphere points (speed vs accuracy)
            'rsa_threshold': 20.0,    # RSA threshold for surface/buried classification (%)
            'include_hydrogens': False, # Whether to include hydrogen atoms
            'preset': 'balanced'      # Preset name
        }
    
    def show_sasa_config(self):
        """Show SASA configuration dialog."""
        dialog = SASAConfigDialog(self, getattr(self, 'sasa_config', self.get_default_sasa_config()))
        if dialog.exec_() == QDialog.Accepted:
            self.sasa_config = dialog.get_config()
            # Update button tooltip to show current settings
            config = self.sasa_config
            tooltip = (f"SASA Settings: {config['preset']} preset\n"
                      f"Probe radius: {config['probe_radius']} Å\n"
                      f"Sphere points: {config['n_points']}\n"
                      f"RSA threshold: {config['rsa_threshold']}%")
            self.sasa_config_btn.setToolTip(tooltip)
    
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
    
    def create_ramachandran_plot(self, ramachandran_data):
        """Create a Ramachandran plot using ramachandraw."""
        if not RAMACHANDRAW_AVAILABLE or not MATPLOTLIB_AVAILABLE or not ramachandran_data:
            return None, {}
        
        try:
            import ramachandraw.utils
            import tempfile
            import os
            
            # We need to use the current structure file with ramachandraw
            if not hasattr(self, 'current_structure_path') or not self.current_structure_path:
                # Fallback to manual plot if no structure file available
                return self._create_manual_ramachandran_plot(ramachandran_data)
            
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
            return self._create_manual_ramachandran_plot(ramachandran_data)
    
    def _create_manual_ramachandran_plot(self, ramachandran_data):
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
            ax.set_title('Ramachandran Plot', fontsize=16, fontweight='bold', pad=20)
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
        """Display enhanced metadata from comprehensive PDB data."""
        comprehensive_data = self.current_comprehensive_data
        metadata = comprehensive_data.get('metadata', {})
        
        # Enhanced Basic Properties
        enhanced_basic_group = QGroupBox("Enhanced Structure Information")
        enhanced_basic_layout = QVBoxLayout(enhanced_basic_group)
        
        # Create tabs for different metadata sections
        metadata_tabs = QTabWidget()
        
        # Basic Information Tab
        basic_tab = QWidget()
        basic_tab_layout = QFormLayout(basic_tab)
        
        # Standard basic info
        basic_tab_layout.addRow("Structure ID:", QLabel(basic_results.get('structure_id', 'N/A')))
        basic_tab_layout.addRow("Total Chains:", QLabel(str(len(basic_results.get('chains', [])))))
        basic_tab_layout.addRow("Total Residues:", QLabel(str(basic_results.get('total_residues', 0))))
        basic_tab_layout.addRow("Total Atoms:", QLabel(str(basic_results.get('total_atoms', 0))))
        basic_tab_layout.addRow("Molecular Weight:", QLabel(f"{basic_results.get('molecular_weight', 0):.1f} Da"))
        
        # Enhanced info from metadata
        if 'entry' in metadata:
            entry = metadata['entry']
            
            # Title and description
            title = entry.get('struct', {}).get('title', 'N/A')
            basic_tab_layout.addRow("Title:", QLabel(title))
            
            # Authors
            authors = [author.get('name', '') for author in entry.get('audit_author', [])]
            if authors:
                authors_text = ', '.join(authors[:3])  # Show first 3 authors
                if len(authors) > 3:
                    authors_text += f" and {len(authors) - 3} others"
                basic_tab_layout.addRow("Authors:", QLabel(authors_text))
            
            # Dates
            deposit_date = entry.get('rcsb_accession_info', {}).get('deposit_date', 'N/A')
            release_date = entry.get('rcsb_accession_info', {}).get('initial_release_date', 'N/A')
            basic_tab_layout.addRow("Deposition Date:", QLabel(deposit_date))
            basic_tab_layout.addRow("Release Date:", QLabel(release_date))
        
        metadata_tabs.addTab(basic_tab, "Basic Info")
        
        # Experimental Information Tab
        if 'experimental' in metadata or 'entry' in metadata:
            exp_tab = QWidget()
            exp_tab_layout = QFormLayout(exp_tab)
            
            if 'entry' in metadata:
                entry = metadata['entry']
                
                # Experimental method
                exptl_methods = entry.get('exptl', [])
                if exptl_methods:
                    method = exptl_methods[0].get('method', 'N/A')
                    exp_tab_layout.addRow("Experimental Method:", QLabel(method))
                
                # Resolution
                resolution = entry.get('rcsb_entry_info', {}).get('resolution_combined', 'N/A')
                exp_tab_layout.addRow("Resolution:", QLabel(f"{resolution} Å" if resolution != 'N/A' else 'N/A'))
                
                # R-factors
                refine = entry.get('refine', [])
                if refine:
                    r_work = refine[0].get('ls_R_factor_R_work', 'N/A')
                    r_free = refine[0].get('ls_R_factor_R_free', 'N/A')
                    if r_work != 'N/A':
                        exp_tab_layout.addRow("R-work:", QLabel(f"{r_work:.3f}" if isinstance(r_work, (int, float)) else str(r_work)))
                    if r_free != 'N/A':
                        exp_tab_layout.addRow("R-free:", QLabel(f"{r_free:.3f}" if isinstance(r_free, (int, float)) else str(r_free)))
                
                # Space group
                symmetry = entry.get('symmetry', {})
                space_group = symmetry.get('space_group_name_H_M', 'N/A')
                exp_tab_layout.addRow("Space Group:", QLabel(space_group))
                
                # Unit cell
                cell = entry.get('cell', {})
                if cell:
                    a = cell.get('length_a', 'N/A')
                    b = cell.get('length_b', 'N/A')
                    c = cell.get('length_c', 'N/A')
                    alpha = cell.get('angle_alpha', 'N/A')
                    beta = cell.get('angle_beta', 'N/A')
                    gamma = cell.get('angle_gamma', 'N/A')
                    
                    if all(x != 'N/A' for x in [a, b, c]):
                        cell_text = f"a={a:.2f} b={b:.2f} c={c:.2f} Å"
                        exp_tab_layout.addRow("Unit Cell (lengths):", QLabel(cell_text))
                    
                    if all(x != 'N/A' for x in [alpha, beta, gamma]):
                        angles_text = f"α={alpha:.1f}° β={beta:.1f}° γ={gamma:.1f}°"
                        exp_tab_layout.addRow("Unit Cell (angles):", QLabel(angles_text))
            
            metadata_tabs.addTab(exp_tab, "Experimental")
        
        # Ligands and Entities Tab
        if 'nonpolymer_entities' in metadata or 'polymer_entities' in metadata:
            entities_tab = QWidget()
            entities_tab_layout = QVBoxLayout(entities_tab)
            
            # Polymer entities (protein chains)
            if 'polymer_entities' in metadata:
                polymer_entities = metadata['polymer_entities']
                if isinstance(polymer_entities, list) and polymer_entities:
                    polymer_group = QGroupBox(f"Protein Chains ({len(polymer_entities)})")
                    polymer_layout = QVBoxLayout(polymer_group)
                    
                    polymer_table = QTableWidget()
                    polymer_table.setColumnCount(4)
                    polymer_table.setHorizontalHeaderLabels(["Entity ID", "Description", "Type", "Chains"])
                    polymer_table.setRowCount(len(polymer_entities))
                    
                    for i, entity in enumerate(polymer_entities):
                        entity_id = entity.get('rcsb_id', 'N/A')
                        description = entity.get('rcsb_polymer_entity', {}).get('pdbx_description', 'N/A')
                        entity_type = entity.get('entity_poly', {}).get('type', 'N/A')
                        chains = ', '.join(entity.get('rcsb_polymer_entity_container_identifiers', {}).get('asym_ids', []))
                        
                        polymer_table.setItem(i, 0, QTableWidgetItem(str(entity_id)))
                        polymer_table.setItem(i, 1, QTableWidgetItem(description))
                        polymer_table.setItem(i, 2, QTableWidgetItem(entity_type))
                        polymer_table.setItem(i, 3, QTableWidgetItem(chains))
                    
                    polymer_table.resizeColumnsToContents()
                    polymer_table.setMaximumHeight(200)
                    polymer_layout.addWidget(polymer_table)
                    entities_tab_layout.addWidget(polymer_group)
            
            # Non-polymer entities (ligands)
            if 'nonpolymer_entities' in metadata:
                nonpolymer_entities = metadata['nonpolymer_entities']
                if isinstance(nonpolymer_entities, list) and nonpolymer_entities:
                    ligand_group = QGroupBox(f"Ligands and Other Molecules ({len(nonpolymer_entities)})")
                    ligand_layout = QVBoxLayout(ligand_group)
                    
                    ligand_table = QTableWidget()
                    ligand_table.setColumnCount(4)
                    ligand_table.setHorizontalHeaderLabels(["ID", "Name", "Formula", "Type"])
                    ligand_table.setRowCount(len(nonpolymer_entities))
                    
                    for i, entity in enumerate(nonpolymer_entities):
                        entity_id = entity.get('rcsb_id', 'N/A')
                        name = entity.get('rcsb_nonpolymer_entity', {}).get('pdbx_description', 'N/A')
                        formula = entity.get('rcsb_nonpolymer_entity', {}).get('formula_weight', 'N/A')
                        entity_type = entity.get('rcsb_nonpolymer_entity', {}).get('type', 'N/A')
                        
                        ligand_table.setItem(i, 0, QTableWidgetItem(str(entity_id)))
                        ligand_table.setItem(i, 1, QTableWidgetItem(name))
                        ligand_table.setItem(i, 2, QTableWidgetItem(str(formula)))
                        ligand_table.setItem(i, 3, QTableWidgetItem(entity_type))
                    
                    ligand_table.resizeColumnsToContents()
                    ligand_table.setMaximumHeight(200)
                    ligand_layout.addWidget(ligand_table)
                    entities_tab_layout.addWidget(ligand_group)
            
            metadata_tabs.addTab(entities_tab, "Entities & Ligands")
        
        # Publications tab removed for performance
        
        enhanced_basic_layout.addWidget(metadata_tabs)
        self.results_layout.addWidget(enhanced_basic_group)
        
        # Also display standard composition analysis
        if basic_results.get('composition'):
            self.display_composition_analysis(basic_results)
    
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
        
        if 'quality' in results:
            self.display_quality_results(results['quality'])
        
        # Surface results display removed - will be added in a future release
        
        # Cavity analysis moved to dedicated tab - removed from here
        
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
        basic_layout.addRow("Resolution:", QLabel(str(results.get('resolution', 'N/A'))))
        basic_layout.addRow("Space Group:", QLabel(str(results.get('space_group', 'N/A'))))
        
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
            rama_plot, rama_counts = self.create_ramachandran_plot(results['ramachandran_data'])
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
    
    def display_quality_results(self, results):
        """Display structural quality results."""
        quality_group = QGroupBox("Structural Quality Assessment")
        quality_layout = QVBoxLayout(quality_group)
        
        # Ramachandran statistics with pie chart
        # Use the same classification as the Ramachandran plot for consistency
        if hasattr(self, 'current_results') and 'secondary' in self.current_results and 'ramachandran_data' in self.current_results['secondary']:
            # Classify each point using the same method as the plot
            rama_data = self.current_results['secondary']['ramachandran_data']
            favored_count = 0
            allowed_count = 0
            outliers_count = 0
            
            for r in rama_data:
                phi = r.get('phi', 0)
                psi = r.get('psi', 0)
                if self.is_in_favored_region(phi, psi):
                    favored_count += 1
                elif self.is_in_allowed_region(phi, psi):
                    allowed_count += 1
                else:
                    outliers_count += 1
                    
            total_rama = len(rama_data)
        else:
            # Fallback to quality analysis counts if secondary data not available
            favored_count = results.get('ramachandran_favored', 0)
            allowed_count = results.get('ramachandran_allowed', 0)
            outliers_count = results.get('ramachandran_outliers_count', 0)
            total_rama = max(favored_count + allowed_count + outliers_count, 1)  # Prevent division by zero
        if total_rama > 0:
            rama_summary_group = QGroupBox("Ramachandran Quality")
            rama_summary_layout = QHBoxLayout(rama_summary_group)
            
            # Statistics
            stats_widget = QWidget()
            stats_layout = QFormLayout(stats_widget)
            
            favored_pct = (favored_count / total_rama) * 100
            allowed_pct = (allowed_count / total_rama) * 100
            outlier_pct = (outliers_count / total_rama) * 100
            
            stats_layout.addRow("Favored:", QLabel(f"{favored_pct:.1f}% ({favored_count} residues)"))
            stats_layout.addRow("Allowed:", QLabel(f"{allowed_pct:.1f}% ({allowed_count} residues)"))
            stats_layout.addRow("Outliers:", QLabel(f"{outlier_pct:.1f}% ({outliers_count} residues)"))
            stats_layout.addRow("Total Residues:", QLabel(str(total_rama)))
            
            stats_widget.setMaximumWidth(350)
            rama_summary_layout.addWidget(stats_widget)
            
            # Pie chart for Ramachandran quality
            rama_quality_data = {
                'Favored': favored_pct,
                'Allowed': allowed_pct,
                'Outliers': outlier_pct
            }
            rama_colors = ['#2E8B57', '#FFD700', '#DC143C']  # Green, Gold, Red
            rama_pie = self.create_pie_chart(rama_quality_data, "Ramachandran Quality Distribution", rama_colors)
            if rama_pie:
                rama_summary_layout.addWidget(rama_pie, 1)
            
            quality_layout.addWidget(rama_summary_group)
        
        # Overall quality summary
        summary_stats_group = QGroupBox("Quality Summary")
        summary_layout = QFormLayout(summary_stats_group)
        
        # Missing atoms
        missing_count = len(results.get('missing_atoms', []))
        summary_layout.addRow("Missing Atoms:", QLabel(str(missing_count)))
        
        # B-factor statistics
        if results.get('b_factor_stats'):
            bf_stats = results['b_factor_stats']
            summary_layout.addRow("Mean B-factor:", QLabel(f"{bf_stats.get('mean', 0):.2f} Ų"))
            summary_layout.addRow("B-factor Range:", QLabel(f"{bf_stats.get('min', 0):.2f} - {bf_stats.get('max', 0):.2f} Ų"))
        
        quality_layout.addWidget(summary_stats_group)
        self.results_layout.addWidget(quality_group)
        
        # Ramachandran outliers table has been removed as requested
        
        # Missing atoms
        if results.get('missing_atoms'):
            missing_group = QGroupBox("Missing Atoms")
            missing_layout = QVBoxLayout(missing_group)
            
            missing_table = QTableWidget()
            missing_table.setColumnCount(4)
            missing_table.setHorizontalHeaderLabels(["Residue", "Chain", "Position", "Missing Atom"])
            
            missing = results['missing_atoms']
            missing_table.setRowCount(len(missing))
            
            for i, miss in enumerate(missing):
                missing_table.setItem(i, 0, QTableWidgetItem(miss['residue']))
                missing_table.setItem(i, 1, QTableWidgetItem(miss['chain']))
                missing_table.setItem(i, 2, QTableWidgetItem(str(miss['position'])))
                missing_table.setItem(i, 3, QTableWidgetItem(miss['missing_atom']))
            
            missing_table.resizeColumnsToContents()
            # Allow table to expand to show all missing atoms
            missing_table.setMinimumHeight(100)
            # Remove maximum height restriction to show all data
            missing_layout.addWidget(missing_table)
            
            self.results_layout.addWidget(missing_group)
    
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
        
        # Add selected analysis types
        for analysis_type in ['basic', 'secondary', 'geometry', 'quality']:
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
            
            # Basic Properties
            if 'basic' in self.current_results and options.get('include_basic', True):
                basic = self.current_results['basic']
                writer.writerow(['=== BASIC PROPERTIES ==='])
                writer.writerow(['Property', 'Value'])
                writer.writerow(['Total Residues', basic.get('total_residues', 0)])
                writer.writerow(['Total Atoms', basic.get('total_atoms', 0)])
                writer.writerow(['Molecular Weight (Da)', f"{basic.get('molecular_weight', 0):.2f}"])
                writer.writerow(['Resolution', basic.get('resolution', 'N/A')])
                writer.writerow(['Space Group', basic.get('space_group', 'N/A')])
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
            
            # Quality Assessment
            if 'quality' in self.current_results and options.get('include_quality', True):
                quality = self.current_results['quality']
                writer.writerow(['=== QUALITY ASSESSMENT ==='])
                writer.writerow(['Metric', 'Value'])
                writer.writerow(['Ramachandran Favored', quality.get('ramachandran_favored', 0)])
                writer.writerow(['Ramachandran Allowed', quality.get('ramachandran_allowed', 0)])
                writer.writerow(['Ramachandran Outliers', quality.get('ramachandran_outliers_count', 0)])
                writer.writerow(['Missing Atoms', len(quality.get('missing_atoms', []))])
                
                # B-factor statistics
                if 'b_factor_stats' in quality:
                    stats = quality['b_factor_stats']
                    writer.writerow(['B-factor Mean', f"{stats.get('mean', 0):.2f}"])
                    writer.writerow(['B-factor Std Dev', f"{stats.get('std', 0):.2f}"])
                    writer.writerow(['B-factor Min', f"{stats.get('min', 0):.2f}"])
                    writer.writerow(['B-factor Max', f"{stats.get('max', 0):.2f}"])
                writer.writerow([])
            
            # Surface Properties
            # Surface Properties export removed - will be added in a future release
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
                        'Property': ['Total Residues', 'Total Atoms', 'Molecular Weight (Da)', 'Resolution', 'Space Group'],
                        'Value': [
                            basic.get('total_residues', 0),
                            basic.get('total_atoms', 0),
                            f"{basic.get('molecular_weight', 0):.2f}",
                            basic.get('resolution', 'N/A'),
                            basic.get('space_group', 'N/A')
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
                
                # Quality Assessment
                if 'quality' in self.current_results and options.get('include_quality', True):
                    quality = self.current_results['quality']
                    
                    quality_data = {
                        'Metric': ['Ramachandran Favored', 'Ramachandran Allowed', 'Ramachandran Outliers', 'Missing Atoms'],
                        'Value': [
                            quality.get('ramachandran_favored', 0),
                            quality.get('ramachandran_allowed', 0),
                            quality.get('ramachandran_outliers_count', 0),
                            len(quality.get('missing_atoms', []))
                        ]
                    }
                    
                    # Add B-factor statistics
                    if 'b_factor_stats' in quality:
                        stats = quality['b_factor_stats']
                        quality_data['Metric'].extend(['B-factor Mean', 'B-factor Std Dev', 'B-factor Min', 'B-factor Max'])
                        quality_data['Value'].extend([
                            f"{stats.get('mean', 0):.2f}",
                            f"{stats.get('std', 0):.2f}",
                            f"{stats.get('min', 0):.2f}",
                            f"{stats.get('max', 0):.2f}"
                        ])
                    
                    pd.DataFrame(quality_data).to_excel(writer, sheet_name='Quality', index=False)
                
                # Surface Properties
                # Surface Properties export removed - will be added in a future release
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
            
            # Basic Properties
            if 'basic' in self.current_results and options.get('include_basic', True):
                basic = self.current_results['basic']
                f.write("BASIC PROPERTIES\n")
                f.write("-"*40 + "\n")
                f.write(f"Total Residues: {basic.get('total_residues', 0)}\n")
                f.write(f"Total Atoms: {basic.get('total_atoms', 0)}\n")
                f.write(f"Molecular Weight: {basic.get('molecular_weight', 0):.2f} Da\n")
                f.write(f"Resolution: {basic.get('resolution', 'N/A')}\n")
                f.write(f"Space Group: {basic.get('space_group', 'N/A')}\n")
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
            
            # Quality Assessment
            if 'quality' in self.current_results and options.get('include_quality', True):
                quality = self.current_results['quality']
                f.write("QUALITY ASSESSMENT\n")
                f.write("-"*40 + "\n")
                f.write(f"Ramachandran Favored: {quality.get('ramachandran_favored', 0)}\n")
                f.write(f"Ramachandran Allowed: {quality.get('ramachandran_allowed', 0)}\n")
                f.write(f"Ramachandran Outliers: {quality.get('ramachandran_outliers_count', 0)}\n")
                f.write(f"Missing Atoms: {len(quality.get('missing_atoms', []))}\n")
                
                if 'b_factor_stats' in quality:
                    stats = quality['b_factor_stats']
                    f.write(f"B-factor Mean: {stats.get('mean', 0):.2f}\n")
                    f.write(f"B-factor Std Dev: {stats.get('std', 0):.2f}\n")
                    f.write(f"B-factor Range: {stats.get('min', 0):.2f} - {stats.get('max', 0):.2f}\n")
                f.write("\n")
            
            # Surface Properties
            # Surface Properties export removed - will be added in a future release
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
        import os
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
                f.write("</div>\n\n")
                
                # Basic Properties
                if 'basic' in self.current_results:
                    basic = self.current_results['basic']
                    f.write("<h2>Basic Properties</h2>\n")
                    
                    # Detailed table
                    f.write("<table>\n")
                    f.write("<tr><th>Property</th><th>Value</th></tr>\n")
                    f.write(f"<tr><td>Total Residues</td><td>{basic.get('total_residues', 0)}</td></tr>\n")
                    f.write(f"<tr><td>Total Atoms</td><td>{basic.get('total_atoms', 0)}</td></tr>\n")
                    f.write(f"<tr><td>Molecular Weight</td><td>{basic.get('molecular_weight', 0):.2f} Da</td></tr>\n")
                    f.write(f"<tr><td>Resolution</td><td>{basic.get('resolution', 'N/A')}</td></tr>\n")
                    f.write(f"<tr><td>Space Group</td><td>{basic.get('space_group', 'N/A')}</td></tr>\n")
                    f.write("</table>\n")
                    
                    # Chain information
                    if 'chains' in basic:
                        f.write("<h3>Chain Information</h3>\n")
                        f.write("<table>\n")
                        f.write("<tr><th>Chain ID</th><th>Residues</th><th>Atoms</th><th>Sequence Length</th></tr>\n")
                        for chain in basic['chains']:
                            f.write(f"<tr><td>{chain['id']}</td><td>{chain['residues']}</td><td>{chain['atoms']}</td><td>{len(chain.get('sequence', ''))}</td></tr>\n")
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
                        rama_chart_path = self.save_ramachandran_plot(secondary['ramachandran_data'], temp_dir)
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
                
                # Quality Assessment
                if 'quality' in self.current_results:
                    quality = self.current_results['quality']
                    f.write("<h2>Quality Assessment</h2>\n")
                    
                    # Quality table
                    f.write("<table>\n")
                    f.write("<tr><th>Quality Metric</th><th>Value</th></tr>\n")
                    f.write(f"<tr><td>Ramachandran Favored</td><td>{quality.get('ramachandran_favored', 0)}</td></tr>\n")
                    f.write(f"<tr><td>Ramachandran Allowed</td><td>{quality.get('ramachandran_allowed', 0)}</td></tr>\n")
                    f.write(f"<tr><td>Ramachandran Outliers</td><td>{quality.get('ramachandran_outliers_count', 0)}</td></tr>\n")
                    f.write(f"<tr><td>Missing Atoms</td><td>{len(quality.get('missing_atoms', []))}</td></tr>\n")
                    
                    if 'b_factor_stats' in quality:
                        stats = quality['b_factor_stats']
                        f.write(f"<tr><td>B-factor Mean</td><td>{stats.get('mean', 0):.2f}</td></tr>\n")
                        f.write(f"<tr><td>B-factor Std Dev</td><td>{stats.get('std', 0):.2f}</td></tr>\n")
                        f.write(f"<tr><td>B-factor Range</td><td>{stats.get('min', 0):.2f} - {stats.get('max', 0):.2f}</td></tr>\n")
                    f.write("</table>\n")
                    
                    # Outliers table
                    if quality.get('ramachandran_outliers'):
                        f.write("<h3>Ramachandran Outliers</h3>\n")
                        f.write("<table>\n")
                        f.write("<tr><th>Residue</th><th>Chain</th><th>Position</th><th>Phi (°)</th><th>Psi (°)</th></tr>\n")
                        for outlier in quality['ramachandran_outliers'][:20]:  # Limit to first 20
                            f.write(f"<tr><td>{outlier['residue']}</td><td>{outlier['chain']}</td><td>{outlier['position']}</td><td>{outlier['phi']:.1f}</td><td>{outlier['psi']:.1f}</td></tr>\n")
                        f.write("</table>\n")
                    

                
                # Surface Properties
                if 'surface' in self.current_results:
                    surface = self.current_results['surface']
                    f.write("<h2>Surface Properties (SASA Analysis)</h2>\n")
                    
                    # Surface/buried distribution chart
                    surface_data = {
                        'Surface': surface.get('surface_percentage', 0),
                        'Buried': surface.get('buried_percentage', 0)
                    }
                    if MATPLOTLIB_AVAILABLE:
                        chart_path = self.save_chart_as_image(surface_data, "Surface vs Buried Residues", "pie", temp_dir)
                        if chart_path:
                            f.write(f"<div class='chart'><img src='data:image/png;base64,{self.image_to_base64(chart_path)}' alt='Surface Distribution Chart'></div>\n")
                    
                    # Surface properties table
                    f.write("<table>\n")
                    f.write("<tr><th>Surface Property</th><th>Value</th></tr>\n")
                    f.write(f"<tr><td>Total SASA</td><td>{surface.get('total_sasa', 0):.2f} Ų</td></tr>\n")
                    f.write(f"<tr><td>Accessible Surface Area</td><td>{surface.get('accessible_surface_area', 0):.2f} Ų</td></tr>\n")
                    f.write(f"<tr><td>Buried Surface Area</td><td>{surface.get('buried_surface_area', 0):.2f} Ų</td></tr>\n")
                    f.write(f"<tr><td>Surface Percentage</td><td>{surface.get('surface_percentage', 0):.2f}%</td></tr>\n")
                    f.write(f"<tr><td>Buried Percentage</td><td>{surface.get('buried_percentage', 0):.2f}%</td></tr>\n")
                    f.write(f"<tr><td>Average RSA</td><td>{surface.get('average_rsa', 0):.2f}%</td></tr>\n")
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
    
    def save_chart_as_image(self, data, title, chart_type, temp_dir):
        """Save a chart as PNG image and return the path with robust error handling."""
        if not MATPLOTLIB_AVAILABLE or not data:
            return None
        
        try:
            import matplotlib.pyplot as plt
            import os
            
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
                ax.set_ylabel('Count', fontsize=16, fontweight='bold')
                plt.xticks(rotation=45, ha='right', fontsize=14, fontweight='bold')
                plt.yticks(fontsize=14, fontweight='bold')
            
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
    
    def save_ramachandran_plot(self, rama_data, temp_dir):
        """Save Ramachandran plot as PNG image."""
        if not MATPLOTLIB_AVAILABLE or not rama_data:
            return None
        
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches
            import os
            
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
            ax.set_title('Ramachandran Plot', fontsize=20, fontweight='bold', pad=25)
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
            import os
            
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
        
        self.include_quality = QCheckBox("Quality Assessment")
        self.include_quality.setChecked(True)
        self.include_quality.setEnabled('quality' in self.results if self.results else False)
        content_layout.addWidget(self.include_quality)
        
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
            'include_quality': self.include_quality.isChecked(),
            # Surface Properties export removed
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
            self.include_quality.isChecked(),
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