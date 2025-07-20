#!/usr/bin/env python3
"""
Structural Analysis Tools for PicoMol.

This module provides comprehensive protein structural analysis including:
- Secondary structure analysis
- Geometric property calculations
- Surface area and volume estimation
- Cavity and binding site detection
- Structural quality assessment
- Ramachandran plot analysis
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
    QGridLayout, QTextBrowser, QApplication
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
        self.surface_method = 'fast'  # Default surface analysis method
    
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
                method_name = "fast neighbor counting" if self.surface_method == 'fast' else "accurate SASA calculation"
                self.progress_update.emit(f"Calculating surface properties using {method_name}...")
                results['surface'] = self.analyze_surface_properties(method=self.surface_method)
            
            if 'quality' in self.analysis_types:
                self.progress_update.emit("Assessing structural quality...")
                results['quality'] = self.analyze_structural_quality()
            
            if 'cavities' in self.analysis_types:
                self.progress_update.emit("Detecting binding sites and cavities...")
                results['cavities'] = self.detect_cavities()
            
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
        """Analyze secondary structure using simple geometric rules."""
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
        
        for model in self.structure:
            for chain in model:
                residues = [res for res in chain if is_aa(res)]
                total_residues += len(residues)
                
                # Calculate phi/psi angles and classify secondary structure
                for i in range(1, len(residues) - 1):
                    try:
                        phi, psi = self.calculate_phi_psi(residues, i)
                        
                        # Store for Ramachandran plot
                        results['ramachandran_data'].append({
                            'residue': residues[i].get_resname(),
                            'phi': phi,
                            'psi': psi,
                            'chain': chain.id,
                            'position': residues[i].id[1]
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
    
    def analyze_surface_properties(self, method='fast'):
        """Analyze surface properties using optimized algorithms.
        
        Args:
            method (str): 'fast' for neighbor counting, 'accurate' for optimized SASA
        
        Returns:
            dict: Surface analysis results
        """
        import time
        start_time = time.time()
        
        if method == 'fast':
            results = self.analyze_surface_properties_fast()
        else:
            results = self.analyze_surface_properties_accurate()
        
        # Add timing information
        calculation_time = time.time() - start_time
        results['calculation_time'] = calculation_time
        results['method_used'] = method
        
        self.progress_update.emit(f"Surface analysis completed in {calculation_time:.2f} seconds using {method} method")
        
        return results
    
    def analyze_surface_properties_fast(self):
        """Fast surface analysis using CA neighbor counting.
        
        This method is ~100x faster than full SASA calculation and provides
        good surface/buried classification for most purposes.
        """
        results = {
            'accessible_surface_area': 0.0,
            'buried_surface_area': 0.0,
            'surface_residues': [],
            'buried_residues': [],
            'surface_percentage': 0.0,
            'buried_percentage': 0.0,
            'total_sasa': 0.0,
            'average_rsa': 0.0,
            'neighbor_threshold': 16  # Threshold used for classification
        }
        
        all_residues = []
        total_estimated_sasa = 0.0
        
        # Get all CA atoms for neighbor search
        ca_atoms = []
        residue_to_ca = {}
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue) and 'CA' in residue:
                        ca_atom = residue['CA']
                        ca_atoms.append(ca_atom)
                        residue_to_ca[residue] = ca_atom
        
        # Use NeighborSearch for efficient neighbor finding
        if ca_atoms:
            try:
                neighbor_search = NeighborSearch(ca_atoms)
                
                for model in self.structure:
                    for chain in model:
                        for residue in chain:
                            if is_aa(residue) and residue in residue_to_ca:
                                ca_atom = residue_to_ca[residue]
                                
                                # Count neighbors within 8Å
                                neighbors = neighbor_search.search(ca_atom.get_coord(), 8.0)
                                neighbor_count = len([n for n in neighbors if n != ca_atom])
                                
                                # Estimate SASA based on neighbor count
                                res_name = residue.get_resname()
                                max_sasa = MAX_SASA_VALUES.get(res_name, 150.0)
                                
                                # Empirical relationship: fewer neighbors = more surface exposed
                                # This is a simplified model but works well in practice
                                if neighbor_count < 16:  # Surface threshold
                                    estimated_rsa = max(20.0, 100.0 - (neighbor_count * 5.0))
                                else:  # Buried
                                    estimated_rsa = max(0.0, 20.0 - ((neighbor_count - 16) * 2.0))
                                
                                estimated_sasa = (estimated_rsa / 100.0) * max_sasa
                                
                                residue_info = {
                                    'residue': res_name,
                                    'chain': chain.id,
                                    'position': residue.id[1],
                                    'sasa': estimated_sasa,
                                    'rsa': estimated_rsa,
                                    'max_sasa': max_sasa,
                                    'neighbors': neighbor_count
                                }
                                
                                all_residues.append(residue_info)
                                total_estimated_sasa += estimated_sasa
                                
                                # Classification based on neighbor count
                                if neighbor_count < 16:  # Surface
                                    results['surface_residues'].append(residue_info)
                                else:  # Buried
                                    results['buried_residues'].append(residue_info)
                                    
            except Exception as e:
                self.progress_update.emit(f"NeighborSearch failed, using fallback method: {e}")
                # Fallback to simple distance calculation
                return self._analyze_surface_fallback()
        
        # Calculate summary statistics
        total_residues = len(all_residues)
        if total_residues > 0:
            surface_count = len(results['surface_residues'])
            buried_count = len(results['buried_residues'])
            
            results['surface_percentage'] = (surface_count / total_residues) * 100
            results['buried_percentage'] = (buried_count / total_residues) * 100
            results['total_sasa'] = total_estimated_sasa
            results['accessible_surface_area'] = sum(r['sasa'] for r in results['surface_residues'])
            results['buried_surface_area'] = sum(r['sasa'] for r in results['buried_residues'])
            results['average_rsa'] = sum(r['rsa'] for r in all_residues) / total_residues
        
        return results
    
    def analyze_surface_properties_accurate(self):
        """Accurate surface analysis using optimized SASA calculation.
        
        Uses an optimized Shrake-Rupley algorithm with NeighborSearch for
        better performance while maintaining accuracy.
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
        
        # Get all atoms for neighbor search
        all_atoms = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        for atom in residue:
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
                            # Calculate optimized SASA for this residue
                            residue_sasa = self.calculate_residue_sasa_optimized(residue, neighbor_search)
                            
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
                            
                            # Classification based on RSA threshold
                            if rsa < 20.0:
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
        """Fallback surface analysis using simple distance calculations."""
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
        
        all_residues = []
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue) and 'CA' in residue:
                        ca_atom = residue['CA']
                        neighbor_count = self.count_neighbors(ca_atom, model, radius=8.0, ca_only=True)
                        
                        res_name = residue.get_resname()
                        max_sasa = MAX_SASA_VALUES.get(res_name, 150.0)
                        
                        # Simple classification
                        if neighbor_count < 16:
                            estimated_rsa = 50.0  # Surface
                            estimated_sasa = 0.5 * max_sasa
                            results['surface_residues'].append({
                                'residue': res_name,
                                'chain': chain.id,
                                'position': residue.id[1],
                                'sasa': estimated_sasa,
                                'rsa': estimated_rsa,
                                'max_sasa': max_sasa,
                                'neighbors': neighbor_count
                            })
                        else:
                            estimated_rsa = 10.0  # Buried
                            estimated_sasa = 0.1 * max_sasa
                            results['buried_residues'].append({
                                'residue': res_name,
                                'chain': chain.id,
                                'position': residue.id[1],
                                'sasa': estimated_sasa,
                                'rsa': estimated_rsa,
                                'max_sasa': max_sasa,
                                'neighbors': neighbor_count
                            })
                        
                        all_residues.append({
                            'residue': res_name,
                            'sasa': estimated_sasa,
                            'rsa': estimated_rsa
                        })
        
        # Calculate summary statistics
        total_residues = len(all_residues)
        if total_residues > 0:
            surface_count = len(results['surface_residues'])
            buried_count = len(results['buried_residues'])
            
            results['surface_percentage'] = (surface_count / total_residues) * 100
            results['buried_percentage'] = (buried_count / total_residues) * 100
            results['total_sasa'] = sum(r['sasa'] for r in all_residues)
            results['accessible_surface_area'] = sum(r['sasa'] for r in results['surface_residues'])
            results['buried_surface_area'] = sum(r['sasa'] for r in results['buried_residues'])
            results['average_rsa'] = sum(r['rsa'] for r in all_residues) / total_residues
        
        return results
    
    def calculate_residue_sasa_optimized(self, residue, neighbor_search, probe_radius=1.4, n_points=30):
        """Calculate SASA for a residue using optimized Shrake-Rupley algorithm.
        
        Args:
            residue: The residue to calculate SASA for
            neighbor_search: Pre-built NeighborSearch object for efficient neighbor finding
            probe_radius: Probe radius in Angstroms (1.4Å for water)
            n_points: Number of points to generate on each atom sphere (reduced for speed)
        
        Returns:
            float: SASA value in Ų
        """
        total_sasa = 0.0
        
        # Calculate SASA for each atom in the residue
        for atom in residue:
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
    
    def calculate_residue_sasa(self, residue, model, probe_radius=1.4, n_points=100):
        """Calculate SASA for a residue using Shrake-Rupley algorithm.
        
        DEPRECATED: Use calculate_residue_sasa_optimized for better performance.
        
        Args:
            residue: The residue to calculate SASA for
            model: The PDB model containing all atoms
            probe_radius: Probe radius in Angstroms (1.4Å for water)
            n_points: Number of points to generate on each atom sphere
        
        Returns:
            float: SASA value in Ų
        """
        total_sasa = 0.0
        
        # Get all atoms in the protein for neighbor checking
        all_atoms = []
        for chain in model:
            for res in chain:
                for atom in res:
                    all_atoms.append(atom)
        
        # Calculate SASA for each atom in the residue
        for atom in residue:
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
        """Analyze structural quality metrics."""
        results = {
            'ramachandran_outliers': [],
            'ramachandran_favored': 0,
            'ramachandran_allowed': 0,
            'ramachandran_outliers_count': 0,
            'clashes': [],
            'missing_atoms': [],
            'b_factor_stats': {}
        }
        
        # Analyze Ramachandran plot data
        ramachandran_data = []
        
        for model in self.structure:
            for chain in model:
                residues = [res for res in chain if is_aa(res)]
                
                for i in range(1, len(residues) - 1):
                    try:
                        phi, psi = self.calculate_phi_psi(residues, i)
                        
                        # Classify Ramachandran regions
                        region = self.classify_ramachandran_region(phi, psi)
                        
                        if region == 'favored':
                            results['ramachandran_favored'] += 1
                        elif region == 'allowed':
                            results['ramachandran_allowed'] += 1
                        else:
                            results['ramachandran_outliers_count'] += 1
                            results['ramachandran_outliers'].append({
                                'residue': residues[i].get_resname(),
                                'chain': chain.id,
                                'position': residues[i].id[1],
                                'phi': phi,
                                'psi': psi
                            })
                        
                        ramachandran_data.append({
                            'phi': phi,
                            'psi': psi,
                            'region': region
                        })
                        
                    except Exception:
                        continue
        
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
    
    def classify_ramachandran_region(self, phi, psi):
        """Classify phi/psi angles into Ramachandran regions."""
        # Simplified Ramachandran classification
        # Favored regions (core)
        if (-180 <= phi <= -30 and -70 <= psi <= 50):  # Alpha helix
            return 'favored'
        elif (-180 <= phi <= -30 and 90 <= psi <= 180):  # Beta sheet
            return 'favored'
        elif (-90 <= phi <= 30 and -180 <= psi <= -90):  # Beta sheet
            return 'favored'
        
        # Allowed regions (expanded)
        elif (-180 <= phi <= -30 and -180 <= psi <= 180):  # Extended allowed
            return 'allowed'
        elif (30 <= phi <= 180 and -180 <= psi <= 180):   # Left-handed helix region
            return 'allowed'
        
        # Outliers
        else:
            return 'outlier'
    
    def detect_cavities(self):
        """Detect potential binding sites and cavities using improved geometric methods.
        
        Uses a combination of grid-based and probe-based approaches inspired by
        conventional cavity detection algorithms like CASTp and fpocket.
        """
        results = {
            'potential_cavities': [],
            'surface_pockets': [],
            'cavity_volume': 0.0,
            'largest_cavity': None,
            'druggable_cavities': []
        }
        
        # Get all protein atoms
        all_atoms = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):  # Only protein atoms
                        for atom in residue:
                            all_atoms.append(atom)
        
        if not all_atoms:
            return results
        
        # Use improved cavity detection with multiple probe sizes
        cavity_points = self.detect_cavity_points_improved(all_atoms)
        
        if cavity_points:
            # Cluster cavity points into distinct cavities
            cavities = self.cluster_cavity_points_improved(cavity_points, all_atoms)
            
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
        
        sample_count = 0
        for x in x_range:
            for y in y_range:
                for z in z_range:
                    point = np.array([x, y, z])
                    
                    # Check if point is in a cavity using improved criteria
                    if self.is_cavity_point_improved(point, all_atoms, probe_radius):
                        cavity_points.append(point)
                    
                    sample_count += 1
                    if sample_count >= max_grid_points:
                        break
                if sample_count >= max_grid_points:
                    break
            if sample_count >= max_grid_points:
                break
        
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
        
        points = np.array(points)
        clusters = []
        used = set()
        
        for i, point in enumerate(points):
            if i in used:
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
        
        for i, cluster in enumerate(clusters):
            cluster_points = cluster['points']
            
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


class StructuralAnalysisTab(QWidget):
    """Tab for structural analysis tools."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_app = parent
        self.analysis_worker = None
        self.current_structure_path = None
        self.pdb_list = PDBList() if BIOPYTHON_AVAILABLE else None
        self.pdb_parser = PDBParser(QUIET=True) if BIOPYTHON_AVAILABLE else None
        
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
        
        self.surface_checkbox = QCheckBox("Surface Properties")
        self.surface_checkbox.setChecked(False)
        self.surface_checkbox.setToolTip("Analyze surface properties")
        options_layout.addWidget(self.surface_checkbox)
        
        # Surface analysis method selection
        self.surface_method_combo = QComboBox()
        self.surface_method_combo.addItems(["Fast (Neighbor Counting)", "Accurate (Optimized SASA)"])
        self.surface_method_combo.setCurrentIndex(0)  # Default to fast
        self.surface_method_combo.setToolTip("Choose surface analysis method: Fast (~1s) vs Accurate (~10s)")
        self.surface_method_combo.setEnabled(False)  # Initially disabled
        options_layout.addWidget(self.surface_method_combo)
        
        # Connect surface checkbox to enable/disable method selection
        self.surface_checkbox.stateChanged.connect(self.on_surface_checkbox_changed)
        
        self.cavities_checkbox = QCheckBox("Binding Sites")
        self.cavities_checkbox.setChecked(False)
        self.cavities_checkbox.setToolTip("Detect potential binding sites (slower)")
        options_layout.addWidget(self.cavities_checkbox)
        
        options_layout.addStretch()
        control_layout.addLayout(options_layout)
        
        # Analysis button
        analyze_btn = QPushButton("Analyze Structure")
        analyze_btn.setToolTip("Perform structural analysis with selected options")
        analyze_btn.clicked.connect(self.analyze_structure)
        control_layout.addWidget(analyze_btn)
        
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
        
        # Set minimum size for results area to accommodate charts
        self.results_area.setMinimumHeight(600)
        
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
            "<li><b>Surface Properties:</b> Fast neighbor counting or accurate SASA calculation</li>"
            "<li><b>Binding Sites:</b> Cavity detection and potential binding sites</li>"
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
                    from .enhanced_pdb_puller_async import OptimizedPDBPuller
                    temp_puller = OptimizedPDBPuller(self.pulled_structures_dir)
                    comprehensive_data = temp_puller.fetch_pdb_data_async(
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
        if self.surface_checkbox.isChecked():
            analysis_types.append('surface')
        if self.cavities_checkbox.isChecked():
            analysis_types.append('cavities')
        
        if not analysis_types:
            QMessageBox.warning(self, "No Analysis Selected", "Please select at least one analysis type.")
            return
        
        # Show progress
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        
        # Get surface analysis method if surface analysis is selected
        surface_method = 'fast'
        if self.surface_checkbox.isChecked():
            if self.surface_method_combo.currentIndex() == 1:
                surface_method = 'accurate'
        
        # Start analysis in worker thread
        self.analysis_worker = StructuralAnalysisWorker(self.current_structure_path, analysis_types)
        self.analysis_worker.surface_method = surface_method  # Pass method to worker
        self.analysis_worker.analysis_complete.connect(self.display_results)
        self.analysis_worker.error_occurred.connect(self.handle_analysis_error)
        self.analysis_worker.progress_update.connect(self.update_progress)
        self.analysis_worker.start()
    
    def on_surface_checkbox_changed(self, state):
        """Enable/disable surface method selection based on checkbox state."""
        self.surface_method_combo.setEnabled(state == Qt.Checked)
    
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
        """Create a pie chart widget."""
        if not MATPLOTLIB_AVAILABLE:
            return None
        
        # Use larger figure for amino acid composition
        if "Amino Acid" in title:
            fig = Figure(figsize=(12, 10), dpi=100)
        else:
            fig = Figure(figsize=(8, 6), dpi=100)
        
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        
        labels = list(data.keys())
        values = list(data.values())
        
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
            canvas.setMinimumSize(900, 750)
        else:
            canvas.setMinimumSize(600, 450)
        
        return canvas
    
    def create_bar_chart(self, data, title, xlabel, ylabel):
        """Create a bar chart widget."""
        if not MATPLOTLIB_AVAILABLE:
            return None
        
        fig = Figure(figsize=(10, 6), dpi=100)
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
        canvas.setMinimumSize(750, 450)
        return canvas
    
    def create_ramachandran_plot(self, rama_data):
        """Create a proper conventional Ramachandran plot."""
        if not MATPLOTLIB_AVAILABLE or not rama_data:
            return None
        
        fig = Figure(figsize=(10, 10), dpi=100)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        
        # Extract phi and psi angles
        phi_angles = [d['phi'] for d in rama_data]
        psi_angles = [d['psi'] for d in rama_data]
        
        # Define conventional Ramachandran regions based on crystallographic data
        # Most favored regions (core)
        favored_regions = [
            # Right-handed alpha helix (αR)
            patches.Polygon([(-180, -70), (-30, -70), (-30, 50), (-180, 50)], 
                          closed=True, facecolor='darkgreen', alpha=0.3, edgecolor='darkgreen', linewidth=2),
            # Beta sheet region
            patches.Polygon([(-180, 90), (-30, 90), (-30, 180), (-180, 180)], 
                          closed=True, facecolor='darkgreen', alpha=0.3, edgecolor='darkgreen', linewidth=2),
            # Extended beta region
            patches.Polygon([(-180, -180), (-30, -180), (-30, -90), (-180, -90)], 
                          closed=True, facecolor='darkgreen', alpha=0.3, edgecolor='darkgreen', linewidth=2),
            # Left-handed alpha helix (αL) - rare but allowed
            patches.Polygon([(30, 30), (90, 30), (90, 90), (30, 90)], 
                          closed=True, facecolor='darkgreen', alpha=0.3, edgecolor='darkgreen', linewidth=2)
        ]
        
        # Additional allowed regions (lighter green)
        allowed_regions = [
            # Extended allowed alpha region
            patches.Polygon([(-180, -90), (-30, -90), (-30, -70), (-180, -70)], 
                          closed=True, facecolor='lightgreen', alpha=0.3, edgecolor='green', linewidth=1),
            # Extended allowed beta region
            patches.Polygon([(-30, 90), (30, 90), (30, 180), (-30, 180)], 
                          closed=True, facecolor='lightgreen', alpha=0.3, edgecolor='green', linewidth=1),
            # Bridge region
            patches.Polygon([(-90, -30), (-30, -30), (-30, 30), (-90, 30)], 
                          closed=True, facecolor='lightgreen', alpha=0.3, edgecolor='green', linewidth=1),
            # Left-handed extended region
            patches.Polygon([(30, -180), (180, -180), (180, -90), (30, -90)], 
                          closed=True, facecolor='lightgreen', alpha=0.3, edgecolor='green', linewidth=1)
        ]
        
        # Add regions to plot
        for region in allowed_regions:
            ax.add_patch(region)
        for region in favored_regions:
            ax.add_patch(region)
        
        # Create scatter plot with different colors for different regions
        # Classify points by region
        favored_phi, favored_psi = [], []
        allowed_phi, allowed_psi = [], []
        outlier_phi, outlier_psi = [], []
        
        for phi, psi in zip(phi_angles, psi_angles):
            if self.is_in_favored_region(phi, psi):
                favored_phi.append(phi)
                favored_psi.append(psi)
            elif self.is_in_allowed_region(phi, psi):
                allowed_phi.append(phi)
                allowed_psi.append(psi)
            else:
                outlier_phi.append(phi)
                outlier_psi.append(psi)
        
        # Plot points with different colors
        if favored_phi:
            ax.scatter(favored_phi, favored_psi, c='darkblue', s=20, alpha=0.7, label=f'Favored ({len(favored_phi)})')
        if allowed_phi:
            ax.scatter(allowed_phi, allowed_psi, c='orange', s=20, alpha=0.7, label=f'Allowed ({len(allowed_phi)})')
        if outlier_phi:
            ax.scatter(outlier_phi, outlier_psi, c='red', s=20, alpha=0.8, label=f'Outliers ({len(outlier_phi)})')
        
        # Set up the plot
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xlabel('φ (phi) degrees', fontsize=12)
        ax.set_ylabel('ψ (psi) degrees', fontsize=12)
        ax.set_title('Ramachandran Plot', fontsize=14, fontweight='bold')
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        ax.set_xticks(range(-180, 181, 60))
        ax.set_yticks(range(-180, 181, 60))
        
        # Add region labels
        ax.text(-105, -20, 'αR', fontsize=16, fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        ax.text(-105, 135, 'β', fontsize=16, fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        ax.text(-105, -135, 'β', fontsize=16, fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        ax.text(60, 60, 'αL', fontsize=16, fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        # Add legend
        legend = ax.legend(loc='upper right', fontsize=10)
        legend.get_frame().set_alpha(0.9)
        
        # Add statistics text
        total_points = len(phi_angles)
        favored_pct = (len(favored_phi) / total_points * 100) if total_points > 0 else 0
        allowed_pct = (len(allowed_phi) / total_points * 100) if total_points > 0 else 0
        outlier_pct = (len(outlier_phi) / total_points * 100) if total_points > 0 else 0
        
        stats_text = f"Favored: {favored_pct:.1f}%\nAllowed: {allowed_pct:.1f}%\nOutliers: {outlier_pct:.1f}%"
        ax.text(-170, -150, stats_text, fontsize=10, 
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.9))
        
        # Improve tick labels
        ax.tick_params(axis='both', which='major', labelsize=10)
        
        fig.tight_layout()
        canvas.setMinimumSize(750, 750)
        return canvas
    
    def is_in_favored_region(self, phi, psi):
        """Check if phi/psi angles are in favored regions."""
        # Right-handed alpha helix
        if -180 <= phi <= -30 and -70 <= psi <= 50:
            return True
        # Beta sheet regions
        if -180 <= phi <= -30 and 90 <= psi <= 180:
            return True
        if -180 <= phi <= -30 and -180 <= psi <= -90:
            return True
        # Left-handed alpha helix (rare)
        if 30 <= phi <= 90 and 30 <= psi <= 90:
            return True
        return False
    
    def is_in_allowed_region(self, phi, psi):
        """Check if phi/psi angles are in allowed regions."""
        # Extended allowed regions around favored ones
        if -180 <= phi <= -30 and -90 <= psi <= -70:
            return True
        if -30 <= phi <= 30 and 90 <= psi <= 180:
            return True
        if -90 <= phi <= -30 and -30 <= psi <= 30:
            return True
        if 30 <= phi <= 180 and -180 <= psi <= -90:
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
        canvas.setMinimumSize(750, 450)
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
        
        if 'surface' in results:
            self.display_surface_results(results['surface'])
        
        if 'cavities' in results:
            self.display_cavity_results(results['cavities'])
        
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
            rama_plot = self.create_ramachandran_plot(results['ramachandran_data'])
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
        total_rama = results.get('ramachandran_favored', 0) + results.get('ramachandran_allowed', 0) + results.get('ramachandran_outliers_count', 0)
        if total_rama > 0:
            rama_summary_group = QGroupBox("Ramachandran Quality")
            rama_summary_layout = QHBoxLayout(rama_summary_group)
            
            # Statistics
            stats_widget = QWidget()
            stats_layout = QFormLayout(stats_widget)
            
            favored_pct = (results.get('ramachandran_favored', 0) / total_rama) * 100
            allowed_pct = (results.get('ramachandran_allowed', 0) / total_rama) * 100
            outlier_pct = (results.get('ramachandran_outliers_count', 0) / total_rama) * 100
            
            stats_layout.addRow("Favored:", QLabel(f"{favored_pct:.1f}% ({results.get('ramachandran_favored', 0)} residues)"))
            stats_layout.addRow("Allowed:", QLabel(f"{allowed_pct:.1f}% ({results.get('ramachandran_allowed', 0)} residues)"))
            stats_layout.addRow("Outliers:", QLabel(f"{outlier_pct:.1f}% ({results.get('ramachandran_outliers_count', 0)} residues)"))
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
        
        # Ramachandran outliers
        if results.get('ramachandran_outliers'):
            outlier_group = QGroupBox("Ramachandran Outliers")
            outlier_layout = QVBoxLayout(outlier_group)
            
            outlier_table = QTableWidget()
            outlier_table.setColumnCount(5)
            outlier_table.setHorizontalHeaderLabels(["Residue", "Chain", "Position", "Phi (°)", "Psi (°)"])
            
            outliers = results['ramachandran_outliers']
            outlier_table.setRowCount(len(outliers))
            
            for i, outlier in enumerate(outliers):
                outlier_table.setItem(i, 0, QTableWidgetItem(outlier['residue']))
                outlier_table.setItem(i, 1, QTableWidgetItem(outlier['chain']))
                outlier_table.setItem(i, 2, QTableWidgetItem(str(outlier['position'])))
                outlier_table.setItem(i, 3, QTableWidgetItem(f"{outlier['phi']:.1f}"))
                outlier_table.setItem(i, 4, QTableWidgetItem(f"{outlier['psi']:.1f}"))
            
            outlier_table.resizeColumnsToContents()
            # Allow table to expand to show all outliers
            outlier_table.setMinimumHeight(100)
            # Remove maximum height restriction to show all data
            outlier_layout.addWidget(outlier_table)
            
            self.results_layout.addWidget(outlier_group)
        
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
        cavity_group = QGroupBox("Binding Sites and Cavities Analysis")
        cavity_layout = QVBoxLayout(cavity_group)
        
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