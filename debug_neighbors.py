#!/usr/bin/env python3
"""
Debug the neighbor counting to understand why counts are so high.
"""

import os
import sys
sys.path.append('src')

try:
    from Bio.PDB import PDBParser
    import numpy as np
    
    # Test with a small structure
    structure_path = "data/pulled_structures/structures/1CRN.pdb"
    
    if os.path.exists(structure_path):
        print(f"Debugging neighbor counting with {os.path.basename(structure_path)}")
        
        # Load structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("test", structure_path)
        
        # Get first model
        model = list(structure)[0]
        
        # Get first few residues to analyze
        residues = []
        for chain in model:
            for residue in chain:
                if residue.has_id('CA'):
                    residues.append(residue)
                if len(residues) >= 3:  # Just test first 3 residues
                    break
            if len(residues) >= 3:
                break
        
        print(f"\\nAnalyzing first {len(residues)} residues:")
        
        for i, residue in enumerate(residues):
            ca_atom = residue['CA']
            ca_coord = ca_atom.get_coord()
            
            print(f"\\nResidue {i+1}: {residue.get_resname()} {residue.get_parent().id}:{residue.id[1]}")
            print(f"CA position: {ca_coord}")
            
            # Count neighbors at different radii
            for radius in [6.0, 8.0, 10.0, 12.0]:
                count = 0
                atom_coord = ca_atom.get_coord()
                
                for chain in model:
                    for res in chain:
                        for other_atom in res:
                            if other_atom != ca_atom:
                                distance = np.linalg.norm(atom_coord - other_atom.get_coord())
                                if distance <= radius:
                                    count += 1
                
                print(f"  Neighbors within {radius}Å: {count}")
        
        # Also check total atoms in structure
        total_atoms = 0
        for chain in model:
            for residue in chain:
                for atom in residue:
                    total_atoms += 1
        
        print(f"\\nTotal atoms in structure: {total_atoms}")
        print(f"Total residues: {len([r for c in model for r in c])}")
        
    else:
        print(f"Test structure not found: {structure_path}")
        
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()