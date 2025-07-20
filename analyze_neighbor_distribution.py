#!/usr/bin/env python3
"""
Analyze the distribution of CA neighbors to find optimal threshold.
"""

import os
import sys
sys.path.append('src')

try:
    from Bio.PDB import PDBParser
    import numpy as np
    
    # Test with hemoglobin (larger protein)
    structure_path = "data/pulled_structures/structures/4HHB.pdb"
    
    if os.path.exists(structure_path):
        print(f"Analyzing neighbor distribution for {os.path.basename(structure_path)}")
        
        # Load structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("test", structure_path)
        model = list(structure)[0]
        
        # Collect all CA neighbor counts
        neighbor_counts = []
        residue_info = []
        
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atom = residue['CA']
                    ca_coord = ca_atom.get_coord()
                    
                    # Count CA neighbors within 8Å
                    count = 0
                    for other_chain in model:
                        for other_residue in other_chain:
                            if 'CA' in other_residue:
                                other_ca = other_residue['CA']
                                if other_ca != ca_atom:
                                    distance = np.linalg.norm(ca_coord - other_ca.get_coord())
                                    if distance <= 8.0:
                                        count += 1
                    
                    neighbor_counts.append(count)
                    residue_info.append({
                        'residue': residue.get_resname(),
                        'chain': chain.id,
                        'position': residue.id[1],
                        'neighbors': count
                    })
        
        # Analyze distribution
        neighbor_counts = np.array(neighbor_counts)
        
        print(f"\\nNeighbor count statistics:")
        print(f"Total residues: {len(neighbor_counts)}")
        print(f"Min neighbors: {np.min(neighbor_counts)}")
        print(f"Max neighbors: {np.max(neighbor_counts)}")
        print(f"Mean neighbors: {np.mean(neighbor_counts):.1f}")
        print(f"Median neighbors: {np.median(neighbor_counts):.1f}")
        print(f"Std deviation: {np.std(neighbor_counts):.1f}")
        
        # Show distribution
        print(f"\\nNeighbor count distribution:")
        for i in range(0, int(np.max(neighbor_counts)) + 1, 2):
            count_in_range = np.sum((neighbor_counts >= i) & (neighbor_counts < i + 2))
            percentage = count_in_range / len(neighbor_counts) * 100
            print(f"  {i:2d}-{i+1:2d} neighbors: {count_in_range:3d} residues ({percentage:5.1f}%)")
        
        # Test different thresholds
        print(f"\\nTesting different thresholds:")
        for threshold in [12, 14, 16, 18, 20]:
            surface_count = np.sum(neighbor_counts < threshold)
            buried_count = np.sum(neighbor_counts >= threshold)
            surface_pct = surface_count / len(neighbor_counts) * 100
            buried_pct = buried_count / len(neighbor_counts) * 100
            print(f"  Threshold {threshold}: {surface_pct:5.1f}% surface, {buried_pct:5.1f}% buried")
        
        # Show some examples of high and low neighbor counts
        print(f"\\nExamples of residues with few neighbors (likely surface):")
        low_neighbor_residues = [r for r in residue_info if r['neighbors'] <= 8]
        for r in low_neighbor_residues[:5]:
            print(f"  {r['residue']} {r['chain']}:{r['position']} ({r['neighbors']} neighbors)")
        
        print(f"\\nExamples of residues with many neighbors (likely buried):")
        high_neighbor_residues = [r for r in residue_info if r['neighbors'] >= 18]
        for r in high_neighbor_residues[:5]:
            print(f"  {r['residue']} {r['chain']}:{r['position']} ({r['neighbors']} neighbors)")
        
    else:
        print(f"Test structure not found: {structure_path}")
        
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()