#!/usr/bin/env python3
"""
Quick test of the surface analysis fix.
"""

import os
import sys
sys.path.append('src')

try:
    from Bio.PDB import PDBParser
    from core.structural_analysis import StructuralAnalysisWorker
    
    # Test with a small structure
    structure_path = "data/pulled_structures/structures/1CRN.pdb"
    
    if os.path.exists(structure_path):
        print(f"Testing surface analysis fix with {os.path.basename(structure_path)}")
        
        # Load structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("test", structure_path)
        
        # Create worker instance
        worker = StructuralAnalysisWorker(structure_path, ['surface'])
        worker.structure = structure
        
        # Run just the surface analysis
        print("Running surface analysis...")
        results = worker.analyze_surface_properties()
        
        if results:
            total_residues = len(results['surface_residues']) + len(results['buried_residues'])
            surface_count = len(results['surface_residues'])
            buried_count = len(results['buried_residues'])
            
            surface_pct = (surface_count / total_residues * 100) if total_residues > 0 else 0
            buried_pct = (buried_count / total_residues * 100) if total_residues > 0 else 0
            
            print(f"\nResults for {os.path.basename(structure_path)}:")
            print(f"Total residues: {total_residues}")
            print(f"Surface residues: {surface_count} ({surface_pct:.1f}%)")
            print(f"Buried residues: {buried_count} ({buried_pct:.1f}%)")
            
            # Check if fix worked
            if buried_pct < 90:
                print("✓ SUCCESS: Fix is working - reasonable surface/buried distribution")
            else:
                print("✗ ISSUE: Still showing mostly buried residues")
                
            # Show some example residues
            print(f"\nExample surface residues (first 5):")
            for i, res in enumerate(results['surface_residues'][:5]):
                print(f"  {res['residue']} {res['chain']}:{res['position']} (RSA: {res['rsa']:.1f}%)")
                
            print(f"\nExample buried residues (first 5):")
            for i, res in enumerate(results['buried_residues'][:5]):
                print(f"  {res['residue']} {res['chain']}:{res['position']} (RSA: {res['rsa']:.1f}%)")
        else:
            print("No results returned")
    else:
        print(f"Test structure not found: {structure_path}")
        
except ImportError as e:
    print(f"Missing dependency: {e}")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()