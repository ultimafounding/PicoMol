#!/usr/bin/env python3
"""
Fixed Enhanced PDB Puller for PicoMol.

This module provides comprehensive PDB data fetching using working RCSB API endpoints:
- Structure files (PDB, mmCIF)
- Comprehensive metadata via RCSB PDB GraphQL API
- Experimental details and basic information
- Sequence information (FASTA)
- Improved error handling and fallback mechanisms
"""

import os
import json
import requests
import time
from urllib.parse import urljoin
from typing import Dict, List, Optional, Any
import warnings

try:
    from Bio.PDB import PDBList, PDBParser
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


class EnhancedPDBPuller:
    """Enhanced PDB data puller with fixed API endpoints."""
    
    def __init__(self, data_dir: str):
        """
        Initialize the enhanced PDB puller.
        
        Args:
            data_dir: Directory to store downloaded PDB data
        """
        self.data_dir = data_dir
        self.pdb_list = PDBList() if BIOPYTHON_AVAILABLE else None
        self.pdb_parser = PDBParser(QUIET=True) if BIOPYTHON_AVAILABLE else None
        
        # Updated RCSB PDB API endpoints
        self.graphql_api = "https://data.rcsb.org/graphql"
        self.files_base = "https://files.rcsb.org/download"
        
        # Create subdirectories for organized storage
        self.structures_dir = os.path.join(data_dir, "structures")
        self.metadata_dir = os.path.join(data_dir, "metadata")
        self.sequences_dir = os.path.join(data_dir, "sequences")
        self.validation_dir = os.path.join(data_dir, "validation")
        
        for directory in [self.structures_dir, self.metadata_dir, 
                         self.sequences_dir, self.validation_dir]:
            os.makedirs(directory, exist_ok=True)
    
    def fetch_comprehensive_pdb_data(self, pdb_id: str, 
                                   include_validation: bool = True,
                                   include_sequences: bool = True,
                                   include_mmcif: bool = False) -> Dict[str, Any]:
        """
        Fetch comprehensive PDB data for a given PDB ID using working endpoints.
        
        Args:
            pdb_id: PDB identifier (e.g., '1CRN')
            include_validation: Whether to fetch validation reports (limited)
            include_sequences: Whether to fetch FASTA sequences
            include_mmcif: Whether to fetch mmCIF format files
            
        Returns:
            Dictionary containing all fetched data and file paths
        """
        pdb_id = pdb_id.upper().strip()
        

        
        result = {
            'pdb_id': pdb_id,
            'files': {},
            'metadata': {},
            'sequences': {},
            'validation': {},
            'errors': []
        }
        
        try:
            # 1. Fetch structure files
            structure_files = self._fetch_structure_files(pdb_id, include_mmcif)
            result['files'].update(structure_files)
            
            # 2. Fetch comprehensive metadata using GraphQL
            metadata = self._fetch_metadata_graphql(pdb_id)
            result['metadata'] = metadata
            
            # 3. Fetch sequence information
            if include_sequences:
                sequences = self._fetch_sequences(pdb_id)
                result['sequences'] = sequences
            
            # 4. Fetch validation data (limited)
            if include_validation:
                validation = self._fetch_validation_data(pdb_id)
                result['validation'] = validation
            
            # 5. Save comprehensive data summary
            self._save_data_summary(pdb_id, result)
            

            
        except Exception as e:
            error_msg = f"Error fetching comprehensive data for {pdb_id}: {str(e)}"
            print(error_msg)
            result['errors'].append(error_msg)
        
        return result
    
    def _fetch_structure_files(self, pdb_id: str, include_mmcif: bool = False) -> Dict[str, str]:
        """Fetch structure files (PDB and optionally mmCIF)."""
        files = {}
        
        # Fetch PDB file
        try:
            pdb_path = os.path.join(self.structures_dir, f"{pdb_id}.pdb")
            
            if not os.path.exists(pdb_path):
                if BIOPYTHON_AVAILABLE:
                    # Use Biopython to download PDB file
                    retrieved_path = self.pdb_list.retrieve_pdb_file(
                        pdb_id, pdir=self.structures_dir, file_format="pdb"
                    )
                    
                    # Rename to consistent format
                    if os.path.exists(retrieved_path) and retrieved_path != pdb_path:
                        if os.path.exists(pdb_path):
                            os.remove(pdb_path)
                        os.rename(retrieved_path, pdb_path)
                else:
                    # Fallback: direct download
                    pdb_url = f"{self.files_base}/{pdb_id}.pdb"
                    self._download_file(pdb_url, pdb_path)
            
            if os.path.exists(pdb_path):
                files['pdb'] = pdb_path

            else:
                raise FileNotFoundError(f"Failed to download PDB file for {pdb_id}")
                
        except Exception as e:
            raise
        
        # Fetch mmCIF file if requested
        if include_mmcif:
            try:
                mmcif_path = os.path.join(self.structures_dir, f"{pdb_id}.cif")
                if not os.path.exists(mmcif_path):
                    mmcif_url = f"{self.files_base}/{pdb_id}.cif"
                    self._download_file(mmcif_url, mmcif_path)
                
                if os.path.exists(mmcif_path):
                    files['mmcif'] = mmcif_path

                    
            except Exception as e:
                pass
        return files
    
    def _fetch_metadata_graphql(self, pdb_id: str) -> Dict[str, Any]:
        """Fetch comprehensive metadata using GraphQL API."""
        metadata = {}
        
        try:
            # Comprehensive entry query with all available metadata
            entry_query = {
                "query": f"""
                {{
                  entry(entry_id: "{pdb_id}") {{
                    rcsb_id
                    struct {{
                      title
                      pdbx_descriptor
                    }}
                    rcsb_entry_info {{
                      resolution_combined
                      experimental_method
                      deposited_polymer_entity_instance_count
                      deposited_nonpolymer_entity_instance_count
                      deposited_atom_count
                      deposited_hydrogen_atom_count
                      molecular_weight
                      polymer_composition
                      selected_polymer_entity_types
                      structure_determination_methodology
                      diffrn_resolution_high {{
                        value
                      }}
                    }}
                    exptl {{
                      method
                      crystals_number
                    }}
                    rcsb_accession_info {{
                      deposit_date
                      initial_release_date
                      revision_date
                      major_revision
                      minor_revision
                    }}
                    audit_author {{
                      name
                      pdbx_ordinal
                    }}
                    cell {{
                      length_a
                      length_b
                      length_c
                      angle_alpha
                      angle_beta
                      angle_gamma
                      volume
                      Z_PDB
                    }}
                    symmetry {{
                      space_group_name_H_M
                      Int_Tables_number
                    }}
                    refine {{
                      ls_R_factor_R_work
                      ls_R_factor_R_free
                      ls_R_factor_obs
                      ls_number_reflns_obs
                      ls_number_reflns_R_free
                      ls_percent_reflns_obs
                      ls_d_res_high
                      ls_d_res_low
                      pdbx_ls_sigma_F
                      pdbx_ls_sigma_I
                      pdbx_R_Free_selection_details
                      pdbx_stereochemistry_target_values
                    }}
                    diffrn {{
                      id
                      ambient_temp
                      crystal_id
                    }}
                    diffrn_detector {{
                      detector
                      type
                    }}
                    diffrn_source {{
                      source
                      type
                      pdbx_synchrotron_site
                      pdbx_synchrotron_beamline
                      pdbx_wavelength_list
                    }}
                    pdbx_database_status {{
                      status_code
                      status_code_cs
                      status_code_mr
                      status_code_sf
                    }}
                    rcsb_primary_citation {{
                      pdbx_database_id_PubMed
                      pdbx_database_id_DOI
                      title
                      journal_abbrev
                      journal_volume
                      page_first
                      page_last
                      year
                      rcsb_authors
                      rcsb_journal_abbrev
                    }}
                    software {{
                      name
                      version
                      classification
                      description
                    }}
                  }}
                }}
                """
            }
            
            response = self._make_graphql_request(entry_query)
            if response and 'data' in response and response['data'] and 'entry' in response['data']:
                entry_data = response['data']['entry']
                if entry_data:
                    metadata['entry'] = entry_data

                else:
                    pass
            else:
                if response and 'errors' in response:
                    pass
            
            # Fetch additional entity details if available
            if 'entry' in metadata:
                entry = metadata['entry']
                
                # Extract polymer entity count
                polymer_count = entry.get('rcsb_entry_info', {}).get('deposited_polymer_entity_instance_count', 0)
                nonpolymer_count = entry.get('rcsb_entry_info', {}).get('deposited_nonpolymer_entity_instance_count', 0)
                
                if polymer_count > 0:
                    metadata['polymer_entity_count'] = polymer_count

                    
                    # Try to fetch detailed polymer entity information
                    polymer_entities = self._fetch_polymer_entities(pdb_id)
                    if polymer_entities:
                        metadata['polymer_entities'] = polymer_entities
                
                if nonpolymer_count > 0:
                    metadata['nonpolymer_entity_count'] = nonpolymer_count

                    
                    # Try to fetch detailed non-polymer entity information
                    nonpolymer_entities = self._fetch_nonpolymer_entities(pdb_id)
                    if nonpolymer_entities:
                        metadata['nonpolymer_entities'] = nonpolymer_entities
                
                # Fetch assembly information
                assemblies = self._fetch_assemblies(pdb_id)
                if assemblies:
                    metadata['assemblies'] = assemblies
            
            # Save metadata to file
            metadata_file = os.path.join(self.metadata_dir, f"{pdb_id}_metadata.json")
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2, default=str)

            
        except Exception as e:
            pass
        return metadata
    
    def _fetch_polymer_entities(self, pdb_id: str) -> List[Dict[str, Any]]:
        """Fetch detailed polymer entity information."""
        try:
            # Query for polymer entities
            polymer_query = {
                "query": f"""
                {{
                  polymer_entities(entry_id: "{pdb_id}") {{
                    rcsb_id
                    rcsb_polymer_entity {{
                      pdbx_description
                      type
                      pdbx_number_of_molecules
                      formula_weight
                      pdbx_mutation
                      pdbx_ec
                    }}
                    entity_poly {{
                      type
                      pdbx_seq_one_letter_code
                      pdbx_seq_one_letter_code_can
                      rcsb_sample_sequence_length
                    }}
                    rcsb_polymer_entity_container_identifiers {{
                      asym_ids
                      auth_asym_ids
                      entity_id
                    }}
                    rcsb_polymer_entity_feature {{
                      type
                      name
                      description
                      feature_positions {{
                        beg_seq_id
                        end_seq_id
                      }}
                    }}
                  }}
                }}
                """
            }
            
            response = self._make_graphql_request(polymer_query)
            if response and 'data' in response and response['data']:
                entities = response['data'].get('polymer_entities', [])
                if entities:

                    return entities
            
        except Exception as e:
            pass
        return []
    
    def _fetch_nonpolymer_entities(self, pdb_id: str) -> List[Dict[str, Any]]:
        """Fetch detailed non-polymer entity information."""
        try:
            # Query for non-polymer entities
            nonpolymer_query = {
                "query": f"""
                {{
                  nonpolymer_entities(entry_id: "{pdb_id}") {{
                    rcsb_id
                    rcsb_nonpolymer_entity {{
                      pdbx_description
                      type
                      formula_weight
                      pdbx_number_of_molecules
                    }}
                    nonpolymer_comp {{
                      chem_comp {{
                        id
                        name
                        type
                        formula
                        formula_weight
                        pdbx_synonyms
                      }}
                    }}
                    rcsb_nonpolymer_entity_container_identifiers {{
                      asym_ids
                      auth_asym_ids
                      entity_id
                    }}
                  }}
                }}
                """
            }
            
            response = self._make_graphql_request(nonpolymer_query)
            if response and 'data' in response and response['data']:
                entities = response['data'].get('nonpolymer_entities', [])
                if entities:

                    return entities
            
        except Exception as e:
            pass
        return []
    
    def _fetch_assemblies(self, pdb_id: str) -> List[Dict[str, Any]]:
        """Fetch assembly information."""
        try:
            # Query for assemblies
            assembly_query = {
                "query": f"""
                {{
                  assemblies(entry_id: "{pdb_id}") {{
                    rcsb_id
                    rcsb_assembly_info {{
                      assembly_id
                      polymer_entity_instance_count
                      nonpolymer_entity_instance_count
                      polymer_monomer_count
                      entry_id
                    }}
                    pdbx_struct_assembly {{
                      details
                      id
                      method_details
                      oligomeric_details
                      oligomeric_count
                    }}
                    pdbx_struct_assembly_gen {{
                      assembly_id
                      asym_id_list
                      oper_expression
                    }}
                  }}
                }}
                """
            }
            
            response = self._make_graphql_request(assembly_query)
            if response and 'data' in response and response['data']:
                assemblies = response['data'].get('assemblies', [])
                if assemblies:

                    return assemblies
            
        except Exception as e:
            pass
        return []
    
    def _fetch_sequences(self, pdb_id: str) -> Dict[str, Any]:
        """Fetch sequence information including FASTA files."""
        sequences = {}
        
        try:
            # Download FASTA file - try multiple URL formats
            fasta_path = os.path.join(self.sequences_dir, f"{pdb_id}.fasta")
            if not os.path.exists(fasta_path):
                # Try different FASTA URL formats
                fasta_urls = [
                    f"{self.files_base}/{pdb_id}.fasta",
                    f"https://www.rcsb.org/fasta/entry/{pdb_id}",
                    f"https://files.rcsb.org/view/{pdb_id}.fasta"
                ]
                
                success = False
                for fasta_url in fasta_urls:
                    success = self._download_file(fasta_url, fasta_path)
                    if success:
                        break
                
                if not success:

                    return sequences
            
            if os.path.exists(fasta_path):
                sequences['fasta_file'] = fasta_path
                
                # Parse FASTA content
                with open(fasta_path, 'r') as f:
                    fasta_content = f.read()
                
                sequences['fasta_content'] = fasta_content
                sequences['chains'] = self._parse_fasta_chains(fasta_content)

            
        except Exception as e:
            pass
        return sequences
    
    def _fetch_validation_data(self, pdb_id: str) -> Dict[str, Any]:
        """Fetch available validation data (limited due to API changes)."""
        validation = {}
        
        try:
            # Try to download validation report PDF if available
            pdf_path = os.path.join(self.validation_dir, f"{pdb_id}_validation.pdf")
            if not os.path.exists(pdf_path):
                pdf_url = f"{self.files_base}/{pdb_id}_validation.pdf"
                success = self._download_file(pdf_url, pdf_path)
                if success:
                    validation['pdf_report'] = pdf_path

                else:
                    pass
        except Exception as e:
            pass
        return validation
    
    def _make_graphql_request(self, query: Dict, timeout: int = 30) -> Optional[Dict]:
        """Make a GraphQL request with improved error handling."""
        try:
            response = requests.post(self.graphql_api, json=query, timeout=timeout)
            
            if response.status_code == 200:
                return response.json()
            else:
                return None
                
        except requests.exceptions.Timeout:
            return None
        except requests.exceptions.RequestException as e:
            return None
        except json.JSONDecodeError as e:
            return None
    
    def _download_file(self, url: str, filepath: str, timeout: int = 60) -> bool:
        """Download a file from URL with improved error handling."""
        try:
            response = requests.get(url, timeout=timeout, stream=True)
            
            if response.status_code == 200:
                with open(filepath, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                return True
            else:
                return False
                
        except Exception as e:
            return False
    
    def _parse_fasta_chains(self, fasta_content: str) -> List[Dict[str, str]]:
        """Parse FASTA content to extract chain information."""
        chains = []
        current_chain = None
        
        for line in fasta_content.split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current_chain:
                    chains.append(current_chain)
                
                # Parse header
                header = line[1:]
                parts = header.split('|')
                
                current_chain = {
                    'header': header,
                    'sequence': '',
                    'chain_id': '',
                    'description': ''
                }
                
                # Extract chain ID and description
                if len(parts) >= 2:
                    current_chain['chain_id'] = parts[0].split('_')[-1] if '_' in parts[0] else ''
                    current_chain['description'] = parts[1] if len(parts) > 1 else ''
                
            elif line and current_chain:
                current_chain['sequence'] += line
        
        if current_chain:
            chains.append(current_chain)
        
        return chains
    
    def _save_data_summary(self, pdb_id: str, data: Dict[str, Any]) -> None:
        """Save a comprehensive data summary."""
        summary_file = os.path.join(self.metadata_dir, f"{pdb_id}_summary.json")
        
        # Create a clean summary without large data
        summary = {
            'pdb_id': pdb_id,
            'fetch_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'files_downloaded': list(data['files'].keys()),
            'metadata_sections': list(data['metadata'].keys()),
            'sequences_available': bool(data['sequences']),
            'validation_available': bool(data['validation']),
            'errors': data['errors']
        }
        
        # Add key metadata if available
        if 'entry' in data['metadata']:
            entry = data['metadata']['entry']
            summary['title'] = entry.get('struct', {}).get('title', 'N/A')
            
            # Handle resolution (can be a list)
            resolution = entry.get('rcsb_entry_info', {}).get('resolution_combined', 'N/A')
            if isinstance(resolution, list) and resolution:
                summary['resolution'] = resolution[0]
            else:
                summary['resolution'] = resolution
            
            summary['experimental_method'] = entry.get('rcsb_entry_info', {}).get('experimental_method', 'N/A')
            summary['deposition_date'] = entry.get('rcsb_accession_info', {}).get('deposit_date', 'N/A')
        
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
    
    def get_structure_info(self, pdb_id: str) -> Dict[str, Any]:
        """Get comprehensive structure information from cached data."""
        pdb_id = pdb_id.upper().strip()
        
        # Check if data exists
        summary_file = os.path.join(self.metadata_dir, f"{pdb_id}_summary.json")
        metadata_file = os.path.join(self.metadata_dir, f"{pdb_id}_metadata.json")
        
        info = {
            'pdb_id': pdb_id,
            'available': False,
            'files': {},
            'summary': {},
            'metadata': {}
        }
        
        # Load summary
        if os.path.exists(summary_file):
            try:
                with open(summary_file, 'r') as f:
                    info['summary'] = json.load(f)
                info['available'] = True
            except Exception as e:
                pass
        # Load full metadata
        if os.path.exists(metadata_file):
            try:
                with open(metadata_file, 'r') as f:
                    info['metadata'] = json.load(f)
            except Exception as e:
                pass
        # Check for files
        pdb_file = os.path.join(self.structures_dir, f"{pdb_id}.pdb")
        if os.path.exists(pdb_file):
            info['files']['pdb'] = pdb_file
        
        fasta_file = os.path.join(self.sequences_dir, f"{pdb_id}.fasta")
        if os.path.exists(fasta_file):
            info['files']['fasta'] = fasta_file
        
        return info
    
    def list_downloaded_structures(self) -> List[str]:
        """List all downloaded PDB structures."""
        structures = []
        
        if os.path.exists(self.structures_dir):
            for filename in os.listdir(self.structures_dir):
                if filename.endswith('.pdb'):
                    pdb_id = filename[:-4].upper()
                    structures.append(pdb_id)
        
        return sorted(structures)


def create_enhanced_pdb_puller(data_dir: str):
    """Create an enhanced PDB puller instance."""
    return EnhancedPDBPuller(data_dir)