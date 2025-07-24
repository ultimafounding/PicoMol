#!/usr/bin/env python3
"""
Enhanced PDB Puller for PicoMol.

This module provides comprehensive PDB data fetching including:
- Structure files (PDB, mmCIF)
- Comprehensive metadata via RCSB PDB REST API
- Experimental details and validation metrics
- Sequence information (FASTA)
- Ligand and heteroatom information
- Related structures and cross-references
- Functional annotations and classifications
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
    """Enhanced PDB data puller with comprehensive information retrieval."""
    
    def __init__(self, data_dir: str):
        """
        Initialize the enhanced PDB puller.
        
        Args:
            data_dir: Directory to store downloaded PDB data
        """
        self.data_dir = data_dir
        self.pdb_list = PDBList() if BIOPYTHON_AVAILABLE else None
        self.pdb_parser = PDBParser(QUIET=True) if BIOPYTHON_AVAILABLE else None
        
        # RCSB PDB API endpoints
        self.api_base = "https://data.rcsb.org/rest/v1/core"
        self.search_api = "https://search.rcsb.org/rcsbsearch/v2/query"
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
        Fetch comprehensive PDB data for a given PDB ID.
        
        Args:
            pdb_id: PDB identifier (e.g., '1CRN')
            include_validation: Whether to fetch validation reports
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
            
            # 2. Fetch comprehensive metadata
            metadata = self._fetch_metadata(pdb_id)
            result['metadata'] = metadata
            
            # 3. Fetch sequence information
            if include_sequences:
                sequences = self._fetch_sequences(pdb_id)
                result['sequences'] = sequences
            
            # 4. Fetch validation data
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
            
            if not os.path.exists(pdb_path) and BIOPYTHON_AVAILABLE:
                # Use Biopython to download PDB file
                retrieved_path = self.pdb_list.retrieve_pdb_file(
                    pdb_id, pdir=self.structures_dir, file_format="pdb"
                )
                
                # Rename to consistent format
                if os.path.exists(retrieved_path) and retrieved_path != pdb_path:
                    if os.path.exists(pdb_path):
                        os.remove(pdb_path)
                    os.rename(retrieved_path, pdb_path)
            
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
    
    def _fetch_metadata(self, pdb_id: str) -> Dict[str, Any]:
        """Fetch comprehensive metadata from RCSB PDB REST API."""
        metadata = {}
        
        try:
            # Fetch entry information
            entry_url = f"{self.api_base}/entry/{pdb_id}"
            entry_data = self._make_api_request(entry_url)
            if entry_data:
                metadata['entry'] = entry_data

            
            # Fetch experimental information
            exptl_url = f"{self.api_base}/experimental_method/{pdb_id}"
            exptl_data = self._make_api_request(exptl_url)
            if exptl_data:
                metadata['experimental'] = exptl_data

            
            # Fetch structure summary
            summary_url = f"{self.api_base}/summary/{pdb_id}"
            summary_data = self._make_api_request(summary_url)
            if summary_data:
                metadata['summary'] = summary_data

            
            # Fetch publication information
            pubmed_url = f"{self.api_base}/pubmed/{pdb_id}"
            pubmed_data = self._make_api_request(pubmed_url)
            if pubmed_data:
                metadata['publications'] = pubmed_data

            
            # Fetch polymer entity information
            polymer_url = f"{self.api_base}/polymer_entity/{pdb_id}"
            polymer_data = self._make_api_request(polymer_url)
            if polymer_data:
                metadata['polymer_entities'] = polymer_data

            
            # Fetch non-polymer entity information (ligands)
            nonpolymer_url = f"{self.api_base}/nonpolymer_entity/{pdb_id}"
            nonpolymer_data = self._make_api_request(nonpolymer_url)
            if nonpolymer_data:
                metadata['nonpolymer_entities'] = nonpolymer_data

            
            # Fetch assembly information
            assembly_url = f"{self.api_base}/assembly/{pdb_id}"
            assembly_data = self._make_api_request(assembly_url)
            if assembly_data:
                metadata['assemblies'] = assembly_data

            
            # Save metadata to file
            metadata_file = os.path.join(self.metadata_dir, f"{pdb_id}_metadata.json")
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2, default=str)

            
        except Exception as e:
            pass
        return metadata
    
    def _fetch_sequences(self, pdb_id: str) -> Dict[str, Any]:
        """Fetch sequence information including FASTA files."""
        sequences = {}
        
        try:
            # Download FASTA file
            fasta_path = os.path.join(self.sequences_dir, f"{pdb_id}.fasta")
            if not os.path.exists(fasta_path):
                fasta_url = f"{self.files_base}/{pdb_id}.fasta"
                self._download_file(fasta_url, fasta_path)
            
            if os.path.exists(fasta_path):
                sequences['fasta_file'] = fasta_path
                
                # Parse FASTA content
                with open(fasta_path, 'r') as f:
                    fasta_content = f.read()
                
                sequences['fasta_content'] = fasta_content
                sequences['chains'] = self._parse_fasta_chains(fasta_content)

            
            # Fetch sequence alignment information
            alignment_url = f"{self.api_base}/sequence_alignment/{pdb_id}"
            alignment_data = self._make_api_request(alignment_url)
            if alignment_data:
                sequences['alignments'] = alignment_data

            
        except Exception as e:
            pass
        return sequences
    
    def _fetch_validation_data(self, pdb_id: str) -> Dict[str, Any]:
        """Fetch validation reports and quality metrics."""
        validation = {}
        
        try:
            # Fetch validation report
            validation_url = f"{self.api_base}/validation/{pdb_id}"
            validation_data = self._make_api_request(validation_url)
            if validation_data:
                validation['report'] = validation_data

            
            # Download validation report PDF if available
            try:
                pdf_path = os.path.join(self.validation_dir, f"{pdb_id}_validation.pdf")
                if not os.path.exists(pdf_path):
                    pdf_url = f"{self.files_base}/{pdb_id}_validation.pdf"
                    self._download_file(pdf_url, pdf_path)
                
                if os.path.exists(pdf_path):
                    validation['pdf_report'] = pdf_path

                    
            except Exception as e:
                pass
            
            # Fetch structure quality metrics
            quality_url = f"{self.api_base}/quality/{pdb_id}"
            quality_data = self._make_api_request(quality_url)
            if quality_data:
                validation['quality_metrics'] = quality_data

            
        except Exception as e:
            pass
        return validation
    
    def _make_api_request(self, url: str, timeout: int = 30) -> Optional[Dict]:
        """Make a request to the RCSB PDB API with error handling."""
        try:
            response = requests.get(url, timeout=timeout)
            
            if response.status_code == 200:
                return response.json()
            elif response.status_code == 404:

                return None
            else:

                return None
                
        except requests.exceptions.Timeout:
            return None
        except requests.exceptions.RequestException as e:
            return None
        except json.JSONDecodeError as e:
            return None
    
    def _download_file(self, url: str, filepath: str, timeout: int = 60) -> bool:
        """Download a file from URL with error handling."""
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
            summary['resolution'] = entry.get('rcsb_entry_info', {}).get('resolution_combined', 'N/A')
            summary['experimental_method'] = entry.get('exptl', [{}])[0].get('method', 'N/A')
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
    
    def get_enhanced_structure_data(self, pdb_id: str) -> Dict[str, Any]:
        """Get enhanced structure data for analysis tools."""
        info = self.get_structure_info(pdb_id)
        
        enhanced_data = {
            'pdb_id': pdb_id,
            'basic_info': {},
            'experimental_info': {},
            'quality_info': {},
            'sequence_info': {},
            'ligand_info': {},
            'files': info['files']
        }
        
        if not info['available']:
            return enhanced_data
        
        metadata = info['metadata']
        
        # Extract basic information
        if 'entry' in metadata:
            entry = metadata['entry']
            enhanced_data['basic_info'] = {
                'title': entry.get('struct', {}).get('title', 'N/A'),
                'authors': [author.get('name', '') for author in entry.get('audit_author', [])],
                'deposition_date': entry.get('rcsb_accession_info', {}).get('deposit_date', 'N/A'),
                'release_date': entry.get('rcsb_accession_info', {}).get('initial_release_date', 'N/A'),
                'revision_date': entry.get('rcsb_accession_info', {}).get('revision_date', 'N/A')
            }
        
        # Extract experimental information
        if 'experimental' in metadata:
            exptl = metadata['experimental']
            enhanced_data['experimental_info'] = {
                'method': exptl.get('method', 'N/A'),
                'resolution': exptl.get('resolution', 'N/A'),
                'r_work': exptl.get('r_work', 'N/A'),
                'r_free': exptl.get('r_free', 'N/A')
            }
        
        # Extract quality information
        if 'validation' in metadata:
            validation = metadata['validation']
            enhanced_data['quality_info'] = validation.get('report', {})
        
        # Extract sequence information
        if 'polymer_entities' in metadata:
            entities = metadata['polymer_entities']
            enhanced_data['sequence_info'] = {
                'entities': entities,
                'chain_count': len(entities) if isinstance(entities, list) else 0
            }
        
        # Extract ligand information
        if 'nonpolymer_entities' in metadata:
            ligands = metadata['nonpolymer_entities']
            enhanced_data['ligand_info'] = {
                'ligands': ligands,
                'ligand_count': len(ligands) if isinstance(ligands, list) else 0
            }
        
        return enhanced_data


def create_enhanced_pdb_puller(data_dir: str) -> EnhancedPDBPuller:
    """Create an enhanced PDB puller instance."""
    return EnhancedPDBPuller(data_dir)