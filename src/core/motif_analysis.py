#!/usr/bin/env python3
"""
Motif and Domain Analysis Tools for PicoMol.

This module provides comprehensive protein motif and domain analysis including:
- InterPro domain annotation
- Pfam domain search
- ExPASy PROSITE motif scanning
- Custom motif pattern matching
- Domain visualization and analysis
"""

import os
import re
import json
import time
import requests
from io import StringIO
from urllib.parse import urlencode, quote
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTabWidget, QGroupBox, QFormLayout,
    QTextEdit, QLineEdit, QPushButton, QLabel, QComboBox, QSpinBox,
    QCheckBox, QTableWidget, QTableWidgetItem, QHeaderView, QScrollArea,
    QMessageBox, QFileDialog, QProgressBar, QSplitter, QDialog, QDialogButtonBox,
    QListWidget, QListWidgetItem, QSplitter as QSplitterWidget, QFrame,
    QGridLayout, QTextBrowser, QApplication
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer, QUrl
from PyQt5.QtGui import QFont, QPixmap, QPainter, QPen, QBrush, QColor, QDesktopServices

try:
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


class InterProSearchWorker(QThread):
    """Worker thread for InterPro domain search."""
    
    search_complete = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)
    progress_update = pyqtSignal(str)
    
    def __init__(self, sequence, email="user@example.com"):
        super().__init__()
        self.sequence = sequence
        self.email = email
        self.base_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
    
    def run(self):
        try:
            self.progress_update.emit("Submitting sequence to InterPro...")
            
            # Submit job to InterPro API
            job_id = self.submit_job()
            if not job_id:
                self.error_occurred.emit("Failed to submit job to InterPro. Please check your internet connection and try again.")
                return
            
            self.progress_update.emit(f"Job submitted (ID: {job_id}). Waiting for results...")
            
            # Poll for results
            results = self.poll_results(job_id)
            if results:
                self.search_complete.emit(results)
            else:
                self.error_occurred.emit("Failed to retrieve results from InterPro. The job may have failed or timed out.")
                
        except Exception as e:
            self.error_occurred.emit(f"InterPro search error: {str(e)}")
    

    
    def submit_job(self):
        """Submit sequence to InterPro."""
        try:
            url = f"{self.base_url}/run"
            
            # Prepare parameters
            params = {
                'email': self.email,
                'sequence': self.sequence,
                'goterms': 'true',
                'pathways': 'true'
            }
            
            print(f"[DEBUG] Submitting to InterPro: {url}")
            print(f"[DEBUG] Sequence length: {len(self.sequence)} amino acids")
            
            response = requests.post(url, data=params, timeout=30)
            
            print(f"[DEBUG] Response status: {response.status_code}")
            if response.status_code == 200:
                job_id = response.text.strip()
                print(f"[DEBUG] Job ID received: {job_id}")
                return job_id
            else:
                print(f"[DEBUG] Failed response: {response.text}")
                return None
                
        except Exception as e:
            print(f"[ERROR] Error submitting InterPro job: {e}")
            return None
    
    def poll_results(self, job_id, max_attempts=60):
        """Poll for job completion and retrieve results."""
        try:
            status_url = f"{self.base_url}/status/{job_id}"
            result_url = f"{self.base_url}/result/{job_id}/json"
            
            for attempt in range(max_attempts):
                # Check status
                status_response = requests.get(status_url, timeout=10)
                if status_response.status_code == 200:
                    status = status_response.text.strip()
                    
                    if status == "FINISHED":
                        self.progress_update.emit("Job completed. Retrieving results...")
                        
                        # Get results
                        result_response = requests.get(result_url, timeout=30)
                        if result_response.status_code == 200:
                            return self.parse_interpro_results(result_response.json(), job_id)
                        else:
                            return None
                    
                    elif status == "RUNNING":
                        self.progress_update.emit(f"Job running... (attempt {attempt + 1}/{max_attempts})")
                        time.sleep(10)  # Wait 10 seconds before next check
                    
                    elif status == "ERROR" or status == "FAILURE":
                        return None
                    
                else:
                    time.sleep(5)
            
            return None
            
        except Exception as e:
            print(f"Error polling InterPro results: {e}")
            return None
    
    def parse_interpro_results(self, json_data, job_id=None):
        """Parse InterPro JSON results."""
        try:
            results = {
                'domains': [],
                'go_terms': [],
                'pathways': [],
                'raw_data': json_data,
                'job_id': job_id  # Add job ID to results
            }
            
            # Debug: Print structure of json_data
            print(f"[DEBUG] InterPro JSON structure: {type(json_data)}")
            if isinstance(json_data, dict):
                print(f"[DEBUG] Top-level keys: {list(json_data.keys())}")
            
            if 'results' in json_data and json_data['results']:
                for result in json_data['results']:
                    if 'matches' in result:
                        for match in result['matches']:
                            domain_info = {
                                'signature': match.get('signature', {}),
                                'locations': match.get('locations', []),
                                'model': match.get('model', ''),
                                'score': match.get('score', 0),
                                'evalue': match.get('evalue', 'N/A')
                            }
                            
                            # Extract signature details
                            sig = domain_info['signature']
                            sig_lib_release = sig.get('signatureLibraryRelease', {}) if sig else {}
                            database_name = sig_lib_release.get('library', '') if sig_lib_release else ''
                            
                            domain_info.update({
                                'accession': sig.get('accession', '') if sig else '',
                                'name': sig.get('name', '') if sig else '',
                                'description': sig.get('description', '') if sig else '',
                                'database': database_name,
                                'interpro': sig.get('entry', {}) if sig else {}
                            })
                            
                            results['domains'].append(domain_info)
                            
                            # Extract GO terms
                            if sig and 'entry' in sig and sig['entry'] and 'goXRefs' in sig['entry']:
                                go_xrefs = sig['entry']['goXRefs']
                                if go_xrefs:  # Check if goXRefs is not None
                                    for go_ref in go_xrefs:
                                        if go_ref:  # Check if go_ref is not None
                                            go_info = {
                                                'id': go_ref.get('identifier', ''),
                                                'name': go_ref.get('name', ''),
                                                'category': go_ref.get('category', '')
                                            }
                                            if go_info not in results['go_terms']:
                                                results['go_terms'].append(go_info)
                            
                            # Extract pathway information
                            if sig and 'entry' in sig and sig['entry'] and 'pathwayXRefs' in sig['entry']:
                                pathway_xrefs = sig['entry']['pathwayXRefs']
                                if pathway_xrefs:  # Check if pathwayXRefs is not None
                                    for pathway_ref in pathway_xrefs:
                                        if pathway_ref:  # Check if pathway_ref is not None
                                            pathway_info = {
                                                'id': pathway_ref.get('identifier', ''),
                                                'name': pathway_ref.get('name', ''),
                                                'database': pathway_ref.get('databaseName', '')
                                            }
                                            if pathway_info not in results['pathways']:
                                                results['pathways'].append(pathway_info)
            
            return results
            
        except Exception as e:
            print(f"Error parsing InterPro results: {e}")
            print(f"[DEBUG] Exception type: {type(e)}")
            print(f"[DEBUG] JSON data type: {type(json_data)}")
            if hasattr(json_data, 'keys'):
                print(f"[DEBUG] JSON keys: {list(json_data.keys())}")
            return None


class PfamSearchWorker(QThread):
    """Worker thread for Pfam domain search using InterPro API."""
    
    search_complete = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)
    progress_update = pyqtSignal(str)
    
    def __init__(self, sequence, email="user@example.com"):
        super().__init__()
        self.sequence = sequence
        self.email = email
        self.base_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
    
    def run(self):
        try:
            self.progress_update.emit("Submitting sequence to Pfam (via InterPro)...")
            
            # Try the main InterPro approach first
            job_id = self.submit_job()
            if not job_id:
                # Try fallback approach
                self.progress_update.emit("Trying alternative Pfam search method...")
                results = self.fallback_pfam_search()
                if results:
                    self.search_complete.emit(results)
                    return
                else:
                    self.error_occurred.emit("Failed to submit job to InterPro for Pfam search")
                    return
            
            self.progress_update.emit(f"Job submitted (ID: {job_id}). Waiting for Pfam results...")
            
            # Poll for results
            results = self.poll_results(job_id)
            if results:
                self.search_complete.emit(results)
            else:
                self.error_occurred.emit("No Pfam domains found or search failed")
                
        except Exception as e:
            self.error_occurred.emit(f"Pfam search error: {str(e)}")
    
    def submit_job(self):
        """Submit sequence to InterPro with Pfam-specific parameters."""
        try:
            url = f"{self.base_url}/run"
            
            # Prepare parameters - search all databases but filter for Pfam later
            # The 'appl' parameter might not work as expected, so we'll filter results
            params = {
                'email': self.email,
                'sequence': self.sequence,
                'goterms': 'false',  # Disable GO terms for faster processing
                'pathways': 'false'  # Disable pathways for faster processing
            }
            
            print(f"[DEBUG] Submitting Pfam job with params: {params}")
            response = requests.post(url, data=params, timeout=30)
            print(f"[DEBUG] Response status: {response.status_code}")
            
            if response.status_code == 200:
                job_id = response.text.strip()
                print(f"[DEBUG] Job ID received: {job_id}")
                return job_id
            else:
                print(f"[DEBUG] Failed response: {response.text}")
                return None
                
        except Exception as e:
            print(f"Error submitting Pfam job: {e}")
            return None
    
    def fallback_pfam_search(self):
        """Fallback Pfam search using a simpler approach."""
        try:
            self.progress_update.emit("Using simplified Pfam search...")
            
            # Create a mock result with a message about the limitation
            results = {
                'domains': [],
                'clans': [],
                'raw_data': {},
                'note': 'Pfam search temporarily unavailable. Please try the InterPro search which includes Pfam domains.'
            }
            
            return results
            
        except Exception as e:
            print(f"Error in fallback Pfam search: {e}")
            return None
    
    def poll_results(self, job_id, max_attempts=60):
        """Poll for job completion and retrieve Pfam results."""
        try:
            status_url = f"{self.base_url}/status/{job_id}"
            result_url = f"{self.base_url}/result/{job_id}/json"
            
            for attempt in range(max_attempts):
                # Check status
                status_response = requests.get(status_url, timeout=10)
                if status_response.status_code == 200:
                    status = status_response.text.strip()
                    
                    if status == "FINISHED":
                        self.progress_update.emit("Job completed. Retrieving Pfam results...")
                        
                        # Get results
                        result_response = requests.get(result_url, timeout=30)
                        if result_response.status_code == 200:
                            return self.parse_pfam_results(result_response.json())
                        else:
                            return None
                    
                    elif status == "RUNNING":
                        self.progress_update.emit(f"Job running... (attempt {attempt + 1}/{max_attempts})")
                        time.sleep(10)  # Wait 10 seconds before next check
                    
                    elif status == "ERROR" or status == "FAILURE":
                        return None
                    
                else:
                    time.sleep(5)
            
            return None
            
        except Exception as e:
            print(f"Error polling Pfam results: {e}")
            return None
    
    def parse_pfam_results(self, json_data):
        """Parse InterPro JSON results and extract Pfam domains only."""
        try:
            results = {
                'domains': [],
                'clans': [],
                'raw_data': json_data
            }
            
            # Debug: Print structure of json_data
            print(f"[DEBUG] Pfam JSON structure: {type(json_data)}")
            if isinstance(json_data, dict):
                print(f"[DEBUG] Top-level keys: {list(json_data.keys())}")
            
            if 'results' in json_data and json_data['results']:
                for result in json_data['results']:
                    if 'matches' in result:
                        for match in result['matches']:
                            signature = match.get('signature', {})
                            
                            # Only process Pfam entries
                            sig_lib = signature.get('signatureLibraryRelease', {})
                            library_name = sig_lib.get('library', '') if sig_lib else ''
                            if library_name and library_name.lower() == 'pfam':
                                
                                domain_info = {
                                    'accession': signature.get('accession', ''),
                                    'name': signature.get('name', ''),
                                    'description': signature.get('description', ''),
                                    'type': signature.get('entry', {}).get('type', 'family'),
                                    'evalue': match.get('evalue', 'N/A'),
                                    'score': match.get('score', 0)
                                }
                                
                                # Extract location information
                                locations = match.get('locations', [])
                                if locations:
                                    location = locations[0]  # Take first location
                                    domain_info['start'] = location.get('start', 0)
                                    domain_info['end'] = location.get('end', 0)
                                else:
                                    domain_info['start'] = 0
                                    domain_info['end'] = 0
                                
                                # Extract clan information if available
                                entry_info = signature.get('entry', {}) if signature else {}
                                if entry_info and 'memberDatabases' in entry_info and entry_info['memberDatabases']:
                                    for member_db in entry_info['memberDatabases']:
                                        if member_db and member_db.get('library') == 'pfam':
                                            clan_info = member_db.get('clan', '')
                                            domain_info['clan'] = clan_info
                                            if clan_info and clan_info not in results['clans']:
                                                results['clans'].append(clan_info)
                                            break
                                else:
                                    domain_info['clan'] = ''
                                
                                results['domains'].append(domain_info)
            
            return results
            
        except Exception as e:
            print(f"Error parsing Pfam results: {e}")
            print(f"[DEBUG] Exception type: {type(e)}")
            print(f"[DEBUG] JSON data type: {type(json_data)}")
            if hasattr(json_data, 'keys'):
                print(f"[DEBUG] JSON keys: {list(json_data.keys())}")
            return None


class PROSITESearchWorker(QThread):
    """Worker thread for PROSITE motif search using official API."""
    
    search_complete = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)
    progress_update = pyqtSignal(str)
    
    def __init__(self, sequence, exclude_high_prob=True):
        super().__init__()
        self.sequence = sequence
        self.exclude_high_prob = exclude_high_prob
    
    def run(self):
        try:
            # Try API first with a short timeout, then fallback
            self.progress_update.emit("Trying PROSITE API...")
            
            try:
                # Quick API attempt with short timeout
                job_id = self.submit_prosite_job_quick()
                if job_id and isinstance(job_id, dict):
                    # Got direct results
                    self.search_complete.emit(job_id)
                    return
                elif job_id:
                    # Got job ID, try to poll briefly
                    self.progress_update.emit(f"Polling PROSITE job {job_id}...")
                    results = self.poll_prosite_results_quick(job_id)
                    if results:
                        self.search_complete.emit(results)
                        return
            except Exception as api_error:
                print(f"[DEBUG] API failed: {str(api_error)}")
            
            # API failed or timed out, use local patterns
            self.progress_update.emit("Using local PROSITE patterns...")
            results = self.fallback_local_search()
            if results:
                self.search_complete.emit(results)
            else:
                self.error_occurred.emit("Both API and local search failed")
                
        except Exception as e:
            print(f"[ERROR] PROSITE search error: {str(e)}")
            # Last resort - try local search
            try:
                results = self.fallback_local_search()
                if results:
                    self.search_complete.emit(results)
                else:
                    self.error_occurred.emit(f"PROSITE search failed: {str(e)}")
            except Exception as fallback_error:
                self.error_occurred.emit(f"All PROSITE search methods failed: {str(e)}")
    
    def submit_prosite_job_quick(self):
        """Submit sequence to PROSITE ScanProsite API."""
        import requests
        
        # PROSITE ScanProsite API endpoint
        url = "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"
        
        # Prepare the sequence in FASTA format
        fasta_sequence = f">Query_Sequence\n{self.sequence}"
        
        # Prepare form data - try simpler approach first
        data = {
            'seq': fasta_sequence,
            'output': 'txt',  # Try text output first
            'skip': '1' if self.exclude_high_prob else '0',  # Skip high-probability motifs
        }
        
        print(f"[DEBUG] Submitting to PROSITE: {url}")
        print(f"[DEBUG] Sequence length: {len(self.sequence)} amino acids")
        print(f"[DEBUG] Skip high-prob: {self.exclude_high_prob}")
        
        try:
            response = requests.post(url, data=data, timeout=15)  # Shorter timeout
            response.raise_for_status()
            
            content = response.text
            print(f"[DEBUG] Response length: {len(content)} characters")
            print(f"[DEBUG] Response preview: {content[:500]}...")
            
            # Check if we got direct results or need to poll
            if 'PS' in content and ('match' in content.lower() or 'found' in content.lower()):
                # Direct results - parse immediately
                print("[DEBUG] Got direct results")
                return self.parse_prosite_response(content)
            
            # Look for job ID in the response
            import re
            job_match = re.search(r'job[_\s]*id[^\w]*([\w\-]+)', content, re.IGNORECASE)
            if job_match:
                job_id = job_match.group(1)
                print(f"[DEBUG] Found job ID: {job_id}")
                return job_id
            
            # If no clear job ID, try to parse as results anyway
            print("[DEBUG] No job ID found, trying to parse as results")
            return self.parse_prosite_response(content)
            
        except requests.exceptions.RequestException as e:
            print(f"[ERROR] Request failed: {str(e)}")
            return None  # Return None instead of raising exception
        except Exception as e:
            print(f"[ERROR] Processing failed: {str(e)}")
            return None  # Return None instead of raising exception
    
    def poll_prosite_results_quick(self, job_id):
        """Quick poll PROSITE API for job results with short timeout."""
        import requests
        import time
        
        # If job_id is actually parsed results, return them
        if isinstance(job_id, dict):
            return job_id
        
        status_url = f"https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"
        
        max_attempts = 5  # Much shorter polling
        attempt = 0
        
        while attempt < max_attempts:
            try:
                # Check job status
                response = requests.get(f"{status_url}?jobid={job_id}&output=txt", timeout=10)
                response.raise_for_status()
                
                content = response.text
                
                # Check if results are ready
                if 'running' in content.lower() or 'pending' in content.lower():
                    attempt += 1
                    time.sleep(1)  # Wait 1 second before next poll
                    continue
                
                # Parse results
                return self.parse_prosite_response(content)
                    
            except requests.exceptions.RequestException as e:
                attempt += 1
                if attempt >= max_attempts:
                    return None
                time.sleep(1)
        
        return None  # Timeout
    
    def poll_prosite_results(self, job_id):
        """Poll PROSITE API for job results."""
        import requests
        import time
        
        # If job_id is actually parsed results, return them
        if isinstance(job_id, dict):
            return job_id
        
        status_url = f"https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"
        
        max_attempts = 30  # Maximum polling attempts
        attempt = 0
        
        while attempt < max_attempts:
            try:
                # Check job status
                response = requests.get(f"{status_url}?jobid={job_id}&output=xml", timeout=30)
                response.raise_for_status()
                
                content = response.text
                
                # Check if results are ready
                if 'running' in content.lower() or 'pending' in content.lower():
                    attempt += 1
                    self.progress_update.emit(f"Waiting for results... (attempt {attempt}/{max_attempts})")
                    time.sleep(2)  # Wait 2 seconds before next poll
                    continue
                
                # Parse results
                if content.startswith('<?xml') or '<prositeresults>' in content.lower():
                    return self.parse_prosite_xml(content)
                else:
                    return self.parse_prosite_response(content)
                    
            except requests.exceptions.RequestException as e:
                attempt += 1
                if attempt >= max_attempts:
                    raise Exception(f"Failed to retrieve results: {str(e)}")
                time.sleep(2)
        
        raise Exception("Timeout waiting for PROSITE results")
    
    def parse_prosite_xml(self, xml_content):
        """Parse PROSITE XML results."""
        try:
            import xml.etree.ElementTree as ET
            
            root = ET.fromstring(xml_content)
            
            results = {
                'motifs': [],
                'pattern_count': 0,
                'categories': {},
                'sequence_length': len(self.sequence),
                'search_date': time.strftime('%Y-%m-%d %H:%M:%S'),
                'source': 'PROSITE API'
            }
            
            # Parse matches from XML
            for match in root.findall('.//match'):
                motif_info = {
                    'id': match.get('id', ''),
                    'name': match.get('name', ''),
                    'description': match.get('description', ''),
                    'start': int(match.get('start', 0)),
                    'end': int(match.get('end', 0)),
                    'sequence': match.get('sequence', ''),
                    'length': int(match.get('end', 0)) - int(match.get('start', 0)) + 1,
                    'category': self.categorize_motif(match.get('id', '')),
                    'function': match.get('function', ''),
                    'skip_flag': self.is_high_probability_motif(match.get('id', ''))
                }
                
                results['motifs'].append(motif_info)
                results['pattern_count'] += 1
                
                # Count by category
                category = motif_info['category']
                if category not in results['categories']:
                    results['categories'][category] = 0
                results['categories'][category] += 1
            
            return results
            
        except Exception as e:
            # Fallback to text parsing
            return self.parse_prosite_response(xml_content)
    
    def parse_prosite_response(self, content):
        """Parse PROSITE text/HTML response."""
        results = {
            'motifs': [],
            'pattern_count': 0,
            'categories': {},
            'sequence_length': len(self.sequence),
            'search_date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'source': 'PROSITE API'
        }
        
        print(f"[DEBUG] Parsing PROSITE response...")
        
        # Check if this is an HTML page with results
        if '<html' in content.lower() or '<body' in content.lower():
            return self.parse_prosite_html(content, results)
        
        # Try to parse as plain text (tab-separated format)
        lines = content.split('\n')
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            # Look for tab-separated PROSITE results
            # Format: Query_Sequence\tstart\tend\tPS_ID\t.\t.\t.\tsequence
            if '\t' in line and 'PS' in line:
                print(f"[DEBUG] Found potential match line: {line}")
                
                # Try to extract information from this tab-separated line
                motif_info = self.extract_motif_from_tab_line(line)
                if motif_info:
                    results['motifs'].append(motif_info)
                    results['pattern_count'] += 1
                    
                    # Count by category
                    category = motif_info['category']
                    if category not in results['categories']:
                        results['categories'][category] = 0
                    results['categories'][category] += 1
            
            # Also check for old format (lines starting with PS)
            elif line.startswith('PS') and any(char.isdigit() for char in line):
                print(f"[DEBUG] Found potential PS line: {line}")
                
                # Try to extract information from this line and surrounding lines
                motif_info = self.extract_motif_from_line(line, lines, i)
                if motif_info:
                    results['motifs'].append(motif_info)
                    results['pattern_count'] += 1
                    
                    # Count by category
                    category = motif_info['category']
                    if category not in results['categories']:
                        results['categories'][category] = 0
                    results['categories'][category] += 1
        
        print(f"[DEBUG] Found {len(results['motifs'])} motifs")
        
        # Debug: Print first few motifs
        for i, motif in enumerate(results['motifs'][:3]):
            print(f"[DEBUG] Motif {i+1}: {motif['id']} - {motif['name']} at {motif['start']}-{motif['end']}")
        
        return results
    
    def parse_prosite_html(self, content, results):
        """Parse HTML response from PROSITE."""
        print(f"[DEBUG] Parsing HTML response")
        
        # Look for table rows or other structured data
        import re
        
        # Try to find PS patterns in the HTML
        ps_pattern = re.compile(r'PS\d{5}')
        ps_matches = ps_pattern.findall(content)
        
        if ps_matches:
            print(f"[DEBUG] Found PS patterns in HTML: {ps_matches}")
            
            for ps_id in set(ps_matches):  # Remove duplicates
                # Create a basic motif entry
                motif_info = {
                    'id': ps_id,
                    'name': f'PROSITE pattern {ps_id}',
                    'description': f'Pattern {ps_id} found in sequence',
                    'start': 1,
                    'end': len(self.sequence),
                    'sequence': 'N/A',
                    'length': 1,
                    'category': self.categorize_motif(ps_id),
                    'function': 'Pattern match',
                    'source': 'PROSITE API'
                }
                
                results['motifs'].append(motif_info)
                results['pattern_count'] += 1
                
                # Count by category
                category = motif_info['category']
                if category not in results['categories']:
                    results['categories'][category] = 0
                results['categories'][category] += 1
        
        return results
    
    def extract_motif_from_line(self, line, lines, line_index):
        """Extract motif information from a line and its context."""
        import re
        
        # Try to extract PS ID
        ps_match = re.search(r'(PS\d{5})', line)
        if not ps_match:
            return None
        
        ps_id = ps_match.group(1)
        
        # Try to extract position information
        pos_match = re.search(r'(\d+)\s*-\s*(\d+)', line)
        if pos_match:
            start = int(pos_match.group(1))
            end = int(pos_match.group(2))
        else:
            start = 1
            end = len(self.sequence)
        
        # Try to extract sequence match
        seq_match = re.search(r'([A-Z]{3,})', line)
        sequence = seq_match.group(1) if seq_match else 'N/A'
        
        motif_info = {
            'id': ps_id,
            'name': f'PROSITE pattern {ps_id}',
            'description': f'Pattern {ps_id} found in sequence',
            'start': start,
            'end': end,
            'sequence': sequence,
            'length': end - start + 1,
            'category': self.categorize_motif(ps_id),
            'function': 'Pattern match',
            'source': 'PROSITE API'
        }
        
        return motif_info
    
    def extract_motif_from_tab_line(self, line):
        """Extract motif information from a tab-separated line."""
        try:
            # Split by tabs
            parts = line.split('\t')
            
            # Expected format: Query_Sequence\tstart\tend\tPS_ID\t.\t.\t.\tsequence
            if len(parts) < 8:
                print(f"[DEBUG] Line has {len(parts)} parts, expected at least 8")
                return None
            
            sequence_name = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            ps_id = parts[3]
            sequence_match = parts[7] if len(parts) > 7 else 'N/A'
            
            # Validate PS ID format
            if not ps_id.startswith('PS') or not any(char.isdigit() for char in ps_id):
                print(f"[DEBUG] Invalid PS ID: {ps_id}")
                return None
            
            # Apply filtering if enabled
            if self.exclude_high_prob and self.is_high_probability_motif(ps_id):
                print(f"[DEBUG] Skipping high-probability motif: {ps_id}")
                return None
            
            motif_info = {
                'id': ps_id,
                'name': self.get_motif_name(ps_id),
                'description': f'PROSITE pattern {ps_id}',
                'start': start,
                'end': end,
                'sequence': sequence_match,
                'length': end - start + 1,
                'category': self.categorize_motif(ps_id),
                'function': 'Pattern match',
                'source': 'PROSITE API'
            }
            
            print(f"[DEBUG] Successfully parsed motif: {ps_id} at {start}-{end}")
            return motif_info
            
        except (ValueError, IndexError) as e:
            print(f"[DEBUG] Error parsing tab line: {e}")
            return None
    
    def get_motif_name(self, ps_id):
        """Get a descriptive name for a PROSITE pattern ID."""
        # Map common PROSITE IDs to their names
        motif_names = {
            'PS00001': 'N-glycosylation site',
            'PS00002': 'cAMP/cGMP-dependent protein kinase',
            'PS00003': 'Casein kinase II phosphorylation site',
            'PS00004': 'Protein kinase C phosphorylation site',
            'PS00005': 'Tyrosine kinase phosphorylation site',
            'PS00006': 'ck2 Phosphorylation site',
            'PS00007': 'Tyrosine sulfation site',
            'PS00008': 'Myristoylation site',
            'PS00009': 'Amidation site',
            'PS00017': 'ATP/GTP-binding site motif A (P-loop)',
            'PS00028': 'Zinc finger C2H2 type',
            'PS00134': 'Serine protease active site',
            'PS00213': 'Leucine zipper pattern',
            'PS00237': 'G-protein coupled receptor signature',
            'PS00239': 'Signal peptide',
            'PS00241': 'EF-hand calcium binding'
        }
        
        return motif_names.get(ps_id, f'PROSITE pattern {ps_id}')
    
    def categorize_motif(self, motif_id):
        """Categorize motif based on ID."""
        # Basic categorization based on known PROSITE patterns
        if motif_id.startswith('PS0000') or motif_id in ['PS00001', 'PS00002', 'PS00003', 'PS00004', 'PS00005', 'PS00006', 'PS00007', 'PS00008', 'PS00009']:
            return 'Post-translational modification'
        elif motif_id.startswith('PS0001') or motif_id.startswith('PS0002'):
            return 'Nucleotide binding'
        elif motif_id.startswith('PS0003'):
            return 'DNA/RNA binding'
        elif motif_id.startswith('PS0013') or motif_id.startswith('PS0014'):
            return 'Enzymatic site'
        elif motif_id.startswith('PS0023'):
            return 'Membrane protein'
        else:
            return 'Other'
    
    def is_high_probability_motif(self, motif_id):
        """Check if motif has high probability of occurrence (PS00001-PS00009)."""
        high_prob_patterns = ['PS00001', 'PS00002', 'PS00003', 'PS00004', 'PS00005', 'PS00006', 'PS00007', 'PS00008', 'PS00009']
        return motif_id in high_prob_patterns
    
    def fallback_local_search(self):
        """Fallback to local pattern matching if API fails."""
        print("[DEBUG] Using fallback local PROSITE search")
        
        # Comprehensive PROSITE patterns for local search
        patterns = {
            # High-probability patterns (PS00001-PS00009) - excluded by default
            'PS00001': {
                'name': 'N-glycosylation site',
                'pattern': r'N[^P][ST][^P]',
                'description': 'Asparagine-linked glycosylation motif',
                'category': 'Post-translational modification',
                'function': 'Protein glycosylation',
                'skip_flag': True
            },
            'PS00002': {
                'name': 'cAMP/cGMP-dependent protein kinase',
                'pattern': r'[RK].{2}[ST]',
                'description': 'Protein kinase phosphorylation site',
                'category': 'Post-translational modification',
                'function': 'Protein phosphorylation',
                'skip_flag': True
            },
            'PS00003': {
                'name': 'Casein kinase II phosphorylation site',
                'pattern': r'[ST].{2}[DE]',
                'description': 'Casein kinase II phosphorylation site',
                'category': 'Post-translational modification',
                'function': 'Protein phosphorylation',
                'skip_flag': True
            },
            'PS00004': {
                'name': 'Protein kinase C phosphorylation site',
                'pattern': r'[ST][RK]',
                'description': 'Protein kinase C phosphorylation site',
                'category': 'Post-translational modification',
                'function': 'Protein phosphorylation',
                'skip_flag': True
            },
            'PS00008': {
                'name': 'Myristoylation site',
                'pattern': r'^MG[^EDHPFYW].{2}[STAGCN][^P]',
                'description': 'N-myristoylation site',
                'category': 'Post-translational modification',
                'function': 'Membrane anchoring',
                'skip_flag': True
            },
            
            # Specific functional patterns (not excluded)
            'PS00017': {
                'name': 'ATP/GTP-binding site motif A (P-loop)',
                'pattern': r'[AG].{4}GK[ST]',
                'description': 'P-loop containing nucleoside triphosphate hydrolases',
                'category': 'Nucleotide binding',
                'function': 'ATP/GTP binding and hydrolysis',
                'skip_flag': False
            },
            'PS00028': {
                'name': 'Zinc finger C2H2 type',
                'pattern': r'C.{2,4}C.{12}H.{3,5}H',
                'description': 'Zinc finger C2H2 type domain',
                'category': 'DNA/RNA binding',
                'function': 'DNA binding',
                'skip_flag': False
            },
            'PS00134': {
                'name': 'Serine protease active site',
                'pattern': r'[LIVM][LIVM]GG[DE][STA][GA][YF]',
                'description': 'Serine protease active site',
                'category': 'Enzymatic site',
                'function': 'Proteolytic activity',
                'skip_flag': False
            },
            'PS00213': {
                'name': 'Leucine zipper pattern',
                'pattern': r'L.{6}L.{6}L.{6}L',
                'description': 'Leucine zipper motif',
                'category': 'Protein-protein interaction',
                'function': 'Protein dimerization',
                'skip_flag': False
            },
            'PS00237': {
                'name': 'G-protein coupled receptor signature',
                'pattern': r'[GSTAIV][NEQHRK][^P][^FYWP][LIVMFYW][^P]DRY',
                'description': 'G-protein coupled receptor signature',
                'category': 'Membrane protein',
                'function': 'Signal transduction',
                'skip_flag': False
            },
            'PS00239': {
                'name': 'Signal peptide',
                'pattern': r'^M[LIVMFYW]{5,15}[LIVMFYWACGS]{5,15}[AGSV][AGSV][AGSV]',
                'description': 'Signal peptide cleavage site',
                'category': 'Signal sequence',
                'function': 'Protein targeting',
                'skip_flag': False
            },
            'PS00241': {
                'name': 'EF-hand calcium binding',
                'pattern': r'D.{1}[DNS][ILVFYW][DNES]G[LIVMFYW][LIVST][LIVMFYW]E',
                'description': 'EF-hand calcium-binding motif',
                'category': 'Structural motif',
                'function': 'Calcium binding',
                'skip_flag': False
            }
        }
        
        results = {
            'motifs': [],
            'pattern_count': 0,
            'categories': {},
            'sequence_length': len(self.sequence),
            'search_date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'source': 'Local PROSITE patterns'
        }
        
        for pattern_id, pattern_info in patterns.items():
            # Apply filtering based on SKIP-FLAG if enabled
            if self.exclude_high_prob and pattern_info.get('skip_flag', False):
                continue
            
            matches = list(re.finditer(pattern_info['pattern'], self.sequence, re.IGNORECASE))
            
            for match in matches:
                motif_info = {
                    'id': pattern_id,
                    'name': pattern_info['name'],
                    'description': pattern_info['description'],
                    'category': pattern_info['category'],
                    'function': pattern_info['function'],
                    'start': match.start() + 1,  # 1-based indexing
                    'end': match.end(),
                    'sequence': match.group(),
                    'length': len(match.group()),
                    'source': 'Local PROSITE patterns'
                }
                results['motifs'].append(motif_info)
                results['pattern_count'] += 1
                
                # Count by category
                category = pattern_info['category']
                if category not in results['categories']:
                    results['categories'][category] = 0
                results['categories'][category] += 1
        
        print(f"[DEBUG] Fallback search found {len(results['motifs'])} motifs")
        return results


class DomainVisualizationWidget(QWidget):
    """Widget for visualizing protein domains and motifs."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.sequence_length = 0
        self.domains = []
        self.motifs = []
        self.setMinimumHeight(200)  # Base minimum height
        self.setMinimumWidth(600)   # Base minimum width
        self.setStyleSheet("background-color: white; border: 1px solid #ccc; border-radius: 5px;")
        self.dynamic_height = 200  # Will be calculated based on content
        self.dynamic_width = 600   # Will be calculated based on content
    
    def set_data(self, sequence_length, domains=None, motifs=None):
        """Set the data for visualization."""
        self.sequence_length = sequence_length
        self.domains = domains or []
        self.motifs = motifs or []
        
        # Calculate required dimensions based on content
        self.calculate_dynamic_dimensions()
        
        # Update the widget size
        self.setMinimumHeight(self.dynamic_height)
        self.setMaximumHeight(self.dynamic_height + 100)  # Allow some flexibility
        self.setMinimumWidth(self.dynamic_width)
        self.setMaximumWidth(self.dynamic_width + 200)    # Allow some flexibility
        
        self.update()
    
    def calculate_dynamic_dimensions(self):
        """Calculate the required width and height based on content amount."""
        # Calculate dynamic width based on sequence length and label needs
        self.calculate_dynamic_width()
        
        # Calculate dynamic height based on number of tracks needed
        self.calculate_dynamic_height()
    
    def calculate_dynamic_width(self):
        """Calculate the required width based on sequence length and content."""
        if self.sequence_length == 0:
            self.dynamic_width = 600
            return
        
        # Base width calculation
        base_width = 200  # Margins and basic elements
        
        # Calculate width needed for sequence visualization
        # Minimum 2 pixels per amino acid, more for shorter sequences
        if self.sequence_length <= 100:
            pixels_per_aa = 8  # Very detailed for short sequences
        elif self.sequence_length <= 300:
            pixels_per_aa = 4  # Good detail for medium sequences
        elif self.sequence_length <= 1000:
            pixels_per_aa = 2  # Reasonable for long sequences
        else:
            pixels_per_aa = 1  # Compressed for very long sequences
        
        sequence_width = self.sequence_length * pixels_per_aa
        
        # Calculate width needed for labels
        max_label_width = 0
        
        # Check domain names
        for domain in self.domains:
            name = domain.get('name', '')
            # Estimate 8 pixels per character for large fonts
            label_width = len(name) * 8
            max_label_width = max(max_label_width, label_width)
        
        # Check motif names
        for motif in self.motifs:
            name = motif.get('name', '')
            label_width = len(name) * 8
            max_label_width = max(max_label_width, label_width)
        
        # Ensure we have enough space for the longest label
        min_width_for_labels = max_label_width + 100  # Extra space for positioning
        
        # Take the maximum of sequence-based width and label-based width
        content_width = max(sequence_width, min_width_for_labels)
        
        self.dynamic_width = max(600, base_width + content_width)
        
        # Cap maximum width to prevent excessive size
        self.dynamic_width = min(self.dynamic_width, 2000)
    
    def calculate_dynamic_height(self):
        """Calculate the required height based on content amount."""
        base_height = 120  # Base space for title, sequence line, scale, legend
        
        # Use the calculated dynamic width for track assignment
        working_width = self.dynamic_width - 160  # Account for margins
        
        # Calculate domain tracks needed
        domain_tracks = 0
        if self.domains:
            temp_tracks = self.assign_to_tracks(self.domains, working_width, 80)
            domain_tracks = len(temp_tracks)
        
        # Calculate motif tracks needed
        motif_tracks = 0
        if self.motifs:
            temp_tracks = self.assign_to_tracks(self.motifs, working_width, 80, min_spacing=40)
            motif_tracks = len(temp_tracks)
        
        # Calculate space needed - MUCH larger for collision-free labels
        domain_space = domain_tracks * 90 + 150   # 90px per track + extra for labels above
        motif_space = motif_tracks * 90 + 150     # 90px per track + extra for labels below
        
        # Add extra space for labels and margins
        label_space = 200  # Much more space for collision-free labels
        margin_space = 200  # Larger top and bottom margins
        
        self.dynamic_height = max(300, base_height + domain_space + motif_space + label_space + margin_space)
        
        # Cap maximum height to prevent excessive size
        self.dynamic_height = min(self.dynamic_height, 1200)
    
    def paintEvent(self, event):
        """Paint the domain visualization with advanced overlap prevention."""
        if self.sequence_length == 0:
            painter = QPainter(self)
            painter.drawText(20, 50, "No sequence data available for visualization")
            return
        
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        
        # Get widget dimensions with better margins
        margin = 80
        width = self.width() - 2 * margin
        height = self.height() - 2 * margin
        
        if width <= 0 or height <= 0:
            return
        
        # Calculate layout positions
        seq_y = margin + height // 2  # Center the sequence line
        domain_area_height = 120  # Space for domains above sequence
        motif_area_height = 120   # Space for motifs below sequence
        
        # Draw title
        painter.setPen(QPen(QColor(0, 0, 0)))
        painter.setFont(QFont("Arial", 14, QFont.Bold))
        title_text = f"Protein Sequence Visualization ({self.sequence_length} amino acids)"
        painter.drawText(margin, margin - 25, title_text)
        
        # Draw sequence line
        painter.setPen(QPen(QColor(0, 0, 0), 4))
        painter.drawLine(margin, seq_y, margin + width, seq_y)
        
        # Draw scale with better spacing
        scale_positions = [0, self.sequence_length // 4, self.sequence_length // 2, 
                          3 * self.sequence_length // 4, self.sequence_length]
        
        painter.setFont(QFont("Arial", 11))
        for pos in scale_positions:
            x = int(margin + (pos / self.sequence_length) * width)
            painter.setPen(QPen(QColor(0, 0, 0), 2))
            painter.drawLine(x, seq_y - 10, x, seq_y + 10)
            
            # Draw scale labels below the line
            painter.setPen(QPen(QColor(80, 80, 80)))
            text_rect = painter.fontMetrics().boundingRect(str(pos))
            painter.drawText(x - text_rect.width() // 2, seq_y + 30, str(pos))
        
        # Draw domains with intelligent layout to prevent overlap
        if self.domains:
            self.draw_domains_with_layout(painter, margin, width, seq_y, domain_area_height)
        
        # Draw motifs with intelligent layout to prevent overlap
        if self.motifs:
            self.draw_motifs_with_layout(painter, margin, width, seq_y, motif_area_height)
        
        # Draw legend
        self.draw_legend(painter, margin)
    
    def draw_domains_with_layout(self, painter, margin, width, seq_y, area_height):
        """Draw domains with intelligent layout to prevent overlap."""
        domain_colors = [
            QColor(255, 120, 120), QColor(120, 255, 120), QColor(120, 120, 255),
            QColor(255, 255, 120), QColor(255, 120, 255), QColor(120, 255, 255),
            QColor(255, 180, 120), QColor(180, 255, 120), QColor(180, 120, 255),
            QColor(200, 150, 255), QColor(255, 200, 150), QColor(150, 255, 200)
        ]
        
        # Sort domains by start position
        sorted_domains = sorted(self.domains, key=lambda d: d.get('start', 0))
        
        # Calculate domain positions and assign to tracks to avoid overlap
        domain_tracks = self.assign_to_tracks(sorted_domains, width, margin)
        
        domain_height = 28  # Much larger domains
        track_spacing = 80  # MUCH more spacing to accommodate labels
        
        # Track used label positions to prevent overlap
        used_label_positions = []
        
        for track_idx, track_domains in enumerate(domain_tracks):
            track_y = seq_y - 80 - (track_idx * track_spacing)  # More space above
            
            for domain_idx, domain in enumerate(track_domains):
                color = domain_colors[domain_idx % len(domain_colors)]
                
                start = domain.get('start', 0)
                end = domain.get('end', 0)
                
                if start > 0 and end > 0:
                    x1 = int(margin + (start / self.sequence_length) * width)
                    x2 = int(margin + (end / self.sequence_length) * width)
                    
                    # Ensure minimum width for visibility
                    if x2 - x1 < 5:
                        x2 = x1 + 5
                    
                    # Draw domain rectangle with border
                    painter.setBrush(QBrush(color))
                    painter.setPen(QPen(color.darker(150), 2))
                    painter.drawRect(x1, track_y, x2 - x1, domain_height)
                    
                    # Draw domain label with collision detection
                    self.draw_domain_label_with_collision_detection(painter, domain, x1, x2, track_y, domain_height, used_label_positions)
    
    def draw_motifs_with_layout(self, painter, margin, width, seq_y, area_height):
        """Draw motifs with intelligent layout to prevent overlap."""
        motif_colors = [
            QColor(255, 0, 0), QColor(0, 150, 0), QColor(0, 0, 255),
            QColor(255, 100, 0), QColor(150, 0, 150), QColor(0, 150, 150),
            QColor(200, 100, 0), QColor(100, 200, 0), QColor(100, 0, 200)
        ]
        
        # Sort motifs by start position
        sorted_motifs = sorted(self.motifs, key=lambda m: m.get('start', 0))
        
        # Calculate motif positions and assign to tracks
        motif_tracks = self.assign_to_tracks(sorted_motifs, width, margin, min_spacing=40)
        
        track_spacing = 80  # MUCH more spacing to accommodate labels
        
        # Track used label positions to prevent overlap
        used_label_positions = []
        
        for track_idx, track_motifs in enumerate(motif_tracks):
            track_y = seq_y + 80 + (track_idx * track_spacing)  # More space below
            
            for motif_idx, motif in enumerate(track_motifs):
                color = motif_colors[motif_idx % len(motif_colors)]
                
                start = motif.get('start', 0)
                end = motif.get('end', 0)
                
                if start > 0 and end > 0:
                    x1 = int(margin + (start / self.sequence_length) * width)
                    x2 = int(margin + (end / self.sequence_length) * width)
                    
                    # Ensure minimum width for visibility
                    if x2 - x1 < 3:
                        x2 = x1 + 3
                    
                    # Draw motif line with thickness
                    painter.setPen(QPen(color, 5))
                    painter.drawLine(x1, track_y, x2, track_y)
                    
                    # Draw motif markers
                    painter.setBrush(QBrush(color))
                    painter.setPen(QPen(color.darker(), 2))
                    painter.drawEllipse(x1 - 3, track_y - 3, 6, 6)
                    painter.drawEllipse(x2 - 3, track_y - 3, 6, 6)
                    
                    # Draw motif label with collision detection
                    self.draw_motif_label_with_collision_detection(painter, motif, x1, x2, track_y, used_label_positions)
    
    def assign_to_tracks(self, elements, width, margin, min_spacing=20):
        """Dynamically assign elements to tracks with intelligent collision detection."""
        if not elements:
            return []
        
        # Create element objects with position info
        positioned_elements = []
        for element in elements:
            start = element.get('start', 0)
            end = element.get('end', 0)
            
            if start <= 0 or end <= 0:
                continue
            
            x1 = int(margin + (start / self.sequence_length) * width)
            x2 = int(margin + (end / self.sequence_length) * width)
            
            # Ensure minimum width
            if x2 - x1 < 10:
                x2 = x1 + 10
            
            positioned_elements.append({
                'element': element,
                'x1': x1,
                'x2': x2,
                'width': x2 - x1,
                'track': -1  # Not assigned yet
            })
        
        # Sort by start position for better packing
        positioned_elements.sort(key=lambda x: x['x1'])
        
        # Dynamic track assignment with collision detection
        tracks = []
        
        for pos_elem in positioned_elements:
            best_track = self.find_best_track(pos_elem, tracks, min_spacing)
            
            if best_track is not None:
                tracks[best_track].append(pos_elem)
                pos_elem['track'] = best_track
            else:
                # Create new track
                tracks.append([pos_elem])
                pos_elem['track'] = len(tracks) - 1
        
        # Convert back to element lists
        result_tracks = []
        for track in tracks:
            result_tracks.append([pos_elem['element'] for pos_elem in track])
        
        return result_tracks
    
    def find_best_track(self, pos_elem, tracks, min_spacing):
        """Find the best track for an element using dynamic collision detection."""
        x1, x2 = pos_elem['x1'], pos_elem['x2']
        
        # Try existing tracks first
        for track_idx, track in enumerate(tracks):
            if self.can_fit_in_track_dynamic(x1, x2, track, min_spacing):
                return track_idx
        
        # No existing track works
        return None
    
    def can_fit_in_track_dynamic(self, x1, x2, track, min_spacing):
        """Dynamic collision detection with gap finding."""
        if not track:
            return True
        
        # Get all occupied ranges in this track
        occupied_ranges = []
        for pos_elem in track:
            occupied_ranges.append((pos_elem['x1'] - min_spacing, pos_elem['x2'] + min_spacing))
        
        # Sort ranges by start position
        occupied_ranges.sort()
        
        # Check if our element fits in any gap
        element_width = x2 - x1
        
        # Check before first element
        if occupied_ranges and x2 + min_spacing <= occupied_ranges[0][0]:
            return True
        
        # Check gaps between elements
        for i in range(len(occupied_ranges) - 1):
            gap_start = occupied_ranges[i][1]
            gap_end = occupied_ranges[i + 1][0]
            gap_width = gap_end - gap_start
            
            if gap_width >= element_width + 2 * min_spacing:
                return True
        
        # Check after last element
        if occupied_ranges and x1 >= occupied_ranges[-1][1] + min_spacing:
            return True
        
        return False
    
    def draw_domain_label(self, painter, domain, x1, x2, y, height):
        """Draw domain label with SIMPLE, RELIABLE text rendering."""
        domain_name = domain.get('name', 'Domain')
        
        # Use large, readable font
        painter.setFont(QFont("Arial", 12, QFont.Bold))
        
        # Get EXACT text dimensions
        metrics = painter.fontMetrics()
        text_width = metrics.width(domain_name)
        text_height = metrics.height()
        
        # Position text ABOVE the domain to avoid any fitting issues
        label_x = x1
        label_y = y - 10
        
        # Draw background box with EXACT dimensions
        padding = 6
        box_x = label_x - padding
        box_y = label_y - text_height - padding
        box_width = text_width + (2 * padding)
        box_height = text_height + (2 * padding)
        
        # Draw background
        painter.setBrush(QBrush(QColor(255, 255, 255, 250)))
        painter.setPen(QPen(QColor(100, 100, 100), 2))
        painter.drawRect(box_x, box_y, box_width, box_height)
        
        # Draw text
        painter.setPen(QPen(QColor(0, 0, 0)))
        painter.drawText(label_x, label_y, domain_name)
        
        # Draw database name if available
        db_name = domain.get('database', '')
        if db_name:
            painter.setFont(QFont("Arial", 10))
            db_metrics = painter.fontMetrics()
            db_width = db_metrics.width(db_name)
            db_height = db_metrics.height()
            
            # Position database label above main label
            db_x = x1
            db_y = box_y - 5
            
            # Draw database background
            db_padding = 4
            db_box_x = db_x - db_padding
            db_box_y = db_y - db_height - db_padding
            db_box_width = db_width + (2 * db_padding)
            db_box_height = db_height + (2 * db_padding)
            
            painter.setBrush(QBrush(QColor(230, 230, 230, 250)))
            painter.setPen(QPen(QColor(150, 150, 150), 1))
            painter.drawRect(db_box_x, db_box_y, db_box_width, db_box_height)
            
            painter.setPen(QPen(QColor(80, 80, 80)))
            painter.drawText(db_x, db_y, db_name)
    
    def draw_domain_label_with_collision_detection(self, painter, domain, x1, x2, y, height, used_positions):
        """Draw domain label with AGGRESSIVE collision avoidance."""
        domain_name = domain.get('name', 'Domain')
        
        # Use large, readable font
        painter.setFont(QFont("Arial", 12, QFont.Bold))
        
        # Get EXACT text dimensions
        metrics = painter.fontMetrics()
        text_width = metrics.width(domain_name)
        text_height = metrics.height()
        
        # MUCH more aggressive position testing - try every 20 pixels
        label_x, label_y = self.find_safe_position_for_domain(
            x1, x2, y, height, text_width, text_height, used_positions
        )
        
        # Draw background box with EXACT dimensions
        padding = 8  # Larger padding
        box_x = label_x - padding
        box_y = label_y - text_height - padding
        box_width = text_width + (2 * padding)
        box_height = text_height + (2 * padding)
        
        # Record this position as used with extra margin
        margin = 15
        used_positions.append((box_x - margin, box_y - margin, 
                             box_width + 2*margin, box_height + 2*margin))
        
        # Draw background
        painter.setBrush(QBrush(QColor(255, 255, 255, 250)))
        painter.setPen(QPen(QColor(100, 100, 100), 2))
        painter.drawRect(box_x, box_y, box_width, box_height)
        
        # Draw text
        painter.setPen(QPen(QColor(0, 0, 0)))
        painter.drawText(label_x, label_y, domain_name)
        
        # Draw database name if available
        db_name = domain.get('database', '')
        if db_name:
            painter.setFont(QFont("Arial", 10))
            db_metrics = painter.fontMetrics()
            db_width = db_metrics.width(db_name)
            db_height = db_metrics.height()
            
            # Find safe position for database label
            db_x, db_y = self.find_safe_position_simple(
                label_x, box_y - 25, db_width, db_height, used_positions
            )
            
            # Draw database background
            db_padding = 4
            db_box_x = db_x - db_padding
            db_box_y = db_y - db_height - db_padding
            db_box_width = db_width + (2 * db_padding)
            db_box_height = db_height + (2 * db_padding)
            
            # Record database label position with margin
            used_positions.append((db_box_x - 10, db_box_y - 10, 
                                 db_box_width + 20, db_box_height + 20))
            
            painter.setBrush(QBrush(QColor(230, 230, 230, 250)))
            painter.setPen(QPen(QColor(150, 150, 150), 1))
            painter.drawRect(db_box_x, db_box_y, db_box_width, db_box_height)
            
            painter.setPen(QPen(QColor(80, 80, 80)))
            painter.drawText(db_x, db_y, db_name)
    
    def draw_motif_label(self, painter, motif, x1, x2, y):
        """Draw motif label with SIMPLE, RELIABLE text rendering."""
        motif_name = motif.get('name', 'Motif')
        
        # Use large, readable font
        painter.setFont(QFont("Arial", 11, QFont.Bold))
        
        # Get EXACT text dimensions
        metrics = painter.fontMetrics()
        text_width = metrics.width(motif_name)
        text_height = metrics.height()
        
        # Position text BELOW the motif line
        label_x = x1
        label_y = y + 25
        
        # Draw background box with EXACT dimensions
        padding = 6
        box_x = label_x - padding
        box_y = label_y - text_height - padding
        box_width = text_width + (2 * padding)
        box_height = text_height + (2 * padding)
        
        # Draw background
        painter.setBrush(QBrush(QColor(255, 255, 255, 250)))
        painter.setPen(QPen(QColor(100, 100, 100), 2))
        painter.drawRect(box_x, box_y, box_width, box_height)
        
        # Draw text
        painter.setPen(QPen(QColor(0, 0, 0)))
        painter.drawText(label_x, label_y, motif_name)
        
        # Draw category if available
        category = motif.get('category', '')
        if category:
            painter.setFont(QFont("Arial", 9))
            cat_metrics = painter.fontMetrics()
            cat_width = cat_metrics.width(category)
            cat_height = cat_metrics.height()
            
            # Position category below main label
            cat_x = x1
            cat_y = box_y + box_height + cat_height + 8
            
            # Draw category background
            cat_padding = 4
            cat_box_x = cat_x - cat_padding
            cat_box_y = cat_y - cat_height - cat_padding
            cat_box_width = cat_width + (2 * cat_padding)
            cat_box_height = cat_height + (2 * cat_padding)
            
            painter.setBrush(QBrush(QColor(240, 240, 240, 250)))
            painter.setPen(QPen(QColor(150, 150, 150), 1))
            painter.drawRect(cat_box_x, cat_box_y, cat_box_width, cat_box_height)
            
            painter.setPen(QPen(QColor(80, 80, 80)))
            painter.drawText(cat_x, cat_y, category)
    
    def draw_motif_label_with_collision_detection(self, painter, motif, x1, x2, y, used_positions):
        """Draw motif label with AGGRESSIVE collision avoidance."""
        motif_name = motif.get('name', 'Motif')
        
        # Use large, readable font
        painter.setFont(QFont("Arial", 11, QFont.Bold))
        
        # Get EXACT text dimensions
        metrics = painter.fontMetrics()
        text_width = metrics.width(motif_name)
        text_height = metrics.height()
        
        # MUCH more aggressive position testing
        label_x, label_y = self.find_safe_position_for_motif(
            x1, x2, y, text_width, text_height, used_positions
        )
        
        # Draw background box with EXACT dimensions
        padding = 8  # Larger padding
        box_x = label_x - padding
        box_y = label_y - text_height - padding
        box_width = text_width + (2 * padding)
        box_height = text_height + (2 * padding)
        
        # Record this position as used with extra margin
        margin = 15
        used_positions.append((box_x - margin, box_y - margin, 
                             box_width + 2*margin, box_height + 2*margin))
        
        # Draw background
        painter.setBrush(QBrush(QColor(255, 255, 255, 250)))
        painter.setPen(QPen(QColor(100, 100, 100), 2))
        painter.drawRect(box_x, box_y, box_width, box_height)
        
        # Draw text
        painter.setPen(QPen(QColor(0, 0, 0)))
        painter.drawText(label_x, label_y, motif_name)
        
        # Draw category if available
        category = motif.get('category', '')
        if category:
            painter.setFont(QFont("Arial", 9))
            cat_metrics = painter.fontMetrics()
            cat_width = cat_metrics.width(category)
            cat_height = cat_metrics.height()
            
            # Find safe position for category
            cat_x, cat_y = self.find_safe_position_simple(
                label_x, box_y + box_height + 20, cat_width, cat_height, used_positions
            )
            
            # Draw category background
            cat_padding = 4
            cat_box_x = cat_x - cat_padding
            cat_box_y = cat_y - cat_height - cat_padding
            cat_box_width = cat_width + (2 * cat_padding)
            cat_box_height = cat_height + (2 * cat_padding)
            
            # Record category label position with margin
            used_positions.append((cat_box_x - 10, cat_box_y - 10, 
                                 cat_box_width + 20, cat_box_height + 20))
            
            painter.setBrush(QBrush(QColor(240, 240, 240, 250)))
            painter.setPen(QPen(QColor(150, 150, 150), 1))
            painter.drawRect(cat_box_x, cat_box_y, cat_box_width, cat_box_height)
            
            painter.setPen(QPen(QColor(80, 80, 80)))
            painter.drawText(cat_x, cat_y, category)
    
    def rect_overlaps_any(self, new_rect, existing_rects):
        """Check if a rectangle overlaps with any existing rectangles with safety margin."""
        x1, y1, w1, h1 = new_rect
        
        # Add safety margin to prevent close proximity
        margin = 10
        x1 -= margin
        y1 -= margin
        w1 += 2 * margin
        h1 += 2 * margin
        
        for existing_rect in existing_rects:
            x2, y2, w2, h2 = existing_rect
            
            # Check if rectangles overlap (with margin)
            if not (x1 + w1 <= x2 or x2 + w2 <= x1 or y1 + h1 <= y2 or y2 + h2 <= y1):
                return True
        
        return False
    
    def find_safe_position_for_domain(self, x1, x2, y, height, text_width, text_height, used_positions):
        """Find a safe position for domain label using grid search."""
        # Try positions in a systematic grid pattern
        search_positions = []
        
        # Above domain - multiple heights
        for offset in [15, 35, 55, 75, 95, 115]:
            search_positions.append((x1, y - offset))
            search_positions.append((x1 + (x2-x1)//2 - text_width//2, y - offset))  # Centered
            search_positions.append((x2 - text_width, y - offset))  # Right-aligned
        
        # To the sides
        for y_offset in [0, -20, -40, 20, 40]:
            search_positions.append((x2 + 15, y + height//2 + y_offset))  # Right
            search_positions.append((x1 - text_width - 15, y + height//2 + y_offset))  # Left
        
        # Test each position
        for pos_x, pos_y in search_positions:
            # Ensure position is within reasonable bounds
            if pos_x < 50 or pos_x + text_width > self.width() - 50:
                continue
            if pos_y < 50 or pos_y > self.height() - 50:
                continue
                
            # Check for collision
            padding = 10
            test_rect = (pos_x - padding, pos_y - text_height - padding,
                        text_width + 2*padding, text_height + 2*padding)
            
            if not self.rect_overlaps_any(test_rect, used_positions):
                return pos_x, pos_y
        
        # Fallback - force position far above
        return x1, y - 120
    
    def find_safe_position_for_motif(self, x1, x2, y, text_width, text_height, used_positions):
        """Find a safe position for motif label using grid search."""
        # Try positions in a systematic grid pattern
        search_positions = []
        
        # Below motif - multiple heights
        for offset in [30, 50, 70, 90, 110, 130]:
            search_positions.append((x1, y + offset))
            search_positions.append((x1 + (x2-x1)//2 - text_width//2, y + offset))  # Centered
            search_positions.append((x2 - text_width, y + offset))  # Right-aligned
        
        # To the sides
        for y_offset in [0, 20, 40, -20, -40]:
            search_positions.append((x2 + 15, y + y_offset))  # Right
            search_positions.append((x1 - text_width - 15, y + y_offset))  # Left
        
        # Above motif (last resort)
        for offset in [20, 40, 60]:
            search_positions.append((x1, y - offset))
        
        # Test each position
        for pos_x, pos_y in search_positions:
            # Ensure position is within reasonable bounds
            if pos_x < 50 or pos_x + text_width > self.width() - 50:
                continue
            if pos_y < 50 or pos_y > self.height() - 50:
                continue
                
            # Check for collision
            padding = 10
            test_rect = (pos_x - padding, pos_y - text_height - padding,
                        text_width + 2*padding, text_height + 2*padding)
            
            if not self.rect_overlaps_any(test_rect, used_positions):
                return pos_x, pos_y
        
        # Fallback - force position far below
        return x1, y + 150
    
    def find_safe_position_simple(self, preferred_x, preferred_y, text_width, text_height, used_positions):
        """Find a safe position near the preferred location."""
        # Try positions in expanding circles around preferred position
        for radius in [0, 20, 40, 60, 80, 100]:
            for x_offset in [-radius, 0, radius]:
                for y_offset in [-radius, 0, radius]:
                    if x_offset == 0 and y_offset == 0 and radius > 0:
                        continue
                    
                    test_x = preferred_x + x_offset
                    test_y = preferred_y + y_offset
                    
                    # Bounds check
                    if test_x < 50 or test_x + text_width > self.width() - 50:
                        continue
                    if test_y < 50 or test_y > self.height() - 50:
                        continue
                    
                    # Collision check
                    padding = 10
                    test_rect = (test_x - padding, test_y - text_height - padding,
                                text_width + 2*padding, text_height + 2*padding)
                    
                    if not self.rect_overlaps_any(test_rect, used_positions):
                        return test_x, test_y
        
        # Fallback
        return preferred_x, preferred_y
    
    def draw_legend(self, painter, margin):
        """Draw legend at the bottom."""
        if not (self.domains or self.motifs):
            return
        
        legend_y = self.height() - 40
        painter.setFont(QFont("Arial", 12))
        painter.setPen(QPen(QColor(80, 80, 80)))
        
        legend_parts = []
        if self.domains:
            legend_parts.append("■ Domains (above sequence)")
        if self.motifs:
            legend_parts.append("● Motifs (below sequence)")
        
        legend_text = "Legend: " + " | ".join(legend_parts)
        painter.drawText(margin, legend_y, legend_text)
        
        # Add track info if there are multiple tracks
        total_tracks = 0
        if self.domains:
            domain_tracks = self.assign_to_tracks(self.domains, self.width() - 160, 80)
            total_tracks += len(domain_tracks)
        if self.motifs:
            motif_tracks = self.assign_to_tracks(self.motifs, self.width() - 160, 80)
            total_tracks += len(motif_tracks)
        
        if total_tracks > 2:
            painter.setFont(QFont("Arial", 10))
            painter.setPen(QPen(QColor(120, 120, 120)))
            painter.drawText(margin, legend_y + 15, f"Elements arranged in {total_tracks} tracks to prevent overlap")


class MotifAnalysisTab(QWidget):
    """Main tab for protein motif and domain analysis."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_app = parent
        self.current_sequence = ""
        self.interpro_worker = None
        self.prosite_worker = None
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)
        
        # Input section
        input_group = QGroupBox("Sequence Input")
        input_layout = QVBoxLayout(input_group)
        
        # Sequence input controls
        controls_layout = QHBoxLayout()
        
        # Load from structure button
        load_structure_btn = QPushButton("Load from Current Structure")
        load_structure_btn.setToolTip("Load sequence from the currently displayed protein structure")
        load_structure_btn.clicked.connect(self.load_from_structure)
        controls_layout.addWidget(load_structure_btn)
        
        # Load from file button
        load_file_btn = QPushButton("Load from File...")
        load_file_btn.setToolTip("Load sequence from a FASTA file")
        load_file_btn.clicked.connect(self.load_from_file)
        controls_layout.addWidget(load_file_btn)
        
        controls_layout.addStretch()
        input_layout.addLayout(controls_layout)
        
        # Sequence input area
        self.sequence_input = QTextEdit()
        self.sequence_input.setPlaceholderText(
            "Enter your protein sequence here or use the buttons above to load from structure/file...\n\n"
            "💡 Supports both plain sequences and FASTA format\n"
            "🌐 Connects to InterPro and PROSITE APIs for comprehensive analysis\n\n"
            "Example protein sequence (Human Serum Albumin):\n"
            "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL"
        )
        self.sequence_input.setMaximumHeight(120)
        self.sequence_input.setFont(QFont("Courier", 10))
        input_layout.addWidget(self.sequence_input)
        
        # Analysis buttons
        button_layout = QHBoxLayout()
        
        self.interpro_btn = QPushButton("Search InterPro")
        self.interpro_btn.setToolTip("Search for protein domains and functional sites using InterPro (includes Pfam, SMART, PROSITE, and more)")
        self.interpro_btn.clicked.connect(self.search_interpro)
        button_layout.addWidget(self.interpro_btn)
        
        self.prosite_btn = QPushButton("Search PROSITE")
        self.prosite_btn.setToolTip("Search for protein motifs using PROSITE API with local pattern fallback")
        self.prosite_btn.clicked.connect(self.search_prosite)
        button_layout.addWidget(self.prosite_btn)
        
        self.search_all_btn = QPushButton("Search All Databases")
        self.search_all_btn.setToolTip("Search InterPro and PROSITE simultaneously for comprehensive analysis")
        self.search_all_btn.clicked.connect(self.search_all)
        self.search_all_btn.setStyleSheet("font-weight: bold; background-color: #4CAF50; color: white;")
        button_layout.addWidget(self.search_all_btn)
        
        input_layout.addLayout(button_layout)
        
        # PROSITE filtering options
        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("PROSITE Options:"))
        
        self.exclude_high_prob_checkbox = QCheckBox("Exclude motifs with a high probability of occurrence")
        self.exclude_high_prob_checkbox.setChecked(True)  # Default to excluding high-probability motifs (matches PROSITE website)
        self.exclude_high_prob_checkbox.setToolTip(
            "🎯 Exclude motifs with a high probability of occurrence from the scan\n\n"
            "These are patterns found in many protein sequences, including:\n"
            "• Common post-translational modifications (N-glycosylation, phosphorylation)\n"
            "• Compositionally biased regions\n"
            "• Frequently occurring sequence patterns\n\n"
            "✅ Recommended: Keep checked for cleaner, more specific results\n"
            "(Matches the PROSITE website's SKIP-FLAG=TRUE filtering)"
        )
        filter_layout.addWidget(self.exclude_high_prob_checkbox)
        
        filter_layout.addStretch()
        input_layout.addLayout(filter_layout)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        input_layout.addWidget(self.progress_bar)
        
        # Status label
        self.status_label = QLabel("")
        self.status_label.setStyleSheet("color: #666; font-style: italic;")
        input_layout.addWidget(self.status_label)
        
        layout.addWidget(input_group)
        
        # Results section with tabs
        self.results_tabs = QTabWidget()
        
        # Save all results section (moved to top)
        save_all_widget = QWidget()
        save_all_layout = QHBoxLayout(save_all_widget)
        
        save_all_layout.addWidget(QLabel("Export Results:"))
        save_all_layout.addStretch()
        
        # Save all results button
        self.save_all_btn = QPushButton("💾 Save All Results")
        self.save_all_btn.setToolTip("Save all analysis results (InterPro + PROSITE) to file")
        self.save_all_btn.clicked.connect(self.save_all_results)
        self.save_all_btn.setStyleSheet("background-color: #FF9800; color: white; padding: 8px 20px; border-radius: 3px; font-weight: bold;")
        self.save_all_btn.setEnabled(False)  # Enable when results are available
        save_all_layout.addWidget(self.save_all_btn)
        
        layout.addWidget(save_all_widget)
        
        # InterPro results tab
        self.interpro_tab = QScrollArea()
        self.interpro_widget = QWidget()
        self.interpro_layout = QVBoxLayout(self.interpro_widget)
        self.interpro_tab.setWidget(self.interpro_widget)
        self.interpro_tab.setWidgetResizable(True)
        self.results_tabs.addTab(self.interpro_tab, "InterPro")
        
        # PROSITE results tab
        self.prosite_tab = QScrollArea()
        self.prosite_widget = QWidget()
        self.prosite_layout = QVBoxLayout(self.prosite_widget)
        self.prosite_tab.setWidget(self.prosite_widget)
        self.prosite_tab.setWidgetResizable(True)
        self.results_tabs.addTab(self.prosite_tab, "PROSITE")
        
        layout.addWidget(self.results_tabs, 1)
        
        # Initially show placeholder
        self.show_placeholder()
    
    def show_placeholder(self):
        """Show placeholder content when no analysis has been performed."""
        placeholder_text = """
        <h3>🎯 Protein Motifs and Domains Analysis</h3>
        <p>Enter a protein sequence above and click one of the search buttons to analyze:</p>
        <ul>
        <li><b>🔬 InterPro:</b> Comprehensive protein domain and functional site annotation using the EBI InterPro API<br>
        <i>Includes Pfam, SMART, PROSITE, CDD, PRINTS, and many more databases</i></li>
        <li><b>🧬 PROSITE:</b> Protein motif and pattern recognition using the official ExPASy PROSITE ScanProsite API<br>
        <i>With intelligent local pattern fallback for reliability</i></li>
        </ul>
        <div style="background-color: #e8f5e8; border: 1px solid #27ae60; border-radius: 6px; padding: 15px; margin: 15px 0;">
        <p style="margin: 0; color: #27ae60;"><b>💡 Pro Tips:</b></p>
        <ul style="margin: 5px 0; color: #27ae60;">
        <li>Use <b>"Search All Databases"</b> for the most comprehensive analysis</li>
        <li>Toggle <b>"Exclude high-probability motifs"</b> to filter common patterns</li>
        <li>Results include interactive visualization and detailed functional annotations</li>
        <li>Export results in multiple formats including HTML reports</li>
        </ul>
        </div>
        <p><i>🌐 This tool connects to live web services for the most up-to-date annotations.</i></p>
        """
        
        for tab_layout in [self.interpro_layout, self.prosite_layout]:
            # Clear existing widgets
            for i in reversed(range(tab_layout.count())):
                item = tab_layout.itemAt(i)
                if item and item.widget():
                    item.widget().setParent(None)
            
            # Add placeholder
            placeholder = QLabel(placeholder_text)
            placeholder.setAlignment(Qt.AlignTop)
            placeholder.setWordWrap(True)
            placeholder.setStyleSheet("color: #666; padding: 20px;")
            tab_layout.addWidget(placeholder)
            tab_layout.addStretch()
    
    def load_from_structure(self):
        """Load sequence from the currently displayed structure."""
        if hasattr(self.parent_app, 'sequence_display'):
            sequence_text = self.parent_app.sequence_display.toPlainText()
            if sequence_text and not sequence_text.startswith("Sequence will appear"):
                # Extract just the sequence part (remove FASTA headers)
                lines = sequence_text.split('\n')
                sequence_lines = [line for line in lines if not line.startswith('>')]
                sequence = ''.join(sequence_lines).replace(' ', '').upper()
                
                if sequence:
                    self.sequence_input.setPlainText(sequence)
                    self.current_sequence = sequence
                    if hasattr(self.parent_app, 'statusBar'):
                        self.parent_app.statusBar().showMessage("Sequence loaded from current structure")
                else:
                    QMessageBox.warning(self, "No Sequence", "No valid sequence found in the current structure.")
            else:
                QMessageBox.warning(self, "No Structure", "No structure is currently loaded.")
        else:
            QMessageBox.warning(self, "Error", "Cannot access sequence display.")
    
    def load_from_file(self):
        """Load sequence from a FASTA file with comprehensive support."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Sequence File", "", 
            "FASTA Files (*.fasta *.fa *.fas *.fna *.faa);;Text Files (*.txt);;All Files (*)"
        )
        
        if file_path:
            try:
                # Parse FASTA file
                sequences = self.parse_fasta_file(file_path)
                
                if not sequences:
                    QMessageBox.warning(self, "No Sequences", "No valid sequences found in the file.")
                    return
                
                # Handle single vs multiple sequences
                if len(sequences) == 1:
                    # Single sequence - load directly
                    seq_id, sequence = sequences[0]
                    self.load_sequence_to_input(sequence, seq_id, file_path)
                else:
                    # Multiple sequences - show selection dialog
                    selected_sequence = self.show_sequence_selection_dialog(sequences, file_path)
                    if selected_sequence:
                        seq_id, sequence = selected_sequence
                        self.load_sequence_to_input(sequence, seq_id, file_path)
                    
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load sequence file:\n{str(e)}")
    
    def parse_fasta_file(self, file_path):
        """Parse FASTA file and return list of (id, sequence) tuples."""
        sequences = []
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Handle different line endings
            content = content.replace('\r\n', '\n').replace('\r', '\n')
            
            # Split into potential FASTA entries
            entries = content.split('>')
            
            for entry in entries:
                entry = entry.strip()
                if not entry:
                    continue
                
                lines = entry.split('\n')
                if not lines:
                    continue
                
                # First line is the header
                header = lines[0].strip()
                
                # Remaining lines are sequence
                sequence_lines = [line.strip() for line in lines[1:] if line.strip()]
                sequence = ''.join(sequence_lines).upper()
                
                # Validate sequence
                if sequence and self.is_valid_protein_sequence(sequence):
                    # Extract meaningful ID from header
                    seq_id = self.extract_sequence_id(header)
                    sequences.append((seq_id, sequence))
                elif sequence:
                    # Check if it might be a nucleotide sequence
                    if self.is_nucleotide_sequence(sequence):
                        QMessageBox.information(
                            self, "Nucleotide Sequence Detected", 
                            f"The sequence '{header[:50]}...' appears to be a nucleotide sequence.\n"
                            "This tool is designed for protein sequences. Please use a protein sequence or "
                            "translate your nucleotide sequence first."
                        )
                    else:
                        print(f"[DEBUG] Invalid sequence found: {header[:50]}...")
            
            return sequences
            
        except Exception as e:
            print(f"[ERROR] Error parsing FASTA file: {e}")
            raise
    
    def extract_sequence_from_text(self, text):
        """Extract protein sequence from text, handling FASTA format."""
        lines = text.split('\n')
        sequence_lines = []
        
        for line in lines:
            line = line.strip()
            # Skip FASTA headers and empty lines
            if line.startswith('>') or not line:
                continue
            # Skip comment lines
            if line.startswith('#') or line.startswith(';'):
                continue
            sequence_lines.append(line)
        
        # Join all sequence lines and remove whitespace
        sequence = ''.join(sequence_lines)
        sequence = ''.join(sequence.split())  # Remove all whitespace
        
        return sequence.upper()
    
    def is_valid_protein_sequence(self, sequence):
        """Check if sequence contains valid protein amino acid codes."""
        if not sequence:
            return False
        
        # Standard amino acid codes + ambiguous codes
        valid_chars = set('ACDEFGHIKLMNPQRSTVWYXBZJUO*-')
        sequence_chars = set(sequence.upper())
        
        # Allow sequence if at least 80% of characters are valid amino acids
        valid_count = len(sequence_chars.intersection(valid_chars))
        total_count = len(sequence_chars)
        
        return (valid_count / total_count) >= 0.8 if total_count > 0 else False
    
    def is_nucleotide_sequence(self, sequence):
        """Check if sequence appears to be nucleotide sequence."""
        if not sequence:
            return False
        
        nucleotide_chars = set('ATCGRYSWKMBDHVN-')
        sequence_chars = set(sequence.upper())
        
        # Consider it nucleotide if >90% are nucleotide characters
        nucleotide_count = len(sequence_chars.intersection(nucleotide_chars))
        total_count = len(sequence_chars)
        
        return (nucleotide_count / total_count) >= 0.9 if total_count > 0 else False
    
    def extract_sequence_id(self, header):
        """Extract a meaningful sequence ID from FASTA header."""
        if not header:
            return "Unknown"
        
        # Try to extract common ID patterns
        import re
        
        # UniProt pattern: >sp|P12345|PROTEIN_HUMAN or >tr|A0A123|PROTEIN_HUMAN
        uniprot_match = re.search(r'>?(?:sp|tr)\|([^|]+)\|([^\s]+)', header)
        if uniprot_match:
            return f"{uniprot_match.group(1)} ({uniprot_match.group(2)})"
        
        # NCBI pattern: >gi|123456|ref|NP_123456.1| or >NP_123456.1
        ncbi_match = re.search(r'>?(?:gi\|\d+\|)?(?:ref\|)?([A-Z]{2}_\d+(?:\.\d+)?)', header)
        if ncbi_match:
            return ncbi_match.group(1)
        
        # PDB pattern: >1ABC_A or >pdb|1ABC|A
        pdb_match = re.search(r'>?(?:pdb\|)?([0-9][A-Z0-9]{3})(?:\|([A-Z]))?', header)
        if pdb_match:
            chain = f"_{pdb_match.group(2)}" if pdb_match.group(2) else ""
            return f"{pdb_match.group(1)}{chain}"
        
        # Generic pattern: take first word after >
        generic_match = re.search(r'>?([^\s|]+)', header)
        if generic_match:
            return generic_match.group(1)
        
        # Fallback: use first 20 characters
        return header[:20].replace('>', '')
    
    def show_sequence_selection_dialog(self, sequences, file_path):
        """Show dialog for selecting sequence from multiple sequences."""
        dialog = QDialog(self)
        dialog.setWindowTitle(f"Select Sequence - {os.path.basename(file_path)}")
        dialog.setMinimumSize(600, 400)
        
        layout = QVBoxLayout(dialog)
        
        # Info label
        info_label = QLabel(f"Found {len(sequences)} sequences in the file. Please select one:")
        layout.addWidget(info_label)
        
        # Sequence list
        sequence_list = QListWidget()
        
        for i, (seq_id, sequence) in enumerate(sequences):
            # Create descriptive item text
            item_text = f"{seq_id} ({len(sequence)} amino acids)"
            if len(sequence) > 50:
                preview = sequence[:50] + "..."
            else:
                preview = sequence
            
            item = QListWidgetItem(f"{item_text}\n{preview}")
            item.setData(Qt.UserRole, (seq_id, sequence))  # Store sequence data
            sequence_list.addItem(item)
        
        sequence_list.setCurrentRow(0)  # Select first item by default
        layout.addWidget(sequence_list)
        
        # Buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)
        
        # Show dialog and get result
        if dialog.exec_() == QDialog.Accepted:
            current_item = sequence_list.currentItem()
            if current_item:
                return current_item.data(Qt.UserRole)
        
        return None
    
    def load_sequence_to_input(self, sequence, seq_id, file_path):
        """Load sequence into the input field with proper formatting."""
        # Create FASTA format for display
        fasta_text = f">{seq_id}\n{sequence}"
        
        self.sequence_input.setPlainText(fasta_text)
        self.current_sequence = sequence
        
        if hasattr(self.parent_app, 'statusBar'):
            self.parent_app.statusBar().showMessage(
                f"Loaded sequence '{seq_id}' ({len(sequence)} amino acids) from {os.path.basename(file_path)}"
            )
    
    def validate_sequence(self):
        """Validate the input sequence with enhanced FASTA support."""
        raw_text = self.sequence_input.toPlainText().strip()
        
        if not raw_text:
            QMessageBox.warning(self, "No Sequence", "Please enter a protein sequence.")
            return None
        
        # Parse sequence (handle FASTA format)
        sequence = self.extract_sequence_from_text(raw_text)
        
        if not sequence:
            QMessageBox.warning(self, "No Sequence", "No valid sequence found in the input.")
            return None
        
        # Check for valid protein sequence
        valid_chars = set('ACDEFGHIKLMNPQRSTVWYXBZJUO*-')
        sequence_chars = set(sequence.upper())
        invalid_chars = sequence_chars - valid_chars
        
        if invalid_chars:
            # Check if it might be a nucleotide sequence
            if self.is_nucleotide_sequence(sequence):
                QMessageBox.warning(
                    self, "Nucleotide Sequence Detected", 
                    "This appears to be a nucleotide sequence. This tool is designed for protein sequences.\n"
                    "Please use a protein sequence or translate your nucleotide sequence first."
                )
            else:
                QMessageBox.warning(
                    self, "Invalid Sequence", 
                    f"Invalid characters for protein sequence: {', '.join(sorted(invalid_chars))}\n\n"
                    "Valid amino acid codes: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y\n"
                    "Ambiguous codes: X, B, Z, J, U, O\n"
                    "Special characters: *, -"
                )
            return None
        
        # Check sequence length
        if len(sequence) < 10:
            QMessageBox.warning(
                self, "Sequence Too Short", 
                f"Sequence is only {len(sequence)} amino acids long. "
                "Please enter a longer sequence for meaningful analysis."
            )
            return None
        
        if len(sequence) > 10000:
            reply = QMessageBox.question(
                self, "Long Sequence", 
                f"This sequence is {len(sequence)} amino acids long. "
                "Analysis may take a long time. Continue?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply != QMessageBox.Yes:
                return None
        
        if len(sequence) < 10:
            QMessageBox.warning(self, "Sequence Too Short", "Sequence must be at least 10 amino acids long.")
            return None
        
        self.current_sequence = sequence
        return sequence
    
    def search_interpro(self):
        """Search InterPro database."""
        sequence = self.validate_sequence()
        if not sequence:
            return
        
        self.set_search_state(True, "Searching InterPro...")
        
        self.interpro_worker = InterProSearchWorker(sequence)
        self.interpro_worker.search_complete.connect(self.handle_interpro_results)
        self.interpro_worker.error_occurred.connect(self.handle_search_error)
        self.interpro_worker.progress_update.connect(self.update_status)
        self.interpro_worker.start()
    

    
    def search_prosite(self):
        """Search PROSITE patterns."""
        sequence = self.validate_sequence()
        if not sequence:
            return
        
        self.set_search_state(True, "Searching PROSITE...")
        
        # Get filtering option
        exclude_high_prob = self.exclude_high_prob_checkbox.isChecked()
        
        self.prosite_worker = PROSITESearchWorker(sequence, exclude_high_prob)
        self.prosite_worker.search_complete.connect(self.handle_prosite_results)
        self.prosite_worker.error_occurred.connect(self.handle_search_error)
        self.prosite_worker.progress_update.connect(self.update_status)
        self.prosite_worker.start()
    
    def search_all(self):
        """Search all databases."""
        sequence = self.validate_sequence()
        if not sequence:
            return
        
        # Start with PROSITE (fastest)
        self.search_prosite()
        
        # Start InterPro search after a short delay
        QTimer.singleShot(1000, self.search_interpro)
    
    def set_search_state(self, searching, message=""):
        """Set the UI state during search."""
        self.interpro_btn.setEnabled(not searching)
        self.prosite_btn.setEnabled(not searching)
        self.search_all_btn.setEnabled(not searching)
        
        self.progress_bar.setVisible(searching)
        if searching:
            self.progress_bar.setRange(0, 0)  # Indeterminate
        
        self.status_label.setText(message)
    
    def update_status(self, message):
        """Update status message."""
        self.status_label.setText(message)
    
    def handle_interpro_results(self, results):
        """Handle InterPro search results."""
        self.set_search_state(False)
        self.display_interpro_results(results)
        self.update_save_button_state()
    

    
    def handle_prosite_results(self, results):
        """Handle PROSITE search results."""
        print(f"[DEBUG] Handling PROSITE results: {type(results)}")
        if results:
            print(f"[DEBUG] Results keys: {list(results.keys()) if isinstance(results, dict) else 'Not a dict'}")
            if isinstance(results, dict) and 'motifs' in results:
                print(f"[DEBUG] Found {len(results['motifs'])} motifs")
                for i, motif in enumerate(results['motifs'][:3]):  # Show first 3
                    print(f"[DEBUG] Motif {i+1}: {motif.get('id', 'No ID')} - {motif.get('name', 'No name')}")
        
        self.set_search_state(False)
        self.display_prosite_results(results)
        self.update_save_button_state()
    
    def handle_search_error(self, error_message):
        """Handle search errors."""
        self.set_search_state(False)
        QMessageBox.critical(self, "Search Error", error_message)
    
    def display_interpro_results(self, results):
        """Display InterPro results in website-like format."""
        # Store results for visualization and saving
        self.current_interpro_results = results
        
        # Clear existing results
        for i in reversed(range(self.interpro_layout.count())):
            item = self.interpro_layout.itemAt(i)
            if item and item.widget():
                item.widget().setParent(None)
        
        if not results or not results.get('domains'):
            no_results = QLabel("No InterPro domains found. This may be due to API connectivity issues or the sequence not containing recognizable domains.")
            no_results.setStyleSheet("color: #666; padding: 20px;")
            self.interpro_layout.addWidget(no_results)
            self.interpro_layout.addStretch()
            return
        
        # Header with job info and controls
        self.create_interpro_header(results)
        
        # Sequence overview section
        self.create_sequence_overview(results)
        
        # InterPro entries section (grouped by InterPro entry)
        self.create_interpro_entries_section(results)
        
        # Member database matches section
        self.create_member_database_section(results)
        
        # GO terms section
        if results.get('go_terms'):
            self.create_go_terms_section(results)
        
        # Pathways section
        if results.get('pathways'):
            self.create_pathways_section(results)
        
        self.interpro_layout.addStretch()
    
    def create_interpro_header(self, results):
        """Create the header section with job info and controls."""
        header_layout = QVBoxLayout()
        
        # Job information section
        if results.get('job_id'):
            job_info_layout = QHBoxLayout()
            
            job_label = QLabel(f"🔍 <b>InterPro Job ID:</b> {results['job_id']}")
            job_label.setStyleSheet("color: #2196F3; padding: 5px; background-color: #e3f2fd; border-radius: 3px;")
            job_info_layout.addWidget(job_label)
            
            # Button to open InterPro website
            view_online_btn = QPushButton("🌐 View on InterPro Website")
            view_online_btn.setToolTip(f"Open InterPro job {results['job_id']} in web browser")
            view_online_btn.clicked.connect(lambda: self.open_interpro_job(results['job_id']))
            view_online_btn.setStyleSheet("background-color: #2196F3; color: white; padding: 5px 15px; border-radius: 3px;")
            job_info_layout.addWidget(view_online_btn)
            
            job_info_layout.addStretch()
            
            job_widget = QWidget()
            job_widget.setLayout(job_info_layout)
            header_layout.addWidget(job_widget)
        
        # Summary and controls
        summary_layout = QHBoxLayout()
        
        # Count different types of matches
        interpro_entries = set()
        member_db_matches = 0
        
        for domain in results['domains']:
            if domain.get('interpro'):
                interpro_entries.add(domain['interpro'].get('accession', ''))
            member_db_matches += 1
        
        summary_text = f"Found {len(interpro_entries)} InterPro entries and {member_db_matches} member database matches"
        summary = QLabel(summary_text)
        summary.setStyleSheet("font-weight: bold; font-size: 14px; padding: 10px;")
        summary_layout.addWidget(summary)
        
        summary_layout.addStretch()
        
        # Save button
        save_interpro_btn = QPushButton("💾 Save Results")
        save_interpro_btn.setToolTip("Save InterPro results to file")
        save_interpro_btn.clicked.connect(self.save_interpro_results)
        save_interpro_btn.setStyleSheet("background-color: #4CAF50; color: white; padding: 5px 15px; border-radius: 3px;")
        summary_layout.addWidget(save_interpro_btn)
        
        summary_widget = QWidget()
        summary_widget.setLayout(summary_layout)
        header_layout.addWidget(summary_widget)
        
        header_widget = QWidget()
        header_widget.setLayout(header_layout)
        self.interpro_layout.addWidget(header_widget)
    
    def create_sequence_overview(self, results):
        """Create sequence overview with basic information."""
        overview_group = QGroupBox("Sequence Overview")
        overview_layout = QVBoxLayout(overview_group)
        
        # Sequence length info
        seq_length = len(self.current_sequence) if hasattr(self, 'current_sequence') else 'Unknown'
        length_label = QLabel(f"Sequence length: {seq_length} amino acids")
        length_label.setStyleSheet("font-weight: bold; padding: 5px;")
        overview_layout.addWidget(length_label)
        
        # Domain count summary
        domain_count = len(results['domains'])
        domain_summary = QLabel(f"Total domains/sites found: {domain_count}")
        domain_summary.setStyleSheet("padding: 5px;")
        overview_layout.addWidget(domain_summary)
        
        self.interpro_layout.addWidget(overview_group)
    
    def create_interpro_entries_section(self, results):
        """Create InterPro entries section grouped by InterPro accession."""
        # Group domains by InterPro entry
        interpro_groups = {}
        for domain in results['domains']:
            interpro_entry = domain.get('interpro')
            if interpro_entry and interpro_entry.get('accession'):
                acc = interpro_entry['accession']
                if acc not in interpro_groups:
                    interpro_groups[acc] = {
                        'entry': interpro_entry,
                        'matches': []
                    }
                interpro_groups[acc]['matches'].append(domain)
        
        if interpro_groups:
            entries_group = QGroupBox(f"InterPro Entries ({len(interpro_groups)})")
            entries_layout = QVBoxLayout(entries_group)
            
            for acc, group_data in sorted(interpro_groups.items()):
                entry_widget = self.create_interpro_entry_widget(acc, group_data)
                entries_layout.addWidget(entry_widget)
            
            self.interpro_layout.addWidget(entries_group)
    
    def create_interpro_entry_widget(self, accession, group_data):
        """Create widget for a single InterPro entry."""
        entry = group_data['entry']
        matches = group_data['matches']
        
        entry_frame = QFrame()
        entry_frame.setFrameStyle(QFrame.Box)
        entry_frame.setStyleSheet("QFrame { border: 1px solid #ddd; border-radius: 5px; margin: 2px; }")
        
        entry_layout = QVBoxLayout(entry_frame)
        
        # Entry header
        header_layout = QHBoxLayout()
        
        # InterPro accession and name
        title_text = f"<b>{accession}</b>"
        if entry.get('name'):
            title_text += f" - {entry['name']}"
        
        title_label = QLabel(title_text)
        title_label.setStyleSheet("font-size: 14px; color: #2196F3; padding: 5px;")
        header_layout.addWidget(title_label)
        
        header_layout.addStretch()
        
        # Entry type
        if entry.get('type'):
            type_label = QLabel(f"Type: {entry['type']}")
            type_label.setStyleSheet("background-color: #f0f0f0; padding: 3px 8px; border-radius: 3px; font-size: 11px;")
            header_layout.addWidget(type_label)
        
        entry_layout.addLayout(header_layout)
        
        # Description
        if entry.get('description'):
            desc_label = QLabel(entry['description'])
            desc_label.setWordWrap(True)
            desc_label.setStyleSheet("color: #666; padding: 5px; font-style: italic;")
            entry_layout.addWidget(desc_label)
        
        # Member database matches
        matches_label = QLabel(f"Member database matches ({len(matches)}):")
        matches_label.setStyleSheet("font-weight: bold; padding: 5px 0px;")
        entry_layout.addWidget(matches_label)
        
        # Matches table
        matches_table = QTableWidget()
        matches_table.setColumnCount(6)
        matches_table.setHorizontalHeaderLabels([
            "Database", "Accession", "Name", "Start", "End", "E-value"
        ])
        matches_table.setRowCount(len(matches))
        
        for i, match in enumerate(matches):
            matches_table.setItem(i, 0, QTableWidgetItem(match.get('database', '')))
            matches_table.setItem(i, 1, QTableWidgetItem(match.get('accession', '')))
            matches_table.setItem(i, 2, QTableWidgetItem(match.get('name', '')))
            
            # Location info
            locations = match.get('locations', [])
            if locations:
                matches_table.setItem(i, 3, QTableWidgetItem(str(locations[0].get('start', ''))))
                matches_table.setItem(i, 4, QTableWidgetItem(str(locations[0].get('end', ''))))
            
            matches_table.setItem(i, 5, QTableWidgetItem(str(match.get('evalue', 'N/A'))))
        
        matches_table.resizeColumnsToContents()
        matches_table.setMaximumHeight(150)
        entry_layout.addWidget(matches_table)
        
        return entry_frame
    
    def create_member_database_section(self, results):
        """Create member database matches section."""
        # Group by database
        db_groups = {}
        for domain in results['domains']:
            db = domain.get('database', 'Unknown')
            if db not in db_groups:
                db_groups[db] = []
            db_groups[db].append(domain)
        
        if db_groups:
            db_group = QGroupBox(f"Member Database Matches ({len(results['domains'])} total)")
            db_layout = QVBoxLayout(db_group)
            
            for db_name, domains in sorted(db_groups.items()):
                db_widget = self.create_database_section_widget(db_name, domains)
                db_layout.addWidget(db_widget)
            
            self.interpro_layout.addWidget(db_group)
    
    def create_database_section_widget(self, db_name, domains):
        """Create widget for matches from a specific database."""
        db_frame = QFrame()
        db_frame.setFrameStyle(QFrame.StyledPanel)
        db_frame.setStyleSheet("QFrame { border: 1px solid #ccc; border-radius: 3px; margin: 1px; }")
        
        db_layout = QVBoxLayout(db_frame)
        
        # Database header
        header_layout = QHBoxLayout()
        
        db_label = QLabel(f"<b>{db_name}</b> ({len(domains)} matches)")
        db_label.setStyleSheet("font-size: 13px; padding: 5px;")
        header_layout.addWidget(db_label)
        
        header_layout.addStretch()
        
        # Database-specific styling
        if db_name.lower() == 'pfam':
            db_label.setStyleSheet("font-size: 13px; padding: 5px; color: #FF6B35;")
        elif db_name.lower() == 'smart':
            db_label.setStyleSheet("font-size: 13px; padding: 5px; color: #4CAF50;")
        elif db_name.lower() == 'prosite':
            db_label.setStyleSheet("font-size: 13px; padding: 5px; color: #2196F3;")
        
        db_layout.addLayout(header_layout)
        
        # Matches table
        table = QTableWidget()
        table.setColumnCount(6)
        table.setHorizontalHeaderLabels([
            "Accession", "Name", "Description", "Start", "End", "E-value"
        ])
        table.setRowCount(len(domains))
        
        for i, domain in enumerate(domains):
            table.setItem(i, 0, QTableWidgetItem(domain.get('accession', '')))
            table.setItem(i, 1, QTableWidgetItem(domain.get('name', '')))
            table.setItem(i, 2, QTableWidgetItem(domain.get('description', '')))
            
            # Location info
            locations = domain.get('locations', [])
            if locations:
                table.setItem(i, 3, QTableWidgetItem(str(locations[0].get('start', ''))))
                table.setItem(i, 4, QTableWidgetItem(str(locations[0].get('end', ''))))
            
            table.setItem(i, 5, QTableWidgetItem(str(domain.get('evalue', 'N/A'))))
            
            # Highlight significant matches
            evalue = domain.get('evalue', 'N/A')
            if evalue != 'N/A':
                try:
                    eval_float = float(evalue)
                    if eval_float < 1e-10:
                        # Very significant - green background
                        for col in range(6):
                            item = table.item(i, col)
                            if item:
                                item.setBackground(QColor(200, 255, 200))
                    elif eval_float < 1e-5:
                        # Significant - light green background
                        for col in range(6):
                            item = table.item(i, col)
                            if item:
                                item.setBackground(QColor(230, 255, 230))
                except ValueError:
                    pass
        
        table.resizeColumnsToContents()
        table.setMaximumHeight(min(200, 30 + len(domains) * 25))
        db_layout.addWidget(table)
        
        return db_frame
    
    def create_go_terms_section(self, results):
        """Create GO terms section."""
        go_group = QGroupBox(f"Gene Ontology Terms ({len(results['go_terms'])})")
        go_layout = QVBoxLayout(go_group)
        
        # Group GO terms by category
        go_categories = {}
        for go_term in results['go_terms']:
            category = go_term.get('category', 'unknown')
            if category not in go_categories:
                go_categories[category] = []
            go_categories[category].append(go_term)
        
        for category, terms in sorted(go_categories.items()):
            category_label = QLabel(f"<b>{category.replace('_', ' ').title()} ({len(terms)})</b>")
            category_label.setStyleSheet("padding: 5px; background-color: #f5f5f5; border-radius: 3px; margin: 2px;")
            go_layout.addWidget(category_label)
            
            for term in terms:
                term_layout = QHBoxLayout()
                
                go_id = QLabel(term.get('id', ''))
                go_id.setStyleSheet("font-family: monospace; color: #2196F3; font-weight: bold;")
                go_id.setMinimumWidth(100)
                term_layout.addWidget(go_id)
                
                go_name = QLabel(term.get('name', ''))
                go_name.setWordWrap(True)
                term_layout.addWidget(go_name)
                
                term_layout.addStretch()
                
                term_widget = QWidget()
                term_widget.setLayout(term_layout)
                go_layout.addWidget(term_widget)
        
        self.interpro_layout.addWidget(go_group)
    
    def create_pathways_section(self, results):
        """Create pathways section."""
        pathways_group = QGroupBox(f"Pathways ({len(results['pathways'])})")
        pathways_layout = QVBoxLayout(pathways_group)
        
        for pathway in results['pathways']:
            pathway_layout = QHBoxLayout()
            
            pathway_id = QLabel(pathway.get('id', ''))
            pathway_id.setStyleSheet("font-family: monospace; color: #FF6B35; font-weight: bold;")
            pathway_id.setMinimumWidth(120)
            pathway_layout.addWidget(pathway_id)
            
            pathway_name = QLabel(pathway.get('name', ''))
            pathway_name.setWordWrap(True)
            pathway_layout.addWidget(pathway_name)
            
            if pathway.get('database'):
                db_label = QLabel(f"[{pathway['database']}]")
                db_label.setStyleSheet("color: #666; font-size: 11px;")
                pathway_layout.addWidget(db_label)
            
            pathway_layout.addStretch()
            
            pathway_widget = QWidget()
            pathway_widget.setLayout(pathway_layout)
            pathways_layout.addWidget(pathway_widget)
        
        self.interpro_layout.addWidget(pathways_group)
    

    
    def display_prosite_results(self, results):
        """Display PROSITE search results."""
        # Store results for saving
        self.current_prosite_results = results
        
        # Clear existing results
        for i in reversed(range(self.prosite_layout.count())):
            item = self.prosite_layout.itemAt(i)
            if item and item.widget():
                item.widget().setParent(None)
        
        if not results or not results.get('motifs'):
            no_results = QLabel("No PROSITE motifs found in this sequence.")
            no_results.setStyleSheet("color: #666; padding: 20px;")
            self.prosite_layout.addWidget(no_results)
            self.prosite_layout.addStretch()
            return
        
        # Header with summary and save button
        header_layout = QHBoxLayout()
        
        # Summary with motif count and categories
        motif_count = results.get('pattern_count', 0)
        category_count = len(results.get('categories', {}))
        
        summary_text = f"Found {motif_count} motif(s) across {category_count} functional categories"
        summary = QLabel(summary_text)
        summary.setStyleSheet("font-weight: bold; font-size: 14px; padding: 10px;")
        header_layout.addWidget(summary)
        
        header_layout.addStretch()
        
        # Save PROSITE results button
        save_prosite_btn = QPushButton("💾 Save PROSITE Results")
        save_prosite_btn.setToolTip("Save PROSITE results to file")
        save_prosite_btn.clicked.connect(self.save_prosite_results)
        save_prosite_btn.setStyleSheet("background-color: #9C27B0; color: white; padding: 5px 15px; border-radius: 3px;")
        header_layout.addWidget(save_prosite_btn)
        
        header_widget = QWidget()
        header_widget.setLayout(header_layout)
        self.prosite_layout.addWidget(header_widget)
        
        # Category statistics
        if results.get('categories'):
            stats_frame = QFrame()
            stats_frame.setStyleSheet("background-color: #f5f5f5; border: 1px solid #ddd; border-radius: 5px; padding: 10px;")
            stats_layout = QVBoxLayout(stats_frame)
            
            stats_title = QLabel("📊 <b>Motif Categories Found:</b>")
            stats_layout.addWidget(stats_title)
            
            categories_layout = QHBoxLayout()
            for category, count in results['categories'].items():
                category_label = QLabel(f"<b>{category}:</b> {count}")
                category_label.setStyleSheet("padding: 3px 8px; margin: 2px; background-color: white; border-radius: 3px;")
                categories_layout.addWidget(category_label)
            
            categories_layout.addStretch()
            stats_layout.addLayout(categories_layout)
            self.prosite_layout.addWidget(stats_frame)
        
        # Motifs table
        motifs_table = QTableWidget()
        motifs_table.setColumnCount(9)
        motifs_table.setHorizontalHeaderLabels([
            "Category", "ID", "Name", "Function", "Start", "End", "Length", "Sequence", "Source"
        ])
        
        # Set table properties
        motifs_table.horizontalHeader().setStretchLastSection(True)
        motifs_table.setAlternatingRowColors(True)
        motifs_table.setSelectionBehavior(QTableWidget.SelectRows)
        motifs_table.setSortingEnabled(True)
        
        # Populate table
        motifs = results.get('motifs', [])
        motifs_table.setRowCount(len(motifs))
        
        # Color coding for categories
        category_colors = {
            'Post-translational modification': '#FF9800',
            'Nucleotide binding': '#4CAF50',
            'DNA/RNA binding': '#2196F3',
            'Protein-protein interaction': '#9C27B0',
            'Enzymatic site': '#F44336',
            'Membrane protein': '#00BCD4',
            'Signal sequence': '#795548',
            'Structural motif': '#607D8B'
        }
        
        for row, motif in enumerate(motifs):
            # Category with color coding
            category = motif.get('category', '')
            category_item = QTableWidgetItem(category)
            if category in category_colors:
                category_item.setBackground(QColor(category_colors[category]))
                category_item.setForeground(QColor('white'))
            motifs_table.setItem(row, 0, category_item)
            
            # Other columns
            motifs_table.setItem(row, 1, QTableWidgetItem(motif.get('id', '')))
            motifs_table.setItem(row, 2, QTableWidgetItem(motif.get('name', '')))
            motifs_table.setItem(row, 3, QTableWidgetItem(motif.get('function', '')))
            motifs_table.setItem(row, 4, QTableWidgetItem(str(motif.get('start', ''))))
            motifs_table.setItem(row, 5, QTableWidgetItem(str(motif.get('end', ''))))
            motifs_table.setItem(row, 6, QTableWidgetItem(str(motif.get('length', ''))))
            
            # Sequence with monospace font
            seq_item = QTableWidgetItem(motif.get('sequence', ''))
            seq_item.setFont(QFont("Courier", 10))
            motifs_table.setItem(row, 7, seq_item)
            
            # Add source column
            source = motif.get('source', 'PROSITE API')
            source_item = QTableWidgetItem(source)
            source_item.setBackground(QColor(230, 240, 255))  # Light blue for API source
            source_item.setToolTip("Results from official PROSITE ScanProsite API")
            motifs_table.setItem(row, 8, source_item)
        
        # Resize columns to content
        motifs_table.resizeColumnsToContents()
        
        self.prosite_layout.addWidget(motifs_table)
        
        # Context viewer for selected motifs
        context_frame = QFrame()
        context_frame.setStyleSheet("background-color: #f9f9f9; border: 1px solid #ddd; border-radius: 5px; padding: 10px;")
        context_layout = QVBoxLayout(context_frame)
        
        context_title = QLabel("🔍 <b>Sequence Context:</b> (Click a motif above to view)")
        context_layout.addWidget(context_title)
        
        self.context_display = QTextEdit()
        self.context_display.setMaximumHeight(100)
        self.context_display.setFont(QFont("Courier", 11))
        self.context_display.setStyleSheet("background-color: white; padding: 10px; border-radius: 3px;")
        self.context_display.setPlainText("Select a motif from the table above to view its sequence context.")
        self.context_display.setReadOnly(True)
        context_layout.addWidget(self.context_display)
        
        self.prosite_layout.addWidget(context_frame)
        
        # Connect table selection to context display
        motifs_table.itemSelectionChanged.connect(lambda: self.show_motif_context(motifs_table, motifs))
        
        self.prosite_layout.addStretch()
    
    def show_motif_context(self, table, motifs):
        """Show sequence context for selected motif."""
        current_row = table.currentRow()
        if current_row >= 0 and current_row < len(motifs):
            motif = motifs[current_row]
            context = motif.get('context', {})
            
            if context and 'sequence' in context:
                context_seq = context['sequence']
                motif_start = context.get('motif_start', 0)
                motif_end = context.get('motif_end', 0)
                
                # Format the context with highlighting
                before = context_seq[:motif_start]
                motif_seq = context_seq[motif_start:motif_end]
                after = context_seq[motif_end:]
                
                formatted_context = f"Motif: {motif.get('name', '')}\nPattern: {motif.get('pattern', '')}\nPosition: {motif.get('start', '')}-{motif.get('end', '')}\nContext: {before}[{motif_seq}]{after}"
                
                self.context_display.setPlainText(formatted_context)
            else:
                self.context_display.setPlainText(f"Context not available for {motif.get('name', 'this motif')}.")
    
    def display_prosite_results(self, results):
        """Display enhanced PROSITE results."""
        # Store results for saving
        self.current_prosite_results = results
        
        # Clear existing results
        for i in reversed(range(self.prosite_layout.count())):
            item = self.prosite_layout.itemAt(i)
            if item and item.widget():
                item.widget().setParent(None)
        
        if not results or not results.get('motifs'):
            no_results = QLabel("No PROSITE motifs found.")
            no_results.setStyleSheet("color: #666; padding: 20px;")
            self.prosite_layout.addWidget(no_results)
            self.prosite_layout.addStretch()
            return
        
        # Header with summary and save button
        header_layout = QHBoxLayout()
        
        # Summary with statistics
        summary_text = f"Found {len(results['motifs'])} motif(s) in {len(results.get('categories', {}))} categories"
        summary = QLabel(summary_text)
        summary.setStyleSheet("font-weight: bold; font-size: 14px; padding: 10px;")
        header_layout.addWidget(summary)
        
        header_layout.addStretch()
        
        # Save button
        save_btn = QPushButton("💾 Save Results")
        save_btn.setToolTip("Save PROSITE results to file")
        save_btn.clicked.connect(self.save_prosite_results)
        save_btn.setStyleSheet("background-color: #2196F3; color: white; padding: 5px 15px; border-radius: 3px;")
        header_layout.addWidget(save_btn)
        
        header_widget = QWidget()
        header_widget.setLayout(header_layout)
        self.prosite_layout.addWidget(header_widget)
        
        # Category statistics
        if results.get('categories'):
            stats_group = QGroupBox("Motif Categories")
            stats_layout = QHBoxLayout(stats_group)
            
            for category, count in results['categories'].items():
                cat_label = QLabel(f"<b>{category}:</b> {count}")
                cat_label.setStyleSheet("padding: 5px; margin: 2px; background-color: #f0f0f0; border-radius: 3px;")
                stats_layout.addWidget(cat_label)
            
            stats_layout.addStretch()
            self.prosite_layout.addWidget(stats_group)
        
        # Motifs table with enhanced columns
        motifs_group = QGroupBox("PROSITE Motifs")
        motifs_layout = QVBoxLayout(motifs_group)
        
        motifs_table = QTableWidget()
        motifs_table.setColumnCount(9)
        motifs_table.setHorizontalHeaderLabels([
            "Category", "ID", "Name", "Function", "Start", "End", "Length", "Sequence", "Source"
        ])
        
        motifs_table.setRowCount(len(results['motifs']))
        
        # Color coding for categories
        category_colors = {
            'Post-translational modification': QColor(255, 230, 230),
            'Nucleotide binding': QColor(230, 255, 230),
            'DNA/RNA binding': QColor(230, 230, 255),
            'Protein-protein interaction': QColor(255, 255, 230),
            'Enzymatic site': QColor(255, 230, 255),
            'Membrane protein': QColor(230, 255, 255),
            'Signal sequence': QColor(255, 245, 230),
            'Structural motif': QColor(245, 230, 255)
        }
        
        for i, motif in enumerate(results['motifs']):
            category = motif.get('category', '')
            
            motifs_table.setItem(i, 0, QTableWidgetItem(category))
            motifs_table.setItem(i, 1, QTableWidgetItem(motif.get('id', '')))
            motifs_table.setItem(i, 2, QTableWidgetItem(motif.get('name', '')))
            motifs_table.setItem(i, 3, QTableWidgetItem(motif.get('function', '')))
            motifs_table.setItem(i, 4, QTableWidgetItem(str(motif.get('start', ''))))
            motifs_table.setItem(i, 5, QTableWidgetItem(str(motif.get('end', ''))))
            motifs_table.setItem(i, 6, QTableWidgetItem(str(motif.get('length', ''))))
            motifs_table.setItem(i, 7, QTableWidgetItem(motif.get('sequence', '')))
            
            # Add source column
            source = motif.get('source', 'PROSITE API')
            source_item = QTableWidgetItem(source)
            source_item.setBackground(QColor(230, 240, 255))  # Light blue for API source
            source_item.setToolTip("Results from official PROSITE ScanProsite API")
            motifs_table.setItem(i, 8, source_item)
            
            # Color code by category (except source column which has its own colors)
            if category in category_colors:
                color = category_colors[category]
                for col in range(8):  # Skip the source column (index 8)
                    item = motifs_table.item(i, col)
                    if item:
                        item.setBackground(color)
        
        motifs_table.resizeColumnsToContents()
        motifs_table.setSortingEnabled(True)
        motifs_layout.addWidget(motifs_table)
        
        # Add context view for selected motifs
        context_group = QGroupBox("Sequence Context (click a motif to view)")
        context_layout = QVBoxLayout(context_group)
        
        self.context_display = QTextEdit()
        self.context_display.setMaximumHeight(100)
        self.context_display.setFont(QFont("Courier", 10))
        self.context_display.setPlainText("Select a motif from the table above to view its sequence context.")
        context_layout.addWidget(self.context_display)
        
        motifs_layout.addWidget(context_group)
        
        # Connect table selection to context display
        motifs_table.itemSelectionChanged.connect(lambda: self.show_motif_context(motifs_table, results['motifs']))
        
        self.prosite_layout.addWidget(motifs_group)
        self.prosite_layout.addStretch()
    
    def show_motif_context(self, table, motifs):
        """Show sequence context for selected motif."""
        current_row = table.currentRow()
        if current_row >= 0 and current_row < len(motifs):
            motif = motifs[current_row]
            
            # Get motif information
            motif_name = motif.get('name', 'Unknown motif')
            motif_id = motif.get('id', '')
            motif_start = motif.get('start', 0)
            motif_end = motif.get('end', 0)
            motif_sequence = motif.get('sequence', '')
            motif_description = motif.get('description', '')
            motif_category = motif.get('category', '')
            motif_function = motif.get('function', '')
            
            # Check if we have the current sequence and valid positions
            if hasattr(self, 'current_sequence') and self.current_sequence and motif_start > 0 and motif_end > 0:
                # Generate sequence context (20 amino acids before and after)
                context_window = 20
                seq_len = len(self.current_sequence)
                
                # Calculate context boundaries
                context_start = max(0, motif_start - 1 - context_window)  # Convert to 0-based indexing
                context_end = min(seq_len, motif_end + context_window)
                
                # Extract context sequence
                full_context = self.current_sequence[context_start:context_end]
                
                # Calculate relative positions within the context
                motif_rel_start = (motif_start - 1) - context_start  # Convert to 0-based
                motif_rel_end = motif_end - context_start
                
                # Split context into before, motif, and after
                before = full_context[:motif_rel_start]
                motif_seq = full_context[motif_rel_start:motif_rel_end]
                after = full_context[motif_rel_end:]
                
                # Format the context display
                context_info = []
                context_info.append(f"🔍 Motif: {motif_name}")
                if motif_id:
                    context_info.append(f"📋 ID: {motif_id}")
                if motif_description:
                    context_info.append(f"📝 Description: {motif_description}")
                if motif_category:
                    context_info.append(f"🏷️ Category: {motif_category}")
                if motif_function:
                    context_info.append(f"⚙️ Function: {motif_function}")
                context_info.append(f"📍 Position: {motif_start}-{motif_end} (length: {motif_end - motif_start + 1})")
                context_info.append("")
                context_info.append("🧬 Sequence Context:")
                
                # Create formatted sequence with position markers
                pos_start = context_start + 1  # Convert back to 1-based for display
                pos_end = context_end
                
                context_info.append(f"Position {pos_start}-{pos_end}:")
                context_info.append(f"{before}[{motif_seq}]{after}")
                context_info.append("")
                context_info.append("Legend: [motif sequence] in brackets")
                
                formatted_context = "\n".join(context_info)
                self.context_display.setPlainText(formatted_context)
            else:
                # Fallback display with available information
                fallback_info = []
                fallback_info.append(f"🔍 Motif: {motif_name}")
                if motif_id:
                    fallback_info.append(f"📋 ID: {motif_id}")
                if motif_description:
                    fallback_info.append(f"📝 Description: {motif_description}")
                if motif_category:
                    fallback_info.append(f"🏷️ Category: {motif_category}")
                if motif_function:
                    fallback_info.append(f"⚙️ Function: {motif_function}")
                if motif_start > 0 and motif_end > 0:
                    fallback_info.append(f"📍 Position: {motif_start}-{motif_end}")
                if motif_sequence:
                    fallback_info.append(f"🧬 Motif Sequence: {motif_sequence}")
                fallback_info.append("")
                fallback_info.append("⚠️ Full sequence context not available.")
                fallback_info.append("Please ensure a sequence is loaded for detailed context view.")
                
                self.context_display.setPlainText("\n".join(fallback_info))
    
    def save_prosite_results(self):
        """Save PROSITE results to file."""
        if not hasattr(self, 'current_prosite_results') or not self.current_prosite_results:
            QMessageBox.warning(self, "No Results", "No PROSITE results to save.")
            return
        
        # Get save location
        file_path, selected_filter = QFileDialog.getSaveFileName(
            self, "Save PROSITE Results", 
            f"prosite_results_{time.strftime('%Y%m%d_%H%M%S')}",
            "CSV Files (*.csv);;JSON Files (*.json);;Text Files (*.txt);;All Files (*)"
        )
        
        if not file_path:
            return
        
        try:
            results = self.current_prosite_results
            
            # Add proper file extension
            if selected_filter.startswith("CSV") and not file_path.endswith('.csv'):
                file_path += '.csv'
            elif selected_filter.startswith("JSON") and not file_path.endswith('.json'):
                file_path += '.json'
            elif selected_filter.startswith("Text") and not file_path.endswith('.txt'):
                file_path += '.txt'
            
            if selected_filter.startswith("CSV"):
                self.save_prosite_csv(file_path, results)
            elif selected_filter.startswith("JSON"):
                self.save_prosite_json(file_path, results)
            else:
                self.save_prosite_text(file_path, results)
            
            QMessageBox.information(self, "Success", f"PROSITE results saved to {file_path}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save results:\n{str(e)}")
    
    def save_prosite_csv(self, file_path, results):
        """Save results as CSV."""
        import csv
        
        with open(file_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Header
            writer.writerow([
                'ID', 'Name', 'Category', 'Function', 'Description', 
                'Start', 'End', 'Length', 'Sequence', 'Context'
            ])
            
            # Data
            for motif in results['motifs']:
                context = motif.get('context', {})
                context_str = context.get('sequence', '') if context else ''
                
                writer.writerow([
                    motif.get('id', ''),
                    motif.get('name', ''),
                    motif.get('category', ''),
                    motif.get('function', ''),
                    motif.get('description', ''),
                    motif.get('start', ''),
                    motif.get('end', ''),
                    motif.get('length', ''),
                    motif.get('sequence', ''),
                    context_str
                ])
    
    def save_prosite_json(self, file_path, results):
        """Save results as JSON."""
        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
    
    def save_prosite_text(self, file_path, results):
        """Save results as formatted text."""
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write("PROSITE Motif Analysis Results\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Analysis Date: {results.get('search_date', 'Unknown')}\n")
            f.write(f"Sequence Length: {results.get('sequence_length', 'Unknown')} amino acids\n")
            f.write(f"Total Motifs Found: {len(results['motifs'])}\n\n")
            
            # Category summary
            if results.get('categories'):
                f.write("Motifs by Category:\n")
                for category, count in results['categories'].items():
                    f.write(f"  {category}: {count}\n")
                f.write("\n")
            
            # Detailed results
            f.write("Detailed Results:\n")
            f.write("-" * 50 + "\n\n")
            
            for i, motif in enumerate(results['motifs'], 1):
                f.write(f"{i}. {motif.get('name', 'Unknown')} ({motif.get('id', '')})\n")
                f.write(f"   Category: {motif.get('category', 'Unknown')}\n")
                f.write(f"   Function: {motif.get('function', 'Unknown')}\n")
                f.write(f"   Position: {motif.get('start', '')}-{motif.get('end', '')}\n")
                f.write(f"   Sequence: {motif.get('sequence', '')}\n")
                f.write(f"   Description: {motif.get('description', '')}\n")
                
                context = motif.get('context', {})
                if context:
                    f.write(f"   Context: {context.get('sequence', '')}\n")
                
                f.write("\n")
    
    def open_interpro_job(self, job_id):
        """Open InterPro job results in web browser."""
        if not job_id:
            QMessageBox.warning(self, "No Job ID", "No InterPro job ID available.")
            return
        
        # Construct InterPro results URL
        interpro_url = f"https://www.ebi.ac.uk/interpro/result/InterProScan/{job_id}/"
        
        try:
            # Open URL in default web browser
            QDesktopServices.openUrl(QUrl(interpro_url))
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open InterPro job in browser:\n{str(e)}")
    
    def save_interpro_results(self):
        """Save InterPro results to file."""
        if not hasattr(self, 'current_interpro_results') or not self.current_interpro_results:
            QMessageBox.warning(self, "No Results", "No InterPro results to save.")
            return
        
        # Get save location
        file_path, selected_filter = QFileDialog.getSaveFileName(
            self, "Save InterPro Results", 
            f"interpro_results_{time.strftime('%Y%m%d_%H%M%S')}",
            "CSV Files (*.csv);;JSON Files (*.json);;Text Files (*.txt);;All Files (*)"
        )
        
        if not file_path:
            return
        
        try:
            results = self.current_interpro_results
            
            # Add proper file extension
            if selected_filter.startswith("CSV") and not file_path.endswith('.csv'):
                file_path += '.csv'
            elif selected_filter.startswith("JSON") and not file_path.endswith('.json'):
                file_path += '.json'
            elif selected_filter.startswith("Text") and not file_path.endswith('.txt'):
                file_path += '.txt'
            
            if selected_filter.startswith("CSV"):
                self.save_interpro_csv(file_path, results)
            elif selected_filter.startswith("JSON"):
                self.save_interpro_json(file_path, results)
            else:
                self.save_interpro_text(file_path, results)
            
            QMessageBox.information(self, "Success", f"InterPro results saved to {file_path}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save results:\n{str(e)}")
    
    def save_interpro_csv(self, file_path, results):
        """Save InterPro results as CSV."""
        import csv
        
        with open(file_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Header
            writer.writerow([
                'Database', 'Accession', 'Name', 'Description', 'Start', 'End', 'E-value'
            ])
            
            # Data
            for domain in results['domains']:
                locations = domain.get('locations', [])
                start = locations[0].get('start', '') if locations else ''
                end = locations[0].get('end', '') if locations else ''
                
                writer.writerow([
                    domain.get('database', ''),
                    domain.get('accession', ''),
                    domain.get('name', ''),
                    domain.get('description', ''),
                    start, end,
                    domain.get('evalue', 'N/A')
                ])
    
    def save_interpro_json(self, file_path, results):
        """Save InterPro results as JSON."""
        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
    
    def save_interpro_text(self, file_path, results):
        """Save InterPro results as formatted text."""
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write("InterPro Domain Analysis Results\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total Domains Found: {len(results['domains'])}\n\n")
            
            # Count by database
            db_counts = {}
            for domain in results['domains']:
                db = domain.get('database', 'Unknown')
                db_counts[db] = db_counts.get(db, 0) + 1
            
            if db_counts:
                f.write("Domains by Database:\n")
                for db, count in db_counts.items():
                    f.write(f"  {db}: {count}\n")
                f.write("\n")
            
            # Detailed results
            f.write("Detailed Results:\n")
            f.write("-" * 50 + "\n\n")
            
            for i, domain in enumerate(results['domains'], 1):
                locations = domain.get('locations', [])
                start = locations[0].get('start', 'N/A') if locations else 'N/A'
                end = locations[0].get('end', 'N/A') if locations else 'N/A'
                
                f.write(f"{i}. {domain.get('name', 'Unknown')} ({domain.get('accession', '')})\n")
                f.write(f"   Database: {domain.get('database', 'Unknown')}\n")
                f.write(f"   Position: {start}-{end}\n")
                f.write(f"   E-value: {domain.get('evalue', 'N/A')}\n")
                f.write(f"   Description: {domain.get('description', 'N/A')}\n\n")
            
            # GO terms if available
            if results.get('go_terms'):
                f.write("\nGene Ontology Terms:\n")
                f.write("-" * 30 + "\n")
                for go_term in results['go_terms']:
                    f.write(f"  {go_term.get('id', '')}: {go_term.get('name', '')} ({go_term.get('category', '')})\n")
    
    def update_visualization(self):
        """Update the domain visualization with current results."""
        if not self.current_sequence:
            return
        
        # Collect domains from InterPro results
        domains = []
        if hasattr(self, 'current_interpro_results') and self.current_interpro_results:
            for domain in self.current_interpro_results.get('domains', []):
                locations = domain.get('locations', [])
                if locations:
                    viz_domain = {
                        'name': domain.get('name', 'Unknown'),
                        'start': locations[0].get('start', 0),
                        'end': locations[0].get('end', 0),
                        'database': domain.get('database', ''),
                        'description': domain.get('description', '')
                    }
                    domains.append(viz_domain)
        
        # Collect motifs from PROSITE results
        motifs = []
        if hasattr(self, 'current_prosite_results') and self.current_prosite_results:
            for motif in self.current_prosite_results.get('motifs', []):
                viz_motif = {
                    'name': motif.get('name', 'Unknown'),
                    'start': motif.get('start', 0),
                    'end': motif.get('end', 0),
                    'category': motif.get('category', ''),
                    'sequence': motif.get('sequence', '')
                }
                motifs.append(viz_motif)
        
        # Update visualization
        self.domain_viz.set_data(len(self.current_sequence), domains, motifs)
    
    def update_save_button_state(self):
        """Enable save all button if any results are available."""
        has_interpro = hasattr(self, 'current_interpro_results') and bool(self.current_interpro_results)
        has_prosite = hasattr(self, 'current_prosite_results') and bool(self.current_prosite_results)
        
        self.save_all_btn.setEnabled(has_interpro or has_prosite)
    
    def save_all_results(self):
        """Save all analysis results to a comprehensive report."""
        has_interpro = hasattr(self, 'current_interpro_results') and bool(self.current_interpro_results)
        has_prosite = hasattr(self, 'current_prosite_results') and bool(self.current_prosite_results)
        
        if not (has_interpro or has_prosite):
            QMessageBox.warning(self, "No Results", "No analysis results to save.")
            return
        
        # Get save location
        file_path, selected_filter = QFileDialog.getSaveFileName(
            self, "Save All Analysis Results", 
            f"motif_analysis_complete_{time.strftime('%Y%m%d_%H%M%S')}",
            "HTML Report (*.html);;JSON Files (*.json);;CSV Files (*.csv);;Text Files (*.txt);;All Files (*)"
        )
        
        if not file_path:
            return
        
        try:
            # Add proper file extension
            if selected_filter.startswith("HTML") and not file_path.endswith('.html'):
                file_path += '.html'
            elif selected_filter.startswith("JSON") and not file_path.endswith('.json'):
                file_path += '.json'
            elif selected_filter.startswith("CSV") and not file_path.endswith('.csv'):
                file_path += '.csv'
            elif selected_filter.startswith("Text") and not file_path.endswith('.txt'):
                file_path += '.txt'
            
            if selected_filter.startswith("HTML"):
                self.save_all_html(file_path)
            elif selected_filter.startswith("JSON"):
                self.save_all_json(file_path)
            elif selected_filter.startswith("CSV"):
                self.save_all_csv(file_path)
            else:
                self.save_all_text(file_path)
            
            QMessageBox.information(self, "Success", f"Complete analysis report saved to {file_path}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save results:\n{str(e)}")
    
    def save_all_html(self, file_path):
        """Save all results as an HTML report."""
        html_content = self.generate_html_report()
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
    
    def save_all_json(self, file_path):
        """Save all results as JSON."""
        combined_results = {
            'analysis_info': {
                'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                'sequence_length': len(self.current_sequence) if self.current_sequence else 0,
                'sequence': self.current_sequence if self.current_sequence else ''
            },
            'interpro_results': getattr(self, 'current_interpro_results', None),
            'prosite_results': getattr(self, 'current_prosite_results', None)
        }
        
        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(combined_results, f, indent=2, ensure_ascii=False)
    
    def save_all_csv(self, file_path):
        """Save all results as CSV (multiple sheets concept in single file)."""
        import csv
        
        with open(file_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Analysis summary
            writer.writerow(['MOTIF ANALYSIS SUMMARY'])
            writer.writerow(['Analysis Date', time.strftime('%Y-%m-%d %H:%M:%S')])
            writer.writerow(['Sequence Length', len(self.current_sequence) if self.current_sequence else 0])
            writer.writerow([])
            
            # InterPro results
            if hasattr(self, 'current_interpro_results') and self.current_interpro_results:
                writer.writerow(['INTERPRO DOMAINS'])
                writer.writerow(['Database', 'Accession', 'Name', 'Description', 'Start', 'End', 'E-value'])
                
                for domain in self.current_interpro_results.get('domains', []):
                    locations = domain.get('locations', [])
                    start = locations[0].get('start', '') if locations else ''
                    end = locations[0].get('end', '') if locations else ''
                    
                    writer.writerow([
                        domain.get('database', ''),
                        domain.get('accession', ''),
                        domain.get('name', ''),
                        domain.get('description', ''),
                        start, end,
                        domain.get('evalue', 'N/A')
                    ])
                writer.writerow([])
            
            # PROSITE results
            if hasattr(self, 'current_prosite_results') and self.current_prosite_results:
                writer.writerow(['PROSITE MOTIFS'])
                writer.writerow(['Category', 'ID', 'Name', 'Function', 'Start', 'End', 'Length', 'Sequence'])
                
                for motif in self.current_prosite_results.get('motifs', []):
                    writer.writerow([
                        motif.get('category', ''),
                        motif.get('id', ''),
                        motif.get('name', ''),
                        motif.get('function', ''),
                        motif.get('start', ''),
                        motif.get('end', ''),
                        motif.get('length', ''),
                        motif.get('sequence', '')
                    ])
    
    def save_all_text(self, file_path):
        """Save all results as formatted text."""
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write("COMPREHENSIVE MOTIF ANALYSIS REPORT\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Sequence Length: {len(self.current_sequence) if self.current_sequence else 0} amino acids\n\n")
            
            # InterPro section
            if hasattr(self, 'current_interpro_results') and self.current_interpro_results:
                interpro = self.current_interpro_results
                f.write("INTERPRO DOMAIN ANALYSIS\n")
                f.write("-" * 40 + "\n\n")
                
                f.write(f"Total domains found: {len(interpro.get('domains', []))}\n")
                
                # Count by database
                db_counts = {}
                for domain in interpro.get('domains', []):
                    db = domain.get('database', 'Unknown')
                    db_counts[db] = db_counts.get(db, 0) + 1
                
                if db_counts:
                    f.write("\nDomains by database:\n")
                    for db, count in db_counts.items():
                        f.write(f"  {db}: {count}\n")
                
                f.write("\nDetailed domain list:\n")
                for i, domain in enumerate(interpro.get('domains', []), 1):
                    locations = domain.get('locations', [])
                    start = locations[0].get('start', 'N/A') if locations else 'N/A'
                    end = locations[0].get('end', 'N/A') if locations else 'N/A'
                    
                    f.write(f"{i}. {domain.get('name', 'Unknown')} ({domain.get('accession', '')})\n")
                    f.write(f"   Database: {domain.get('database', 'Unknown')}\n")
                    f.write(f"   Position: {start}-{end}\n")
                    f.write(f"   Description: {domain.get('description', 'N/A')}\n\n")
            
            # PROSITE section
            if hasattr(self, 'current_prosite_results') and self.current_prosite_results:
                prosite = self.current_prosite_results
                f.write("\nPROSITE MOTIF ANALYSIS\n")
                f.write("-" * 40 + "\n\n")
                
                f.write(f"Total motifs found: {len(prosite.get('motifs', []))}\n")
                
                # Category summary
                if prosite.get('categories'):
                    f.write("\nMotifs by category:\n")
                    for category, count in prosite['categories'].items():
                        f.write(f"  {category}: {count}\n")
                
                f.write("\nDetailed motif list:\n")
                for i, motif in enumerate(prosite.get('motifs', []), 1):
                    f.write(f"{i}. {motif.get('name', 'Unknown')} ({motif.get('id', '')})\n")
                    f.write(f"   Category: {motif.get('category', 'Unknown')}\n")
                    f.write(f"   Function: {motif.get('function', 'Unknown')}\n")
                    f.write(f"   Position: {motif.get('start', '')}-{motif.get('end', '')}\n")
                    f.write(f"   Sequence: {motif.get('sequence', '')}\n")
                    f.write(f"   Description: {motif.get('description', '')}\n\n")
    
    def generate_html_report(self):
        """Generate a comprehensive HTML report."""
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PicoMol Motif Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }}
        .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
        .section {{ margin-bottom: 30px; }}
        .section h2 {{ color: #333; border-bottom: 2px solid #ddd; padding-bottom: 5px; }}
        .section h3 {{ color: #555; margin-top: 20px; margin-bottom: 10px; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 10px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; font-weight: bold; }}
        .summary {{ background-color: #f9f9f9; padding: 15px; border-radius: 5px; margin-bottom: 15px; }}
        .interpro-entry {{ background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; margin: 10px 0; padding: 15px; }}
        .interpro-header {{ color: #2196F3; font-weight: bold; margin-bottom: 10px; }}
        .database-section {{ margin: 15px 0; }}
        .database-header {{ font-weight: bold; padding: 5px 0; }}
        .pfam {{ color: #FF6B35; }}
        .smart {{ color: #4CAF50; }}
        .prosite {{ color: #2196F3; }}
        .go-term {{ margin: 5px 0; padding: 5px; background-color: #f5f5f5; border-radius: 3px; }}
        .go-id {{ font-family: monospace; color: #2196F3; font-weight: bold; }}
        .pathway {{ margin: 5px 0; padding: 5px; background-color: #fff5f5; border-radius: 3px; }}
        .pathway-id {{ font-family: monospace; color: #FF6B35; font-weight: bold; }}
        .motif-category {{ padding: 3px 8px; border-radius: 3px; font-size: 0.9em; }}
        .ptm {{ background-color: #ffe6e6; }}
        .nucleotide {{ background-color: #e6ffe6; }}
        .dna {{ background-color: #e6e6ff; }}
        .protein {{ background-color: #ffffe6; }}
        .enzyme {{ background-color: #ffe6ff; }}
        .membrane {{ background-color: #e6ffff; }}
        .signal {{ background-color: #fff5e6; }}
        .structural {{ background-color: #f5e6ff; }}
        .significant {{ background-color: #d4edda; }}
        .very-significant {{ background-color: #c8e6c9; }}
        code {{ background-color: #f8f9fa; padding: 2px 4px; border-radius: 3px; font-family: monospace; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 PicoMol Motif Analysis Report</h1>
        <p><strong>Analysis Date:</strong> {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Sequence Length:</strong> {len(self.current_sequence) if self.current_sequence else 0} amino acids</p>
    </div>
"""
        
        # InterPro section
        if hasattr(self, 'current_interpro_results') and self.current_interpro_results:
            interpro = self.current_interpro_results
            
            # Count different types of matches
            interpro_entries = set()
            member_db_matches = 0
            db_groups = {}
            
            for domain in interpro.get('domains', []):
                if domain.get('interpro'):
                    interpro_entries.add(domain['interpro'].get('accession', ''))
                member_db_matches += 1
                
                db = domain.get('database', 'Unknown')
                if db not in db_groups:
                    db_groups[db] = []
                db_groups[db].append(domain)
            
            html += f"""
    <div class="section">
        <h2>🔍 InterPro Analysis Results</h2>
        <div class="summary">
            <p><strong>Job ID:</strong> {interpro.get('job_id', 'N/A')}</p>
            <p><strong>InterPro Entries:</strong> {len(interpro_entries)}</p>
            <p><strong>Member Database Matches:</strong> {member_db_matches}</p>
            <p><strong>Databases Searched:</strong> {', '.join(sorted(db_groups.keys()))}</p>
        </div>
"""
            
            # InterPro Entries section
            if interpro_entries:
                html += "<h3>InterPro Entries</h3>"
                interpro_groups = {}
                for domain in interpro.get('domains', []):
                    interpro_entry = domain.get('interpro')
                    if interpro_entry and interpro_entry.get('accession'):
                        acc = interpro_entry['accession']
                        if acc not in interpro_groups:
                            interpro_groups[acc] = {
                                'entry': interpro_entry,
                                'matches': []
                            }
                        interpro_groups[acc]['matches'].append(domain)
                
                for acc, group_data in sorted(interpro_groups.items()):
                    entry = group_data['entry']
                    matches = group_data['matches']
                    
                    html += f"""
        <div class="interpro-entry">
            <div class="interpro-header">{acc} - {entry.get('name', '')}</div>
            <p><em>{entry.get('description', '')}</em></p>
            <p><strong>Type:</strong> {entry.get('type', 'N/A')} | <strong>Member matches:</strong> {len(matches)}</p>
            <table>
                <tr>
                    <th>Database</th>
                    <th>Accession</th>
                    <th>Name</th>
                    <th>Start</th>
                    <th>End</th>
                    <th>E-value</th>
                </tr>
"""
                    
                    for match in matches:
                        locations = match.get('locations', [])
                        start = locations[0].get('start', 'N/A') if locations else 'N/A'
                        end = locations[0].get('end', 'N/A') if locations else 'N/A'
                        
                        # Add significance highlighting
                        row_class = ""
                        evalue = match.get('evalue', 'N/A')
                        if evalue != 'N/A':
                            try:
                                eval_float = float(evalue)
                                if eval_float < 1e-10:
                                    row_class = ' class="very-significant"'
                                elif eval_float < 1e-5:
                                    row_class = ' class="significant"'
                            except ValueError:
                                pass
                        
                        html += f"""
                <tr{row_class}>
                    <td><span class="{match.get('database', '').lower()}">{match.get('database', '')}</span></td>
                    <td>{match.get('accession', '')}</td>
                    <td>{match.get('name', '')}</td>
                    <td>{start}</td>
                    <td>{end}</td>
                    <td>{evalue}</td>
                </tr>
"""
                    
                    html += "            </table>\n        </div>\n"
            
            # Member Database Matches section
            html += "<h3>Member Database Matches by Database</h3>"
            for db_name, domains in sorted(db_groups.items()):
                html += f"""
        <div class="database-section">
            <div class="database-header {db_name.lower()}">{db_name} ({len(domains)} matches)</div>
            <table>
                <tr>
                    <th>Accession</th>
                    <th>Name</th>
                    <th>Description</th>
                    <th>Start</th>
                    <th>End</th>
                    <th>E-value</th>
                </tr>
"""
                
                for domain in domains:
                    locations = domain.get('locations', [])
                    start = locations[0].get('start', 'N/A') if locations else 'N/A'
                    end = locations[0].get('end', 'N/A') if locations else 'N/A'
                    
                    # Add significance highlighting
                    row_class = ""
                    evalue = domain.get('evalue', 'N/A')
                    if evalue != 'N/A':
                        try:
                            eval_float = float(evalue)
                            if eval_float < 1e-10:
                                row_class = ' class="very-significant"'
                            elif eval_float < 1e-5:
                                row_class = ' class="significant"'
                        except ValueError:
                            pass
                    
                    html += f"""
                <tr{row_class}>
                    <td>{domain.get('accession', '')}</td>
                    <td>{domain.get('name', '')}</td>
                    <td>{domain.get('description', '')}</td>
                    <td>{start}</td>
                    <td>{end}</td>
                    <td>{evalue}</td>
                </tr>
"""
                
                html += "            </table>\n        </div>\n"
            
            # GO Terms section
            if interpro.get('go_terms'):
                html += "<h3>Gene Ontology Terms</h3>"
                go_categories = {}
                for go_term in interpro['go_terms']:
                    category = go_term.get('category', 'unknown')
                    if category not in go_categories:
                        go_categories[category] = []
                    go_categories[category].append(go_term)
                
                for category, terms in sorted(go_categories.items()):
                    html += f"<h4>{category.replace('_', ' ').title()} ({len(terms)})</h4>"
                    for term in terms:
                        html += f"""
        <div class="go-term">
            <span class="go-id">{term.get('id', '')}</span> - {term.get('name', '')}
        </div>
"""
            
            # Pathways section
            if interpro.get('pathways'):
                html += "<h3>Pathways</h3>"
                for pathway in interpro['pathways']:
                    html += f"""
        <div class="pathway">
            <span class="pathway-id">{pathway.get('id', '')}</span> - {pathway.get('name', '')}
            {f' [{pathway["database"]}]' if pathway.get('database') else ''}
        </div>
"""
            
            html += "    </div>\n"
        
        # PROSITE section
        if hasattr(self, 'current_prosite_results') and self.current_prosite_results:
            prosite = self.current_prosite_results
            html += f"""
    <div class="section">
        <h2>🎯 PROSITE Motif Analysis</h2>
        <div class="summary">
            <p><strong>Total motifs found:</strong> {len(prosite.get('motifs', []))}</p>
            <p><strong>Functional categories:</strong> {len(prosite.get('categories', {}))}</p>
        </div>
"""
            
            if prosite.get('motifs'):
                # Add API source info
                html += """
        <div style="margin: 10px 0; padding: 10px; background-color: #e3f2fd; border-radius: 5px; border-left: 4px solid #2196f3;">
            <strong>🔗 Data Source:</strong> Official PROSITE ScanProsite API<br>
            <small>Results obtained directly from the PROSITE database at ExPASy, ensuring accuracy and up-to-date pattern matching.</small>
        </div>
"""
                
                # Group motifs by category
                category_groups = {}
                for motif in prosite['motifs']:
                    category = motif.get('category', 'Unknown')
                    if category not in category_groups:
                        category_groups[category] = []
                    category_groups[category].append(motif)
                
                category_classes = {
                    'Post-translational modification': 'ptm',
                    'Nucleotide binding': 'nucleotide',
                    'DNA/RNA binding': 'dna',
                    'Protein-protein interaction': 'protein',
                    'Enzymatic site': 'enzyme',
                    'Membrane protein': 'membrane',
                    'Signal sequence': 'signal',
                    'Structural motif': 'structural'
                }
                
                for category, motifs in sorted(category_groups.items()):
                    css_class = category_classes.get(category, '')
                    html += f"""
        <h3><span class="motif-category {css_class}">{category}</span> ({len(motifs)} motifs)</h3>
        <table>
            <tr>
                <th>ID</th>
                <th>Name</th>
                <th>Pattern</th>
                <th>Function</th>
                <th>Position</th>
                <th>Sequence</th>
            </tr>
"""
                    
                    for motif in motifs:
                        # Use consistent styling for API results
                        row_color = '#f0f8ff'  # Light blue for API results
                        
                        html += f"""
            <tr style="background-color: {row_color};">
                <td>{motif.get('id', '')}</td>
                <td>{motif.get('name', '')}</td>
                <td><code>{motif.get('pattern', motif.get('sequence', ''))}</code></td>
                <td>{motif.get('function', '')}</td>
                <td>{motif.get('start', '')}-{motif.get('end', '')}</td>
                <td><code>{motif.get('sequence', '')}</code></td>
            </tr>
"""
                    
                    html += "        </table>\n"
            
            html += "    </div>\n"
        
        html += """
    <div class="section">
        <p><em>Report generated by PicoMol - Molecular Visualization and Bioinformatics Suite</em></p>
        <p><small>InterPro data from EMBL-EBI | PROSITE patterns from SIB Swiss Institute of Bioinformatics</small></p>
    </div>
</body>
</html>
"""
        
        return html


def create_motif_analysis_tab(parent=None):
    """Create the protein motif and domain analysis tab."""
    return MotifAnalysisTab(parent)


if __name__ == "__main__":
    # Test the motif analysis tools
    import sys
    from PyQt5.QtWidgets import QApplication
    
    app = QApplication(sys.argv)
    
    widget = create_motif_analysis_tab()
    widget.setWindowTitle("PicoMol Motif Analysis Tools")
    widget.resize(1000, 700)
    widget.show()
    
    sys.exit(app.exec_())