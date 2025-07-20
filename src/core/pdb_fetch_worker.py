#!/usr/bin/env python3
"""
PDB Fetch Worker for non-blocking PDB data retrieval.

This module provides a QThread-based worker for fetching PDB data
without blocking the main UI thread.
"""

from PyQt5.QtCore import QThread, pyqtSignal
from typing import Dict, Any, Optional
import traceback

try:
    from .enhanced_pdb_puller_async import OptimizedPDBPuller
    OPTIMIZED_PULLER_AVAILABLE = True
except ImportError:
    OPTIMIZED_PULLER_AVAILABLE = False


class PDBFetchWorker(QThread):
    """Worker thread for fetching PDB data without blocking the UI."""
    
    # Signals
    progress_update = pyqtSignal(str)  # Progress message
    data_fetched = pyqtSignal(dict)    # Fetched data
    error_occurred = pyqtSignal(str)   # Error message
    finished = pyqtSignal()            # Completion signal
    
    def __init__(self, puller: OptimizedPDBPuller, pdb_id: str, 
                 include_validation: bool = False,
                 include_sequences: bool = True,
                 include_mmcif: bool = False):
        """
        Initialize the PDB fetch worker.
        
        Args:
            puller: OptimizedPDBPuller instance
            pdb_id: PDB ID to fetch
            include_validation: Whether to fetch validation data
            include_sequences: Whether to fetch sequences
            include_mmcif: Whether to fetch mmCIF files
        """
        super().__init__()
        self.puller = puller
        self.pdb_id = pdb_id
        self.include_validation = include_validation
        self.include_sequences = include_sequences
        self.include_mmcif = include_mmcif
        self._is_cancelled = False
    
    def run(self):
        """Run the PDB fetch operation."""
        try:
            # Define progress callback
            def progress_callback(message: str):
                if not self._is_cancelled:
                    self.progress_update.emit(message)
            
            # Fetch the data
            data = self.puller.fetch_pdb_data_async(
                self.pdb_id,
                progress_callback=progress_callback,
                include_validation=self.include_validation,
                include_sequences=self.include_sequences,
                include_mmcif=self.include_mmcif
            )
            
            if not self._is_cancelled:
                if data.get('errors'):
                    error_msg = "\\n".join(data['errors'])
                    self.error_occurred.emit(f"Errors occurred:\\n{error_msg}")
                else:
                    self.data_fetched.emit(data)
            
        except Exception as e:
            if not self._is_cancelled:
                error_msg = f"Failed to fetch PDB data: {str(e)}"
                # Add traceback for debugging
                traceback_str = traceback.format_exc()
                self.error_occurred.emit(f"{error_msg}\\n\\nDetails:\\n{traceback_str}")
        
        finally:
            if not self._is_cancelled:
                self.finished.emit()
    
    def cancel(self):
        """Cancel the fetch operation."""
        self._is_cancelled = True
        self.quit()
        self.wait(1000)  # Wait up to 1 second for thread to finish


class PDBFetchManager:
    """Manager for PDB fetch operations with UI integration."""
    
    def __init__(self, data_dir: str):
        """
        Initialize the PDB fetch manager.
        
        Args:
            data_dir: Directory for storing PDB data
        """
        if not OPTIMIZED_PULLER_AVAILABLE:
            raise ImportError("Optimized PDB puller not available")
        
        self.puller = OptimizedPDBPuller(data_dir)
        self.current_worker = None
    
    def fetch_pdb_async(self, pdb_id: str,
                       progress_callback: Optional[callable] = None,
                       success_callback: Optional[callable] = None,
                       error_callback: Optional[callable] = None,
                       finished_callback: Optional[callable] = None,
                       include_validation: bool = False,
                       include_sequences: bool = True,
                       include_mmcif: bool = False):
        """
        Fetch PDB data asynchronously with callbacks.
        
        Args:
            pdb_id: PDB ID to fetch
            progress_callback: Called with progress messages
            success_callback: Called with fetched data on success
            error_callback: Called with error message on failure
            finished_callback: Called when operation completes
            include_validation: Whether to fetch validation data
            include_sequences: Whether to fetch sequences
            include_mmcif: Whether to fetch mmCIF files
        """
        # Cancel any existing operation
        self.cancel_current_fetch()
        
        # Create new worker
        self.current_worker = PDBFetchWorker(
            self.puller, pdb_id,
            include_validation=include_validation,
            include_sequences=include_sequences,
            include_mmcif=include_mmcif
        )
        
        # Connect callbacks
        if progress_callback:
            self.current_worker.progress_update.connect(progress_callback)
        
        if success_callback:
            self.current_worker.data_fetched.connect(success_callback)
        
        if error_callback:
            self.current_worker.error_occurred.connect(error_callback)
        
        if finished_callback:
            self.current_worker.finished.connect(finished_callback)
        
        # Clean up worker when finished
        self.current_worker.finished.connect(self._cleanup_worker)
        
        # Start the worker
        self.current_worker.start()
    
    def cancel_current_fetch(self):
        """Cancel the current fetch operation if any."""
        if self.current_worker and self.current_worker.isRunning():
            self.current_worker.cancel()
            self.current_worker = None
    
    def _cleanup_worker(self):
        """Clean up the worker thread."""
        if self.current_worker:
            self.current_worker.deleteLater()
            self.current_worker = None
    
    def is_fetching(self) -> bool:
        """Check if a fetch operation is currently running."""
        return self.current_worker is not None and self.current_worker.isRunning()
    
    def get_structure_info(self, pdb_id: str) -> Dict[str, Any]:
        """Get cached structure information."""
        return self.puller.get_structure_info(pdb_id)
    
    def list_downloaded_structures(self) -> list:
        """List all downloaded structures."""
        return self.puller.list_downloaded_structures()
    
    def clear_cache(self):
        """Clear the API response cache."""
        self.puller.clear_cache()
    
    def get_cache_size(self) -> int:
        """Get the current cache size."""
        return self.puller.get_cache_size()