import functools
import os
import sys
import threading
import socketserver
import http.server
import time
import shutil
import subprocess
import sys
import warnings
import webbrowser

from src.blast_utils import (
    create_ncbi_style_blastp_tab, create_ncbi_style_blastn_tab,
    create_ncbi_style_blastx_tab, create_ncbi_style_tblastn_tab,
    create_ncbi_style_tblastx_tab
)

from src.core.preferences import PreferencesDialog, PreferencesManager
from src.gui.theme_manager import apply_theme
from src.core.bioinformatics_tools import create_bioinformatics_tab
from src.gui.welcome_dialog import WelcomeDialog

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QFileDialog, QMessageBox,
    QComboBox, QCheckBox, QGroupBox, QTextEdit, QDialog, QDialogButtonBox, QAction,
    QTabWidget, QSizePolicy, QColorDialog, QFormLayout, QScrollArea
)


class AboutDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About PicoMol")
        self.setMinimumWidth(600)
        self.setMinimumHeight(700)
        
        # Create scroll area to handle overflow
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        # Create the widget that contains the content
        content = QWidget()
        layout = QVBoxLayout(content)
        
        # Add logo at the top
        try:
            project_root = os.path.dirname(os.path.abspath(__file__))
            logo_path = os.path.join(project_root, "PicoMol.png")
            if os.path.exists(logo_path):
                logo_label = QLabel()
                pixmap = QPixmap(logo_path)
                if not pixmap.isNull():
                    # Scale the logo to a reasonable size for the dialog
                    scaled_pixmap = pixmap.scaled(120, 120, Qt.KeepAspectRatio, Qt.SmoothTransformation)
                    logo_label.setPixmap(scaled_pixmap)
                    logo_label.setAlignment(Qt.AlignCenter)
                    layout.addWidget(logo_label)
        except Exception as e:
            print(f"Error adding logo to About dialog: {e}")
        
        about_text = QLabel(
            """
            <style>
                body { font-family: Arial, sans-serif; line-height: 1.5; }
                h2 { color: #2c3e50; margin-top: 0.5em; }
                h3 { color: #3498db; margin-top: 1em; }
                a { color: #2980b9; text-decoration: none; }
                a:hover { text-decoration: underline; }
                ul { margin: 0.5em 0; padding-left: 1.5em; }
                li { margin-bottom: 0.3em; }
                .citation { font-size: 0.9em; margin: 0.5em 0; }
            </style>
            
            <h2>🧬 PicoMol</h2>
            <p><b>Version:</b> 0.1.0 (2025-07-24)</p>
            <p>A comprehensive molecular visualization and bioinformatics suite for protein structures, sequence analysis, and functional annotation.</p>
            
            <h3>🔬 Core Features</h3>
            <ul>
              <li><b>3D Molecular Visualization:</b> Interactive protein structure viewing with multiple representations</li>
              <li><b>BLAST Integration:</b> Complete BLAST suite (BLASTP, BLASTN, BLASTX, TBLASTN, TBLASTX)</li>
              <li><b>Motif & Domain Analysis:</b> InterPro and PROSITE integration for comprehensive functional annotation</li>
              <li><b>PDB Support:</b> Fetch from PDB database or load local files</li>
              <li><b>Sequence Analysis:</b> View and analyze protein/nucleotide sequences</li>
              <li><b>Comprehensive Data Export:</b> Export analysis results in multiple formats (HTML, JSON, CSV, Excel, Text)</li>
            </ul>
            
            <h3>🛠️ Technology Stack</h3>
            <ul>
              <li><b>GUI Framework:</b> <a href='https://riverbankcomputing.com/software/pyqt/intro'>PyQt5</a> with QWebEngine</li>
              <li><b>3D Visualization:</b> <a href='https://github.com/arose/ngl'>NGL.js</a> (MIT License)</li>
              <li><b>Bioinformatics:</b> <a href='https://biopython.org/'>Biopython</a></li>
              <li><b>BLAST:</b> NCBI BLAST API integration</li>
              <li><b>InterPro:</b> EBI InterPro web service API</li>
              <li><b>PROSITE:</b> ExPASy PROSITE ScanProsite API</li>
            </ul>
            
            <p><b>Developer:</b> Jack Magson<br>
            <b>License:</b> GNU GPL v3.0 (see LICENSE file)<br>
            <b>Repository:</b> <a href='https://github.com/ultimafounding/PicoMol'>GitHub</a></p>
            
            <h3>📚 Citations</h3>
            
            <p><b>NGL.js:</b></p>
            <div class="citation">
            Rose, A. S., Bradley, A. R., Valasatava, Y., Duarte, J. M., Prlić, A., & Rose, P. W. (2018). NGL viewer: web-based molecular graphics for large complexes. 
            <i>Bioinformatics</i>, 34(21), 3755-3758. 
            <a href='https://doi.org/10.1093/bioinformatics/bty419'>doi:10.1093/bioinformatics/bty419</a>
            </div>
            
            <div class="citation">
            Rose, A. S., & Hildebrand, P. W. (2015). NGL Viewer: a web application for molecular visualization. 
            <i>Nucleic Acids Research</i>, 43(W1), W576-W579. 
            <a href='https://doi.org/10.1093/nar/gkv402'>doi:10.1093/nar/gkv402</a>
            </div>
            
            <p><b>BLAST:</b></p>
            <div class="citation">
            Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. 
            <i>Journal of Molecular Biology</i>, 215(3), 403-410. 
            <a href='https://doi.org/10.1016/S0022-2836(05)80360-2'>doi:10.1016/S0022-2836(05)80360-2</a>
            </div>
            
            <div class="citation">
            Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J., Zhang, Z., Miller, W., & Lipman, D. J. (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. 
            <i>Nucleic Acids Research</i>, 25(17), 3389-3402. 
            <a href='https://doi.org/10.1093/nar/25.17.3389'>doi:10.1093/nar/25.17.3389</a>
            </div>
            
            <p><b>InterPro:</b></p>
            <div class="citation">
            Apweiler, R., Attwood, T. K., Bairoch, A., Bateman, A., Birney, E., Biswas, M., ... & Mulder, N. J. (2001). The InterPro database, an integrated documentation resource for protein families, domains and functional sites. 
            <i>Nucleic Acids Research</i>, 29(1), 37-40. 
            <a href='https://doi.org/10.1093/nar/29.1.37'>doi:10.1093/nar/29.1.37</a>
            </div>
            
            <div class="citation">
            Paysan-Lafosse, T., Blum, M., Chuguransky, S., Grego, T., Pinto, B. L., Salazar, G. A., ... & Bridge, A. (2023). InterPro in 2022. 
            <i>Nucleic Acids Research</i>, 51(D1), D418-D427. 
            <a href='https://doi.org/10.1093/nar/gkac993'>doi:10.1093/nar/gkac993</a>
            </div>
            
            <p><b>Pfam:</b></p>
            <div class="citation">
            Bateman, A., Birney, E., Durbin, R., Eddy, S. R., Howe, K. L., & Sonnhammer, E. L. (2000). The Pfam protein families database. 
            <i>Nucleic Acids Research</i>, 28(1), 263-266. 
            <a href='https://doi.org/10.1093/nar/28.1.263'>doi:10.1093/nar/28.1.263</a>
            </div>
            
            <div class="citation">
            Mistry, J., Chuguransky, S., Williams, L., Qureshi, M., Salazar, G. A., Sonnhammer, E. L., ... & Bateman, A. (2021). Pfam: The protein families database in 2021. 
            <i>Nucleic Acids Research</i>, 49(D1), D412-D419. 
            <a href='https://doi.org/10.1093/nar/gkaa913'>doi:10.1093/nar/gkaa913</a>
            </div>
            
            <p><b>PROSITE:</b></p>
            <div class="citation">
            Bairoch, A. (1991). PROSITE: a dictionary of sites and patterns in proteins. 
            <i>Nucleic Acids Research</i>, 19(suppl), 2241-2245. 
            <a href='https://doi.org/10.1093/nar/19.suppl.2241'>doi:10.1093/nar/19.suppl.2241</a>
            </div>
            
            <div class="citation">
            Sigrist, C. J., de Castro, E., Cerutti, L., Cuche, B. A., Hulo, N., Bridge, A., ... & Xenarios, I. (2013). New and continuing developments at PROSITE. 
            <i>Nucleic Acids Research</i>, 41(D1), D344-D347. 
            <a href='https://doi.org/10.1093/nar/gks1067'>doi:10.1093/nar/gks1067</a>
            </div>
            
            <p><b>PDB Format:</b></p>
            <div class="citation">
            Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N., Weissig, H., ... & Bourne, P. E. (2000). The Protein Data Bank. 
            <i>Nucleic Acids Research</i>, 28(1), 235-242. 
            <a href='https://doi.org/10.1093/nar/28.1.235'>doi:10.1093/nar/28.1.235</a>
            </div>
            
            <p><b>ramachandraw:</b></p>
            <div class="citation">
            Alexandre Cirilo. (2024). ramachandraw: A Ramachandran plotting tool (v1.0.1). Zenodo. 
            <a href='https://doi.org/10.5281/zenodo.10585423'>doi:10.5281/zenodo.10585423</a>
            </div>
            
            
            """
        )
        about_text.setWordWrap(True)
        about_text.setOpenExternalLinks(True)
        about_text.setTextFormat(Qt.RichText)
        about_text.setTextInteractionFlags(Qt.TextBrowserInteraction)
        
        layout.addWidget(about_text)
        
        # Set the scroll area's widget
        scroll.setWidget(content)
        
        # Create main layout with scroll area and button box
        main_layout = QVBoxLayout(self)
        main_layout.addWidget(scroll)
        
        # Add OK button at the bottom
        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(self.accept)
        main_layout.addWidget(button_box)

from PyQt5.QtCore import QSettings
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEngineSettings
from PyQt5.QtCore import QUrl, Qt
from PyQt5.QtGui import QFont, QPixmap, QIcon

from Bio.PDB import PDBList, PDBParser, PDBIO
from Bio.SeqIO.PdbIO import BiopythonParserWarning

try:
    from src.core.pdb_fetch_worker import PDBFetchManager
    OPTIMIZED_FETCH_AVAILABLE = True
except ImportError:
    from src.core.enhanced_pdb_puller_fixed import EnhancedPDBPuller
    OPTIMIZED_FETCH_AVAILABLE = False


class ServerThread(threading.Thread):
    def __init__(self, directory, start_port=8000):
        super().__init__(daemon=True)
        self.port = start_port
        self.directory = directory
        self.httpd = None

    def run(self):
        handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=self.directory)
        for _ in range(20):
            try:
                self.httpd = socketserver.TCPServer(('', self.port), handler)
                print(f"Serving on port {self.port} from directory: {self.directory}")
                self.httpd.serve_forever()
                break
            except OSError:
                print(f"Port {self.port} busy, trying next...")
                self.port += 1
        else:
            raise RuntimeError("Could not bind to any port.")

    def shutdown(self):
        if self.httpd:
            self.httpd.shutdown()
            self.httpd.server_close()
            print("Server stopped.")



class ProteinViewerApp(QMainWindow):
    def __init__(self, port):
        super().__init__()
        self.port = port
        self._undo_stack = []
        self._redo_stack = []
        self._is_restoring_state = False
        self.setWindowTitle("PicoMol - Molecular Visualization Suite")
        self.setAcceptDrops(True)
        
        # Set application icon
        self.set_application_icon()
        
        # Set responsive window size based on screen dimensions
        self.setup_responsive_window()
        
        # Initialize preferences manager
        self.preferences_manager = PreferencesManager()

        # Show welcome screen if user wants it
        settings = QSettings("PicoMolApp", "PicoMol")
        show_welcome = settings.value("show_welcome", False, type=bool)
        if show_welcome:
            self._welcome_dialog = WelcomeDialog(self)
            self._welcome_dialog.setModal(True)
            def handle_close():
                settings.setValue("show_welcome", self._welcome_dialog.should_show_next_time())
                self._welcome_dialog.deleteLater()
            self._welcome_dialog.finished.connect(handle_close)
            self._welcome_dialog.show()

        # Create directories if not exist
        project_root = os.path.dirname(os.path.abspath(__file__))
        self.pulled_structures_dir = os.path.join(project_root, "data", "pulled_structures")
        os.makedirs(self.pulled_structures_dir, exist_ok=True)
        self.ngl_assets_dir = os.path.join(project_root, "assets", "ngl_assets")
        os.makedirs(self.ngl_assets_dir, exist_ok=True)

        ngl_min_js_path = os.path.join(self.ngl_assets_dir, "ngl.min.js")
        
        # Check if file exists
        if not os.path.exists(ngl_min_js_path):
            QMessageBox.information(self, "NGL.js Missing", "NGL.js not found. Attempting to download...")
            try:
                # Try to import the setup_ngl function directly
                sys.path.append(os.path.dirname(os.path.abspath(__file__)))
                from setup_ngl import setup_ngl
                
                if not setup_ngl():
                    raise Exception("Failed to download NGL.js")
                    
                if not os.path.exists(ngl_min_js_path):
                    raise Exception("Downloaded file not found")
                    
                QMessageBox.information(self, "Success", "Successfully downloaded NGL.js")
                    
            except Exception as e:
                self.show_error_dialog(
                    "NGL.js Setup Failed",
                    "Failed to download NGL.js.",
                    suggestion="Please check your internet connection and try again.",
                    details=str(e)
                )
                sys.exit(1)

        # Copy ngl.min.js to pulled_structures_dir if it doesn't exist there
        dest_path = os.path.join(self.pulled_structures_dir, "ngl.min.js")
        if not os.path.exists(dest_path):
            try:
                shutil.copy(ngl_min_js_path, self.pulled_structures_dir)
            except Exception as e:
                print(f"Warning: Could not copy ngl.min.js to {self.pulled_structures_dir}: {e}")

        self.pdb_parser = PDBParser()
        self.pdb_list = PDBList()
        
        # Initialize PDB fetch manager (optimized or fallback)
        if OPTIMIZED_FETCH_AVAILABLE:
            self.pdb_fetch_manager = PDBFetchManager(self.pulled_structures_dir)
            self.enhanced_pdb_puller = None  # For compatibility
        else:
            self.enhanced_pdb_puller = EnhancedPDBPuller(self.pulled_structures_dir)
            self.pdb_fetch_manager = None
        
        self.init_ui()
    
    def setup_responsive_window(self):
        """Set up responsive window sizing based on screen dimensions."""
        from PyQt5.QtWidgets import QDesktopWidget
        
        # Get screen geometry
        desktop = QDesktopWidget()
        screen_rect = desktop.screenGeometry()
        screen_width = screen_rect.width()
        screen_height = screen_rect.height()
        
        # Calculate responsive window size (80% of screen, with reasonable limits)
        min_width = 900
        min_height = 700
        max_width = 1600
        max_height = 1200
        
        # Calculate preferred size as 80% of screen
        preferred_width = int(screen_width * 0.8)
        preferred_height = int(screen_height * 0.8)
        
        # Apply limits
        window_width = max(min_width, min(preferred_width, max_width))
        window_height = max(min_height, min(preferred_height, max_height))
        
        # Center the window on screen
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2
        
        # Try to restore saved window geometry first
        settings = QSettings("PicoMolApp", "PicoMol")
        saved_geometry = settings.value("window_geometry")
        
        if saved_geometry:
            # Restore saved geometry
            self.restoreGeometry(saved_geometry)
            print("Restored saved window geometry")
        else:
            # Set calculated geometry
            self.setGeometry(x, y, window_width, window_height)
            print(f"Screen: {screen_width}x{screen_height}, Window: {window_width}x{window_height} at ({x}, {y})")
        
        # Set minimum size to ensure usability
        self.setMinimumSize(min_width, min_height)
    
    def set_application_icon(self):
        """Set the application icon from the logo file."""
        try:
            # Get the path to the logo file
            project_root = os.path.dirname(os.path.abspath(__file__))
            logo_path = os.path.join(project_root, "PicoMol.png")
            
            if os.path.exists(logo_path):
                # Load the logo and set as window icon
                pixmap = QPixmap(logo_path)
                if not pixmap.isNull():
                    icon = QIcon(pixmap)
                    self.setWindowIcon(icon)
                    
                    # Also set as application icon
                    QApplication.instance().setWindowIcon(icon)
                    print(f"Successfully set application icon from {logo_path}")
                else:
                    print(f"Failed to load logo from {logo_path}")
            else:
                print(f"Logo file not found at {logo_path}")
        except Exception as e:
            print(f"Error setting application icon: {e}")

    def show_error_dialog(self, title, summary, suggestion=None, details=None):
        dlg = QDialog(self)
        dlg.setWindowTitle(title)
        dlg.setMinimumWidth(420)
        layout = QVBoxLayout(dlg)
        summary_label = QLabel(f"<b>{summary}</b>")
        summary_label.setWordWrap(True)
        layout.addWidget(summary_label)
        if suggestion:
            suggestion_label = QLabel(suggestion)
            suggestion_label.setStyleSheet("color: #0077cc;")
            suggestion_label.setWordWrap(True)
            layout.addWidget(suggestion_label)
        if details:
            details_box = QTextEdit()
            details_box.setReadOnly(True)
            details_box.setPlainText(details)
            details_box.setMaximumHeight(80)
            details_box.setVisible(False)
            toggle_btn = QPushButton("Show Details")
            def toggle():
                if details_box.isVisible():
                    details_box.setVisible(False)
                    toggle_btn.setText("Show Details")
                else:
                    details_box.setVisible(True)
                    toggle_btn.setText("Hide Details")
            toggle_btn.clicked.connect(toggle)
            layout.addWidget(toggle_btn)
            layout.addWidget(details_box)
        btn_box = QDialogButtonBox(QDialogButtonBox.Ok)
        btn_box.accepted.connect(dlg.accept)
        layout.addWidget(btn_box)
        dlg.exec_()



    def capture_state(self):
        """Capture the current user-facing state for undo/redo."""
        state = {
            'representation': self.representation_combo.currentText() if hasattr(self, 'representation_combo') else None,
            'color_scheme': self.color_combo.currentText() if hasattr(self, 'color_combo') else None,
            'spin': self.spin_checkbox.isChecked() if hasattr(self, 'spin_checkbox') else None,
    
            'background_color': self.background_color_entry.text() if hasattr(self, 'background_color_entry') else None,
            'custom_color': self.custom_color_entry.text() if hasattr(self, 'custom_color_entry') else None,
            'structure_id': getattr(self, 'current_structure_id', None),
        }
        return state

    def restore_state(self, state):
        """Restore the viewer state from a snapshot."""
        self._is_restoring_state = True
        try:
            if state.get('representation') and hasattr(self, 'representation_combo'):
                self.representation_combo.setCurrentText(state['representation'])
            if state.get('color_scheme') and hasattr(self, 'color_combo'):
                self.color_combo.setCurrentText(state['color_scheme'])
            if state.get('spin') is not None and hasattr(self, 'spin_checkbox'):
                self.spin_checkbox.setChecked(state['spin'])
                self.toggle_spin()

            if state.get('background_color') and hasattr(self, 'background_color_entry'):
                self.background_color_entry.setText(state['background_color'])
                self.update_background_color()
            if hasattr(self, 'custom_color_entry'):
                val = state.get('custom_color') or ''
                print(f"[UNDO] restore_state: setting custom_color to {val!r}")
                self._is_restoring_state = True
                try:
                    self.custom_color_entry.blockSignals(True)
                    self.custom_color_entry.setText(val)
                    self._last_custom_color = val
                finally:
                    self.custom_color_entry.blockSignals(False)
                    self._is_restoring_state = False
                if hasattr(self, 'current_structure_id') and self.current_structure_id:
                    self.web_view.page().runJavaScript(f"if (typeof setCustomColor === 'function') setCustomColor('{val}');")
                else:
                    print(f"[DEBUG] Restore state: setCustomColor called but no structure loaded")
            # Structure reload if needed
            if state.get('structure_id') and hasattr(self, 'pulled_structures_dir'):
                pdb_path = os.path.join(self.pulled_structures_dir, f"{state['structure_id']}.pdb")
                if os.path.exists(pdb_path):
                    structure = self.pdb_parser.get_structure(state['structure_id'], pdb_path)
                    self.display_structure(structure)
                    self.current_structure_id = state['structure_id']
        finally:
            self._is_restoring_state = False

    def push_undo(self):
        if self._is_restoring_state:
            return
        state = self.capture_state()
        if self._undo_stack:
            last_state = self._undo_stack[-1]
            if state == last_state:
                print("[UNDO] push_undo: No change, not pushing duplicate state.")
                return
        self._undo_stack.append(state)
        self._redo_stack.clear()
        print(f"[UNDO] push_undo: Stack size: {len(self._undo_stack)}")
        print(f"[UNDO] Stack: {[s['custom_color'] for s in self._undo_stack]}")

    def undo(self):
        if len(self._undo_stack) > 1:
            prev_state = self._undo_stack.pop()
            self._redo_stack.append(prev_state)
            print(f"[UNDO] undo: Stack size: {len(self._undo_stack)}")
            print(f"[UNDO] Stack: {[s['custom_color'] for s in self._undo_stack]}")
            state = self._undo_stack[-1]
            self.restore_state(state)

    def redo(self):
        if self._redo_stack:
            state = self._redo_stack.pop()
            self._undo_stack.append(state)
            print(f"[UNDO] redo: Stack size: {len(self._undo_stack)}")
            print(f"[UNDO] Stack: {[s['custom_color'] for s in self._undo_stack]}")
            self.restore_state(state)


    def init_ui(self):
        # Create the web view first
        self.web_view = QWebEngineView()
        self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Menu bar with File, Edit, Recent Files
        menubar = self.menuBar()
        file_menu = menubar.addMenu("File")

        save_as_action = QAction("Save Structure As...", self)
        save_as_action.setShortcut("Ctrl+Shift+S")
        save_as_action.triggered.connect(self.save_structure_as)
        file_menu.addAction(save_as_action)

        self.recent_files_menu = file_menu.addMenu("Recent Files")
        self.update_recent_files_menu()

        edit_menu = menubar.addMenu("Edit")
        undo_action = QAction("Undo", self)
        undo_action.setShortcut("Ctrl+Z")
        undo_action.triggered.connect(self.undo)
        edit_menu.addAction(undo_action)
        
        redo_action = QAction("Redo", self)
        redo_action.setShortcut("Ctrl+Y")
        redo_action.triggered.connect(self.redo)
        edit_menu.addAction(redo_action)

        # Reset View Action
        reset_view_action = QAction("Reset View to Defaults", self)
        reset_view_action.setShortcut("Ctrl+R")
        reset_view_action.setToolTip("Reset all viewer settings to their default values.")
        reset_view_action.triggered.connect(self.reset_view_to_defaults)
        edit_menu.addAction(reset_view_action)

        # Tools menu with Preferences
        tools_menu = menubar.addMenu("Tools")
        preferences_action = QAction("Preferences...", self)
        preferences_action.setShortcut("Ctrl+,")
        preferences_action.setToolTip("Open application preferences")
        preferences_action.triggered.connect(self.show_preferences_dialog)
        tools_menu.addAction(preferences_action)
        
        tools_menu.addSeparator()
        
        # PDB Information action
        pdb_info_action = QAction("Show PDB Information...", self)
        pdb_info_action.setShortcut("Ctrl+I")
        pdb_info_action.setToolTip("Show comprehensive information about the current PDB structure")
        pdb_info_action.triggered.connect(self.show_pdb_information)
        tools_menu.addAction(pdb_info_action)
        
        # Help menu with About and Welcome
        help_menu = menubar.addMenu("Help")
        
        welcome_action = QAction("Show Welcome Screen", self)
        welcome_action.setToolTip("Show the welcome screen with feature overview and quick start guide.")
        welcome_action.triggered.connect(self.show_welcome_dialog)
        help_menu.addAction(welcome_action)
        
        help_menu.addSeparator()
        
        about_action = QAction("About PicoMol", self)
        about_action.setToolTip("Show version info, credits, and license.")
        about_action.triggered.connect(self.show_about_dialog)
        help_menu.addAction(about_action)

        # Create main widget and layout
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(5, 5, 5, 5)
        main_layout.setSpacing(5)
        
        # Create main tab widget
        self.main_tabs = QTabWidget()
        main_layout.addWidget(self.main_tabs)

        # Drag-and-drop overlay label (hidden by default)
        self.drag_overlay = QLabel("\n\nDrop PDB or ENT files here to open", self)
        self.drag_overlay.setAlignment(Qt.AlignCenter)
        self.drag_overlay.setStyleSheet("""
            background: rgba(30, 144, 255, 0.7); 
            color: white; 
            font-size: 28px; 
            font-weight: bold;
            border: 3px dashed white; 
            border-radius: 24px;
            padding: 20px;
        """)
        self.drag_overlay.setVisible(False)
        self.drag_overlay.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.drag_overlay.setGeometry(0, 0, self.width(), self.height())

        # Ensure overlay resizes with window
        self.resizeEvent = self._resizeEventWithOverlay

        # Status bar
        self.statusBar().showMessage("Ready")

        # Create visualization tab
        visualization_tab = QWidget()
        visualization_layout = QVBoxLayout(visualization_tab)  # Changed to QVBoxLayout
        visualization_layout.setContentsMargins(5, 5, 5, 5)
        visualization_layout.setSpacing(5)
        
        # Add sequence display at the top
        self.sequence_display = QTextEdit()
        self.sequence_display.setReadOnly(True)
        self.sequence_display.setMaximumHeight(100)  # Use maximum instead of fixed for responsiveness
        self.sequence_display.setMinimumHeight(60)   # Set minimum height
        self.sequence_display.setPlaceholderText("Sequence will appear here when a structure is loaded...")
        self.sequence_display.setStyleSheet("""
            QTextEdit {
                font-family: monospace;
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 4px;
                padding: 8px;
                margin: 5px;
                font-size: 12px;
            }
        """)
        visualization_layout.addWidget(self.sequence_display)
        
        # Create horizontal layout for control panel and viewer
        viewer_container = QWidget()
        viewer_layout = QHBoxLayout(viewer_container)
        viewer_layout.setContentsMargins(0, 0, 0, 0)
        viewer_layout.setSpacing(5)
        
        # Control Panel
        control_panel = QWidget()
        control_layout = QVBoxLayout(control_panel)
        control_layout.setContentsMargins(5, 5, 5, 5)
        
        # PDB Input Section
        pdb_group = QGroupBox("PDB Input")
        pdb_layout = QVBoxLayout()
        
        # PDB ID Input
        pdb_id_layout = QHBoxLayout()
        pdb_id_layout.addWidget(QLabel("PDB ID:"))
        self.pdb_id_entry = QLineEdit()
        self.pdb_id_entry.setPlaceholderText("e.g., 1CRN, 4HHB")
        self.pdb_id_entry.setToolTip("Enter a valid PDB ID to fetch a protein structure from the PDB database.")
        pdb_id_layout.addWidget(self.pdb_id_entry)
        
        fetch_button = QPushButton("Fetch")
        fetch_button.setToolTip("Download and visualize a protein structure using the entered PDB ID.")
        fetch_button.clicked.connect(self.fetch_pdb_id)
        pdb_id_layout.addWidget(fetch_button)
        pdb_layout.addLayout(pdb_id_layout)
        
        # File Open Button
        open_button = QPushButton("Open Local PDB File...")
        open_button.setToolTip("Open and visualize a local .pdb or .ent file from your computer.")
        open_button.clicked.connect(self.open_local_pdb)
        pdb_layout.addWidget(open_button)
        
        # Screenshot Button
        screenshot_button = QPushButton("Save Screenshot...")
        screenshot_button.setToolTip("Save a screenshot of the current protein structure view as a PNG image.")
        screenshot_button.clicked.connect(self.save_screenshot)
        pdb_layout.addWidget(screenshot_button)
        
        # Feedback Button
        feedback_button = QPushButton("Send Feedback...")
        feedback_button.setToolTip("Report a bug or suggest an improvement.")
        feedback_button.clicked.connect(self.show_feedback_dialog)
        pdb_layout.addWidget(feedback_button)
        
        pdb_group.setLayout(pdb_layout)
        control_layout.addWidget(pdb_group)
        
        # Visualization Options
        vis_group = QGroupBox("Visualization Options")
        vis_layout = QVBoxLayout()
        
        # Representation
        rep_layout = QHBoxLayout()
        rep_layout.addWidget(QLabel("Representation:"))
        self.representation_combo = QComboBox()
        self.representation_combo.setToolTip("Select the 3D representation style for the protein structure.")
        self.representation_combo.addItems(["cartoon", "ball+stick", "spacefill", "surface", "ribbon", "licorice", "tube"])
        self.representation_combo.setCurrentText("cartoon")
        self.representation_combo.currentIndexChanged.connect(self.update_representation)
        rep_layout.addWidget(self.representation_combo)
        vis_layout.addLayout(rep_layout)
        
        # Color Scheme
        color_layout = QHBoxLayout()
        color_layout.addWidget(QLabel("Color:"))
        self.color_combo = QComboBox()
        self.color_combo.setToolTip("Choose a color scheme for the structure.")
        self.color_combo.addItems(["chainid", "residueindex", "sstruc", "resname", "element", "uniform"])
        self.color_combo.setCurrentText("chainid")
        self.color_combo.currentIndexChanged.connect(self.update_color_scheme)
        color_layout.addWidget(self.color_combo)
        
        # Uniform color picker (initially hidden)
        self.uniform_color_button = QPushButton()
        self.uniform_color_button.setFixedSize(24, 24)
        self.uniform_color_button.setStyleSheet("background-color: #FF0000; border: 1px solid #999;")
        self.uniform_color_button.setToolTip("Select color for uniform coloring")
        self.uniform_color_button.clicked.connect(self.pick_uniform_color)
        self.uniform_color_button.hide()  # Initially hidden
        color_layout.addWidget(self.uniform_color_button)
        
        vis_layout.addLayout(color_layout)
        
        # Background Color
        bg_layout = QHBoxLayout()
        bg_layout.addWidget(QLabel("Background:"))
        self.background_color_entry = QLineEdit("black")
        self.background_color_entry.setToolTip("Enter a color name or hex code (e.g., 'white', '#000000').")
        self.background_color_entry.setMaximumWidth(100)
        bg_layout.addWidget(self.background_color_entry)
        
        apply_bg_button = QPushButton("Apply")
        apply_bg_button.setToolTip("Apply background color")
        apply_bg_button.clicked.connect(self.update_background_color)
        bg_layout.addWidget(apply_bg_button)
        vis_layout.addLayout(bg_layout)
        
        # Spin Toggle
        self.spin_checkbox = QCheckBox("Auto-rotate")
        self.spin_checkbox.setToolTip("Toggle automatic rotation of the 3D view.")
        self.spin_checkbox.setChecked(False)
        self.spin_checkbox.stateChanged.connect(self.toggle_spin)
        vis_layout.addWidget(self.spin_checkbox)
        
        vis_group.setLayout(vis_layout)
        control_layout.addWidget(vis_group)
        
        # Add stretch to push everything up
        control_layout.addStretch()
        
        # Add control panel and viewer to horizontal layout
        viewer_layout.addWidget(control_panel, 0)
        viewer_layout.addWidget(self.web_view, 1)
        
        # Add the viewer container to the main layout
        visualization_layout.addWidget(viewer_container, 1)
        
        # Add visualization tab to main tabs
        self.main_tabs.addTab(visualization_tab, "3D Viewer")
        
        # Create bioinformatics tab with comprehensive tools
        bioinformatics_tab = create_bioinformatics_tab(self)
        
        # Add bioinformatics tab to main tabs
        self.main_tabs.addTab(bioinformatics_tab, "Bioinformatics")

        # Create Plasmid Viewer tab
        from src.gui.plasmid_viewer import create_plasmid_viewer_tab
        plasmid_viewer_tab = create_plasmid_viewer_tab(self)
        self.main_tabs.addTab(plasmid_viewer_tab, "Plasmid Viewer")
        
        # Create BLAST tab
        blast_tab = QWidget()
        blast_layout = QVBoxLayout(blast_tab)
        
        # Create a tab widget for the different BLAST types
        blast_type_tabs = QTabWidget()
        blast_layout.addWidget(blast_type_tabs)
        
        # Add tabs for each BLAST type with NCBI-style layouts
        blastn_tab = create_ncbi_style_blastn_tab(self)
        blast_type_tabs.addTab(blastn_tab, "blastn")
        
        blastp_tab = create_ncbi_style_blastp_tab(self)
        blast_type_tabs.addTab(blastp_tab, "blastp")
        
        # Use consistent NCBI-style BLASTX interface
        blastx_tab = create_ncbi_style_blastx_tab(self)
        blast_type_tabs.addTab(blastx_tab, "blastx")
        
        tblastn_tab = create_ncbi_style_tblastn_tab(self)
        blast_type_tabs.addTab(tblastn_tab, "tblastn")
        
        tblastx_tab = create_ncbi_style_tblastx_tab(self)
        blast_type_tabs.addTab(tblastx_tab, "tblastx")
        
        self.main_tabs.addTab(blast_tab, "BLAST")
        
        # Initialize undo stack
        self._undo_stack = []
        self._redo_stack = []
        self.push_undo()
        
        # Apply preferences on startup
        self.apply_preferences()

    # BLAST methods are now handled by blast_utils module with online functionality

    def reset_view_to_defaults(self):
        """Reset all visible viewer settings to their default values and update the viewer."""
        # Store current state for undo
        self.push_undo()
        # Set defaults
        if hasattr(self, 'representation_combo'):
            self.representation_combo.setCurrentText("cartoon")
        if hasattr(self, 'color_combo'):
            self.color_combo.setCurrentIndex(0)  # atomindex or first in list
        if hasattr(self, 'spin_checkbox'):
            self.spin_checkbox.setChecked(False)
            self.toggle_spin()
        if hasattr(self, 'background_color_entry'):
            self.background_color_entry.setText("black")
            self.update_background_color()
        if hasattr(self, 'custom_color_entry'):
            self.custom_color_entry.setText("")
            self.update_custom_color()
        self.statusBar().showMessage("Viewer settings reset to defaults.")

    class FeedbackDialog(QDialog):
        def __init__(self, parent=None):
            super().__init__(parent)
            self.setWindowTitle("Send Feedback")
            self.setMinimumWidth(400)
            layout = QVBoxLayout(self)
            label = QLabel("<b>We value your feedback!</b><br>Describe any issues or suggestions below:")
            label.setWordWrap(True)
            layout.addWidget(label)
            self.text_edit = QTextEdit()
            self.text_edit.setPlaceholderText("Describe your bug, suggestion, or general feedback here...")
            layout.addWidget(self.text_edit)

            # GitHub Issues Button
            github_button = QPushButton("Open GitHub Issues Page")
            github_button.setToolTip("Report issues or suggestions directly on GitHub.")
            github_button.clicked.connect(self.open_github_issues)
            layout.addWidget(github_button)

            btn_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            btn_box.accepted.connect(self.accept)
            btn_box.rejected.connect(self.reject)
            layout.addWidget(btn_box)

        def get_feedback(self):
            return self.text_edit.toPlainText().strip()
        def open_github_issues(self):
            webbrowser.open_new_tab("https://github.com/ultimafounding/PicoMol/issues")

    def show_feedback_dialog(self):
        import urllib.parse
        dlg = self.FeedbackDialog(self)
        if dlg.exec_() == QDialog.Accepted:
            feedback = dlg.get_feedback()
            if feedback:
                # Prepare GitHub new issue URL
                base_url = "https://github.com/ultimafounding/PicoMol/issues/new"
                title = urllib.parse.quote("Feedback from PicoMol User")
                body = urllib.parse.quote(feedback)
                url = f"{base_url}?title={title}&body={body}"
                webbrowser.open_new_tab(url)
                self.statusBar().showMessage("Opening GitHub to submit your feedback...")
            else:
                self.statusBar().showMessage("No feedback entered.")

    def show_about_dialog(self):
        dlg = AboutDialog(self)
        dlg.exec_()
    
    def show_welcome_dialog(self):
        """Show the welcome dialog."""
        dlg = WelcomeDialog(self)
        dlg.exec_()
    
    def show_preferences_dialog(self):
        """Show the preferences dialog."""
        dlg = PreferencesDialog(self)
        dlg.preferences_applied.connect(self.apply_preferences)
        dlg.exec_()
    
    def notify_structure_loaded(self, structure_id, structure_path):
        """Notify components that a new structure has been loaded.
        
        This method is called when a new protein structure is loaded and allows
        various components (like the structural analysis tab) to update accordingly.
        
        Args:
            structure_id (str): The PDB ID or identifier of the loaded structure
            structure_path (str): Path to the loaded structure file
        """
        # Store current structure information
        self.current_structure_id = structure_id
        self.current_structure_path = structure_path
        
        # Notify the structural analysis tab if it exists
        if hasattr(self, 'main_tabs'):
            # Find the bioinformatics tab
            for i in range(self.main_tabs.count()):
                if self.main_tabs.tabText(i) == "Bioinformatics":
                    bio_tab = self.main_tabs.widget(i)
                    if hasattr(bio_tab, 'findChild'):
                        # Look for the structural analysis tab within the bioinformatics tab
                        from PyQt5.QtWidgets import QTabWidget
                        bio_tabs = bio_tab.findChild(QTabWidget)
                        if bio_tabs:
                            for j in range(bio_tabs.count()):
                                if bio_tabs.tabText(j) == "Structure":
                                    structure_tab = bio_tabs.widget(j)
                                    if hasattr(structure_tab, 'on_structure_loaded'):
                                        structure_tab.on_structure_loaded(structure_id, structure_path)
                                    break
                    break
        
        # Update status bar
        self.statusBar().showMessage(f"Structure {structure_id} loaded and components notified", 2000)
    
    def apply_preferences(self):
        """Apply preferences to the current application state."""
        # Get visualization defaults and apply them
        viz_defaults = self.preferences_manager.get_visualization_defaults()
        
        # Apply default representation if no structure is loaded
        if hasattr(self, 'representation_combo'):
            current_rep = self.representation_combo.currentText()
            default_rep = viz_defaults['representation']
            if current_rep == 'cartoon':  # Only change if still at default
                self.representation_combo.setCurrentText(default_rep)
        
        # Apply default color scheme
        if hasattr(self, 'color_combo'):
            current_color = self.color_combo.currentText()
            default_color = viz_defaults['color_scheme']
            if current_color == 'chainid':  # Only change if still at default
                self.color_combo.setCurrentText(default_color)
        
        # Apply default background color
        if hasattr(self, 'background_color_entry'):
            default_bg = viz_defaults['background_color']
            self.background_color_entry.setText(default_bg)
            self.update_background_color()
        
        # Apply auto-spin setting
        if hasattr(self, 'spin_checkbox'):
            default_spin = viz_defaults['auto_spin']
            self.spin_checkbox.setChecked(default_spin)
            self.toggle_spin()
        
        # Apply interface settings
        interface_settings = self.preferences_manager.get_interface_settings()
        
        # Apply theme
        theme = interface_settings['theme']
        apply_theme(theme)
        
        # Apply font if specified
        font_family = interface_settings['font_family']
        font_size = interface_settings['font_size']
        if font_family and font_size:
            font = QFont(font_family, font_size)
            QApplication.instance().setFont(font)
        
        # Update tooltips visibility
        show_tooltips = interface_settings['show_tooltips']
        if not show_tooltips:
            # Disable tooltips for all widgets
            for widget in self.findChildren(QWidget):
                widget.setToolTip("")
        
        self.statusBar().showMessage("Preferences applied successfully.", 3000)

    def update_custom_color(self):
        if hasattr(self, 'custom_color_entry'):
            color = self.custom_color_entry.text()
            print(f"[UNDO] update_custom_color: Applying color {color}")
            # Only call JavaScript if a structure is loaded
            if hasattr(self, 'current_structure_id') and self.current_structure_id:
                self.web_view.page().runJavaScript(f"if (typeof setCustomColor === 'function') setCustomColor('{color}');")
            if not self._is_restoring_state:
                prev_color = None
                if self._undo_stack:
                    prev_color = self._undo_stack[-1].get('custom_color', None)
                if color != prev_color:
                    self.push_undo()
                    self._last_custom_color = color

    def update_background_color(self):
        color = self.background_color_entry.text()
        print(f"[UNDO] update_background_color: Applying color {color}")
        # Only call JavaScript if a structure is loaded
        if hasattr(self, 'current_structure_id') and self.current_structure_id:
            self.web_view.page().runJavaScript(f"if (typeof setBackgroundColor === 'function') setBackgroundColor('{color}');")
        if not self._is_restoring_state:
            prev_color = None
        if self._undo_stack:
            prev_color = self._undo_stack[-1].get('background_color', None)
        if color != prev_color:
            self.push_undo()

    def update_representation(self):
        self.push_undo()
        representation = self.representation_combo.currentText()
        # Only call JavaScript if a structure is loaded
        if hasattr(self, 'current_structure_id') and self.current_structure_id:
            self.web_view.page().runJavaScript(f"if (typeof setRepresentation === 'function') setRepresentation('{representation}');")
        else:
            print(f"[DEBUG] Representation changed to {representation} but no structure loaded")

    def clear_all_representations(self):
        # Only call JavaScript if a structure is loaded
        if hasattr(self, 'current_structure_id') and self.current_structure_id:
            self.web_view.page().runJavaScript("if (typeof clearAllRepresentations === 'function') clearAllRepresentations();")
        else:
            print("[DEBUG] Clear representations called but no structure loaded")

    def pick_uniform_color(self):
        """Open a color dialog to pick a uniform color"""
        color = QColorDialog.getColor()
        if color.isValid():
            # Convert QColor to hex string
            color_hex = color.name()
            # Update button color
            self.uniform_color_button.setStyleSheet(f"background-color: {color_hex}; border: 1px solid #999;")
            # Apply the color to the structure only if loaded
            if hasattr(self, 'current_structure_id') and self.current_structure_id:
                self.web_view.page().runJavaScript(f"if (typeof setUniformColor === 'function') setUniformColor('{color_hex}');")
            else:
                print(f"[DEBUG] Uniform color changed to {color_hex} but no structure loaded")
    
    def update_color_scheme(self):
        self.push_undo()
        color_scheme = self.color_combo.currentText()
        
        # Show/hide color picker based on selection
        self.uniform_color_button.setVisible(color_scheme == "uniform")
        
        # Only call JavaScript if a structure is loaded
        if hasattr(self, 'current_structure_id') and self.current_structure_id:
            self.web_view.page().runJavaScript(f"if (typeof setColorScheme === 'function') setColorScheme('{color_scheme}');")
            
            # If switching to uniform, apply the current button color
            if color_scheme == "uniform" and self.uniform_color_button.styleSheet():
                # Extract current color from button style
                style = self.uniform_color_button.styleSheet()
                if 'background-color:' in style:
                    color = style.split('background-color:')[1].split(';')[0].strip()
                    self.web_view.page().runJavaScript(f"if (typeof setUniformColor === 'function') setUniformColor('{color}');")
        else:
            print(f"[DEBUG] Color scheme changed to {color_scheme} but no structure loaded")

    def toggle_spin(self):
        self.push_undo()
        spin_enabled = self.spin_checkbox.isChecked()
        # Only call JavaScript if a structure is loaded
        if hasattr(self, 'current_structure_id') and self.current_structure_id:
            self.web_view.page().runJavaScript(f"if (typeof setSpin === 'function') setSpin({str(spin_enabled).lower()});")
        else:
            print(f"[DEBUG] Spin changed to {spin_enabled} but no structure loaded")

    def fetch_pdb_id(self):
        pdb_id = self.pdb_id_entry.text().strip().upper()
        if not pdb_id:
            self.statusBar().showMessage("Please enter a PDB ID.")
            return
        
        # Use optimized fetch manager if available
        if OPTIMIZED_FETCH_AVAILABLE and self.pdb_fetch_manager:
            self._fetch_pdb_optimized(pdb_id)
        else:
            self._fetch_pdb_fallback(pdb_id)
    
    def _fetch_pdb_optimized(self, pdb_id: str):
        """Fetch PDB data using the optimized, non-blocking approach."""
        # Check if already fetching
        if self.pdb_fetch_manager.is_fetching():
            self.statusBar().showMessage("Another fetch operation is in progress...")
            return
        
        # Disable fetch button during operation
        fetch_button = self.sender() if hasattr(self, 'sender') else None
        if fetch_button:
            fetch_button.setEnabled(False)
        
        # Define callbacks
        def on_progress(message: str):
            self.statusBar().showMessage(message)
        
        def on_success(comprehensive_data: dict):
            try:
                # Get the PDB file path
                pdb_path = comprehensive_data['files'].get('pdb')
                if not pdb_path or not os.path.exists(pdb_path):
                    raise FileNotFoundError(f"PDB file not found for {pdb_id}")

                # Load and display the structure
                structure = self.pdb_parser.get_structure(pdb_id, pdb_path)
                self.current_structure_id = pdb_id
                self.display_structure(structure)
                
                # Store comprehensive data for later use
                self.current_pdb_data = comprehensive_data
                
                # Update status with comprehensive info
                metadata = comprehensive_data.get('metadata', {})
                if 'entry' in metadata:
                    entry = metadata['entry']
                    
                    # Extract title with fallback
                    title = 'Unknown'
                    if 'struct' in entry and 'title' in entry['struct']:
                        title = entry['struct']['title']
                    
                    # Extract method with fallback
                    method = 'Unknown'
                    if 'exptl' in entry and isinstance(entry['exptl'], list) and entry['exptl']:
                        method = entry['exptl'][0].get('method', 'Unknown')
                    elif 'rcsb_entry_info' in entry and 'experimental_method' in entry['rcsb_entry_info']:
                        method = entry['rcsb_entry_info']['experimental_method']
                    
                    self.statusBar().showMessage(f"Loaded {pdb_id}: {title[:50]}... ({method})")
                else:
                    self.statusBar().showMessage(f"Loaded {pdb_id} with comprehensive data")
                
                self.add_to_recent_files(pdb_path)
                
                # Notify structural analysis tab that a new structure is loaded
                self.notify_structure_loaded(pdb_id, pdb_path)
                
                # Show summary of fetched data (non-blocking)
                from PyQt5.QtCore import QTimer
                QTimer.singleShot(100, lambda: self.show_fetch_summary(comprehensive_data))
                
            except Exception as e:
                self.statusBar().showMessage(f"Error loading structure {pdb_id}")
                self.show_error_dialog(
                    "Error Loading Structure",
                    f"Could not load structure for PDB ID {pdb_id}.",
                    suggestion="The file may be corrupted or invalid.",
                    details=str(e)
                )
        
        def on_error(error_message: str):
            self.statusBar().showMessage(f"Error fetching PDB ID {pdb_id}")
            self.show_error_dialog(
                "Error Fetching PDB",
                f"Could not fetch data for PDB ID {pdb_id}.",
                suggestion="Please check your internet connection and verify the PDB ID is correct.",
                details=error_message
            )
        
        def on_finished():
            # Re-enable fetch button
            if fetch_button:
                fetch_button.setEnabled(True)
        
        # Start the async fetch
        self.pdb_fetch_manager.fetch_pdb_async(
            pdb_id,
            progress_callback=on_progress,
            success_callback=on_success,
            error_callback=on_error,
            finished_callback=on_finished,
            include_validation=False,  # Skip validation for speed
            include_sequences=True,
            include_mmcif=False
        )
    
    def _fetch_pdb_fallback(self, pdb_id: str):
        """Fallback PDB fetch using the original enhanced puller."""
        try:
            self.statusBar().showMessage(f"Fetching comprehensive data for PDB ID {pdb_id}...")
            
            # Use enhanced PDB puller to get comprehensive data
            comprehensive_data = self.enhanced_pdb_puller.fetch_comprehensive_pdb_data(
                pdb_id, 
                include_validation=False,  # Skip validation for speed
                include_sequences=True,
                include_mmcif=False
            )
            
            if comprehensive_data['errors']:
                error_details = "\n".join(comprehensive_data['errors'])
                raise Exception(f"Errors occurred during data fetching:\n{error_details}")
            
            # Get the PDB file path
            pdb_path = comprehensive_data['files'].get('pdb')
            if not pdb_path or not os.path.exists(pdb_path):
                raise FileNotFoundError(f"PDB file not found for {pdb_id}")

            # Load and display the structure
            structure = self.pdb_parser.get_structure(pdb_id, pdb_path)
            self.current_structure_id = pdb_id
            self.display_structure(structure)
            
            # Store comprehensive data for later use
            self.current_pdb_data = comprehensive_data
            
            # Update status with comprehensive info
            metadata = comprehensive_data.get('metadata', {})
            if 'entry' in metadata:
                entry = metadata['entry']
                
                # Extract title with fallback
                title = 'Unknown'
                if 'struct' in entry and 'title' in entry['struct']:
                    title = entry['struct']['title']
                
                # Extract method with fallback
                method = 'Unknown'
                if 'exptl' in entry and isinstance(entry['exptl'], list) and entry['exptl']:
                    method = entry['exptl'][0].get('method', 'Unknown')
                elif 'rcsb_entry_info' in entry and 'experimental_method' in entry['rcsb_entry_info']:
                    method = entry['rcsb_entry_info']['experimental_method']
                
                self.statusBar().showMessage(f"Loaded {pdb_id}: {title[:50]}... ({method})")
            else:
                self.statusBar().showMessage(f"Loaded {pdb_id} with comprehensive data")
            
            self.add_to_recent_files(pdb_path)
            
            # Notify structural analysis tab that a new structure is loaded
            self.notify_structure_loaded(pdb_id, pdb_path)
            
            # Show summary of fetched data
            self.show_fetch_summary(comprehensive_data)

        except Exception as e:
            self.statusBar().showMessage(f"Error fetching PDB ID {pdb_id}")
            self.show_error_dialog(
                "Error Fetching PDB",
                f"Could not fetch comprehensive data for PDB ID {pdb_id}.",
                suggestion="Please check your internet connection and verify the PDB ID is correct.",
                details=str(e)
            )

    def show_fetch_summary(self, comprehensive_data):
        """Show a summary of the fetched comprehensive data with advanced metadata details."""
        pdb_id = comprehensive_data['pdb_id']
        files = comprehensive_data['files']
        metadata = comprehensive_data['metadata']
        
        summary_parts = [f"Successfully fetched comprehensive data for {pdb_id}:"]
        
        # Files downloaded
        if files:
            summary_parts.append(f"\n📁 Files downloaded: {', '.join(files.keys())}")
        
        # Advanced metadata sections
        if metadata:
            summary_parts.append(f"\n🔬 Advanced metadata sections: {', '.join(metadata.keys())}")
        
        # Detailed info if available
        if 'entry' in metadata:
            entry = metadata['entry']
            
            # Basic structure info
            title = entry.get('struct', {}).get('title', 'N/A')
            rcsb_info = entry.get('rcsb_entry_info', {})
            method = rcsb_info.get('experimental_method', 'N/A')
            resolution = rcsb_info.get('resolution_combined', 'N/A')
            if isinstance(resolution, list) and resolution:
                resolution = resolution[0]
            
            summary_parts.extend([
                f"\n📋 Title: {title}",
                f"🧪 Method: {method}",
                f"🔍 Resolution: {resolution} Å" if resolution != 'N/A' else "🔍 Resolution: N/A"
            ])
            
            # Advanced experimental details
            mol_weight = rcsb_info.get('molecular_weight')
            atom_count = rcsb_info.get('deposited_atom_count')
            if mol_weight is not None and isinstance(mol_weight, (int, float)):
                summary_parts.append(f"⚖️ Molecular Weight: {mol_weight:,.0f} Da")
            if atom_count is not None and isinstance(atom_count, (int, float)):
                summary_parts.append(f"⚛️ Total Atoms: {atom_count:,}")
            
            # Refinement statistics
            if 'refine' in entry and entry['refine']:
                refine = entry['refine'][0] if isinstance(entry['refine'], list) else entry['refine']
                r_work = refine.get('ls_R_factor_R_work')
                r_free = refine.get('ls_R_factor_R_free')
                # Check for both None and valid numeric values
                if r_work is not None and r_free is not None and isinstance(r_work, (int, float)) and isinstance(r_free, (int, float)):
                    summary_parts.append(f"📊 R-factors: R-work={r_work:.3f}, R-free={r_free:.3f}")
                elif r_work is not None and isinstance(r_work, (int, float)):
                    summary_parts.append(f"📊 R-work: {r_work:.3f}")
                elif r_free is not None and isinstance(r_free, (int, float)):
                    summary_parts.append(f"📊 R-free: {r_free:.3f}")
            
            # Crystal parameters
            if 'cell' in entry and entry['cell']:
                cell = entry['cell']
                volume = cell.get('volume')
                if volume is not None and isinstance(volume, (int, float)):
                    summary_parts.append(f"💎 Unit Cell Volume: {volume:,.0f} Å³")
            
            # Space group
            if 'symmetry' in entry and entry['symmetry']:
                space_group = entry['symmetry'].get('space_group_name_H_M', 'N/A')
                if space_group != 'N/A':
                    summary_parts.append(f"🔷 Space Group: {space_group}")
        
        # Entity information
        if 'polymer_entities' in metadata and metadata['polymer_entities']:
            entity_count = len(metadata['polymer_entities'])
            summary_parts.append(f"🧬 Polymer Entities: {entity_count} detailed entities")
        
        if 'nonpolymer_entities' in metadata and metadata['nonpolymer_entities']:
            ligand_count = len(metadata['nonpolymer_entities'])
            summary_parts.append(f"💊 Ligands/Small Molecules: {ligand_count} entities")
        
        # Assembly information
        if 'assemblies' in metadata and metadata['assemblies']:
            assembly_count = len(metadata['assemblies'])
            summary_parts.append(f"🏗️ Biological Assemblies: {assembly_count} assemblies")
        
        # Sequences info
        sequences = comprehensive_data.get('sequences', {})
        if sequences.get('chains'):
            chain_count = len(sequences['chains'])
            summary_parts.append(f"🔗 Sequences: {chain_count} protein chains")
        
        # Publication info
        if 'entry' in metadata and 'rcsb_primary_citation' in metadata['entry']:
            citation = metadata['entry']['rcsb_primary_citation']
            if citation and citation.get('title'):
                summary_parts.append(f"📚 Primary Citation: Available")
        
        summary_text = "\n".join(summary_parts)
        summary_text += "\n\n✨ This comprehensive dataset includes experimental details, refinement statistics, crystal parameters, entity information, and much more!"
        
        # Show in a message box
        from PyQt5.QtWidgets import QMessageBox, QTextEdit
        msg = QMessageBox(self)
        msg.setWindowTitle("🎉 Advanced Metadata Fetched Successfully")
        msg.setText(f"Successfully fetched comprehensive PDB data for {pdb_id}")
        msg.setDetailedText(summary_text)
        msg.setIcon(QMessageBox.Information)
        msg.exec_()
    
    def show_pdb_information(self):
        """Show comprehensive PDB information dialog."""
        if not hasattr(self, 'current_structure_id') or not self.current_structure_id:
            QMessageBox.warning(self, "No Structure", "No PDB structure is currently loaded.")
            return
        
        pdb_id = self.current_structure_id
        
        # Get comprehensive data if available
        if hasattr(self, 'current_pdb_data') and self.current_pdb_data:
            self.show_comprehensive_pdb_dialog(self.current_pdb_data)
        else:
            # Try to get cached data from appropriate puller
            if OPTIMIZED_FETCH_AVAILABLE and self.pdb_fetch_manager:
                info = self.pdb_fetch_manager.get_structure_info(pdb_id)
            elif self.enhanced_pdb_puller:
                info = self.enhanced_pdb_puller.get_structure_info(pdb_id)
            else:
                QMessageBox.warning(self, "No Puller Available", "No PDB puller is available.")
                return
            
            if info['available']:
                # Create a dialog with available information
                self.show_basic_pdb_info_dialog(info)
            else:
                QMessageBox.information(self, "No Enhanced Data", 
                    f"No comprehensive data available for {pdb_id}.\n\n"
                    f"Use 'Fetch PDB' to download comprehensive information.")
    
    def show_comprehensive_pdb_dialog(self, comprehensive_data):
        """Show comprehensive PDB information in a dialog."""
        from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTabWidget, QTextBrowser, QDialogButtonBox
        
        dialog = QDialog(self)
        dialog.setWindowTitle(f"PDB Information: {comprehensive_data['pdb_id']}")
        dialog.setMinimumSize(800, 600)
        
        layout = QVBoxLayout(dialog)
        
        # Create tab widget
        tabs = QTabWidget()
        
        # Summary tab
        summary_browser = QTextBrowser()
        summary_text = self.format_pdb_summary(comprehensive_data)
        summary_browser.setHtml(summary_text)
        tabs.addTab(summary_browser, "Summary")
        
        # Metadata tab
        metadata_browser = QTextBrowser()
        metadata_text = self.format_metadata_display(comprehensive_data['metadata'])
        metadata_browser.setHtml(metadata_text)
        tabs.addTab(metadata_browser, "Detailed Metadata")
        
        # Sequences tab
        if comprehensive_data.get('sequences'):
            seq_browser = QTextBrowser()
            seq_text = self.format_sequences_display(comprehensive_data['sequences'])
            seq_browser.setHtml(seq_text)
            tabs.addTab(seq_browser, "Sequences")
        
        # Files tab
        files_browser = QTextBrowser()
        files_text = self.format_files_display(comprehensive_data['files'])
        files_browser.setHtml(files_text)
        tabs.addTab(files_browser, "Downloaded Files")
        
        layout.addWidget(tabs)
        
        # Button box
        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(dialog.accept)
        layout.addWidget(button_box)
        
        dialog.exec_()
    
    def format_pdb_summary(self, comprehensive_data):
        """Format PDB summary for display with comprehensive metadata."""
        pdb_id = comprehensive_data['pdb_id']
        metadata = comprehensive_data.get('metadata', {})
        
        html = f"<h2>PDB Structure: {pdb_id}</h2>"
        
        if 'entry' in metadata:
            entry = metadata['entry']
            
            # Title
            title = entry.get('struct', {}).get('title', 'N/A')
            html += f"<h3>{title}</h3>"
            
            # Basic information
            html += "<h4>Basic Information</h4><ul>"
            
            # Method and resolution
            rcsb_info = entry.get('rcsb_entry_info', {})
            method = rcsb_info.get('experimental_method', 'N/A')
            html += f"<li><b>Experimental Method:</b> {method}</li>"
            
            resolution = rcsb_info.get('resolution_combined', 'N/A')
            if resolution != 'N/A':
                if isinstance(resolution, list) and resolution:
                    resolution = resolution[0]
                html += f"<li><b>Resolution:</b> {resolution} Å</li>"
            
            # Molecular weight and composition
            mol_weight = rcsb_info.get('molecular_weight')
            if mol_weight is not None and isinstance(mol_weight, (int, float)):
                html += f"<li><b>Molecular Weight:</b> {mol_weight:,.0f} Da</li>"
            
            # Atom and residue counts
            atom_count = rcsb_info.get('deposited_atom_count')
            residue_count = rcsb_info.get('deposited_residue_count')
            if atom_count is not None and isinstance(atom_count, (int, float)):
                html += f"<li><b>Total Atoms:</b> {atom_count:,}</li>"
            if residue_count is not None and isinstance(residue_count, (int, float)):
                html += f"<li><b>Total Residues:</b> {residue_count:,}</li>"
            
            # Entity counts
            polymer_count = rcsb_info.get('deposited_polymer_entity_instance_count', 0)
            nonpolymer_count = rcsb_info.get('deposited_nonpolymer_entity_instance_count', 0)
            if polymer_count > 0:
                html += f"<li><b>Polymer Entities:</b> {polymer_count}</li>"
            if nonpolymer_count > 0:
                html += f"<li><b>Non-polymer Entities (Ligands):</b> {nonpolymer_count}</li>"
            
            # Dates
            accession_info = entry.get('rcsb_accession_info', {})
            deposit_date = accession_info.get('deposit_date', 'N/A')
            release_date = accession_info.get('initial_release_date', 'N/A')
            revision_date = accession_info.get('revision_date', 'N/A')
            html += f"<li><b>Deposition Date:</b> {deposit_date}</li>"
            html += f"<li><b>Release Date:</b> {release_date}</li>"
            if revision_date != 'N/A':
                html += f"<li><b>Last Revision:</b> {revision_date}</li>"
            
            html += "</ul>"
            
            # Experimental details
            if 'refine' in entry and entry['refine']:
                refine = entry['refine'][0] if isinstance(entry['refine'], list) else entry['refine']
                html += "<h4>Refinement Statistics</h4><ul>"
                
                r_work = refine.get('ls_R_factor_R_work')
                r_free = refine.get('ls_R_factor_R_free')
                if r_work is not None and isinstance(r_work, (int, float)):
                    html += f"<li><b>R-work:</b> {r_work:.3f}</li>"
                if r_free is not None and isinstance(r_free, (int, float)):
                    html += f"<li><b>R-free:</b> {r_free:.3f}</li>"
                
                res_high = refine.get('ls_d_res_high')
                res_low = refine.get('ls_d_res_low')
                if res_high is not None and isinstance(res_high, (int, float)):
                    html += f"<li><b>High Resolution Limit:</b> {res_high} Å</li>"
                if res_low is not None and isinstance(res_low, (int, float)):
                    html += f"<li><b>Low Resolution Limit:</b> {res_low} Å</li>"
                
                html += "</ul>"
            
            # Crystal information
            if 'cell' in entry and entry['cell']:
                cell = entry['cell']
                html += "<h4>Crystal Parameters</h4><ul>"
                html += f"<li><b>Unit Cell:</b> a={cell.get('length_a', 'N/A')} Å, "
                html += f"b={cell.get('length_b', 'N/A')} Å, c={cell.get('length_c', 'N/A')} Å</li>"
                html += f"<li><b>Angles:</b> α={cell.get('angle_alpha', 'N/A')}°, "
                html += f"β={cell.get('angle_beta', 'N/A')}°, γ={cell.get('angle_gamma', 'N/A')}°</li>"
                
                volume = cell.get('volume')
                if volume is not None and isinstance(volume, (int, float)):
                    html += f"<li><b>Volume:</b> {volume:,.0f} Å³</li>"
                html += "</ul>"
            
            if 'symmetry' in entry and entry['symmetry']:
                symmetry = entry['symmetry']
                space_group = symmetry.get('space_group_name_H_M', 'N/A')
                if space_group != 'N/A':
                    html += f"<h4>Space Group</h4><p>{space_group}</p>"
            
            # Authors
            authors = [author.get('name', '') for author in entry.get('audit_author', [])]
            if authors:
                html += "<h4>Authors</h4><p>"
                html += ', '.join(authors[:5])
                if len(authors) > 5:
                    html += f" and {len(authors) - 5} others"
                html += "</p>"
            
            # Publication information
            if 'rcsb_primary_citation' in entry and entry['rcsb_primary_citation']:
                citation = entry['rcsb_primary_citation']
                html += "<h4>Primary Citation</h4>"
                title = citation.get('title', 'N/A')
                journal = citation.get('journal_abbrev', 'N/A')
                year = citation.get('year', 'N/A')
                doi = citation.get('pdbx_database_id_DOI', 'N/A')
                
                html += f"<p><b>Title:</b> {title}</p>"
                html += f"<p><b>Journal:</b> {journal} ({year})</p>"
                if doi != 'N/A':
                    html += f"<p><b>DOI:</b> <a href='https://doi.org/{doi}' target='_blank'>{doi}</a></p>"
        
        # Enhanced entity information
        if 'polymer_entities' in metadata and metadata['polymer_entities']:
            entities = metadata['polymer_entities']
            html += f"<h4>Polymer Entities ({len(entities)})</h4><ul>"
            for entity in entities[:3]:  # Show first 3
                entity_info = entity.get('rcsb_polymer_entity', {})
                description = entity_info.get('pdbx_description', 'Unknown')
                entity_type = entity_info.get('type', 'Unknown')
                html += f"<li><b>{entity_type}:</b> {description}</li>"
            if len(entities) > 3:
                html += f"<li><i>... and {len(entities) - 3} more entities</i></li>"
            html += "</ul>"
        
        if 'nonpolymer_entities' in metadata and metadata['nonpolymer_entities']:
            ligands = metadata['nonpolymer_entities']
            html += f"<h4>Ligands and Small Molecules ({len(ligands)})</h4><ul>"
            for ligand in ligands[:5]:  # Show first 5
                ligand_info = ligand.get('rcsb_nonpolymer_entity', {})
                description = ligand_info.get('pdbx_description', 'Unknown')
                html += f"<li>{description}</li>"
            if len(ligands) > 5:
                html += f"<li><i>... and {len(ligands) - 5} more ligands</i></li>"
            html += "</ul>"
        
        # Files information
        files = comprehensive_data.get('files', {})
        if files:
            html += "<h4>Downloaded Files</h4><ul>"
            for file_type, file_path in files.items():
                html += f"<li><b>{file_type.upper()}:</b> {os.path.basename(file_path)}</li>"
            html += "</ul>"
        
        # Sequences information
        sequences = comprehensive_data.get('sequences', {})
        if sequences.get('chains'):
            html += f"<h4>Sequences</h4><p>{len(sequences['chains'])} protein chains available</p>"
        
        return html
    
    def format_metadata_display(self, metadata):
        """Format metadata for detailed display."""
        html = "<h3>Detailed Metadata</h3>"
        
        for section, data in metadata.items():
            html += f"<h4>{section.replace('_', ' ').title()}</h4>"
            html += "<pre style='background-color: #f5f5f5; padding: 10px; font-size: 12px;'>"
            
            # Format the data nicely
            if isinstance(data, dict):
                html += self.format_dict_for_display(data, indent=0)
            elif isinstance(data, list):
                for i, item in enumerate(data[:3]):  # Show first 3 items
                    html += f"Item {i+1}:\n"
                    if isinstance(item, dict):
                        html += self.format_dict_for_display(item, indent=2)
                    else:
                        html += f"  {item}\n"
                if len(data) > 3:
                    html += f"... and {len(data) - 3} more items\n"
            else:
                html += str(data)
            
            html += "</pre>"
        
        return html
    
    def format_dict_for_display(self, data, indent=0):
        """Format dictionary data for display."""
        result = ""
        prefix = "  " * indent
        
        for key, value in data.items():
            if isinstance(value, dict):
                result += f"{prefix}{key}:\n"
                result += self.format_dict_for_display(value, indent + 1)
            elif isinstance(value, list):
                result += f"{prefix}{key}: [{len(value)} items]\n"
                if value and len(value) <= 3:
                    for item in value:
                        result += f"{prefix}  - {item}\n"
            else:
                result += f"{prefix}{key}: {value}\n"
        
        return result
    
    def format_sequences_display(self, sequences):
        """Format sequences for display."""
        html = "<h3>Protein Sequences</h3>"
        
        chains = sequences.get('chains', [])
        for chain in chains:
            html += f"<h4>Chain: {chain.get('chain_id', 'Unknown')}</h4>"
            html += f"<p><b>Description:</b> {chain.get('description', 'N/A')}</p>"
            html += f"<p><b>Length:</b> {len(chain.get('sequence', ''))} residues</p>"
            
            sequence = chain.get('sequence', '')
            if sequence:
                html += "<p><b>Sequence:</b></p>"
                html += "<pre style='background-color: #f5f5f5; padding: 10px; font-family: monospace; font-size: 12px; word-wrap: break-word;'>"
                # Format sequence with line breaks every 80 characters
                for i in range(0, len(sequence), 80):
                    html += sequence[i:i+80] + "\n"
                html += "</pre>"
        
        return html
    
    def format_files_display(self, files):
        """Format files information for display."""
        html = "<h3>Downloaded Files</h3>"
        
        for file_type, file_path in files.items():
            html += f"<h4>{file_type.upper()} File</h4>"
            html += f"<p><b>Path:</b> {file_path}</p>"
            
            if os.path.exists(file_path):
                file_size = os.path.getsize(file_path)
                html += f"<p><b>Size:</b> {file_size:,} bytes</p>"
                
                # Show modification time
                import time
                mod_time = os.path.getmtime(file_path)
                mod_time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(mod_time))
                html += f"<p><b>Downloaded:</b> {mod_time_str}</p>"
            else:
                html += "<p><i>File not found</i></p>"
        
        return html
    
    def show_basic_pdb_info_dialog(self, info):
        """Show basic PDB information dialog."""
        from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTextBrowser, QDialogButtonBox
        
        dialog = QDialog(self)
        dialog.setWindowTitle(f"PDB Information: {info['pdb_id']}")
        dialog.setMinimumSize(600, 400)
        
        layout = QVBoxLayout(dialog)
        
        browser = QTextBrowser()
        
        html = f"<h2>PDB Structure: {info['pdb_id']}</h2>"
        
        if info.get('summary'):
            summary = info['summary']
            html += f"<h3>{summary.get('title', 'Unknown Title')}</h3>"
            html += f"<p><b>Method:</b> {summary.get('experimental_method', 'N/A')}</p>"
            html += f"<p><b>Resolution:</b> {summary.get('resolution', 'N/A')}</p>"
            html += f"<p><b>Fetched:</b> {summary.get('fetch_timestamp', 'N/A')}</p>"
        
        if info.get('files'):
            html += "<h4>Available Files</h4><ul>"
            for file_type, file_path in info['files'].items():
                html += f"<li><b>{file_type.upper()}:</b> {os.path.basename(file_path)}</li>"
            html += "</ul>"
        
        html += "<p><i>For comprehensive information, use 'Fetch PDB' to download enhanced data.</i></p>"
        
        browser.setHtml(html)
        layout.addWidget(browser)
        
        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(dialog.accept)
        layout.addWidget(button_box)
        
        dialog.exec_()
    
    def open_local_pdb(self, file_path=None):
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Open Local PDB File", "", "PDB Files (*.pdb *.ent);;All Files (*)"
            )
        if not file_path:
            return
        try:
            self.statusBar().showMessage(f"Opening local PDB file: {os.path.basename(file_path)}...")
            structure_id = os.path.basename(file_path).split(".")[0]
            pdb_filename = f"{structure_id}.pdb"
            target_path = os.path.join(self.pulled_structures_dir, pdb_filename)
            
            # Only copy if the file is not already in the target location
            if os.path.abspath(file_path) != os.path.abspath(target_path):
                import shutil
                shutil.copy(file_path, target_path)
                
            self.statusBar().showMessage(f"Loading structure from {os.path.basename(file_path)}...")
            structure = self.pdb_parser.get_structure(structure_id, target_path)
            self.current_structure_id = structure_id  # Set the current structure ID
            
            # Clear comprehensive data since we're loading a local file without metadata
            self.current_pdb_data = None
            
            self.display_structure(structure)
            
            # Notify structural analysis tab that a new structure is loaded
            self.notify_structure_loaded(structure_id, target_path)
            
            self.statusBar().showMessage(f"Displayed {structure_id}")
            self.add_to_recent_files(target_path)
            
        except Exception as e:
            self.statusBar().showMessage(f"Error opening file: {os.path.basename(file_path)}")
            self.show_error_dialog(
                "Error Opening File",
                f"Could not open the selected file: {os.path.basename(file_path)}",
                suggestion="Make sure the file exists and is a valid PDB or ENT file.",
                details=str(e)
            )

    def display_structure(self, structure):
        from pathlib import Path
        io = PDBIO()
        io.set_structure(structure)

        # Save structure PDB file
        pdb_filename = f"{structure.id}.pdb"
        pdb_path = os.path.join(self.pulled_structures_dir, pdb_filename)
        try:
            io.save(pdb_path)
        except Exception as e:
            self.show_error_dialog(
                "Error Writing PDB",
                f"Could not write PDB file for structure {structure.id}.",
                suggestion="Check directory permissions or disk space.",
                details=str(e)
            )
            print(f"Error writing PDB file: {e}")
            return

        print(f"Structure ID for sequence extraction: {structure.id}")

        # Extract and display sequence using Bio.PDB structure iteration
        try:
            from Bio.PDB.Polypeptide import protein_letters_3to1, is_aa
            
            sequences = []
            
            for model in structure:
                for chain in model:
                    chain_sequence = ""
                    residue_count = 0
                    
                    # Build sequence from amino acid residues
                    for residue in chain:
                        if is_aa(residue):
                            try:
                                res_name = residue.get_resname()
                                chain_sequence += protein_letters_3to1.get(res_name, 'X')
                                residue_count += 1
                            except:
                                chain_sequence += 'X'  # Unknown amino acid
                                residue_count += 1
                    
                    # Only add chains that have amino acid residues
                    if chain_sequence:
                        chain_id = chain.id if chain.id.strip() else 'A'  # Default to 'A' if empty
                        display_id = f"{structure.id}:Chain_{chain_id}"
                        sequence_header = f">{display_id} | {residue_count} residues"
                        
                        # Format sequence with line breaks every 80 characters for readability
                        formatted_sequence = ""
                        for i in range(0, len(chain_sequence), 80):
                            formatted_sequence += chain_sequence[i:i+80] + "\n"
                        
                        sequences.append(f"{sequence_header}\n{formatted_sequence.rstrip()}")
            
            if sequences:
                self.sequence_display.setText("\n\n".join(sequences))
            else:
                self.sequence_display.setText("No amino acid sequences found in this structure.")
                
        except Exception as e:
            self.sequence_display.setText(f"Error extracting sequence: {e}")
            print(f"Error extracting sequence: {e}")

        # Prepare HTML file path
        html_filename = f"{structure.id}.html"
        html_path = os.path.join(self.pulled_structures_dir, html_filename)

        # Generate the HTML for NGL
        html_content = f'''
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Protein Structure: {structure.id}</title>
    <script src="ngl.min.js"></script>
    <style>
        html, body, #viewport {{
            width: 100%;
            height: 100%;
            margin: 0;
            padding: 0;
            overflow: hidden;
        }}
        #loading {{
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            font-family: Arial, sans-serif;
            font-size: 18px;
            color: #666;
        }}
    </style>
</head>
<body>
    <div id="viewport"></div>
    <div id="loading">Loading structure...</div>
    <script>
        // Global variables
        var stage = null;
        var currentComponent = null;
        var isLoaded = false;
        
        document.addEventListener('DOMContentLoaded', function() {{
            try {{
                stage = new NGL.Stage("viewport");
                
                function loadAndRepresent(pdbFilename, representation) {{
                    stage.removeAllComponents();
                    return stage.loadFile(pdbFilename, {{ defaultRepresentation: false }}).then(function (comp) {{
                        currentComponent = comp;
                        comp.addRepresentation(representation);
                        stage.setSpin(false);
                        stage.autoView();
                        isLoaded = true;
                        
                        // Hide loading indicator
                        var loading = document.getElementById('loading');
                        if (loading) loading.style.display = 'none';
                        
                        // Set color scheme after component is loaded
                        if (window.initialColorScheme) {{
                            currentComponent.eachRepresentation(function (repr) {{
                                repr.setColor(window.initialColorScheme);
                            }});
                        }}
                        
                        console.log('Structure loaded successfully');
                        return comp;
                    }}).catch(function(error) {{
                        console.error('Error loading structure:', error);
                        var loading = document.getElementById('loading');
                        if (loading) loading.textContent = 'Error loading structure';
                    }});
                }}

                // Initial load
                loadAndRepresent("{pdb_filename}", "cartoon");

                // Define global functions for external control
                window.setRepresentation = function(representation) {{
                    if (!isLoaded || !currentComponent) {{
                        console.warn('Structure not loaded yet, cannot set representation');
                        return;
                    }}
                    try {{
                        currentComponent.removeAllRepresentations();
                        currentComponent.addRepresentation(representation);
                        // Update color scheme if we have one
                        if (window.initialColorScheme) {{
                            currentComponent.eachRepresentation(function (repr) {{
                                repr.setColor(window.initialColorScheme);
                            }});
                        }}
                        console.log('Representation set to:', representation);
                    }} catch (error) {{
                        console.error('Error setting representation:', error);
                    }}
                }};

                window.setColorScheme = function(colorScheme) {{
                    if (!isLoaded || !currentComponent) {{
                        console.warn('Structure not loaded yet, cannot set color scheme');
                        return;
                    }}
                    try {{
                        currentComponent.eachRepresentation(function (repr) {{
                            repr.setColor(colorScheme);
                        }});
                        console.log('Color scheme set to:', colorScheme);
                    }} catch (error) {{
                        console.error('Error setting color scheme:', error);
                    }}
                }};

                window.setUniformColor = function(color) {{
                    if (!isLoaded || !currentComponent) {{
                        console.warn('Structure not loaded yet, cannot set uniform color');
                        return;
                    }}
                    try {{
                        currentComponent.eachRepresentation(function (repr) {{
                            repr.setColor(color);
                        }});
                        console.log('Uniform color set to:', color);
                    }} catch (error) {{
                        console.error('Error setting uniform color:', error);
                    }}
                }};

                window.setCustomColor = function(color) {{
                    if (!isLoaded || !currentComponent) {{
                        console.warn('Structure not loaded yet, cannot set custom color');
                        return;
                    }}
                    try {{
                        currentComponent.eachRepresentation(function (repr) {{
                            repr.setColor(color);
                        }});
                        console.log('Custom color set to:', color);
                    }} catch (error) {{
                        console.error('Error setting custom color:', error);
                    }}
                }};

                window.setSpin = function(spinEnabled) {{
                    if (!stage) {{
                        console.warn('Stage not initialized yet, cannot set spin');
                        return;
                    }}
                    try {{
                        stage.setSpin(spinEnabled);
                        console.log('Spin set to:', spinEnabled);
                    }} catch (error) {{
                        console.error('Error setting spin:', error);
                    }}
                }};

                window.setBackgroundColor = function(color) {{
                    if (!stage) {{
                        console.warn('Stage not initialized yet, cannot set background color');
                        return;
                    }}
                    try {{
                        stage.viewer.setBackground(color);
                        console.log('Background color set to:', color);
                    }} catch (error) {{
                        console.error('Error setting background color:', error);
                    }}
                }};
                
                window.clearAllRepresentations = function() {{
                    if (!isLoaded || !currentComponent) {{
                        console.warn('Structure not loaded yet, cannot clear representations');
                        return;
                    }}
                    try {{
                        currentComponent.removeAllRepresentations();
                        console.log('All representations cleared');
                    }} catch (error) {{
                        console.error('Error clearing representations:', error);
                    }}
                }};

                window.addEventListener('resize', function(event) {{
                    if (stage) {{
                        stage.handleResize();
                    }}
                }}, false);
                
            }} catch (error) {{
                console.error('Error initializing NGL viewer:', error);
                var loading = document.getElementById('loading');
                if (loading) loading.textContent = 'Error initializing viewer';
            }}
        }});
    </script>
</body>
</html>
'''

        # Write the HTML file, with error handling
        try:
            with open(html_path, "w", encoding='utf-8') as f:
                f.write(html_content)
        except Exception as e:
            self.show_error_dialog(
                "Error Writing HTML",
                f"Could not write HTML viewer file for structure {structure.id}.",
                suggestion="Check directory permissions or disk space.",
                details=str(e)
            )
            print(f"Error writing HTML file: {e}")
            return

        # Load the HTML in the QtWebEngineView
        color_scheme = self.color_combo.currentText() if hasattr(self, 'color_combo') else 'atomindex'
        html_content = html_content.replace("</head>", f"<script>window.initialColorScheme = '{color_scheme}';</script></head>")
        
        with open(html_path, "w", encoding='utf-8') as f:
            f.write(html_content)
        
        self.web_view.setUrl(QUrl(f"http://localhost:{self.port}/data/pulled_structures/{html_filename}"))
        # Wait for page to load before updating color scheme
        def on_load_finished(success):
            if success:
                # Add a small delay to ensure NGL is fully initialized
                def apply_initial_settings():
                    # Apply current settings from UI
                    if hasattr(self, 'color_combo'):
                        current_scheme = self.color_combo.currentText()
                        if current_scheme != color_scheme:
                            self.web_view.page().runJavaScript(f"if (typeof setColorScheme === 'function') setColorScheme('{current_scheme}');")
                    
                    if hasattr(self, 'background_color_entry'):
                        bg_color = self.background_color_entry.text()
                        if bg_color:
                            self.web_view.page().runJavaScript(f"if (typeof setBackgroundColor === 'function') setBackgroundColor('{bg_color}');")
                    
                    if hasattr(self, 'spin_checkbox'):
                        spin_enabled = self.spin_checkbox.isChecked()
                        self.web_view.page().runJavaScript(f"if (typeof setSpin === 'function') setSpin({str(spin_enabled).lower()});")
                    
                    if hasattr(self, 'representation_combo'):
                        representation = self.representation_combo.currentText()
                        if representation != 'cartoon':  # Only change if different from default
                            self.web_view.page().runJavaScript(f"if (typeof setRepresentation === 'function') setRepresentation('{representation}');")
                
                # Use QTimer to delay the application of settings
                from PyQt5.QtCore import QTimer
                QTimer.singleShot(1000, apply_initial_settings)  # 1 second delay
        
        # Disconnect any existing signal handler first
        try:
            self.web_view.page().loadFinished.disconnect()
        except TypeError:
            pass  # No existing connection
        self.web_view.page().loadFinished.connect(on_load_finished)

    def save_screenshot(self):
        # Grab the current view of the QWebEngineView
        pixmap = self.web_view.grab()
        if pixmap.isNull():
            QMessageBox.warning(self, "Screenshot Failed", "Could not capture the current view.")
            return
        # Prompt for file location
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Screenshot", "screenshot.png", "PNG Files (*.png);;All Files (*)")
        if not file_path:
            return
        # Ensure the file has a .png extension
        if not (file_path.lower().endswith('.png')):
            file_path += '.png'
        # Save the screenshot
        if pixmap.save(file_path, 'PNG'):
            QMessageBox.information(self, "Screenshot Saved", f"Screenshot saved to: {file_path}")

    def closeEvent(self, event):
        # Aggressively cleanup QWebEngineView and its page
        if hasattr(self, 'web_view'):
            try:
                # Remove from parent layout if present
                parent_widget = self.web_view.parentWidget()
                if parent_widget is not None:
                    layout = parent_widget.layout()
                    if layout is not None:
                        layout.removeWidget(self.web_view)
                self.web_view.setParent(None)
                self.web_view.close()
                page = self.web_view.page()
                if page is not None:
                    page.deleteLater()
                self.web_view.deleteLater()
                self.web_view = None  # Remove lingering reference
            except Exception as e:
                print(f"Error during web_view cleanup: {e}")
        if hasattr(self, '_server_thread'):
            self._server_thread.shutdown()
        # Process events to ensure deletions are handled
        QApplication.processEvents()
        super().closeEvent(event)

    def add_to_recent_files(self, file_path):
        settings = QSettings("PicoMolApp", "PicoMol")
        recents = settings.value("recent_files", [], type=list)
        if file_path in recents:
            recents.remove(file_path)
        recents.insert(0, file_path)
        recents = recents[:10]
        settings.setValue("recent_files", recents)
        self.update_recent_files_menu()

    def update_recent_files_menu(self):
        self.recent_files_menu.clear()
        settings = QSettings("PicoMolApp", "PicoMol")
        recents = settings.value("recent_files", [], type=list)
        for path in recents:
            action = QAction(path, self)
            action.triggered.connect(lambda checked, p=path: self.open_local_pdb(p))
            self.recent_files_menu.addAction(action)

    def save_structure_as(self):
        if not hasattr(self, 'current_structure_id') or not self.current_structure_id:
            self.statusBar().showMessage("No structure loaded to save.")
            return
        pdb_path = os.path.join(self.pulled_structures_dir, f"{self.current_structure_id}.pdb")
        if not os.path.exists(pdb_path):
            self.statusBar().showMessage("Current structure file not found.")
            return
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Structure As", f"{self.current_structure_id}.pdb", "PDB Files (*.pdb);;All Files (*)")
        if not file_path:
            return
        if not (file_path.lower().endswith('.pdb')):
            file_path += '.pdb'
        try:
            import shutil
            shutil.copy(pdb_path, file_path)
            self.statusBar().showMessage(f"Structure saved to: {file_path}")
            self.add_to_recent_files(file_path)
        except Exception as e:
            self.statusBar().showMessage(f"Failed to save structure: {e}")

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                if url.toLocalFile().lower().endswith(('.pdb', '.ent')):
                    self.drag_overlay.setVisible(True)
                    event.acceptProposedAction()
                    return
        self.drag_overlay.setVisible(False)
        event.ignore()

    def dragLeaveEvent(self, event):
        self.drag_overlay.setVisible(False)
        event.accept()

    def dropEvent(self, event):
        self.drag_overlay.setVisible(False)
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                file_path = url.toLocalFile()
                if file_path.lower().endswith(('.pdb', '.ent')):
                    self.open_local_pdb(file_path)
                    event.acceptProposedAction()
                    return
        event.ignore()

    def setAcceptDrops(self, accept):
        super().setAcceptDrops(accept)

    def _resizeEventWithOverlay(self, event):
        # Ensure drag overlay always covers the full window
        self.drag_overlay.setGeometry(0, 0, self.width(), self.height())
        if hasattr(super(), 'resizeEvent'):
            super().resizeEvent(event)


def main():
    os.environ["QTWEBENGINE_REMOTE_DEBUGGING"] = "9222"

    # Suppress BiopythonParserWarning about missing HEADER line
    warnings.filterwarnings("ignore", message="'HEADER' line not found; can't determine PDB ID.", category=BiopythonParserWarning)

    app = QApplication(sys.argv)

    settings = QWebEngineSettings.defaultSettings()

    # Use the project root directory for serving files
    project_root = os.path.dirname(os.path.abspath(__file__))
    server_thread = ServerThread(directory=project_root)
    server_thread.start()

    # Wait for server to start and get port
    while server_thread.httpd is None:
        time.sleep(0.1)

    viewer = ProteinViewerApp(server_thread.port)
    viewer._server_thread = server_thread
    viewer.show()

    exit_code = app.exec_()
    server_thread.shutdown()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()