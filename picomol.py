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
            <p><b>Version:</b> 0.0.3 (2025-07-18)</p>
            <p>A comprehensive molecular visualization and bioinformatics suite for protein structures and sequence analysis.</p>
            
            <h3>🔬 Core Features</h3>
            <ul>
              <li><b>3D Molecular Visualization:</b> Interactive protein structure viewing with multiple representations</li>
              <li><b>BLAST Integration:</b> Complete BLAST suite (BLASTP, BLASTN, BLASTX, TBLASTN, TBLASTX)</li>
              <li><b>PDB Support:</b> Fetch from PDB database or load local files</li>
              <li><b>Sequence Analysis:</b> View and analyze protein/nucleotide sequences</li>
              <li><b>Export Tools:</b> Save screenshots and structure files</li>
            </ul>
            
            <h3>🛠️ Technology Stack</h3>
            <ul>
              <li><b>GUI Framework:</b> <a href='https://riverbankcomputing.com/software/pyqt/intro'>PyQt5</a> with QWebEngine</li>
              <li><b>3D Visualization:</b> <a href='https://github.com/arose/ngl'>NGL.js</a> (MIT License)</li>
              <li><b>Bioinformatics:</b> <a href='https://biopython.org/'>Biopython</a></li>
              <li><b>BLAST:</b> NCBI BLAST API integration</li>
            </ul>
            
            <p><b>Developer:</b> Jack Magson<br>
            <b>License:</b> GNU GPL v3.0 (see LICENSE file)<br>
            <b>Repository:</b> <a href='https://github.com/ultimafounding/PicoMol'>GitHub</a></p>
            
            <h3>📚 Citations</h3>
            
            <p><b>NGL.js:</b></p>
            <div class="citation">
            AS Rose, AR Bradley, Y Valasatava, JM Duarte, A Prlić and PW Rose. 
            NGL viewer: web-based molecular graphics for large complexes. 
            <i>Bioinformatics</i> 34(21): 3755-3758, 2018. 
            <a href='https://doi.org/10.1093/bioinformatics/bty419'>doi:10.1093/bioinformatics/bty419</a>
            </div>
            
            <div class="citation">
            AS Rose and PW Hildebrand. NGL Viewer: a web application for molecular visualization. 
            <i>Nucleic Acids Research</i> 43(W1): W576-W579, 2015. 
            <a href='https://doi.org/10.1093/nar/gkv402'>doi:10.1093/nar/gkv402</a>
            </div>
            
            <p><b>BLAST:</b></p>
            <div class="citation">
            Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 
            Basic local alignment search tool. 
            <i>Journal of Molecular Biology</i> 215(3): 403-410, 1990.
            <a href='https://doi.org/10.1016/S0022-2836(05)80360-2'>doi:10.1016/S0022-2836(05)80360-2</a>
            </div>
            
            <div class="citation">
            Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. 
            Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. 
            <i>Nucleic Acids Research</i> 25(17): 3389-3402, 1997.
            <a href='https://doi.org/10.1093/nar/25.17.3389'>doi:10.1093/nar/25.17.3389</a>
            </div>
            
            <p><b>PDB Format:</b></p>
            <div class="citation">
            Berman HM, Westbrook J, Feng Z, Gilliland G, Bhat TN, Weissig H, 
            Shindyalov IN, Bourne PE. 
            The Protein Data Bank. 
            <i>Nucleic Acids Research</i> 28(1): 235-242, 2000.
            <a href='https://doi.org/10.1093/nar/28.1.235'>doi:10.1093/nar/28.1.235</a>
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
        about_text.setOpenExternalLinks(True)
        about_text.setWordWrap(True)
        layout.addWidget(about_text)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(self.accept)
        layout.addWidget(button_box)

from PyQt5.QtCore import QSettings
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEngineSettings
from PyQt5.QtCore import QUrl, Qt
from PyQt5.QtGui import QFont

from Bio.PDB import PDBList, PDBParser, PDBIO
from Bio.SeqIO.PdbIO import BiopythonParserWarning


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
        self.setGeometry(100, 100, 1000, 800)
        self.setAcceptDrops(True)
        
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
        self.pulled_structures_dir = os.path.join(os.getcwd(), "data", "pulled_structures")
        os.makedirs(self.pulled_structures_dir, exist_ok=True)
        self.ngl_assets_dir = os.path.join(os.getcwd(), "assets", "ngl_assets")
        os.makedirs(self.ngl_assets_dir, exist_ok=True)

        ngl_min_js_path = os.path.join(self.ngl_assets_dir, "ngl.min.js")
        if not os.path.exists(ngl_min_js_path):
            QMessageBox.information(self, "NGL.js Missing", "ngl.min.js not found. Attempting to run setup_ngl.py...")
            try:
                # Check if requests is installed
                try:
                    import requests
                except ImportError:
                    self.show_error_dialog(
                        "Dependency Missing",
                        "The 'requests' library is required to download ngl.js.",
                        suggestion="Please install it using 'pip install requests' and restart the application."
                    )
                    sys.exit(1)

                subprocess.run([sys.executable, os.path.join(os.getcwd(), "setup_ngl.py")], check=True)
                QMessageBox.information(self, "NGL.js Setup", "ngl.min.js downloaded successfully. Please restart the application.")
                sys.exit(0)
            except subprocess.CalledProcessError as e:
                self.show_error_dialog(
                    "NGL.js Setup Failed",
                    "Failed to run setup_ngl.py.",
                    suggestion="Please download ngl.min.js manually or check your internet connection.",
                    details=str(e)
                )
                sys.exit(1)
            except Exception as e:
                self.show_error_dialog(
                    "Error",
                    "An unexpected error occurred during NGL.js setup.",
                    details=str(e)
                )
                sys.exit(1)

        # Copy ngl.min.js to pulled_structures_dir
        shutil.copy(ngl_min_js_path, self.pulled_structures_dir)

        self.pdb_parser = PDBParser()
        self.pdb_list = PDBList()
        self.init_ui()

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

        # Undo/redo stacks
        self._undo_stack = []
        self._redo_stack = []
        self._is_restoring_state = False

        super().__init__()
        self.setWindowTitle("Basic Protein Structure Viewer")
        self.setGeometry(100, 100, 1000, 800)
        self.setAcceptDrops(True)

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
        self.pulled_structures_dir = os.path.join(os.getcwd(), "data", "pulled_structures")
        os.makedirs(self.pulled_structures_dir, exist_ok=True)
        self.ngl_assets_dir = os.path.join(os.getcwd(), "assets", "ngl_assets")
        os.makedirs(self.ngl_assets_dir, exist_ok=True)

        ngl_min_js_path = os.path.join(self.ngl_assets_dir, "ngl.min.js")
        if not os.path.exists(ngl_min_js_path):
            QMessageBox.information(self, "NGL.js Missing", "ngl.min.js not found. Attempting to run setup_ngl.py...")
            try:
                # Check if requests is installed
                try:
                    import requests
                except ImportError:
                    self.show_error_dialog(
    "Dependency Missing",
    "The 'requests' library is required to download ngl.js.",
    suggestion="Please install it using 'pip install requests' and restart the application."
)
                    sys.exit(1)

                subprocess.run([sys.executable, os.path.join(os.getcwd(), "setup_ngl.py")], check=True)
                QMessageBox.information(self, "NGL.js Setup", "ngl.min.js downloaded successfully. Please restart the application.")
                sys.exit(0)
            except subprocess.CalledProcessError as e:
                self.show_error_dialog(
    "NGL.js Setup Failed",
    "Failed to run setup_ngl.py.",
    suggestion="Please download ngl.min.js manually or check your internet connection.",
    details=str(e)
)
                sys.exit(1)
            except Exception as e:
                self.show_error_dialog(
    "Error",
    "An unexpected error occurred during NGL.js setup.",
    details=str(e)
)
                sys.exit(1)

        # Copy ngl.min.js to pulled_structures_dir
        shutil.copy(ngl_min_js_path, self.pulled_structures_dir)

        self.pdb_parser = PDBParser()
        self.pdb_list = PDBList()

        self.init_ui()

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
        self.sequence_display.setFixedHeight(75)  # Set fixed height to 50px
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
        
        # Create Sequence Analysis tab (empty for now, will be implemented later)
        seq_tab = QWidget()
        seq_layout = QVBoxLayout(seq_tab)
        seq_label = QLabel("Sequence Analysis Tools (Coming Soon)")
        seq_label.setAlignment(Qt.AlignCenter)
        seq_label.setStyleSheet("font-size: 16px; font-weight: bold; color: #666; margin-top: 20px;")
        seq_layout.addWidget(seq_label)
        seq_layout.addStretch()
        self.main_tabs.addTab(seq_tab, "Sequence Tools")
        
        # Create Structural Analysis tab
        struct_tab = QWidget()
        struct_layout = QVBoxLayout(struct_tab)
        struct_label = QLabel("Structural Analysis (Coming in Next Release)")
        struct_label.setAlignment(Qt.AlignCenter)
        struct_label.setStyleSheet("font-size: 16px; font-weight: bold; color: #666; margin-top: 20px;")
        struct_layout.addWidget(struct_label)
        struct_layout.addStretch()
        self.main_tabs.addTab(struct_tab, "Structure")
        
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
        try:
            self.statusBar().showMessage(f"Fetching PDB ID {pdb_id}...")
            # Define a consistent path for the PDB file.
            pdb_path = os.path.join(self.pulled_structures_dir, f"{pdb_id}.pdb")

            # Retrieve PDB file only if it doesn't exist
            if not os.path.exists(pdb_path):
                try:
                    retrieved_file_path = self.pdb_list.retrieve_pdb_file(
                        pdb_id, pdir=self.pulled_structures_dir, file_format="pdb"
                    )
                    # If we get here but no file was actually downloaded
                    if not os.path.exists(retrieved_file_path):
                        raise FileNotFoundError(f"PDB file for {pdb_id} was not downloaded")
                    # Rename the downloaded file to our consistent path, overwriting if necessary.
                    if os.path.exists(pdb_path):
                        os.remove(pdb_path)
                    os.rename(retrieved_file_path, pdb_path)
                except Exception as e:
                    # If there's any error during download, ensure we don't have a partial file
                    if os.path.exists(pdb_path):
                        try:
                            os.remove(pdb_path)
                        except:
                            pass
                    raise

            # If we get here, we should have a valid PDB file
            if not os.path.exists(pdb_path):
                raise FileNotFoundError(f"Failed to create PDB file for {pdb_id}")

            structure = self.pdb_parser.get_structure(pdb_id, pdb_path)
            self.current_structure_id = pdb_id  # Set the current structure ID
            self.display_structure(structure)
            self.statusBar().showMessage(f"Displayed {pdb_id}")
            self.add_to_recent_files(pdb_path)

        except Exception as e:
            self.statusBar().showMessage(f"Error fetching PDB ID {pdb_id}")
            self.show_error_dialog(
                "Error Fetching PDB",
                f"Could not fetch PDB ID {pdb_id}.",
                suggestion="Please check your internet connection and verify the PDB ID is correct. The application will now close due to an unresolved error.",
                details=str(e)
            )
            # Close the application after showing the error
            QApplication.quit()
            sys.exit(1)

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
            self.display_structure(structure)
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
        from Bio import SeqIO
        from io import StringIO
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

        # Extract and display sequence
        try:
            # Read the PDB file content into a string buffer
            with open(pdb_path, 'r') as f:
                pdb_content = f.read()
            pdb_buffer = StringIO(pdb_content)

            # Parse the PDB file to extract sequence information
            sequences = []
            for i, record in enumerate(SeqIO.parse(pdb_buffer, "pdb-atom")): # Use "pdb-atom" for parsing PDB files
                # Construct a more informative ID using the structure.id and chain ID
                display_id = f"{structure.id}:{record.id.split(':')[-1]}" if ':' in record.id else f"{structure.id}:{record.id}"
                sequences.append(f">{display_id}\n{record.seq}")
            self.sequence_display.setText("\n".join(sequences))
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
            with open(html_path, "w") as f:
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
        
        with open(html_path, "w") as f:
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

    server_thread = ServerThread(directory=os.getcwd())
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