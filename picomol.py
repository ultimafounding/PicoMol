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

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QFileDialog, QMessageBox,
    QComboBox, QCheckBox, QGroupBox, QTextEdit, QDialog, QDialogButtonBox, QAction
)


class AboutDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About PicoMol")
        self.setMinimumWidth(420)
        layout = QVBoxLayout(self)
        about_text = QLabel(
            """
            <h2>PicoMol</h2>
            <p><b>Version:</b> 0.0.2 (2025-07-15)</p>
            <p>A simple, user-friendly molecular visualization tool for protein structures.<br>
            Built with PyQt5 and NGL.js, powered by Biopython.</p>
            <ul>
              <li>Fetch and visualize PDB structures</li>
              <li>Open local PDB files</li>
              <li>Undo/redo, screenshots, color schemes, and more</li>
            </ul>
            <p><b>Credits:</b><br>
            Developed by Jack Magson.<br>
            Uses <a href='https://github.com/arose/ngl'>NGL.js</a> (MIT License) for 3D visualization.<br>
            Powered by <a href='https://biopython.org/'>Biopython</a> and <a href='https://riverbankcomputing.com/software/pyqt/intro'>PyQt5</a>.<br>
            </p>
            <p><b>License:</b> GNU GPL v3.0<br>
            See the LICENSE file for details.</p>
            <p><b>NGL.js citation:</b><br>
            AS Rose, AR Bradley, Y Valasatava, JM Duarte, A Prlić and PW Rose. NGL viewer: web-based molecular graphics for large complexes. <i>Bioinformatics</i>: bty419, 2018. <a href='https://doi.org/10.1093/bioinformatics/bty419'>doi:10.1093/bioinformatics/bty419</a><br>
            AS Rose and PW Hildebrand. NGL Viewer: a web application for molecular visualization. <i>Nucleic Acids Res</i> (1 July 2015) 43 (W1): W576-W579 first published online April 29, 2015. <a href='https://doi.org/10.1093/nar/gkv402'>doi:10.1093/nar/gkv402</a>
            </p>
            """
        )
        about_text.setOpenExternalLinks(True)
        about_text.setWordWrap(True)
        layout.addWidget(about_text)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(self.accept)
        layout.addWidget(button_box)

from PyQt5.QtCore import QSettings
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEngineSettings
from PyQt5.QtCore import QUrl, Qt

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


class WelcomeDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Welcome to PicoMol!")
        self.setMinimumWidth(400)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        welcome_label = QLabel(
            "<h2>Welcome to PicoMol!</h2>"
            "<p>This is a simple molecular visualization tool for protein structures.</p>"
            "<ul>"
            "<li>Fetch PDB structures using their ID</li>"
            "<li>Open local PDB files</li>"
            "<li>Adjust visualization and color schemes</li>"
            "<li>Take screenshots and more!</li>"
            "</ul>"
            "<p>Hover over any control for tooltips and inline help.</p>"
        )
        welcome_label.setWordWrap(True)
        layout.addWidget(welcome_label)

        self.checkbox = QCheckBox("Show this welcome screen on startup")
        self.checkbox.setChecked(True)
        layout.addWidget(self.checkbox)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(self.accept)
        layout.addWidget(button_box)

    def should_show_next_time(self):
        return self.checkbox.isChecked()



class ProteinViewerApp(QMainWindow):
    def __init__(self, port):
        super().__init__()
        self.port = port
        self._undo_stack = []
        self._redo_stack = []
        self._is_restoring_state = False
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
        self.pulled_structures_dir = os.path.join(os.getcwd(), "pulled_structures")
        os.makedirs(self.pulled_structures_dir, exist_ok=True)
        self.ngl_assets_dir = os.path.join(os.getcwd(), "ngl_assets")
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
        self.pulled_structures_dir = os.path.join(os.getcwd(), "pulled_structures")
        os.makedirs(self.pulled_structures_dir, exist_ok=True)
        self.ngl_assets_dir = os.path.join(os.getcwd(), "ngl_assets")
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
                self.web_view.page().runJavaScript(f"setCustomColor('{val}');")
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

        # Help menu with About
        help_menu = menubar.addMenu("Help")
        about_action = QAction("About PicoMol", self)
        about_action.setToolTip("Show version info, credits, and license.")
        about_action.triggered.connect(self.show_about_dialog)
        help_menu.addAction(about_action)

        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)

        # Drag-and-drop overlay label (hidden by default)
        self.drag_overlay = QLabel("\n\nDrop PDB or ENT files here to open", self)
        self.drag_overlay.setAlignment(Qt.AlignCenter)
        self.drag_overlay.setStyleSheet("background: rgba(30, 144, 255, 0.7); color: white; font-size: 28px; border: 3px dashed white; border-radius: 24px;")
        self.drag_overlay.setVisible(False)
        self.drag_overlay.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.drag_overlay.setGeometry(0, 0, self.width(), self.height())

        # Ensure overlay resizes with window
        self.resizeEvent = self._resizeEventWithOverlay

        # Status bar
        self.statusBar().showMessage("Ready")

        # Control Panel
        control_panel_widget = QWidget()
        control_panel_layout = QVBoxLayout(control_panel_widget)

        control_panel_layout.addWidget(QLabel("Enter PDB ID:"))

        self.pdb_id_entry = QLineEdit()
        self.pdb_id_entry.setToolTip("Enter a valid PDB ID (e.g., 1CRN, 4HHB) to fetch a protein structure from the PDB database.")
        control_panel_layout.addWidget(self.pdb_id_entry)

        fetch_button = QPushButton("Fetch PDB")
        fetch_button.setToolTip("Download and visualize a protein structure using the entered PDB ID.")
        fetch_button.clicked.connect(self.fetch_pdb_id)
        control_panel_layout.addWidget(fetch_button)

        open_button = QPushButton("Open Local PDB File")
        open_button.setToolTip("Open and visualize a local .pdb or .ent file from your computer.")
        open_button.clicked.connect(self.open_local_pdb)
        control_panel_layout.addWidget(open_button)

        # Screenshot Button
        screenshot_button = QPushButton("Save Screenshot")
        screenshot_button.setToolTip("Save a screenshot of the current protein structure view as a PNG image.")
        screenshot_button.clicked.connect(self.save_screenshot)
        control_panel_layout.addWidget(screenshot_button)

        # NGL.js Options Group
        ngl_options_group = QGroupBox("NGL.js Display Options")
        ngl_layout = QVBoxLayout()

        representation_label = QLabel("Representation:")
        ngl_layout.addWidget(representation_label)
        self.representation_combo = QComboBox()
        self.representation_combo.setToolTip("Select the 3D representation style for the protein structure (e.g., cartoon, surface, ball+stick, etc.).")
        self.representation_combo.addItems(["axes", "backbone", "ball+stick", "base", "cartoon", "contact", "distance", "helixorient", "hyperball", "label", "licorice", "line", "point", "ribbon", "rocket", "rope", "spacefill", "surface", "trace", "tube", "unitcell", "validation"])
        self.representation_combo.setCurrentText("cartoon") # Set default to cartoon
        self.representation_combo.currentIndexChanged.connect(self.update_representation)
        ngl_layout.addWidget(self.representation_combo)

        color_label = QLabel("Color Scheme:")
        ngl_layout.addWidget(color_label)
        self.color_combo = QComboBox()
        self.color_combo.setToolTip("Choose a color scheme for the structure (e.g., by atom, chain, residue, etc.).")
        self.color_combo.addItems(["atomindex", "bfactor", "chainid", "chainindex", "chainname", "densityfit", "electrostatic", "element", "entityindex", "entitytype", "geoquality", "hydrophobicity", "modelindex", "moleculetype", "occupancy", "random", "residueindex", "resname", "sstruc", "uniform", "value", "volume"])
        self.color_combo.currentIndexChanged.connect(self.update_color_scheme)
        ngl_layout.addWidget(self.color_combo)

        self.spin_checkbox = QCheckBox("Spin")
        self.spin_checkbox.setToolTip("Toggle automatic rotation (spin) of the 3D protein structure.")
        self.spin_checkbox.setChecked(False)
        self.spin_checkbox.stateChanged.connect(self.toggle_spin)
        ngl_layout.addWidget(self.spin_checkbox)



        # Background Color Option
        background_color_label = QLabel("Background Color (hex or name):")
        ngl_layout.addWidget(background_color_label)
        self.background_color_entry = QLineEdit("black") # Default to black
        self.background_color_entry.setToolTip("Enter a color name or hex code for the background (e.g., 'black', '#ffffff').")
        ngl_layout.addWidget(self.background_color_entry)
        apply_bg_color_button = QPushButton("Apply Background Color")
        apply_bg_color_button.setToolTip("Apply the chosen background color to the viewer.")
        apply_bg_color_button.clicked.connect(self.update_background_color)
        ngl_layout.addWidget(apply_bg_color_button)

        # Custom Color Option
        custom_color_label = QLabel("Custom Color (hex or name):")
        ngl_layout.addWidget(custom_color_label)
        self.custom_color_entry = QLineEdit()
        self.custom_color_entry.setToolTip("Enter a custom color (name or hex code) to apply to the protein structure.")
        self._last_custom_color = self.custom_color_entry.text()
        ngl_layout.addWidget(self.custom_color_entry)
        apply_custom_color_button = QPushButton("Apply Custom Color")

        # Initial undo stack state after UI setup
        self._undo_stack.clear()
        self._redo_stack.clear()
        self._undo_stack.append(self.capture_state())
        print(f"[UNDO] Initial stack: {[s['custom_color'] for s in self._undo_stack]}")
        apply_custom_color_button.setToolTip("Apply the custom color to the protein structure.")
        apply_custom_color_button.clicked.connect(self.update_custom_color)
        ngl_layout.addWidget(apply_custom_color_button)

        ngl_options_group.setLayout(ngl_layout)
        control_panel_layout.addWidget(ngl_options_group)
        control_panel_layout.addStretch(1) # Push everything to the top

        # Feedback Button
        feedback_button = QPushButton("Send Feedback")
        feedback_button.setToolTip("Report a bug or suggest an improvement.")
        feedback_button.clicked.connect(self.show_feedback_dialog)
        control_panel_layout.addWidget(feedback_button)

        # Sequence Display
        sequence_group = QGroupBox("Sequence Data")
        sequence_layout = QVBoxLayout()
        self.sequence_display = QTextEdit()
        self.sequence_display.setReadOnly(True)
        self.sequence_display.setToolTip("Displays the amino acid sequence of the loaded protein structure.")
        sequence_layout.addWidget(self.sequence_display)
        sequence_group.setLayout(sequence_layout)
        sequence_group.setMaximumHeight(100) # Set a maximum height for the sequence box

        self.web_view = QWebEngineView()

        # Create a new vertical layout for sequence and web view
        right_panel_layout = QVBoxLayout()
        right_panel_layout.addWidget(sequence_group, 1) # Give sequence_group a smaller stretch factor
        right_panel_layout.addWidget(self.web_view, 3) # Give web_view a larger stretch factor

        main_layout.addWidget(control_panel_widget)
        main_layout.addLayout(right_panel_layout, 3) # Give right_panel_layout a stretch factor of 3

        # Initialize undo stack with initial state
        self._undo_stack.clear()
        self._redo_stack.clear()
        self._undo_stack.append(self.capture_state())

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

    def update_custom_color(self):
        color = self.custom_color_entry.text()
        print(f"[UNDO] update_custom_color: Applying color {color}")
        self.web_view.page().runJavaScript(f"setCustomColor('{color}');")
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
        self.web_view.page().runJavaScript(f"setBackgroundColor('{color}');")
        if not self._is_restoring_state:
            prev_color = None
        if self._undo_stack:
            prev_color = self._undo_stack[-1].get('background_color', None)
        if color != prev_color:
            self.push_undo()

    def update_representation(self):
        self.push_undo()
        representation = self.representation_combo.currentText()
        # This will call a JavaScript function in the web view
        self.web_view.page().runJavaScript(f"setRepresentation('{representation}');")

    def clear_all_representations(self):
        self.web_view.page().runJavaScript("clearAllRepresentations();")

    def update_color_scheme(self):
        self.push_undo()
        color_scheme = self.color_combo.currentText()
        # This will call a JavaScript function in the web view
        self.web_view.page().runJavaScript(f"setColorScheme('{color_scheme}');")

    def toggle_spin(self):
        self.push_undo()
        spin_enabled = self.spin_checkbox.isChecked()
        # This will call a JavaScript function in the web view
        self.web_view.page().runJavaScript(f"setSpin({str(spin_enabled).lower()});")

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
    </style>
</head>
<body>
    <div id="viewport"></div>
    <script>
        document.addEventListener('DOMContentLoaded', function() {{
            var stage = new NGL.Stage( "viewport" );
            var currentComponent = null;

            function loadAndRepresent(pdbFilename, representation) {{
                stage.removeAllComponents();
                stage.loadFile( pdbFilename, {{ defaultRepresentation: false }} ).then(function (comp) {{
                    currentComponent = comp;
                    comp.addRepresentation(representation);
                    stage.setSpin(false);
                    stage.autoView();
                    // Set color scheme after component is loaded
                    if (window.initialColorScheme) {{
                        currentComponent.eachRepresentation(function (repr) {{
                            repr.setColor(window.initialColorScheme);
                        }});
                    }}
                }});
            }}

            // Initial load
            loadAndRepresent("{pdb_filename}", "cartoon");

            window.setRepresentation = function(representation) {{
                if (currentComponent) {{
                    currentComponent.removeAllRepresentations(); // Clear existing representations
                    currentComponent.addRepresentation(representation);
                    // Update color scheme if we have one
                    if (window.initialColorScheme) {{
                        currentComponent.eachRepresentation(function (repr) {{
                            repr.setColor(window.initialColorScheme);
                        }});
                    }}
                }} else {{
                    loadAndRepresent("{pdb_filename}", representation);
                }}
            }};

            window.setColorScheme = function(colorScheme) {{
                if (currentComponent) {{
                    currentComponent.eachRepresentation(function (repr) {{
                        repr.setColor(colorScheme);
                    }});
                }}
            }}

            window.setCustomColor = function(color) {{
                if (currentComponent) {{
                    currentComponent.eachRepresentation(function (repr) {{
                        repr.setColor(color);
                    }});
                }}
            }};

            window.setSpin = function(spinEnabled) {{
                stage.setSpin(spinEnabled);
            }};

            window.setBackgroundColor = function(color) {{
                stage.viewer.setBackground(color);
            }};

            window.addEventListener( "resize", function( event ) {{
                stage.handleResize();
            }}, false );
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
        
        self.web_view.setUrl(QUrl(f"http://localhost:{self.port}/pulled_structures/{html_filename}"))
        # Wait for page to load before updating color scheme
        def on_load_finished(success):
            if success:
                # Ensure we have a color combo
                if hasattr(self, 'color_combo'):
                    current_scheme = self.color_combo.currentText()
                    # Only update if it's different from initial
                    if current_scheme != color_scheme:
                        self.web_view.page().runJavaScript(f"setColorScheme('{current_scheme}');")
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