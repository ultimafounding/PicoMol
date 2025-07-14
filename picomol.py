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

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QFileDialog, QMessageBox,
    QComboBox, QCheckBox, QGroupBox, QTextEdit
)
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEngineSettings
from PyQt5.QtCore import QUrl

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
        self.setWindowTitle("Basic Protein Structure Viewer")
        self.setGeometry(100, 100, 1000, 800)

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
                    QMessageBox.critical(self, "Dependency Missing", "The 'requests' library is required to download ngl.js. Please install it using 'pip install requests' and restart the application.")
                    sys.exit(1)

                subprocess.run([sys.executable, os.path.join(os.getcwd(), "setup_ngl.py")], check=True)
                QMessageBox.information(self, "NGL.js Setup", "ngl.min.js downloaded successfully. Please restart the application.")
                sys.exit(0)
            except subprocess.CalledProcessError as e:
                QMessageBox.critical(self, "NGL.js Setup Failed", f"Failed to run setup_ngl.py: {e}. Please download ngl.min.js manually or check your internet connection.")
                sys.exit(1)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"An unexpected error occurred during NGL.js setup: {e}")
                sys.exit(1)

        # Copy ngl.min.js to pulled_structures_dir
        shutil.copy(ngl_min_js_path, self.pulled_structures_dir)

        self.port = port
        self.pdb_parser = PDBParser()
        self.pdb_list = PDBList()

        self.init_ui()

    def init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)

        # Status bar
        self.statusBar = self.statusBar()
        self.statusBar.showMessage("Ready")

        # Control Panel
        control_panel_widget = QWidget()
        control_panel_layout = QVBoxLayout(control_panel_widget)

        control_panel_layout.addWidget(QLabel("Enter PDB ID:"))

        self.pdb_id_entry = QLineEdit()
        control_panel_layout.addWidget(self.pdb_id_entry)

        fetch_button = QPushButton("Fetch PDB")
        fetch_button.clicked.connect(self.fetch_pdb_id)
        control_panel_layout.addWidget(fetch_button)

        open_button = QPushButton("Open Local PDB File")
        open_button.clicked.connect(self.open_local_pdb)
        control_panel_layout.addWidget(open_button)

        # NGL.js Options Group
        ngl_options_group = QGroupBox("NGL.js Display Options")
        ngl_layout = QVBoxLayout()

        representation_label = QLabel("Representation:")
        ngl_layout.addWidget(representation_label)
        self.representation_combo = QComboBox()
        self.representation_combo.addItems(["axes", "backbone", "ball+stick", "base", "cartoon", "contact", "distance", "helixorient", "hyperball", "label", "licorice", "line", "point", "ribbon", "rocket", "rope", "spacefill", "surface", "trace", "tube", "unitcell", "validation"])
        self.representation_combo.setCurrentText("cartoon") # Set default to cartoon
        self.representation_combo.currentIndexChanged.connect(self.update_representation)
        ngl_layout.addWidget(self.representation_combo)

        color_label = QLabel("Color Scheme:")
        ngl_layout.addWidget(color_label)
        self.color_combo = QComboBox()
        self.color_combo.addItems(["atomindex", "bfactor", "chainid", "chainindex", "chainname", "densityfit", "electrostatic", "element", "entityindex", "entitytype", "geoquality", "hydrophobicity", "modelindex", "moleculetype", "occupancy", "random", "residueindex", "resname", "sstruc", "uniform", "value", "volume"])
        self.color_combo.currentIndexChanged.connect(self.update_color_scheme)
        ngl_layout.addWidget(self.color_combo)

        self.spin_checkbox = QCheckBox("Spin")
        self.spin_checkbox.setChecked(False)
        self.spin_checkbox.stateChanged.connect(self.toggle_spin)
        ngl_layout.addWidget(self.spin_checkbox)

        self.remove_waters_checkbox = QCheckBox("Remove Waters")
        self.remove_waters_checkbox.setChecked(False)
        self.remove_waters_checkbox.stateChanged.connect(self.toggle_remove_waters)
        ngl_layout.addWidget(self.remove_waters_checkbox)

        # Background Color Option
        background_color_label = QLabel("Background Color (hex or name):")
        ngl_layout.addWidget(background_color_label)
        self.background_color_entry = QLineEdit("black") # Default to black
        ngl_layout.addWidget(self.background_color_entry)
        apply_bg_color_button = QPushButton("Apply Background Color")
        apply_bg_color_button.clicked.connect(self.update_background_color)
        ngl_layout.addWidget(apply_bg_color_button)

        # Custom Color Option
        custom_color_label = QLabel("Custom Color (hex or name):")
        ngl_layout.addWidget(custom_color_label)
        self.custom_color_entry = QLineEdit()
        ngl_layout.addWidget(self.custom_color_entry)
        apply_custom_color_button = QPushButton("Apply Custom Color")
        apply_custom_color_button.clicked.connect(self.update_custom_color)
        ngl_layout.addWidget(apply_custom_color_button)

        ngl_options_group.setLayout(ngl_layout)
        control_panel_layout.addWidget(ngl_options_group)
        control_panel_layout.addStretch(1) # Push everything to the top

        # Sequence Display
        sequence_group = QGroupBox("Sequence Data")
        sequence_layout = QVBoxLayout()
        self.sequence_display = QTextEdit()
        self.sequence_display.setReadOnly(True)
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

    def update_custom_color(self):
        color = self.custom_color_entry.text()
        self.web_view.page().runJavaScript(f"setCustomColor('{color}');")

    def update_background_color(self):
        color = self.background_color_entry.text()
        self.web_view.page().runJavaScript(f"setBackgroundColor('{color}');")

    def update_representation(self):
        representation = self.representation_combo.currentText()
        # This will call a JavaScript function in the web view
        self.web_view.page().runJavaScript(f"setRepresentation('{representation}');")

    def clear_all_representations(self):
        self.web_view.page().runJavaScript("clearAllRepresentations();")

    def update_color_scheme(self):
        color_scheme = self.color_combo.currentText()
        # This will call a JavaScript function in the web view
        self.web_view.page().runJavaScript(f"setColorScheme('{color_scheme}');")

    def toggle_spin(self):
        spin_enabled = self.spin_checkbox.isChecked()
        # This will call a JavaScript function in the web view
        self.web_view.page().runJavaScript(f"setSpin({str(spin_enabled).lower()});")

    def toggle_remove_waters(self):
        remove_waters = self.remove_waters_checkbox.isChecked()
        self.web_view.page().runJavaScript(f"setRemoveWaters({str(remove_waters).lower()});")

    def update_background_color(self):
        color = self.background_color_entry.text()
        self.web_view.page().runJavaScript(f"setBackgroundColor('{color}');")

    def fetch_pdb_id(self):
        pdb_id = self.pdb_id_entry.text().strip().upper()
        if not pdb_id:
            self.statusBar.showMessage("Please enter a PDB ID.")
            return
        try:
            self.statusBar.showMessage(f"Fetching PDB ID {pdb_id}...")
            # Define a consistent path for the PDB file.
            pdb_path = os.path.join(self.pulled_structures_dir, f"{pdb_id}.pdb")

            # Retrieve PDB file only if it doesn't exist
            if not os.path.exists(pdb_path):
                retrieved_file_path = self.pdb_list.retrieve_pdb_file(
                    pdb_id, pdir=self.pulled_structures_dir, file_format="pdb"
                )
                # Rename the downloaded file to our consistent path, overwriting if necessary.
                if os.path.exists(pdb_path):
                    os.remove(pdb_path)
                os.rename(retrieved_file_path, pdb_path)

            self.statusBar.showMessage(f"Loading structure {pdb_id}...")
            structure = self.pdb_parser.get_structure(pdb_id, pdb_path)
            self.display_structure(structure)
            self.statusBar.showMessage(f"Displayed {pdb_id}")
        except Exception as e:
            self.statusBar.showMessage(f"Error fetching PDB ID {pdb_id}")
            QMessageBox.critical(self, "Error", f"Could not fetch PDB ID {pdb_id}:\n{e}")

    def open_local_pdb(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Local PDB File", "", "PDB Files (*.pdb *.ent);;All Files (*)"
        )
        if not file_path:
            return
        try:
            self.statusBar.showMessage(f"Opening local PDB file: {os.path.basename(file_path)}...")
            structure_id = os.path.basename(file_path).split(".")[0]
            # To display correctly, copy the file to the served directory
            pdb_filename = f"{structure_id}.pdb"
            target_path = os.path.join(self.pulled_structures_dir, pdb_filename)

            # If the file is not already in the correct location, copy it.
            if os.path.abspath(file_path) != os.path.abspath(target_path):
                import shutil
                shutil.copy(file_path, target_path)

            self.statusBar.showMessage(f"Loading structure from {os.path.basename(file_path)}...")
            structure = self.pdb_parser.get_structure(structure_id, target_path)
            self.display_structure(structure)
            self.statusBar.showMessage(f"Displayed {structure_id}")
        except Exception as e:
            self.statusBar.showMessage(f"Error opening file: {os.path.basename(file_path)}")
            QMessageBox.critical(self, "Error", f"Could not open file:\n{e}")

    def display_structure(self, structure):
        from pathlib import Path
        from Bio import SeqIO
        from io import StringIO
        io = PDBIO()
        io.set_structure(structure)

        # Save structure PDB file
        pdb_filename = f"{structure.id}.pdb"
        pdb_path = os.path.join(self.pulled_structures_dir, pdb_filename)
        io.save(pdb_path)

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
                }});
            }}

            // Initial load
            loadAndRepresent("{pdb_filename}", "cartoon");

            window.setRepresentation = function(representation) {{
                if (currentComponent) {{
                    currentComponent.removeAllRepresentations(); // Clear existing representations
                    currentComponent.addRepresentation(representation);
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
            }};

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

        with open(html_path, "w") as f:
            f.write(html_content)

        # Load the HTML in the QtWebEngineView
        self.web_view.setUrl(QUrl(f"http://localhost:{self.port}/pulled_structures/{html_filename}"))


    def closeEvent(self, event):
        if hasattr(self, '_server_thread'):
            self._server_thread.shutdown()
        super().closeEvent(event)


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