from PyQt6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QFileDialog, QTextEdit, QGraphicsView, QGraphicsScene
from PyQt6.QtGui import QPainter, QPen, QBrush, QColor, QFont
from PyQt6.QtCore import Qt
from Bio import SeqIO
import math

class PlasmidViewer(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        
        self.load_button = QPushButton("Load Plasmid File")
        self.load_button.clicked.connect(self.load_plasmid_file)
        self.layout.addWidget(self.load_button)
        
        self.plasmid_view = QGraphicsView()
        self.plasmid_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.layout.addWidget(self.plasmid_view)

        self.plasmid_details = QTextEdit()
        self.plasmid_details.setReadOnly(True)
        self.layout.addWidget(self.plasmid_details)

    def load_plasmid_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Plasmid File", "", "GenBank Files (*.gb *.gbk);;FASTA Files (*.fa *.fasta);;All Files (*)"
        )
        if not file_path:
            return

        try:
            record = SeqIO.read(file_path, "genbank") # Assuming GenBank for now
            self.display_plasmid(record)
        except Exception as e:
            self.plasmid_details.setText(f"Error parsing file: {e}")

    def display_plasmid(self, record):
        self.plasmid_details.clear()
        self.plasmid_details.append(f"ID: {record.id}")
        self.plasmid_details.append(f"Name: {record.name}")
        self.plasmid_details.append(f"Description: {record.description}")
        self.plasmid_details.append(f"Length: {len(record.seq)} bp")
        self.plasmid_details.append("\nFeatures:")
        for feature in record.features:
            self.plasmid_details.append(f"- {feature.type} ({feature.location})")

        scene = QGraphicsScene()
        self.plasmid_view.setScene(scene)

        plasmid_len = len(record.seq)
        radius = 100
        center_x, center_y = 0, 0

        # Draw plasmid backbone
        pen = QPen(Qt.GlobalColor.black, 2)
        scene.addEllipse(center_x - radius, center_y - radius, radius * 2, radius * 2, pen)

        # Ensure the full plasmid is visible and centered
        scene.setSceneRect(scene.itemsBoundingRect().adjusted(-20, -20, 20, 20))
        self.plasmid_view.fitInView(scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)

        # Draw features
        for feature in record.features:
            if feature.type in ["gene", "promoter", "CDS"]:
                start = feature.location.start
                end = feature.location.end
                
                start_angle = (start / plasmid_len) * 360
                span_angle = ((end - start) / plasmid_len) * 360

                # Use QColor constants
                pen_color = QColor(0, 0, 255)  # blue
                if feature.type == "promoter":
                    pen_color = QColor(255, 0, 0)  # red
                elif feature.type == "CDS":
                    pen_color = QColor(0, 128, 0)  # green

                pen = QPen(pen_color, 10)
                # QGraphicsScene doesn't have addArc; draw using addEllipse with a QPainterPath or use addPath
                # For now add a simple arc placeholder using addEllipse with pen width to represent feature
                scene.addEllipse(center_x - radius, center_y - radius, radius * 2, radius * 2, pen)

def create_plasmid_viewer_tab(parent):
    """Creates the Plasmid Viewer tab."""
    return PlasmidViewer(parent)

