from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
                             QComboBox, QSpinBox, QCheckBox, QGroupBox, QFormLayout,
                             QFileDialog, QMessageBox, QProgressDialog, QColorDialog,
                             QSlider, QTabWidget, QWidget, QTextEdit)
from PyQt5.QtGui import QPixmap, QPainter, QColor, QPen, QBrush, QFont
from PyQt5.QtCore import Qt, QSize, QThread, pyqtSignal
from PyQt5.QtSvg import QSvgGenerator
import os

class ExportWorker(QThread):
    """Worker thread for export operations"""
    progress = pyqtSignal(int)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    
    def __init__(self, export_func, *args, **kwargs):
        super().__init__()
        self.export_func = export_func
        self.args = args
        self.kwargs = kwargs
    
    def run(self):
        try:
            result = self.export_func(*self.args, **self.kwargs)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))

class ExportDialog(QDialog):
    """Enhanced export dialog for plasmid maps and sequences"""
    
    def __init__(self, sequence_viewer, parent=None):
        super().__init__(parent)
        self.sequence_viewer = sequence_viewer
        self.setWindowTitle("Export Sequence and Maps")
        self.setMinimumSize(500, 600)
        self.setup_ui()
    
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Create tab widget for different export types
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)
        
        # Map export tab
        self.create_map_export_tab()
        
        # Sequence export tab
        self.create_sequence_export_tab()
        
        # Advanced options tab
        self.create_advanced_tab()
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.preview_button = QPushButton("Preview")
        self.preview_button.clicked.connect(self.show_preview)
        button_layout.addWidget(self.preview_button)
        
        button_layout.addStretch()
        
        self.export_button = QPushButton("Export")
        self.export_button.clicked.connect(self.export)
        button_layout.addWidget(self.export_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
    
    def create_map_export_tab(self):
        """Create map export options tab"""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Format selection
        format_group = QGroupBox("Export Format")
        format_layout = QFormLayout(format_group)
        
        self.format_combo = QComboBox()
        self.format_combo.addItems(["PNG (Raster)", "SVG (Vector)", "PDF (Vector)"])
        format_layout.addRow("Format:", self.format_combo)
        
        layout.addWidget(format_group)
        
        # View selection
        view_group = QGroupBox("View Type")
        view_layout = QFormLayout(view_group)
        
        self.view_combo = QComboBox()
        self.view_combo.addItems(["Circular Map", "Linear Map", "Both Views"])
        view_layout.addRow("View:", self.view_combo)
        
        layout.addWidget(view_group)
        
        # Size and quality
        size_group = QGroupBox("Size and Quality")
        size_layout = QFormLayout(size_group)
        
        self.width_spin = QSpinBox()
        self.width_spin.setRange(400, 4000)
        self.width_spin.setValue(1200)
        self.width_spin.setSuffix(" px")
        size_layout.addRow("Width:", self.width_spin)
        
        self.height_spin = QSpinBox()
        self.height_spin.setRange(400, 4000)
        self.height_spin.setValue(800)
        self.height_spin.setSuffix(" px")
        size_layout.addRow("Height:", self.height_spin)
        
        self.dpi_spin = QSpinBox()
        self.dpi_spin.setRange(72, 600)
        self.dpi_spin.setValue(300)
        size_layout.addRow("DPI:", self.dpi_spin)
        
        layout.addWidget(size_group)
        
        # Elements to include
        elements_group = QGroupBox("Elements to Include")
        elements_layout = QVBoxLayout(elements_group)
        
        self.include_features = QCheckBox("Features")
        self.include_features.setChecked(True)
        elements_layout.addWidget(self.include_features)
        
        self.include_enzymes = QCheckBox("Restriction Sites")
        self.include_enzymes.setChecked(True)
        elements_layout.addWidget(self.include_enzymes)
        
        self.include_labels = QCheckBox("Feature Labels")
        self.include_labels.setChecked(True)
        elements_layout.addWidget(self.include_labels)
        
        self.include_ruler = QCheckBox("Scale/Ruler")
        self.include_ruler.setChecked(True)
        elements_layout.addWidget(self.include_ruler)
        
        self.include_title = QCheckBox("Title and Info")
        self.include_title.setChecked(True)
        elements_layout.addWidget(self.include_title)
        
        layout.addWidget(elements_group)
        
        layout.addStretch()
        self.tab_widget.addTab(tab, "Map Export")
    
    def create_sequence_export_tab(self):
        """Create sequence export options tab"""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Format selection
        format_group = QGroupBox("Sequence Format")
        format_layout = QFormLayout(format_group)
        
        self.seq_format_combo = QComboBox()
        self.seq_format_combo.addItems([
            "FASTA", "GenBank", "EMBL", "Plain Text", 
            "SnapGene-style (HTML)", "Annotated Text"
        ])
        format_layout.addRow("Format:", self.seq_format_combo)
        
        layout.addWidget(format_group)
        
        # Sequence options
        seq_group = QGroupBox("Sequence Options")
        seq_layout = QVBoxLayout(seq_group)
        
        self.include_complement_seq = QCheckBox("Include Complement Strand")
        seq_layout.addWidget(self.include_complement_seq)
        
        self.include_translation_seq = QCheckBox("Include Translation")
        seq_layout.addWidget(self.include_translation_seq)
        
        self.include_features_seq = QCheckBox("Include Feature Annotations")
        self.include_features_seq.setChecked(True)
        seq_layout.addWidget(self.include_features_seq)
        
        self.include_enzymes_seq = QCheckBox("Include Restriction Sites")
        seq_layout.addWidget(self.include_enzymes_seq)
        
        layout.addWidget(seq_group)
        
        # Formatting options
        format_group = QGroupBox("Text Formatting")
        format_layout = QFormLayout(format_group)
        
        self.bases_per_line_spin = QSpinBox()
        self.bases_per_line_spin.setRange(20, 120)
        self.bases_per_line_spin.setValue(60)
        format_layout.addRow("Bases per line:", self.bases_per_line_spin)
        
        self.line_numbers = QCheckBox("Include line numbers")
        self.line_numbers.setChecked(True)
        format_layout.addRow("", self.line_numbers)
        
        layout.addWidget(format_group)
        
        layout.addStretch()
        self.tab_widget.addTab(tab, "Sequence Export")
    
    def create_advanced_tab(self):
        """Create advanced options tab"""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Color customization
        color_group = QGroupBox("Color Customization")
        color_layout = QFormLayout(color_group)
        
        self.background_color = QColor(255, 255, 255)
        self.bg_color_button = QPushButton()
        self.bg_color_button.setStyleSheet(f"background-color: {self.background_color.name()}")
        self.bg_color_button.clicked.connect(self.choose_background_color)
        color_layout.addRow("Background:", self.bg_color_button)
        
        layout.addWidget(color_group)
        
        # Font options
        font_group = QGroupBox("Font Options")
        font_layout = QFormLayout(font_group)
        
        self.font_size_spin = QSpinBox()
        self.font_size_spin.setRange(8, 24)
        self.font_size_spin.setValue(12)
        font_layout.addRow("Font Size:", self.font_size_spin)
        
        self.font_family_combo = QComboBox()
        self.font_family_combo.addItems(["Consolas", "Courier New", "Monaco", "Menlo"])
        font_layout.addRow("Font Family:", self.font_family_combo)
        
        layout.addWidget(font_group)
        
        # Metadata
        metadata_group = QGroupBox("Metadata")
        metadata_layout = QVBoxLayout(metadata_group)
        
        self.custom_title = QTextEdit()
        self.custom_title.setMaximumHeight(60)
        self.custom_title.setPlaceholderText("Custom title (leave empty for default)")
        metadata_layout.addWidget(QLabel("Custom Title:"))
        metadata_layout.addWidget(self.custom_title)
        
        self.author_field = QTextEdit()
        self.author_field.setMaximumHeight(40)
        self.author_field.setPlaceholderText("Author/Creator")
        metadata_layout.addWidget(QLabel("Author:"))
        metadata_layout.addWidget(self.author_field)
        
        layout.addWidget(metadata_group)
        
        layout.addStretch()
        self.tab_widget.addTab(tab, "Advanced")
    
    def choose_background_color(self):
        """Choose background color"""
        color = QColorDialog.getColor(self.background_color, self)
        if color.isValid():
            self.background_color = color
            self.bg_color_button.setStyleSheet(f"background-color: {color.name()}")
    
    def show_preview(self):
        """Show export preview"""
        # Create preview dialog
        preview_dialog = QDialog(self)
        preview_dialog.setWindowTitle("Export Preview")
        preview_dialog.setMinimumSize(600, 400)
        
        layout = QVBoxLayout(preview_dialog)
        
        # Generate preview
        preview_pixmap = self.generate_preview()
        
        if preview_pixmap:
            preview_label = QLabel()
            # Scale preview to fit dialog
            scaled_pixmap = preview_pixmap.scaled(
                550, 350, Qt.KeepAspectRatio, Qt.SmoothTransformation
            )
            preview_label.setPixmap(scaled_pixmap)
            preview_label.setAlignment(Qt.AlignCenter)
            layout.addWidget(preview_label)
        else:
            layout.addWidget(QLabel("Preview not available"))
        
        # Close button
        close_button = QPushButton("Close")
        close_button.clicked.connect(preview_dialog.accept)
        layout.addWidget(close_button)
        
        preview_dialog.exec_()
    
    def generate_preview(self):
        """Generate a preview of the export"""
        try:
            if self.tab_widget.currentIndex() == 0:  # Map export
                return self.generate_map_preview()
            else:  # Sequence export
                return self.generate_sequence_preview()
        except Exception as e:
            QMessageBox.warning(self, "Preview Error", f"Could not generate preview: {str(e)}")
            return None
    
    def generate_map_preview(self):
        """Generate map preview"""
        if not self.sequence_viewer.record:
            return None
        
        # Create a small version for preview
        pixmap = QPixmap(400, 300)
        pixmap.fill(self.background_color)
        
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.Antialiasing)
        
        # Draw a simplified version of the map
        if self.view_combo.currentText() == "Circular Map":
            self.draw_circular_preview(painter, pixmap.size())
        elif self.view_combo.currentText() == "Linear Map":
            self.draw_linear_preview(painter, pixmap.size())
        else:  # Both views
            # Draw both in split view
            self.draw_circular_preview(painter, QSize(200, 300))
            painter.translate(200, 0)
            self.draw_linear_preview(painter, QSize(200, 300))
        
        painter.end()
        return pixmap
    
    def draw_circular_preview(self, painter, size):
        """Draw circular map preview"""
        center_x, center_y = size.width() // 2, size.height() // 2
        radius = min(size.width(), size.height()) // 3
        
        # Draw backbone
        painter.setPen(QPen(Qt.black, 2))
        painter.drawEllipse(center_x - radius, center_y - radius, radius * 2, radius * 2)
        
        # Draw features if enabled
        if self.include_features.isChecked() and self.sequence_viewer.record:
            self.draw_features_preview(painter, center_x, center_y, radius, circular=True)
        
        # Add title if enabled
        if self.include_title.isChecked():
            painter.setPen(QPen(Qt.black))
            painter.setFont(QFont("Arial", 8))
            title = self.custom_title.toPlainText() or f"{self.sequence_viewer.record.id}"
            painter.drawText(10, 20, title)
    
    def draw_linear_preview(self, painter, size):
        """Draw linear map preview"""
        y_center = size.height() // 2
        line_length = size.width() - 40
        
        # Draw backbone
        painter.setPen(QPen(Qt.black, 2))
        painter.drawLine(20, y_center, 20 + line_length, y_center)
        
        # Draw features if enabled
        if self.include_features.isChecked() and self.sequence_viewer.record:
            self.draw_features_preview(painter, 20, y_center, line_length, circular=False)
        
        # Add title if enabled
        if self.include_title.isChecked():
            painter.setPen(QPen(Qt.black))
            painter.setFont(QFont("Arial", 8))
            title = self.custom_title.toPlainText() or f"{self.sequence_viewer.record.id}"
            painter.drawText(10, 20, title)
    
    def draw_features_preview(self, painter, center_x, center_y, radius_or_length, circular=True):
        """Draw simplified features for preview"""
        if not hasattr(self.sequence_viewer.record, 'features'):
            return
        
        seq_len = len(self.sequence_viewer.record.seq)
        colors = {
            "gene": QColor(255, 255, 0),
            "CDS": QColor(255, 165, 0),
            "promoter": QColor(255, 0, 0),
            "terminator": QColor(128, 0, 128),
            "rep_origin": QColor(0, 255, 0)
        }
        
        for feature in self.sequence_viewer.record.features[:5]:  # Limit for preview
            if feature.type in colors:
                color = colors[feature.type]
                painter.setPen(QPen(color, 3))
                
                start = int(feature.location.start)
                end = int(feature.location.end)
                
                if circular:
                    # Draw arc for circular
                    start_angle = (start / seq_len) * 360
                    span_angle = ((end - start) / seq_len) * 360
                    painter.drawArc(
                        center_x - radius_or_length, center_y - radius_or_length,
                        radius_or_length * 2, radius_or_length * 2,
                        int(-start_angle * 16), int(-span_angle * 16)
                    )
                else:
                    # Draw line for linear
                    start_x = center_x + (start / seq_len) * radius_or_length
                    end_x = center_x + (end / seq_len) * radius_or_length
                    painter.drawLine(int(start_x), center_y - 10, int(end_x), center_y - 10)
    
    def generate_sequence_preview(self):
        """Generate sequence preview"""
        if not self.sequence_viewer.record:
            return None
        
        # Create text preview
        pixmap = QPixmap(400, 300)
        pixmap.fill(Qt.white)
        
        painter = QPainter(pixmap)
        painter.setFont(QFont("Consolas", 8))
        painter.setPen(QPen(Qt.black))
        
        # Show first few lines of sequence
        sequence = str(self.sequence_viewer.record.seq)
        bases_per_line = self.bases_per_line_spin.value()
        
        y = 20
        for i in range(0, min(len(sequence), bases_per_line * 8), bases_per_line):
            line = sequence[i:i + bases_per_line]
            if self.line_numbers.isChecked():
                text = f"{i+1:>6} {line}"
            else:
                text = line
            painter.drawText(10, y, text)
            y += 15
            if y > 280:
                break
        
        painter.end()
        return pixmap
    
    def export(self):
        """Perform the export"""
        if not self.sequence_viewer.record:
            QMessageBox.warning(self, "No Sequence", "No sequence loaded to export.")
            return
        
        # Get file path
        if self.tab_widget.currentIndex() == 0:  # Map export
            format_text = self.format_combo.currentText()
            if "PNG" in format_text:
                filter_str = "PNG Images (*.png)"
                default_ext = ".png"
            elif "SVG" in format_text:
                filter_str = "SVG Images (*.svg)"
                default_ext = ".svg"
            else:  # PDF
                filter_str = "PDF Documents (*.pdf)"
                default_ext = ".pdf"
        else:  # Sequence export
            format_text = self.seq_format_combo.currentText()
            if "FASTA" in format_text:
                filter_str = "FASTA Files (*.fasta *.fa)"
                default_ext = ".fasta"
            elif "GenBank" in format_text:
                filter_str = "GenBank Files (*.gb *.gbk)"
                default_ext = ".gb"
            elif "HTML" in format_text:
                filter_str = "HTML Files (*.html)"
                default_ext = ".html"
            else:
                filter_str = "Text Files (*.txt)"
                default_ext = ".txt"
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export File", 
            f"{self.sequence_viewer.record.id}{default_ext}",
            filter_str
        )
        
        if not file_path:
            return
        
        # Show progress dialog
        progress = QProgressDialog("Exporting...", "Cancel", 0, 100, self)
        progress.setWindowModality(Qt.WindowModal)
        progress.show()
        
        try:
            if self.tab_widget.currentIndex() == 0:  # Map export
                self.export_map(file_path, progress)
            else:  # Sequence export
                self.export_sequence(file_path, progress)
            
            progress.setValue(100)
            QMessageBox.information(self, "Export Complete", f"Successfully exported to:\n{file_path}")
            self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Export failed:\n{str(e)}")
        finally:
            progress.close()
    
    def export_map(self, file_path, progress):
        """Export map to file"""
        progress.setValue(20)
        
        width = self.width_spin.value()
        height = self.height_spin.value()
        
        if file_path.endswith('.svg'):
            # SVG export
            generator = QSvgGenerator()
            generator.setFileName(file_path)
            generator.setSize(QSize(width, height))
            generator.setViewBox(QRect(0, 0, width, height))
            generator.setTitle(self.custom_title.toPlainText() or f"{self.sequence_viewer.record.id}")
            generator.setDescription("Plasmid map exported from PicoMol")
            
            painter = QPainter(generator)
        else:
            # Raster export (PNG/PDF)
            pixmap = QPixmap(width, height)
            pixmap.fill(self.background_color)
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.Antialiasing)
        
        progress.setValue(40)
        
        # Draw the map
        self.draw_full_map(painter, QSize(width, height))
        
        progress.setValue(80)
        
        painter.end()
        
        if not file_path.endswith('.svg'):
            # Save raster image
            pixmap.save(file_path, "PNG" if file_path.endswith('.png') else "PDF")
    
    def draw_full_map(self, painter, size):
        """Draw the full quality map"""
        # This would contain the full implementation of map drawing
        # For now, use the preview drawing as a placeholder
        if self.view_combo.currentText() == "Circular Map":
            self.draw_circular_preview(painter, size)
        elif self.view_combo.currentText() == "Linear Map":
            self.draw_linear_preview(painter, size)
        else:  # Both views
            half_width = size.width() // 2
            self.draw_circular_preview(painter, QSize(half_width, size.height()))
            painter.translate(half_width, 0)
            self.draw_linear_preview(painter, QSize(half_width, size.height()))
    
    def export_sequence(self, file_path, progress):
        """Export sequence to file"""
        progress.setValue(20)
        
        format_text = self.seq_format_combo.currentText()
        
        if "FASTA" in format_text:
            self.export_fasta(file_path)
        elif "GenBank" in format_text:
            self.export_genbank(file_path)
        elif "HTML" in format_text:
            self.export_html(file_path)
        else:
            self.export_text(file_path)
        
        progress.setValue(100)
    
    def export_fasta(self, file_path):
        """Export as FASTA"""
        from Bio import SeqIO
        with open(file_path, 'w') as f:
            SeqIO.write(self.sequence_viewer.record, f, "fasta")
    
    def export_genbank(self, file_path):
        """Export as GenBank"""
        from Bio import SeqIO
        with open(file_path, 'w') as f:
            SeqIO.write(self.sequence_viewer.record, f, "genbank")
    
    def export_html(self, file_path):
        """Export as HTML with SnapGene-style formatting"""
        sequence = str(self.sequence_viewer.record.seq)
        bases_per_line = self.bases_per_line_spin.value()
        
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>{self.sequence_viewer.record.id}</title>
            <style>
                body {{ font-family: {self.font_family_combo.currentText()}, monospace; }}
                .sequence {{ font-size: {self.font_size_spin.value()}px; }}
                .A {{ color: red; }}
                .T {{ color: blue; }}
                .C {{ color: green; }}
                .G {{ color: orange; }}
                .line-number {{ color: gray; margin-right: 10px; }}
            </style>
        </head>
        <body>
            <h1>{self.custom_title.toPlainText() or self.sequence_viewer.record.id}</h1>
            <p>Length: {len(sequence)} bp</p>
            <div class="sequence">
        """
        
        for i in range(0, len(sequence), bases_per_line):
            line = sequence[i:i + bases_per_line]
            if self.line_numbers.isChecked():
                html += f'<span class="line-number">{i+1:>6}</span>'
            
            for base in line:
                html += f'<span class="{base}">{base}</span>'
            html += '<br>\n'
        
        html += """
            </div>
        </body>
        </html>
        """
        
        with open(file_path, 'w') as f:
            f.write(html)
    
    def export_text(self, file_path):
        """Export as plain text"""
        sequence = str(self.sequence_viewer.record.seq)
        bases_per_line = self.bases_per_line_spin.value()
        
        with open(file_path, 'w') as f:
            f.write(f"{self.custom_title.toPlainText() or self.sequence_viewer.record.id}\n")
            f.write(f"Length: {len(sequence)} bp\n\n")
            
            for i in range(0, len(sequence), bases_per_line):
                line = sequence[i:i + bases_per_line]
                if self.line_numbers.isChecked():
                    f.write(f"{i+1:>6} {line}\n")
                else:
                    f.write(f"{line}\n")