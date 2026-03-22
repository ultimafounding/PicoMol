from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
                             QComboBox, QSpinBox, QCheckBox, QGroupBox, QFormLayout,
                             QFileDialog, QMessageBox, QProgressDialog, QColorDialog,
                             QSlider, QTabWidget, QWidget, QTextEdit, QGraphicsView)
from PyQt6.QtGui import QPixmap, QPainter, QColor, QPen, QBrush, QFont
from PyQt6.QtCore import Qt, QSize, QThread, pyqtSignal, QRectF
from PyQt6.QtSvg import QSvgGenerator
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
        format_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 13px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        format_layout = QFormLayout(format_group)
        format_layout.setContentsMargins(10, 15, 10, 10)
        format_layout.setSpacing(8)
        
        self.format_combo = QComboBox()
        self.format_combo.addItems(["PNG (Raster)", "SVG (Vector)", "PDF (Vector)"])
        format_layout.addRow("Format:", self.format_combo)
        
        layout.addWidget(format_group)
        
        # View selection
        view_group = QGroupBox("View Type")
        view_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 13px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        view_layout = QFormLayout(view_group)
        view_layout.setContentsMargins(10, 15, 10, 10)
        view_layout.setSpacing(8)
        
        self.view_combo = QComboBox()
        self.view_combo.addItems(["Circular Map", "Linear Map", "Both Views"])
        view_layout.addRow("View:", self.view_combo)
        
        layout.addWidget(view_group)
        
        # Size and quality
        size_group = QGroupBox("Size and Quality")
        size_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 13px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        size_layout = QFormLayout(size_group)
        size_layout.setContentsMargins(10, 15, 10, 10)
        size_layout.setSpacing(8)
        
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
        self.dpi_spin.setRange(72, 1200)
        self.dpi_spin.setValue(600)
        size_layout.addRow("DPI:", self.dpi_spin)
        
        layout.addWidget(size_group)
        
        layout.addStretch()
        self.tab_widget.addTab(tab, "Map Export")
    
    def create_advanced_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Color customization
        color_group = QGroupBox("Color Customization")
        color_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 13px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        color_layout = QFormLayout(color_group)
        color_layout.setContentsMargins(10, 15, 10, 10)
        color_layout.setSpacing(8)
        
        self.background_color = QColor(255, 255, 255)
        self.bg_color_button = QPushButton()
        self.bg_color_button.setStyleSheet(f"background-color: {self.background_color.name()}")
        self.bg_color_button.clicked.connect(self.choose_background_color)
        color_layout.addRow("Background:", self.bg_color_button)
        
        layout.addWidget(color_group)
        
        # Font options
        font_group = QGroupBox("Font Options")
        font_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 13px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        font_layout = QFormLayout(font_group)
        font_layout.setContentsMargins(10, 15, 10, 10)
        font_layout.setSpacing(8)
        
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
        metadata_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 13px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        metadata_layout = QVBoxLayout(metadata_group)
        metadata_layout.setContentsMargins(10, 15, 10, 10)
        metadata_layout.setSpacing(8)
        
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
                550, 350, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
            )
            preview_label.setPixmap(scaled_pixmap)
            preview_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(preview_label)
        else:
            layout.addWidget(QLabel("Preview not available"))
        
        # Close button
        close_button = QPushButton("Close")
        close_button.clicked.connect(preview_dialog.accept)
        layout.addWidget(close_button)
        
        preview_dialog.exec()
    
    def generate_preview(self):
        """Generate a preview of the export"""
        try:
            # Only map export is available now
            return self.generate_map_preview()
        except Exception as e:
            QMessageBox.warning(self, "Preview Error", f"Could not generate preview: {str(e)}")
            return None
    
    def generate_map_preview(self):
        """Generate map preview using actual scene rendering"""
        if not self.sequence_viewer.record:
            return None
        
        # Force the sequence viewer to update the current view first
        if self.view_combo.currentText() == "Circular Map":
            self.sequence_viewer.show_circular_view()
        elif self.view_combo.currentText() == "Linear Map":
            self.sequence_viewer.show_linear_view()
        else:  # Both views
            self.sequence_viewer.show_circular_view()  # Start with circular for preview
        
        # Create a small version for preview
        pixmap = QPixmap(400, 300)
        pixmap.fill(self.background_color)
        
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        
        # Use the same rendering logic as the full export
        self.draw_full_map(painter, pixmap.size())
        
        painter.end()
        return pixmap
    
    def draw_circular_preview(self, painter, size):
        """Draw circular map preview"""
        center_x, center_y = size.width() // 2, size.height() // 2
        radius = min(size.width(), size.height()) // 3
        
        # Draw backbone
        painter.setPen(QPen(Qt.GlobalColor.black, 2))
        painter.drawEllipse(center_x - radius, center_y - radius, radius * 2, radius * 2)
        
        # Draw features if enabled
        if self.include_features.isChecked() and self.sequence_viewer.record:
            self.draw_features_preview(painter, center_x, center_y, radius, circular=True)
        
        # Add title if enabled
        if self.include_title.isChecked():
            painter.setPen(QPen(Qt.GlobalColor.black))
            painter.setFont(QFont("Arial", 8))
            title = self.custom_title.toPlainText() or f"{self.sequence_viewer.record.id}"
            painter.drawText(10, 20, title)
    
    def draw_linear_preview(self, painter, size):
        """Draw linear map preview"""
        y_center = size.height() // 2
        line_length = size.width() - 40
        
        # Draw backbone
        painter.setPen(QPen(Qt.GlobalColor.black, 2))
        painter.drawLine(20, y_center, 20 + line_length, y_center)
        
        # Draw features if enabled
        if self.include_features.isChecked() and self.sequence_viewer.record:
            self.draw_features_preview(painter, 20, y_center, line_length, circular=False)
        
        # Add title if enabled
        if self.include_title.isChecked():
            painter.setPen(QPen(Qt.GlobalColor.black))
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
        pixmap.fill(Qt.GlobalColor.white)
        
        painter = QPainter(pixmap)
        painter.setFont(QFont("Consolas", 8))
        painter.setPen(QPen(Qt.GlobalColor.black))
        
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
        """Perform the map export"""
        if not self.sequence_viewer.record:
            QMessageBox.warning(self, "No Sequence", "No sequence loaded to export.")
            return
        
        # Get file path for map export
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
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Map", 
            f"{self.sequence_viewer.record.id}{default_ext}",
            filter_str
        )
        
        if not file_path:
            return
        
        # Show progress dialog
        progress = QProgressDialog("Exporting map...", "Cancel", 0, 100, self)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.show()
        
        try:
            self.export_map(file_path, progress)
            progress.setValue(100)
            QMessageBox.information(self, "Export Complete", f"Successfully exported map to:\n{file_path}")
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
        dpi = self.dpi_spin.value()
        
        if file_path.endswith('.svg'):
            # SVG export
            generator = QSvgGenerator()
            generator.setFileName(file_path)
            generator.setSize(QSize(width, height))
            generator.setViewBox(QRect(0, 0, width, height))
            generator.setTitle(self.custom_title.toPlainText() or f"{self.sequence_viewer.record.id}")
            generator.setDescription("Plasmid map exported from PicoMol")
            
            painter = QPainter(generator)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
            
            # Draw the map at target size
            self.draw_full_map(painter, QSize(width, height))
            
            painter.end()
        else:
            # Raster export (PNG/PDF) - scale pixel dimensions based on DPI
            # Base DPI is 96, so scale factor = target_dpi / 96
            scale_factor = dpi / 96.0
            pixel_width = int(width * scale_factor)
            pixel_height = int(height * scale_factor)
            
            # Create high-resolution pixmap
            pixmap = QPixmap(pixel_width, pixel_height)
            pixmap.fill(self.background_color)
            
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
            
            # Scale painter to work in logical units (width x height)
            painter.scale(scale_factor, scale_factor)
            
            # Draw the map at logical size (width x height)
            self.draw_full_map(painter, QSize(width, height))
            
            painter.end()
            
            # Save with maximum quality
            if file_path.endswith('.png'):
                pixmap.save(file_path, "PNG", quality=100)
            elif file_path.endswith('.pdf'):
                pixmap.save(file_path, "PDF", quality=100)
        
        progress.setValue(100)
    
    def draw_full_map(self, painter, size):
        """Draw the full quality map using the actual sequence viewer rendering"""
        # Get the current scene from the sequence viewer
        scene = self.sequence_viewer.scene
        if not scene:
            # Fallback to preview if no scene
            if self.view_combo.currentText() == "Circular Map":
                self.draw_circular_preview(painter, size)
            elif self.view_combo.currentText() == "Linear Map":
                self.draw_linear_preview(painter, size)
            else:  # Both views
                half_width = size.width() // 2
                self.draw_circular_preview(painter, QSize(half_width, size.height()))
                painter.translate(half_width, 0)
                self.draw_linear_preview(painter, QSize(half_width, size.height()))
            return
        
        # Force view update and save current state
        current_view = self.sequence_viewer.current_view if hasattr(self.sequence_viewer, 'current_view') else 'circular'
        
        if self.view_combo.currentText() == "Circular Map":
            # Force circular view update
            self.sequence_viewer.show_circular_view()
            # Render circular view with dynamic sizing
            self.render_view_to_painter_dynamic(self.sequence_viewer, painter, size)
            
        elif self.view_combo.currentText() == "Linear Map":
            # Force linear view update
            self.sequence_viewer.show_linear_view()
            # Render linear view with dynamic sizing
            self.render_view_to_painter_dynamic(self.sequence_viewer, painter, size)
            
        else:  # Both views
            # Render both views side by side with dynamic sizing
            half_width = size.width() // 2
            
            # Render circular view on left half
            self.sequence_viewer.show_circular_view()
            
            # Render circular view in left half with offset
            left_size = QSize(half_width, size.height())
            self.render_view_to_painter_dynamic(self.sequence_viewer, painter, left_size, offset_x=0, offset_y=0)
            
            # Render linear view on right half
            self.sequence_viewer.show_linear_view()
            
            # Render linear view in right half with offset
            right_size = QSize(half_width, size.height())
            self.render_view_to_painter_dynamic(self.sequence_viewer, painter, right_size, offset_x=half_width, offset_y=0)
        
        # Restore original view
        if current_view == 'circular':
            self.sequence_viewer.show_circular_view()
        else:
            self.sequence_viewer.show_linear_view()
    
    def render_view_to_painter_dynamic(self, sequence_viewer, painter, target_size, offset_x=0, offset_y=0):
        """Render a QGraphicsView to a painter with proper sizing"""
        # Get the scene and calculate proper bounds
        scene = sequence_viewer.scene
        scene_rect = scene.itemsBoundingRect()
        
        if scene_rect.isEmpty():
            return
        
        # Add padding to prevent cutting off edges
        padding = 60
        scene_rect = scene_rect.adjusted(-padding, -padding, padding, padding)
        
        # Calculate scale to fit the content in the target area
        margin = 30
        available_width = target_size.width() - (2 * margin)
        available_height = target_size.height() - (2 * margin)
        
        scale_x = available_width / scene_rect.width()
        scale_y = available_height / scene_rect.height()
        scale = min(scale_x, scale_y)
        
        # Calculate the scaled content size
        scaled_width = scene_rect.width() * scale
        scaled_height = scene_rect.height() * scale
        
        # Center the content in the target area
        content_x = offset_x + (target_size.width() - scaled_width) / 2
        content_y = offset_y + (target_size.height() - scaled_height) / 2
        
        # Create a target rectangle for rendering
        target_rect = QRectF(content_x, content_y, scaled_width, scaled_height)
        
        # Render the scene directly to the target rectangle
        scene.render(painter, target_rect, scene_rect)
        
        # Draw a border for debugging (remove this later)
        painter.setPen(QPen(Qt.GlobalColor.red, 1))
        painter.drawRect(target_rect)