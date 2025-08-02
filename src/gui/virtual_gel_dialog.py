from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QGraphicsView, QGraphicsScene, 
                             QHBoxLayout, QLabel, QComboBox, QPushButton, QGroupBox,
                             QSpinBox, QFormLayout, QTextEdit, QSplitter)
from PyQt5.QtGui import QPen, QBrush, QColor, QFont, QLinearGradient
from PyQt5.QtCore import Qt
import math

class VirtualGelDialog(QDialog):
    def __init__(self, fragments, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Virtual Agarose Gel Electrophoresis")
        self.setMinimumSize(700, 600)
        self.fragments = fragments
        
        # Main layout
        main_layout = QHBoxLayout(self)
        
        # Left panel for controls
        controls_panel = QGroupBox("Gel Parameters")
        controls_layout = QFormLayout(controls_panel)
        
        # Gel concentration
        self.gel_concentration = QComboBox()
        self.gel_concentration.addItems(["0.8%", "1.0%", "1.2%", "1.5%", "2.0%"])
        self.gel_concentration.setCurrentText("1.0%")
        self.gel_concentration.currentTextChanged.connect(self.update_gel)
        controls_layout.addRow("Agarose concentration:", self.gel_concentration)
        
        # Voltage
        self.voltage = QSpinBox()
        self.voltage.setRange(50, 200)
        self.voltage.setValue(100)
        self.voltage.setSuffix(" V")
        self.voltage.valueChanged.connect(self.update_gel)
        controls_layout.addRow("Voltage:", self.voltage)
        
        # Run time
        self.run_time = QSpinBox()
        self.run_time.setRange(30, 180)
        self.run_time.setValue(60)
        self.run_time.setSuffix(" min")
        self.run_time.valueChanged.connect(self.update_gel)
        controls_layout.addRow("Run time:", self.run_time)
        
        # Marker selection
        self.marker_combo = QComboBox()
        self.marker_combo.addItems(["1kb DNA Ladder", "100bp DNA Ladder", "λ DNA/HindIII", "No Marker"])
        self.marker_combo.currentTextChanged.connect(self.update_gel)
        controls_layout.addRow("DNA Marker:", self.marker_combo)
        
        # Update button
        update_button = QPushButton("Update Gel")
        update_button.clicked.connect(self.update_gel)
        controls_layout.addRow(update_button)
        
        main_layout.addWidget(controls_panel)
        
        # Right panel with gel view and fragment info
        right_panel = QSplitter(Qt.Vertical)
        
        # Gel view
        self.view = QGraphicsView()
        self.scene = QGraphicsScene()
        self.view.setScene(self.scene)
        self.view.setMinimumHeight(400)
        right_panel.addWidget(self.view)
        
        # Fragment information
        info_group = QGroupBox("Fragment Information")
        info_layout = QVBoxLayout(info_group)
        self.fragment_info = QTextEdit()
        self.fragment_info.setReadOnly(True)
        self.fragment_info.setMaximumHeight(150)
        info_layout.addWidget(self.fragment_info)
        right_panel.addWidget(info_group)
        
        main_layout.addWidget(right_panel)
        
        # Set initial sizes
        main_layout.setStretch(0, 0)  # Controls panel fixed width
        main_layout.setStretch(1, 1)  # Gel view expandable
        
        self.draw_gel()
        self.update_fragment_info()

    def get_marker_fragments(self, marker_type):
        """Get standard DNA marker fragment sizes"""
        markers = {
            "1kb DNA Ladder": [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250],
            "100bp DNA Ladder": [1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
            "λ DNA/HindIII": [23130, 9416, 6557, 4361, 2322, 2027, 564, 125],
            "No Marker": []
        }
        return markers.get(marker_type, [])
    
    def calculate_migration_distance(self, fragment_size, gel_concentration, voltage, time):
        """Calculate migration distance based on fragment size and gel parameters"""
        # More realistic migration model based on molecular sieving
        # Smaller fragments migrate further through gel matrix
        
        base_distance = 300  # Maximum migration distance in pixels
        
        # Handle edge cases
        if fragment_size <= 0:
            return 0
        
        # Size factor - logarithmic relationship with better scaling
        # Fragments 100-10000 bp should have good separation
        log_size = math.log10(max(fragment_size, 50))
        
        # Normalize to 100-10000 bp range
        log_min = math.log10(100)   # ~2.0
        log_max = math.log10(10000) # ~4.0
        
        # Invert so smaller fragments migrate further
        size_factor = max(0, (log_max - log_size) / (log_max - log_min))
        
        # Gel concentration factor (higher concentration = less migration)
        concentration = float(gel_concentration.rstrip('%'))
        conc_factor = 1.5 / concentration  # 1.5% is reference
        
        # Voltage and time factor (linear relationship)
        vt_factor = (voltage / 100) * (time / 60)
        
        # Apply molecular sieving effects
        # Very large fragments (>20kb) migrate very little
        if fragment_size > 20000:
            size_factor *= 0.3
        elif fragment_size > 10000:
            size_factor *= 0.6
        
        # Very small fragments (<100bp) may run off gel
        if fragment_size < 100:
            size_factor = min(1.2, size_factor * 1.5)
        
        # Calculate final distance
        distance = base_distance * size_factor * conc_factor * vt_factor * 0.4
        
        # Add some realistic variation
        import random
        variation = random.uniform(0.95, 1.05)  # ±5% variation
        distance *= variation
        
        return min(distance, base_distance * 0.95)  # Cap at 95% of gel length
    
    def draw_gel(self):
        """Draw the virtual gel with fragments and markers"""
        self.scene.clear()
        
        if not self.fragments:
            # Show empty gel with instructions
            self.scene.addText("No fragments to display.\nLoad sequences and perform restriction digest.", 
                             QFont("Arial", 12))
            return
        
        # Gel parameters
        gel_width = 400
        gel_height = 350
        well_width = 30
        well_height = 15
        
        # Draw gel background with gradient
        gel_rect = self.scene.addRect(50, 50, gel_width, gel_height, 
                                     QPen(Qt.black, 2), QBrush(QColor(240, 245, 255)))
        
        # Draw buffer chambers
        # Top buffer (negative electrode)
        top_buffer = self.scene.addRect(50, 30, gel_width, 20,
                                       QPen(Qt.darkGray), QBrush(QColor(200, 200, 255)))
        # Bottom buffer (positive electrode)
        bottom_buffer = self.scene.addRect(50, 400, gel_width, 20,
                                          QPen(Qt.darkGray), QBrush(QColor(255, 200, 200)))
        
        # Calculate lanes
        marker_type = self.marker_combo.currentText()
        has_marker = marker_type != "No Marker"
        num_sample_lanes = min(len(self.fragments), 8)  # Max 8 sample lanes
        num_lanes = num_sample_lanes + (1 if has_marker else 0)
        
        if num_lanes == 0:
            return
            
        lane_spacing = gel_width / num_lanes
        well_y = 40
        
        # Draw wells
        for i in range(num_lanes):
            well_x = 50 + (i + 0.5) * lane_spacing - well_width/2
            well = self.scene.addRect(well_x, well_y, well_width, well_height,
                                    QPen(Qt.black), QBrush(Qt.black))
        
        # Get gel parameters
        gel_conc = self.gel_concentration.currentText()
        voltage = self.voltage.value()
        time = self.run_time.value()
        
        # Draw marker lane (if selected)
        lane_index = 0
        if has_marker:
            marker_fragments = self.get_marker_fragments(marker_type)
            lane_x = 50 + (lane_index + 0.5) * lane_spacing
            
            for fragment_size in marker_fragments:
                distance = self.calculate_migration_distance(fragment_size, gel_conc, voltage, time)
                band_y = 50 + 20 + distance
                
                # Skip bands that would run off the gel
                if band_y > 50 + gel_height - 10:
                    continue
                
                # Draw band with realistic appearance
                band_width = 25
                band_height = 2 if fragment_size > 5000 else 3
                
                band = self.scene.addRect(lane_x - band_width/2, band_y, band_width, band_height,
                                        QPen(Qt.darkBlue), QBrush(QColor(0, 100, 255)))
                
                # Add size label (only for major bands)
                if fragment_size in [10000, 5000, 3000, 1500, 1000, 500, 250] or fragment_size > 15000:
                    if fragment_size >= 1000:
                        label_text = f"{fragment_size/1000:.1f}kb"
                    else:
                        label_text = f"{fragment_size}bp"
                    
                    label = self.scene.addText(label_text, QFont("Arial", 7))
                    label.setPos(lane_x + 18, band_y - 3)
            
            # Lane label
            lane_label = self.scene.addText("Marker", QFont("Arial", 9, QFont.Bold))
            lane_label.setPos(lane_x - 20, 15)
            lane_index += 1
        
        # Draw sample lanes
        # Group fragments by size for better lane organization
        fragment_groups = self.group_similar_fragments(self.fragments)
        
        for group_idx, fragment_group in enumerate(fragment_groups[:num_sample_lanes]):
            lane_x = 50 + (lane_index + 0.5) * lane_spacing
            
            # Draw all fragments in this group
            for frag_idx, fragment in enumerate(fragment_group):
                fragment_size = len(fragment)
                distance = self.calculate_migration_distance(fragment_size, gel_conc, voltage, time)
                band_y = 50 + 20 + distance
                
                # Skip bands that would run off the gel
                if band_y > 50 + gel_height - 10:
                    continue
                
                # Band appearance based on fragment properties
                band_width = 25
                band_height = 3
                
                # Intensity based on fragment size (larger = brighter)
                base_intensity = min(255, max(80, 120 + fragment_size // 100))
                
                # Multiple fragments of same size = brighter band
                if len(fragment_group) > 1:
                    base_intensity = min(255, base_intensity + 30)
                
                band_color = QColor(base_intensity, 0, base_intensity // 2)
                
                # Slight horizontal offset for multiple fragments
                x_offset = (frag_idx - len(fragment_group)/2) * 2
                
                band = self.scene.addRect(lane_x - band_width/2 + x_offset, band_y, 
                                        band_width, band_height,
                                        QPen(band_color), QBrush(band_color))
            
            # Add size label for the group
            representative_size = len(fragment_group[0])
            distance = self.calculate_migration_distance(representative_size, gel_conc, voltage, time)
            band_y = 50 + 20 + distance
            
            if band_y <= 50 + gel_height - 10:  # Only if visible
                if representative_size >= 1000:
                    label_text = f"{representative_size/1000:.1f}kb"
                else:
                    label_text = f"{representative_size}bp"
                
                if len(fragment_group) > 1:
                    label_text += f" (×{len(fragment_group)})"
                
                label = self.scene.addText(label_text, QFont("Arial", 7))
                label.setPos(lane_x + 18, band_y - 3)
            
            # Lane label
            lane_label = self.scene.addText(f"Sample {group_idx+1}", QFont("Arial", 9))
            lane_label.setPos(lane_x - 25, 15)
            lane_index += 1
        
        # Add gel information
        gel_info = f"Agarose: {gel_conc}, {voltage}V, {time}min"
        info_label = self.scene.addText(gel_info, QFont("Arial", 10))
        info_label.setPos(50, gel_height + 70)
        
        # Add direction indicator with electrode labels
        direction_label = self.scene.addText("(-) Cathode ← Migration → Anode (+)", QFont("Arial", 9))
        direction_label.setPos(150, gel_height + 90)
        
        # Add subtle migration direction indicators
        self.add_migration_arrows(gel_width, gel_height)
        
        # Add scale reference
        scale_label = self.scene.addText(f"Migration distance scale: 0-{int(gel_height*0.8)} pixels", QFont("Arial", 8))
        scale_label.setPos(50, gel_height + 110)
    
    def group_similar_fragments(self, fragments):
        """Group fragments of similar sizes together"""
        if not fragments:
            return []
        
        # Sort fragments by size
        sorted_fragments = sorted(fragments, key=len, reverse=True)
        
        groups = []
        current_group = [sorted_fragments[0]]
        
        for fragment in sorted_fragments[1:]:
            # Group fragments within 10% size difference
            current_size = len(fragment)
            group_size = len(current_group[0])
            
            size_diff = abs(current_size - group_size) / group_size
            
            if size_diff < 0.1:  # Within 10%
                current_group.append(fragment)
            else:
                groups.append(current_group)
                current_group = [fragment]
        
        groups.append(current_group)
        return groups
    
    def add_migration_arrows(self, gel_width, gel_height):
        """Add subtle migration direction indicators"""
        try:
            from PyQt5.QtGui import QPolygonF
            from PyQt5.QtCore import QPointF
            
            # Add subtle side indicators
            arrow_color = QColor(120, 120, 120, 80)  # Very subtle gray
            
            # Left side - just a few small arrows
            for y in [100, 200, 300]:
                if y < gel_height + 50:
                    arrow_x = 25
                    arrow_y = y
                    
                    # Small downward chevron
                    line_length = 4
                    self.scene.addLine(arrow_x - 2, arrow_y, arrow_x, arrow_y + line_length, 
                                     QPen(arrow_color, 1))
                    self.scene.addLine(arrow_x + 2, arrow_y, arrow_x, arrow_y + line_length, 
                                     QPen(arrow_color, 1))
            
            # Right side - matching arrows
            for y in [100, 200, 300]:
                if y < gel_height + 50:
                    arrow_x = 50 + gel_width + 25
                    arrow_y = y
                    
                    # Small downward chevron
                    line_length = 4
                    self.scene.addLine(arrow_x - 2, arrow_y, arrow_x, arrow_y + line_length, 
                                     QPen(arrow_color, 1))
                    self.scene.addLine(arrow_x + 2, arrow_y, arrow_x, arrow_y + line_length, 
                                     QPen(arrow_color, 1))
            
        except Exception as e:
            print(f"Error adding migration arrows: {e}")
    
    def update_gel(self):
        """Update the gel visualization"""
        self.draw_gel()
    
    def update_fragment_info(self):
        """Update fragment information display"""
        if not self.fragments:
            self.fragment_info.setText("No fragments to display.")
            return
        
        # Sort fragments by size
        sorted_fragments = sorted(self.fragments, key=len, reverse=True)
        
        info_text = f"Total fragments: {len(self.fragments)}\n\n"
        info_text += "Fragment sizes (bp):\n"
        
        for i, fragment in enumerate(sorted_fragments, 1):
            size = len(fragment)
            if size >= 1000:
                size_str = f"{size/1000:.1f} kb"
            else:
                size_str = f"{size} bp"
            
            info_text += f"{i:2d}. {size_str:>8s} ({size:,} bp)\n"
        
        # Add statistics
        sizes = [len(f) for f in self.fragments]
        if sizes:
            info_text += f"\nStatistics:\n"
            info_text += f"Largest:  {max(sizes):,} bp\n"
            info_text += f"Smallest: {min(sizes):,} bp\n"
            info_text += f"Average:  {sum(sizes)//len(sizes):,} bp\n"
            info_text += f"Total:    {sum(sizes):,} bp"
        
        self.fragment_info.setText(info_text)
