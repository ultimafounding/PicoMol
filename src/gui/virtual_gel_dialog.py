from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QGraphicsView, QGraphicsScene, 
                             QHBoxLayout, QLabel, QComboBox, QPushButton, QGroupBox,
                             QSpinBox, QFormLayout, QTextEdit, QSplitter)
from PyQt5.QtGui import QPen, QBrush, QColor, QFont, QLinearGradient
from PyQt5.QtCore import Qt
import math
from collections import defaultdict

class VirtualGelDialog(QDialog):
    def __init__(self, fragments, parent=None, enzyme_fragments=None):
        super().__init__(parent)
        self.setWindowTitle("Virtual Agarose Gel Electrophoresis")
        self.setMinimumSize(700, 600)
        self.fragments = fragments
        self.enzyme_fragments = enzyme_fragments
        
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
        
        # Migration scaling control
        self.migration_scaling = QSpinBox()
        self.migration_scaling.setRange(50, 200)
        self.migration_scaling.setValue(100)
        self.migration_scaling.setSuffix("%")
        self.migration_scaling.setToolTip("Adjust migration distance scaling for better separation")
        self.migration_scaling.valueChanged.connect(self.update_gel)
        controls_layout.addRow("Migration scaling:", self.migration_scaling)
        
        # Image size scaling control
        self.image_scale = QSpinBox()
        self.image_scale.setRange(25, 300)
        self.image_scale.setValue(100)
        self.image_scale.setSuffix("%")
        self.image_scale.setToolTip("Adjust the overall size of the gel image")
        self.image_scale.valueChanged.connect(self.update_gel)
        controls_layout.addRow("Image size:", self.image_scale)
        
        # Text labels control
        from PyQt5.QtWidgets import QCheckBox
        self.show_labels = QCheckBox("Show fragment labels")
        self.show_labels.setChecked(True)
        self.show_labels.setToolTip("Show/hide fragment size labels on bands")
        self.show_labels.toggled.connect(self.update_gel)
        controls_layout.addRow(self.show_labels)
        
        # Label density control
        self.label_density = QComboBox()
        self.label_density.addItems(["Minimal", "Moderate", "All visible"])
        self.label_density.setCurrentText("Moderate")
        self.label_density.setToolTip("Control how many labels are shown to reduce clutter")
        self.label_density.currentTextChanged.connect(self.update_gel)
        controls_layout.addRow("Label density:", self.label_density)
        
        # Label style control
        self.label_style = QComboBox()
        self.label_style.addItems(["Clean (recommended)", "Detailed", "Minimal"])
        self.label_style.setCurrentText("Clean (recommended)")
        self.label_style.setToolTip("Choose label style for best appearance")
        self.label_style.currentTextChanged.connect(self.update_gel)
        controls_layout.addRow("Label style:", self.label_style)
        
        # Marker selection
        self.marker_combo = QComboBox()
        self.marker_combo.addItems([
            "1kb DNA Ladder", 
            "100bp DNA Ladder", 
            "GeneRuler 1kb",
            "GeneRuler 100bp Plus",
            "NEB 1kb Plus",
            "Quick-Load 1kb",
            "Quick-Load 100bp",
            "λ DNA/HindIII",
            "No Marker"
        ])
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
        
        self.update_fragment_info()
        self.draw_gel()

    def get_marker_fragments(self, marker_type):
        """Get standard DNA marker fragment sizes"""
        markers = {
            "1kb DNA Ladder": [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250],
            "100bp DNA Ladder": [1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
            "λ DNA/HindIII": [23130, 9416, 6557, 4361, 2322, 2027, 564, 125],
            "GeneRuler 1kb": [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 750, 500, 250],
            "GeneRuler 100bp Plus": [3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
            "NEB 1kb Plus": [12000, 11000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 650, 400, 300, 200, 100],
            "Quick-Load 1kb": [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 500],
            "Quick-Load 100bp": [1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
            "No Marker": []
        }
        return markers.get(marker_type, [])
    
    def calculate_migration_distance(self, fragment_size, gel_concentration, voltage, time, all_fragment_sizes=None):
        """Calculate migration distance based on fragment size and gel parameters with dynamic scaling"""
        # Handle edge cases
        if fragment_size <= 0:
            return 0
        
        # Get gel height for scaling (use base height, scaling applied later)
        base_gel_height = 350
        usable_distance = base_gel_height * 0.8  # Use 80% of gel height for migration
        
        # Dynamic scaling based on actual fragment sizes present
        if all_fragment_sizes and len(all_fragment_sizes) > 1:
            # Use actual size range for optimal scaling
            min_size = min(all_fragment_sizes)
            max_size = max(all_fragment_sizes)
            
            # Ensure reasonable range for logarithmic scaling
            if max_size / min_size < 2:  # If range is too small, expand it
                center = (min_size + max_size) / 2
                min_size = max(50, center * 0.5)
                max_size = center * 2
        else:
            # Default range for single fragments or fallback
            min_size = max(50, fragment_size * 0.1)
            max_size = fragment_size * 10
        
        # Logarithmic scaling for better separation
        log_fragment = math.log10(max(fragment_size, 50))
        log_min = math.log10(max(min_size, 50))
        log_max = math.log10(max_size)
        
        # Prevent division by zero
        if log_max == log_min:
            size_factor = 0.5
        else:
            # Invert so smaller fragments migrate further (0 = largest, 1 = smallest)
            size_factor = max(0, min(1, (log_max - log_fragment) / (log_max - log_min)))
        
        # Gel concentration factor (higher concentration = less migration)
        concentration = float(gel_concentration.rstrip('%'))
        conc_factor = 1.5 / concentration  # 1.5% is reference
        
        # Voltage and time factor (linear relationship)
        vt_factor = (voltage / 100) * (time / 60)
        
        # Apply molecular sieving effects for very large or small fragments
        if fragment_size > 20000:
            size_factor *= 0.4  # Very large fragments barely migrate
        elif fragment_size > 10000:
            size_factor *= 0.7  # Large fragments migrate less
        elif fragment_size < 100:
            # Small fragments migrate more but may run off gel
            size_factor = min(1.0, size_factor * 1.3)
        
        # Get migration scaling factor from UI control
        migration_scaling = getattr(self, 'migration_scaling', None)
        if migration_scaling:
            scale_multiplier = migration_scaling.value() / 100.0
        else:
            scale_multiplier = 1.0
        
        # Calculate final distance (deterministic, no random variation)
        distance = usable_distance * size_factor * conc_factor * vt_factor * 0.6 * scale_multiplier
        
        # Ensure minimum separation between different sized fragments
        if all_fragment_sizes and len(set(all_fragment_sizes)) > 1:
            # Add small offset based on fragment size for better visual separation
            size_offset = (fragment_size % 100) * 0.1  # Small deterministic offset
            distance += size_offset
        
        return min(distance, usable_distance * 0.95)  # Cap at 95% of usable distance
    
    def get_all_fragment_sizes(self):
        """Get all fragment sizes for dynamic scaling"""
        all_sizes = []
        
        # Collect sizes from enzyme fragments
        if self.enzyme_fragments:
            for fragments in self.enzyme_fragments.values():
                all_sizes.extend([len(f) for f in fragments])
        
        # Collect sizes from regular fragments
        if self.fragments:
            all_sizes.extend([len(f) for f in self.fragments])
        
        # Add marker sizes if present
        marker_type = self.marker_combo.currentText()
        if marker_type != "No Marker":
            marker_fragments = self.get_marker_fragments(marker_type)
            all_sizes.extend(marker_fragments)
        
        return list(set(all_sizes))  # Remove duplicates
    
    def draw_gel(self):
        """Draw the virtual gel with fragments and markers"""
        self.scene.clear()
        
        # Initialize simple text layout manager
        self.used_label_positions = []  # Simple list of used Y positions for labels
        self.label_side_counter = 0  # Counter for alternating label sides
        
        if not self.fragments:
            # Show empty gel with instructions
            self.scene.addText("No fragments to display.\nLoad sequences and perform restriction digest.", 
                             QFont("Arial", 12))
            return
        
        # Get all fragment sizes for dynamic scaling
        all_fragment_sizes = self.get_all_fragment_sizes()
        
        # Calculate lanes first
        marker_type = self.marker_combo.currentText()
        has_marker = marker_type != "No Marker"
        # Calculate number of sample lanes based on available data
        if self.enzyme_fragments:
            # Prioritize important combinations and limit total lanes
            total_combinations = len(self.enzyme_fragments)
            if total_combinations <= 8:
                num_sample_lanes = total_combinations
            else:
                # Show most important combinations: uncut, singles, doubles, complete
                num_sample_lanes = min(12, total_combinations)  # Allow more lanes for combinations
        else:
            # Fallback to fragment-based calculation
            num_sample_lanes = min(len(self.fragments), 8)
        num_lanes = num_sample_lanes + (1 if has_marker else 0)
        
        # Get image scaling factor
        image_scale = getattr(self, 'image_scale', None)
        if image_scale:
            image_scale_factor = image_scale.value() / 100.0
        else:
            image_scale_factor = 1.0
        
        # Gel parameters (now that we know num_lanes) - scaled by image scale factor
        base_gel_width = max(400, num_lanes * 50)
        base_gel_height = 350
        
        gel_width = int(base_gel_width * image_scale_factor)
        gel_height = int(base_gel_height * image_scale_factor)
        well_width = int(30 * image_scale_factor)
        well_height = int(15 * image_scale_factor)
        
        # Scale all positions and sizes
        margin = int(50 * image_scale_factor)
        buffer_height = int(20 * image_scale_factor)
        pen_width = max(1, int(2 * image_scale_factor))
        
        # Draw gel background with gradient
        gel_rect = self.scene.addRect(margin, margin, gel_width, gel_height, 
                                     QPen(Qt.black, pen_width), QBrush(QColor(240, 245, 255)))
        
        # Draw buffer chambers
        # Top buffer (negative electrode)
        top_buffer_y = margin - buffer_height - int(10 * image_scale_factor)
        top_buffer = self.scene.addRect(margin, top_buffer_y, gel_width, buffer_height,
                                       QPen(Qt.darkGray), QBrush(QColor(200, 200, 255)))
        # Bottom buffer (positive electrode)
        bottom_buffer_y = margin + gel_height + int(10 * image_scale_factor)
        bottom_buffer = self.scene.addRect(margin, bottom_buffer_y, gel_width, buffer_height,
                                          QPen(Qt.darkGray), QBrush(QColor(255, 200, 200)))
        
        if num_lanes == 0:
            return
            
        lane_spacing = gel_width / num_lanes
        well_y = int(40 * image_scale_factor)
        
        # Draw wells
        for i in range(num_lanes):
            well_x = margin + (i + 0.5) * lane_spacing - well_width/2
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
            lane_x = margin + (lane_index + 0.5) * lane_spacing
            
            for fragment_size in marker_fragments:
                distance = self.calculate_migration_distance(fragment_size, gel_conc, voltage, time, all_fragment_sizes)
                # Scale the distance and position
                scaled_distance = distance * image_scale_factor
                band_y = margin + int(20 * image_scale_factor) + scaled_distance
                
                # Skip bands that would run off the gel
                if band_y > margin + gel_height - int(10 * image_scale_factor):
                    continue
                
                # Draw band with realistic appearance (scaled)
                band_width = int(25 * image_scale_factor)
                band_height = int((2 if fragment_size > 5000 else 3) * image_scale_factor)
                
                band = self.scene.addRect(lane_x - band_width/2, band_y, band_width, band_height,
                                        QPen(Qt.darkBlue), QBrush(QColor(0, 100, 255)))
                
                # Add clean marker labels
                if self.show_labels.isChecked() and self.should_show_marker_label(fragment_size):
                    self.add_clean_label(fragment_size, lane_x, band_y, image_scale_factor, is_marker=True)
            
            # Lane label
            font_size = max(7, int(9 * image_scale_factor))
            lane_label = self.scene.addText("Marker", QFont("Arial", font_size, QFont.Bold))
            lane_label.setPos(lane_x - int(20 * image_scale_factor), int(15 * image_scale_factor))
            lane_index += 1
        
        # Draw sample lanes
        # Use enzyme-specific lanes if available, otherwise group by size
        if self.enzyme_fragments:
            # Prioritize and organize enzyme combinations for display
            enzyme_items = list(self.enzyme_fragments.items())
            
            # Sort combinations by priority: uncut, singles, doubles, triples, complete
            def get_priority(item):
                name, fragments = item
                if name == "Uncut":
                    return (0, name)  # Highest priority
                elif " + " not in name:
                    return (1, name)  # Single enzymes
                elif name.startswith("Complete digest"):
                    return (9, name)  # Complete digest at end
                else:
                    # Count number of enzymes in combination
                    enzyme_count = name.count(" + ") + 1
                    return (enzyme_count + 1, name)
            
            sorted_items = sorted(enzyme_items, key=get_priority)
            
            # Limit to displayable number of lanes
            displayed_items = sorted_items[:num_sample_lanes]
            fragment_groups = [fragments for enzyme, fragments in displayed_items]
            enzyme_names = [enzyme for enzyme, fragments in displayed_items]
        else:
            # Fallback to size-based grouping
            fragment_groups = self.group_similar_fragments(self.fragments)
            enzyme_names = None
        
        for group_idx, fragment_group in enumerate(fragment_groups[:num_sample_lanes]):
            lane_x = margin + (lane_index + 0.5) * lane_spacing
            
            # Draw all fragments in this group
            # Only draw fragments if the enzyme cuts
            if fragment_group:
                for frag_idx, fragment in enumerate(fragment_group):
                    fragment_size = len(fragment)
                    distance = self.calculate_migration_distance(fragment_size, gel_conc, voltage, time, all_fragment_sizes)
                    # Scale the distance and position
                    scaled_distance = distance * image_scale_factor
                    band_y = margin + int(20 * image_scale_factor) + scaled_distance
                    
                    # Skip bands that would run off the gel
                    if band_y > margin + gel_height - int(10 * image_scale_factor):
                        continue
                    
                    # Band appearance based on fragment properties (scaled)
                    band_width = int(25 * image_scale_factor)
                    band_height = int(3 * image_scale_factor)
                    
                    # Intensity based on fragment size (larger = brighter)
                    base_intensity = min(255, max(80, 120 + fragment_size // 100))
                    
                    # Multiple fragments of same size = brighter band
                    if len(fragment_group) > 1:
                        base_intensity = min(255, base_intensity + 30)
                    
                    band_color = QColor(base_intensity, 0, base_intensity // 2)
                    
                    # Slight horizontal offset for multiple fragments (scaled)
                    x_offset = (frag_idx - len(fragment_group)/2) * int(2 * image_scale_factor)
                    
                    band = self.scene.addRect(lane_x - band_width/2 + x_offset, band_y, 
                                            band_width, band_height,
                                            QPen(band_color), QBrush(band_color))
                
                # Add size label for the group
                representative_size = len(fragment_group[0])
                distance = self.calculate_migration_distance(representative_size, gel_conc, voltage, time, all_fragment_sizes)
                # Scale the distance and position
                scaled_distance = distance * image_scale_factor
                band_y = margin + int(20 * image_scale_factor) + scaled_distance
                
                # Add clean fragment labels
                if (self.show_labels.isChecked() and 
                    band_y <= margin + gel_height - int(10 * image_scale_factor) and
                    self.should_show_fragment_label(representative_size, group_idx)):
                    
                    # Add multiplicity info if needed
                    multiplicity = len(fragment_group) if len(fragment_group) > 1 else None
                    self.add_clean_label(representative_size, lane_x, band_y, image_scale_factor, 
                                       is_marker=False, multiplicity=multiplicity)
            
            # Lane label
            # Create lane name - use enzyme name if available
            if enzyme_names:
                lane_name = str(enzyme_names[group_idx])
                # Truncate very long combination names
                if len(lane_name) > 20:
                    # For long names, show abbreviated version
                    if "Complete digest" in lane_name:
                        lane_name = "Complete"
                    elif " + " in lane_name:
                        # Show first and last enzyme with count
                        parts = lane_name.split(" + ")
                        if len(parts) > 2:
                            lane_name = f"{parts[0]} + ... + {parts[-1]} ({len(parts)})"
            else:
                lane_name = f"Sample {group_idx+1}"
            
            # Use smaller font for long names (scaled)
            base_font_size = 7 if len(lane_name) > 15 else 9
            font_size = max(6, int(base_font_size * image_scale_factor))
            
            # Position lane labels with better spacing
            label_x = lane_x - int(30 * image_scale_factor)
            label_y = int(15 * image_scale_factor)
            
            # For lane labels, we don't need overlap detection as they're at fixed positions
            lane_label = self.scene.addText(lane_name, QFont("Arial", font_size))
            lane_label.setPos(label_x, label_y)
            lane_index += 1
        
        # Add gel information (scaled)
        gel_info = f"Agarose: {gel_conc}, {voltage}V, {time}min"
        info_font_size = max(8, int(10 * image_scale_factor))
        info_label = self.scene.addText(gel_info, QFont("Arial", info_font_size))
        info_label.setPos(margin, margin + gel_height + int(70 * image_scale_factor))
        
        # Add direction indicator with electrode labels (scaled)
        direction_font_size = max(7, int(9 * image_scale_factor))
        direction_label = self.scene.addText("(-) Cathode ← Migration → Anode (+)", QFont("Arial", direction_font_size))
        direction_label.setPos(margin + int(100 * image_scale_factor), margin + gel_height + int(90 * image_scale_factor))
        
        # Add subtle migration direction indicators
        self.add_migration_arrows(gel_width, gel_height, image_scale_factor, margin)
        
        # Add scale reference (scaled)
        scale_font_size = max(6, int(8 * image_scale_factor))
        scale_label = self.scene.addText(f"Image scale: {int(image_scale_factor*100)}%", QFont("Arial", scale_font_size))
        scale_label.setPos(margin, margin + gel_height + int(110 * image_scale_factor))
    
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
    
    def add_migration_arrows(self, gel_width, gel_height, image_scale_factor=1.0, margin=50):
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
        
        # Update the view to show the scaled content
        if hasattr(self, 'scene') and hasattr(self, 'view'):
            # Get the scene bounds and add some padding
            scene_rect = self.scene.itemsBoundingRect()
            padding = 20
            scene_rect.adjust(-padding, -padding, padding, padding)
            self.scene.setSceneRect(scene_rect)
            self.view.fitInView(scene_rect, Qt.KeepAspectRatio)
    
    def update_fragment_info(self):
        """Update fragment information display"""
        if not self.fragments and not self.enzyme_fragments:
            self.fragment_info.setText("No fragments to display.")
            return
        
        info_text = ""
        
        if self.enzyme_fragments:
            # Show information about enzyme combinations
            info_text += f"Enzyme Combinations: {len(self.enzyme_fragments)}\n\n"
            
            # Categorize combinations
            uncut = []
            singles = []
            doubles = []
            triples = []
            complete = []
            
            for name, fragments in self.enzyme_fragments.items():
                if name == "Uncut":
                    uncut.append((name, fragments))
                elif name.startswith("Complete digest"):
                    complete.append((name, fragments))
                elif " + " not in name:
                    singles.append((name, fragments))
                elif name.count(" + ") == 1:
                    doubles.append((name, fragments))
                else:
                    triples.append((name, fragments))
            
            # Display summary by category
            if uncut:
                info_text += f"Uncut: {len(uncut[0][1])} fragment(s)\n"
            if singles:
                info_text += f"Single digests: {len(singles)}\n"
            if doubles:
                info_text += f"Double digests: {len(doubles)}\n"
            if triples:
                info_text += f"Triple+ digests: {len(triples)}\n"
            if complete:
                info_text += f"Complete digest: {len(complete[0][1])} fragment(s)\n"
            
            info_text += "\nFragment size ranges:\n"
            
            # Show size ranges for each category
            for category_name, category_list in [("Uncut", uncut), ("Singles", singles), 
                                                ("Doubles", doubles), ("Triples+", triples), 
                                                ("Complete", complete)]:
                if category_list:
                    all_sizes = []
                    for _, fragments in category_list:
                        all_sizes.extend([len(f) for f in fragments])
                    
                    if all_sizes:
                        min_size = min(all_sizes)
                        max_size = max(all_sizes)
                        avg_size = sum(all_sizes) // len(all_sizes)
                        
                        info_text += f"{category_name}: {min_size}-{max_size} bp (avg: {avg_size} bp)\n"
        
        else:
            # Fallback to original fragment display
            sorted_fragments = sorted(self.fragments, key=len, reverse=True)
            
            info_text = f"Total fragments: {len(self.fragments)}\n\n"
            info_text += "Fragment sizes (bp):\n"
            
            for i, fragment in enumerate(sorted_fragments[:10], 1):  # Limit to first 10
                size = len(fragment)
                if size >= 1000:
                    size_str = f"{size/1000:.1f} kb"
                else:
                    size_str = f"{size} bp"
                
                info_text += f"{i:2d}. {size_str:>8s} ({size:,} bp)\n"
            
            if len(sorted_fragments) > 10:
                info_text += f"... and {len(sorted_fragments) - 10} more\n"
            
            # Add statistics
            sizes = [len(f) for f in self.fragments]
            if sizes:
                info_text += f"\nStatistics:\n"
                info_text += f"Largest:  {max(sizes):,} bp\n"
                info_text += f"Smallest: {min(sizes):,} bp\n"
                info_text += f"Average:  {sum(sizes)//len(sizes):,} bp\n"
                info_text += f"Total:    {sum(sizes):,} bp"
        
        self.fragment_info.setText(info_text)
    
    def find_available_y_position(self, preferred_y, font_size):
        """Find the next available Y position that doesn't overlap with existing labels"""
        text_height = font_size * 1.5  # Height including padding
        min_spacing = text_height + 8  # Increased minimum space between labels
        
        # Check if preferred position is available
        if not self.position_conflicts(preferred_y, min_spacing):
            return preferred_y
        
        # Try positions above and below the preferred position
        for offset in range(1, 20):  # Try up to 20 positions
            # Try above
            test_y_above = preferred_y - (offset * min_spacing)
            if not self.position_conflicts(test_y_above, min_spacing):
                return test_y_above
            
            # Try below
            test_y_below = preferred_y + (offset * min_spacing)
            if not self.position_conflicts(test_y_below, min_spacing):
                return test_y_below
        
        # If no position found, use preferred position anyway
        return preferred_y
    
    def position_conflicts(self, y_position, min_spacing):
        """Check if a Y position conflicts with existing labels"""
        for used_y in self.used_label_positions:
            if abs(y_position - used_y) < min_spacing:
                return True
        return False
    
    def should_show_marker_label(self, fragment_size):
        """Determine if a marker fragment should have a label based on density setting"""
        density = self.label_density.currentText()
        style = self.label_style.currentText()
        
        # Adjust based on style - minimal style shows fewer labels
        if style == "Minimal":
            if density == "Minimal":
                return fragment_size in [10000, 5000, 1000, 500]
            elif density == "Moderate":
                return fragment_size in [10000, 5000, 3000, 1000, 500]
            else:
                return fragment_size in [10000, 8000, 5000, 3000, 1000, 500, 250]
        else:
            # Standard logic for other styles
            if density == "Minimal":
                return fragment_size in [10000, 5000, 3000, 1000, 500]
            elif density == "Moderate":
                return (fragment_size in [10000, 8000, 5000, 3000, 2000, 1500, 1000, 500, 250] or 
                       fragment_size > 15000)
            else:  # "All visible"
                return (fragment_size in [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250] or 
                       fragment_size > 15000)
    
    def should_show_fragment_label(self, fragment_size, group_index):
        """Determine if a fragment should have a label based on density setting"""
        density = self.label_density.currentText()
        style = self.label_style.currentText()
        
        # Adjust based on style
        if style == "Minimal":
            if density == "Minimal":
                return group_index % 4 == 0 and fragment_size > 1000
            elif density == "Moderate":
                return group_index % 3 == 0 and fragment_size > 500
            else:
                return group_index % 2 == 0
        else:
            # Standard logic for other styles
            if density == "Minimal":
                return (group_index % 3 == 0 and 
                       (fragment_size > 2000 or fragment_size in [1000, 1500, 500]))
            elif density == "Moderate":
                return (group_index % 2 == 0 or fragment_size > 5000 or 
                       fragment_size in [1000, 1500, 2000, 3000, 500])
            else:  # "All visible"
                return True
    
    def add_clean_label(self, fragment_size, lane_x, band_y, image_scale_factor, is_marker=False, multiplicity=None):
        """Add a clean, well-positioned label for a fragment"""
        style = self.label_style.currentText()
        
        # Format the text based on style
        if style == "Minimal":
            if fragment_size >= 1000:
                label_text = f"{fragment_size/1000:.0f}k"
            else:
                label_text = f"{fragment_size}"
        elif style == "Clean (recommended)":
            if fragment_size >= 1000:
                label_text = f"{fragment_size/1000:.1f}kb"
            else:
                label_text = f"{fragment_size}bp"
        else:  # "Detailed"
            if fragment_size >= 1000:
                label_text = f"{fragment_size/1000:.1f}kb ({fragment_size:,}bp)"
            else:
                label_text = f"{fragment_size}bp"
            
            if multiplicity and multiplicity > 1:
                label_text += f" ×{multiplicity}"
        
        # Only add multiplicity for clean style if it's significant
        if style == "Clean (recommended)" and multiplicity and multiplicity > 1:
            label_text += f" ×{multiplicity}"
        
        # Choose font size based on style and scale
        if style == "Minimal":
            font_size = max(5, int(6 * image_scale_factor))
        else:
            font_size = max(6, int(8 * image_scale_factor))
        
        # Position the label cleanly
        if is_marker:
            # Markers go on the right side, well spaced
            label_x = lane_x + int(35 * image_scale_factor)
            label_y = band_y + int(2 * image_scale_factor)  # Slightly below band center
        else:
            # Fragment labels alternate sides for better distribution
            if self.label_side_counter % 2 == 0:
                label_x = lane_x + int(35 * image_scale_factor)  # Right side
            else:
                label_x = lane_x - int(80 * image_scale_factor)  # Left side with more space
            label_y = band_y + int(2 * image_scale_factor)
        
        # Find available Y position to avoid overlap
        final_y = self.find_available_y_position(label_y, font_size)
        
        # Create the label with better styling
        label = self.scene.addText(label_text, QFont("Arial", font_size))
        
        # Style the label based on type
        if is_marker:
            label.setDefaultTextColor(QColor(0, 0, 150))  # Dark blue for markers
        else:
            label.setDefaultTextColor(QColor(80, 80, 80))  # Dark gray for fragments
        
        label.setPos(label_x, final_y)
        
        # Record position and increment counter
        self.used_label_positions.append(final_y)
        if not is_marker:
            self.label_side_counter += 1
