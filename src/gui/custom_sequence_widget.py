from PyQt5.QtWidgets import QWidget, QScrollArea, QVBoxLayout
from PyQt5.QtGui import QPainter, QPen, QBrush, QFont, QColor, QFontMetrics
from PyQt5.QtCore import Qt, QRect, QPoint, QSize, pyqtSignal
from Bio.Seq import Seq
import math

class CustomSequenceWidget(QWidget):
    """Custom painted sequence widget using QPainter for SnapGene-style display"""
    
    positionChanged = pyqtSignal(int)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.record = None
        self.bases_per_line = 60
        self.show_features = True
        self.show_translation = True
        self.show_complement = False
        self.show_enzymes = False
        self.enzyme_sites = {}
        
        # Visual settings
        self.base_font = QFont("Consolas", 11)
        self.ruler_font = QFont("Consolas", 9)
        self.number_font = QFont("Consolas", 10)
        
        # Measurements - use precise font metrics
        fm = QFontMetrics(self.base_font)
        self.char_width = fm.horizontalAdvance("A")  # Consolas is monospace
        self.line_height = fm.height() + 6  # Extra spacing for readability
        self.margin = 20
        self.number_width = 80  # More space for line numbers
        self.group_spacing = self.char_width // 2  # Space between 10-base groups
        
        # Colors - SnapGene style
        self.colors = {
            'A': QColor(255, 0, 0),      # Red
            'T': QColor(0, 0, 255),      # Blue  
            'C': QColor(0, 150, 0),      # Green
            'G': QColor(255, 140, 0),    # Orange
            'background': QColor(255, 255, 255),
            'ruler': QColor(120, 120, 120),
            'numbers': QColor(80, 80, 80),
            'translation': QColor(0, 0, 150),
            'complement': QColor(150, 0, 150),
            'enzyme': QColor(200, 0, 0)
        }
        
        # Feature colors - SnapGene style
        self.feature_colors = {
            "gene": QColor(255, 255, 0, 80),        # Yellow
            "CDS": QColor(255, 165, 0, 80),         # Orange  
            "promoter": QColor(255, 0, 0, 80),      # Red
            "terminator": QColor(128, 0, 128, 80),  # Purple
            "rep_origin": QColor(0, 255, 0, 80),    # Green
            "misc_feature": QColor(192, 192, 192, 80), # Gray
            "regulatory": QColor(255, 192, 203, 80), # Pink
            "enhancer": QColor(255, 255, 224, 80),   # Light yellow
            "primer_bind": QColor(0, 255, 255, 80),  # Cyan
            "protein_bind": QColor(144, 238, 144, 80) # Light green
        }
        
        self.setMinimumSize(800, 400)
        self.setMouseTracking(True)
        
        # Current position tracking
        self.current_position = 0
    
    def get_x_position(self, base_index):
        """Calculate precise x position for a base at given index"""
        x = self.number_width + base_index * self.char_width
        # Add spacing for every complete group of 10 bases
        if base_index > 0:
            x += (base_index // 10) * self.group_spacing
        return x
    
    def set_sequence(self, record):
        """Set the sequence to display"""
        self.record = record
        self.update_size()
        self.update()
    
    def update_size(self):
        """Update widget size based on sequence length"""
        if not self.record:
            self.setMinimumSize(800, 400)
            return
        
        sequence = str(self.record.seq)
        num_lines = (len(sequence) + self.bases_per_line - 1) // self.bases_per_line
        
        # Calculate lines per block
        lines_per_block = 2  # ruler + sequence
        if self.show_features:
            lines_per_block += 1
        if self.show_translation:
            lines_per_block += 1
        if self.show_complement:
            lines_per_block += 1
        if self.show_enzymes:
            lines_per_block += 1
        lines_per_block += 1  # spacing
        
        height = num_lines * lines_per_block * self.line_height + self.margin * 2 + 50
        
        # Calculate width using precise positioning
        if self.bases_per_line > 0:
            last_base_x = self.get_x_position(self.bases_per_line - 1)
            width = last_base_x + self.char_width + self.margin * 2
        else:
            width = self.number_width + self.margin * 2
        
        self.setMinimumSize(width, height)
        self.resize(width, height)
    
    def paintEvent(self, event):
        """Custom paint event"""
        if not self.record:
            return
        
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        
        # Fill background
        painter.fillRect(self.rect(), self.colors['background'])
        
        sequence = str(self.record.seq).upper()
        y_pos = self.margin
        
        # Draw header
        painter.setFont(QFont("Consolas", 12, QFont.Bold))
        painter.setPen(QPen(Qt.black))
        painter.drawText(self.margin, y_pos, f"Sequence: {self.record.id} ({len(sequence)} bp)")
        y_pos += self.line_height * 2
        
        # Draw sequence in blocks
        for i in range(0, len(sequence), self.bases_per_line):
            line_start = i + 1
            line_end = min(i + self.bases_per_line, len(sequence))
            line_seq = sequence[i:line_end]
            
            block_start_y = y_pos
            
            # Draw feature backgrounds first (behind everything)
            if self.show_features and hasattr(self.record, 'features'):
                self.draw_feature_backgrounds(painter, i, line_seq, block_start_y)
            
            # Draw feature annotation line
            if self.show_features:
                self.draw_feature_line(painter, i, line_seq, y_pos)
                y_pos += self.line_height
            
            # Draw ruler
            self.draw_ruler(painter, len(line_seq), y_pos)
            y_pos += self.line_height
            
            # Draw sequence line with line number
            self.draw_sequence_line(painter, line_seq, line_start, y_pos)
            y_pos += self.line_height
            
            # Draw translation if enabled
            if self.show_translation:
                self.draw_translation(painter, line_seq, line_start, y_pos)
                y_pos += self.line_height
            
            # Draw complement if enabled
            if self.show_complement:
                self.draw_complement(painter, line_seq, line_start, y_pos)
                y_pos += self.line_height
            
            # Draw enzyme sites if enabled
            if self.show_enzymes:
                self.draw_enzyme_line(painter, i, line_seq, y_pos)
                y_pos += self.line_height
            
            y_pos += self.line_height  # spacing between blocks
    
    def draw_ruler(self, painter, length, y_pos):
        """Draw position ruler"""
        painter.setFont(self.ruler_font)
        painter.setPen(QPen(self.colors['ruler']))
        
        for i in range(length):
            x = self.get_x_position(i)
            
            if i % 10 == 0:
                painter.drawText(x, y_pos, "|")
            elif i % 5 == 0:
                painter.drawText(x, y_pos, ".")
    
    def draw_feature_line(self, painter, start_pos, sequence, y_pos):
        """Draw feature annotation line"""
        if not hasattr(self.record, 'features'):
            return
        
        painter.setFont(self.base_font)
        
        feature_chars = [" "] * len(sequence)
        
        for feature in self.record.features:
            if feature.type == 'source':
                continue
                
            feat_start = int(feature.location.start)
            feat_end = int(feature.location.end)
            line_end = start_pos + len(sequence)
            
            if feat_start < line_end and feat_end > start_pos:
                overlap_start = max(feat_start - start_pos, 0)
                overlap_end = min(feat_end - start_pos, len(sequence))
                
                char = self.get_feature_char(feature.type)
                for pos in range(overlap_start, overlap_end):
                    if 0 <= pos < len(feature_chars):
                        feature_chars[pos] = char
        
        # Draw feature characters
        painter.setPen(QPen(Qt.darkGray))
        for i, char in enumerate(feature_chars):
            x = self.get_x_position(i)
            painter.drawText(x, y_pos, char)
    
    def get_feature_char(self, feature_type):
        """Get character for feature type"""
        chars = {
            "gene": "=",
            "CDS": "-", 
            "promoter": ">",
            "terminator": "|",
            "rep_origin": "o",
            "misc_feature": ".",
            "regulatory": "~",
            "enhancer": "+",
            "primer_bind": "^",
            "protein_bind": "#"
        }
        return chars.get(feature_type, ".")
    
    def draw_sequence_line(self, painter, sequence, line_number, y_pos):
        """Draw sequence line with colored bases"""
        painter.setFont(self.base_font)
        
        # Draw line number and direction indicator
        painter.setPen(QPen(self.colors['numbers']))
        painter.drawText(10, y_pos, f"{line_number:>6}")
        painter.drawText(self.number_width - 20, y_pos, "5'->")
        
        # Draw sequence with colored bases
        for i, base in enumerate(sequence):
            x = self.get_x_position(i)
            
            # Set color for base
            if base in self.colors:
                painter.setPen(QPen(self.colors[base]))
            else:
                painter.setPen(QPen(Qt.black))
            
            painter.drawText(x, y_pos, base)
    
    def draw_translation(self, painter, sequence, start_pos, y_pos):
        """Draw translation line"""
        painter.setFont(self.base_font)
        painter.setPen(QPen(self.colors['translation']))
        
        # Calculate reading frame offset
        frame_offset = (start_pos - 1) % 3
        
        # Translate sequence in proper reading frame
        for i in range(0, len(sequence) - 2, 3):
            # Check if we're in the correct reading frame
            if (i + frame_offset) % 3 == 0:
                codon = sequence[i:i+3]
                if len(codon) == 3 and all(c in 'ATCG' for c in codon):
                    try:
                        aa = str(Seq(codon).translate())
                        
                        # Position amino acid at the middle of the codon
                        middle_pos = i + 1  # Middle of the 3-base codon
                        x = self.get_x_position(middle_pos)
                        
                        painter.drawText(x, y_pos, aa)
                    except:
                        pass
    
    def draw_complement(self, painter, line_sequence, line_start_pos, y_pos):
        """Draw complement line (antiparallel strand, 3' to 5' direction)"""
        painter.setFont(self.base_font)
        painter.setPen(QPen(self.colors['complement']))
        
        # Draw label for complement (3' to 5' direction)
        painter.drawText(10, y_pos, "<-3'")
        
        comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        
        # For the antiparallel complementary strand:
        # Top strand:    5'-ATCG-3'
        # Bottom strand: 3'-TAGC-5'
        # 
        # When displaying the bottom strand from 3' to 5' (left to right),
        # we show the simple complement of the line sequence
        line_complement = "".join(comp_map.get(base, base) for base in line_sequence)
        
        # Draw the complement bases for this line
        for i, base in enumerate(line_complement):
            x = self.get_x_position(i)
            
            # Color the complement bases
            if base in self.colors:
                painter.setPen(QPen(self.colors[base]))
            else:
                painter.setPen(QPen(self.colors['complement']))
            
            painter.drawText(x, y_pos, base)
    
    def draw_enzyme_line(self, painter, start_pos, sequence, y_pos):
        """Draw enzyme cut sites"""
        if not self.enzyme_sites:
            return
        
        painter.setFont(self.base_font)
        painter.setPen(QPen(self.colors['enzyme']))
        
        enzyme_chars = [" "] * len(sequence)
        line_end = start_pos + len(sequence)
        
        # Mark enzyme cut sites
        for enzyme_name, sites in self.enzyme_sites.items():
            for site in sites:
                cut_pos = site - 1  # Convert to 0-based
                if start_pos <= cut_pos < line_end:
                    relative_pos = cut_pos - start_pos
                    if 0 <= relative_pos < len(enzyme_chars):
                        enzyme_chars[relative_pos] = "|"
        
        # Draw enzyme characters
        for i, char in enumerate(enzyme_chars):
            x = self.get_x_position(i)
            painter.drawText(x, y_pos, char)
    
    def draw_feature_backgrounds(self, painter, start_pos, sequence, y_pos):
        """Draw feature background highlights"""
        if not hasattr(self.record, 'features'):
            return
        
        for feature in self.record.features:
            if feature.type == 'source':
                continue
            
            feat_start = int(feature.location.start)
            feat_end = int(feature.location.end)
            line_end = start_pos + len(sequence)
            
            # Check if feature overlaps with current line
            if feat_start < line_end and feat_end > start_pos:
                overlap_start = max(feat_start - start_pos, 0)
                overlap_end = min(feat_end - start_pos, len(sequence))
                
                if feature.type in self.feature_colors:
                    color = self.feature_colors[feature.type]
                    
                    # Calculate rectangle position using precise positioning
                    x_start = self.get_x_position(overlap_start)
                    x_end = self.get_x_position(overlap_end - 1) + self.char_width
                    width = x_end - x_start
                    
                    # Draw background for multiple lines if needed
                    lines_to_highlight = 1
                    if self.show_features:
                        lines_to_highlight += 1
                    if self.show_translation:
                        lines_to_highlight += 1
                    if self.show_complement:
                        lines_to_highlight += 1
                    if self.show_enzymes:
                        lines_to_highlight += 1
                    
                    painter.fillRect(
                        x_start, y_pos - 5,
                        width, lines_to_highlight * self.line_height,
                        color
                    )
    
    def mousePressEvent(self, event):
        """Handle mouse clicks to track position"""
        if not self.record:
            return
        
        # Convert mouse position to sequence position
        seq_pos = self.mouse_to_sequence_position(event.pos())
        if seq_pos >= 0:
            self.current_position = seq_pos
            self.positionChanged.emit(seq_pos)
    
    def mouse_to_sequence_position(self, pos):
        """Convert mouse position to sequence position"""
        if not self.record:
            return -1
        
        x = pos.x()
        y = pos.y() - self.margin - self.line_height * 2  # Account for header
        
        if x < self.number_width or y < 0:
            return -1
        
        # Calculate which line block we're in
        lines_per_block = 2  # ruler + sequence
        if self.show_features:
            lines_per_block += 1
        if self.show_translation:
            lines_per_block += 1
        if self.show_complement:
            lines_per_block += 1
        if self.show_enzymes:
            lines_per_block += 1
        lines_per_block += 1  # spacing
        
        block_index = y // (lines_per_block * self.line_height)
        
        # Find the closest base position
        closest_base = -1
        min_distance = float('inf')
        
        for i in range(self.bases_per_line):
            base_x = self.get_x_position(i)
            distance = abs(x - base_x)
            if distance < min_distance:
                min_distance = distance
                closest_base = i
        
        if closest_base >= 0:
            seq_pos = block_index * self.bases_per_line + closest_base
            return min(seq_pos, len(self.record.seq) - 1)
        
        return -1
    
    # Methods to match the original interface
    def set_show_features(self, show):
        self.show_features = show
        self.update_size()
        self.update()
    
    def set_show_translation(self, show):
        self.show_translation = show
        self.update_size()
        self.update()
    
    def set_show_complement(self, show):
        self.show_complement = show
        self.update_size()
        self.update()
    
    def set_show_enzymes(self, show):
        self.show_enzymes = show
        self.update_size()
        self.update()
    
    def set_bases_per_line(self, bases):
        self.bases_per_line = bases
        self.update_size()
        self.update()
    
    def display_sequence(self, record):
        self.set_sequence(record)