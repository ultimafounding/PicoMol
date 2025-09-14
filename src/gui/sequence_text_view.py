from PyQt5.QtWidgets import QTextEdit, QWidget, QPlainTextEdit, QGraphicsView, QGraphicsScene
from PyQt5.QtGui import QColor, QTextCharFormat, QFont, QPainter, QTextCursor
from PyQt5.QtCore import Qt, QRect, QSize

class LineNumberArea(QWidget):
    def __init__(self, editor):
        super().__init__(editor)
        self.editor = editor

    def sizeHint(self):
        return QSize(self.editor.lineNumberAreaWidth(), 0)

    def paintEvent(self, event):
        self.editor.lineNumberAreaPaintEvent(event)

class SequenceTextView(QPlainTextEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setFont(QFont("monospace", 10))
        self.lineNumberArea = LineNumberArea(self)

        self.blockCountChanged.connect(self.updateLineNumberAreaWidth)
        self.updateRequest.connect(self.updateLineNumberArea)
        self.cursorPositionChanged.connect(self.highlightCurrentLine)

        self.updateLineNumberAreaWidth(0)

    def lineNumberAreaWidth(self):
        digits = 1
        count = max(1, self.blockCount())
        while count >= 10:
            count /= 10
            digits += 1
        space = 3 + self.fontMetrics().width('9') * digits
        return space

    def updateLineNumberAreaWidth(self, _):
        self.setViewportMargins(self.lineNumberAreaWidth(), 0, 0, 0)

    def updateLineNumberArea(self, rect, dy):
        if dy:
            self.lineNumberArea.scroll(0, dy)
        else:
            self.lineNumberArea.update(0, rect.y(), self.lineNumberArea.width(), rect.height())

        if rect.contains(self.viewport().rect()):
            self.updateLineNumberAreaWidth(0)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        cr = self.contentsRect()
        self.lineNumberArea.setGeometry(QRect(cr.left(), cr.top(), self.lineNumberAreaWidth(), cr.height()))

    def lineNumberAreaPaintEvent(self, event):
        painter = QPainter(self.lineNumberArea)
        painter.fillRect(event.rect(), Qt.lightGray)

        block = self.firstVisibleBlock()
        blockNumber = block.blockNumber()
        top = self.blockBoundingGeometry(block).translated(self.contentOffset()).top()
        bottom = top + self.blockBoundingRect(block).height()

        while block.isValid() and top <= event.rect().bottom():
            if block.isVisible() and bottom >= event.rect().top():
                number = str(blockNumber + 1)
                painter.setPen(Qt.black)
                painter.drawText(0, int(top), self.lineNumberArea.width(), self.fontMetrics().height(),
                                 Qt.AlignRight, number)

            block = block.next()
            top = bottom
            bottom = top + self.blockBoundingRect(block).height()
            blockNumber += 1

    def highlightCurrentLine(self):
        extraSelections = []
        if not self.isReadOnly():
            selection = QTextEdit.ExtraSelection()
            lineColor = QColor(Qt.yellow).lighter(160)
            selection.format.setBackground(lineColor)
            selection.format.setProperty(QTextCharFormat.FullWidthSelection, True)
            selection.cursor = self.textCursor()
            selection.cursor.clearSelection()
            extraSelections.append(selection)
        self.setExtraSelections(extraSelections)

    def display_sequence(self, record, restriction_batch=None):
        self.clear()
        self.record = record
        self.restriction_batch = restriction_batch
        self.setPlainText(str(self.record.seq))

        # Apply feature highlighting
        self.highlight_features()
        
        # Apply restriction site highlighting
        if self.restriction_batch:
            self.highlight_restriction_sites()

    def highlight_features(self):
        cursor = self.textCursor()
        for feature in self.record.features:
            if feature.type in ["gene", "CDS", "promoter"]:
                color = self.get_feature_color(feature.type)
                char_format = self.create_char_format(color)
                
                cursor.setPosition(feature.location.start)
                cursor.setPosition(feature.location.end, QTextCursor.KeepAnchor)
                cursor.setCharFormat(char_format)

    def create_char_format(self, color):
        char_format = QTextCharFormat()
        char_format.setBackground(QColor(color))
        return char_format

    def get_feature_color(self, feature_type):
        colors = {
            "gene": "#3498db",      # Blue
            "CDS": "#2ecc71",       # Green
            "promoter": "#e74c3c",   # Red
        }
        return colors.get(feature_type, "#bdc3c7") # Gray for others
    
    def highlight_restriction_sites(self):
        """Highlight restriction enzyme cut sites exactly like SnapGene"""
        if not self.restriction_batch or not self.record:
            return
        
        try:
            # Get restriction analysis
            analysis = self.restriction_batch.search(self.record.seq)
            
            # Store cut site information
            cut_sites_info = []
            
            # SnapGene-style colors (more professional)
            enzyme_colors = [
                QColor(255, 102, 102),  # Light red
                QColor(102, 255, 102),  # Light green  
                QColor(102, 178, 255),  # Light blue
                QColor(255, 255, 102),  # Light yellow
                QColor(255, 178, 102),  # Light orange
                QColor(178, 102, 255),  # Light purple
                QColor(102, 255, 255),  # Light cyan
                QColor(255, 102, 255),  # Light magenta
            ]
            
            color_index = 0
            
            # Collect cut site information without modifying text yet
            for enzyme, sites in analysis.items():
                if not sites:
                    continue
                    
                enzyme_obj = enzyme
                enzyme_name = str(enzyme_obj)
                recognition_site = str(enzyme_obj.site)
                site_length = len(recognition_site)
                
                # Get enzyme color
                enzyme_color = enzyme_colors[color_index % len(enzyme_colors)]
                color_index += 1
                
                # Get cut position
                try:
                    cut_position = enzyme_obj.fst5
                    if cut_position > 0:
                        cut_position = cut_position - 1  # Convert to 0-based
                    else:
                        cut_position = None
                except:
                    cut_position = None
                
                for site_pos in sites:
                    cut_sites_info.append({
                        'enzyme': enzyme_name,
                        'position': site_pos,
                        'recognition_site': recognition_site,
                        'sequence': str(self.record.seq)[site_pos:site_pos + site_length],
                        'cut_within_site': cut_position,
                        'color': enzyme_color,
                        'site_length': site_length
                    })
            
            # Store cut sites info
            self.cut_sites_info = cut_sites_info
            
            # Create SnapGene-style display
            self.create_snapgene_display(cut_sites_info)
                    
        except Exception as e:
            print(f"Error highlighting restriction sites: {e}")
    
    def create_snapgene_display(self, cut_sites_info):
        """Create SnapGene-style restriction site display"""
        if not cut_sites_info:
            return
        
        try:
            # Get original sequence
            original_sequence = str(self.record.seq)
            
            # Create formatted text with SnapGene-style layout
            formatted_text = self.format_sequence_with_enzymes(original_sequence, cut_sites_info)
            
            # Set the formatted text
            self.setPlainText(formatted_text)
            
            # Apply SnapGene-style formatting
            self.apply_snapgene_formatting(cut_sites_info)
            
        except Exception as e:
            print(f"Error creating SnapGene display: {e}")
    
    def format_sequence_with_enzymes(self, sequence, cut_sites_info):
        """Format sequence text with SnapGene-style enzyme boxes above"""
        lines = []
        bases_per_line = 60  # SnapGene typically uses 60 bases per line
        
        for line_start in range(0, len(sequence), bases_per_line):
            line_end = min(line_start + bases_per_line, len(sequence))
            line_seq = sequence[line_start:line_end]
            
            # Create enzyme annotation line for this sequence line
            enzyme_line = self.create_enzyme_line(line_start, line_end, cut_sites_info)
            
            # Add position numbers
            pos_line = f"{line_start + 1:>6}  "
            
            # Add the lines in SnapGene order: enzymes, position, sequence
            if enzyme_line.strip():  # Only add enzyme line if there are enzymes
                lines.append(" " * 8 + enzyme_line)
            lines.append(pos_line + line_seq)
            lines.append("")  # Empty line for spacing
        
        return "\n".join(lines)
    
    def create_enzyme_line(self, line_start, line_end, cut_sites_info):
        """Create enzyme annotation line for a sequence line"""
        enzyme_line = [" "] * (line_end - line_start)
        
        # Find enzymes that cut in this line
        for site_info in cut_sites_info:
            site_pos = site_info['position']
            enzyme_name = site_info['enzyme']
            site_length = site_info['site_length']
            
            # Check if enzyme site overlaps with this line
            if line_start <= site_pos < line_end:
                relative_pos = site_pos - line_start
                
                # Place enzyme name above the recognition site
                enzyme_text = f"[{enzyme_name}]"
                
                # Try to center the enzyme name over the recognition site
                start_pos = max(0, relative_pos - len(enzyme_text) // 2)
                end_pos = min(len(enzyme_line), start_pos + len(enzyme_text))
                
                # Place the enzyme name
                for i, char in enumerate(enzyme_text):
                    if start_pos + i < len(enzyme_line):
                        enzyme_line[start_pos + i] = char
        
        return "".join(enzyme_line)
    
    def apply_snapgene_formatting(self, cut_sites_info):
        """Apply SnapGene-style formatting to the text"""
        try:
            cursor = self.textCursor()
            text = self.toPlainText()
            lines = text.split('\n')
            
            # Apply formatting line by line
            current_pos = 0
            
            for line in lines:
                if not line.strip():
                    current_pos += len(line) + 1  # +1 for newline
                    continue
                
                # Check if this is an enzyme line (contains [enzyme_name])
                if '[' in line and ']' in line:
                    # Format enzyme annotations
                    import re
                    enzyme_matches = list(re.finditer(r'\[([A-Za-z0-9]+)\]', line))
                    
                    for match in enzyme_matches:
                        enzyme_name = match.group(1)
                        
                        # Find corresponding color
                        enzyme_color = QColor(255, 100, 100)  # Default
                        for site_info in cut_sites_info:
                            if site_info['enzyme'] == enzyme_name:
                                enzyme_color = site_info['color']
                                break
                        
                        # Create SnapGene-style enzyme box format
                        enzyme_format = QTextCharFormat()
                        enzyme_format.setBackground(enzyme_color)
                        enzyme_format.setForeground(QColor(255, 255, 255))  # White text
                        enzyme_format.setFontWeight(QFont.Bold)
                        enzyme_format.setFontPointSize(9)
                        
                        # Apply formatting to the enzyme name
                        start_pos = current_pos + match.start()
                        end_pos = current_pos + match.end()
                        
                        cursor.setPosition(start_pos)
                        cursor.setPosition(end_pos, QTextCursor.KeepAnchor)
                        cursor.setCharFormat(enzyme_format)
                
                elif any(c in 'ATCG' for c in line):  # This is a sequence line
                    # Apply base coloring
                    for i, char in enumerate(line):
                        if char in 'ATCG':
                            char_pos = current_pos + i
                            
                            # Base colors (SnapGene style)
                            base_colors = {
                                'A': QColor(255, 0, 0),    # Red
                                'T': QColor(0, 0, 255),    # Blue
                                'C': QColor(0, 150, 0),    # Green
                                'G': QColor(255, 140, 0)   # Orange
                            }
                            
                            base_format = QTextCharFormat()
                            base_format.setForeground(base_colors.get(char, QColor(0, 0, 0)))
                            base_format.setFontWeight(QFont.Bold)
                            
                            cursor.setPosition(char_pos)
                            cursor.setPosition(char_pos + 1, QTextCursor.KeepAnchor)
                            cursor.setCharFormat(base_format)
                
                current_pos += len(line) + 1  # +1 for newline
                
        except Exception as e:
            print(f"Error applying SnapGene formatting: {e}")
    
    def update_restriction_highlighting(self, restriction_batch):
        """Update restriction site highlighting with new enzyme selection"""
        self.restriction_batch = restriction_batch
        
        # Re-display the sequence to refresh highlighting
        if hasattr(self, 'record') and self.record:
            self.display_sequence(self.record, restriction_batch)
