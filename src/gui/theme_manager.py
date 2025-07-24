#!/usr/bin/env python3
"""
Theme Manager for PicoMol.

This module provides theme management functionality for the PicoMol application,
including predefined themes and the ability to apply them to the application.
"""

from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPalette, QColor


class ThemeManager:
    """Manages application themes for PicoMol."""
    
    def __init__(self):
        self.themes = {
            "System Default": self._get_system_theme,
            "Light": self._get_light_theme,
            "Dark": self._get_dark_theme,
            "Blue": self._get_blue_theme,
            "Green": self._get_green_theme
        }
    
    def apply_theme(self, theme_name):
        """Apply a theme to the application."""
        if theme_name in self.themes:
            app = QApplication.instance()
            if app:
                if theme_name == "System Default":
                    # Reset to system default
                    app.setPalette(app.style().standardPalette())
                    app.setStyleSheet("")
                else:
                    palette = self.themes[theme_name]()
                    app.setPalette(palette)
                    # Apply additional stylesheet if needed
                    stylesheet = self._get_theme_stylesheet(theme_name)
                    app.setStyleSheet(stylesheet)
                return True
        return False
    
    def get_available_themes(self):
        """Get list of available theme names."""
        return list(self.themes.keys())
    
    def _get_system_theme(self):
        """Get system default theme."""
        app = QApplication.instance()
        return app.style().standardPalette() if app else QPalette()
    
    def _get_light_theme(self):
        """Get light theme palette."""
        palette = QPalette()
        
        # Window colors
        palette.setColor(QPalette.Window, QColor(240, 240, 240))
        palette.setColor(QPalette.WindowText, QColor(0, 0, 0))
        
        # Base colors (for input fields)
        palette.setColor(QPalette.Base, QColor(255, 255, 255))
        palette.setColor(QPalette.AlternateBase, QColor(245, 245, 245))
        
        # Text colors
        palette.setColor(QPalette.Text, QColor(0, 0, 0))
        palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
        
        # Button colors
        palette.setColor(QPalette.Button, QColor(225, 225, 225))
        palette.setColor(QPalette.ButtonText, QColor(0, 0, 0))
        
        # Highlight colors
        palette.setColor(QPalette.Highlight, QColor(76, 163, 224))
        palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))
        
        # Disabled colors
        palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(120, 120, 120))
        palette.setColor(QPalette.Disabled, QPalette.Text, QColor(120, 120, 120))
        palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(120, 120, 120))
        
        return palette
    
    def _get_dark_theme(self):
        """Get dark theme palette with improved contrast and readability."""
        palette = QPalette()
        
        # Window colors
        palette.setColor(QPalette.Window, QColor(45, 45, 48))  # Slightly lighter than before
        palette.setColor(QPalette.WindowText, QColor(240, 240, 240))  # Off-white for better readability
        
        # Base colors (for input fields)
        palette.setColor(QPalette.Base, QColor(30, 30, 30))  # Slightly lighter than before
        palette.setColor(QPalette.AlternateBase, QColor(45, 45, 48))
        
        # Text colors
        palette.setColor(QPalette.Text, QColor(240, 240, 240))  # Brighter text for better contrast
        palette.setColor(QPalette.BrightText, QColor(255, 180, 180))  # Softer red for bright text
        
        # Button colors
        palette.setColor(QPalette.Button, QColor(63, 63, 70))  # Lighter button background
        palette.setColor(QPalette.ButtonText, QColor(240, 240, 240))  # Brighter button text
        
        # Highlight colors
        palette.setColor(QPalette.Highlight, QColor(0, 122, 204))  # Brighter blue for better visibility
        palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))  # White text on highlight
        
        # Disabled colors
        palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(160, 160, 160))  # Lighter gray for disabled
        palette.setColor(QPalette.Disabled, QPalette.Text, QColor(160, 160, 160))
        palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(140, 140, 140))  # Slightly darker for buttons
        
        # Additional UI elements
        palette.setColor(QPalette.Link, QColor(86, 156, 214))  # Brighter link color
        palette.setColor(QPalette.LinkVisited, QColor(156, 86, 214))  # Purple for visited links
        palette.setColor(QPalette.ToolTipBase, QColor(45, 45, 48))  # Match window color
        palette.setColor(QPalette.ToolTipText, QColor(240, 240, 240))  # Match window text
        palette.setColor(QPalette.PlaceholderText, QColor(120, 120, 120))  # Placeholder text color
        
        return palette
    
    def _get_blue_theme(self):
        """Get blue theme palette."""
        palette = QPalette()
        
        # Window colors
        palette.setColor(QPalette.Window, QColor(230, 240, 250))
        palette.setColor(QPalette.WindowText, QColor(0, 0, 0))
        
        # Base colors (for input fields)
        palette.setColor(QPalette.Base, QColor(255, 255, 255))
        palette.setColor(QPalette.AlternateBase, QColor(240, 248, 255))
        
        # Text colors
        palette.setColor(QPalette.Text, QColor(0, 0, 0))
        palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
        
        # Button colors
        palette.setColor(QPalette.Button, QColor(200, 220, 240))
        palette.setColor(QPalette.ButtonText, QColor(0, 0, 0))
        
        # Highlight colors
        palette.setColor(QPalette.Highlight, QColor(0, 120, 215))
        palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))
        
        # Disabled colors
        palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(120, 120, 120))
        palette.setColor(QPalette.Disabled, QPalette.Text, QColor(120, 120, 120))
        palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(120, 120, 120))
        
        return palette
    
    def _get_green_theme(self):
        """Get green theme palette."""
        palette = QPalette()
        
        # Window colors
        palette.setColor(QPalette.Window, QColor(240, 250, 240))
        palette.setColor(QPalette.WindowText, QColor(0, 0, 0))
        
        # Base colors (for input fields)
        palette.setColor(QPalette.Base, QColor(255, 255, 255))
        palette.setColor(QPalette.AlternateBase, QColor(248, 255, 248))
        
        # Text colors
        palette.setColor(QPalette.Text, QColor(0, 0, 0))
        palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
        
        # Button colors
        palette.setColor(QPalette.Button, QColor(220, 240, 220))
        palette.setColor(QPalette.ButtonText, QColor(0, 0, 0))
        
        # Highlight colors
        palette.setColor(QPalette.Highlight, QColor(34, 139, 34))
        palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))
        
        # Disabled colors
        palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(120, 120, 120))
        palette.setColor(QPalette.Disabled, QPalette.Text, QColor(120, 120, 120))
        palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(120, 120, 120))
        
        return palette
    
    def _get_theme_stylesheet(self, theme_name):
        """Get additional stylesheet for specific themes."""
        stylesheets = {
            "Dark": """
                /* Base styles */
                QWidget {
                    color: #f0f0f0;
                    font-size: 9pt;
                }
                
                /* Tooltip */
                QToolTip {
                    color: #f0f0f0;
                    background-color: #2d2d30;
                    border: 1px solid #3f3f46;
                    padding: 4px;
                    border-radius: 3px;
                }
                
                /* Menu Bar */
                QMenuBar {
                    background-color: #2d2d30;
                    color: #f0f0f0;
                    border-bottom: 1px solid #3f3f46;
                }
                QMenuBar::item {
                    background-color: transparent;
                    padding: 4px 8px;
                    border-radius: 3px;
                }
                QMenuBar::item:selected {
                    background-color: #3f3f46;
                }
                QMenuBar::item:pressed {
                    background-color: #007acc;
                }
                
                /* Menus */
                QMenu {
                    background-color: #2d2d30;
                    color: #f0f0f0;
                    border: 1px solid #3f3f46;
                    padding: 4px;
                }
                QMenu::item {
                    padding: 6px 24px 6px 28px;
                    border-radius: 3px;
                }
                QMenu::item:selected {
                    background-color: #3f3f46;
                }
                QMenu::item:pressed {
                    background-color: #007acc;
                }
                QMenu::icon {
                    left: 5px;
                }
                QMenu::separator {
                    height: 1px;
                    background: #3f3f46;
                    margin: 4px 8px;
                }
                
                /* Tab Widget */
                QTabWidget::pane {
                    border: 1px solid #3f3f46;
                    background-color: #2d2d30;
                    border-radius: 4px;
                    margin: 0px;
                    padding: 0px;
                }
                QTabBar::tab {
                    background-color: #2d2d30;
                    color: #f0f0f0;
                    padding: 8px 16px;
                    margin-right: 2px;
                    border: 1px solid #3f3f46;
                    border-bottom: none;
                    border-top-left-radius: 4px;
                    border-top-right-radius: 4px;
                }
                QTabBar::tab:selected, QTabBar::tab:hover {
                    background-color: #1e1e1e;
                    border-color: #3f3f46;
                    border-bottom-color: #1e1e1e;
                }
                QTabBar::tab:!selected {
                    margin-top: 2px;
                }
                
                /* Group Box */
                QGroupBox {
                    color: #f0f0f0;
                    border: 1px solid #3f3f46;
                    border-radius: 4px;
                    margin-top: 1.5em;
                    padding-top: 10px;
                }
                QGroupBox::title {
                    subcontrol-origin: margin;
                    left: 12px;
                    padding: 0 5px 0 5px;
                }
                
                /* Scroll Bars */
                QScrollBar:vertical {
                    background-color: #2d2d30;
                    width: 12px;
                    border-radius: 6px;
                }
                QScrollBar::handle:vertical {
                    background-color: #3f3f46;
                    min-height: 30px;
                    border-radius: 6px;
                }
                QScrollBar::handle:vertical:hover {
                    background-color: #5f5f66;
                }
                QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
                    height: 0px;
                }
                
                /* Buttons */
                QPushButton {
                    background-color: #3f3f46;
                    color: #f0f0f0;
                    border: 1px solid #3f3f46;
                    border-radius: 4px;
                    padding: 5px 12px;
                    min-width: 80px;
                }
                QPushButton:hover {
                    background-color: #5f5f66;
                    border-color: #5f5f66;
                }
                QPushButton:pressed {
                    background-color: #007acc;
                    border-color: #007acc;
                }
                QPushButton:disabled {
                    background-color: #2d2d30;
                    color: #6d6d6d;
                    border: 1px solid #3f3f46;
                }
                
                /* Line Edits */
                QLineEdit, QTextEdit, QPlainTextEdit, QComboBox, QSpinBox, QDoubleSpinBox {
                    background-color: #1e1e1e;
                    color: #f0f0f0;
                    border: 1px solid #3f3f46;
                    border-radius: 4px;
                    padding: 4px 8px;
                    selection-background-color: #007acc;
                }
                QLineEdit:focus, QTextEdit:focus, QPlainTextEdit:focus {
                    border: 1px solid #007acc;
                }
                
                /* Checkboxes and Radio Buttons */
                QCheckBox::indicator, QRadioButton::indicator {
                    width: 16px;
                    height: 16px;
                }
                QCheckBox::indicator:unchecked, QRadioButton::indicator:unchecked {
                    border: 1px solid #3f3f46;
                    background-color: #2d2d30;
                }
                QCheckBox::indicator:checked, QRadioButton::indicator:checked {
                    border: 1px solid #007acc;
                    background-color: #007acc;
                }
                QRadioButton::indicator {
                    border-radius: 8px;
                }
                QRadioButton::indicator:checked {
                    background-color: #2d2d30;
                    border: 4px solid #007acc;
                }
                
                /* Progress Bar */
                QProgressBar {
                    border: 1px solid #3f3f46;
                    border-radius: 4px;
                    text-align: center;
                    background-color: #2d2d30;
                }
                QProgressBar::chunk {
                    background-color: #007acc;
                    width: 10px;
                    margin: 0.5px;
                }
                
                /* Slider */
                QSlider::groove:horizontal {
                    border: 1px solid #3f3f46;
                    height: 4px;
                    background: #2d2d30;
                    margin: 2px 0;
                    border-radius: 2px;
                }
                QSlider::handle:horizontal {
                    background: #f0f0f0;
                    border: 1px solid #3f3f46;
                    width: 12px;
                    margin: -6px 0;
                    border-radius: 6px;
                }
                QSlider::handle:horizontal:hover {
                    background: #ffffff;
                    border-color: #007acc;
                }
            """,
            "Blue": """
                QToolTip {
                    color: #000000;
                    background-color: #e6f3ff;
                    border: 1px solid #0078d4;
                }
                QMenuBar {
                    background-color: #cce7ff;
                }
                QMenuBar::item:selected {
                    background-color: #0078d4;
                    color: #ffffff;
                }
                QMenu::item:selected {
                    background-color: #0078d4;
                    color: #ffffff;
                }
                QTabWidget::pane {
                    border: 1px solid #0078d4;
                }
                QTabBar::tab:selected {
                    background-color: #e6f3ff;
                }
                QGroupBox {
                    border: 2px solid #0078d4;
                    border-radius: 5px;
                    margin-top: 1ex;
                }
            """,
            "Green": """
                QToolTip {
                    color: #000000;
                    background-color: #f0fff0;
                    border: 1px solid #228b22;
                }
                QMenuBar {
                    background-color: #e6ffe6;
                }
                QMenuBar::item:selected {
                    background-color: #228b22;
                    color: #ffffff;
                }
                QMenu::item:selected {
                    background-color: #228b22;
                    color: #ffffff;
                }
                QTabWidget::pane {
                    border: 1px solid #228b22;
                }
                QTabBar::tab:selected {
                    background-color: #f0fff0;
                }
                QGroupBox {
                    border: 2px solid #228b22;
                    border-radius: 5px;
                    margin-top: 1ex;
                }
            """
        }
        return stylesheets.get(theme_name, "")


# Global theme manager instance
theme_manager = ThemeManager()


def apply_theme(theme_name):
    """Convenience function to apply a theme."""
    return theme_manager.apply_theme(theme_name)


def get_available_themes():
    """Convenience function to get available themes."""
    return theme_manager.get_available_themes()


if __name__ == "__main__":
    # Test the theme manager
    import sys
    from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QPushButton, QLabel
    
    app = QApplication(sys.argv)
    
    class TestWindow(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("Theme Manager Test")
            self.setGeometry(100, 100, 400, 300)
            
            central = QWidget()
            self.setCentralWidget(central)
            layout = QVBoxLayout(central)
            
            layout.addWidget(QLabel("Theme Manager Test"))
            
            for theme in get_available_themes():
                btn = QPushButton(f"Apply {theme}")
                btn.clicked.connect(lambda checked, t=theme: apply_theme(t))
                layout.addWidget(btn)
    
    window = TestWindow()
    window.show()
    
    sys.exit(app.exec_())