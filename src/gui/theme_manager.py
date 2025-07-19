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
        """Get dark theme palette."""
        palette = QPalette()
        
        # Window colors
        palette.setColor(QPalette.Window, QColor(53, 53, 53))
        palette.setColor(QPalette.WindowText, QColor(255, 255, 255))
        
        # Base colors (for input fields)
        palette.setColor(QPalette.Base, QColor(25, 25, 25))
        palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
        
        # Text colors
        palette.setColor(QPalette.Text, QColor(255, 255, 255))
        palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
        
        # Button colors
        palette.setColor(QPalette.Button, QColor(53, 53, 53))
        palette.setColor(QPalette.ButtonText, QColor(255, 255, 255))
        
        # Highlight colors
        palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        palette.setColor(QPalette.HighlightedText, QColor(0, 0, 0))
        
        # Disabled colors
        palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(127, 127, 127))
        palette.setColor(QPalette.Disabled, QPalette.Text, QColor(127, 127, 127))
        palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(127, 127, 127))
        
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
                QToolTip {
                    color: #ffffff;
                    background-color: #2a2a2a;
                    border: 1px solid #555555;
                }
                QMenuBar {
                    background-color: #353535;
                    color: #ffffff;
                }
                QMenuBar::item {
                    background-color: transparent;
                }
                QMenuBar::item:selected {
                    background-color: #2a82da;
                }
                QMenu {
                    background-color: #353535;
                    color: #ffffff;
                    border: 1px solid #555555;
                }
                QMenu::item:selected {
                    background-color: #2a82da;
                }
                QTabWidget::pane {
                    border: 1px solid #555555;
                    background-color: #353535;
                }
                QTabBar::tab {
                    background-color: #2a2a2a;
                    color: #ffffff;
                    padding: 8px 16px;
                    margin-right: 2px;
                }
                QTabBar::tab:selected {
                    background-color: #353535;
                }
                QGroupBox {
                    color: #ffffff;
                    border: 2px solid #555555;
                    border-radius: 5px;
                    margin-top: 1ex;
                }
                QGroupBox::title {
                    subcontrol-origin: margin;
                    left: 10px;
                    padding: 0 5px 0 5px;
                }
                QScrollBar:vertical {
                    background-color: #2a2a2a;
                    width: 15px;
                    border-radius: 0px;
                }
                QScrollBar::handle:vertical {
                    background-color: #555555;
                    min-height: 5px;
                    border-radius: 4px;
                }
                QScrollBar::handle:vertical:hover {
                    background-color: #666666;
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