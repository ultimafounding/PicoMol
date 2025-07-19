#!/usr/bin/env python3
"""
Preferences system for PicoMol.

This module provides a comprehensive preferences dialog and settings management
for customizing various aspects of the PicoMol application.
"""

import os
import json
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTabWidget, QWidget, QLabel,
    QLineEdit, QPushButton, QCheckBox, QComboBox, QSpinBox, QDoubleSpinBox,
    QGroupBox, QFormLayout, QColorDialog, QFileDialog, QMessageBox,
    QDialogButtonBox, QSlider, QTextEdit, QScrollArea, QFrame, QButtonGroup,
    QRadioButton, QGridLayout, QSizePolicy, QApplication, QFontDialog
)
from PyQt5.QtCore import QSettings, Qt, pyqtSignal
from PyQt5.QtGui import QFont, QColor, QPalette


class PreferencesDialog(QDialog):
    """Main preferences dialog for PicoMol."""
    
    # Signal emitted when preferences are applied
    preferences_applied = pyqtSignal()
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("PicoMol Preferences")
        self.setMinimumSize(800, 600)
        self.setModal(True)
        
        # Initialize settings
        self.settings = QSettings("PicoMolApp", "PicoMol")
        
        # Create the UI
        self.init_ui()
        
        # Load current settings
        self.load_preferences()
    
    def init_ui(self):
        """Initialize the preferences dialog UI."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)
        
        # Create tab widget
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)
        
        # Create preference tabs
        self.create_general_tab()
        self.create_visualization_tab()
        self.create_interface_tab()
        
        # Button box
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel | 
            QDialogButtonBox.Apply | QDialogButtonBox.RestoreDefaults
        )
        button_box.accepted.connect(self.accept_preferences)
        button_box.rejected.connect(self.reject)
        button_box.button(QDialogButtonBox.Apply).clicked.connect(self.apply_preferences)
        button_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self.restore_defaults)
        layout.addWidget(button_box)
    
    def create_general_tab(self):
        """Create the General preferences tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(15)
        
        # Startup group
        startup_group = QGroupBox("Startup")
        startup_layout = QFormLayout(startup_group)
        
        self.show_welcome_cb = QCheckBox("Show welcome dialog on startup")
        self.show_welcome_cb.setToolTip("Display the welcome screen when PicoMol starts")
        startup_layout.addRow(self.show_welcome_cb)
        
        layout.addWidget(startup_group)
        
        # File handling group
        file_group = QGroupBox("File Handling")
        file_layout = QFormLayout(file_group)
        
        self.max_recent_files = QSpinBox()
        self.max_recent_files.setRange(5, 50)
        self.max_recent_files.setValue(10)
        self.max_recent_files.setToolTip("Maximum number of recent files to remember")
        file_layout.addRow("Max recent files:", self.max_recent_files)
        
        layout.addWidget(file_group)
        
        layout.addStretch()
        self.tab_widget.addTab(tab, "General")
    
    def create_visualization_tab(self):
        """Create the Visualization preferences tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(15)
        
        # Default visualization settings
        viz_group = QGroupBox("Default Visualization Settings")
        viz_layout = QFormLayout(viz_group)
        
        self.default_representation = QComboBox()
        self.default_representation.addItems([
            "cartoon", "ball+stick", "spacefill", "surface", 
            "ribbon", "licorice", "tube"
        ])
        self.default_representation.setCurrentText("cartoon")
        viz_layout.addRow("Default representation:", self.default_representation)
        
        self.default_color_scheme = QComboBox()
        self.default_color_scheme.addItems([
            "chainid", "residueindex", "sstruc", "resname", 
            "element", "uniform"
        ])
        self.default_color_scheme.setCurrentText("chainid")
        viz_layout.addRow("Default color scheme:", self.default_color_scheme)
        
        self.default_background_color = QPushButton()
        self.default_background_color.setFixedSize(50, 30)
        self.default_background_color.setStyleSheet("background-color: black; border: 1px solid gray;")
        self.default_background_color.clicked.connect(self.pick_background_color)
        viz_layout.addRow("Default background:", self.default_background_color)
        
        self.auto_spin_cb = QCheckBox("Enable auto-rotation by default")
        viz_layout.addRow(self.auto_spin_cb)
        
        layout.addWidget(viz_group)
        
        layout.addStretch()
        self.tab_widget.addTab(tab, "Visualization")
    

    
    def create_interface_tab(self):
        """Create the Interface preferences tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(15)
        
        # Appearance
        appearance_group = QGroupBox("Appearance")
        appearance_layout = QFormLayout(appearance_group)
        
        # Theme selection
        self.theme_combo = QComboBox()
        self.theme_combo.addItems(["System Default", "Light", "Dark", "Blue", "Green"])
        self.theme_combo.setToolTip("Choose the application theme")
        appearance_layout.addRow("Theme:", self.theme_combo)
        
        self.font_button = QPushButton("Choose Font...")
        self.font_button.clicked.connect(self.choose_font)
        self.current_font = QApplication.font()
        self.update_font_button_text()
        appearance_layout.addRow("Application font:", self.font_button)
        
        layout.addWidget(appearance_group)
        
        # Behavior
        behavior_group = QGroupBox("Interface Behavior")
        behavior_layout = QFormLayout(behavior_group)
        
        self.tooltips_cb = QCheckBox("Show tooltips")
        self.tooltips_cb.setChecked(True)
        behavior_layout.addRow(self.tooltips_cb)
        
        layout.addWidget(behavior_group)
        
        layout.addStretch()
        self.tab_widget.addTab(tab, "Interface")
    

    

    
    def pick_background_color(self):
        """Pick default background color."""
        color = QColorDialog.getColor()
        if color.isValid():
            self.default_background_color.setStyleSheet(
                f"background-color: {color.name()}; border: 1px solid gray;"
            )
            self.default_background_color.setProperty("color", color.name())
    
    def choose_font(self):
        """Choose application font."""
        font, ok = QFontDialog.getFont(self.current_font, self)
        if ok:
            self.current_font = font
            self.update_font_button_text()
    
    def update_font_button_text(self):
        """Update font button text to show current font."""
        font_text = f"{self.current_font.family()}, {self.current_font.pointSize()}pt"
        if self.current_font.bold():
            font_text += ", Bold"
        if self.current_font.italic():
            font_text += ", Italic"
        self.font_button.setText(font_text)
    
    def load_preferences(self):
        """Load preferences from settings."""
        # General tab
        self.show_welcome_cb.setChecked(
            self.settings.value("show_welcome", True, type=bool)
        )
        self.max_recent_files.setValue(
            self.settings.value("max_recent_files", 10, type=int)
        )
        
        # Visualization tab
        self.default_representation.setCurrentText(
            self.settings.value("default_representation", "cartoon", type=str)
        )
        self.default_color_scheme.setCurrentText(
            self.settings.value("default_color_scheme", "chainid", type=str)
        )
        bg_color = self.settings.value("default_background_color", "black", type=str)
        self.default_background_color.setStyleSheet(
            f"background-color: {bg_color}; border: 1px solid gray;"
        )
        self.default_background_color.setProperty("color", bg_color)
        
        self.auto_spin_cb.setChecked(
            self.settings.value("auto_spin", False, type=bool)
        )
        
        # Interface tab
        self.theme_combo.setCurrentText(
            self.settings.value("theme", "System Default", type=str)
        )
        font_family = self.settings.value("font_family", "", type=str)
        font_size = self.settings.value("font_size", 0, type=int)
        if font_family and font_size:
            self.current_font = QFont(font_family, font_size)
            self.update_font_button_text()
        
        self.tooltips_cb.setChecked(
            self.settings.value("show_tooltips", True, type=bool)
        )
    
    def save_preferences(self):
        """Save preferences to settings."""
        # General tab
        self.settings.setValue("show_welcome", self.show_welcome_cb.isChecked())
        self.settings.setValue("max_recent_files", self.max_recent_files.value())
        
        # Visualization tab
        self.settings.setValue("default_representation", self.default_representation.currentText())
        self.settings.setValue("default_color_scheme", self.default_color_scheme.currentText())
        bg_color = self.default_background_color.property("color") or "black"
        self.settings.setValue("default_background_color", bg_color)
        self.settings.setValue("auto_spin", self.auto_spin_cb.isChecked())
        
        # Interface tab
        self.settings.setValue("theme", self.theme_combo.currentText())
        self.settings.setValue("font_family", self.current_font.family())
        self.settings.setValue("font_size", self.current_font.pointSize())
        self.settings.setValue("show_tooltips", self.tooltips_cb.isChecked())
    
    def apply_preferences(self):
        """Apply preferences without closing dialog."""
        self.save_preferences()
        self.preferences_applied.emit()
        QMessageBox.information(self, "Preferences Applied", 
                              "Preferences have been applied successfully.")
    
    def accept_preferences(self):
        """Accept and apply preferences, then close dialog."""
        self.save_preferences()
        self.preferences_applied.emit()
        self.accept()
    
    def restore_defaults(self):
        """Restore all preferences to default values."""
        reply = QMessageBox.question(
            self, "Restore Defaults",
            "Are you sure you want to restore all preferences to their default values?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            # Clear all settings
            self.settings.clear()
            # Reload default values
            self.load_preferences()
            QMessageBox.information(self, "Defaults Restored", 
                                  "All preferences have been restored to default values.")
    



class PreferencesManager:
    """Manager class for handling application preferences."""
    
    def __init__(self):
        self.settings = QSettings("PicoMolApp", "PicoMol")
    
    def get_preference(self, key, default=None, type_hint=None):
        """Get a preference value."""
        if type_hint:
            return self.settings.value(key, default, type=type_hint)
        return self.settings.value(key, default)
    
    def set_preference(self, key, value):
        """Set a preference value."""
        self.settings.setValue(key, value)
    

    
    def get_visualization_defaults(self):
        """Get default visualization settings."""
        return {
            'representation': self.get_preference('default_representation', 'cartoon', str),
            'color_scheme': self.get_preference('default_color_scheme', 'chainid', str),
            'background_color': self.get_preference('default_background_color', 'black', str),
            'auto_spin': self.get_preference('auto_spin', False, bool)
        }
    
    def get_interface_settings(self):
        """Get interface settings."""
        return {
            'theme': self.get_preference('theme', 'System Default', str),
            'font_family': self.get_preference('font_family', '', str),
            'font_size': self.get_preference('font_size', 0, int),
            'show_tooltips': self.get_preference('show_tooltips', True, bool)
        }


def show_preferences_dialog(parent=None):
    """Show the preferences dialog."""
    dialog = PreferencesDialog(parent)
    return dialog.exec_()


if __name__ == "__main__":
    # Test the preferences dialog
    import sys
    app = QApplication(sys.argv)
    
    dialog = PreferencesDialog()
    dialog.show()
    
    sys.exit(app.exec_())