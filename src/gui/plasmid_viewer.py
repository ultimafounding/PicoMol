from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton

def create_plasmid_viewer_tab(parent):
    """Creates the Plasmid Viewer tab."""
    tab = QWidget()
    layout = QVBoxLayout(tab)
    
    load_button = QPushButton("Load Plasmid File")
    layout.addWidget(load_button)
    
    return tab
