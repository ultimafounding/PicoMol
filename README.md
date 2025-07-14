# PicoMol: The Miniature Molecular Visualisation System

PicoMol is a simple desktop application built with PyQt5 and NGL.js for visualizing protein structures. It allows users to fetch protein data from the PDB database or load local PDB files, and then display them with various representations and coloring schemes.

## Features

*   **Fetch PDB Structures:** Download protein data directly from the Protein Data Bank (PDB) using a PDB ID.
*   **Open Local PDB Files:** Load and visualize protein structures from `.pdb` or `.ent` files stored on your local machine.

*   **Improved User Interface:** Features a cleaner, more organized layout with controls in a sidebar and a larger viewing area for protein structures.
*   **Interactive 3D Visualization:** Utilizes NGL.js for high-quality, interactive 3D rendering of protein structures.
*   **Representation Options:** Switch between different molecular representations (e.g., cartoon, ball+stick, surface, licorice).
*   **Coloring Schemes:** Apply various coloring schemes based on atom properties, residue types, chain IDs, and more.
*   **Spin Control:** Toggle automatic rotation of the protein structure.
*   **Background Color Customization:** Change the background color of the viewer.
*   **Water Molecule Removal:** Option to hide water molecules from the visualization.

## Requirements

*   Python 3.x
*   `PyQt5`
*   `Bio.PDB` (part of Biopython)
*   `requests` (for downloading NGL.js)

## Installation

1.  **Clone the repository (or download the files):**

    ```bash
    git clone https://github.com/your-username/PicoMol.git
    cd PicoMol
    ```

2.  **Install dependencies:**

    ```bash
    pip install PyQt5 biopython requests
    ```

3.  **Download NGL.js:**

    The application will attempt to download `ngl.min.js` automatically when first run if it's not found. Alternatively, you can run the setup script manually:

    ```bash
    python setup_ngl.py
    ```

    This will download `ngl.min.js` into the `ngl_assets` directory.

## Usage

To run the application, execute `picomol.py`:

```bash
python picomol.py
```

Once the application starts:

1.  **Fetch PDB:** Enter a PDB ID (e.g., `1CRN`, `4HHB`) into the text field and click "Fetch PDB". The application will download the structure and display it.
2.  **Open Local PDB:** Click "Open Local PDB File" to browse and select a `.pdb` or `.ent` file from your computer.
3.  **Adjust Visualization:** Use the dropdown menus and checkboxes to change the representation, coloring, spin, and background color.

## How it Works

PicoMol combines a PyQt5 graphical user interface with a local HTTP server and the NGL.js library for 3D visualization.

*   **Backend (Python):**
    *   `picomol.py`: Manages the PyQt5 GUI, handles user input, fetches PDB files using Biopython's `PDBList`, and serves the NGL.js viewer via a `SimpleHTTPRequestHandler` in a separate thread.
    *   `setup_ngl.py`: A utility script to download the `ngl.min.js` library from a CDN.
*   **Frontend (HTML/JavaScript/NGL.js):**
    *   When a protein structure is loaded, `picomol.py` generates a dynamic HTML file that embeds NGL.js. This HTML file is then loaded into a `QWebEngineView` (a web browser component) within the PyQt5 application.
    *   NGL.js (MIT licensed) handles the parsing of the PDB data and the 3D rendering. JavaScript functions are exposed to the Python side, allowing the PyQt5 application to control NGL.js visualization parameters (representation, color, spin, etc.). The MIT license is compatible with the project's GPLv3 license, allowing for its inclusion.

## Project Structure

```
PicoMol/
├── picomol.py              # Main application script
├── setup_ngl.py            # Script to download ngl.min.js
├── ngl_assets/             # Directory for NGL.js library
│   └── ngl.min.js          # NGL.js library (downloaded by setup_ngl.py)
└── pulled_structures/      # Directory for downloaded PDB files and generated HTML viewers
    ├── example.html        # Example generated HTML viewer
    ├── example.pdb         # Example downloaded PDB file
    └── ngl.min.js          # Copy of ngl.min.js for local serving
```

## License

This project is open-source and available under the GNU General Public License v3. See the `LICENSE` file for more details.

## Citing NGL.js

When using NGL.js, please cite:

*   AS Rose, AR Bradley, Y Valasatava, JM Duarte, A Prlić and PW Rose. NGL viewer: web-based molecular graphics for large complexes. Bioinformatics: bty419, 2018. doi:10.1093/bioinformatics/bty419
*   AS Rose and PW Hildebrand. NGL Viewer: a web application for molecular visualization. Nucl Acids Res (1 July 2015) 43 (W1): W576-W579 first published online April 29, 2015. doi:10.1093/nar/gkv402
