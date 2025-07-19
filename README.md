# PicoMol: The Miniature Molecular Visualization and Bioinformatics Suite

PicoMol is a powerful, user-friendly desktop application for molecular visualization, structural analysis, and bioinformatics. Built with PyQt5 and NGL.js, it provides a comprehensive suite of tools for researchers, students, and developers working with protein structures and biomolecular data.

## Features

### 🧬 Core Molecular Visualization
- **Fetch PDB Structures:** Download protein data directly from the Protein Data Bank (PDB) using a PDB ID
- **Open Local PDB Files:** Load and visualize protein structures from `.pdb` or `.ent` files
- **Drag-and-Drop Support:** Instantly open files by dragging them onto the PicoMol window
- **Interactive 3D Visualization:** High-quality, interactive 3D rendering powered by NGL.js
- **Multiple Representations:** Switch between cartoon, ball+stick, surface, licorice, and more
- **Flexible Coloring:** Apply various coloring schemes (chain ID, residue type, secondary structure, etc.)
- **Screenshot Export:** Save high-quality images of your molecular structures

### 🔬 BLAST Integration
- **Complete BLAST Suite:** Full implementation of BLASTP, BLASTN, BLASTX, TBLASTN, and TBLASTX
- **NCBI-Style Interface:** Authentic user experience matching the official NCBI BLAST interface
- **Online BLAST Searches:** Connect directly to NCBI's BLAST servers for real-time searches
- **Comprehensive Results:** Detailed results with alignment tables, statistics, and downloadable outputs
- **Sequence Validation:** Built-in sequence validation and formatting
- **Multiple Database Support:** Access to all major NCBI databases (nr, nt, RefSeq, etc.)

### 📊 Sequence Analysis
- **Sequence Display:** View amino acid or nucleotide sequences extracted from PDB structures
- **FASTA Support:** Import and export sequences in FASTA format
- **Sequence Validation:** Automatic validation for protein and nucleotide sequences

### 🎨 User Experience
- **Tabbed Interface:** Organized workspace with separate tabs for visualization, BLAST, and analysis tools
- **Robust Error Handling:** User-friendly error dialogs with actionable suggestions
- **Undo/Redo System:** Full undo/redo support for visualization changes
- **Recent Files:** Quick access to recently opened structures
- **Customizable Settings:** Adjustable background colors, spin controls, and more
- **Welcome Screen:** Optional welcome dialog for new users

## Installation

### Prerequisites
- Python 3.7 or higher
- Internet connection (for downloading NGL.js and BLAST functionality)

### Quick Install

1. **Clone the repository:**
   ```bash
   git clone https://github.com/ultimafounding/PicoMol.git
   cd PicoMol
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Run PicoMol:**
   ```bash
   python picomol.py
   ```

The application will automatically download NGL.js on first run if needed.

### Manual Dependency Installation
If you prefer to install dependencies manually:
```bash
pip install PyQt5 PyQtWebEngine biopython requests
```

## Usage

### Getting Started
1. **Launch PicoMol:** Run `python picomol.py`
2. **Load a Structure:**
   - Enter a PDB ID (e.g., `1CRN`, `4HHB`) and click "Fetch"
   - Use "Open Local PDB File" to load your own structures
   - Drag and drop `.pdb` or `.ent` files directly onto the window

### 3D Visualization
- **Change Representation:** Use the dropdown to switch between cartoon, ball+stick, surface, etc.
- **Adjust Colors:** Select different coloring schemes or use uniform colors
- **Control View:** Toggle auto-rotation, change background color
- **Save Images:** Click "Save Screenshot" to export the current view

### BLAST Searches
1. **Navigate to BLAST Tab:** Click on the "BLAST" tab in the main interface
2. **Choose BLAST Type:** Select from BLASTP, BLASTN, BLASTX, TBLASTN, or TBLASTX
3. **Enter Query:** Paste your sequence or upload a FASTA file
4. **Configure Search:** Select database, adjust parameters as needed
5. **Run Search:** Click "BLAST" to submit your search to NCBI
6. **View Results:** Comprehensive results with alignments and statistics

## Project Structure

```
PicoMol/
├── picomol.py                  # Main application
├── setup_ngl.py               # NGL.js setup utility
├── requirements.txt           # Python dependencies
├── blast_utils/               # BLAST functionality package
│   ├── __init__.py           # Package initialization
│   ├── blast_utils.py        # Core BLAST functionality
│   ├── blast_utils_ncbi_layout.py  # NCBI-style UI layouts
│   └── blast_results_parser.py     # Results parsing and formatting
├── ngl_assets/               # NGL.js library storage
├── pulled_structures/        # Downloaded PDB files and generated viewers
└── release_notes/           # Version history and release notes
```

## Technical Details

### Architecture
## 🧩 Technology Stack

PicoMol combines several powerful technologies:

- **PyQt5:** Cross-platform GUI framework providing the main interface
- **QWebEngineView:** Embedded web browser for 3D visualization
- **NGL.js:** High-performance molecular graphics library (MIT licensed)
- **Biopython:** Bioinformatics tools for sequence and structure analysis
- **NCBI BLAST+:** Local BLAST searches and sequence analysis
- **PDBx/mmCIF:** Structure file format support
- **Biopython:** Molecular biology tools for PDB parsing and sequence handling
- **Local HTTP Server:** Serves molecular data to the embedded browser

### BLAST Implementation
- **Online Integration:** Direct connection to NCBI BLAST servers
- **Asynchronous Processing:** Non-blocking searches with progress indicators
- **Result Parsing:** Comprehensive parsing of XML and text BLAST outputs
- **NCBI Compliance:** Interface design matches official NCBI BLAST pages

## Requirements

- **Python:** 3.7+
- **PyQt5:** 5.15.0+
- **PyQtWebEngine:** 5.15.0+
- **Biopython:** 1.79+
- **Requests:** 2.25.0+

See `requirements.txt` for complete dependency list.

## Contributing

We welcome contributions! Please feel free to:
- Report bugs or request features via [GitHub Issues](https://github.com/ultimafounding/PicoMol/issues)
- Submit pull requests for improvements
- Share feedback and suggestions

## 📚 Citations

If you use PicoMol in your research, please cite the following:

### NGL.js
```
Rose AS, Bradley AR, Valasatava Y, Duarte JM, Prlić A, Rose PW.
NGL viewer: web-based molecular graphics for large complexes.
Bioinformatics 34(21): 3755-3758, 2018.
doi:10.1093/bioinformatics/bty419

Rose AS, Hildebrand PW.
NGL Viewer: a web application for molecular visualization.
Nucleic Acids Research 43(W1): W576-W579, 2015.
doi:10.1093/nar/gkv402
```

### BLAST
```
Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ.
Basic local alignment search tool.
Journal of Molecular Biology 215(3): 403-410, 1990.
doi:10.1016/S0022-2836(05)80360-2

Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ.
Gapped BLAST and PSI-BLAST: a new generation of protein database search programs.
Nucleic Acids Research 25(17): 3389-3402, 1997.
doi:10.1093/nar/25.17.3389
```

### PDB Format
```
Berman HM, Westbrook J, Feng Z, Gilliland G, Bhat TN, Weissig H, 
Shindyalov IN, Bourne PE.
The Protein Data Bank.
Nucleic Acids Research 28(1): 235-242, 2000.
doi:10.1093/nar/28.1.235
```

## 📄 License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

### Third-Party Licenses
- **NGL.js:** MIT License (compatible with GPL v3.0)
- **Biopython:** Biopython License (compatible with GPL v3.0)


When using PicoMol in your research, please cite:

**For NGL.js visualization:**
- AS Rose, AR Bradley, Y Valasatava, JM Duarte, A Prlić and PW Rose. NGL viewer: web-based molecular graphics for large complexes. *Bioinformatics*: bty419, 2018. [doi:10.1093/bioinformatics/bty419](https://doi.org/10.1093/bioinformatics/bty419)
- AS Rose and PW Hildebrand. NGL Viewer: a web application for molecular visualization. *Nucleic Acids Research* 43 (W1): W576-W579, 2015. [doi:10.1093/nar/gkv402](https://doi.org/10.1093/nar/gkv402)

**For BLAST functionality:**
- Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. Basic local alignment search tool. *Journal of Molecular Biology* 215, 403-410, 1990.

## Support

- **Documentation:** Check this README and inline help tooltips
- **Issues:** Report problems on [GitHub Issues](https://github.com/ultimafounding/PicoMol/issues)
- **Feedback:** Use the built-in feedback dialog (Help → Send Feedback)

---

*PicoMol: Making molecular visualization and bioinformatics accessible to everyone.*
