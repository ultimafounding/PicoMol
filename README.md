# PicoMol: The Miniature Molecular Visualization and Bioinformatics Suite

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

PicoMol is a comprehensive desktop application for molecular visualization, structural analysis, and bioinformatics. Built with PyQt5 and NGL.js, it provides a complete suite of tools for researchers, students, and developers working with protein structures and biomolecular data.

## Features

### 🧬 Core Visualization
- **3D Molecular Viewer**: Interactive visualization powered by NGL.js
- **Multiple Representations**: Cartoon, ball+stick, surface, and more
- **Flexible Coloring**: By chain, residue, secondary structure, etc.
- **PDB Integration**: Direct download and local file support
- **Screenshot Export**: Save high-quality images of your structures

### 🔬 Bioinformatics Tools
- **Sequence Analysis**:
  - Protein properties (MW, pI, hydrophobicity, aromaticity, instability index)
  - Nucleic acid analysis (GC content, translation in all 6 frames, complements)
  - Secondary structure prediction (helix, turn, sheet content)
  - Amino acid and nucleotide composition analysis
- **Motif & Domain Analysis**:
  - InterPro integration for comprehensive domain annotation
  - PROSITE pattern matching with official API
  - Interactive visualization of domains and motifs
  - Sequence context display with highlighting
  - Export capabilities (CSV, Excel, HTML, PDF)
- **Structural Analysis**:
  - Bond angle and distance calculations
  - Secondary structure analysis
  - Ramachandran plot generation
  - Contact map visualization
  - Advanced structural metrics

### 🧪 BLAST Integration
- **All BLAST Variants**: BLASTP, BLASTN, BLASTX, TBLASTN, TBLASTX
- **NCBIBLAST-like Interface**: Familiar experience for users
- **Multiple Database Support**: Access to major NCBI databases
- **Comprehensive Results**: Detailed alignments and statistics

### 🎨 User Experience
- **Modern Interface**: Clean, tabbed layout with dark/light themes
- **Cross-Platform**: Works on Windows, macOS, and Linux
- **Performance Optimized**: Handles large structures efficiently
- **Drag & Drop**: Easy file loading with drag-and-drop support

## 🚀 Installation

### Prerequisites
- Python 3.7 or higher (Python 3.8+ recommended)
- Internet connection (for downloading NGL.js and online BLAST functionality)
- At least 100MB of free disk space

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

**Note:** PicoMol features an organized project structure with source code in `src/`, assets in `assets/`, data files in `data/`, and documentation in `docs/`. The main application file (`picomol.py`) remains in the root directory for easy execution.

### Manual Dependency Installation
If you prefer to install dependencies manually:
```bash
pip install PyQt5 PyQtWebEngine biopython requests
```

## Usage

### Getting Started
1. **Launch PicoMol:** Run `python picomol.py`
2. **Load a Structure:**
   - Enter a PDB ID (e.g., `1CRN`, `4HHB`) and click \"Fetch\"
   - Use \"Open Local PDB File\" to load your own structures
   - Drag and drop `.pdb` or `.ent` files directly onto the window

### 3D Visualization
- **Change Representation:** Use the dropdown to switch between cartoon, ball+stick, surface, etc.
- **Adjust Colors:** Select different coloring schemes or use uniform colors
- **Control View:** Toggle auto-rotation, change background color
- **Save Images:** Click \"Save Screenshot\" to export the current view

### BLAST Searches
1. **Navigate to BLAST Tab:** Click on the \"BLAST\" tab in the main interface
2. **Choose BLAST Type:** Select from BLASTP, BLASTN, BLASTX, TBLASTN, or TBLASTX
3. **Enter Query:** Paste your sequence or upload a FASTA file
4. **Configure Search:** Select database, adjust parameters as needed
5. **Run Search:** Click \"BLAST\" to submit your search to NCBI
6. **View Results:** Comprehensive results with alignments and statistics

### Bioinformatics Analysis
1. **Navigate to Bioinformatics Tab:** Click on the "Bioinformatics" tab
2. **Choose Analysis Type:** Select from the available sub-tabs:
   - **Sequence Analysis:** Comprehensive sequence analysis tools
   - **Protein Motifs and Domains:** Motif and domain analysis tools
   - **Structure:** Structural analysis tools (coming soon)

#### Sequence Analysis
1. **Input Methods:**
   - **Direct Entry:** Paste sequences (plain text or FASTA format)
   - **Load from Structure:** Use sequence from currently loaded PDB structure
   - **Load from File:** Import from FASTA files with multi-sequence support
2. **Analysis Options:**
   - **Automatic Type Detection:** Sequences are automatically classified
   - **Manual Override:** Change sequence type if needed
   - **Comprehensive Results:** Detailed analysis with exportable data

#### Motif and Domain Analysis
1. **Sequence Input:**
   - **Load from Current Structure:** Extract sequence from loaded PDB structure
   - **Load from File:** Import FASTA sequences
   - **Direct Entry:** Paste protein sequences with example provided
2. **Search Options:**
   - **InterPro Search:** Comprehensive domain annotation (Pfam, SMART, PROSITE, CDD)
   - **PROSITE Search:** Motif pattern recognition with official API
   - **Search All Databases:** Combined search for complete analysis
3. **Advanced Features:**
   - **Filtering:** Exclude high-probability motifs for cleaner results
   - **Visualization:** Interactive domain and motif visualization
   - **Export:** Save results and generate detailed HTML reports
   - **Context Display:** Click any motif in results to view sequence context with highlighting
   - **Progress Tracking:** Progress bar remains visible until all searches complete

## Project Structure

```
PicoMol/
├── picomol.py                     # Main application entry point
├── setup_ngl.py                   # NGL.js setup utility
├── requirements.txt               # Python dependencies
├── README.md                      # This file
├── LICENSE                        # GPL v3.0 license
├── .gitignore                     # Git ignore rules
├── src/                          # Source code modules
│   ├── core/                     # Core functionality
│   │   ├── bioinformatics_tools.py  # Comprehensive bioinformatics analysis suite
│   │   ├── motif_analysis.py        # Protein motif and domain analysis tools
│   │   ├── structural_analysis.py   # Structural analysis and calculations
│   │   ├── enhanced_pdb_puller.py   # Enhanced PDB fetching and metadata
│   │   ├── pdb_fetch_worker.py      # Optimized PDB fetch worker
│   │   └── preferences.py           # Settings and preferences management
│   ├── gui/                      # User interface components
│   │   ├── theme_manager.py         # Theme system and UI customization
│   │   └── welcome_dialog.py        # Welcome screen dialog
│   └── blast_utils/              # BLAST functionality package
│       ├── __init__.py              # Package initialization
│       ├── blast_utils.py           # Core BLAST functionality
│       ├── blast_utils_ncbi_layout.py  # NCBI-style UI layouts
│       ├── blast_results_parser.py     # Results parsing and formatting
│       └── accession_validator.py      # Sequence and accession validation
├── assets/                       # Static resources
│   └── ngl_assets/              # NGL.js library storage
├── data/                         # Data files and runtime storage
│   └── pulled_structures/       # Downloaded PDB files and generated viewers
└── docs/                         # Documentation
    ├── CHANGELOG.md             # Version history and changes
    ├── KNOWN_BUGS.md            # Known issues and limitations
    ├── next_release.md          # Upcoming features and changes
    └── release_notes/           # Detailed release notes
```

## Technical Details

### Architecture
PicoMol combines several powerful technologies:

**Core Technologies:**
- **PyQt5:** Cross-platform GUI framework providing the main interface
- **QWebEngineView:** Embedded web browser for 3D visualization
- **NGL.js:** High-performance molecular graphics library (MIT licensed)
- **Biopython:** Comprehensive bioinformatics tools for sequence and structure analysis
- **NCBI BLAST+:** Online BLAST searches and sequence analysis

**Enhanced Features (Optional Dependencies):**
- **NumPy:** Scientific computing for structural calculations
- **Matplotlib:** Advanced plotting and visualization capabilities
- **Pandas:** Data manipulation and export functionality
- **OpenPyXL:** Excel file format support
- **WeasyPrint:** PDF generation capabilities

**Design Principles:**
- **Modular Architecture:** Clean separation of concerns with organized source structure
- **Graceful Degradation:** Core functionality works without optional dependencies
- **Theme System:** Customizable appearance with multiple built-in themes
- **Asynchronous Processing:** Non-blocking operations for better user experience

### Bioinformatics Implementation
- **Dual-Mode Analysis:** Full-featured analysis with Biopython, basic analysis without
- **Smart FASTA Parsing:** Handles both Biopython and manual parsing methods
- **Sequence Validation:** Real-time validation with detailed error reporting
- **Multi-threading:** Background analysis prevents UI freezing
- **API Integration:** Direct connection to InterPro and PROSITE web services
- **Fallback Systems:** Local pattern matching when APIs are unavailable
- **Intelligent Parsing:** Robust parsing of various API response formats
- **Extensible Design:** Easy to add new analysis tools and features

### BLAST Implementation
- **Online Integration:** Direct connection to NCBI BLAST servers
- **Asynchronous Processing:** Non-blocking searches with progress indicators
- **Result Parsing:** Comprehensive parsing of XML and text BLAST outputs
- **NCBI Compliance:** Interface design matches official NCBI BLAST pages

### Motif and Domain Analysis Implementation
- **InterPro API:** Official EBI InterPro web service integration
- **PROSITE API:** Direct connection to ExPASy PROSITE ScanProsite service
- **Hybrid Approach:** API-first with local pattern fallback for reliability
- **Sequence Context Display:** Interactive motif context viewer with highlighting
- **Progress Management:** Smart progress tracking for concurrent searches
- **Result Caching:** Efficient storage and retrieval of analysis results
- **Export System:** Multiple output formats including HTML reports
- **Error Handling:** Robust error recovery and user feedback

## Contributing

We welcome contributions! Please feel free to:
- Report bugs or request features via [GitHub Issues](https://github.com/ultimafounding/PicoMol/issues)
- Submit pull requests for improvements
- Share feedback and suggestions

## 📚 Citations

Here are the citations for projects and tools included in PicoMol:

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

### InterPro
```
Apweiler R, Attwood TK, Bairoch A, Bateman A, Birney E, Biswas M, Bucher P, Cerutti L, Corpet F, Croning MD, et al.
The InterPro database, an integrated documentation resource for protein families, domains and functional sites.
Nucleic Acids Research 29(1): 37-40, 2001.
doi:10.1093/nar/29.1.37

Paysan-Lafosse T, Blum M, Chuguransky S, Grego T, Pinto BL, Salazar GA, Bileschi ML, Bork P, Bridge A, Colwell L, et al.
InterPro in 2022.
Nucleic Acids Research 51(D1): D418-D427, 2023.
doi:10.1093/nar/gkac993
```

### Pfam
```
Bateman A, Birney E, Durbin R, Eddy SR, Howe KL, Sonnhammer EL.
The Pfam protein families database.
Nucleic Acids Research 28(1): 263-266, 2000.
doi:10.1093/nar/28.1.263

Mistry J, Chuguransky S, Williams L, Qureshi M, Salazar GA, Sonnhammer ELL, Tosatto SCE, Paladin L, Raj S, Richardson LJ, et al.
Pfam: The protein families database in 2021.
Nucleic Acids Research 49(D1): D412-D419, 2021.
doi:10.1093/nar/gkaa913
```

### PROSITE
```
Bairoch A.
PROSITE: a dictionary of sites and patterns in proteins.
Nucleic Acids Research 19(Suppl): 2241-2245, 1991.
doi:10.1093/nar/19.suppl.2241

Bairoch A, Bucher P, Hofmann K.
The PROSITE database, its status in 1997.
Nucleic Acids Research 25(1): 217-221, 1997.
doi:10.1093/nar/25.1.217

Sigrist CJA, de Castro E, Cerutti L, Cuche BA, Hulo N, Bridge A, Bougueleret L, Xenarios I.
New and continuing developments at PROSITE.
Nucleic Acids Research 41(D1): D344-D347, 2013.
doi:10.1093/nar/gks1067
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

## Support

- **Documentation:** Check this README and inline help tooltips
- **Issues:** Report problems on [GitHub Issues](https://github.com/ultimafounding/PicoMol/issues)
- **Feedback:** Use the built-in feedback dialog (Help → Send Feedback)

---

*PicoMol: Making molecular visualization and bioinformatics accessible to everyone.*
