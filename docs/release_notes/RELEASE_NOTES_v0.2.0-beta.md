# PicoMol v0.2.0-beta Release Notes
## Advanced Sequence Analysis Suite

**Release Date:** March 22, 2026  
**Version:** 0.2.0-beta  
**Previous Version:** 0.1.0 (Production Ready)

---

## 🎉 Major Feature Release

Version 0.2.0-beta represents a significant milestone in PicoMol's evolution, introducing professional-grade advanced sequence analysis capabilities that transform the application from a molecular visualization tool into a comprehensive bioinformatics platform. This release builds upon the solid foundation of v0.1.0 with cutting-edge sequence analysis tools designed for researchers, educators, and molecular biologists.

---

## 🆕 What's New in v0.2.0-beta

### 🧬 **Advanced Sequence Analysis Suite**

#### **ORF (Open Reading Frame) Analysis**
- **Comprehensive ORF Detection:** Find open reading frames with full customization
  - Configurable minimum length (30-3000 bp)
  - Custom start and stop codons (default: ATG start, TAA/TAG/TGA stop)
  - Analysis of all 6 reading frames (3 forward, 3 reverse)
- **Detailed Results Display:**
  - Interactive results table with frame, position, and sequence information
  - Statistics panel with total ORFs, strand distribution, and length analysis
  - Sequence preview for each identified ORF
- **Research-Ready Output:** Export ORF data for downstream analysis

#### **Professional Primer Design Tool**
- **Target Region Selection:**
  - Interactive amplicon region definition
  - Real-time amplicon length calculation
  - Support for custom target regions
- **Advanced Primer Scoring Algorithm:**
  - Melting temperature (Tm) calculation with salt-adjusted formula
  - GC content optimization (40-60% preferred range)
  - Secondary structure analysis and avoidance
  - Self-complementarity and hairpin detection
- **Comprehensive Validation:**
  - SnapGene-style primer validation with detailed analysis
  - Primer pair analysis with dimer prediction
  - Quality scoring system (0-100 scale)
- **Export Capabilities:**
  - Copy primer sequences to clipboard
  - Export primer data to CSV format
  - Detailed validation reports

#### **Sequence Comparison Tools**
- **Identity Analysis:** Compare multiple sequences with precision
  - Longest common substring detection
  - Percentage identity calculation
  - Alignment visualization
- **Flexible Input Options:**
  - Load sequences from FASTA files
  - Compare with currently loaded sequences
  - Support for multiple sequence comparison

### 🎨 **Enhanced Visualization & Export**

#### **Professional Export System**
- **Multi-Format Support:**
  - **Map Export:** PNG, SVG, PDF with customizable resolution
  - **Sequence Export:** FASTA, GenBank, HTML, annotated text
- **Advanced Export Options:**
  - Customizable size and DPI settings (up to 1200 DPI)
  - Element selection (features, enzymes, scale, labels)
  - Color scheme customization
- **Live Preview:** Real-time preview before saving
- **Publication Quality:** High-resolution outputs suitable for academic publications

#### **SnapGene-Style Sequence Viewer**
- **Enhanced Enzyme Visualization:**
  - Colored enzyme boxes with professional styling
  - Connection lines for enzyme-site relationships
  - Interactive enzyme information display
- **Improved Feature Rendering:**
  - Background highlighting for genetic features
  - Professional labeling and annotation
  - Optimized for large sequences
- **Performance Optimizations:**
  - Smooth scrolling for long sequences
  - Efficient memory usage
  - Responsive interface during complex operations

### 🔧 **Technical Improvements**

#### **PyQt6 Migration**
- **Complete Framework Update:** Migrated from PyQt5 to PyQt6
  - Enhanced performance on modern systems
  - Improved stability and compatibility
  - Future-proof codebase for continued development
- **API Updates:** All GUI components updated to PyQt6 standards
- **Cross-Platform Enhancement:** Better compatibility across Windows, macOS, and Linux

#### **Architecture Enhancements**
- **Modular Design:** Separate analysis widgets for maintainability
- **Background Processing:** Threading for complex calculations prevents UI freezing
- **Error Handling:** Comprehensive validation and user feedback
- **Extensibility:** Plugin-ready architecture for future tools

---

## 🔄 Major Changes Since v0.1.0

### **New Features Added:**
1. **Advanced Sequence Analysis Suite** - ORF detection, primer design, sequence comparison
2. **Professional Export System** - Multi-format export with live preview
3. **SnapGene-Style Visualization** - Enhanced sequence viewer with professional styling
4. **PyQt6 Migration** - Complete framework update for modern compatibility
5. **Cloning Workflow Tools** - Gibson assembly, Golden Gate, restriction cloning
6. **Virtual Gel Electrophoresis** - Simulate gel results for planning
7. **Enhanced Enzyme Tools** - Comprehensive enzyme selection and analysis

### **Technical Improvements:**
1. **Performance Optimization** - Better handling of large sequences
2. **Memory Management** - Efficient memory usage for complex analyses
3. **User Interface** - Responsive design with modern styling
4. **Error Handling** - Comprehensive validation and user feedback
5. **Code Quality** - Modular architecture with improved maintainability

---

## 📁 New Project Structure

The project has been expanded to accommodate new features:

```
src/gui/
├── advanced_analysis.py          # ORF, primer design, sequence comparison
├── export_dialog.py             # Professional export system
├── snapgene_sequence_view.py    # Enhanced sequence visualization
├── enzyme_selection_dialog.py    # Enzyme selection tools
├── restriction_cloning_dialog.py # Restriction cloning workflow
├── gibson_assembly_dialog.py   # Gibson assembly workflow
├── golden_gate_dialog.py        # Golden Gate cloning workflow
├── virtual_gel_dialog.py       # Virtual gel simulation
└── [additional enhanced components...]
```

---

## 🚀 Installation & Upgrade

### **New Installation:**
```bash
git clone https://github.com/ultimafounding/PicoMol.git
cd PicoMol
pip install -r requirements.txt
python picomol.py
```

### **Upgrading from v0.1.0:**
1. **Backup your data:** Save any important sequences or analysis results
2. **Pull latest changes:** `git pull origin main`
3. **Update dependencies:** `pip install -r requirements.txt` (PyQt6 will be installed)
4. **Run PicoMol:** `python picomol.py`

**Note:** PyQt6 is now required. The upgrade process will automatically replace PyQt5 with PyQt6.

---

## 🎯 Target Users

Version 0.2.0-beta is designed for:

- **Molecular Biologists:** Professional primer design and cloning workflows
- **Researchers:** Advanced sequence analysis and comparison tools
- **Educators:** Comprehensive teaching tools with professional visualizations
- **Students:** User-friendly interface with advanced analysis capabilities
- **Bioinformaticians:** Export capabilities for downstream analysis

---

## 🔧 System Requirements

### **Minimum Requirements:**
- Python 3.7 or higher (Python 3.8+ recommended)
- 4GB RAM (8GB+ recommended for large sequences)
- 200MB free disk space
- Internet connection (for online features and API access)

### **Recommended Setup:**
- Python 3.9+
- 8GB+ RAM
- 1GB+ free disk space
- Stable internet connection
- Modern graphics card for optimal visualization

---

## 🐛 Beta Status & Known Issues

**This is a beta release** with advanced features that may contain:

- **New Features:** Advanced sequence analysis tools are in beta testing
- **Performance:** Large sequences (>100kb) may experience slower performance
- **Export:** Some export formats may need refinement

**Feedback Welcome:** Please report issues and suggestions on [GitHub Issues](https://github.com/ultimafounding/PicoMol/issues).

---

## 🔮 Future Roadmap

Planned for upcoming releases:
- **Batch Processing:** Analyze multiple sequences simultaneously
- **Advanced Cloning:** Additional cloning workflows and simulations
- **Plugin System:** Extensible architecture for custom analysis tools
- **Cloud Integration:** Save and share analyses online
- **Mobile Support:** Tablet compatibility for field work

---

## 🙏 Acknowledgments

### **Core Technologies:**
- **NGL.js:** High-performance molecular graphics (MIT License)
- **Biopython:** Comprehensive bioinformatics tools
- **PyQt6:** Cross-platform GUI framework
- **Matplotlib:** Advanced plotting and visualization
- **NCBI BLAST:** Sequence similarity search

### **Development:**
- **Developer:** Jack Magson
- **License:** GNU GPL v3.0
- **Repository:** [GitHub](https://github.com/ultimafounding/PicoMol)

---

## 📞 Support & Feedback

- **Documentation:** Complete user guides and technical documentation
- **Issues:** Report problems on [GitHub Issues](https://github.com/ultimafounding/PicoMol/issues)
- **Feedback:** Use the built-in feedback dialog (Help → Send Feedback)
- **Community:** Join discussions and share your experiences

---

## 📚 Citations

When using PicoMol v0.2.0-beta in your research, please cite:

**PicoMol:**
```
PicoMol: A comprehensive molecular visualization and bioinformatics suite.
Version 0.2.0-beta (2026). Available at: https://github.com/ultimafounding/PicoMol
```

**New Features in v0.2.0-beta:**
- Advanced sequence analysis suite (ORF detection, primer design, sequence comparison)
- Professional export system with multi-format support
- SnapGene-style sequence visualization
- PyQt6 migration for enhanced compatibility

---

**🎆 Welcome to the advanced sequence analysis era of PicoMol! 🎆**

*Making molecular visualization, bioinformatics, and sequence analysis accessible to everyone.*
