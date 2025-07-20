# Changelog

All notable changes to PicoMol will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.4] - 2025-07-20
### Added
- **Comprehensive Motif and Domain Analysis:** New dedicated tab for protein functional annotation
  - **InterPro Integration:** Official EBI InterPro API integration for comprehensive domain annotation
    - Includes Pfam, SMART, PROSITE, CDD, PRINTS, and many more databases
    - GO terms and pathway annotation
    - Detailed functional site identification
  - **PROSITE Motif Search:** Official ExPASy PROSITE ScanProsite API integration
    - Real-time motif pattern recognition
    - Intelligent local pattern fallback for reliability
    - Filtering options for high-probability motifs
    - Comprehensive motif categorization
  - **Advanced Visualization:** Interactive domain and motif visualization
    - Dynamic layout with intelligent collision detection
    - Multi-track display for complex proteins
    - Scalable visualization adapting to sequence length
  - **Export Capabilities:** Multiple output formats including detailed HTML reports
  - **Search Options:** Individual database searches or comprehensive "Search All" functionality

### Enhanced
- **Welcome Dialog:** Updated with new motif analysis features and improved tips
- **About Dialog:** Added InterPro and PROSITE citations and updated version information
- **Documentation:** Comprehensive updates to README with motif analysis features
- **User Interface:** Enhanced tooltips and help text throughout the application
- **API Integration:** Robust web service connections with intelligent fallback systems

### Technical Improvements
- **Hybrid API Approach:** API-first with local pattern fallback for maximum reliability
- **Advanced Parsing:** Robust parsing of various API response formats (XML, tab-separated, HTML)
- **Multi-threading:** Background analysis prevents UI freezing during searches
- **Error Handling:** Comprehensive error recovery and user feedback
- **Result Caching:** Efficient storage and retrieval of analysis results

### Fixed
- **PROSITE API Integration:** Resolved parsing issues with tab-separated response format
- **Result Display:** Fixed table population and result visualization
- **Sequence Input:** Enhanced placeholder text and input validation
- **UI Responsiveness:** Improved handling of long-running API requests

## [1.0.0] - 2025-01-15 (Planned)
### Added
- **Comprehensive Bioinformatics Suite:** Complete sequence analysis tools for proteins, DNA, and RNA
  - Protein analysis: molecular weight, isoelectric point, hydropathy (GRAVY), aromaticity, instability index
  - Secondary structure prediction with helix, turn, and sheet content
  - Charge analysis at different pH values (5, 7, 9)
  - Detailed amino acid composition with counts and percentages
  - Nucleic acid analysis: GC content, complement sequences, translation frames
  - Nucleotide composition in conventional order (A, T/U, G, C, N)
- **Enhanced FASTA Support:** Smart parsing and handling of sequence files
  - Multi-sequence FASTA file support with interactive selection dialog
  - Automatic sequence type detection (protein/DNA/RNA)
  - Direct FASTA text pasting with intelligent parsing
  - Sequence validation with helpful error messages
  - Cross-platform integration between structure viewer and analysis tools
- **Advanced User Interface:** Modern, customizable interface with comprehensive features
  - Theme system with multiple built-in themes (System Default, Light, Dark, Blue, Green)
  - Comprehensive preferences dialog for appearance and behavior customization
  - Enhanced tabbed interface with organized tool sections
  - Improved error handling with actionable suggestions
  - Responsive design that adapts to different screen sizes
- **Robust Analysis Engine:** Dual-mode analysis system for maximum compatibility
  - Full-featured analysis with Biopython integration
  - Fallback basic analysis mode for systems without Biopython
  - Multi-threaded analysis to prevent UI freezing
  - Intelligent partial codon handling with user notifications
  - Warning suppression for better user experience

### Enhanced
- **Complete BLAST Integration:** Full NCBI-style interface with all BLAST types
  - BLASTP, BLASTN, BLASTX, TBLASTN, and TBLASTX support
  - Online BLAST searches connecting directly to NCBI servers
  - Comprehensive results parsing and formatting
  - Multiple database options (nr, nt, RefSeq, etc.)
- **3D Molecular Visualization:** Improved structure viewing and interaction
  - Enhanced drag-and-drop support for PDB and ENT files
  - Better sequence display integration with structure viewer
  - Improved screenshot functionality
  - Enhanced undo/redo system for visualization changes
- **Project Organization:** Better code structure and maintainability
  - Organized `blast_utils` package for BLAST functionality
  - Modular `bioinformatics_tools` for sequence analysis
  - Separate `theme_manager` and `preferences` modules
  - Comprehensive documentation and inline help

### Changed
- **Version Numbering:** Updated to semantic versioning (1.0.0)
- **Interface Layout:** Reorganized tabs for better workflow
  - 3D Viewer tab for molecular visualization
  - Bioinformatics tab with comprehensive analysis tools
  - BLAST tab with all search types
  - Removed placeholder tabs, replaced with functional tools
- **Documentation:** Comprehensive updates to README, about dialog, and help text
- **Welcome Dialog:** Enhanced with feature overview and modern design
- **Error Handling:** Improved user feedback with detailed error messages

### Fixed
- **Biopython Compatibility:** Resolved deprecation warnings and version issues
  - Fixed `get_amino_acids_percent` deprecation warning
  - Improved GC content function compatibility across Biopython versions
  - Better handling of partial codon translation warnings
- **UI Stability:** Enhanced interface reliability and responsiveness
  - Fixed layout clearing issues in results display
  - Improved widget cleanup and memory management
  - Better handling of theme switching and preferences
- **Sequence Processing:** More robust sequence handling and validation
  - Enhanced FASTA parsing with better error recovery
  - Improved sequence type detection algorithms
  - Better handling of edge cases in sequence analysis

### Technical Improvements
- **Code Architecture:** Modular design with clear separation of concerns
- **Performance:** Optimized analysis algorithms and UI responsiveness
- **Compatibility:** Better cross-platform support and dependency handling
- **Testing:** Enhanced error handling and validation throughout the application

## [0.0.3] - 2025-07-18
### Added
- Complete BLAST integration with full NCBI-style interface
  - BLASTP, BLASTN, BLASTX, TBLASTN, and TBLASTX support
  - Online BLAST searches connecting directly to NCBI servers
  - Comprehensive results parsing and formatting
  - Sequence validation and FASTA file support
  - Multiple database options (nr, nt, RefSeq, etc.)
- Organized project structure with `blast_utils` package
- Enhanced tabbed interface for better feature organization
- Sequence display above 3D viewer for better visibility
- Color picker for uniform color scheme selection
- Requirements.txt and setup.py for easier installation
- Comprehensive documentation updates

### Changed
- Restructured codebase with organized `blast_utils` package
- Updated import statements for new package structure
- Improved README with detailed feature descriptions
- Enhanced error handling and user feedback
- Better project organization and file structure

### Removed
- Legacy BLAST files and duplicate code
- Large HTML reference files (blastn.html, etc.)
- Unnecessary development documentation files
- Python cache files and temporary artifacts

### Fixed
- Resolved BLAST function import errors
- Fixed sequence display issues
- Improved tab switching functionality
- Enhanced undo/redo system stability

## [0.0.2] - 2025-07-15
### Added
- Initial release of PicoMol
- Basic protein structure visualization
- PDB file loading and fetching
- Basic visualization controls (color schemes, representations)
- Screenshot functionality
- Undo/redo support