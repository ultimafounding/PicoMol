# Changelog

All notable changes to PicoMol will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2025-09-14
### Added
- **Advanced Sequence Analysis Suite:** Comprehensive new analysis tools accessible via the sequence viewer
  - **ORF (Open Reading Frame) Analysis:** Find and analyze open reading frames with customizable parameters
    - Configurable minimum length, start/stop codons
    - Analysis of all 6 reading frames (forward and reverse)
    - Detailed statistics and results table with sequence preview
  - **Primer Design Tool:** Professional-grade PCR primer design with enhanced features
    - Target region selection with amplicon length calculation
    - Advanced primer scoring algorithm considering Tm, GC content, secondary structures
    - SnapGene-style primer validation with detailed analysis
    - Export capabilities (clipboard, CSV) and primer pair analysis
    - Comprehensive primer validation including self-complementarity analysis
  - **Sequence Comparison:** Tools for comparing multiple sequences
    - Identity analysis with longest common substring detection
    - Support for loading sequences from files or current viewer
    - Dot plot visualization (planned)

- **Enhanced Export System:** Professional export capabilities for sequences and maps
  - **Multi-format Export:** PNG, SVG, PDF for maps; FASTA, GenBank, HTML for sequences
  - **Customizable Export Options:** Size, DPI, elements to include, color schemes
  - **Preview Functionality:** Live preview of exports before saving
  - **Advanced Formatting:** Custom titles, metadata, font options
  - **Professional Quality:** High-resolution outputs suitable for publications

- **Improved Sequence Visualization:** Enhanced sequence viewer with SnapGene-style features
  - **Restriction Site Highlighting:** Visual enzyme boxes with colored backgrounds
  - **Enhanced Enzyme Display:** Professional enzyme labeling with connection lines
  - **Better Feature Rendering:** Improved feature background highlighting
  - **Optimized Performance:** Better handling of large sequences and complex annotations

### Enhanced
- **User Interface:** Advanced analysis tools integrated seamlessly into the sequence viewer
- **Code Architecture:** Modular design with separate analysis widgets for maintainability
- **Documentation:** Updated to reflect new advanced analysis capabilities

### Technical Improvements
- **Threading:** Background analysis workers prevent UI freezing during complex calculations
- **Error Handling:** Comprehensive validation and user feedback for all analysis tools
- **Performance:** Optimized algorithms for large sequence analysis
- **Extensibility:** Plugin-ready architecture for future analysis tools

## [0.1.0] - 2025-01-20
### Major Milestone Release - Production Ready

### Added
- **Production Quality Codebase:** Comprehensive testing and error handling throughout the application
- **Stable API Integration:** Robust InterPro and PROSITE integration with intelligent fallback systems
- **Performance Optimizations:** Enhanced handling of large protein structures and complex analyses
- **Comprehensive Documentation:** Updated README, technical documentation, and user guides
- **Version Management:** Proper semantic versioning and release management

### Enhanced
- **User Experience:** Improved interface responsiveness and error messaging
- **Data Export:** Refined export functionality with consistent formatting across all output types
- **Visualization:** Optimized 3D rendering and graph generation performance
- **Code Quality:** Refactored codebase with improved modularity and maintainability

### Fixed
- **Stability Issues:** Resolved various edge cases and improved application stability
- **Memory Management:** Optimized memory usage for large structure analysis
- **Cross-Platform Compatibility:** Enhanced compatibility across Windows, macOS, and Linux
- **API Reliability:** Improved error handling for network requests and API timeouts

### Technical Improvements
- **Modular Architecture:** Better separation of concerns and code organization
- **Error Recovery:** Comprehensive error handling with user-friendly messages
- **Performance Monitoring:** Improved performance tracking and optimization
- **Documentation:** Complete technical documentation and API references

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

## [0.0.5] - 2025-07-24
### Added
- **Comprehensive Data Export:** All export formats (HTML, JSON, CSV, Excel, and Text) now include detailed analysis results, including all tables and graphs.
- **Enhanced HTML Reports:** HTML exports now feature improved layouts, higher-resolution graphs, and better styling to prevent overlapping elements.

### Changed
- **Default Export Format:** The default export option is now HTML, providing a more user-friendly and comprehensive report out of the box.

### Fixed
- **Graph Rendering:** Fixed issues with overlapping elements in exported graphs, particularly the amino acid composition pie chart.
- **Data Consistency:** Ensured that all calculated data, including detailed Ramachandran, bond length, bond angle, and B-factor tables, is consistently included in all export formats.

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