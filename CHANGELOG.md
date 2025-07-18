# Changelog

All notable changes to PicoMol will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
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
