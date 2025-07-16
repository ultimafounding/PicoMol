# Changelog

All notable changes to PicoMol will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Tabbed interface for better organization of features
  - 3D Viewer: Main protein visualization with sequence display
  - Sequence Tools: Placeholder for future sequence analysis features
  - BLAST: Placeholder for BLAST search functionality
  - Structure: Placeholder for structural analysis tools
- Sequence display moved above the 3D viewer for better visibility
- Improved UI with better organization and visual hierarchy
- Color picker for uniform color scheme selection
  - Appears only when "uniform" color scheme is selected
  - Click to choose from system color dialog
  - Visual feedback of current color

### Changed
- Restructured main window layout to accommodate new tabbed interface
- Updated sequence display styling for better readability
- Adjusted sequence display height to 75px for better space efficiency
- Improved window resizing behavior

### Fixed
- Fixed sequence display not showing up after structure loading
- Resolved tab switching issues
- Fixed undo/redo functionality with new UI elements

## [0.0.2] - 2025-07-15
### Added
- Initial release of PicoMol
- Basic protein structure visualization
- PDB file loading and fetching
- Basic visualization controls (color schemes, representations)
- Screenshot functionality
- Undo/redo support
