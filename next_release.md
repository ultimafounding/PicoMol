# Next Release: Planned Features and Improvements

This file tracks ongoing and planned changes for future releases of PicoMol.

## Planned for Next Release

- **Preferences/settings dialog**
    - Allow users to customize appearance and behavior.
    - Some improvements are underway, but full support is planned for the next release.
    - Feedback on accessibility is welcome!
- **Accessibility improvements**
    - Keyboard navigation, colorblind-friendly palettes, and better focus indicators.
    - Some improvements are underway, but full support is planned for the next release.
    - Feedback on accessibility is welcome!

- **Advanced Visualization Features**
    - Support for multiple viewports (split-screen viewing)
    - Customizable camera controls (zoom, pan, rotate)
    - Support for different lighting schemes
    - Ability to add labels to specific atoms/residues
    - Support for custom color gradients

- **Structural Analysis Tools**
    - Distance measurement between atoms/residues
    - Angle measurement between bonds
    - Hydrogen bond detection and visualization
    - Secondary structure identification and highlighting
    - RMSD calculation between structures

- **Sequence and Structure Alignment**
    - Sequence viewer synchronized with 3D structure
    - Structure alignment viewer
    - Support for multiple structure comparison
    - Highlight conserved regions

- **Export and Sharing**
    - Enhanced screenshot options (different resolutions, formats)
    - Video recording of molecular animations
    - Export structure in multiple formats (PDB, mmCIF, etc.)
    - Share current view state via URL

- **Performance Improvements**
    - Optimized loading of large structures
    - Better memory management
    - Improved rendering performance
    - Support for GPU acceleration

- **User Experience Enhancements**
    - Recent structures list with thumbnails
    - Favorite structures list
    - Quick access to frequently used features
    - Improved error handling and user feedback
    - Better support for high DPI displays

- **Data Integration**
    - Integration with UniProt data
    - Support for additional structure databases
    - Support for custom structure databases
    - Integration with sequence analysis tools

---

## Ideas/Requests

- _(Add your ideas, user requests, or feature suggestions here!)_

---

## Bug Fixes

- **Partially Resolved:** Loading a file after inputting a non-valid file no longer results in a white screen. The application now closes with an error message, with future improvements planned for better error recovery.
- **Resolved:** The default color now correctly matches the color scheme box setting, as color scheme initialization has been moved to component load callbacks and representation changes.

---

Feel free to update this list as development progresses or new ideas arise.
