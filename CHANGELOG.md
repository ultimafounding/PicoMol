# Changelog

## [Unreleased]
### Added
- Reusable, user-friendly error dialog via `show_error_dialog` in `ProteinViewerApp`.
- Improved error messages for missing dependencies, failed downloads, and file errors.

### Changed
- Status bar usage refactored to use `self.statusBar().showMessage(...)` throughout the app, preventing crashes.
- `self.port` is now reliably set in the `ProteinViewerApp` constructor, fixing attribute errors.
- All error dialogs now prevent the app from crashing and provide actionable feedback.
- Refactored and consolidated class structure, removed duplicate and misplaced assignments.

### Fixed
- Fixed `TypeError: 'QStatusBar' object is not callable` by removing incorrect assignment.
- Fixed `NameError: name 'port' is not defined` and `'ProteinViewerApp' object has no attribute 'port'` crashes.
- App no longer crashes on error dialog or status bar usage.

## [Older versions]
- See git history for prior changes.
