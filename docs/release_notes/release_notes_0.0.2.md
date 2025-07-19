# Release Notes: PicoMol v0.0.2

## New Features & Improvements

1. **Add a screenshot feature** to capture the current view of the protein structure.
2. **Provide tooltips and inline help** for all controls and options.
    - Tooltips have been added to all major controls in the UI, providing concise descriptions for each feature. Users now receive inline help by hovering over buttons, fields, and options.
3. **Implement a welcome screen or onboarding tutorial** for first-time users.
    - Welcome dialog now appears only if the user wants it, and remembers their preference. Minimal and non-intrusive.
4. **Add undo/redo support** for user actions.
    - Undo/redo now tracks all view-related changes, including custom color and background color, with a smart stack that only records real changes.
    - The stack is minimal and avoids redundant entries, providing a clean user experience.
    - **Note:** Undo/redo is still a work in progress. Minor issues may remain and feedback is welcome.
5. **Make file loading and saving as simple as possible** (with clear dialogs, recent files list, and drag-and-drop support).
    - Added drag-and-drop file opening with a visible overlay.
    - Added "Save Structure As..." and Recent Files menu for quick access.
6. **Improve error messages** to be actionable and user-friendly.
    - Added a reusable, user-friendly error dialog for all major error scenarios (missing files, failed downloads, dependency errors, etc.).
    - Status bar usage is now robust and no longer causes runtime crashes.
    - App is much more stable and provides actionable feedback instead of crashing.
9. **Allow users to reset the view or return to defaults easily.**
    - Added a "Reset View to Defaults" action in the Edit menu (Ctrl+R). This resets all visible viewer settings (representation, color scheme, spin, background color, custom color) to their defaults, with undo support.
10. **Add a “feedback” button** for users to quickly report issues or suggest improvements.
    - The feedback button now opens a dialog for user input. When submitted, it opens a pre-filled GitHub issue with the user's feedback in the body, making it easy for users to report bugs or suggestions directly to the developer. No local file is saved.
    - **Usage:** Click "Send Feedback" in the control panel, enter your feedback, and submit. Your browser will open a new GitHub issue page with your feedback ready to send.
11. **Provide a simple “about” dialog** with version info and credits.
    - About dialog added under the Help menu. Shows version (0.0.2), credits (Jack Magson), license, and full NGL.js citations. Styled consistently with the rest of the app.

---

## Bug Fixes

- **Partially Resolved:** Loading a file after inputting a non-valid file no longer results in a white screen. The application now closes with an error message, with future improvements planned for better error recovery.
- **Resolved:** The default color now correctly matches the color scheme box setting, as color scheme initialization has been moved to component load callbacks and representation changes.

---

## Known Issues

- Undo/redo functionality is still being refined and may have minor issues
- Some error states may still cause the application to close unexpectedly
- Performance may degrade with very large structures

---

For a full list of changes, see the changelog or previous ongoing_changes.md.
