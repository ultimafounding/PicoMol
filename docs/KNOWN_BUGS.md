PARTIALLY RESOLVED:
- Loading a file after inputting a non-valid file results in a failure and a white screen.
  - Issue: Application would show a white screen when an invalid PDB code was entered
  - Fix: Application now closes with an error message when an invalid PDB code is entered
  - Status: Partially Resolved - The application no longer shows a white screen, but closes instead. Future improvements could include better error recovery.

FIXED:
- Default colour is not what is set in the colour scheme box (Fixed in commit [commit_hash])
  - Issue: Color scheme wasn't properly initialized after loading new structures
  - Fix: Moved color scheme initialization to component load callback and representation changes
  - Status: Resolved