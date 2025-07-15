- Loading a file after inputting a non-valid file results in a failure and a white screen.

FIXED:
- Default colour is not what is set in the colour scheme box (Fixed in commit [commit_hash])
  - Issue: Color scheme wasn't properly initialized after loading new structures
  - Fix: Moved color scheme initialization to component load callback and representation changes
  - Status: Resolved