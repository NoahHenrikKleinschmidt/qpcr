# New Features in Next Version

This is a more cleaned up version of the previous TODO.md

### New `qpcr.stats` module
The new `qpcr.stats` module defines the `qpcr.Evaluator` as well as some direct API functions to
perform statistical evaluations on the data. Currently supported actions are:
- `multiple pair-wise T-Tests` either group-wise (per assay / per column) or assay-wise (per group).

### `qpcr.Results.setup_cols` automated
It is not necessary anymore to call `setup_cols` manually before adding data using `add_Ct`, `add_dCt`, `add_ddCt`, they will do so if necessary.

### Bugfixes and refactorings...
The main functionality has not been altered.