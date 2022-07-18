# New Features in Next Version

This is a more cleaned up version of the previous TODO.md

### `pipe` method comprehension
The `pipe` methods of `qpcr.Analyser`, `qpcr.Normaliser`, and the `qpcr.Filters`, as well as the `read` method of the `qpcr.DataReader` can now directly be fed with a `list`. 

### Extended interface for `qpcr.Assay` and `qpcr.Results` 
`qpcr.Assay` and `qpcr.Results` objects now allow direct item setting, getting and deleting on their dataframes. `qpcr.Results` can now be merged together using the `+` operator.

### Function API
There are now a number of functions that wrap the `qpcr` classes to allow even easier
data handling. Among them are:
- `qpcr.read(...)` to wrap a standard `qpcr.DataReader` and its `read` method
- `qpcr.analyse(...)` to wrap a standard `qpcr.Analyser` and its `pipe` method
- `qpcr.normalise(...)` to wrap a standard `qpcr.Normalaliser` and its `pipe method.

### New stats
The `qpcr.Results.stats()` dataframe now also includes `IQR` (by default, but adjustable to any two quantiles), and `CI` (assuming a normal distribution, default at 95% but also adjustable).

### New `__str__` representations
The main `qpcr` classes now all have a `__str__` method to allow easier user-interaction, where they display their dataframes as well as other information. 

### SampleReader drop
The `qpcr.SampleReader` was dropped from `qpcr`.

### `qpcr.Results.split` drop
The `split` method that previously generated a number of `qpcr.Results` objects from a single one based on its `_rel_` columns has been dropped. 

### Method changes
Some attributes such as the `qpcr.Assay.dCt` or `qpcr.Results.is_empty` are now properties and now longer callable methods.

### Bugfixes and Code Refactorizations
The huge `__init__` method was refactored into a proper `main` submodule. Also the `defaults` were refactored into their own submodule instead of being part of `_auxiliary`. Also the `id` policy has been changed to allow repetitive `.id(...)` calling without requiring the use of `.id_reset()`. 

