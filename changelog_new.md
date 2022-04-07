# New Features in Next Version

This is a more cleaned up version of the previous TODO.md

### New `qpcr.Calibrator` class
The new `qpcr.Calibrator` allows calculation of qPCR primer efficiencies based 
a dilution series. It is able to work either with an entire assay as input or with
just a subset of specially labeled replicates and can compute, store, save, and load and later assign pre-computed efficiencies from and to `qpcr.Assay` objects. 

### New `EfficiencyLines` Plotter 
Accompanying the `Calibrator` is a new plotter that can visualise the linear regression performed during efficiency computation, similar to how the `FilterSummary` plotter can summarise the workings of the `qpcr.Filters`. 

### Assays and Efficiency
Efficiency is now directly handled by the `qpcr.Assay` objects (`qpcr.Analyser` can still receive `.efficiency()` as before, but this will be dropped at some point...)

### New Plotting shortcuts
A `Plotters.PreviewResults` can now directly be called from `qpcr.Results` without having to import and setting up manually. The new `qpcr.Results.preview()` method calls on `PreviewResults` object to plot the stored data and returns the produced figure.

The same goes for a `Plotters.ReplicateBoxplot` which can be called from a `qpcr.Assay` directly via `qpcr.Assay.boxplot()`. 

Also, the `qpcr.Calibrator` is equipped with a `qpcr.Calibrator.plot()` method to 
call `Plotters.EfficiencyLines`. 