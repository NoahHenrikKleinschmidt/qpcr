# <img src="https://user-images.githubusercontent.com/89252165/153070064-4d3fb42e-a5f9-40fd-b856-755d58a52687.svg" width="30"> qpcr

### A python module to analyse qPCR data on single-datasets or high-throughput

[![DOI](https://zenodo.org/badge/398244987.svg)](https://zenodo.org/badge/latestdoi/398244987)


This project presents a python module designed to facilitate the analysis of qPCR data through established `Delta-Delta-Ct` analysis. To that end, this module provides a set of processing classes that may be assembled into a fully-fledged analysis pipeline, starting from raw `Ct` values stored in `csv` files all the way to finished visualisations of the results. 

User-friendliness and quick and easy workflows were of primary concern during development, so that this module is suitable for both non-experienced users as well as veteran coders. Virtually all steps within the analysis workflow are customizible to allow a streamlined analysis of any given dataset.

The exported results are formatted to be readily imported into graphing software. However, this module also includes methods to readily generate both static and interactive figures for rapid data exploration. 

### Installation
This module can directly be installed via `pip`.
```
pip install qpcr
```

### Example usage
To facilitate data analysis, common workflows have been implemented in pre-defined `pipelines` that allow for quick data analysis with minimal user effort. An example analysis using the pre-defined `BasicPlus` pipeline:

```python

from qpcr.Pipes import BasicPlus
from qpcr.Plotters import PreviewResults

# get our datafiles
normaliser_files = [
    "../Example Data/28S.csv",
    "../Example Data/actin.csv"
]

sample_files = [
    "../Example Data/HNRNPL_nmd.csv",
    "../Example Data/HNRNPL_prot.csv",
    "../Example Data/SRSF11_nmd.csv",
    "../Example Data/SRSF11_prot.csv",
]

# define our experimental parameters
reps = 6
group_names = ["WT-", "WT+", "KO-", "KO+"] 

# setting up the pipeline
pipeline = BasicPlus()
pipeline.save_to("./Example Results")

pipeline.replicates(reps)
pipeline.names(group_names)

# set up a preview of results
preview = PreviewResults(mode = "static")
preview.params( # setting some custom style
                color = "xkcd:sandy yellow", 
                edgecolor = "black", 
                edgewidth = 1, 
                figsize = (8,4)
            )
pipeline.add_plotters(preview)

# feed in our data
pipeline.add_assays(sample_files)
pipeline.add_normalisers(normaliser_files)

# run the pipeline
pipeline.run()

# and at this point the results are already saved...
```

![](https://github.com/NoahHenrikKleinschmidt/qpcr/blob/main/Examples/Example%20Results/PreviewResults_1.jpg)


### Getting started
A set of basic introductory tutorials is available as `jupyter notebooks` in the [Examples](https://github.com/NoahHenrikKleinschmidt/qpcr/tree/main/Examples) directory in this repository. For more information about the API, checkout the documentation at the [github-pages](https://noahhenrikkleinschmidt.github.io/qpcr/index.html).


#### Note
This is the second version of this module but the first to be properly released. Architecture and functionality have been fundamentally re-developed since Version 1 (which can still be accessed via the `legacy_V1` branch of the project repository).
Development of Version 1 is discontinued.

#### Citation
Kleinschmidt, N. (2022). qpcr - a python module for easy and versatile qPCR data analysis for small-scale datasets and high-throughput (Version 2.1.0) [Computer software]. https://github.com/NoahHenrikKleinschmidt/qpcr.git

