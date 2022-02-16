# <img src="https://user-images.githubusercontent.com/89252165/153070064-4d3fb42e-a5f9-40fd-b856-755d58a52687.svg" width="32"> qpcr

### A python module to analyse qPCR data on single-datasets or high-throughput

[![DOI](https://zenodo.org/badge/398244987.svg)](https://zenodo.org/badge/latestdoi/398244987)

---

> NOTE: This is the development branch of the next planned release `v3.0.1`.
> Check out the `TODO.md` to check on new features and progress...
---


This project presents a python module designed to facilitate the analysis of qPCR data through established `Delta-Delta-Ct` analysis. To that end, this module provides a set of processing classes that may be assembled into a fully-fledged analysis pipeline, starting from raw `Ct` values stored in `csv` or `excel` files all the way to finished visualisations of the results. 

User-friendliness and quick and easy workflows were of primary concern during development, so that this module is suitable for both non-experienced users as well as veteran coders. Virtually all steps within the analysis workflow are customizible to allow a streamlined analysis of any given dataset.

The exported results are formatted to be readily imported into other analysis or graphing software.

### Installation
This module can directly be installed via `pip`.

```
pip install qpcr
```

### What does `qpcr` do?
The "core business" that `qpcr` was designed for is `Delta-Delta-Ct` analysis starting from raw Ct values. It offers automated processes to read datafiles, filter out outlier Ct values, compute Delta-Ct, normalise assays against one another, and visualise the results. Hence, `qpcr` offers a full suite for automated Delta-Delta-Ct analyses.
However, even if your analysis is not going to be Delta-Delta-Ct, you may wish to check out how the `qpcr.Readers` might help you faciliate your workflow, as they are streamlined to read diversely structured datafiles contianing both a single and multiple qPCR datasets. 

#### Customisibility
A technical note at this point. At it's core `qpcr` offers very versatile data manipulation through two processing classes called the `qpcr.Analyser` and the `qpcr.Normaliser`. The `qpcr.Analyser` performs actions on a single qPCR datasets / assay stored in a `qpcr.Assay` object (the central data storage unit of the `qpcr` module). 
It is used to perform `Delta-Ct` computation. However, the precise function that it _applies_ to the single Assay is costumisible, so there is no restriction as such to what the Analyser will do to an Assay. 
On the other hand the `qpcr.Normaliser` will perform actions on a single Assay using data from a second Assay. It is used to perform normalisation of assays-of-interst against normaliser-assays. However, here too the precise function that is applied to the Assay is costumisible.  

### Example usage
To facilitate data analysis, common workflows have been implemented in pre-defined `pipelines` that allow for quick data analysis with minimal user effort. An example analysis using the pre-defined `BasicPlus` pipeline:

```python

from qpcr.Pipes import BasicPlus
from qpcr.Plotters import PreviewResults

# get our datafiles
normaliser_files = [
                    "./Examples/Example Data/28S.csv",
                    "./Examples/Example Data/actin.csv"
                ]

assay_files = [
                "./Examples/Example Data/HNRNPL_nmd.csv",
                "./Examples/Example Data/HNRNPL_prot.csv",
                "./Examples/Example Data/SRSF11_nmd.csv",
                "./Examples/Example Data/SRSF11_prot.csv",
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
colors = ["xkcd:pastel blue", "xkcd:sapphire", "xkcd:rose pink", "xkcd:raspberry"]
preview.params( # setting some custom style
                color = colors, 
                edgecolor = "black", 
                edgewidth = 1, 
                figsize = (8,4)
            )
pipeline.add_plotters(preview)

# feed in our data
pipeline.add_assays(assay_files)
pipeline.add_normalisers(normaliser_files)

# run the pipeline
pipeline.run()

# and at this point the results are already saved...
```

![](https://github.com/NoahHenrikKleinschmidt/qpcr/blob/dev_next/Examples/Example%20Results/colorful.jpg)


### Getting started
A set of basic introductory tutorials is available as `jupyter notebooks` in the [Examples](https://github.com/NoahHenrikKleinschmidt/qpcr/tree/main/Examples) directory in this repository. For more information about the API, checkout the documentation at the [github-pages](https://noahhenrikkleinschmidt.github.io/qpcr/index.html).


#### Note
This is the second version of this module but the first to be properly released. Architecture and functionality have been fundamentally re-developed since Version 1 (which can still be accessed via the `legacy_V1` branch of the project repository).
Development of Version 1 is discontinued.

#### Citation
Kleinschmidt, N. (2022). qpcr - a python module for easy and versatile qPCR data analysis for small-scale datasets and high-throughput (Version 3.0.1) [Computer software]. https://github.com/NoahHenrikKleinschmidt/qpcr.git

