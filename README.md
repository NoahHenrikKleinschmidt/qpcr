<!-- # <img src="https://user-images.githubusercontent.com/89252165/153070064-4d3fb42e-a5f9-40fd-b856-755d58a52687.svg" width="32"> qpcr -->
# <img src="./docs/source/qpcr_tiny.svg" width="25"> qpcr

### A python module to analyse qPCR data on single-datasets or high-throughput

[![DOI](https://zenodo.org/badge/398244987.svg)](https://zenodo.org/badge/latestdoi/398244987)
[![Generic badge](https://img.shields.io/badge/made_for-qPCR-yellow.svg)](https://shields.io/)
[![PyPI version](https://badge.fury.io/py/qpcr.svg)](https://badge.fury.io/py/qpcr)
[![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](http://qpcr.readthedocs.io/?badge=latest)
[![Downloads](https://static.pepy.tech/personalized-badge/qpcr?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=Downloads)](https://pepy.tech/project/qpcr)



This project presents a python package designed to facilitate the analysis of qPCR data through established `Delta-Delta-Ct` analysis. To that end, this module provides a set of processing classes that may be assembled into a fully-fledged analysis pipeline, starting from raw `Ct` values stored in `csv` or `excel` files all the way to finished visualisations of the results. 

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

#### Customisability
A technical note at this point. At it's core `qpcr` offers very versatile data manipulation through two processing classes called the `qpcr.Analyser` and the `qpcr.Normaliser`. The `qpcr.Analyser` performs actions on a single qPCR datasets / assay stored in a `qpcr.Assay` object (the central data storage unit of the `qpcr` module). 
It is used to perform `Delta-Ct` computation. However, the precise function that it _applies_ to the single Assay is costumisable, so there is no restriction as such to what the Analyser will do to an Assay. 
On the other hand the `qpcr.Normaliser` will perform actions on a single Assay using data from a second Assay. It is used to perform normalisation of assays-of-interst against normaliser-assays. However, here too the precise function that is applied to the Assay is costumisable.  

### Example usage
A very simple analysis may start from a single `excel` file containing data of several qPCR assays, some of which are "assays-of-interest" and some of which are "normalisers" such as 28S rRNA. If we have performed a tiny bit of pre-processing on our datafile, our analysis may look like this:

```python
import qpcr
from qpcr.Pipes import ddCt

myfile = "todays_qPCR_run.xlsx"

# read the datafile
assays, normalisers = qpcr.read(  myfile, 
                                  multi_assay = True, 
                                  decorator = True 
                                )

# we can use a pre-defined default Delta-Delta-Ct pipeline
pipe = ddCt()
pipe.link( assays = assays, normalisers = normalisers )
pipe.run()

# at this point we can save our results to a file
results = pipe.results()
results.save( "./MyResults/" )

# and generate a preview using some custom settings
colors = ["xkcd:pastel blue", "xkcd:sapphire", "xkcd:rose pink", "xkcd:raspberry"]

fig = results.preview( color = colors, edgecolor = "black" )
```
![colorful](https://user-images.githubusercontent.com/89252165/158015384-d26fcfec-0ad6-44bc-a771-35a5dd43a380.png)


### Getting started
A set of basic introductory tutorials is available as `jupyter notebooks` in the [Examples](https://github.com/NoahHenrikKleinschmidt/qpcr/tree/main/Examples) directory in this repository. For more information about the API, checkout the documentation on [read-the-docs](https://qpcr.readthedocs.io/en/latest/).

### Qupid Web App
In case you prefer a graphical-user-interface, `qpcr` offers the *Qupid* web app built with `streamlit`. Qupid offers the bulk of qpcr's main features with costumizability to some degree. Qupid is openly available [via streamlit](https://share.streamlit.io/noahhenrikkleinschmidt/qupid/main/src/main.py). 

#### Citation
Kleinschmidt, N. (2022). qpcr - a python module for easy and versatile qPCR data analysis for small-scale datasets and high-throughput (Version 4.0.0) [Computer software]. https://github.com/NoahHenrikKleinschmidt/qpcr.git

