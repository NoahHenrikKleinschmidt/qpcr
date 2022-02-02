# qpcr
### A python module to analyse qPCR data on single-datasets or high-throughput

This project represents a python package that includes a number of functions useful for the analysis of qPCR data generated generated by Qiagen RotorGene® 
and is taylored to work with the Excel Spreadhseet exported from this device (or, more precisely, a `csv`-copy of the same). However, any `csv`-formatted data should be valid input.

User friendliness and quick and easy workflows were of primary concern during development. The exported results are formatted to be readily imported into graphing software. However, this module also includes methods to readily generate static or interactive figures of analysed results. 

### Version 2 is here
This is the second version of this module. Architecture and functionality have been fundamentally re-developed since the first version (which can still be accessed via the `legacy_V1` branch of this repository). The current Version 2 implements a series of data-processing classes that interact to form an analysis-pipeline from start to finish. 

### Getting started...
A set of basic introductory tutorials is available as `jupyter notebooks` in the `Examples` directory in this repository. 



<!-- 

## An Example: __qPCR Analysis__ can be so quick 'n easy :-)
```python
import qpcr.Analysis as qA
hnrnpl_nmd = "Example Data/HNRNPL_nmd.csv"
hnrnpl_prot = "Example Data/HNRNPL_prot.csv"
s28 = "Example Data/28S.csv"

groups = ["wt-", "wt+", "ko-", "ko+"]
result = qA.delta_deltaCt([s28, hnrnpl_nmd, hnrnpl_prot], 
                        replicates=6, normaliser="28S",
                        anchor="first", group_names=groups)

qA.preview_results(result)
```
![image](https://user-images.githubusercontent.com/89252165/130353861-98ca0083-383b-45e1-87e1-99f1ab99daf7.png)

> Do these figures not look like what you want? No problem, they're just instant previews, you can quickly generate csv files with your results that you can import into Graphpad Prism or whatever your favourite graphing software happens to be ;-)

### For a detailed introduction, check out the wiki page!

## Now available as a web-app that facilitates quick and versatile (and high-throughput) qPCR data analysis without any coding: <a href = "https://share.streamlit.io/noahhenrikkleinschmidt/qpcr-analyser/main/main.py"> qPCR Analysis Tool </a> -->
