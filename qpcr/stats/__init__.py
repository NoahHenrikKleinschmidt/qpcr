"""
The stats module provides a set of functions for statistical analysis of qPCR data.
This module defines stand-alone functions for performing statistical analysis as well as the `Evaluator` class
that handles linking to results and statistical evaluation.

Supported analysis include _ANOVA_ analysis over entire datasets or _multiple pairwise T-Tests_. 
For each analysis two modes are supported: ``groupwise`` and ``assaywise``. In `assaywise` mode, the tests are performed to compare the groups within each data column (assay) with each other.
In `groupwise` mode, the tests are performed to compare the data columns within each group of the dataframe overall.
In case of multiple T-Tests p-values are corrected for multiple testing using the `Benjamini-Hochberg` procedure. 
"""

from .Evaluator import Evaluator
from .func_api import *

__default_Evaluator__ = Evaluator()
"""The default Evaluator"""
