"""
This module contains a set of common lightweight wrappers to perform simple and non-customised DeltaDeltaCt analyses.
This is designed for lightweight users that do not wish to employ specialised pipelines.
"""

import __init__ as qpcr
import matplotlib.pyplot as plt
import pandas as pd 
import statistics as stats
# import qpcr.auxiliary as aux
# import qpcr.auxiliary.graphical.auxiliaries as gx
import difflib


class DeltaDeltaCt(qpcr.SampleReader):
    """
    Performs simple standardized DeltaDeltaCt analysis 
    based on two lists of files, one for normalisers, one for Sample Assays...
    This makes use of standard settings for qpcr.Analyser() and qpcr.Assay()
    which cannot be customized! For customization generate your own pipeline!
    """
    def __init__(self):
        super().__init__()
        self._Normalisers = []
        self._Assays = []
        self._save_to = None
        self._df = None
        self._stats_df = None
        self._Results = None

    def add_normalisers(self, normalisers:list):
        """
        Links normalisers (filepaths)
        """
        self._Normalisers.extend(normalisers)
    
    def add_assays(self, samples:list):
        """
        Links sample assays (filepaths)
        """
        self._Assays.extend(samples)

    def save_to(self, directory:str):
        """
        Set the location where to save result files
        """
        self._save_to = directory

    def get(self, kind="stats"):
        """
        Returns a pandas dataframe either in replicate version kind="df"
        or in stats version kind="stats" (default).
        """
        if kind == "stats":
            return self._stats_df
        elif kind == "df":
            return self._df
        elif kind == "obj":
            return self._Results
    

    def run(self):
        """
        The automated standard DeltaDeltaCt pipeline. 
        Returns a pandas dataframe and can save results to csv file 
        (if save_to() location has been defined. )
        """
        reader = qpcr.SampleReader()
        reader.replicates(self._replicates)
        if self._names is not None:
            reader.names(self._names)
        
        analyser = qpcr.Analyser()
        normaliser = qpcr.Normaliser()

        normalisers = []
        samples = []

        # analyse normalisers:
        for _normaliser in self._Normalisers:
            norm = reader.read(_normaliser)
            norm = analyser.pipe(norm)
            normalisers.append(norm)
        normaliser.link(normalisers = normalisers)

        # analyse sample assays
        for sample in self._Assays:
            _sample = reader.read(sample)
            _sample = analyser.pipe(_sample)
            samples.append(_sample)
        normaliser.link(samples = samples)

        normaliser.normalise()
        results = normaliser.get()

        self._Results = results
        self._df = results.get()
        self._stats_df = results.stats()

        if self._save_to is not None:
            results.save(self._save_to)


if __name__ == "__main__":

    norm_files = ["Example Data/28S.csv", "Example Data/actin.csv"]
    sample_files = ["Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]
    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    analysis = DeltaDeltaCt()
    analysis.add_assays(sample_files)
    analysis.add_normalisers(norm_files)
    analysis.replicates(6)
    analysis.names(groupnames)

    analysis.run()

    results = analysis.get()
    print(results)
    

    exit(0)
