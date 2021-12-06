"""
This module contains a set of common lightweight wrappers to perform simple and non-customised DeltaDeltaCt analyses.
This is designed for lightweight users that do not wish to employ specialised pipelines.
"""

import __init__ as qpcr
import matplotlib.pyplot as plt
import pandas as pd 
import statistics as stats
import auxiliary.warnings as wa
import auxiliary as aux
import Plotters
import re
import os 
import difflib

class Pipeline(qpcr.SampleReader):
    """
    This is the basic template class for qpcr Pipelines. 
    It contains a set of basic preliminary methods
    that ensure that elementary required inputs are provided.
    """
    def __init__(self):
        super().__init__()
        self._Normalisers = []
        self._Assays = []
        self._save_to = None
        self._df = None
        self._stats_df = None
        self._Results = None
        self._replicates = None
        self._names = None
        self._softlink = True

    def run(self, **kwargs):
        """
        Will run the pipeline provided that at least minimal inputs have been provided.
        This is a wrapper, the actual pipeline is defined in the function self._run(). 
        To implement your own pipeline, make sure to define your own _run(), redefine run()
        as you require...
        """
        if self._replicates is None:
            wa.HardWarning("Pipeline:no_reps")
        elif self._Normalisers == [] or self._Assays == []:
            wa.HardWarning("Pipeline:no_data")

        self._run(**kwargs)
    
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
    
    def link(self, samples:list):
        """
        Links new sample assays to the pipline, either replacing old ones or keeping them, 
        depending on self.softlink()
        """
        if self._softlink:
            self.prune()
        self.add_assays(samples)

    def prune(self, assays = True, results = True, normalisers = False):
        """
        Will clear assays, results, and/or normalisers
        """
        if assays: self._Assays = []
        if normalisers: self._Normalisers = []
        if results: 
            self._df = None
            self._stats_df = None
            self._Results = None

    def add_normalisers(self, normalisers:list):
        """
        Adds normalisers (filepaths) (keeping any already present)
        """
        self._Normalisers.extend(normalisers)
    
    def add_assays(self, samples:list):
        """
        Adds sample assays (filepaths) (keeping any already present)
        """
        self._Assays.extend(samples)

    def softlink(self, bool = None):
        """
        If softlink = True, then .link_assays() will 
        prune any previous assays. Otherwise, it will 
        add new ones and keep old ones.
        """
        if bool is None:
            return self._softlink
        else:
            self._softlink = bool



    def _run(self, **kwargs):
        """
        This is the actual function that each custom Pipeline has to define...
        """
        print("self._run() is the actual function that each custom Pipeline has to define...")


class Basic(Pipeline):
    """
    Performs simple standardized DeltaDeltaCt analysis 
    based on two lists of files, one for normalisers, one for Sample Assays...
    This makes use of standard settings for qpcr.Analyser() and qpcr.Assay()
    which cannot be customized! For customization generate your own pipeline!
    """
    def __init__(self):
        super().__init__()
    
    def _run(self):
        """
        The automated standard DeltaDeltaCt pipeline. 
        Using default settings of qpcr.Analyser()
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

class BasicPlus(Basic):
    """
    The same as the Basic Pipeline, but has the option to integrate Plotters.
    """
    def __init__(self):
        super().__init__()
        self._Plotters = []
        self._Figures = []
    
    def add_plotters(self, *Plotters:object):
        """
        Adds already specified qpcr.Plotter instances to the Pipeline.
        """
        self._Plotters.extend(Plotters)
    
    def Figures(self):
        """
        Returns a list of all Figures generated
        """
        return self._Figures

    def _run(self):
        """
        The automated standard DeltaDeltaCt pipeline. 
        Also produces applies and Plotters to produce figures...
        """

        super()._run()
        
        for plotter in self._Plotters:
            
            plotter.link(self._Results)
            fig = plotter.plot()
            self._Figures.append(fig)

            if self._save_to is not None:
                filename = self._make_figure_filename(plotter)
                fig.savefig(
                   filename, res = 500
                )

    def _make_figure_filename(self, plotter):
        """
        Increments a filename with a numeric counter...
        """
        filename = os.path.join(self._save_to, f"{plotter.id()}_1.png")
        while os.path.exists(filename):
            old_num  = re.search(f"_(\d).png", filename).group(1)
            new_num = int(old_num) + 1
            filename = filename.replace(f"_{old_num}", f"_{new_num}")
        return filename
        

if __name__ == "__main__":

    norm_files = ["Example Data/28S.csv", "Example Data/actin.csv"]
    sample_files = ["Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]
    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    analysis = BasicPlus()
    analysis.link(sample_files)
    analysis.add_normalisers(norm_files)
    analysis.replicates(6)
    analysis.names(groupnames)
    analysis.save_to("qpcr/data")

    preview = Plotters.PreviewResults(mode = "static")
    preview.params(
        headers = ["NMD", "Prot"], 
        frame = False, 
        color = "green",
        show = False
    )

    analysis.add_plotters(preview)

    # now that pipeline is ready, we can run!
    analysis.run()

    # now we can get results!
    results = analysis.get()
    print(results)
    

    exit(0)
