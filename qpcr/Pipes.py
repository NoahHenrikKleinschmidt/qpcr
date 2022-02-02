"""
This module contains a set of common lightweight wrappers to perform simple and non-customised DeltaDeltaCt analyses.
This is designed for lightweight users that do not wish to employ specialised pipelines.
"""

import qpcr.__init__ as qpcr
import matplotlib.pyplot as plt
import pandas as pd 
import statistics as stats
import qpcr._auxiliary.warnings as wa
import qpcr._auxiliary as aux
import qpcr.Plotters as Plotters
import qpcr.Filters as Filters
import re
import os 
import difflib

class Pipeline(qpcr.SampleReader):
    """
    This is the basic template class for qpcr Pipelines. 
    It contains a set of basic preliminary methods
    that ensure that elementary required inputs are provided.

    Note
    ----
    The simplest implementation of this `Pipeline` template is the `Basic` pipeline.
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
        Will run the pipeline provided that at least minimal inputs have been provided (i.e. Assays and Normalisers, as well as replicate specifics have been provided).
        This is a wrapper, the actual functional core is defined in the method `self._run()`. 
        To implement your own pipeline, make sure to define your own `_run()`, redefine `run()`
        as you require.

        Parameters
        ----------
        **kwargs
            Any additional keyword arguments that will be passed to the actual `_run()` method.
        """
        if self._replicates is None:
            wa.HardWarning("Pipeline:no_reps")
        elif self._Normalisers == [] or self._Assays == []:
            wa.HardWarning("Pipeline:no_data")

        self._run(**kwargs)
    
    def save_to(self, directory:str):
        """
        Set the location where to save result files

        Parameters
        ----------
        directory : str
            A directory to save the results to.
        """
        self._save_to = directory
        # if the directory does not yet exist, we make it
        if not os.path.exists(self._save_to):
            os.mkdir(self._save_to)

    def get(self, kind="stats"):
        """
        Returns
        -------
        data 
            A pandas dataframe either in replicate version `kind="df"`
            or in stats version `kind="stats"` (default). Alternatively,
            a qpcr.Results object can be returned using `kind="obj"`.
        """
        if kind == "stats":
            return self._stats_df
        elif kind == "df":
            return self._df
        elif kind == "obj":
            return self._Results
    
    def link(self, samples:(list or str)):
        """
        Links new sample assays to the pipline, either replacing old ones or keeping them, 
        depending on `softlink()` settings.

        Parameters
        ----------
        samples : list or str
            A `list` of filepaths to raw datafiles of sample assays, or a directory (`str`) where these are stored.
        """
        if self._softlink:
            self.prune()
        self.add_assays(samples)

    def prune(self, assays = True, results = True, normalisers = False):
        """
        Will clear assays, results, and/or normalisers

        Parameters
        ----------
        assays : bool
            Will clear any sample assays in the pipeline if True (default).
        
        results : bool
            Will clear any computed results in the pipline if True (default).
        
        normalisers : bool
            Will clear any normalisers in the pipline if True (default is False).
        """
        if assays: self._Assays = []
        if normalisers: self._Normalisers = []
        if results: 
            self._df = None
            self._stats_df = None
            self._Results = None

    def add_normalisers(self, normalisers):
        """
        Adds normalisers (filepaths) (keeping any already present)
        
        Parameters
        ----------
        normalisers : list or str
            A `list` of filepaths to raw datafiles of normaliser assays, or a directory (`str`) where these are stored.
        """
        normalisers = self._from_directory(normalisers)
        self._Normalisers.extend(normalisers)
    
    def add_assays(self, samples):
        """
        Adds sample assays (filepaths) (keeping any already present)

        Parameters
        ----------
        samples : list or str
            A `list` of filepaths to raw datafiles of sample assays, or a directory (`str`) where these are stored.
        """
        samples = self._from_directory(samples)
        self._Assays.extend(samples)

    def softlink(self, bool = None):
        """
        If `softlink = True`, then `link_assays()` will 
        prune any previous assays. Otherwise, it will 
        add new ones and keep old ones.

        Parameters
        ----------
        bool
            Set to False to disable `softlinking` (default is True).
        """
        if bool is None:
            return self._softlink
        else:
            self._softlink = bool

    def _from_directory(self, files):
        """
        Checks if a directory was provided for assays / normalisers and returns a list of all contained files if so.
        Otherwise it just returns the list of files.
        This is used for add_assays and add_normalisers
        """
        if isinstance(files, str):
            if os.path.isdir(files):
                norms = os.listdir(files)
                norms = [os.path.join(files,n) for n in norms]
                return norms
        else:
            return files

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
    The same as the Basic Pipeline, but has the option to integrate Plotters and Filters.
    """
    def __init__(self):
        super().__init__()
        self._Plotters = []
        self._Figures = []
        self._Filters = []
    
    def add_plotters(self, *Plotters:object):
        """
        Adds already specified qpcr.Plotter objects to the Pipeline.
        
        Parameters
        ----------
        *Plotters
            Any number of qpcr.Plotters.Plotter objects.
        """
        self._Plotters.extend(Plotters)
    
    def add_filters(self, *Filters:object):
        """
        Adds already specified qpcr.Filter instances to the Pipeline.
        
        Parameters
        ----------
        *Filters
            Any number of `qpcr.Filters.Filter` objects.
        """
        self._Filters.extend(Filters)

    def Figures(self):
        """
        Returns
        -------
        list
            A list of all Figures generated by the pipeline's Plotters.
        """
        return self._Figures

    def _run(self):
        """
        The automated standard DeltaDeltaCt pipeline. 
        Also produces applies and Plotters to produce figures...
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

            for filter in self._Filters:
                norm = filter.pipe(norm)

            norm = analyser.pipe(norm)
            normalisers.append(norm)
        normaliser.link(normalisers = normalisers)

        # analyse sample assays
        for sample in self._Assays:
            _sample = reader.read(sample)

            for filter in self._Filters:
                _sample = filter.pipe(_sample)

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

        # plot replicate overview
        if self._save_to is not None:
            for filter in self._Filters:

                # add report location if none was specified...
                if filter.report() is None: 
                    filter.report(self._save_to)
                
                figs = filter.plot()
                self._Figures.extend(figs)

        # plot results
        for plotter in self._Plotters:
            plotter.link(self._Results)
            fig = plotter.plot()
            self._Figures.append(fig)

            if self._save_to is not None:
                filename = self._make_figure_filename(plotter)
                plotter.save(filename)

    def _make_figure_filename(self, plotter):
        """
        Increments a filename with a numeric counter...
        """
        num = 1
        suffix = plotter.suffix()
        while True:
            filename = os.path.join(self._save_to, f"{plotter.id()}_{num}.{suffix}")
            if not os.path.exists(filename):
                break
            num+=1
        return filename



if __name__ == "__main__":

    norm_files = ["Example Data/28S.csv", "Example Data/actin.csv"]
    sample_files = ["Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]

    norm_folder = "Example Data 2/normalisers"
    sample_folder = "Example Data 2/samples"

    groupnames = ["wt-", "wt+", "ko-", "ko+"]
    
    analysis = BasicPlus()
    analysis.save_to("Example Data 2")
    analysis.add_assays(sample_folder) # alternative: link() for iteratively linking new assays...
    analysis.add_normalisers(norm_folder)
    
    analysis.replicates(6)
    analysis.names(groupnames)

    iqr_filter = Filters.RangeFilter()
    # iqr_filter.report("Example Data 2")
    analysis.add_filters(iqr_filter)

    preview = Plotters.PreviewResults(mode = "interactive")
    analysis.add_plotters(preview)

    # now that pipeline is ready, we can run!
    analysis.run()

    # now we can get results!
    results = analysis.get(kind="df")
    print(results)
    
    
    # norm_folder = "Example Data 3/normalisers"
    # exp_folder = "Example Data 3/experimental/samples"
    # ctr_folder = "Example Data 3/experimental/samples"

    # groupnames = ["wt-", "wt+", "ko-", "ko+"]
    
    # analysis = BasicPlus(index = "Example Data 2/index.csv")
    # analysis.save_to("Example Data ")
    
    # print(analysis._index)
    # # analysis.set_experimental(assays = exp_folder, normalisers = norm_folder)
    # # analysis.set_control(assays = ctr_folder, normalisers = norm_folder)

    # # analysis.replicates(6)
    # # analysis.names(groupnames)

    # iqr_filter = Filters.RangeFilter()
    # # iqr_filter.report("Example Data 2")
    # analysis.add_filters(iqr_filter)

    # # preview = Plotters.PreviewResults(mode = "interactive")
    # # analysis.add_plotters(preview)

    # # # now that pipeline is ready, we can run!
    # # analysis.run()

    # # # now we can get results!
    # # results = analysis.get(kind="df")
    # # print(results)

    # # exit(0)
