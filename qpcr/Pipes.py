"""
This module contains a set of common lightweight wrappers to perform simple and non-customised DeltaDeltaCt analyses.
This is designed for lightweight users that do not wish to employ specialised pipelines.
"""

import qpcr.__init__ as qpcr
import matplotlib.pyplot as plt
import pandas as pd 
import statistics as stats
import qpcr.auxiliary.warnings as wa
import qpcr.auxiliary as aux
import qpcr.Plotters
import qpcr.Filters
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

    def add_normalisers(self, normalisers):
        """
        Adds normalisers (filepaths) (keeping any already present)
        """
        normalisers = self._from_directory(normalisers)
        self._Normalisers.extend(normalisers)
    
    def add_assays(self, samples):
        """
        Adds sample assays (filepaths) (keeping any already present)
        """
        samples = self._from_directory(samples)
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
        Adds already specified qpcr.Plotter instances to the Pipeline.
        """
        self._Plotters.extend(Plotters)
    
    def add_filters(self, *Filters:object):
        """
        Adds already specified qpcr.Filter instances to the Pipeline.
        """
        self._Filters.extend(Filters)

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
    

# Alright we tackle this as:
# - we need a list to store experimental runs (as folders or filelists...)
# - we need a list to store corresponding control runs ()

class DoubleNorm(BasicPlus):
    """
    This pipeline allows to normalise entire qpcr runs against each other. 
    It adds a second normaliser to the BasicPlus pipeline. Assays will be analysed by default DeltaDeltaCt,
    normalised run-internally against a normaliser (like Actin), 
    and eventually normalised against their corresponding assay from a control run. 

    NOTE: For this to work paired assays (such as experimental and ctrl) have 
          to have the same name (= id, which is by default the filename)!

    NOTE: THIS PIPELINE IS KINDA FLAWED; ALTHOUGH THIS IS PROBABLY DUE TO A BUG WITHIN THE qpcr.Normaliser() !!! 
    """

    def __init__(self):
        super().__init__()
        self._experimental_assays = [] 
        self._experimental_normalisers = []
        self._control_assays = []
        self._control_normalisers = []
        self._exp_id = "exp"
        self._ctrl_id = "ctrl"

    def link(self, f:dict, kind = "exp", **kwargs):
        """
        Link either experimental (kind = "exp") or control (kind = "ctrl") assays and/or normalisers
        using a dictionary.
        """
        
        if kind == "exp" and self._softlink:
            self.prune(**kwargs)
        elif kind == "ctrl" and self._softlink:
            self.prune(experimental=False, control=True, **kwargs)

        if kind == "exp":
            self.set_experimental(f["assay"], f["norm"])
        elif kind == "ctrl":
            self.set_control(f["assay"], f["norm"])
        

    def prune(self, experimental=True, control=False, assays=True, normalisers=True, **kwargs):
        """
        Clears experimental and/or control assays and/or normalisers
        """
        if experimental:
            if assays: self._experimental_assays = []
            if normalisers: self._experimental_normalisers = []
        
        if control: 
            if assays: self._control_assays  = []
            if normalisers: self._control_normalisers = []
    

    def set_experimental(self, assays:(list or str), normalisers:(list or str)):
        """
        Set experimental run files (list) or folder (string)
        """
        assays = self._from_directory(assays)
        self._experimental_assays.extend(assays)
        normalisers = self._from_directory(normalisers)
        self._experimental_normalisers.extend(normalisers)

    def set_control(self, assays:(list or str), normalisers:(list or str)):
        """        
        Sets control run files (list) or folder (string)
        """
        assays = self._from_directory(assays)
        self._control_assays.extend(assays)
        normalisers = self._from_directory(normalisers)
        self._control_normalisers.extend(normalisers)

    def ids(self, exp=None, ctrl=None):
        """
        Setup optional ids for experimental and/or control runs
        """
        if exp is not None: self._exp_id = exp
        if ctrl is not None: self._ctrl_id = ctrl
        return self._exp_id, self._ctrl_id


    def run(self):
        """
        The automated standard DeltaDeltaCt pipeline. 
        Using default settings of qpcr.Analyser()
        """
        reader = qpcr.SampleReader()
        reader.replicates(self._replicates)
        if self._names is not None:
            reader.names(self._names)
        
        analyser = qpcr.Analyser()
        first_normaliser = qpcr.Normaliser()
        second_normaliser = qpcr.Normaliser()


        # run controls
        control_results = self.__run(reader, analyser, first_normaliser, self._control_normalisers, self._control_assays)
        control_results.id(self.ids()[0])
        control_results.drop_rel()

        # run experimentals
        experimental_results = self.__run(reader, analyser, first_normaliser, self._experimental_normalisers, self._experimental_assays)
        experimental_results.id(self.ids()[1])
        # experimental_results.drop_rel()
        experimental_results = experimental_results.split()
        print(experimental_results)
        # now normalise experimentals against control
        second_normaliser.link(
                                samples = experimental_results, 
                                normalisers = [control_results]
                            )
        
        second_normaliser.normalise(dCt_col = "named", norm_col = "same")

        # alright, we next need to forward the information about dCt_col to the normaliser
        # but this needs to be done iteratively, meaning we have to allow each column to be used as dCt...

        # We can extract all the experimental results into separate dataframes within the normalisers self._Assay list
        # like this we can have the one normaliser df with the corresponding columns and a number of Assays (experimental) that we can use for the iteration within normaliser.normalise()
        # for this we need to implement a split function for Results() that will return a list of Results instances with only one single dCt column...

        results = second_normaliser.get()

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

    def __run(self, reader, analyser, first_normaliser, normalisers, assays):
        """
        Core of _run
        """
        normalisers_list = []
        assays_list = []

        # analyse normalisers:
        for _normaliser in normalisers:
            norm = reader.read(_normaliser)

            for filter in self._Filters:
                norm = filter.pipe(norm)

            norm = analyser.pipe(norm)
            normalisers_list.append(norm)
        first_normaliser.link(normalisers = normalisers_list)

        # analyse sample assays
        for sample in assays:
            _sample = reader.read(sample)

            for filter in self._Filters:
                _sample = filter.pipe(_sample)

            _sample = analyser.pipe(_sample)
            assays_list.append(_sample)
        first_normaliser.link(samples = assays_list)

        first_normaliser.normalise()
        results = first_normaliser.get(copy=True)
        return results


class DoubleNorm2(BasicPlus):
    """
    This pipeline allows to normalise entire qpcr runs against each other. 
    It adds a second normaliser to the BasicPlus pipeline. Assays will be analysed by default DeltaDeltaCt,
    normalised run-internally against a normaliser (like Actin), 
    and eventually normalised against their corresponding assay from a control run. 

    NOTE: This pipeline works with an index table (csv) where filepaths alongside with an ID and corresponding experimental condition (such as CTRL, KD, whatever have to be specified...)
    """
    def __init__(self, index:str):
        super().__init__()
        self._index = self._read_index(index)
        self._Assays = {}
        self._Normalisers = {}

    def _read_index(self, file):
        """
        Reads the index table
        """
        index = pd.read_csv(file)
        all_good = all([i in index.columns for i in ["id", "condition", "normaliser", "path"]])
        if not all_good:
            wa.HardWarning("Pipeline:faulty_index")
        index["normaliser"] = index["normaliser"].astype(bool)
        
        return index

    def prune(self, results = True ):
        """
        Will clear results
        """
        if results: 
            self._df = None
            self._stats_df = None
            self._Results = None

    def _run(self):
        """
        Runs the pipeline
        """
        reader = qpcr.SampleReader()
        reader.replicates(self._replicates)
        if self._names is not None:
            reader.names(self._names)
        
        analyser = qpcr.Analyser()
        normaliser = qpcr.Normaliser()




if __name__ == "__main__":

    # norm_files = ["Example Data/28S.csv", "Example Data/actin.csv"]
    # sample_files = ["Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]

    # norm_folder = "Example Data 2/normalisers"
    # sample_folder = "Example Data 2/samples"

    # groupnames = ["wt-", "wt+", "ko-", "ko+"]
    
    # analysis = BasicPlus()
    # analysis.save_to("Example Data 2")
    # analysis.add_assays(sample_folder) # alternative: link() for iteratively linking new assays...
    # analysis.add_normalisers(norm_folder)
    
    # analysis.replicates(6)
    # analysis.names(groupnames)

    # iqr_filter = Filters.IQRFilter()
    # # iqr_filter.report("Example Data 2")
    # analysis.add_filters(iqr_filter)

    # preview = Plotters.PreviewResults(mode = "interactive")
    # analysis.add_plotters(preview)

    # # now that pipeline is ready, we can run!
    # analysis.run()

    # # now we can get results!
    # results = analysis.get(kind="df")
    # print(results)
    
    
    norm_folder = "Example Data 3/normalisers"
    exp_folder = "Example Data 3/experimental/samples"
    ctr_folder = "Example Data 3/experimental/samples"

    groupnames = ["wt-", "wt+", "ko-", "ko+"]
    
    analysis = DoubleNorm2(index = "Example Data 2/index.csv")
    analysis.save_to("Example Data 3")
    
    print(analysis._index)
    # analysis.set_experimental(assays = exp_folder, normalisers = norm_folder)
    # analysis.set_control(assays = ctr_folder, normalisers = norm_folder)

    # analysis.replicates(6)
    # analysis.names(groupnames)

    # iqr_filter = Filters.IQRFilter()
    # # iqr_filter.report("Example Data 2")
    # analysis.add_filters(iqr_filter)

    # preview = Plotters.PreviewResults(mode = "interactive")
    # analysis.add_plotters(preview)

    # # now that pipeline is ready, we can run!
    # analysis.run()

    # # now we can get results!
    # results = analysis.get(kind="df")
    # print(results)

    # exit(0)
