"""
This module contains a set of common lightweight wrappers to perform simple and non-customised DeltaDeltaCt analyses.
This is designed for lightweight users that do not wish to employ specialised pipelines.

> ### A Word on Plotters 
> Please, note that pipelines fully support mixing "static" and "interactive" Plotters, 
> but static figures will not stay open if interactive plotters are called to plot after them! 
> Because `qpcr.Filters` are always called to plot *before* any other `qpcr.Plotters`, this will mainly
> affect visualising the `qpcr.Plotters.ReplicateBoxplots` generated as Filter-Summaries.
"""

import qpcr
import matplotlib.pyplot as plt
import pandas as pd 
import statistics as stats
import qpcr._auxiliary.warnings as aw
import qpcr._auxiliary as aux
import qpcr.Plotters as Plotters
import qpcr.Filters as Filters
import qpcr.Readers as Readers
import re
import os 
import difflib

class Pipeline:
    """
    This is the basic template class for qpcr Pipelines. 
    It contains a set of basic preliminary methods
    that ensure that elementary required inputs are provided.

    Note
    ----
    The simplest implementation of this `Pipeline` template is the `Basic` pipeline.
    """
    def __init__(self):
        # super().__init__()
        self._Normalisers = []
        self._Assays = []
        self._save_to = None
        self._df = None
        self._stats_df = None
        self._Results = None
        self._replicates = None
        self._names = None
        self._softlink = True

    def assays(self):
        """
        Returns
        -------
        assays : list
            The linked `qpcr.Assay` objects for assays-of-interest.
        """
        assays = self._Assays
        return assays
    
    def normalisers(self):
        """
        Returns
        -------
        normalisers : list
            The linked `qpcr.Assay` objects for normaliser-assays.
        """
        normalisers = self._Normalisers
        return normalisers

    def replicates(self, replicates:(int or tuple)):
        """
        Set the replicates specifics to use for grouping.

        Parameters
        ----------
        replicates : int or tuple
            Can be an `integer` (equal group sizes, e.g. `3` for triplicates), 
            or a `tuple` (uneven group sizes, e.g. `(3,2,3)` if the second group is only a duplicate). 
        """
        self._replicates = replicates

    def names(self, names:(list or dict)):
        """
        Set names for replicates groups.

        Parameters
        ----------
        names : list or dict
            Either a `list` (new names without repetitions) or `dict` (key = old name, value = new name) specifying new group names. 
            Group names only need to be specified once, and are applied to all replicate entries.
        """
        self._names = names


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

        # vet if there are at least one normaliser and assay present
        if self._Normalisers == [] or self._Assays == []:
            aw.HardWarning("Pipeline:no_data")

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
            df = self._df
            if "assay" in df.columns: 
                df = df.drop(columns = ["assay"])
            return df
        elif kind == "obj":
            return self._Results
    
    def link(self, assays:(list or str) = None, normalisers:(list or str) = None):
        """
        Links new assays-of-interest / sample assays and/or normaliser assays 
        to the pipline, either replacing old ones or keeping them, 
        depending on `softlink()` settings.

        Parameters
        ----------
        assays : list or str
            A `list` of filepaths to raw datafiles of assays-of-interest, or a directory (`str`) where these are stored.
        normalisers : list or str
            A `list` of filepaths to raw datafiles of normaliser assays, or a directory (`str`) where these are stored.
        """
        if self._softlink:
            self.prune()
        self.add_assays(assays)
        self.add_normalisers(normalisers)


    def prune(self, assays = True, results = True, normalisers = True):
        """
        Will clear assays, results, and/or normalisers

        Parameters
        ----------
        assays : bool
            Will clear any sample assays in the pipeline if True (default).
        
        results : bool
            Will clear any computed results in the pipline if True (default).
        
        normalisers : bool
            Will clear any normalisers in the pipline if True (default).
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
    
    def add_assays(self, assays):
        """
        Adds assays-of-interest / sample assays (filepaths) (keeping any already present)

        Parameters
        ----------
        assays : list or str
            A `list` of filepaths to raw datafiles of sample assays, or a directory (`str`) where these are stored.
        """
        assays = self._from_directory(assays)
        self._Assays.extend(assays)

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

            # if inputs are a directory
            if os.path.isdir(files):
            
                # get files from the directory
                datafiles = os.listdir(files)
            
                # if no files are found, raise error
                if len(datafiles) == 0:
                    aw.HardWarning("Pipeline:no_data_input", file = files, traceback = False)
            
                # combine paths with parent directory
                datafiles = [os.path.join(files,n) for n in datafiles]
                return datafiles

            # else check if inputs are a single file
            # then just put it into a list to be valid input for .extend()
            elif os.path.isfile(files):
                
                return [files]

            else: 
                aw.HardWarning("Pipeline:no_data_input", file = files, traceback = False)
        
        # else check if we got a list or tuple of files
        elif isinstance(files, (list or tuple)):
            return files
        
        # raise error for anything else...
        else: 
            aw.HardWarning("Pipeline:no_data_input", file = files, traceback = False)


    def _run(self, **kwargs):
        """
        This is the actual function that each custom Pipeline has to define...
        """
        print("self._run() is the actual function that each custom Pipeline has to define...")


class Basic(Pipeline):
    """
    Performs simple standardized DeltaDeltaCt analysis 
    based on two lists of files, one for normaliser assays and one for sample assays.
    This makes use of standard settings for `qpcr.Analyser` and `qpcr.Normaliser`
    which cannot be customized! 
    For customization check out the `Blueprint` pipeline or generate your own.
    """
    def __init__(self):
        super().__init__()
    
    def _run(self, **kwargs):
        """
        The automated standard DeltaDeltaCt pipeline. 
        Using default settings of `qpcr.Analyser` and `qpcr.Normaliser`
        """
        reader = qpcr.DataReader()
        analyser = qpcr.Analyser()
        normaliser = qpcr.Normaliser()

        # add replicate and names info to kwargs
        kwargs = dict(kwargs, replicates = self._replicates, names = self._names)

        # analyse normalisers
        normalisers = [  reader.read(i, **kwargs) for i in self._Normalisers  ]     
        normalisers = [  analyser.pipe(i) for i in normalisers  ]
       
        # analyse sample assays
        assays = [  reader.read(i, **kwargs) for i in self._Assays  ]   
        assays = [  analyser.pipe(i) for i in assays  ]

        # now normalise
        normaliser.link( assays = assays, normalisers = normalisers )
        normaliser.normalise()

        # get and store results
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
    
    def prune(self, assays = True, results = True, normalisers = True, figures = True, plotters = False, filters = False):
        """
        Will clear the pipeline.

        Parameters
        ----------
        assays : bool
            Will clear any sample assays in the pipeline if True (default).
        
        results : bool
            Will clear any computed results in the pipline if True (default).
        
        normalisers : bool
            Will clear any normalisers in the pipline if True (default).
        figures : bool
            Will clear any figures in the pipline if True (default).
        plotters : bool
            Will clear any Plotters from the pipeline if True (default False).
        filters : bool
            Will clear any Filters in the pipline if True (default False).
        """
        super().prune(assays = assays, results = results, normalisers = normalisers)
        if figures: 
            self._Figures = []
        if filters:
            self._Filters = []
        if plotters: 
            self._Plotters = []

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

    def _run(self, **kwargs):
        """
        The automated standard DeltaDeltaCt pipeline. 
        Also produces applies and Plotters to produce figures...
        """

        # simply set up the Blueprint pipeline 
        # and run with default settings
        pipeline = Blueprint()

        pipeline.replicates(self._replicates)
        pipeline.names(self._names)

        if self._save_to is not None:
            pipeline.save_to(self._save_to)

        pipeline.add_assays(self._Assays)
        pipeline.add_normalisers(self._Normalisers)

        pipeline.add_plotters(*self._Plotters)
        pipeline.add_filters(*self._Filters)

        pipeline.run()

        # and store the results
        self._Results = pipeline.get(kind = "obj")
        self._df = pipeline.get(kind = "df")
        self._stats_df = pipeline.get(kind = "stats")




class Blueprint(BasicPlus):
    """
    Performs simple Delta-Delta-Ct analysis based on the same workflow as the `BasicPlus` pipeline, but allows full costumization of SampleReader, Analyser, and Normaliser objects.
    Optionally, `qpcr.SampleReader`, `qpcr.Analyser`, and `qpcr.Normaliser` may be set up externally and linked into the pipeline. Any non-linked processing classes will be set up using defaults.
    """
    def __init__(self):
        super().__init__()
        self._Reader = None
        self._Analyser = None
        self._Normaliser = None

    def prune(self, assays = True, results = True, normalisers = True, figures = True, plotters = False, filters = False, cores = False, reader = False, analyser = False, normaliser = False):
        """
        Will clear the pipeline.

        Parameters
        ----------
        assays : bool
            Will clear any sample assays in the pipeline if True (default).
        
        results : bool
            Will clear any computed results in the pipline if True (default).
        
        normalisers : bool
            Will clear any normalisers in the pipline if True (default).
        figures : bool
            Will clear any figures in the pipline if True (default).
        plotters : bool
            Will clear any Plotters from the pipeline if True (default False).
        filters : bool
            Will clear any Filters in the pipline if True (default False).
        cores : bool
            Will clear Reader, Analyser, and Normaliser if True (default False).
        reader : bool
            Will only clear the Reader if True (default False).
        analyser : bool
            Will only clear the Analyser if True (default False).
        normaliser : bool
            Will only clear the Normaliser if True (default False).
            Note, that clearing `results` will also clear the `Normaliser`'s results,
            but keep the Normaliser itself!
        """
        super().prune(
                        assays = assays, 
                        results = results, 
                        normalisers = normalisers,
                        figures = figures, 
                        plotters = plotters, 
                        filters = filters
                    )
        
        self._Normaliser.prune(
                                assays = assays, 
                                normalisers = normalisers, 
                                results = results
                            )
        if cores: 
            self._Reader = None
            self._Analyser = None
            self._Normaliser = None
        else:
            if reader:
                self._Reader = None
            if analyser:
                self._Analyser = None
            if normaliser:
                self._Normaliser = None

    def Reader(self, Reader : qpcr.SampleReader = None):
        """
        Links a `qpcr.SampleReader` object to the pipeline.

        Parameters
        ----------
        Reader : qpcr.SampleReader
            A `qpcr.SampleReader` object
        """
        if Reader is not None: 
            self._Reader = Reader
        return self._Reader

    def Analyser(self, Analyser : qpcr.Analyser = None):
        """
        Links a `qpcr.Analyser` object to the pipeline.

        Parameters
        ----------
        Analyser : qpcr.Analyser
            A `qpcr.Analyser` object
        """
        if Analyser is not None:
            self._Analyser = Analyser
        return self._Analyser

    def Normaliser(self, Normaliser : qpcr.Normaliser = None):
        """
        Links a `qpcr.Normaliser` object to the pipeline.

        Parameters
        ----------
        Normaliser : qpcr.Normaliser
            A `qpcr.Normaliser` object
        """
        if Normaliser is not None:
            self._Normaliser = Normaliser
        return self._Normaliser
    
    def _setup_cores(self):
        """
        Sets SampleReader, Analyser, and Normaliser to defaults, if no external ones were provided...
        """
        if self.Reader() is None: 
            self.Reader(qpcr.DataReader())
        if self.Analyser() is None: 
            self.Analyser(qpcr.Analyser())
        if self.Normaliser() is None:
            self.Normaliser(qpcr.Normaliser())

    def _run(self, **kwargs):
        """
        The automated standard DeltaDeltaCt pipeline. 
        Also produces applies and Plotters to produce figures...
        """
        # setup Analyser, and Normaliser (if none were provided)
        self._setup_cores()
        
        reader = self.Reader()
        analyser = self.Analyser()
        normaliser = self.Normaliser()

        # check if we have plotters or filters
        have_plotters = len(self._Plotters) != 0 
        have_filters = len(self._Filters) != 0

        # add replicate and names info to kwargs
        kwargs = dict(kwargs, replicates = self._replicates, names = self._names)

        # analyse normalisers
        normalisers = [  reader.read(i, **kwargs) for i in self._Normalisers  ]

        if have_filters:
            for filter in self._Filters:
                normalisers = [  filter.pipe(i) for i in normalisers  ]
            
        normalisers = [  analyser.pipe(i) for i in normalisers  ]
       
        # analyse sample assays
        assays = [  reader.read(i, **kwargs) for i in self._Assays  ]

        if have_filters:
            for filter in self._Filters:
                assays = [  filter.pipe(i) for i in assays  ]
        
        assays = [  analyser.pipe(i) for i in assays  ]

        # now normalise
        normaliser.link( assays = assays, normalisers = normalisers )
        normaliser.normalise()

        # get and store results
        results = normaliser.get()
        self._Results = results
        self._df = results.get()
        self._stats_df = results.stats()

        if self._save_to is not None:
            results.save(self._save_to)

        # plot filtering report
        if have_filters:
            for filter in self._Filters:
                if self._save_to is not None or filter.report() is not None:
                    
                    # add report location if none was specified...
                    if filter.report() is None: 
                        filter.report(self._save_to)

                # save filter summary fig
                figs = filter.plot()
                self._Figures.extend(figs)

        # plot results
        if have_plotters:
            for plotter in self._Plotters:

                # generate plotter fig and save
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


class _Qupid_Blueprint(Blueprint):
    """
    This is the implementation of the Plueprint pipeline for the Qupid webapp.
    It redefines the file-reading functions with methods compatible with the UploadedFile objects from streamlit.
    """
    def __init__(self):
        super().__init__()

    def _setup_cores(self):
        """
        Sets SampleReader, Analyser, and Normaliser to defaults, if no external ones were provided...
        """
        if self.Reader() is None: 
            self.Reader(qpcr._Qupid_SampleReader())
        if self.Analyser() is None: 
            self.Analyser(qpcr.Analyser())
        if self.Normaliser() is None:
            self.Normaliser(qpcr.Normaliser())


# We push this one behind in the release...
# Since we now have the ddCt pipeline we should be able to slim down most of the pipelines here anyway...

# class MultiAssay(Blueprint):
#     """
#     Performs Delta-Delta-Ct based on data from a single multi-assay datafile.
#     Datasets within this datafile must be decorated to identify them as assays-of-interest or normalisers.
#     Check out the documentation of `qpcr.Parsers` for more information on decorators.

#     Note
#     -------
#     This class relies on `qpcr.Parsers` and decorated assays to get its input data. 
#     If your file does not offer this kind of architecture, choose another pipeline.

#     """
#     def __init__(self):
#         super().__init__()
#         self._src = None
    
#     def link(self, filename : str, **kwargs):
#         """
#         Reads a datafile in csv or excel format containing 
#         multiple decorated datasets and extracts sample and normalisers assays. 

#         Parameters
#         ----------
#         filename : str
#             A filepath to a raw data file.
#         **kwargs
#             Any additional keyword arguments that shall be passed to the `qpcr.MultiReader`'s `pipe` method.
#         """
#         self._src = filename

#         # setup Reader, Analyser, and Normaliser (if none were provided)
#         self._setup_cores(**kwargs)
#         reader = self.Reader()

#         # read the multi-assay datafile
#         self._Assays, self._Normalisers = reader.pipe(self._src, **kwargs)
        
        
#     def _run(self, **kwargs):
#         """
#         The main workflow
#         """
#         # setup Reader, Analyser, and Normaliser (if none were provided)
#         self._setup_cores()
        
#         analyser = self.Analyser()
#         normaliser = self.Normaliser()

#         normalisers = []
#         samples = []

#         # analyse normalisers:
#         for norm in self._Normalisers:
#             for filter in self._Filters:
#                 norm = filter.pipe(norm)
#             norm = analyser.pipe(norm)
#             normalisers.append(norm)
#         normaliser.link(normalisers = normalisers)

#         # analyse sample assays
#         for sample in self._Assays:
#             for filter in self._Filters:
#                 sample = filter.pipe(sample)
#             sample = analyser.pipe(sample)
#             samples.append(sample)
#         normaliser.link(assays = samples)

#         normaliser.normalise()
#         results = normaliser.get()

#         self._Results = results
#         self._df = results.get()
#         self._stats_df = results.stats()

#         if self._save_to is not None:
#             results.save(self._save_to)

#         # plot filtering report
#         for filter in self._Filters:
#             if self._save_to is not None or filter.report() is not None:
#                 # add report location if none was specified...
#                 if filter.report() is None: 
#                     filter.report(self._save_to)

#             figs = filter.plot()
#             self._Figures.extend(figs)

#         # plot results
#         for plotter in self._Plotters:
#             plotter.link(self._Results)
#             fig = plotter.plot()
#             self._Figures.append(fig)

#             if self._save_to is not None:
#                 filename = self._make_figure_filename(plotter)
#                 plotter.save(filename)

#     def _setup_cores(self, **kwargs):
#         """
#         Sets Reader, Analyser, and Normaliser to defaults, if no external ones were provided...
#         """
#         # check if a sheet_name was specified in case of multi-sheet files...
#         use_multi_sheet = "sheet_name" not in kwargs
#         if self.Reader() is None: 
#             if self._is_multisheet() and use_multi_sheet:
#                 self.Reader(Readers.MultiSheetReader())
#             else: 
#                 self.Reader(Readers.MultiReader())
#         if self.Analyser() is None: 
#             self.Analyser(qpcr.Analyser())
#         if self.Normaliser() is None:
#             self.Normaliser(qpcr.Normaliser())
    
#     def _is_multisheet(self):
#         """
#         Checks if a provided excel file contains multiple sheets
#         """
#         verdict = False
#         if self._src.endswith("xlsx"):
#             data = pd.read_excel(self._src, sheet_name = None)
#             verdict = len(data.keys()) > 1
#         return verdict

#     # add_assays is disabled
#     def add_assays(self):
#         print("To provide a data input file use link()!\nIf you wish to supply separate files for assays-of-interest and normalisers, checkout another Pipeline as this one only works with a single Mlit-Assay file!")

#     # add_normalisers is disabled
#     def add_assays(self):
#         print("To provide a data input file use link()!\nIf you wish to supply separate files for assays-of-interest and normalisers, checkout another Pipeline as this one only works with a single Mlit-Assay file!")
    
class ddCt(Blueprint):
    """
    Performs only Delta-Delta-Ct and requires `qpcr.Assay` objects as inputs.
    Hence, this pipeline does NOT read any files!

    It follows the default workflow of the `BasicPlus`pipeline and is based on the `BluePrint`
    pipeline to allow customisation.

    Note
    ----
    As the pipeline inherits from the `Blueprint` pipeline it does have a `Reader` method (which won't do anything though!). 
    """
    def __init__(self):
        super().__init__()
        
    def _setup_cores(self):
        """
        Sets Analyser, and Normaliser to defaults, if no external ones were provided...
        """
        if self.Analyser() is None: 
            self.Analyser(qpcr.Analyser())
        if self.Normaliser() is None:
            self.Normaliser(qpcr.Normaliser())
        
    def _run(self, **kwargs):
        # setup Analyser, and Normaliser (if none were provided)
        self._setup_cores()
        
        analyser = self.Analyser()
        normaliser = self.Normaliser()

        # check if we have plotters or filters
        have_plotters = len(self._Plotters) != 0 
        have_filters = len(self._Filters) != 0

        # analyse normalisers
        normalisers = [  i for i in self._Normalisers  ]

        if have_filters:
            for filter in self._Filters:
                normalisers = [  filter.pipe(i) for i in normalisers  ]
            
        normalisers = [  analyser.pipe(i) for i in normalisers  ]
       
        # analyse sample assays
        assays = [  i for i in self._Assays  ]

        if have_filters:
            for filter in self._Filters:
                assays = [  filter.pipe(i) for i in assays  ]
        
        assays = [  analyser.pipe(i) for i in assays  ]

        # now normalise
        normaliser.link( assays = assays, normalisers = normalisers )
        normaliser.normalise()

        # get and store results
        results = normaliser.get()
        self._Results = results
        self._df = results.get()
        self._stats_df = results.stats()

        if self._save_to is not None:
            results.save(self._save_to)

        # plot filtering report
        if have_filters:
            for filter in self._Filters:
                if self._save_to is not None or filter.report() is not None:
                    
                    # add report location if none was specified...
                    if filter.report() is None: 
                        filter.report(self._save_to)

                # save filter summary fig
                figs = filter.plot()
                self._Figures.extend(figs)

        # plot results
        if have_plotters:
            for plotter in self._Plotters:

                # generate plotter fig and save
                plotter.link(self._Results)
                fig = plotter.plot()
                self._Figures.append(fig)

                if self._save_to is not None:
                    filename = self._make_figure_filename(plotter)
                    plotter.save(filename)

if __name__ == "__main__":

    norm_files = ["./Examples/Example Data/28S.csv", "./Examples/Example Data/actin.csv"]
    sample_files = ["./Examples/Example Data/HNRNPL_nmd.csv", "./Examples/Example Data/HNRNPL_prot.csv"]

    norm_folder = "./Example Data 3/normalisers/"
    sample_folder = "./Example Data 3/samples/"

    groupnames = ["wt-", "wt+", "ko-", "ko+"]
    
    analysis = Blueprint()
    
    # grouped_analyser = qpcr.Analyser()
    # grouped_analyser.anchor("first")
    # analysis.Analyser(grouped_analyser)

    # analysis.save_to("Example Data 2")
    # analysis.add_assays(sample_files) # alternative: link() for iteratively linking new assays...
    # analysis.add_normalisers(norm_files)

    analysis.add_assays(sample_folder)
    analysis.add_normalisers(norm_folder)
    
    # print("No Reps specified, all inferred!")
    analysis.replicates(6)
    analysis.names(groupnames)

    iqr_filter = Filters.RangeFilter()
    iqr_filter.plotmode("interactive")

    iqr_filter.report("Example Data 4")
    analysis.add_filters(iqr_filter)

    preview = Plotters.PreviewResults(mode = "static")
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
