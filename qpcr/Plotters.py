"""
This module is designed for streamlined data visualisation of the qpcr generated results
It is designed to work directly with qpcr.Results() instances. 
"""

import __init__ as qpcr
import Pipes
import auxiliary.graphical as gx
import auxiliary as aux 
import auxiliary.warnings as wa
import pandas as pd 
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.express as px
import numpy as np 


# setup some default settings for charts

_default_static_PreviewResults = dict(
                                        color = "Lightgray",
                                        rot = 0,
                                        legend = False, 
                                        title = "Preview of Results"
                                    )

_default_interactive_PreviewResults = dict(
                                            template = "plotly",
                                            title = "Preview of Results"
                                        )


# Concent: 
# Plotter: 
#       We define a superclass Plotter that will handle linking to data, 
#       linking to default parameters, and setting up default parameters.
#       It also provides a generic plot() method that will call on a FigureClass specific _plot() method...
#       Wether or not to use static or interactive plots is also handled by this class...
#
# FigureClass:
#       A parent class for each type of figure. It contains two potential _plot() 
#       methods, one for interactive one for static plotting... Which one to use is decided based on the plotting mode...
#       Each FigureClass has therefore to link default parameters, then init Plotter (superclass), 
#       and define their own _static_plot() and _interactive_plot() methods! Any additionally required methods can be written as well...


class Plotter:
    """
    A superclass that handles Data Linking and Parameter setup for FigureClasses
    """
    def __init__(self, mode = "interactive"):
        self._default_params = None
        self._PARAMS = {}
        self._Results = None
        self._data = None
        self._id = type(self).__name__
        self._MODE = mode
        self._set_plot()

    def link(self, Results:(qpcr.Results or pd.DataFrame)):
        """
        Link a Results object or pandas DataFrame of the same architecture
        as one handled by a Results object. 
        Note, that this will replace any previously linked data!
        """
        if isinstance(Results, qpcr.Results):
            self._Results = Results
            self._data = self._Results.stats()
        elif isinstance(Results, pd.DataFrame):
            self._Results = None
            self._data = Results
        else:
            wa.HardWarning("Plotter:unknown_data")

    def plot(self, **kwargs):
        """
        Generate Figure
        """
        total_kwargs = self.update_params(kwargs)
        self._plot(**total_kwargs)

    def id(self, id:str = None):
        """
        Set a unique Id, by default the Classname will be used.
        """
        if id is not None:
            self._id = id
        else:
            return self._id

    def get(self):
        """
        Returns the DataFrame used for plotting
        """
        return self._data

    def params(self, **params):
        """
        Set default parameters for plotting (will be forwarded to **kwargs)
        Returns default parameters if no new parameters are added.
        """
        if params != {} and self._PARAMS == {}:
            self._PARAMS = params
        elif params != {}:
            self.update_params(params, store = True)
        return self._PARAMS

    def update_params(self, kwargs, supersede = True, store = False):
        """
        Appends pre-set parameters to kwargs. 
        In case of key duplications: 
        It will either replace old ones with new ones (supersede = True, default), 
        or keep old ones (supersede = False)
        """
        if supersede:
            kwargs = dict(self._PARAMS, **kwargs)
        else:
            kwargs = kwargs = dict(kwargs, **self._PARAMS)
        
        if store:
            self._PARAMS = kwargs

        return kwargs

    def _setup_default_params(self, static:dict, interactive:dict):
        """
        Setup and set default parameters for a FigureClass
        """
        self._static_default = static
        self._interactive_default = interactive

    def _set_plot(self):
        """
        Sets self._plot either interactive or static depending on MODe
        """
        self._plot  = self._static_plot if self._MODE == "static" else self._interactive_plot
        prev_default = self._default_params
        self._default_params = self._static_default if self._MODE == "static" else self._interactive_default

        if self.params() == {} or self.params() == prev_default:
            self.params(**self._default_params)
        else:
            self.update_params(self._default_params, store = True, supersede = False)

    def _static_plot(self, **kwargs):
        """
        The plot function that will handle static plotting (will be redefined for each FigureClass)
        """
        print("The plot function that will handle static plotting...")
        print("Surprised to see this?, \nPerhaps your desired plotting methods are not named properly. Make sure to name your method _static_plot()!")

    def _interactive_plot(self, **kwargs):
        """
        The plot function that will handle interactive plotting (will be redefined for each FigureClass)
        """
        print("The plot function that will handle interactive plotting...")
        print("Surprised to see this?, \nPerhaps your desired plotting methods are not named properly. Make sure to name your method _interactive_plot()!")

    def _prep_properties(self):
        """
        Setup ncols, nrows, subplot titles (headers), x, y, and sterr, columns for figure...
        """
        self._setup_default_plot_cols()
        kwargs = self.params()
        data = self.get()

        # setup reference column and figure subplots
        ref_col = aux.from_kwargs("key", "assay", kwargs, rm=True)
        ncols, nrows = aux.from_kwargs("subplots", 
                                        gx.make_layout(data, ref_col), 
                                        kwargs, rm=True
                                       )

        # transposing currently not supported...
        # if aux.from_kwargs("transpose", False, kwargs, rm = True):
        #     ncols, nrows = nrows, ncols

        headers = aux.from_kwargs("headers", aux.sorted_set(data[ref_col]), kwargs, rm = True)

        x = aux.from_kwargs("x", self._default_x, kwargs, rm = True)
        y = aux.from_kwargs("y", self._default_y, kwargs, rm = True)
        sterr = aux.from_kwargs("sterr", self._default_sterr, kwargs, rm = True)
        
        # the query to be used if group or group_name is ref_col
        query = "{ref_col} == '{q}'" if isinstance(headers[0], str) else "{ref_col} == {q}"

        return ref_col,ncols,nrows, headers, x, y, sterr, query

    def _setup_default_plot_cols(self):
        """
        Sets default columns for x and y in barcharts,
        default x
        - "group_name" column if present 
        - "group" Otherwise
        default y
        - "mean"
        default sterr
        - "stdev"
        """
        columns = self._data.columns
        self._default_x = "group_name" if "group_name" in columns else "group"
        self._default_y = "mean"
        self._default_sterr = "stdev"

class PreviewResults(Plotter):
    """
    Generate a Preview of all results from all Assays in subplots.
    """
    def __init__(self, mode:str):
        self._setup_default_params(
                                    static = _default_static_PreviewResults, 
                                    interactive = _default_interactive_PreviewResults
                                )
        # __init__ is going to require default_params to be already set!
        super().__init__(mode)

    def _static_plot(self, **kwargs):
        """
        The static Preview Results Figure
        """
        kwargs = self.update_params(kwargs)
        data = self.get()

        ref_col, ncols, nrows, headings, x, y, sterr, query = self._prep_properties()

        headers = aux.from_kwargs("headers", None, kwargs, rm = True)
        label_subplots = aux.from_kwargs("label_subplots", True, kwargs, rm = True)
        start_character = aux.from_kwargs("labeltype", "A", kwargs, rm = True)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        show = aux.from_kwargs("show", True, kwargs, rm = True)

        figsize = aux.from_kwargs("figsize", None, kwargs, rm = True)
        fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize)

        # setup Coords
        Coords = gx.AxesCoords(fig, axs, (nrows, ncols))

        idx = 0
        for assay in headings: 
            try: 
                tmp_df = data.query(query.format(ref_col = ref_col, q = assay))
                
                # now plot a new bar chart 
                subplot = Coords.subplot()

                tmp_df.plot.bar(
                                x = x, y = y, 
                                ax = subplot,

                                edgecolor = aux.from_kwargs("edgecolor", "black", kwargs),
                                linewidth = aux.from_kwargs("borderwidth", 0.3, kwargs, rm = True),
                                **kwargs
                            )

                subplot.errorbar(
                                x = tmp_df[x], y = tmp_df[y], 
                                yerr = tmp_df[sterr], 
                                fmt = ".", markersize = 0, capsize = 3, 
                                ecolor = aux.from_kwargs("edgecolor", "black", kwargs),
                            )

                subplot.set(
                            title = assay if headers is None else headers[idx],
                            xlabel = aux.from_kwargs("xlabel", None, kwargs),
                            ylabel = aux.from_kwargs("ylabel", "$\Delta\Delta$Ct", kwargs),        
                        )

                if not show_spines:
                    subplot.spines["right"].set_visible(False)
                    subplot.spines["top"].set_visible(False)
                    subplot.spines["left"].set_linewidth(1.05)
                    subplot.spines["bottom"].set_linewidth(1.05)

                # add ABCD... label to subplot
                if label_subplots:
                    self._add_subplot_label(idx, subplot, start_character)         

                Coords.increment()
                idx += 1
            except Exception as e:
                print(e) 
                break
        
        plt.tight_layout()

        
        if show:
            plt.show()

        # print("<<< Static Check >>>")

        return fig 

    def _add_subplot_label(self, idx, subplot, start_character):
        """
        Adds A B C ... to upper left corner of a subplot...
        """
        subplot_label = chr(ord(start_character)+idx)
        subplot.annotate(
                    xy = (-0.1,1.03), 
                    text = subplot_label, 
                    xycoords = "axes fraction",
                    weight = "bold", fontsize = 12
                )

    def _interactive_plot(self, **kwargs):
        """
        The interactive Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self.get()

        try: 
            # setup figure framework variables that have to be removed from 
            ref_col, ncols, nrows, headings, x, y, sterr, query = self._prep_properties()
            headers = aux.from_kwargs("headers", headings, kwargs, rm = True)
            show = aux.from_kwargs("show", True, kwargs, rm = True)

            speclist = gx.make_speclist(nrows, ncols, "xy")
            fig = make_subplots(rows = nrows, cols = ncols, specs = speclist, subplot_titles = headers)


            Coords = gx.AxesCoords(fig, [], (ncols, nrows))

            fig.update_layout(
                        title = aux.from_kwargs("title", "Results Preview", kwargs, rm = True), 
                        height = aux.from_kwargs("height", None, kwargs, rm = True), 
                        width = aux.from_kwargs("width", None, kwargs, rm = True),
                        margin = {"autoexpand": True, "pad" : 0, "b": 1, "t": 50}, 
                        autosize = True, 
                        template = aux.from_kwargs("template", "plotly_white", kwargs, rm = True),
                        legend = {"title" : aux.from_kwargs("legend_title", ref_col, kwargs, rm = True),
                        },
                    )

            idx = 0
            for assay in headings:
                row, col = Coords.get()
                tmp_df = data.query(query.format(ref_col = ref_col, q = assay))
                # now plot a new bar chart 
                fig.add_trace(

                    go.Bar(
                        name = assay,
                        y = tmp_df[y], x = tmp_df[x], 
                        error_y=dict(type='data', array = tmp_df[sterr]),
                        hoverinfo = aux.from_kwargs("hoverinfo", "y", kwargs, rm = True), 
                        **kwargs
                    ), 

                    row, col
                )
                Coords.increment()

                idx += 1

            if show:
                fig.show()

            # print("<<< Interactive Check >>>")
            return fig 
        except Exception as e: 
            raise e 

if __name__ == '__main__':

    files = ["Example Data/28S.csv", "Example Data/28S_again.csv", "Example Data/actin.csv",
    "Example Data/actin2.csv",
    "Example Data/actin3.csv",
    "Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]
    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    # setup figures
    a = PreviewResults(mode = "static")
    b = PreviewResults(mode = "interactive")

    a.params(
          frame = True, labeltype = "a", show = True
    )
    b.params(
        show = True, template = "plotly"  
    )

    # predefined pipeline use
    pipe = Pipes.Basic()
    pipe.link(files)
    pipe.add_normalisers(files[:2])
    pipe.replicates(6)
    pipe.names(groupnames)


    print()
    print("==== Multiple (ALL) Samples ====")
    pipe.run()
    result = pipe.get(kind = "obj")
    try:
        a.link(result)
        a.plot(color = "green")
    except: 
        print("Static Failed!")

    try:
        b.link(result)
        b.plot()
    except: print("Interactive Failed!")

    print()
    print("==== Multiple (7) Samples ====")
    pipe.link(files[:7])
    pipe.run()
    result = pipe.get(kind = "obj")
    try:
        a.link(result)
        a.plot(show = True)
    except: 
        print("Static Failed!")

    try:
        b.link(result)
        b.plot(template = "plotly")
    except: print("Interactive Failed!")

    print()
    print("==== Multiple (6) Samples ====")
    pipe.link(files[:6])
    pipe.run()
    result = pipe.get(kind = "obj")
    try:
        a.link(result)
        a.plot()
    except: 
        print("Static Failed!")

    try:
        b.link(result)
        b.plot()
    except: print("Interactive Failed!")

    print()
    print("==== Multiple (4) Samples ====")
    pipe.link(files[:4])
    pipe.run()
    result = pipe.get(kind = "obj")
    try:
        a.link(result)
        a.plot()
    except: 
        print("Static Failed!")

    try:
        b.link(result)
        b.plot()
    except: print("Interactive Failed!")

    print()
    print("==== Multiple (3) Samples ====")
    pipe.link(files[:3])
    pipe.run()
    result = pipe.get(kind = "obj")
    try:
        a.link(result)
        a.plot()
    except Exception as e: 
        print(e)
        print("Static Failed!")

    try:
        b.link(result)
        b.plot()
    except Exception as e: 
        print("Interactive Failed!")
        raise e

    print()
    print("==== One (1) Sample ====")
    pipe.link([files[0]])
    pipe.run()
    result = pipe.get(kind = "obj")
    try:
        a.link(result)
        a.plot()
    except: 
        print("Static Failed!")

    try:
        b.link(result)
        b.plot()
    except: print("Interactive Failed!")
    exit(0)