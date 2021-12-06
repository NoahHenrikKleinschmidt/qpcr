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

class Chart:
    """
    Setup charts to be either interactive (produced with plotly) or static (matplotlib).
    mode = "interactive" or "static" can be set.
    """
    def __init__(self, mode = "interactive"):
        self._MODE = mode if mode is not None else "interactive"
        self._data = None
        self._PARAMS = {}
        self._id = type(self).__name__

    def id(self, id:str = None):
        """
        Set a unique Id, by default simply the Classname will be used...
        """
        if id is not None:
            self._id = id
        else:
            return self._id

    def get(self):
        """
        Returns the stats results dataframe
        """
        return self._data
    
    def link(self, data:(qpcr.Results or pd.DataFrame) = None):
        """
        Link a qpcr.Results instance or a pd.DataFrame 
        Note a pd.DataFrame has to have the same architecture 
        as one obtained from qpcr.Results.stats()!
        Can also return the linked dataframe
        """
        if isinstance(data, qpcr.Results):
            self._data = data.stats()
        elif isinstance(data, pd.DataFrame):
            self._data = data
        else:
            wa.HardWarning("Chart:unknown_data", obj = data)
        
        self._setup_default_plot_cols()

    def params(self, **params):
        """
        Set default parameters for plotting (will be forwarded to **kwargs)
        """
        self._PARAMS = params

    def update_params(self, kwargs):
        """
        Appends default parameters to kwargs
        """
        kwargs = dict(kwargs, **self._PARAMS)
        return kwargs

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

    

class PreviewResults(Chart):
    """
    Generates a bar chart for each Assay separately, and produces a subplots figure.
    """
    def __init__(self, mode:str = None):
        super().__init__(mode)
        if self._MODE == "interactive":
            self._CORE = _plotly_PreviewResults
        else:
            self._CORE = _matplotlib_PreviewResults
    

    def plot(self, **kwargs):
        """
        Plots the data into figure
        """
        kwargs = self.update_params(kwargs)
        plotter = self._CORE(self)
        fig = plotter._plot(**kwargs)
        return fig 


    def _prep_properties(self, kwargs):
        """
        Setup ncols, nrows, heading titles, x, y, and sterr, columns for figure...
        """
        # setup reference column and figure subplots
        ref_col = aux.from_kwargs("key", "assay", kwargs, rm=True)
        ncols, nrows = aux.from_kwargs("subplots", 
                                        gx.make_layout(self._data, ref_col), 
                                        kwargs, rm=True
                                       )

        if aux.from_kwargs("transpose", False, kwargs):
            ncols, nrows = nrows, ncols

        headings = aux.sorted_set(self._data[ref_col])
        x = aux.from_kwargs("x", self._MASTER._default_x, kwargs, rm = True)
        y = aux.from_kwargs("y", self._MASTER._default_y, kwargs, rm = True)
        sterr = aux.from_kwargs("sterr", self._MASTER._default_sterr, kwargs, rm = True)
        
        # the query to be used if group or group_name is ref_col
        query = "{ref_col} == '{q}'" if isinstance(headings[0], str) else "{ref_col} == {q}"

        return ref_col,ncols,nrows,headings, x, y, sterr, query
    


class _plotly_PreviewResults(PreviewResults):
    """
    Interactive PreviewChart
    """
    def __init__(self, master):
        self._MASTER = master
        self._data = self._MASTER.get()
        

    def _plot(self, **kwargs):
        """
        Generates a PreviewResults Figure
        """
        try: 
            # setup figure framework variables that have to be removed from 
            # kwargs before passing them to df.plot()
            ref_col, ncols, nrows, headings, x, y, sterr, query = self._prep_properties(kwargs)

            headers = aux.from_kwargs("headers", headings, kwargs, rm = True)

            

            speclist = gx.make_speclist(nrows, ncols, "xy")
            fig = make_subplots(rows = nrows, cols = ncols, specs = speclist, subplot_titles = headers)

            # fig.print_grid()

            Coords = gx.AxesCoords(fig, [], (ncols, nrows))
            # Coords.transpose(not aux.from_kwargs("transpose", False, kwargs, rm = True))


            idx = 0
            for assay in headings:
                
                row, col = Coords.get()
                
                tmp_df = self._data.query(query.format(ref_col = ref_col, q = assay))
                # print("cooreds: ", row, col)
                # now plot a new bar chart 
                fig.add_trace(

                    go.Bar(
                        name = assay,
                        y = tmp_df[y], x = tmp_df[x], 
                        error_y=dict(type='data', array = tmp_df[sterr]),
                        hoverinfo = aux.from_kwargs("hoverinfo", "y", kwargs), 
                    ), 
                    row, col
                )
                Coords.increment()

                idx += 1

            fig.update_layout(
                        title = aux.from_kwargs("title", "Results Preview", kwargs), 
                        height = aux.from_kwargs("height", None, kwargs), 
                        width = aux.from_kwargs("width", None, kwargs),
                        margin = {"autoexpand": True, "pad" : 0, "b": 1, "t": 50}, 
                        autosize = True, 
                        template = aux.from_kwargs("template", "plotly_white", kwargs),
                        legend = {"title" : aux.from_kwargs("legend_title", ref_col, kwargs),
                        },
                    )

            if aux.from_kwargs("show", True, kwargs):
                fig.show()

            print("<<< Interactive Check >>>")
            return fig 
        except Exception as e: 
            raise e 


class _matplotlib_PreviewResults(PreviewResults):
    """
    Static PreviewChart
    """
    def __init__(self, master:Chart):
        self._MASTER = master
        self._data = self._MASTER.get()

    def _plot(self, **kwargs):
        """
        Generates a PreviewResults Figure
        """
        # setup figure framework variables that have to be removed from 
        # kwargs before passing them to df.plot()
        ref_col, ncols, nrows, headings, x, y, sterr, query = self._prep_properties(kwargs)

        headers = aux.from_kwargs("headers", None, kwargs, rm = True)
        label_subplots = aux.from_kwargs("label_subplots", True, kwargs, rm = True)
        start_character = aux.from_kwargs("labeltype", "A", kwargs, rm = True)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        color = aux.from_kwargs("color", "Lightgray", kwargs, rm = True)
        show = aux.from_kwargs("show", True, kwargs, rm = True)

        figsize = aux.from_kwargs("figsize", None, kwargs, rm = True)
        fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize)


        # setup Coords
        Coords = gx.AxesCoords(fig, axs, (nrows, ncols))

        idx = 0
        for assay in headings: 
            try: 
                tmp_df = self._data.query(query.format(ref_col = ref_col, q = assay))
                # now plot a new bar chart 
                subplot = Coords.subplot()

                tmp_df.plot.bar(
                                x = x, y = y, 
                                ax = subplot,

                                edgecolor = aux.from_kwargs("edgecolor", "black", kwargs),
                                linewidth = aux.from_kwargs("borderwidth", 0.3, kwargs, rm = True),
                                color = color,
                                rot = aux.from_kwargs("rot", 0, kwargs, rm=True),
                                legend = False,
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

        print("<<< Static Check >>>")

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

if __name__ == '__main__':

    files = ["Example Data/28S.csv", "Example Data/28S_again.csv", "Example Data/actin.csv",
    "Example Data/actin2.csv",
    "Example Data/actin3.csv",
     "Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]
    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    # now generate figure
    a = PreviewResults("static")
    b = PreviewResults()

    # predefined pipeline use
    pipe = Pipes.Basic()
    pipe.link(files)
    pipe.add_normalisers(files[:2])
    pipe.replicates(6)
    pipe.names(groupnames)

    a.params(
          frame = False, labeltype = "A", show = True
    )
    b.params(
        show = True      
    )

    print()
    print("==== Multiple (ALL) Samples ====")
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
    print("==== Multiple (7) Samples ====")
    pipe.link(files[:7])
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
    print("==== Multiple (2) Samples ====")
    pipe.link(files[:2])
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





