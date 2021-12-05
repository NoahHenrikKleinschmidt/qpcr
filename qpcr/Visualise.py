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
import numpy as np 

class Chart:
    """
    Setup charts to be either interactive (produced with plotly) or static (matplotlib).
    mode = "interactive" or "static" can be set.
    """
    def __init__(self, mode = "interactive"):
        self._MODE = mode
        self._data = None
        self._PARAMS = {}

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

    def params(self, params:dict):
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

class PreviewResults(Chart):
    """
    Generates a bar chart for each Assay separately, and produces a subplots figure.
    """
    def __init__(self):
        super().__init__()
        if self._MODE == "interactive":
            self._CORE = _plotly_PreviewResults
        else:
            self._CORE = _matplotlib_PreviewResults
    

    def plot(self, **kwargs):
        """
        Plots the data into figure
        """
        kwargs = self.update_params(kwargs)
        self._CORE = self._CORE(self)
        fig = self._CORE.plot(**kwargs)
        return fig 

    


class _plotly_PreviewResults(PreviewResults):
    """
    Interactive PreviewChart
    """
    def __init__(self, master):
        self._MASTER = master
        self._data = self._MASTER.get()
        

    def plot(self, **kwargs):
        """
        Generates a PreviewResults Figure
        """
        # setup reference column and figure subplots
        ref_col = aux.from_kwargs("key", "assay", kwargs, rm=True)
        ncols, nrows = aux.from_kwargs("subplots", 
                                        gx.make_layout(self._data, ref_col), 
                                        kwargs, rm=True
                                       )
        
        headings = aux.sorted_set(self._data[ref_col])
        speclist = gx.make_speclist(nrows, ncols, "xy")
        fig = make_subplots(nrows, ncols, specs = speclist, subplot_titles = headings)

        x = aux.from_kwargs("x", self._MASTER._default_x, kwargs, rm = True)
        y = aux.from_kwargs("y", self._MASTER._default_y, kwargs, rm = True)
        sterr = aux.from_kwargs("sterr", self._MASTER._default_sterr, kwargs, rm = True)

        query = "{ref_col} == '{q}'" if isinstance(headings[0], str) else "{ref_col} == {q}"

        row, col = 1, 1
        for assay in headings:
            tmp_df = self._data.query(query.format(ref_col = ref_col, q = assay))

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

            if col == ncols:
                col = 1
                row += 1
            else:
                col += 1
        
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

        return fig 


class _matplotlib_PreviewResults(PreviewResults):
    """
    Static PreviewChart
    """
    def __init__(self, master:Chart):
        super().__init__()
        self._MASTER = master




if __name__ == '__main__':

    files = ["Example Data/28S.csv", "Example Data/actin.csv", "Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]
    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    # manual use
    # analysers = []

    # reader = qpcr.SampleReader()
    # reader.replicates(6)
    # reader.names(groupnames)

    # analyser = qpcr.Analyser()

    # for file in files: 

    #     sample = reader.read(file)
    #     res = analyser.pipe(sample)

    #     analysers.append(res)

    # normaliser = qpcr.Normaliser()
    # normaliser.link(normalisers = analysers[:2])
    # normaliser.link(samples = analysers[2:])

    # normaliser.normalise()
    
    # result = normaliser.get()

    # predefined pipeline use
    pipe = Pipes.DeltaDeltaCt()
    pipe.add_assays(files[2:])
    pipe.add_normalisers(files[:2])
    pipe.replicates(6)
    pipe.names(groupnames)

    pipe.run()

    result = pipe.get(kind = "obj")


    # now generate figure
    a = PreviewResults()
    a.params(
        dict(
            template = "plotly_dark"
        )
    )
    a.link(result)

    a.plot()



