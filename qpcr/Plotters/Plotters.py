"""
This module is designed for streamlined data visualisation of the qpcr generated results
It is designed to work directly with `qpcr.Results` objects. 

## `Static` vs `Interactive` Figures
----

The `Plotters` are designed to produce two kinds of figures each, either a `"static"` or  an `"interactive"` figure. 
The type of figure a specific Plotter should produce has to be specified using the `mode` argument. 

### Static Figures
_Static_  figures are made using `matplotlib` and they will open through whatever backend your matplotlib configuration 
as specified. Static figures are primarily designed for printing into labjournals and offer a greater flexibility with 
style customizibility (you can use `seaborn` styles for instance, or the matplotlib `rcparams` to style your figures). 

### Interactive Figures
_Interactive_  figures are made using `plotly` and they will open in your browser. Interactive Figures are primarily designed
for cases where your figures contain a lot of data so having a static view on them might be insufficient. Interactive figures
offer `plotly`'s native features like zooming, cropping, size-adjustments and so forth. It comes at the price of less flexibility
with regard to styling. You can set plotly `templates`, but that's about it. However, also interactive figures are perfectly 
adequate for your labjournal, and you may prefer using these for their dynamic figure size adjustments directly from your browser. 

### Plotting `kwargs` 
Both Static and Interactive Figures support a variety of keyword arguments that allow you to customise many of their 
characteristics and their underlying data handling. You can check which kwargs are passable to each type of figure in 
the documentation of each Plotter. 
"""

import qpcr.__init__ as qpcr
import qpcr.Pipes
import qpcr._auxiliary.graphical as gx
import qpcr._auxiliary as aux 
import qpcr._auxiliary.defaults as defaults
import qpcr._auxiliary.warnings as wa
import pandas as pd 
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly
import seaborn as sns 

# setup default settings for charts

_default_static_PreviewResults = defaults.static_PreviewResults
_default_interactive_PreviewResults = defaults.interactive_PreviewResults

_default_static_ReplicateBoxPlot = defaults.static_ReplicateBoxPlot
_default_interactive_ReplicateBoxPlot = defaults.interactive_ReplicateBoxPlot


# Concept: 
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
    (not for End-User usage)

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).    
    """
    def __init__(self, mode = None):
        self._default_params = None
        self._PARAMS = {} #: _PARAMS will store the plotting parameter kwargs from both _default_params and any additional user specified parameters
        self._Results = None
        self._data = None
        self._id = type(self).__name__
        self._MODE = "interactive" if mode is None else mode
        self._fig = None
        self._set_plot()

    def link(self, Results:(qpcr.Results or pd.DataFrame)):
        """
        Links a Results object or pandas DataFrame of the same architecture
        as one handled by a Results object to the Plotter. 
        Note, that this will replace any previously linked data!

        Parameters
        ----------
        Results : qpcr.Results or pd.DataFrame
            A qpcr.Results Object or a pandas DataFrame of the same architecture.
        """
        
        if aux.same_type(Results, qpcr.Results()):
            self._Results = Results
            self._data = self._Results.stats()
        elif isinstance(Results, pd.DataFrame):
            self._Results = None
            self._data = Results
        else:
            wa.HardWarning("Plotter:unknown_data", obj = Results)

    def plot(self, **kwargs):
        """
        Generate Figure

        Parameters
        ----------
        **kwargs
            Any arbitrary keyword arguments to be passed to the plotting method
        
        Returns
        -------
        fig
            A figure object, either from matplotlib or plotly
        """
        total_kwargs = self.update_params(kwargs)
        fig = self._plot(**total_kwargs)
        self._fig = fig
        return fig

    def id(self, id:str = None):
        """
        Set a unique Id, by default the Classname will be used.

        Parameters
        ----------
        id : str
            A unique identifier for the Plotter object
        
        Returns
        -------
        id 
            The Plotter Id (if no id keyword was entered)
        """
        if id is not None:
            self._id = id
        else:
            return self._id

    def get(self):
        """
        Returns the DataFrame used for plotting

        Returns
        -------
        data
            A pandas DataFrame containing the data underlying the plot.
        """
        return self._data

    def params(self, **params):
        """
        Set default parameters for plotting (will be forwarded to **kwargs)
        Returns default parameters if no new parameters are added.

        Parameters
        ----------
        **params
            Any arbitrary keyword arguments to be passed on to the pre-set plotting kwargs.
            This will either set the default parameters (if none has been specified up to the function call)
            or will update the parameters. In case of keyword duplications any OLD key : value pairs will be 
            OVERWRITTEN by new ones. 
        
        Returns
        -------
        params : dict
            A dictionary of the new pre-set plotting parameters
        """
        if params != {} and self._PARAMS == {}:
            self._PARAMS = params
        elif params != {}:
            self.update_params(params, store = True)
        return self._PARAMS

    def update_params(self, kwargs, supersede = True, store = False):
        """
        Appends pre-set parameters to kwargs. 
        
        Parameters
        ----------
        kwargs : dict
            A dictionary of arbitrary keywords for plotting 
        supersede : bool 
            In case of key duplications: Old values will be replaced with new ones (if supersede = True, default), 
            or keep old ones (supersede = False).
        store : bool
            The newly set plotting parameter dictionary will be stored by the plotting object as new default (if store = True, default is False).
        
        Returns
        -------
        params : dict
            A dictionary of updated plotting parameters
        """
        if supersede:
            kwargs = dict(self._PARAMS, **kwargs)
        else:
            kwargs = dict(kwargs, **self._PARAMS)
        
        if store:
            self._PARAMS = kwargs

        return kwargs

    def save(self, filename, **kwargs):
        """
        Saves the figure to a file. 
        Figures are either saved using `plt.savefig` or `plotly.offline.plot`.

        Parameters
        ----------
        filename : str 
            A filename to save the figure to.
        
        **kwargs
            Any arbitrary keyword arguments for the respective figure saving method. 
        """
        if self._fig is None:
            wa.SoftWarning("Plotter:no_fig_yet")
        else:
            if self._MODE == "static":
                self._fig.savefig(filename, bbox_inches = 'tight', **kwargs)
            elif self._MODE == "interactive":
                plotly.offline.plot(self._fig, filename=filename, **kwargs)

    def suffix(self):
        """
        Returns
        ------- 
        suffix : str
            The appropriate file suffix for save (html for "interactive" or jpg for "static" figures)
        """
        suffix = "jpg" if self._MODE == "static" else "html"
        return suffix

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

        headers = aux.sorted_set(data[ref_col])

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

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).


    Plotting Kwargs
    ----
    
    #### `"static"` Kwargs
    Static PreviewResults figures accept the following kwargs:
    
    |   Argument  |  Description    |  Example    |
    | ---- | ---- | ---- |
    |  show : `bool`    |  Whether or not to show the figure    |  `show = True` (default)   |
    |   figsize : `tuple`   |  The figure size    | `figsize = (10, 4)`     |
    | title : `str`   |  The overall figure title   | `title = "Today's results"    |
    | xlabel : `str`| The x-axis label of each subplot | `xlabel = "Conditions"` |
    | ylabel : `str`| The y-axis label of each subplot | `ylabel = "Mean $\Delta\Delta Ct$"` |
    |    headers : `list` |  A list of titles for each subplot in the preview figure    | `headers = ["transcript A", "transcript B"]`     |
    |  label_subplots  : `bool`   |   Add each subplot with A, B, C ... (if True, default)   | `label_subplots = True` (default)     |
    | labeltype : `str`| The starting character for subplot labelling. By default an `"A"`. | `labeltype = "a"` |
    |   frame   : `bool` |  Show left and top spines of subplots (if True)    | `frame = False` (default)     |
    |   color : `str or list`   | The fillcolor for the individual bars   | `color = "yellow"`     |
    |   edgecolor : `str or list`   | The edgecolor for the individual bars   | `edgecolor = "black"`     |
    |   edgewidth : `float`   |  The width of the edge of individual bars  | `edgewidth = 0.5`     |
    |  ecolor : `str`    |   The color of errorbars   |  `ecolor = "orange"`    |
    |  **kwargs    | Any additional kwargs that can be passed to the `matplotlib`-backend pandas `.plot.bar()` API.     |      |

    <br></br>
    #### `"interactive"` Kwargs
    Interactive PreviewResults figures accept the following kwargs:

    |   Argument  |  Description    |  Example    |
    | ---- | ---- | ---- |
    |  show : `bool`    |  Whether or not to show the figure    |  `show = True` (default)   |
    | title : `str`   |  The overall figure title   | `title = "Today's results"`    |
    | xlabel : `str`   |  The x axis label   | `xlabel = "My super qPCR samples"`    |
    | ylabel : `str`   |  The y axis label   | `ylabel = "Mean of ddCt"`    |
    |  height : `int`   |   Height of the figure   | `height = 50`    |
    |  width : `int`   |   Width of the figure   | `width = 50`    |
    |  padding : `float or tuple`   |   Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.   | `padding = 0.2`    |
    |  template : `str`   | The `plotly` template to use. Check out available templates [here](https://plotly.com/python/templates/).     | `template = "plotly_dark"`    |
    |    headers : `list` |  A list of titles for each subplot in the preview figure    | `headers = ["transcript A", "transcript B"]`     |
    | legend_title : `str`    | The title to be displayed above the legend   |  `legend_title = "my assays"`   |
    |  hoverinfo : `str`   | The type of hoverinfo to display. By default just `"y"`. Learn more about plotly hoverinfo [here](https://plotly.com/python/hover-text-and-formatting/). Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`    |
    |  **kwargs    | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Bar()`.     |      |


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

        edgecolor = aux.from_kwargs("edgecolor", "black", kwargs, rm = True)
        edgewidth = aux.from_kwargs("edgewidth", 0.3, kwargs, rm = True)
        figsize = aux.from_kwargs("figsize", None, kwargs, rm = True)
        fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize)

        # setup Coords
        Coords = gx.AxesCoords(fig, axs, (nrows, ncols))
        idx = 0
        for assay in headings: 
            try: 
                q = query.format(ref_col = ref_col, q = assay)
                tmp_df = data.query(q)
                tmp_df = tmp_df.sort_values("group")


                # now plot a new bar chart 
                subplot = Coords.subplot()

                tmp_df.plot.bar(
                                x = x, y = y, 
                                ax = subplot,

                                edgecolor = edgecolor,
                                linewidth = edgewidth,
                                **kwargs
                            )

                subplot.errorbar(
                                x = tmp_df[x], y = tmp_df[y], 
                                yerr = tmp_df[sterr], 
                                fmt = ".", markersize = 0, capsize = 3, 
                                ecolor = aux.from_kwargs("ecolor", "black", kwargs),
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
                raise e 
                break
        
        plt.tight_layout()

        
        if show:
            plt.show()

        # print("<<< Static Check >>>")

        return fig 

    def _add_subplot_label(self, idx, subplot, start_character):
        """
        Adds A B C ... to upper left corner of a subplot...
        It will start labelling at any arbitrary start_character
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
            padding = aux.from_kwargs("padding", 0.1, kwargs, rm = True)

            if isinstance(padding, (list, tuple)):
                hpad, vpad = padding
            else: 
                hpad, vpad = padding, None

            speclist = gx.make_speclist(nrows, ncols, "xy")
            fig = make_subplots(
                                    rows = nrows, cols = ncols, 
                                    specs = speclist, 
                                    subplot_titles = headers,
                                    horizontal_spacing = hpad,
                                    vertical_spacing = vpad
                            )


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

            # set default axes labels
            xlabel = aux.from_kwargs("xlabel", "", kwargs, rm = True) 
            fig.update_xaxes(title_text = xlabel)
            
            ylabel = aux.from_kwargs("ylabel", "DeltaDeltaCt", kwargs, rm = True) 
            fig.update_yaxes(title_text = ylabel)

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


class ReplicateBoxPlot(Plotter):
    """
    Generate a boxplot figure summary for the input sample replicates.
    

    Parameters
    ----------
    Filter : qpcr.Filters.Filter
        A qpcr.Filters.Filter object. The Filter can also be set using the `filter()` method.

    mode : str
        The plotting mode (either `"interactive"` or `"static"`).


    
    Plotting Kwargs
    ----
    
    #### `"static"` Kwargs
    Static ReplicateBoxPlot figures accept the following kwargs:
    
    |   Argument  |  Description    |  Example    |
    | ---- | ---- | ---- |
    |  show : `bool`    |  Whether or not to show the figure    |  `show = True` (default)   |
    |   figsize : `tuple`   |  The figure size    | `figsize = (10, 4)`     |
    |   subplots : `tuple`   |  A tuple specifying the number of colums and rows (in that order) for the figure    | `subplots = (2, 3)`     |
    | title : `str`   |  The overall figure title   | `title = "My assays"    |
    |  **kwargs    | Any additional kwargs that can be passed to the `seaborn`'s `boxplot()`.     |      |

    <br></br>
    #### `"interactive"` Kwargs
    Interactive ReplicateBoxPlot figures accept the following kwargs:

    |   Argument  |  Description    |  Example    |
    | ---- | ---- | ---- |
    |  show : `bool`    |  Whether or not to show the figure    |  `show = True` (default)   |
    | title : `str`   |  The overall figure title   | `title = "My run"`    |
    |  height : `int`   |   Height of the figure   | `height = 50`    |
    |  width : `int`   |   Width of the figure   | `width = 50`    |
    |  template : `str`   | The `plotly` template to use. Check out available templates [here](https://plotly.com/python/templates/).     | `template = "plotly_dark"`    |
    |  **kwargs    | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Box()`.     |      |



    """
    def __init__(self, Filter = None, mode="interactive"):
        self._setup_default_params(
                                    static = _default_static_ReplicateBoxPlot, 
                                    interactive = _default_interactive_ReplicateBoxPlot
                                )
        # __init__ is going to require default_params to be already set!
        super().__init__(mode=mode)
        self._data = None
        self._Filter = Filter
        self._filter_stats = None
    
    def filter(self, Filter):
        """
        Links a Filter object for which a report figure shall be generated

        Parameters
        ----------
        Filter : qpcr.Filters.Filter
            A qpcr.Filters.Filter object.

        """
        self._Filter = Filter

    def link(self, Assay:qpcr.Assay):
        """
        Links an Assay object to the BoxPlotter.
        This will simply add the Ct column to the current overall data!

        Parameters
        ----------
        Assay : qpcr.Assay
            A qpcr.Assay object.
        """
        if aux.same_type(Assay, qpcr.Assay()):
            self._Results = Assay
            data = self._Results.get()
        else:
            wa.HardWarning("Plotter:unknown_data", obj = Assay)

        # add itentifier column
        data["assay"] = [self._Results.id() for i in range(len(data))]

        if self._data is None: 
            self._data = data
        else:
            self._data = pd.concat([self._data, data], ignore_index=True)
        

    def _interactive_plot(self, **kwargs):
        """
        Generates an interactive Boxplot summary of the input Ct values
        """

        show = aux.from_kwargs("show", True, kwargs, rm = True)
        template = aux.from_kwargs("template", None, kwargs, rm=True)
        title = aux.from_kwargs("title", None, kwargs, rm=True)

        data = self._data
        
        # currently unused (needs to be integrated in the future..., see TODO below)
        # filter_stats = self._Filter.get_stats()

        groups = aux.sorted_set(data["group"])
        group_names = aux.sorted_set(data["group_name"])

        fig = go.Figure()

        fig.update_layout(
                            boxmode="group", 
                            height = aux.from_kwargs("height", None, kwargs, rm = True), 
                            width = aux.from_kwargs("width", None, kwargs, rm = True),
                            template = template, 
                            title = title,
                        )

        for group, name in zip(groups, group_names):
            tmp_df = data.query(f"group == {group}")
            
            # future TODO: 
            # currently the inclusion-range is NOT yet represented.
            # It appears as if boxplot cannot use precomputed mean-sd when x and y are specified... (hypothesis of mine)
            # so we'd need an alternative way of plotting this stuff so we can include it...
             
            fig.add_trace(
                            go.Box(
                                x = tmp_df["assay"],
                                y = tmp_df["Ct"],
                                name = name,
                                hoverinfo = "y+name",
                                **kwargs,
                            ),
                        )

        if show:
            fig.show()

        return fig 


    def _static_plot(self, **kwargs):
        """
        Generates a static Boxplot summary of the input Ct values
        """

        data = self._data
        
        show = aux.from_kwargs("show", True, kwargs, rm = True)
        figsize = aux.from_kwargs("figsize", (7,5), kwargs, rm=True)
        title = aux.from_kwargs("title", None, kwargs, rm=True)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        ncols, nrows = aux.from_kwargs("subplots", gx.make_layout(data, "assay"), kwargs, rm=True)

        # future TODO: 
        # currently the inclusion-range is NOT yet represented.
        # It appears as if boxplot cannot use precomputed mean-sd when x and y are specified... (hypothesis of mine)
        # so we'd need an alternative way of plotting this stuff so we can include it...
            
        # filter_stats = self._Filter.get_stats()

        # sns.set_style("whitegrid")

        fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize=figsize)
        
        # possibly nrows and ncols should be switched...
        Coords = gx.AxesCoords(fig, axs, (nrows, ncols))
        
        for assay in aux.sorted_set(data["assay"]):
            
            try: 
                ax = Coords.subplot()
            except: 
                print("We ran out of axes ... ")
                break
            tmp = data.query(f"assay == '{assay}'")

            sns.boxplot(
                        data=tmp, 
                        x = "assay", y = "Ct", hue="group_name", 
                        ax = ax, 
                        **kwargs
                    )

            ax.legend(bbox_to_anchor=(1,1), loc = None).remove()
            ax.set(title=assay, xlabel = "", xticklabels = [],)
            
            if not show_spines:
                sns.despine()

            Coords.increment()

        # add one single legend to the last plot
        ax.legend(bbox_to_anchor=(1,1), loc = None)

        fig.suptitle(title)
        plt.tight_layout()

        if show:
            fig.show()

        return fig 


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
          frame = True, labeltype = "a", show = False
    )
    b.params(
        show = False, template = "plotly"  
    )

    # predefined pipeline use
    pipe = qpcr.Pipes.Basic()
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
        a.plot(show = False)
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
    print("==== Multiple (2) Samples ====")
    pipe.link(files[:2])
    pipe.run()
    result = pipe.get(kind = "obj")
    try:
        a.link(result)
        a.plot(show = False)
    except Exception as e: 
        print(e)
        print("Static Failed!")

    try:
        b.link(result)
        b.plot(show = False)
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