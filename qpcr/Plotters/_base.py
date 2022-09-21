"""
The base classes for the Plotters (FigureClasses) and Wrappers.
"""

import pandas as pd
import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import qpcr._auxiliary.warnings as aw
import qpcr._auxiliary.graphical as gx
import qpcr.main as main

import matplotlib.pyplot as plt
import plotly 
from plotly.subplots import make_subplots


logger = aux.default_logger()

class _Base(aux._ID):
    """
    The Base class for Plotters and Wrappers
    """
    __slots__ = ["_MODE"]
    def __init__(self, mode = None ):
        super().__init__()
        self._MODE = defaults.plotmode if mode is None else mode
        self._id = self.__class__.__name__
    
    def suffix(self):
        """
        Returns
        ------- 
        suffix : str
            The appropriate file suffix for save (html for "interactive" or jpg for "static" figures)
        """
        suffix = "jpg" if self._MODE == "static" else "html"
        return suffix

class Plotter(_Base):
    """
    A superclass that handles Data Linking and Parameter setup for FigureClasses
    (not for End-User usage)

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).    
    """
    __slots__ = "_default_params", "_data", "_PARAMS", "_Results", "_data", "_fig", "_static_default", "_interactive_default", "_default_x", "_default_y", "_default_sterr"
        
    def __init__(self, mode = None):
        super().__init__(mode)
        
        self._default_params = None
        self._PARAMS = {}

        self._obj = None # the data source object
        self._data = None # the data to plot (df)
        
        self._fig = None
        self._set_plot()

    def clear(self):
        """
        Will clear the currently stored data and figure
        """
        self._data = None
        self._obj = None
        self._fig = None

    def link(self, obj):
        """
        Links a qpcr object as data source.

        Parameters
        ----------
        obj : some data source object
            A qpcr object.
        """
        print( "Define here a dedicated method to linking the respective data source object(s)" )
        self._obj = obj

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
            self._PARAMS = dict(params)
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
    
    def reset_params(self):
        """
        Reset to the pre-set default plotting parameters.
        """
        params = { "static": self._static_default, "interactive": self._interactive_default }
        self._PARAMS = params[self._MODE]

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
            e = aw.PlotterError( "no_fig_yet" )
            logger.info( e )
        else:
            if self._MODE == "static":
                self._fig.savefig(filename, bbox_inches = 'tight', **kwargs)
            elif self._MODE == "interactive":
                plotly.offline.plot(self._fig, filename=filename, **kwargs)

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
        options = { 
                        "static": (self._static_plot, self._static_default), 
                        "interactive": (self._interactive_plot,self._interactive_default) 
                    }
        prev_default = self._default_params
        self._plot, self._default_params = options[self._MODE]
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


    @staticmethod
    def _add_subplot_label(idx, subplot, start_character):
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

    @staticmethod
    def _axeslabels(kwargs):
        """
        Pops the kwargs: xlabel, ylabel.
        """
        xlabel = kwargs.pop("xlabel", None)
        ylabel = kwargs.pop("ylabel", None)
        return xlabel, ylabel

    @staticmethod
    def _get_labels_and_spines_and_rot(kwargs):
        """
        Pops the kwargs: label_subplots, start_character, frame, and rot.
        """
        label_subplots = kwargs.pop("label_subplots", True)
        start_character = kwargs.pop("labeltype", "A")
        show_spines = kwargs.pop("frame", True)
        rot = kwargs.pop("rot", None)
        return label_subplots,start_character,show_spines, rot

    @staticmethod
    def _get_title_and_show(kwargs):
        """
        Pops the kwargs: title and show.
        """
        title = kwargs.pop("title", "Results Preview")
        show = kwargs.pop("show", False)
        return title,show

    def _get_interactive_setup_kwargs(self, kwargs):
        """
        Pops the kwargs: padding (is unpacked into hpad,vpad), height, width, template, hoverinfo.
        """
        hpad, vpad = self._get_padding(kwargs)
        height, width = self._get_height_and_width(kwargs)
        template = kwargs.pop("template", "plotly_white")
        hoverinfo = kwargs.pop("hoverinfo", "y")
        return hpad, vpad, height, width, template, hoverinfo

    def _get_height_and_width(self, kwargs):
        """
        Pops the kwargs: height, width
        """
        height = kwargs.pop("height", None)
        width = kwargs.pop("width", None)
        return height, width

    @staticmethod
    def _set_xtick_rotation( subplot, rot ):
        """
        Sets the rotation of x-axis ticks.
        """

        # adjust anchor location
        if rot > 0:
            align = "left"
        elif rot < 0:
            align = "right"
        else:
            align = "center"

        # set rotation
        plt.setp( subplot.xaxis.get_majorticklabels(), rotation = -rot, ha=align, rotation_mode="anchor") 

    @staticmethod
    def _despine( ax ):
        """
        Performs despining but makes remaining x and y axes a bit thicker...
        """
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_linewidth(1.05)
        ax.spines["bottom"].set_linewidth(1.05)

    @staticmethod
    def _get_groups_and_names(data):
        """
        Gets unique groups and names from the data's correspondingly named columns.
        """
        groups = data["group"].unique()
        group_names = data["group_name"].unique()
        return groups, group_names

    @staticmethod
    def _get_padding(kwargs):
        """
        Pops the padding for interactive figures and unravels it to hpad, vpad
        """
        padding = kwargs.pop("padding", 0.1)
        if isinstance(padding, (list, tuple)):
            hpad, vpad = padding
        else: 
            hpad, vpad = padding, None
        return hpad,vpad

    @staticmethod
    def _setup_static_figure(ncols, nrows, title, kwargs):
        """
        Set up a static figure
        """
        figsize = kwargs.pop("figsize", None)
        fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize)
        fig.suptitle( title )

        # setup Coords
        Coords = gx.AxesCoords(fig, axs, (nrows, ncols))
        Coords.autoincrement()
        return fig,Coords

    @staticmethod
    def _setup_interactive_figure(ncols, nrows, xlabel, ylabel, title, headers, hpad, vpad, height, width, template):
        """
        Set up an interactive figure
        """
        speclist = gx.make_speclist(nrows, ncols, "xy")
        fig = make_subplots(
                                rows = nrows, cols = ncols, 
                                specs = speclist, 
                                subplot_titles = headers,
                                horizontal_spacing = hpad,
                                vertical_spacing = vpad
                        )


        Coords = gx.AxesCoords(fig, [], (ncols, nrows))
        Coords.autoincrement()

        fig.update_layout(
                            title = title, 
                            height = height, 
                            width = width,
                            margin = {"autoexpand": True, "pad" : 0, "b": 1, "t": 50}, 
                            autosize = True, 
                            template = template,
                            # legend = {"title" : kwargs.pop("legend_title", ref_col),
                            # },
                        )

        # set default axes labels
        fig.update_xaxes(title_text = xlabel)
        fig.update_yaxes(title_text = ylabel)
        return fig,Coords


# The Wrapper essentially defines public methods from Plotter but 
# replaces self with self._Plotter instead...
class Wrapper(_Base):
    """
    A superclass that allows to make wrappers for multiple Plotters 
    (not for End-User usage).

    Parameters
    ----------
    kind : str
        The kind of Plotter to call. This can be any of the wrapped 
        Plotters, e.g. `kind = "GroupBars"`.   
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly). 

    """
    __slots__ = "_Plotter"

    def __init__(self, kind : str, mode : str):
        self._Plotter = None
        self._id = type(self).__name__

        # Application of this would work as:
        # plotters = {
        #                 "some_Plotter" : some_Plotter, 
        #         }
            
        # self._Plotter = plotters[ kind ]
        # self._Plotter = self._Plotter( mode = mode )
    
    def plotter(self):
        """
        Returns
        -------
        plotter
            The currently used core Plotter
        """
        plotter = self._Plotter
        return plotter


    def link(self, obj):
        """
        Links a data obj as data source.

        Parameters
        ----------
        obj 
            A qpcr obj to use as data source.
        """
        plotter = self._Plotter
        plotter.link( obj )


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
        fig = self._Plotter.plot( **kwargs )
        return fig 

    def get(self):
        """
        Returns the DataFrame used for plotting

        Returns
        -------
        data
            A pandas DataFrame containing the data underlying the plot.
        """
        return self._Plotter.get()

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
        params = self._Plotter.params( **params )
        return params

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
        params = self._Plotter.update_params( kwargs, supersede = supersede, store = store )
        return params

    def reset_params(self):
        """
        Reset to the pre-set default plotting parameters.
        """
        self._Plotter.reset_params()

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
        self._Plotter.save( filename = filename, **kwargs )
    
    def suffix(self):
        """
        Returns
        ------- 
        suffix : str
            The appropriate file suffix for save (html for "interactive" or jpg for "static" figures)
        """
        suffix = self._Plotter.suffix()
        return suffix



# =============================================================================
# Define some super classes for qpcr-object specific FigureClasses
# =============================================================================

class ResultsPlotter(Plotter):
    """
    A super class for FigureClasses that work with Results objects.
    """
    __slots__ = ["_rep_data"]
    def __init__(self, mode : str ):
        super().__init__(mode)
        self._rep_data = None # replicate data 

    def link(self, obj:(main.Results or pd.DataFrame)):
        """
        Links a qpcr object or pandas DataFrame of the same architecture
        as one handled by a Results object to the Plotter. 
        Note, that this will replace any previously linked data!

        Parameters
        ----------
        obj : qpcr.Results or pd.DataFrame
            A qpcr.Results object or a pandas DataFrame of the same architecture.
        """
        
        if isinstance(obj, main.Results):
            self._obj = obj
            self._data = self._obj.stats()
            self._rep_data = self._obj.get()
        elif isinstance(obj, pd.DataFrame):
            self._obj = None
            self._data = obj
        else:
            e = aw.PlotterError( "unknown_data", obj = obj )
            logger.critical( e )
            raise e 

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

    def _prep_properties(self, kwargs):
        """
        Setup ncols, nrows, subplot titles (headers), x, y, and sterr, columns for figure...
        """
        self._setup_default_plot_cols()
        data = self.get()

        # setup reference column and figure subplots
        ref_col = kwargs.pop( "key", "assay" )
        ncols, nrows = kwargs.pop("subplots", gx.make_layout(data, ref_col) )

        # transposing currently not supported...
        # if kwargs.pop("transpose", False, kwargs):
        #     ncols, nrows = nrows, ncols

        headers = aux.sorted_set( data[ref_col] )

        x, y, sterr = self._get_xysterr(kwargs)

        xlabel, ylabel = self._axeslabels(kwargs)

        # the query to be used if group or group_name is ref_col
        query = "{ref_col} == '{q}'" if isinstance(headers[0], str) else "{ref_col} == {q}"

        return ref_col,ncols,nrows, headers, x, y, sterr, query, xlabel, ylabel

    def _get_xysterr(self, kwargs):
        """
        Pops the kwargs: x, y, and sterr.
        """
        self._setup_default_plot_cols()
        x = kwargs.pop( "x", self._default_x )
        y = kwargs.pop( "y", self._default_y )
        sterr = kwargs.pop( "sterr", self._default_sterr )
        return x,y,sterr


class AssayPlotter(Plotter):
    """
    A superclass for all FigureClasses that work with Assay objects.
    """
    def __init__(self, mode : str ):
        super().__init__(mode)

    def link(self, obj:main.Assay):
        """
        link a qpcr.Assay object to the Plotter.

        Parameters
        ----------
        obj : qpcr.Assay
            A qpcr.Assay object.
        """
        if isinstance( obj, main.Assay) :
            self._obj = obj
            self._data = self._obj.get( copy = True )
        else:
            e = aw.PlotterError("unknown_data", obj = obj )
            logger.critical( e )
            raise e
        