"""
.. _qpcr.Plotters:

This module is designed for streamlined data visualisation directly from ``qpcr`` classes. 
To that end, it defines specific ``Plotter`` classes that are dedicated to producing one specific kind of figure each.
All ``qpcr`` classes that support visualisation from their data have a built-in method to call their appropriate Plotter directly
(such as the ``qpcr.Results.preview`` method that calls the ``PreviewResults`` wrapper). At this point: a "wrapper" is a wrapper for multiple Plotters.
You may choose to setup your own Plotter instead of relying on the built-in methods if you have to visualize many different objects, 
because like this you save the repeated plotter creation by the built-in methods. Also, manual setup is required for Plotter adding to the :ref:`qpcr.Pipes <Pipes>`.


`Static` vs `Interactive` Figures
=================================

The `Plotters` are designed to produce two kinds of figures each, either a ``"static"`` or  an ``"interactive"`` figure. 
The type of figure a specific Plotter should produce has to be specified using the ``mode`` argument. While the two 
visualizations aim to be as close as possible to one another, they are naturally not identical, however. 
Also they offer different kinds of costumization options.


Static Figures
--------------

**Static**  figures are made using ``matplotlib`` and they will open through whatever backend your matplotlib configuration 
as specified. Static figures are primarily designed for printing into labjournals and offer a greater flexibility with 
style customizibility (you can use ``seaborn`` styles for instance, or the matplotlib ``rcparams`` to style your figures). 


Interactive Figures
-------------------

**Interactive** figures are made using ``plotly`` and they will open in your browser. Interactive Figures are primarily designed
for cases where your figures contain a lot of data so having a static view on them might be insufficient. Interactive figures
offer ``plotly``'s native features like zooming, cropping, size-adjustments and so forth. It comes at the price of less flexibility
with regard to styling. You can set plotly ``templates``, but that is about the extend of it. However, also interactive figures are perfectly 
adequate for your labjournal, and you may prefer using these for their dynamic figure size adjustments directly from your browser. 

You can trigger interactive figures either using the ``mode`` argument (``mode = "interactive"``) when creating a Plotter or calling a ``plot`` method of one of the non-Plotter classes.
Or you can set the plotting mode globally to interactive using:

.. code-block:: python

    qpcr.Plotters.interactive()


Plotting ``kwargs``
===================

Both Static and Interactive Figures support a variety of keyword arguments that allow you to customise many of their 
characteristics and their underlying data handling. In fact, the entire styling is carried out via these arguments (a relic from earlier code phases).
You can check which kwargs are passable to each type of figure in 
the documentation of each Plotter. 


Visualizing your data
=====================

There are a number of ways to visualize data in ``qpcr``.


Manual setup (the classic way)
--------------------------------

Just as with the ``qpcr.Analyser`` we can *set up* a ``qpcr.Plotter``, then ``link`` the object we want to visualize from, and call the Plotter's ``plot`` method.
We can pass all additional styling arguments either to the ``plot`` method or use the ``params`` method to set them up beforehand.

Let us for instance visualise the raw Ct values from a ``qpcr.Assay`` we just loaded from a datafile. 
The dedicated Plotter for this is the ``ReplicateBoxPlot``.
We can set it up like this:

.. code-block:: python

    from qpcr.Plotters import ReplicateBoxPlot 

    # by default it will produce a "static" figure by default
    myboxplot = ReplicateBoxPlot()

    # now we link the assay
    myboxplot.link( myassay )

    # and add some graphing paramters
    my_params = { 
                    "color" : "Reds", # we want to change the default colormap to red tones
                    "title" : "my customized figure",
                    "style" : "darkgrid", # and use the seaborn darkgrid style 
                }
    myboxplot.params( my_params )

    # and now visualise the figure
    fig = myboxplot.plot()

    # alternatively: fig = myboxplot.plot( **my_params ) 


This is a rather tedious way of setting up for just a quick look at our Assay, but it may be worth it if we have to repeat this many times...

The Built-in shortuct
----------------------
For more convenience, however, the ``qpcr.Assay`` lets us directly visualize its Ct values through the ``boxplot`` method. We can call it easily like:

.. code-block:: python

    fig = myassay.boxplot( **my_params )


This will perform the setup and linking to a new instance of ``ReplicateBoxPlot``. 



The ``qpcr.plot`` function
--------------------------

As a generic way to visualize data, ``qpcr`` offers the stand-alone ``plot`` function that accepts a data-containing object and will call it's built-in plotting method.
Hence, this is just another way of calling the object's built-in plotting method. It may be convenient for you if you are used to the R plot API (and was inspired by it).
However, it is important to note that you cannot add a `Plotter` to the ``plot`` function, only data-containing classes like the ``qpcr.Assay`` or ``qpcr.Results`` support the generic ``plot`` function!
Thus, a final way of obtaining the exact same figure as before, is using:

.. code-block:: python

    fig = qpcr.plot( myassay, **my_params )


What does the output look like? Here is a possible figure output for a `ReplicateBoxPlot` in static mode:

.. image:: ../../docs/source/resources/boxplot.png
    :align: center

Customizing your Figures
-------------------------

Let us work through a quick example of creating a customized figure from a ``qpcr.Results`` object. 
The ``qpcr.Results`` class works uses the ``PreviewResults`` class to visualise its data.``PreviewResults`` is a wrapper for `AssayBars`, `AssayDots`, `GroupBars`, and `GroupDots`, which are the four native visualizations supported for ``qpcr.Results``.
You can easily call on a `PreviewResults` wrapper using the ``qpcr.Results.preview`` method and specify which kind of visualisation to perform using the ``kind`` argument.

For instance, we might want to visualize our data by having each assay in a separate subplot but visualizing all computed values. This would be a job for the ``AssayDots`` Plotter, so we can set up.

.. code-block:: python

    fig = myresults.preview( kind = "AssayDots" )


Which will return the following figure 


.. image:: ../../docs/source/resources/preview_AssayDots.png
    :align: center


Naturally, we can edit the default settings if we like. Maybe we do not like the *_rel_28S+Actin* in the assay headers. 
We could specify our own headers using the ``headers`` arguments, or simply call the ``drop_rel`` method on our ``qpcr.Results`` object to reduce the assay ids back to their original states.

With a few lines of code we can adjust the figure to something like this:

.. code-block:: python

    # get rid of the _rel_28S+Actin
    result.drop_rel()

    # specify some fancy colors
    mycolors = ["crimson", "black" ]

    fig = result.preview( kind = "AssayDots", 
                          title = "Normalized to 28S+Actin", 
                          style = "ticks",
                          color = mycolors 
                        )


.. image:: ../../docs/source/resources/preview_AssayDots_modified.png
    :align: center


If you wish to edit settings more globally you can either set up your Plotters directly and call their ``params`` method to store dedicated plotting parameters, or you may choose to directly edit the
default plotting settings in the ``qpcr.defaults`` submodule. All the plotting settings are stored as dictionaries there, so you can freely edit these as you like. 
To achieve the same figure as before, we might also do:

.. code-block:: python

    # change default plot settings

    qpcr.defaults.default_preview = "AssayDots"
    qpcr.defaults.static_PreviewDots["style"] = "ticks"
    qpcr.defaults.static_PreviewDots["color"] = ["crimson", "black" ]

    fig = results.preview( title = "Normalized to 28S+Actin" )


However, this will now cause **all** of your Results to generate a black and red `AssayDots` figure in `ticks` style when calling their ``preview`` method (unless you specify new parameters manually).
As you will have noticed, the dictionary we edited was not a "qpcr default plotting parameters" dictionary but rather specific for the `static AssayDots` Plotters. All Plotters have two dedicated dictionaries to store their default parameters in, one for static mode, one for interactive.
You can therefore selectively adjust your favourite default settings without having to worry to skew with other figure types. However, `rcParams` will of course, also work just fine globally. 
"""

import qpcr.main as main

import qpcr._auxiliary.graphical as gx
import qpcr.defaults as defaults
import qpcr._auxiliary as aux 
import qpcr._auxiliary.warnings as wa
import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly
import seaborn as sns 
import numpy as np 
import logging

logger = aux.default_logger()

# get default colnames Id + Ct
raw_col_names = defaults.raw_col_names

# Concept: 
# Plotter: 
#       We define a superclass Plotter that will handle linking to data, 
#       linking to default parameters, and setting up default parameters.
#       It also provides a generic plot() method that will call on a FigureClass specific _plot() method...
#       Wether or not to use static or interactive plots is also handled by this class...
#
# Wrapper: 
#       The Wrapper is an equivalent to the Plotter superclass, but allows a facilitated user interface to 
#       several different plotters. It re-defines the public methods of the Plotter superclass and creates 
#       an instance of a specific Plotter and forwards data to that instance. 
#
# FigureClass:
#       A parent class for each type of figure. It contains two potential _plot() 
#       methods, one for interactive one for static plotting... Which one to use is decided based on the plotting mode...
#       Each FigureClass has therefore to link default parameters, then init Plotter (superclass), 
#       and define their own _static_plot() and _interactive_plot() methods! Any additionally required methods can be written as well...
#
# plot function API
#       Different classes can define their own __qplot__ method which makes them plottable via the qpcr.plot function.
#       __qplot__ must return the method to call for plotting (but not call the method itself)!
#       Since the Plotters themselves are at the basis of plotting they do NOT have a __qplot__ method...
#       This is primarily because the Plotters must be set up to be either static or interactive, while 
#       the other classes can decide which plotmode to call on.

class Plotter:
    """
    A superclass that handles Data Linking and Parameter setup for FigureClasses
    (not for End-User usage)

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).    
    """
    __slots__ = "_default_params", "_data", "_PARAMS", "_Results", "_data", "_rep_data", "_id", "_MODE", "_fig", "_static_default", "_interactive_default", "_default_x", "_default_y", "_default_sterr", "__dict__"
        
    def __init__(self, mode = None):
        self._default_params = None
        self._PARAMS = {} #: _PARAMS will store the plotting parameter kwargs from both _default_params and any additional user specified parameters
        self._Results = None
        self._data = None
        self._rep_data = None
        self._id = type(self).__name__
        self._MODE = defaults.plotmode if mode is None else mode
        self._fig = None
        self._set_plot()

    def link(self, Results:(main.Results or pd.DataFrame)):
        """
        Links a Results object or pandas DataFrame of the same architecture
        as one handled by a Results object to the Plotter. 
        Note, that this will replace any previously linked data!

        Parameters
        ----------
        Results : qpcr.Results or pd.DataFrame
            A qpcr.Results Object or a pandas DataFrame of the same architecture.
        """
        
        if aux.same_type(Results, main.Results()):
            self._Results = Results
            self._data = self._Results.stats()
            self._rep_data = self._Results.get()
        elif isinstance(Results, pd.DataFrame):
            self._Results = None
            self._data = Results
        else:
            e = wa.PlotterError( "unknown_data", obj = Results )
            logger.critical( e )
            raise e 

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

    def _prep_properties(self, kwargs):
        """
        Setup ncols, nrows, subplot titles (headers), x, y, and sterr, columns for figure...
        """
        self._setup_default_plot_cols()
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

# The Wrapper essentially defines public methods from Plotter but 
# replaces self with self._Plotter instead...
class Wrapper:
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
    __slots__ = "_Plotter", "_id", "__dict__"

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


    def link(self, Results:(main.Results or pd.DataFrame)):
        """
        Links a Results object or pandas DataFrame of the same architecture
        as one handled by a Results object to the Plotter. 
        Note, that this will replace any previously linked data!

        Parameters
        ----------
        Results : qpcr.Results or pd.DataFrame
            A qpcr.Results Object or a pandas DataFrame of the same architecture.
        """
        plotter = self._Plotter
        plotter.link( Results )


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


class PreviewResults(Wrapper):
    """
    Generate a Preview of all results from all Assays in subplots.

    This is a wrapper for the Plotters: 

    - `AssayBars` ( which was previously `PreviewResults`, and is now the default setting )
    - `GroupBars`
    - `AssayDots`
    - `GroupDots` 

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).

    kind : str
        The kind of Plotter to call. This can be any of the four wrapped 
        Plotters, e.g. `kind = "GroupBars"`.
    """
    def __init__(self, mode : str = None, kind : str = None ):

        super().__init__( kind = kind, mode = mode )
        
        if kind is None: kind = defaults.default_preview
        plotters = {
                        "AssayBars" : AssayBars,
                        "GroupBars" : GroupBars,
                        "AssayDots" : AssayDots,
                        "GroupDots" : GroupDots,  
                }
        
        self._Plotter = plotters[ kind ]
        self._Plotter = self._Plotter( mode = mode )


class AssayBars(Plotter):
    """
    Generate a Preview of all results from all Assays in subplots as Bar Charts.
    Each assay will be shown in a separate subplot.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).


    Plotting Kwargs
    ===============
    
    `"static"` Kwargs
    
        Static AssayBars figures accept the following kwargs:

        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                  | Description                                                                                                            | Example                                       |
        +===========================+========================================================================================================================+===============================================+
        | show: `bool`              | Whether or not to show the figure                                                                                      | `show = True` (default)                       |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | figsize: `tuple`          | The figure size                                                                                                        | `figsize = (10, 4)`                           |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title: `str`              | The overall figure title                                                                                               | `title = "Today's results"                    |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel: `str`             | The x-axis label of each subplot                                                                                       | `xlabel = "Conditions"`                       |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel: `str`             | The y-axis label of each subplot                                                                                       | `ylabel = "Mean Fold Change"`                 |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | rot: `float`              | The rotation of x-axis labels                                                                                          | `rot = 0.3`                                   |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers: `list`           | A list of titles for each subplot in the preview figure                                                                | `headers = ["transcript A", "transcript B"]`  |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | label_subplots: `bool`    | Add each subplot with A, B, C ... (if True, default)                                                                   | `label_subplots = True` (default)             |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | labeltype: `str`          | The starting character for subplot labelling. By default an `"A"`.                                                     | `labeltype = "a"`                             |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | frame: `bool`             | Show left and top spines of subplots (if True)                                                                         | `frame = False` (default)                     |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | color: `str or list`      | The fillcolor for the individual bars                                                                                  | `color = "yellow"`                            |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | style: `str`              | A `seaborn` style to set. Check out available styles [1].                                                              | `style = "darkgrid"`                          |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgecolor: `str or list`  | The edgecolor for the individual bars                                                                                  | `edgecolor = "black"`                         |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgewidth: `float`        | The width of the edge of individual bars                                                                               | `edgewidth = 0.5`                             |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ecolor: `str`             | The color of errorbars                                                                                                 | `ecolor = "orange"`                           |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                | Any additional kwargs that can be passed to the `matplotlib`-backend pandas `.plot.bar()` API.                         |                                               |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+


    
    `"interactive"` Kwargs
    
        Interactive AssayBars figures accept the following kwargs:

        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                    | Example                                       |
        +=============================+================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                              | `show = True` (default)                       |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                       | `title = "Today's results"`                   |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`              | The x axis label                                                                                                                                               | `xlabel = "My super qPCR samples"`            |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                               | `ylabel = "Mean of ddCt"`                     |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                           | `height = 50`                                 |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                            | `width = 50`                                  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                       | `padding = 0.2`                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                               | `template = "plotly_dark"`                    |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                        | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | legend_title : `str`        | The title to be displayed above the legend                                                                                                                     | `legend_title = "my assays"`                  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default just `"y"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                  | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Bar()`.                                                                                     |                                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+

    
    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """
    def __init__(self, mode:str = None):
        self._setup_default_params(
                                    static = defaults.static_PreviewBars, 
                                    interactive = defaults.interactive_PreviewBars
                                )
        # __init__ is going to require default_params to be already set!
        super().__init__(mode)

    def _static_plot(self, **kwargs):
        """
        The static Preview Results Figure
        """
        kwargs = self.update_params(kwargs)
        data = self.get()

        ref_col, ncols, nrows, headings, x, y, sterr, query = self._prep_properties(kwargs)       

        title = aux.from_kwargs("title", "Results Preview", kwargs, rm = True)
        headers = aux.from_kwargs("headers", None, kwargs, rm = True)
        label_subplots = aux.from_kwargs("label_subplots", True, kwargs, rm = True)
        start_character = aux.from_kwargs("labeltype", "A", kwargs, rm = True)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        show = aux.from_kwargs("show", True, kwargs, rm = True)
        rot = aux.from_kwargs("rot", None, kwargs, rm = True)

        palette = gx.generate_palette(kwargs)

        # set a seaborn style
        style = aux.from_kwargs("style", "ticks", kwargs, rm = True)
        sns.set_style( style )

        edgecolor = aux.from_kwargs("edgecolor", "white", kwargs, rm = True)
        edgewidth = aux.from_kwargs("edgewidth", None, kwargs, rm = True)

        if edgewidth is None: 
            edgewidth = aux.from_kwargs("linewidth", 1, kwargs, rm = True) 
        else: 
            aux.from_kwargs("linewidth", None, kwargs, rm = True )

        figsize = aux.from_kwargs("figsize", None, kwargs, rm = True)
        fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize)
        fig.suptitle( title )

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
                                color = palette,
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

                if rot is not None: 
                    align = "center" if rot == 0 else "left"
                    plt.setp( subplot.xaxis.get_majorticklabels(), rotation = -rot, ha=align, rotation_mode="anchor") 

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

    # def _add_subplot_label(self, idx, subplot, start_character):
    #     """
    #     Adds A B C ... to upper left corner of a subplot...
    #     It will start labelling at any arbitrary start_character
    #     """
    #     subplot_label = chr(ord(start_character)+idx)
    #     subplot.annotate(
    #                 xy = (-0.1,1.03), 
    #                 text = subplot_label, 
    #                 xycoords = "axes fraction",
    #                 weight = "bold", fontsize = 12
    #             )

    def _interactive_plot(self, **kwargs):
        """
        The interactive Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self.get()

        try: 
            # setup figure framework variables that have to be removed from 
            ref_col, ncols, nrows, headings, x, y, sterr, query = self._prep_properties(kwargs)
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
                tmp_df = tmp_df.sort_values("group")

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
    

    Note
    ------
    This used to be the core of making filter summary figures. However, this has now
    been replaced by a dedicated `FilterSummary` figure class. Hence, support for linking
    `qpcr.Filters` directly to this figure class will be dropped in a future release!

    Parameters
    ----------
    Filter : qpcr.Filters.Filter
        A qpcr.Filters.Filter object. The Filter can also be set using the `filter()` method.

    mode : str
        The plotting mode (either `"interactive"` or `"static"`).



    Plotting Kwargs
    ===============
    
    `"static"` Kwargs
        
        Static ReplicateBoxPlot figures accept the following kwargs:
    
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | Argument               | Description                                                                      | Example                  |
        +========================+==================================================================================+==========================+
        | show : `bool`          | Whether or not to show the figure                                                | `show = True` (default)  |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | figsize : `tuple`      | The figure size                                                                  | `figsize = (10, 4)`      |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | subplots : `tuple`     | A tuple specifying the number of colums and rows (in that order) for the figure  | `subplots = (2, 3)`      |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | title : `str`          | The overall figure title                                                         | `title = "My assays"     |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | style : `str`          | A `seaborn` style to set. Check out available styles [1].                        | `style = "darkgrid"`     |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | color : `str or list`  | The fillcolor for the boxes.                                                     | `color = "yellow"`       |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | \*\*kwargs             | Any additional kwargs that can be passed to the `seaborn`'s `boxplot()`.         |                          |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+

    
    `"interactive"` Kwargs
        
        Interactive ReplicateBoxPlot figures accept the following kwargs:

        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | Argument          | Description                                                                 | Example                    |
        +===================+=============================================================================+============================+
        | show : `bool`     | Whether or not to show the figure                                           | `show = True` (default)    |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | title : `str`     | The overall figure title                                                    | `title = "My run"`         |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | ylabel : `str`    | The y-axis title                                                            | `ylabel = "Raw Ct value"`  |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | height : `int`    | Height of the figure                                                        | `height = 50`              |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | width : `int`     | Width of the figure                                                         | `width = 50`               |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | template : `str`  | The `plotly` template to use. Check out available templates [2].            | `template = "plotly_dark"` |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | \*\*kwargs        | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Box()`.  |                            |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+


    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    """
    def __init__(self, mode : str = None):
        self._setup_default_params(
                                    static = defaults.static_ReplicateBoxPlot, 
                                    interactive = defaults.interactive_ReplicateBoxPlot
                                )
        # __init__ is going to require default_params to be already set!
        super().__init__(mode=mode)
        self._data = None

    def clear(self):
        """
        Will clear the currently stored Assay data
        """
        self._data = None

    def link(self, Assay:main.Assay):
        """
        Links an Assay object to the BoxPlotter.
        This will simply add the Ct column to the current overall data!
        Hence, repeated linking of Assay objects will add data and NOT replace any existing.

        Parameters
        ----------
        Assay : qpcr.Assay
            A qpcr.Assay object.
        """
        if type(Assay).__name__ == "Assay" :
            self._Results = Assay
            data = self._Results.get( copy = True )
        else:
            e = wa.PlotterError("unknown_data", obj = Assay )
            logger.critical( e )
            raise e

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
        kwargs = self.update_params(kwargs)

        show = aux.from_kwargs("show", True, kwargs, rm = True)
        template = aux.from_kwargs("template", None, kwargs, rm=True)
        title = aux.from_kwargs("title", None, kwargs, rm=True)

        data = self._data
    
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

        # add default ylabel
        ylabel = aux.from_kwargs("ylabel", "Ct", kwargs, rm = True) 
        fig.update_yaxes(title_text = ylabel)

        fig.update_xaxes(showgrid=True)

        for group, name in zip(groups, group_names):
            tmp_df = data.query(f"group == {group}")
            
            
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
        kwargs = self.update_params(kwargs)
        data = self._data
        
        # set a seaborn style
        style = aux.from_kwargs("style", "ticks", kwargs, rm = True)
        sns.set_style( style )

        show = aux.from_kwargs("show", True, kwargs, rm = True)
        figsize = aux.from_kwargs("figsize", (7,5), kwargs, rm=True)
        title = aux.from_kwargs("title", None, kwargs, rm=True)
        palette = gx.generate_palette(kwargs)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        ncols, nrows = aux.from_kwargs("subplots", gx.make_layout(data, "assay"), kwargs, rm=True)

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
                        palette = palette, 
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


class FilterSummary(Plotter):
    """
    Generates a summary figure of the replicate
    Ct values for each assay before and after filtering.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).


    
    `"static"` Kwargs

	    Static FilterSummary figures accept the following kwargs:
    
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | Argument               | Description                                                                      | Example                        |
        +========================+==================================================================================+================================+
        | show : `bool`          | Whether or not to show the figure                                                | `show = True` (default)        |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | figsize : `tuple`      | The figure size                                                                  | `figsize = (10, 4)`            |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | subplots : `tuple`     | A tuple specifying the number of colums and rows (in that order) for the figure  | `subplots = (2, 3)`            |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | title : `str`          | The overall figure title                                                         | `title = "Today's run"         |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | xlabel : `str`         | The x-axis label of each subplot                                                 | `xlabel = "Conditions"`        |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | ylabel : `str`         | The y-axis label of each subplot                                                 | `ylabel = "My Ct values"`      |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | rot : `float`          | The rotation of x-axis labels                                                    | `rot = 0.3`                    |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | frame   : `bool`       | Show left and top spines of subplots (if True)                                   | `frame = False` (default)      |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | style : `str`          | A `seaborn` style to set. Check out available styles [1].                        | `style = "darkgrid"`           |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | color : `str or list`  | The color for the boxes.                                                         | `color = ["yellow", "green"]`  |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | \*\*kwargs             | Any additional kwargs that can be passed to the `seaborn`'s `boxplot()`.         |                                |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+

  
    
    `"interactive"` Kwargs

	    Interactive FilterSummary figures accept the following kwargs:

        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                      | Example                                       |
        +=============================+==================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                                | `show = True` (default)                       |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                         | `title = "Today's run"`                       |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                                 | `ylabel = "Raw Ct values"`                    |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                             | `height = 50`                                 |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                              | `width = 50`                                  |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                         | `padding = 0.2`                               |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | colors : `list`             | List of two colors for before and after filtering boxes.                                                                                                         | `colors = ["red", "green"]`                   |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                                 | `template = "plotly_dark"`                    |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                          | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default `"y+x+name"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                  | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Box()`.                                                                                       |                                               |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+

    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """
    def __init__(self, mode = None ):
        self._setup_default_params(
                                    static = defaults.static_FilterSummary,
                                    interactive = defaults.interactive_FilterSummary
                                )
        super().__init__( mode = mode)
        self._before = main.Results()
        self._after = main.Results()

    def clear(self):
        """
        Clears the pre- and post-filtering Ct value records.
        """
        self._before = main.Results()
        self._after = main.Results()

    def add_before(self, assay : main.Assay ):
        """
        Add a pre-filtered set of Ct values
        from an `qpcr.Assay` object. 

        Parameters
        ----------
        assay : qpcr.Assay
            An `qpcr.Assay` object
        """
        # setup groups and stuff information 
        if self._before.is_empty:
            self._before.setup_cols( assay )
            self._after.setup_cols( assay )

        self._before.add_Ct( assay )
    
    def add_after(self, assay : main.Assay ):
        """
        Add a post-filtered set of Ct values
        from an `qpcr.Assay` object. 

        Parameters
        ----------
        assay : qpcr.Assay
            An `qpcr.Assay` object
        """
        self._after.add_Ct( assay )


    def _interactive_plot(self, **kwargs):
        """
        Generates an interactive Filter Summary Figure
        """
        kwargs = self.update_params(kwargs)

        show = aux.from_kwargs("show", True, kwargs, rm = True)
        template = aux.from_kwargs("template", None, kwargs, rm=True)
        title = aux.from_kwargs("title", None, kwargs, rm=True)
        padding = aux.from_kwargs("padding", 0.1, kwargs, rm = True)
        
        # setup color scheme
        colors = aux.from_kwargs("colors", ["blue", "crimson"], kwargs, rm = True)
        before_color, after_color = colors

        if isinstance(padding, (list, tuple)):
            hpad, vpad = padding
        else: 
            hpad, vpad = padding, None

        # get the data
        before = self._before.get()
        after = self._after.get()

        # get the number of assays to visualse
        assays = [ i for i in before.columns if i not in [ "group", "group_name", raw_col_names[0]] ]
        headers = aux.from_kwargs("headers", assays, kwargs, rm = True)

        # each assay will be on its own plot,
        ncols, nrows = gx.make_layout_from_list( assays )
        
        speclist = gx.make_speclist(nrows, ncols, "xy")
        fig = make_subplots(
                                rows = nrows, cols = ncols, 
                                specs = speclist, 
                                subplot_titles = headers,
                                horizontal_spacing = hpad,
                                vertical_spacing = vpad,
                                # sharex = True
                        )
        
        fig.update_layout(
                            boxmode="group", 
                            height = aux.from_kwargs("height", None, kwargs, rm = True), 
                            width = aux.from_kwargs("width", None, kwargs, rm = True),
                            template = template, 
                            title = title,
                        )

        # add default ylabel
        ylabel = aux.from_kwargs("ylabel", "Ct", kwargs, rm = True) 
        fig.update_yaxes(title_text = ylabel)
        
        Coords = gx.AxesCoords(fig, [], (ncols, nrows))

        show_legend = True
        # now iterate over each assay and make a BoxPlot
        for assay in assays : 
            
            # get before and after filter datasets for each assay
            pre_Cts = before[ ["group", "group_name", assay] ]
            post_Cts = after[ ["group", "group_name", assay] ]
            
            row, col = Coords.get()

            # pre-filter boxes 
            fig.add_trace(
                            go.Box(
                                x = pre_Cts["group_name"],
                                y = pre_Cts[ assay ],
                                name = "Before",
                                hoverinfo = "y+x+name",
                                marker_color = before_color,
                                legendgroup='group1',
                                showlegend = show_legend,
                                **kwargs,
                            ),
                            row, col
                        )
            # post-filter boxes
            fig.add_trace(
                            go.Box(
                                x = post_Cts["group_name"],
                                y = post_Cts[ assay ],
                                name = "After",
                                hoverinfo = "y+x+name",
                                marker_color = after_color,
                                legendgroup='group2',
                                showlegend = show_legend,
                                **kwargs,
                            ),
                            row, col
                        )

            Coords.increment()

            # set show_legend to False after the first assay,
            # to avoid repetitive legends...
            show_legend = False

        fig.update_layout( showlegend = True )
        if show:
            fig.show()

        return fig 


    def _static_plot(self, **kwargs):
        """
        Generates a static Filter Summary Fig
        """
        kwargs = self.update_params(kwargs)

        show = aux.from_kwargs("show", True, kwargs, rm = True)
        figsize = aux.from_kwargs("figsize", (7,5), kwargs, rm=True)
        title = aux.from_kwargs("title", None, kwargs, rm=True)
        ylabel = aux.from_kwargs("ylabel", "Ct", kwargs, rm = True)
        xlabel = aux.from_kwargs("xlabel", "", kwargs, rm = True)
        rot = aux.from_kwargs("rot",None, kwargs, rm = True)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        palette = gx.generate_palette(kwargs)

        # set a seaborn style
        style = aux.from_kwargs("style", "ticks", kwargs, rm = True)
        sns.set_style( style )

        # get the data
        before = self._before.get()
        after = self._after.get()

        # get the number of assays to visualse
        assays = [ i for i in before.columns if i not in [ "group", "group_name", raw_col_names[0]] ]
        headers = aux.from_kwargs("headers", assays, kwargs, rm = True)

        # each assay will be on its own plot
        ncols, nrows = gx.make_layout_from_list( assays )

        ncols, nrows = aux.from_kwargs("subplots", (ncols, nrows), kwargs, rm=True)


        fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize=figsize)
        
        # possibly nrows and ncols should be switched...
        Coords = gx.AxesCoords(fig, axs, (nrows, ncols))

        # set a kind column for grouping in the boxplots
        before[ "kind" ] = "before" 
        after[ "kind" ] = "after"

        show_legend = True
        for assay, header in zip( assays, headers ): 

            # get before and after filter datasets for each assay
            pre_Cts = before[ ["group", "group_name", "kind", assay] ]
            post_Cts = after[ ["group", "group_name", "kind", assay] ]

            # assemble data to one dataframe
            df = pd.concat( (pre_Cts, post_Cts) )

            # get the new subplot
            try: 
                ax = Coords.subplot()
            except: 
                print("We ran out of axes ... ")
                break

            # main box plot
            sns.boxplot(
                        data = df, 
                        x = "group_name",
                        y = assay,
                        hue = "kind",
                        palette = palette, 
                        ax = ax,
                        **kwargs
                    )
            
            # format the subplot
            ax.set(title = header, ylabel = ylabel, xlabel = xlabel)
            if rot is not None and rot != 0 :
                plt.setp( ax.xaxis.get_majorticklabels(), rotation = -rot, ha="left", rotation_mode="anchor") 
            if not show_legend: 
                ax.legend().remove()
            else:
                ax.legend( bbox_to_anchor=(1,1), loc = None, frameon = False )
            show_legend = False

            Coords.increment()

        if not show_spines: 
            sns.despine()

        fig.suptitle(title)

        plt.tight_layout()
        if show: 
            plt.show()

        return fig


class AssayDots(Plotter):
    """
    Generate a Preview of all results from all Assays in subplots, plotting individual values
    in a Dot Plot rather than Bar Plot. Each assay will be shown in a separate subplot.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).


    Plotting Kwargs
    ================
    
    `"static"` Kwargs

	    Static AssayDots figures accept the following kwargs:
    
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                  | Description                                                               | Example                                       |
        +===========================+===========================================================================+===============================================+
        | show : `bool`             | Whether or not to show the figure                                         | `show = True` (default)                       |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | figsize : `tuple`         | The figure size                                                           | `figsize = (10, 4)`                           |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`             | The overall figure title                                                  | `title = "Today's results"                    |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`            | The x-axis label of each subplot                                          | `xlabel = "Conditions"`                       |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`            | The y-axis label of each subplot                                          | `ylabel = "Mean Fold Change"`                 |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | rot : `float`             | The rotation of x-axis labels                                             | `rot = 0.3`                                   |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`          | A list of titles for each subplot in the preview figure                   | `headers = ["transcript A", "transcript B"]`  |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | label_subplots  : `bool`  | Add each subplot with A, B, C ... (if True, default)                      | `label_subplots = True` (default)             |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | labeltype : `str`         | The starting character for subplot labelling. By default an `"A"`.        | `labeltype = "a"`                             |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | frame   : `bool`          | Show left and top spines of subplots (if True)                            | `frame = False` (default)                     |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | violin   : `bool`         | Show symmetric kde of the dots of each group (if True).                   | `violin = False`                              |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | color : `str or list`     | The fillcolor for the individual dots from replicate groups               | `color = "yellow"`                            |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | style : `str`             | A `seaborn` style to set. Check out available styles [1].                 | `style = "darkgrid"`                          |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                | Any additional kwargs that can be passed to the `seaborn`'s `stripplot`.  |                                               |
        +---------------------------+---------------------------------------------------------------------------+-----------------------------------------------+


    
    `"interactive"` Kwargs

	    Interactive AssayDots figures accept the following kwargs:

        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                    | Example                                       |
        +=============================+================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                              | `show = True` (default)                       |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                       | `title = "Today's results"`                   |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`              | The x axis label                                                                                                                                               | `xlabel = "My super qPCR samples"`            |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                               | `ylabel = "Mean of ddCt"`                     |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                           | `height = 50`                                 |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                            | `width = 50`                                  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                       | `padding = 0.2`                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                               | `template = "plotly_dark"`                    |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                        | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | violin   : `bool`           | Show symmetric kde of the dots of each group (if True).                                                                                                        | `violin = False`                              |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default just `"y"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                  | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Violin()`.                                                                                  |                                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+

    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """
    def __init__(self, mode:str = None ):
        self._setup_default_params(
                                    static = defaults.static_PreviewDots, 
                                    interactive = defaults.interactive_PreviewDots
                                )
        # __init__ is going to require default_params to be already set!
        super().__init__(mode)

    def _static_plot(self, **kwargs):
        """
        The static Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self._rep_data

        # get assays to plot
        assays = [ i for i in data.columns if i not in [ "group", "group_name", raw_col_names[0] ] ]

        # generate subplot layout
        ncols, nrows = gx.make_layout_from_list( assays )

        # get kwargs incompatible with the main plotting method
        title = aux.from_kwargs("title", None, kwargs, rm = True)
        headers = aux.from_kwargs("headers", assays, kwargs, rm = True)
        label_subplots = aux.from_kwargs("label_subplots", True, kwargs, rm = True)
        start_character = aux.from_kwargs("labeltype", "A", kwargs, rm = True)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        show = aux.from_kwargs("show", True, kwargs, rm = True)
        show_violins = aux.from_kwargs("violin", True, kwargs, rm = True)
        rot = aux.from_kwargs("rot", None, kwargs, rm = True)
        alpha = aux.from_kwargs("alpha", 1, kwargs, rm = True)
        xlabel = aux.from_kwargs("xlabel", None, kwargs, rm = True)
        ylabel = aux.from_kwargs("ylabel", None, kwargs, rm = True)

        # generate a custom color palette in case color kwarg is provided
        palette = gx.generate_palette(kwargs)

        # set a seaborn style
        style = aux.from_kwargs("style", "dark", kwargs, rm = True)
        sns.set_style( style )

        # make figure
        figsize = aux.from_kwargs("figsize", None, kwargs, rm = True)
        fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize)
        fig.suptitle( title )

        # setup Coords
        Coords = gx.AxesCoords(fig, axs, (nrows, ncols))
        idx = 0
        for assay in assays: 

            try: 
                tmp_df = data[ ["group", "group_name", assay] ]
                tmp_df = tmp_df.sort_values("group")

                # now plot a new violin chart 
                subplot = Coords.subplot()

                if show_violins:
                    sns.violinplot(
                                    x = tmp_df["group_name"],
                                    y = tmp_df[ assay ],
                                    color = None,
                                    inner = None, 
                                    palette = palette,
                                    ax = subplot,

                                )
                    for i in subplot.collections:
                        i.set_alpha( alpha * 0.3 )

                sns.stripplot(
                                x = tmp_df["group_name"],
                                y = tmp_df[ assay ],
                                palette = palette,
                                alpha = alpha,
                                ax = subplot,
                                **kwargs
                            )
            
                subplot.set(
                            title = assay if headers is None else headers[idx],
                            xlabel = xlabel,
                            ylabel = ylabel,        
                        )

                # adjust xtick rotation   
                if rot is not None: 
                    align = "center" if rot == 0 else "left"
                    plt.setp( subplot.xaxis.get_majorticklabels(), rotation = -rot, ha=align, rotation_mode="anchor") 

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
        data = self._rep_data

        # get assays to plot
        assays = [ i for i in data.columns if i not in [ "group", "group_name", raw_col_names[0] ] ]
        
        # get the info of groups present in the data
        # for later use as xtick labels
        group_names = data[ ["group", "group_name"] ]
        group_names = group_names.sort_values( "group" )
        group_names = aux.sorted_set( group_names[ "group_name"] )
        ticks = np.arange( len(group_names) )

        # make subplot layout
        ncols, nrows = gx.make_layout_from_list( assays )

        try: 
            # get incompaltible kwargs 
            headers = aux.from_kwargs("headers", assays, kwargs, rm = True)
            show = aux.from_kwargs("show", True, kwargs, rm = True)
            show_violins = aux.from_kwargs("violin", True, kwargs, rm = True)

            # setup padding
            padding = aux.from_kwargs("padding", 0.1, kwargs, rm = True)
            if isinstance(padding, (list, tuple)):
                hpad, vpad = padding
            else: 
                hpad, vpad = padding, None

            # make figure
            speclist = gx.make_speclist(nrows, ncols, "xy")
            fig = make_subplots(
                                    rows = nrows, cols = ncols, 
                                    specs = speclist, 
                                    subplot_titles = headers,
                                    horizontal_spacing = hpad,
                                    vertical_spacing = vpad
                            )

            # setup subplot coords handling
            Coords = gx.AxesCoords(fig, [], (ncols, nrows))

            # setup figure
            fig.update_layout(
                                title = aux.from_kwargs("title", "Results Preview", kwargs, rm = True), 
                                height = aux.from_kwargs("height", None, kwargs, rm = True), 
                                width = aux.from_kwargs("width", None, kwargs, rm = True),
                                margin = {"autoexpand": True, "pad" : 0, "b": 1, "t": 50}, 
                                autosize = True, 
                                template = aux.from_kwargs("template", "plotly_white", kwargs, rm = True),
                                # legend = {  "title" : aux.from_kwargs("legend_title", "Assay", kwargs, rm = True),  },
                                showlegend = False
                            )

            # set default axes labels
            xlabel = aux.from_kwargs("xlabel", "", kwargs, rm = True) 
            fig.update_xaxes(title_text = xlabel)
            
            ylabel = aux.from_kwargs("ylabel", "", kwargs, rm = True) 
            fig.update_yaxes(title_text = ylabel)

            idx = 0
            for assay, header in zip( assays, headers ):
                
                row, col = Coords.get()

                tmp_df = data[ ["group", "group_name", assay] ]
                tmp_df = tmp_df.sort_values("group")

                # now plot a new violin chart 
                fig.add_trace(

                    go.Violin(
                                    name = header,
                                    y = tmp_df[ assay ], x = tmp_df[ "group" ], 
                                    points = "all",
                                    pointpos = 0,
                                    hoverinfo = aux.from_kwargs("hoverinfo", "y", kwargs, rm = True), 
                                    **kwargs
                            ), 
                                    row, col
                        )

                # remove violins if not desired (leaving only dot plot)
                if not show_violins: 
                    fig.update_traces(
                                        fillcolor="rgba(0,0,0,0)", 
                                        line_width = 0, 
                                        selector=dict(type='violin'),
                                    )


                # update x axis to categorical group names
                fig.update_layout(
                                    { 
                                        f"xaxis{idx+1}" : dict(
                                                                tickmode = 'array',
                                                                tickvals = ticks,
                                                                ticktext = group_names,
                                                        )         
                                    }
                                )

                Coords.increment()

                idx += 1

            if show:
                fig.show()

            # print("<<< Interactive Check >>>")
            return fig 
        except Exception as e: 
            raise e 


class GroupBars(Plotter):
    """
    Generates a Bar plot figure with a separate subplot for each group. 

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).



    Plotting Kwargs
    ================
    
    `"static"` Kwargs

	    Static GroupBars figures accept the following kwargs:
    
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                  | Description                                                                                     | Example                                       |
        +===========================+=================================================================================================+===============================================+
        | show : `bool`             | Whether or not to show the figure                                                               | `show = True` (default)                       |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | figsize : `tuple`         | The figure size                                                                                 | `figsize = (10, 4)`                           |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`             | The overall figure title                                                                        | `title = "Today's results"                    |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`            | The x-axis label of each subplot                                                                | `xlabel = "Conditions"`                       |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`            | The y-axis label of each subplot                                                                | `ylabel = "Mean Mean Fold Change"`            |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | rot : `float`             | The rotation of x-axis labels                                                                   | `rot = 0.3`                                   |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`          | A list of titles for each subplot in the preview figure                                         | `headers = ["transcript A", "transcript B"]`  |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | label_subplots  : `bool`  | Add each subplot with A, B, C ... (if True, default)                                            | `label_subplots = True` (default)             |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | labeltype : `str`         | The starting character for subplot labelling. By default an `"A"`.                              | `labeltype = "a"`                             |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | frame   : `bool`          | Show left and top spines of subplots (if True)                                                  | `frame = False` (default)                     |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | color : `str or list`     | The fillcolor for the individual bars                                                           | `color = "yellow"`                            |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | style : `str`             | A `seaborn` style to set. Check out available styles [1].                                       | `style = "darkgrid"`                          |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgecolor : `str or list` | The edgecolor for the individual bars                                                           | `edgecolor = "black"`                         |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgewidth : `float`       | The width of the edge of individual bars                                                        | `edgewidth = 0.5`                             |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ecolor : `str`            | The color of errorbars                                                                          | `ecolor = "orange"`                           |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | **kwargs                  | Any additional kwargs that can be passed to the `matplotlib`-backend pandas `.plot.bar()` API.  |                                               |
        +---------------------------+-------------------------------------------------------------------------------------------------+-----------------------------------------------+


    
    `"interactive"` Kwargs

	    Interactive GroupBars figures accept the following kwargs:

        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                    | Example                                       |
        +=============================+================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                              | `show = True` (default)                       |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                       | `title = "Today's results"`                   |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`              | The x axis label                                                                                                                                               | `xlabel = "My super qPCR samples"`            |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                               | `ylabel = "Mean of ddCt"`                     |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                           | `height = 50`                                 |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                            | `width = 50`                                  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                       | `padding = 0.2`                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                               | `template = "plotly_dark"`                    |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                        | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | legend_title : `str`        | The title to be displayed above the legend                                                                                                                     | `legend_title = "my assays"`                  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default just `"y"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                  | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Bar()`.                                                                                     |                                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+


    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """
    def __init__(self, mode : str = None ):
        self._setup_default_params(
                                    static = defaults.static_PreviewBars, 
                                    interactive = defaults.interactive_PreviewBars
                                )
        super().__init__( mode = mode )
    
    
    def _static_plot(self, **kwargs):
        """
        Generate a static figure
        """

        kwargs = self.update_params(kwargs)
        data = self.get()
        groups = data["group"].unique()
        names = data["group_name"].unique()

        ref_col, ncols, nrows, headings, x, y, sterr, query = self._prep_properties(kwargs)
        x = "assay" 

        title = aux.from_kwargs("title", "Results Preview", kwargs, rm = True)
        headers = aux.from_kwargs("headers", names, kwargs, rm = True)
        label_subplots = aux.from_kwargs("label_subplots", True, kwargs, rm = True)
        start_character = aux.from_kwargs("labeltype", "A", kwargs, rm = True)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        show = aux.from_kwargs("show", True, kwargs, rm = True)
        rot = aux.from_kwargs("rot", None, kwargs, rm = True)

        palette = gx.generate_palette(kwargs)

        xlabel = aux.from_kwargs("xlabel", None, kwargs, rm = True)
        ylabel = aux.from_kwargs("ylabel", "$\Delta\Delta$Ct", kwargs, rm = True)
        # set a seaborn style
        style = aux.from_kwargs("style", "ticks", kwargs, rm = True)
        sns.set_style( style )

        edgecolor = aux.from_kwargs("edgecolor", "white", kwargs, rm = True)
        edgewidth = aux.from_kwargs("linewidth", 1, kwargs, rm = True)
        edgewidth = aux.from_kwargs("edgewidth", edgewidth, kwargs, rm = True)
        figsize = aux.from_kwargs("figsize", None, kwargs, rm = True)

        groups = data["group"].unique()
        nrows, ncols = gx.make_layout_from_list( groups )
        
        fig, axs = plt.subplots( ncols, nrows, figsize = figsize )
        fig.suptitle( title )

        coords = gx.AxesCoords( fig, axs, (ncols, nrows) )

        idx = 0 
        for group, name in zip(groups, headers): 

            ax = coords.subplot()
            
            tmp_df = data.query( f"group == {group}" )
            tmp_df = tmp_df.sort_values( x )

            tmp_df.plot.bar(
                        x = x, 
                        y = y,
                        # yerr = "stdev",
                        color = palette,
                        edgecolor = edgecolor,
                        linewidth = edgewidth,
                        ax = ax,
                        **kwargs
                    )

            ax.errorbar(
                            x = tmp_df[x], y = tmp_df[y], 
                            yerr = tmp_df[sterr], 
                            fmt = ".", markersize = 0, capsize = 3, 
                            ecolor = aux.from_kwargs("ecolor", "black", kwargs),
                    )
            
            ax.set(
                            title = name,
                            xlabel = xlabel,
                            ylabel = ylabel,        
                )

            if rot is not None: 
                align = "center" if rot == 0 else "left"
                plt.setp( ax.xaxis.get_majorticklabels(), rotation = -rot, ha=align, rotation_mode="anchor") 


            if not show_spines:
                ax.spines["right"].set_visible(False)
                ax.spines["top"].set_visible(False)
                ax.spines["left"].set_linewidth(1.05)
                ax.spines["bottom"].set_linewidth(1.05)

                # add ABCD... label to subplot
                if label_subplots:
                    self._add_subplot_label(idx, ax, start_character)         

            coords.increment()
            idx += 1

        plt.tight_layout()

        if show:
            plt.show()
        
        return fig 


    def _interactive_plot(self, **kwargs):
        """
        The interactive Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self.get()
        groups = data["group"].unique()
        names = data["group_name"].unique()

        try: 
            # setup figure framework variables that have to be removed from 
            ref_col, ncols, nrows, headings, x, y, sterr, query = self._prep_properties(kwargs)
            x = "assay"
            headers = aux.from_kwargs("headers", names, kwargs, rm = True)
            show = aux.from_kwargs("show", True, kwargs, rm = True)
            padding = aux.from_kwargs("padding", 0.1, kwargs, rm = True)

            if isinstance(padding, (list, tuple)):
                hpad, vpad = padding
            else: 
                hpad, vpad = padding, None

            # make new ncols and nrows based on the groups list
            nrows, ncols = gx.make_layout_from_list(groups)

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
                        legend = {"title" : aux.from_kwargs("legend_title", "group name", kwargs, rm = True),
                        },
                    )

            # set default axes labels
            xlabel = aux.from_kwargs("xlabel", "", kwargs, rm = True) 
            fig.update_xaxes(title_text = xlabel)
            
            ylabel = aux.from_kwargs("ylabel", "DeltaDeltaCt", kwargs, rm = True) 
            fig.update_yaxes(title_text = ylabel)

            idx = 0
            for group, name in zip( groups, names ):
                row, col = Coords.get()
                tmp_df = data.query(f"group == {group}")
                tmp_df = tmp_df.sort_values( x )

                # now plot a new bar chart 
                fig.add_trace(

                    go.Bar(
                        name = name,
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





class GroupDots(Plotter):
    """
    Generate a Preview of all results from all Assays in subplots, plotting individual values
    in a Dot Plot rather than Bar Plot.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).


    Plotting Kwargs
    ================
    
    `"static"` Kwargs

	    Static GroupDots figures accept the following kwargs:

        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                  | Description                                                                 | Example                                       |
        +===========================+=============================================================================+===============================================+
        | show : `bool`             | Whether or not to show the figure                                           | `show = True` (default)                       |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | figsize : `tuple`         | The figure size                                                             | `figsize = (10, 4)`                           |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`             | The overall figure title                                                    | `title = "Today's results"                    |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`            | The x-axis label of each subplot                                            | `xlabel = "Conditions"`                       |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`            | The y-axis label of each subplot                                            | `ylabel = "Mean Fold Change"`                 |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | rot : `float`             | The rotation of x-axis labels                                               | `rot = 0.3`                                   |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`          | A list of titles for each subplot in the preview figure                     | `headers = ["transcript A", "transcript B"]`  |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | label_subplots  : `bool`  | Add each subplot with A, B, C ... (if True, default)                        | `label_subplots = True` (default)             |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | labeltype : `str`         | The starting character for subplot labelling. By default an `"A"`.          | `labeltype = "a"`                             |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | frame   : `bool`          | Show left and top spines of subplots (if True)                              | `frame = False` (default)                     |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | violin   : `bool`         | Show symmetric kde of the dots of each group (if True).                     | `violin = False`                              |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | color : `str or list`     | The fillcolor for the individual dots from replicate groups                 | `color = "yellow"`                            |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | style : `str`             | A `seaborn` style to set. Check out available styles [1].                   | `style = "darkgrid"`                          |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                | Any additional kwargs that can be passed to the `seaborn`'s `stripplot()`.  |                                               |
        +---------------------------+-----------------------------------------------------------------------------+-----------------------------------------------+

    
    `"interactive"` Kwargs

	    Interactive GroupDots figures accept the following kwargs:

        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                    | Example                                       |
        +=============================+================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                              | `show = True` (default)                       |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                       | `title = "Today's results"`                   |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`              | The x axis label                                                                                                                                               | `xlabel = "My super qPCR samples"`            |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                               | `ylabel = "Mean of ddCt"`                     |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                           | `height = 50`                                 |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                            | `width = 50`                                  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                       | `padding = 0.2`                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                               | `template = "plotly_dark"`                    |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                        | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | violin   : `bool`           | Show symmetric kde of the dots of each group (if True).                                                                                                        | `violin = False`                              |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default just `"y"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\kwargs                   | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Violin()`.                                                                                  |                                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+

    
    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """
    def __init__(self, mode:str = None):
        self._setup_default_params(
                                    static = defaults.static_PreviewDots, 
                                    interactive = defaults.interactive_PreviewDots
                                )
        # __init__ is going to require default_params to be already set!
        super().__init__(mode)

    def _static_plot(self, **kwargs):
        """
        The static Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self._rep_data

        # get assays to plot
        assays = [ i for i in data.columns if i not in [ "group", "group_name", raw_col_names[0] ] ]
        assays = sorted( assays )

        # get the groups
        groups = data["group"].unique()
        names = data["group_name"].unique()
        
        # generate subplot layout
        ncols, nrows = gx.make_layout_from_list( groups )

        # get kwargs incompatible with the main plotting method
        title = aux.from_kwargs("title", None, kwargs, rm = True)
        headers = aux.from_kwargs("headers", names, kwargs, rm = True)
        label_subplots = aux.from_kwargs("label_subplots", True, kwargs, rm = True)
        start_character = aux.from_kwargs("labeltype", "A", kwargs, rm = True)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        show = aux.from_kwargs("show", True, kwargs, rm = True)
        show_violins = aux.from_kwargs("violin", True, kwargs, rm = True)
        rot = aux.from_kwargs("rot", None, kwargs, rm = True)
        alpha = aux.from_kwargs("alpha", 1, kwargs, rm = True)
        xlabel = aux.from_kwargs("xlabel", None, kwargs, rm = True)
        ylabel = aux.from_kwargs("ylabel", None, kwargs, rm = True)

        # generate a custom color palette in case color kwarg is provided
        palette = gx.generate_palette(kwargs)

        # set a seaborn style
        style = aux.from_kwargs("style", "dark", kwargs, rm = True)
        sns.set_style( style )

        # make figure
        figsize = aux.from_kwargs("figsize", None, kwargs, rm = True)
        fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize)
        fig.suptitle( title )

        # setup Coords
        Coords = gx.AxesCoords(fig, axs, (nrows, ncols))
        idx = 0
        for group, name in zip( groups, headers ): 

            try:
                
                # prepare a vertical bigtable format dataframe 
                tmp_df = self._prepare_df(data, assays, group)
                
                # now plot a new violin chart 
                subplot = Coords.subplot()

                if show_violins:
                    sns.violinplot(
                                    x = tmp_df[ "assay" ],
                                    y = tmp_df[ "value" ],
                                    color = None,
                                    inner = None, 
                                    palette = palette,
                                    ax = subplot,

                                )
                    for i in subplot.collections:
                        i.set_alpha( alpha * 0.3 )

                sns.stripplot(
                                x = tmp_df[ "assay" ],
                                y = tmp_df[ "value" ],
                                palette = palette,
                                alpha = alpha,
                                ax = subplot,
                                **kwargs
                            )
            
                subplot.set(
                            title = name,
                            xlabel = xlabel,
                            ylabel = ylabel,        
                        )

                # adjust xtick rotation   
                if rot is not None: 
                    align = "center" if rot == 0 else "left"
                    plt.setp( subplot.xaxis.get_majorticklabels(), rotation = -rot, ha=align, rotation_mode="anchor") 

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

    def _prepare_df(self, data, assays, group):
        """
        Concatenates the different assay colums into a single 
        set of two columns, one for all "value"s and one for 
        the "assay" identifiers. 
        """
        # get a group subset 
        group_df = data.query( f"group == {group}" )

        # concatenate all separate assay columns into a set of 
        # 'value' and 'assay' columns.

        # setup the _assays concatenated dataframe   
        _assays = group_df[ [ "group", assays[0] ] ]
        _assays = _assays.rename( columns = { assays[0] : "value" } )
        _assays["assay"] = [ assays[0] for i in _assays["group"] ]

        # now iteratively add all remaining assays
        for assay in assays[1:]:
            tmp_df = group_df[ [ "group", assay ] ]
            tmp_df = tmp_df.rename( columns = { assay : "value" } )
            tmp_df["assay"] = [ assay for i in tmp_df["group"] ]
            _assays = pd.concat( [_assays, tmp_df], ignore_index = True)
                
        # remove the group col, and sort
        tmp_df = _assays.drop( columns = ["group"] )
        tmp_df = tmp_df.sort_values( "assay" )

        return tmp_df

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
        data = self._rep_data

        # get assays to plot
        assays = [ i for i in data.columns if i not in [ "group", "group_name", raw_col_names[0] ] ]
        assays = sorted( assays )

        # get the groups
        groups = data["group"].unique()
        names = data["group_name"].unique()
        ticks = np.arange( len(assays) )

        # make subplot layout
        ncols, nrows = gx.make_layout_from_list( groups )

        try: 
            # get incompaltible kwargs 
            headers = aux.from_kwargs("headers", names, kwargs, rm = True)
            show = aux.from_kwargs("show", True, kwargs, rm = True)
            show_violins = aux.from_kwargs("violin", True, kwargs, rm = True)

            # setup padding
            padding = aux.from_kwargs("padding", 0.1, kwargs, rm = True)
            if isinstance(padding, (list, tuple)):
                hpad, vpad = padding
            else: 
                hpad, vpad = padding, None

            # make figure
            speclist = gx.make_speclist(nrows, ncols, "xy")
            fig = make_subplots(
                                    rows = nrows, cols = ncols, 
                                    specs = speclist, 
                                    subplot_titles = headers,
                                    horizontal_spacing = hpad,
                                    vertical_spacing = vpad
                            )

            # setup subplot coords handling
            Coords = gx.AxesCoords(fig, [], (ncols, nrows))

            # setup figure
            fig.update_layout(
                                title = aux.from_kwargs("title", "Results Preview", kwargs, rm = True), 
                                height = aux.from_kwargs("height", None, kwargs, rm = True), 
                                width = aux.from_kwargs("width", None, kwargs, rm = True),
                                margin = {"autoexpand": True, "pad" : 0, "b": 1, "t": 50}, 
                                autosize = True, 
                                template = aux.from_kwargs("template", "plotly_white", kwargs, rm = True),
                                # legend = {  "title" : aux.from_kwargs("legend_title", "Assay", kwargs, rm = True),  },
                                showlegend = False
                            )

            # set default axes labels
            xlabel = aux.from_kwargs("xlabel", "", kwargs, rm = True) 
            fig.update_xaxes(title_text = xlabel)
            
            ylabel = aux.from_kwargs("ylabel", "", kwargs, rm = True) 
            fig.update_yaxes(title_text = ylabel)

            idx = 0
            for group, name in zip( groups, headers ):
                
                # prepare a vertical bigtable format dataframe 
                tmp_df = self._prepare_df(data, assays, group)

                row, col = Coords.get()

                # now plot a new violin chart 
                fig.add_trace(

                    go.Violin(
                                    name = name,
                                    x = tmp_df[ "assay" ], 
                                    y = tmp_df[ "value" ], 
                                    points = "all",
                                    pointpos = 0,
                                    hoverinfo = aux.from_kwargs("hoverinfo", "y", kwargs, rm = True), 
                                    **kwargs
                            ), 
                                    row, col
                        )

                # remove violins if not desired (leaving only dot plot)
                if not show_violins: 
                    fig.update_traces(
                                        fillcolor="rgba(0,0,0,0)", 
                                        line_width = 0, 
                                        selector=dict(type='violin'),
                                    )


                # update x axis to categorical group names
                fig.update_layout(
                                    { 
                                        f"xaxis{idx+1}" : dict(
                                                                tickmode = 'array',
                                                                tickvals = ticks,
                                                                ticktext = assays,
                                                        )         
                                    }
                                )

                Coords.increment()

                idx += 1

            if show:
                fig.show()

            # print("<<< Interactive Check >>>")
            return fig 
        except Exception as e: 
            raise e 


class EfficiencyLines(Plotter):
    """
    Generates a Figure for the linear regressions used for Assay efficiency 
    calculations. This FigureClass specifically works with the `qpcr.Calibrator` 
    class. 

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).
    
    

    Plotting Kwargs
    ================

    `"static"` Kwargs

	    Static EfficiencyLines figures accept the following kwargs:

        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                  | Description                                                                                    | Example                                       |
        +===========================+================================================================================================+===============================================+
        | show : `bool`             | Whether or not to show the figure                                                              | `show = True` (default)                       |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | figsize : `tuple`         | The figure size                                                                                | `figsize = (10, 4)`                           |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`             | The overall figure title                                                                       | `title = "Today's efficiencies"               |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`            | The x-axis label of each subplot                                                               | `xlabel = "Log Dilution"`                     |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`            | The y-axis label of each subplot                                                               | `ylabel = "My Ct Values"`                     |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | rot : `float`             | The rotation of x-axis labels                                                                  | `rot = 0.3`                                   |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`          | A list of titles for each subplot in the preview figure                                        | `headers = ["transcript A", "transcript B"]`  |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | label_subplots  : `bool`  | Add each subplot with A, B, C ... (if True, default)                                           | `label_subplots = True` (default)             |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | labeltype : `str`         | The starting character for subplot labelling. By default an `"A"`.                             | `labeltype = "a"`                             |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | frame   : `bool`          | Show left and top spines of subplots (if True)                                                 | `frame = False` (default)                     |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | color : `str or list`     | The fillcolor for the individual dots                                                          | `color = "yellow"`                            |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | style : `str`             | A `seaborn` style to set. Check out available styles [1].                                      | `style = "darkgrid"`                          |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgecolor : `str or list` | The edgecolor for the individual dots for the datapoints.                                      | `edgecolor = "black"`                         |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgewidth : `float`       | The width of the edge of individual dots.                                                      | `edgewidth = 0.5`                             |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | linecolor : `str or list` | The color for regression line.                                                                 | `linecolor = "crimson"`                       |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | linewidth : `float`       | The width of the regression line.                                                              | `edgewidth = 0.5`                             |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\kwargs                 | Any additional kwargs that can be passed to `seaborn`'s `scatterplot` and `lineplot` (both!).  |                                               |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+


    
    `"interactive"` Kwargs

	    Interactive EfficiencyLines figures accept the following kwargs:

        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                      | Example                                       |
        +=============================+==================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                                | `show = True` (default)                       |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                         | `title = "Today's efficiencies"`              |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`              | The x axis label                                                                                                                                                 | `xlabel = "Log Dilution"`                     |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                                 | `ylabel = My super Ct values"`                |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                             | `height = 50`                                 |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                              | `width = 50`                                  |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                         | `padding = 0.2`                               |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                                 | `template = "plotly_dark"`                    |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                          | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default just `"y+x"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\kwargs                   | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Scatter()`.                                                                                   |                                               |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+


    
    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """
    def __init__(self, mode : str = None ):
        self._setup_default_params(
                                    static = defaults.static_EfficiencyLines, 
                                    interactive = defaults.interactive_EfficiencyLines
                                )
        super().__init__( mode = mode )
        self._Calibrator = None

    def link( self, calibrator : main.Calibrator ):
        """
        Links a `qpcr.Calibrator` object to source data from.

        Note
        -----
        Only a `qpcr.Calibrator` that has actually *de novo*
        computed efficiencies will have data to plot, and only for 
        the newly computed effiencies! 

        Parameters
        ----------
        calibrator : qpcr.Calibrator
            A `qpcr.Calibrator` object that has computed new efficiencies.
        """
        self._Calibrator = calibrator
    
    def _static_plot(self, **kwargs):
        """
        Generates a static Effiency Regression Figure
        """
        kwargs = self.update_params(kwargs)
        data = self._Calibrator._computed_values
        headers = list(  data.keys()  )

        title = aux.from_kwargs("title", None, kwargs, rm = True)
        headers = aux.from_kwargs("headers", headers, kwargs, rm = True)
        label_subplots = aux.from_kwargs("label_subplots", True, kwargs, rm = True)
        start_character = aux.from_kwargs("labeltype", "A", kwargs, rm = True)
        show_spines = aux.from_kwargs("frame", True, kwargs, rm = True)
        show = aux.from_kwargs("show", True, kwargs, rm = True)
        rot = aux.from_kwargs("rot", None, kwargs, rm = True)
        
        palette, is_palette = gx.generate_palette(kwargs, True )

        xlabel = aux.from_kwargs("xlabel", None, kwargs, rm = True)
        ylabel = aux.from_kwargs("ylabel", None, kwargs, rm = True)

        # set a seaborn style
        style = aux.from_kwargs("style", "ticks", kwargs, rm = True)
        sns.set_style( style )

        edgecolor = aux.from_kwargs("edgecolor", "white", kwargs, rm = True)
        edgewidth = aux.from_kwargs("edgewidth", None, kwargs, rm = True)

        linecolor = aux.from_kwargs("linecolor", "black", kwargs, rm = True)
        linewidth = aux.from_kwargs("linewidth", 1, kwargs, rm = True)

        ncols, nrows = gx.make_layout_from_list( headers )

        fig, axs = plt.subplots( nrows = nrows, ncols = ncols)
        fig.suptitle( title )
        Coords = gx.AxesCoords( fig, axs, subplots = (nrows, ncols) )

        idx = 0 
        for id, obj in data.items(): 
            
            ax = Coords.subplot()

            # get data to plot              
            dilutions, cts = obj.values()
            eff = obj.efficiency()
            model = obj.model()
            Rsquare = model.rvalue ** 2

            # plot the regression line
            yvals = model.slope * dilutions + model.intercept
            sns.lineplot(
                                x = dilutions,
                                y = yvals,
                                color = linecolor,
                                linewidth = linewidth,
                                ax = ax,
                                **kwargs
                    )

            # plot the original data
            sns.scatterplot(
                                x = dilutions,
                                y = cts,

                                # this is somewhat hacky but it allows us to use 
                                # either a colormap or a single color
                                hue = dilutions if is_palette else None,
                                color = None if is_palette else palette,
                                palette = palette if is_palette else None,

                                edgecolor = edgecolor,
                                linewidth = edgewidth,
                                ax = ax,
                                **kwargs
                            )

            # add additional info to a legend
            ax.legend( handles = [
                                    # R^2 Value
                                    Line2D(  [0], [0], 
                                            color = "black", 
                                            visible = False, 
                                            label = f"$R^2$ \t = {Rsquare:.4f}"
                                        ),
                                    # Efficiency Value
                                    Line2D(  [0], [0], 
                                            color = "black", 
                                            visible = False, 
                                            label = f"$eff.$\t= {eff:.4f}"
                                        )
                                ], 
                                loc = "upper right",
                                frameon = False 
                    )

            # some formatting...
            ax.set(
                    title = id, 
                    xlabel = xlabel,
                    ylabel = ylabel,
                )

            if rot is not None: 
                align = "center" if rot == 0 else "left"
                plt.setp( ax.xaxis.get_majorticklabels(), rotation = -rot, ha=align, rotation_mode="anchor") 

            # add ABCD... label to subplot
            if label_subplots:
                self._add_subplot_label(idx, ax, start_character)         

            if not show_spines:
                sns.despine()
            
            idx += 1
            Coords.increment()

        plt.tight_layout()
        if show: 
            plt.show()

        return fig 

    def _interactive_plot(self, **kwargs):
        """
        Generates an interactive EfficiencyLines figure
        """
        kwargs = self.update_params(kwargs)
        data = self._Calibrator._computed_values
        headers = list(  data.keys()  )

        title = aux.from_kwargs("title", None, kwargs, rm = True)
        headers = aux.from_kwargs("headers", headers, kwargs, rm = True)
        show = aux.from_kwargs("show", True, kwargs, rm = True)
        padding = aux.from_kwargs("padding", 0.1, kwargs, rm = True)

        if isinstance(padding, (list, tuple)):
            hpad, vpad = padding
        else: 
            hpad, vpad = padding, None

        # make new ncols and nrows based on the headers list
        nrows, ncols = gx.make_layout_from_list(headers)

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
                    title = title, 
                    height = aux.from_kwargs("height", None, kwargs, rm = True), 
                    width = aux.from_kwargs("width", None, kwargs, rm = True),
                    margin = {"autoexpand": True, "pad" : 0, "b": 1, "t": 50}, 
                    autosize = True, 
                    template = aux.from_kwargs("template", "plotly_white", kwargs, rm = True),
                    # legend = {"title" : aux.from_kwargs("legend_title", "Assay", kwargs, rm = True), },
                )

        # set default axes labels
        xlabel = aux.from_kwargs("xlabel", None, kwargs, rm = True) 
        fig.update_xaxes(title_text = xlabel)
        
        ylabel = aux.from_kwargs("ylabel", None, kwargs, rm = True) 
        fig.update_yaxes(title_text = ylabel)

        Coords = gx.AxesCoords(fig, [], (ncols, nrows))

        # iterate over the data...
        idx = 0 
        for id, obj in data.items():

            row, col = Coords.get()

            # get data to plot              
            dilutions, cts = obj.values()
            eff = obj.efficiency()
            model = obj.model()
            Rsquare = model.rvalue ** 2

            # plot the original data
            fig.add_trace(
                            go.Scatter(
                                name = id,
                                x = dilutions,
                                y = cts,  
                                mode = "markers",
                                showlegend = False,
                                hoverinfo = aux.from_kwargs("hoverinfo", "y+x", kwargs, rm = True), 
                                **kwargs
                            ), 
                            row, col
                        )
            
            # plot the regression line
            yvals = model.slope * dilutions + model.intercept

            fig.add_trace(
                            go.Scatter(
                                name = id,
                                x = dilutions,
                                y = yvals,  
                                mode = "lines",
                                showlegend = False,
                                **kwargs
                            ), 
                            row, col
                        )


            # add infos as annotations (because custom legends don't work here)
            # NOTE: The xref x{idx} is a pretty nice lifehack to essentially 
            #       emulate the ax behaviour from matplotib.subplots... 
            xref = "x" if idx == 0 else f"x{idx+1}"
            yref = "y" if idx == 0 else f"y{idx+1}"
            info_text = f"R^2 = {Rsquare:.4f}<br>eff. = {eff:.4f}"

            fig.add_annotation(
                            dict(
                                    x = 0, 
                                    y = np.max( cts ), 
                                    xref = xref, 
                                    yref = yref, 
                                    text = info_text, 
                                    showarrow = False
                                )
                            )

            Coords.increment()
            idx += 1

        if show: 
            fig.show()

        return fig 


def plot( obj, mode = None, **kwargs ):
    """
    A generic plotting shortcut to visualise the data from a `qpcr` class object, if it supports visualisation.

    Parameters
    ----------
    obj 
        The object to visualise from.
    mode : str
        The plotting mode to use. This can be either "static" (matplotlib) or "interactive" (plotly)
    **kwargs
        Any additional keyword arguments.

    Returns
    -------
    fig
        The figure produced.
    """
    kwargs["mode"] = mode
    func = obj.__qplot__( **kwargs )
    logger.debug( kwargs )
    fig = func( **kwargs )
    return fig

def interactive():
    """Set the default plotting mode to ``interactive``"""
    defaults.plotmode = "interactive"
    return defaults.plotmode

def static():
    """Set the default plotting mode to ``static``"""
    defaults.plotmode = "static"
    return defaults.plotmode

# if __name__ == '__main__':

    # files = ["Example Data/28S.csv", "Example Data/28S_again.csv", "Example Data/actin.csv",
    # "Example Data/actin2.csv",
    # "Example Data/actin3.csv",
    # "Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]
    # groupnames = ["wt-", "wt+", "ko-", "ko+"]

    # # setup figures
    # a = PreviewResults(mode = "static")
    # b = PreviewResults(mode = "interactive")

    # a.params(
    #       frame = True, labeltype = "a", show = False
    # )
    # b.params(
    #     show = False, template = "plotly"  
    # )

    # # predefined pipeline use
    # pipe = qpcr.Pipes.Basic()
    # pipe.link(files)
    # pipe.add_normalisers(files[:2])
    # pipe.replicates(6)
    # pipe.names(groupnames)


    # print()
    # print("==== Multiple (ALL) Samples ====")
    # pipe.run()
    # result = pipe.get(kind = "obj")
    # try:
    #     a.link(result)
    #     a.plot(color = "green")
    # except: 
    #     print("Static Failed!")

    # try:
    #     b.link(result)
    #     b.plot()
    # except: print("Interactive Failed!")

    # print()
    # print("==== Multiple (7) Samples ====")
    # pipe.link(files[:7])
    # pipe.run()
    # result = pipe.get(kind = "obj")
    # try:
    #     a.link(result)
    #     a.plot(show = False)
    # except: 
    #     print("Static Failed!")

    # try:
    #     b.link(result)
    #     b.plot(template = "plotly")
    # except: print("Interactive Failed!")

    # print()
    # print("==== Multiple (6) Samples ====")
    # pipe.link(files[:6])
    # pipe.run()
    # result = pipe.get(kind = "obj")
    # try:
    #     a.link(result)
    #     a.plot()
    # except: 
    #     print("Static Failed!")

    # try:
    #     b.link(result)
    #     b.plot()
    # except: print("Interactive Failed!")

    # print()
    # print("==== Multiple (4) Samples ====")
    # pipe.link(files[:4])
    # pipe.run()
    # result = pipe.get(kind = "obj")
    # try:
    #     a.link(result)
    #     a.plot()
    # except: 
    #     print("Static Failed!")

    # try:
    #     b.link(result)
    #     b.plot()
    # except: print("Interactive Failed!")

    # print()
    # print("==== Multiple (3) Samples ====")
    # pipe.link(files[:3])
    # pipe.run()
    # result = pipe.get(kind = "obj")
    # try:
    #     a.link(result)
    #     a.plot()
    # except Exception as e: 
    #     print(e)
    #     print("Static Failed!")

    # try:
    #     b.link(result)
    #     b.plot()
    # except Exception as e: 
    #     print("Interactive Failed!")
    #     raise e

    # print()
    # print("==== Multiple (2) Samples ====")
    # pipe.link(files[:2])
    # pipe.run()
    # result = pipe.get(kind = "obj")
    # try:
    #     a.link(result)
    #     a.plot(show = False)
    # except Exception as e: 
    #     print(e)
    #     print("Static Failed!")

    # try:
    #     b.link(result)
    #     b.plot(show = False)
    # except Exception as e: 
    #     print("Interactive Failed!")
    #     raise e


    # print()
    # print("==== One (1) Sample ====")
    # pipe.link([files[0]])
    # pipe.run()
    # result = pipe.get(kind = "obj")
    # try:
    #     a.link(result)
    #     a.plot()
    # except: 
    #     print("Static Failed!")

    # try:
    #     b.link(result)
    #     b.plot()
    # except: print("Interactive Failed!")
    # exit(0)