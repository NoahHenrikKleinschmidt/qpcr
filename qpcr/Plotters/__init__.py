"""
.. _qpcrPlotters:

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
    myboxplot.link(myassay)

    # and add some graphing paramters
    my_params = { 
                    "color" : "Reds", # we want to change the default colormap to red tones
                    "title" : "my customized figure",
                    "style" : "darkgrid", # and use the seaborn darkgrid style 
                }
    myboxplot.params(my_params)

    # and now visualise the figure
    fig = myboxplot.plot()

    # alternatively: fig = myboxplot.plot(**my_params) 


This is a rather tedious way of setting up for just a quick look at our Assay, but it may be worth it if we have to repeat this many times...

The Built-in shortuct
----------------------
For more convenience, however, the ``qpcr.Assay`` lets us directly visualize its Ct values through the ``boxplot`` method. We can call it easily like:

.. code-block:: python

    fig = myassay.boxplot(**my_params)


This will perform the setup and linking to a new instance of ``ReplicateBoxPlot``. 



The ``qpcr.plot`` function
--------------------------

As a generic way to visualize data, ``qpcr`` offers the stand-alone ``plot`` function that accepts a data-containing object and will call it's built-in plotting method.
Hence, this is just another way of calling the object's built-in plotting method. It may be convenient for you if you are used to the R plot API (and was inspired by it).
However, it is important to note that you cannot add a `Plotter` to the ``plot`` function, only data-containing classes like the ``qpcr.Assay`` or ``qpcr.Results`` support the generic ``plot`` function!
Thus, a final way of obtaining the exact same figure as before, is using:

.. code-block:: python

    fig = qpcr.plot(myassay, **my_params)


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

    fig = myresults.preview(kind = "AssayDots")


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
    mycolors = ["crimson", "black"]

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

    fig = results.preview(title = "Normalized to 28S+Actin")

However, this will now cause **all** of your Results to generate a black and red `AssayDots` figure in `ticks` style when calling their ``preview`` method (unless you specify new parameters manually).
As you will have noticed, the dictionary we edited was not a "qpcr default plotting parameters" dictionary but rather specific for the `static AssayDots` Plotters. All Plotters have two dedicated dictionaries to store their default parameters in, one for static mode, one for interactive.
You can therefore selectively adjust your favourite default settings without having to worry to skew with other figure types. However, `rcParams` will of course, also work just fine globally. 

Statistical Evaluations
-----------------------

The ``qpcr`` package offers a few statistical evaluations for the ``qpcr.Results`` class. For plotting, multiple _T-Tests_ are of interest because they can be integrated into the preview figures.
In fact, if t-tests are performed on a Results object, the significance bars are automatically integrated into the preview. 

.. code-block:: python

    # perform a t-test
    qpcr.stats.assaywise_ttests(results, pairs=[("WT-", "WT+"), ("KO-", "KO+")])

    # plot the results
    fig = results.preview()

.. image:: ../../docs/source/resources/preview_with_stats.png
    :align: center

We can further style the representation using the `pval_kws` kwarg which accepts a dictionary containing any of:

.. code-block:: 

    style : str
        The style in which to encode the p-values. Available are 
        - `"p<"`, style by pvalue levels (e.g. `p < 0.05`)
        - `"p="`, style by actual pvalues if they are below a threshold (e.g. `p = 0.000124`)
        - `"*"`, style using asterisks by 

    threshold : float 
        The significance p-value threshold.

    step : float
        The step between each significance level. I.e. a step of `0.1` means that the threshold for
        level 2 is 10 times smaller than for level 1 (e.g. from 0.05, to 0.005, to 0.0005 etc.).

    max_levels : int
        The maximal number of levels to include. This is only used for styles with levels (i.e. `"p<"` and "*"` ). 

    levels : tuple 
        A custom iterable of p-value level threshold. Only used for style `"p<"`.

    ns_default : str
        The default string to use for non-significant p-values.

    asterisk : str
        The character to use for encoding in the `"*"` style. 

    fmt : str
        A string formatter for the numeric p-values. By default these will be: 
        - `".2g"` for style `"p<"`, 
        - `".3g"` for style `"p="`,
        - not applicable for style `"*"`

"""

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
# FigureClasses:
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


from .func_api import *
from .PreviewResults import PreviewResults
from .ReplicateBoxPlot import ReplicateBoxPlot
from .FilterSummary import FilterSummary
from .EfficiencyLines import EfficiencyLines
from .AssaySubplotResults import AssayBars, AssayDots
from .GroupSubplotResults import GroupBars, GroupDots
