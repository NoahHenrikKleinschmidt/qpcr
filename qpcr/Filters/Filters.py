"""
This submodule defines a number of filters that can be used to 
remove faulty replicates from a ``qpcr.Assay`` before before passing it to a ``qpcr.Analyser``.


RangeFilter 
===========

The `RangeFilter` filters out all *raw Ct values* that do not comply to a user-specified range (default is `+/-1` around replicate group median)
The user has the option of specifying another anchor and limits for the inclusion range. The limits may be asymmetrical (e.g. 1.5 above but only 0.7 below if this is desired).
The limits are *static* and hence only the anchor is adjusted to each replicate group within each assay (unless a custom anchor is provided). 


As an example, we might wish to filter out any Ct values that are outside of *+/- 0.5* around the *group mean*. To do this we need to set up a ``RangeFilter`` and adjust both the limits and the anchor.

.. code-block:: python

    myfilter = RangeFilter()

    # set up the limits
    myfilter.set_lim( 0.5 )
    
    # set a new anchor
    def mean_anchor( values ):
        return mean( values )
    
    # or as: mean_anchor = lambda values: mean( values )

    myfilter.set_anchor( mean_anchor )

    # and now we can filter our assays

    myassay = myfilter.pipe( myassay )


By default, any Ct values that are filtered out are set to ``np.nan`` and not actually removed from the dataframes. This is because ``qpcr`` tries to retain dimensionality between the Assays.
You can choose to drop faulty outliers by using ``drop_outliers``. 


IQRFilter
=========

The `IQRFilter` filters out any outliers by `n x IQR`, where `n` is a scaling factor (default `n = 1.5`) around the replicate group median (anchor). 
Here the anchor *cannot* be re-set freely, but the limits can be adjusted symmetrically or asymmetrically as desired.

For instance, we might want to filter out any values that are `+2 IQR` but retain all lower ones (maybe because we looked at their actual values and found that they were all acceptable).

.. code-block:: python

    myfilter = IQRFilter()

    # we set some big value for retaining lower values
    myfilter.set_lim( upper = 2, lower = 100 )

    # and now we can filter our assays
    myassay = myfilter.pipe( myassay )

Filter Summary
=====================

The Filters offer a summary of their activity. This is primarily through a *summary figure* that can be called via the ``plot`` method or the ``qpcr.plot`` function.

.. code-block:: python

    fig = myfilter.plot( mode = "interactive" )


# .. raw:: html
#     :file: ../../docs/source/resources/filter_summary.html


The Filters can also export a ``txt``-file summary of the indices they filtered out for each assay. This file is designed as a human-readable quick-check and does not offer a lot of information, however.
The text summary is generated automatically if the ``report`` method is called, where you may specify a `directory` to store the filter reports to. Each filtered assay will be given a single report file.


.. code-block:: python

    # the report setup has to be done BEFORE the actual filtering!
    myfilter.report( "./filter_reports_experiment1" )

    # now we can filter
    myassay = myfilter.pipe( myassay )

Note
-----
These files are just named `filter_{assay name}.txt` so if you intend to save these reports from multiple experiments, make sure to store them in separate directories!
"""

from re import L
import qpcr
import pandas as pd
import numpy as np
import qpcr._auxiliary.warnings as aw
import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import os 
import qpcr.Plotters as Plotters
import logging

logger = aux.default_logger()
class Filter(aux._ID):
    """
    The super class of the Filters that takes in a `qpcr.Assay` object and updates its dataframe to a filtered version.
    """
    def __init__(self):
        super().__init__()
        self._Assay = None
        self._report_loc = None
        self._id = type(self).__name__
        
        # by default we ignore groups with NaN anchors
        # like this we avoid Errors and don't filter out 
        # any unicate groups like the diluent sample...
        self._ignore_nan = True

        # by default we set outliers to NaN and 
        # don't drop them anymore
        self._drop_outliers = False

        self._boxplot_mode = defaults.plotmode
        self._BoxPlotter = Plotters.FilterSummary( mode = self._boxplot_mode )

        self._filter_stats = pd.DataFrame({
                                            "assay" : [], "group" : [], 
                                            "anchor" : [], "upper" : [], "lower" : []
                                        })
    
    def plot_params(self, **params):
        """
        Allows to pre-specify plotting parameters of the FilterSummary Figure.
        This can also be passed directly while calling `Filter.plot`.

        Parameters
        ----------
        **kwargs
            Any accepted additional keyword arguments. 
        """
        self._BoxPlotter.params(**params)
            

    def get_stats(self):
        """
        Returns 
        -------
        stats : pd.DataFrame
            The filtering statistics dataframe (a summary of filtering parameters used)
        """
        return self._filter_stats

    def ignore_nan(self, bool):
        """
        Set a policy for how to deal with groups that have a `NaN` anchor.
        If set to `True` such groups will be ignored and filtering will proceed.
        If set to `False` the Filter will raise an Error!
        """
        self._ignore_nan = bool

    def drop_outliers(self, bool):
        """
        If True, will completely remove outlier Ct values from the dataset.
        If False, will set outliers to NaN
        """
        self._drop_outliers = bool

    def plotmode(self, mode ):
        """
        Set graph mode if a summary Boxplot shall be made

        Parameters
        ----------
        mode : str
            Can be either "interactive" (plotly) or "static" (matplotlib), or None to disable plotting.
        """
        if mode is None: mode = defaults.plotmode
        self._boxplot_mode = mode
        self._BoxPlotter = Plotters.FilterSummary(mode = self._boxplot_mode)
        self._BoxPlotter.params(title = "Filter Summary")

    def plot(self, **kwargs):
        """
        Generates a boxplot summary plot. 

        Note
        ----
        This is designed to be done AFTER all samples have passed the filter.

        Parameters
        ----------
        **kwargs
            Any keyword arguments that should be passed to the plotting method.
        """
        plotter = self._BoxPlotter
        fig = plotter.plot(**kwargs)
        if self._report_loc is not None and self._boxplot_mode is not None: 
            filename = f"{self.id()}_summary"
            suffix = plotter.suffix()
            plotter.save(os.path.join(self._report_loc, f"{filename}.{suffix}"))
        return fig

    def __qplot__( self, **kwargs ):
        return self.plot

    def link(self, Assay:qpcr.Assay):
        """
        Links a `qpcr.Assay` to be filtered
        
        Parameters
        ----------
            Assay : `qpcr.Assay`
                A `qpcr.Assay` object to be filtered.
        """
        self._Assay = Assay

    def pipe(self, Assay:qpcr.Assay, **kwargs):
        """
        A shortcut for link+filter.
        This is the suggested usage for Filters.
        
        Parameters
        ----------
        Assay : `qpcr.Assay`
            A `qpcr.Assay` object to be filtered.
        **kwargs
            Any keyword arguments that should be passed to the plotting method.
        
        Returns
        -------
        Assay : `qpcr.Assay`
            A `qpcr.Assay` object containing only entries that passed the filter.

        """
        if isinstance( Assay, list ):
            return [ self.pipe( assay ) for assay in Assay ]
        else:
            self.link(Assay)
            self.filter(**kwargs)
            return self._Assay    

    def filter(self, **kwargs):
        """
        Applies the filter 
        
        Parameters
        ----------
        **kwargs
            Any keyword arguments that should be passed to the filtering method.
        
        Returns
        -------
        Assay : `qpcr.Assay`
            An updated `qpcr.Assay` object containing only entries that passed the filter.
        """
        if self._Assay is not None:
            self._filter(**kwargs)
            return self._Assay
        else: 
            e = aw.FilterError( "no_assay" )
            logger.critical( e )
            raise e 

    def report(self, directory = None):
        """
        Sets up a location to store a report of any replicates that were filtered out.

        Parameters
        ----------
        directory : str
            A directory where to store report text-files and summary boxplots.

        Returns
        -------
        location : str
            If no new directory is provided, it returns the current report location.
        """
        if directory is not None:
            self._report_loc = directory

            # if the directory does not yet exist, we make it
            if not os.path.exists(self._report_loc):
                os.mkdir(self._report_loc)

        else: 
            return self._report_loc

    def reset(self):
        """
        Resets the excluded indices
        """
        self._faulty_indices = []

    def set_lim(self, lim = None, upper = None, lower = None):
        """
        Sets the range limits for the inclusion range.
        Limits can be either specified symmetrically using `lim` or asymmetrically, using `upper` and `lower`.
        
        Parameters
        ----------
        lim : float
            Sets symmetric upper and lower bounds. 
            Default settings are `lim = 1` setting both `upper` and `lower` to `1`.
        upper : float
            Sets the upper inclusion-range boundary.
        lower : float
            Sets the lower inclusion-range boundary.
        """
        if lim is not None: self._upper, self._lower = lim, lim
        if upper is not None: self._upper = upper
        if lower is not None: self._lower = lower


    def _filter(self, **kwargs):
        """
        The actual filtering function that each FilterObject will define.
        """
        print("The actual filtering function that each FilterObject will define")
        # do stuff
        return self._Assay

    def _write_report(self, faulty_indices, details={}):
        """
        Generates a filtering report file
        """
        filename = "filter_" + self._Assay.id() + ".txt"
        filename = os.path.join(self._report_loc, filename)

        report_string = f"""
Filtering Report

Filter: 
{self._id}
Assay: 
{self._Assay.id()}
Found faulty Replicates: 
{len(faulty_indices)}
Found Indices: 
{faulty_indices}
Details: 
{details}
        """
        report_string = report_string.strip()
        if os.path.exists(filename):
            with open(filename, "a") as f:
                f.write(report_string.replace("Filtering Report", ""))
        else:
            with open(filename, "w") as f:
                f.write(report_string)

    def _filter_out(self, faulty_indices):
        """
        Removes any faulty replicates based on their indices
        """
        # exclude faulty entries
        if len(faulty_indices) > 0:
            self._Assay.ignore(faulty_indices, drop = self._drop_outliers)
    
    def _save_stats(self, assay, group, anchor, upper, lower):
        """
        Saves filtering stats for a given group to self._filter_stats
        """
        new_stats = pd.DataFrame({
                                    "assay" : [assay], "group" : [group], 
                                    "anchor" : [anchor], "upper" : [upper], "lower" : [lower]
                                })
        self._filter_stats = self._filter_stats.append(new_stats, ignore_index=True)

class RangeFilter(Filter):
    """
    Filters out any replicate that lie outside a user-specified range.
    Default are `+/- 1` around the replicate-group median. 
    """
    def __init__(self):
        super().__init__()
        self._upper = 1
        self._lower = 1
        self._anchor = None

    
    def set_anchor(self, anchor):
        """
        Set the range anchor (center of inclusion range)

        Parameters
        ----------
        anchor 
            Supported types for `anchor` are: a numeric value (`int or float`),
            an `iterable` of same length as groups in the dataframe, 
            a `dict` where keys must be numeric group identifiers (starting from 0) and values are numeric values to be used as anchor (`int or float`),
            or a `function` that works with a pandas dataframe as stored by `qpcr.Assay` objects, 
            which must return a single numeric value for the anchor (it will be applied to replicate-grouped subsets of the total dataframe).
        """
        self._anchor = anchor

    def _filter(self, **kwargs):
        """
        Filters out any replicates that are out of range and updates the Assay's dataframe.
        """

        plotter = self._BoxPlotter

        plotter.add_before( self._Assay )

        df = self._Assay.get()
        groups = self._Assay.groups()

        faulty_indices = []
        for group in groups:
            tmp = df.query(f"group == {group}")

            # get anchor and check if its nan
            anchor = self._get_anchor(kwargs, group, tmp)
            if self._ignore_nan and anchor != anchor: 
                continue 

            # generate inclusion range boundries
            upper, lower = self._set_bounds(anchor)

            # get faulty indices
            Ct = defaults.raw_col_names[1]
            faulty_replicates = tmp.query(f"{Ct} < {lower} or {Ct} > {upper}")
            faulty_indices.extend(list(faulty_replicates.index))

            self._save_stats(self._Assay.id(), group, anchor, upper, lower)
        
        # remove faulty indices
        self._filter_out(faulty_indices)

        plotter.add_after( self._Assay )


        if self._report_loc is not None: 
            self._write_report(faulty_indices, details = {
                                                            "anchor" : "group median" if self._anchor is None else self._anchor,
                                                            "upper_bound" : str(self._upper),
                                                            "lower_bound" : str(self._lower), 
                                                        }
                                                    )

        return self._Assay

    
    def _set_bounds(self, anchor):
        """
        Set upper and lower boundaries of inclusion_range
        """
        upper, lower = anchor + self._upper, anchor - self._lower
        return upper,lower

    def _get_anchor(self, kwargs, group, tmp):
        """
        Set anchor for inclusion range
        """
        Ct = defaults.raw_col_names[1]
        if self._anchor is None:
            anchor = np.median(tmp[Ct])
        elif type(self._anchor) == type(print):
            anchor = self._anchor(tmp, **kwargs)
        elif isinstance(self._anchor, (list, tuple, dict)):
            anchor = self._anchor[group]
        elif isinstance(self._anchor, (int, float)):
            anchor = self._anchor
        return anchor

class IQRFilter(Filter):
    """
    Filters out outliers based on the classical n x IQR (with n = 1.5 by default) approach.
    """
    def __init__(self):
        super().__init__()
        self._upper = 1.5
        self._lower = 1.5

    def _filter(self, **kwargs):
        """
        Gets IQR for each group and finds outliers based on self._upper / lower
        """

        plotter = self._BoxPlotter

        plotter.add_before( self._Assay )
    
        df = self._Assay.get()
        groups = self._Assay.groups()
        Ct = defaults.raw_col_names[1]
        
        faulty_indices = []
        for group in groups:
            tmp = df.query(f"group == {group}")

            # get anchor
            anchor = np.nanmedian(tmp[Ct])
            # ignore Nan if so specified
            if self._ignore_nan and anchor != anchor: 
                continue

            # generate inclusion range boundries
            first, third = np.nanquantile(tmp[Ct], 0.26), np.nanquantile(tmp[Ct], 0.76)
            upper, lower = self._set_bounds(anchor, first, third)
            
            # get faulty replicates
            faulty_replicates = tmp.query(f"{Ct} < {lower} or {Ct} > {upper}")
            faulty_indices.extend(list(faulty_replicates.index))

            self._save_stats(self._Assay.id(), group, anchor, upper, lower)        
        
        self._filter_out(faulty_indices)

        plotter.add_after( self._Assay )


        if self._report_loc is not None: 
            self._write_report(faulty_indices, details = {
                                                            "upper_max" : str(self._upper),
                                                            "lower_max" : str(self._lower), 
                                                        }
                                                    )

        return self._Assay
    
    def _set_bounds(self, anchor, first, third):
        """
        Set upper and lower boundaries of inclusion_range
        """
        iqr = third - first
        upper, lower = anchor + iqr * self._upper, anchor - iqr * self._lower
        return upper,lower




def filter( assay, mode: str = "range", lim: (float or tuple) = None, anchor = None, ignore_nan : bool = True, drop_outliers : bool = False ):
    """
    Filter a single or multiple `qpcr.Assay` objects using default `Filter` setups.

    Parameters
    ----------
    assay : qpcr.Assay or list
        A single `qpcr.Assay` object or a list thereof.
    
    mode : str
        Either `"range"` to call a `RangeFilter` that uses a static range to filter values, or `"iqr"` to call a `IQRFilter` that uses the Interquantile Range
        to filter values. 
    
    lim : float or tuple
        The filtering limits for the inclusion range. Any values outside of these will be filtered out. For the `RangeFilter` these are absolute values around the `anchor` (+/- 1 by default)
        while for the `IQRFilter` these are scalars (+- 1.5 by default) for the IQR. If a single `float` is supplied the limits are interpreted symmetrically, while a `tuple` is read as first lower bound, then upper bound.
    
    anchor : 
            Only used for `RangeFilters`. Supported types for `anchor` are: a numeric value (`int or float`),
            an `iterable` of same length as groups in the dataframe, 
            a `dict` where keys must be numeric group identifiers (starting from 0) and values are numeric values to be used as anchor (`int or float`),
            or a `function` that works with a pandas dataframe as stored by `qpcr.Assay` objects, 
            which must return a single numeric value for the anchor (it will be applied to replicate-grouped subsets of the total dataframe).
    
    ignore_nan : bool
        Ignore NaN values when computing inclusion ranges. If set to `False` a single NaN value will render the entire replicate group unfilterable!
    
    drop_outliers : bool
        If `True` entries are actually dropped from the dataframe. By default any entries that do not match the inclusion range are set to NaN.


    Parameters
    ----------
    assay : qpcr.Assay or list
        The same as input but with filtered dataframes.
    """
    f = RangeFilter() if mode == "range" else IQRFilter

    if lim is not None:
        if isinstance(lim, (tuple, list, np.ndarray)): 
            f.set_lim( upper = lim[0], lower = lim[1] )
        else: 
            f.set_lim( lim = lim )
    if mode == "range" and anchor is not None:
        f.set_anchor( anchor = anchor )

    f.ignore_nan( ignore_nan )
    f.drop_outliers( drop_outliers ) 

    return f.pipe( assay )


if __name__ == "__main__":
    
    normalisers = ["./Examples/Example Data/28S.csv", "./Examples/Example Data/actin.csv"]
    assays = ["./Examples/Example Data/HNRNPL_nmd.csv", "./Examples/Example Data/HNRNPL_prot.csv"]

    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    reader = qpcr.DataReader()
    assays = [ reader.read(i, replicates = 6) for i in assays ]
    normalisers = [ reader.read(i, replicates = 6) for i in normalisers ]

    filter = RangeFilter()
    # filter.plotmode("static")
    assays = [ filter.pipe(i) for i in assays ]
    normalisers = [ filter.pipe(i) for i in normalisers ]

    analyser = qpcr.Analyser()
    assays = [ analyser.pipe(i) for i in assays ]
    normalisers = [ analyser.pipe(i) for i in normalisers ]

    normaliser = qpcr.Normaliser()
    normaliser.link(assays, normalisers)
    normaliser.normalise()

    fig = filter.plot(show = False)

    print( type( fig ))

    # prev = Plotters.PreviewResults("interactive")
    # prev.link( normaliser.get() )
    # prev.plot()
