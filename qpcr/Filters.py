"""
This submodule defines a number of filters that can be used to 
remove faulty reads prior to deltact analysis. 

Filters implemented are currently
---------------------------------
RangeFilter 
    Filters out all raw Ct values that do not comply to a user-specified range (default is `+/-1` around replicate group median)
    The user has the option of specifying another anchor and limits for the inclusion range.

IQRFilter
    Filters out any outliers by `n x IQR`, where `n` is a scaling factor (default `n = 1.5`) around the replicate group median (anchor).
"""

from re import L
import qpcr.__init__ as qpcr
import pandas as pd
import numpy as np
import qpcr.auxiliary.warnings as aw
import qpcr.auxiliary as aux
import os 
import qpcr.Plotters

class Filter(aux._ID):
    """
    The super filtering class that takes in a qpcr.Assay object and updates its dataframe to a filtered version.
    """
    def __init__(self):
        super().__init__()
        self._Assay = None
        self._report_loc = None
        self._id = type(self).__name__
        self._boxplot_mode = "interactive"
        self._before_BoxPlotter = Plotters.ReplicateBoxPlot(Filter = self, mode = self._boxplot_mode)
        self._after_BoxPlotter = Plotters.ReplicateBoxPlot(Filter = self, mode = self._boxplot_mode)
        self._filter_stats = pd.DataFrame({
                                            "assay" : [], "group" : [], 
                                            "anchor" : [], "upper" : [], "lower" : []
                                        })
    
    def get_stats(self):
        """
        Returns 
        -------
        stats : pd.DataFrame
            The filtering statistics dataframe (a summary of filtering parameters used)
        """
        return self._filter_stats

    def plotmode(self, mode = "interactive"):
        """
        Set graph mode if a summary Boxplot shall be made

        Parameters
        ----------
        mode : str
            Can be either "interactive" (plotly) or "static" (matplotlib), or None to disable plotting.
        """
        self._boxplot_mode = mode
        self._before_BoxPlotter = Plotters.ReplicateBoxPlot(Filter = self, mode = self._boxplot_mode)
        self._after_BoxPlotter = Plotters.ReplicateBoxPlot(Filter = self, mode = self._boxplot_mode)
    
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
        figs = []
        filenames = [f"{self.id()}_before", f"{self.id()}_after"]
        for plotter, filename in zip([self._before_BoxPlotter, self._after_BoxPlotter], filenames):
            fig = plotter.plot(**kwargs)
            figs.append(fig)
            if self._report_loc is not None and self._boxplot_mode is not None: 
                suffix = plotter.suffix()
                plotter.save(os.path.join(self._report_loc, f"{filename}.{suffix}"))
        return figs

    def link(self, Assay:qpcr.Assay):
        """
        Links a qpcr.Assay to be filtered
        
        Parameters
        ----------
            Assay : qpcr.Assay
                A qpcr.Assay object to be filtered.
        """
        self._Assay = Assay
        self._before_BoxPlotter.link(self._Assay)

    def pipe(self, Assay:qpcr.Assay, **kwargs):
        """
        A shortcut for link+filter.
        This is the suggested usage for Filters.
        
        Parameters
        ----------
        Assay : qpcr.Assay
            A qpcr.Assay object to be filtered.
        **kwargs
            Any keyword arguments that should be passed to the plotting method.
        
        Returns
        -------
        Assay : qpcr.Assay
            A qpcr.Assay object containing only entries that passed the filter.

        """
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
        Assay : qpcr.Assay
            An updated qpcr.Assay object containing only entries that passed the filter.
        """
        if self._Assay is not None:
            self._filter(**kwargs)
            self._after_BoxPlotter.link(self._Assay)
            return self._Assay
        else: 
            aw.HardWarning("Filter:no_assay")

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
            self._Assay.ignore(faulty_indices)
    
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
            a `dict` where keys must be group indices (starting from 0) and values are a numeric value to be used as anchor (`int or float`),
            or a `function` that works with a pandas dataframe as stored by qpcr.Assay objects, 
            which must return a single numeric value for the anchor (it will be applied to replicate-grouped subsets of the total dataframe).
        """
        self._anchor = anchor

    def _filter(self, **kwargs):
        """
        Filters out any replicates that are out of range and updates the Assay's dataframe.
        It stores the original in a new attribute called _orig_df.
        """
    
        df = self._Assay.get()
        groups = self._Assay.groups()

        faulty_indices = []
        for group in groups:
            tmp = df.query(f"group == {group}")
            anchor = self._get_anchor(kwargs, group, tmp)
            upper, lower = self._set_bounds(anchor)
            faulty_replicates = tmp.query(f"Ct < {lower} or Ct > {upper}")
            faulty_indices.extend(list(faulty_replicates.index))

            self._save_stats(self._Assay.id(), group, anchor, upper, lower)

        self._filter_out(faulty_indices)

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
        if self._anchor is None:
            anchor = np.median(tmp["Ct"])
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

    def _filter(self):
        """
        Gets IQR for each group and finds outliers based on self._upper / lower
        """
    
        df = self._Assay.get()
        groups = self._Assay.groups()

        faulty_indices = []
        for group in groups:
            tmp = df.query(f"group == {group}")

            anchor = np.nanmedian(tmp["Ct"])
            first, third = np.nanquantile(tmp["Ct"], 0.26), np.nanquantile(tmp["Ct"], 0.76)
            upper, lower = self._set_bounds(anchor, first, third)
            
            faulty_replicates = tmp.query(f"Ct < {lower} or Ct > {upper}")
            faulty_indices.extend(list(faulty_replicates.index))

            self._save_stats(self._Assay.id(), group, anchor, upper, lower)        
        
        self._filter_out(faulty_indices)

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



if __name__ == "__main__":
    
    files = ["Example Data/28S.csv", "Example Data/actin.csv", "Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]
    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    analysers = []

    reader = qpcr.SampleReader()
    reader.replicates(6)
    reader.names(groupnames)

    analyser = qpcr.Analyser()
    analyser.anchor("first")

    iqr_filter = RangeFilter()
    iqr_filter.plotmode("static")
    iqr_filter.report(".")
    # iqr_filter.set_lim(1.6)

    for file in files: 
        
        sample = reader.read(file)
        
        # here comes in the filter...
        sample = iqr_filter.pipe(sample)
        # print(sample.get())

        res = analyser.pipe(sample)
        analysers.append(res)

    iqr_filter.plot(show = False)

    normaliser = qpcr.Normaliser()
    normaliser.link(normalisers = analysers[:2])
    normaliser.link(samples = analysers[2:])

    normaliser.normalise()
    
    result = normaliser.get()

    print("first result...")
    print(result.get())
    # print(result.stats())
    
    # just for fun, let's try to normalise nmd against prot...

    nmd = normaliser.get(copy=True)
    prot = normaliser.get(copy=True)

    nmd.drop_cols("HNRNPL_prot_rel_28S+actin", "assay")
    nmd.rename_cols(
        {"HNRNPL_nmd_rel_28S+actin" : "dCt"}
        )
    
    nmd.id("nmd")

    prot.drop_cols("HNRNPL_nmd_rel_28S+actin", "assay")
    prot.rename_cols(
        {"HNRNPL_prot_rel_28S+actin" : "dCt"}
    
        )

    prot.id("prot")

    # Alright: At the moment we cannot use a second_normaliser, 
    # as it for whatever reason overwrites ALL Ct containing 
    # columns with the last one??? >> happy bug hunting...

    # UDPATE: It seems like nmd is used for normaliser instead of prot??
    # 
    # SOLUTION: Solution is: use copy.deepcopy() otherwise both nmd 
    #           and prot are just referencing the same Results instance!

    second_normaliser = qpcr.Normaliser()
    second_normaliser.link(normalisers = [prot])
    second_normaliser.link(samples = [nmd, prot])

    second_normaliser.normalise()
    
    result = second_normaliser.get()

    print("second result...")
    print(result.get())

    # p = Plotters.PreviewResults(mode = "interactive")
    # p.link(result)
    # p.plot()

    p1 = Plotters.PreviewResults(mode = "static")
    p1.link(result)
    p1.plot()
    
    #  
    # result.save("..")
    
    # result.add_names(samples)


    exit(0)