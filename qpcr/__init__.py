"""
This module is designed to provide functions to analyse qPCR data. 
It is designed for maximal user-friendliness and streamlined data-visualisation.
"""

import pandas as pd
import qpcr._auxiliary as aux
from qpcr._auxiliary import warnings as aw
import os
import numpy as np 
from copy import deepcopy 


# default column names for raw Ct data files
RAW_COL_NAMES = ["Sample", "Ct"]

class Reader(aux._ID):
    """
    Reads qpcr raw data files in csv format. 

    Parameters
    ----------
    filename : str
        A filename to a csv containing Ct values. 
        The file has to have two named columns; one for sample names, one for Ct values. 
        Both csv (, spearated) and csv2 (; separated) are accepted.
    """
    def __init__(self, filename:str) -> pd.DataFrame: 
        super().__init__()
        self._src = filename
        self._delimiter = ";" if self._is_csv2() else ","
        self.read()

    def get(self):
        """
        Returns
        -------
        data : pd.DataFrame
            The samples dataframe.
        """
        return self._df

    def n(self):
        """
        Returns
        -------
        n : int
            The number of samples in the dataframe.
        """
        return len(self._df["Sample"])

    def read(self):
        """
        Reads the given data file
        """
        self._df = pd.read_csv(
                                self._src, 
                                sep = self._delimiter, 
                                header = self._has_header(), 
                                names = RAW_COL_NAMES
                            )
        # self._df["_index"] = list(self._df.index)

    def _is_csv2(self):
        """
        Tests if csv file is ; delimited (True) or common , (False)
        """
        with open(self._src, "r") as openfile: 
            content = openfile.read()
        if ";" in content: 
            return True
        return False

    def _has_header(self):
        """
        Checks if column headers are provided in the data file
        It does so by checking if the second element in the first row is numeric
        if it is numeric (returns None << False) no headers are presumed. Otherwise
        it returns 0 (as in first row has headers)...
        """
        with open(self._src, "r") as openfile: 
            content = openfile.read().split("\n")[0]
            content = content.split(self._delimiter)
        try: 
            second_col = content[1]
            second_col = float(second_col)
        except ValueError:
            return 0 # Headers in row 0
        return None  # no headers


class Assay(aux._ID):
    """
    Places a set of samples into groups of replicates as specified by the user.
    Also adds a "group" (numeric) column and "group_name" (string) to the Reader dataframe that specifies the replicate groups. 
    Optionally, users may re-name the groups manually (otherwise Group1,... will be used by default).

    Parameters
    ----------
    Reader : qpcr.Reader
        A qpcr.Reader object (optional). 
        Reader objects can also be linked to the assay after setup iteratively.

    """
    def __init__(self, Reader:Reader = None) -> dict:
        super().__init__()
        self._Reader = Reader
        if Reader is not None:
            self.adopt_id(Reader)
        self._df = None
        self._replicates = None
        self._renamed = False

    def get(self):
        """
        Returns
        -------
        data : pd.DataFrame
            The stored dataframe
        """
        return self._df

    def link(self, Reader:Reader):
        """
        Links a qpcr.Reader object to the Assay.

        Parameters
        ----------
        Reader : qpcr.Reader
            A qpcr.Reader object.
        """
        self._Reader = Reader
        self.adopt_id(Reader)

    def names(self, as_set = True):
        """
        Returns a set of sample group names (maintaing group order).

        Parameters
        ----------
        as_set : bool
            If `as_set = True` (default) it returns a set (as list without duplicates) 
            of assigned group names for replicate groups.
            If `as_set = False` it returns the full group_name column (including all repeated entries).
        
        Returns
        -------
        names : list or pd.Series
            The given group names of all replicate groups.
        """
        if as_set:
            return aux.sorted_set(list(self._df["group_name"]))
        else: 
            return self._df["group_name"]
    
    def is_named(self): # not used so far...
        """
        Returns 
        -------
        bool
            `True` if .rename() was performed and custom group names are provided
            else `False`.
        """
        return self._renamed

    def groups(self):
        """
        Returns a set of sample groups (numeric).

        Returns
        -------
        groups : list
            The given numeric group identifiers of all replicate groups.
        """
        return sorted(list(set(self._df["group"])))

    def replicates(self, replicates : (int or tuple) = None):
        """
        Either sets or gets the replicates to be used for grouping the samples
        Before they are assigned, replicates are vetted to ensure they cover all data entries.

        Parameters
        ----------
        replicates : int or tuple
            Can be an `integer` (equal group sizes, e.g. `3` for triplicates), 
            or a `tuple` (uneven group sizes, e.g. `(3,2,3)` if the second group is only a duplicate). 
        """
        if replicates is None:
            return self._replicates
        else: 
            if self._vet_replicates(replicates):
                self._replicates = replicates
            else: 
                aw.HardWarning("Assay:reps_dont_cover", n_samples = self._Reader.n(), reps = replicates)

    def group(self):
        """
        Groups the samples according to replicates specified.
        """
        df = self._Reader.get()
        
        # generate group and group_names columns
        if isinstance(self._replicates, int):
            samples = self._Reader.n()
            groups, group_names = self._make_equal_groups(samples)            
        elif isinstance(self._replicates, tuple):
            groups, group_names = self._make_unequal_groups()
        else:
            aw.HardWarning("Assay:no_reps_yet")

        df["group"], df["group_name"] = groups, group_names
        self._df = df

    def rename(self, names:(list or dict)):
        """
        Replaces the generic Group0,... in the "group_name" column.

        Parameters
        ----------
        names : list or dict
            Either a `list` (new names without repetitions) or `dict` (key = old name, value = new name) specifying new group names. 
            Group names only need to be specified once, and are applied to all replicate entries.
        """
        # get new group names based on list (index) or dict (key)
        if isinstance(names, (list, tuple, set)):
            new_names = self._rename_per_index(names)       
        elif isinstance(names, dict):
            new_names = self._rename_per_key(names)
        else:
            aw.HardWarning("Assay:no_groupname_assignment", names = names)

        # update "group_name"
        self._df["group_name"] = new_names
        self._renamed = True

    def ignore(self, entries:tuple):
        """
        Remove lines based on index from the dataframe.
        This is useful when removing corrupted data entries.

        Parameters
        ----------
        entries : tuple
            Tuple of row indices from the dataframe to drop.
        """
        self._df = self._df.drop(index = list(entries))
        
    def _rename_per_key(self, names):
        """
        Generates new name list based on current names in "group_name" and uses string.replace()
        to update groupnames, based on key (old name) : value (new name) indexing. 
        Before applying it checks if all groups are covered by new names
        """
        current_names = aux.sorted_set(self._df["group_name"])
        all_groups_covered = len(names) == len(current_names)
        if all_groups_covered:
            current_names = list(self._df["group_name"])
            new_names = "$".join(current_names)
            for old_name, new_name in names.items():
                new_names = new_names.replace(old_name, new_name)
            new_names = new_names.split("$")       
            return new_names
        else:
            aw.HardWarning("Assay:groupnames_dont_colver", current_groups = current_names, new_received = names)

    def _rename_per_index(self, names):
        """
        Generates new name list based on current names in "group_names" and uses string.replace()
        to update groupnames to new names based on index (using a the order 
        of groups as is currently present in "group_name"). 
        """
        current_names_set = aux.sorted_set(self._df["group_name"])
        all_groups_covered = len(names) == len(current_names_set)
        if all_groups_covered:
            current_names = list(self._df["group_name"])
            new_names = "$".join(current_names)
            names = list(names)
            for old_name, new_name in zip(current_names_set, names):
                new_names = new_names.replace(old_name, new_name)
            new_names = new_names.split("$")
            return new_names
        else:
            aw.HardWarning("Assay:groupnames_dont_colver", current_groups = current_names_set, new_received = names)

    def _make_unequal_groups(self):
        """
        Returns two lists of [0,0,0,1,1,1] and 
        [Group0, Group0, Group0, Group1,...] 
        to cover all sample entries.
        (this function works with a tuple for replicate group sizes)
        """
        groups = []
        group_names = []
        for rep, idx in zip(self._replicates, range(len(self._replicates))): 
            groups.extend([idx] * rep)
            group_names.extend([f"Group{idx}"] * rep)
        return groups, group_names

    def _make_equal_groups(self, samples):
        """
        Returns two lists of [0,0,0,1,1,1] and 
        [Group0, Group0, Group0, Group1,...] 
        to cover all sample entries.
        (this function works with an integer group size, 
        assuming all groups have the same size)
        """
        groups = []
        group_names = []
        slices = range(int(samples / self._replicates))
        for i in slices:
            groups.extend([i] * self._replicates)
            group_names.extend([f"Group{i}"] * self._replicates)
        return groups, group_names

    def _vet_replicates(self, replicates : (int or tuple)):
        """
        Checks if provided replicates will place all sample entries into a group
        returns True if all samples are covered, False if not...
        """
        samples = self._Reader.n()

        # for INT -> modulo will be 0 if all samples are covered
        # for TUPLE -> sum(replicates) should cover all samples...

        if isinstance(replicates, int):
            verdict = True if samples % replicates == 0 else False
        elif isinstance(replicates, tuple): 
            verdict = True if sum(replicates) == samples else False
        return verdict

class SampleReader(Assay):
    """
    Sets up a Reader+Assay pipeline that reads in a sample file and handles the 
    stored raw data in a pandas dataframe. 
    Its `read()` method directly returns a `qpcr.Assay` object that can be piped to Analyser. 
    Note
    ----
    This is the suggested to read in data, instead of manually setting up Reader and Assay objects.
    """
    def __init__(self):
        super().__init__()
        self._replicates = None
        self._names = None
        self._Reader = None
        self._Assay = None

    def replicates(self, replicates:(int or tuple)):
        """
        Set the replicates to group samples.

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
        
    def read(self, filename):
        """
        Reads one raw datafile (csv format).

        Parameters
        ----------
        filename : str
            A filename to a csv containing Ct values. 
            The file has to have two named columns; one for sample names, one for Ct values. 
            Both csv (, spearated) and csv2 (; separated) are accepted.

        Returns
        -------
        Assay : qpcr.Assay
            A `qpcr.Assay` object containing the grouped and renamed data.
        """
        self._Reader = Reader(filename)
        self._Reader.id(aux.fileID(filename))

        self._Assay = Assay(self._Reader)
        self._Assay.adopt_id(self._Reader)

        if self._replicates is not None:
            self._Assay.replicates(self._replicates)
            self._Assay.group()
        
        if self._names is not None:
            self._Assay.rename(self._names)

        return self._Assay


class Results(aux._ID):
    """
    Handles a pandas dataframe for the results from qpcr.Analyser.
    """
    def __init__(self):
        super().__init__()
        self._df = None
        self._Assay = None
        self._stats_results = {"group" : [], "assay" : [], "mean" : [], "stdev" : [], "median" : []}
        self._stats_df = None

    def adopt_names(self, Assay:Assay):
        """
        Links an instance of Assay to be used as reference for group_names
        It copies the group_name column to the results storing dataframe.
        This step can only be performed once!

        Parameters
        ----------
        Assay : qpcr.Assay
            A `qpcr.Assay` object whose group_name column will be copied.
        """
        self._Assay = Assay
        if self.is_empty():
            self._df = self._Assay.get()
            self._drop_setup_cols()
        else:
            aw.SoftWarning("Results:cannot_link")
    
    def is_named(self):
        """
        Note
        ----
        This is primarily a legacy function from a development-version that discarded `group_name` if no custom names were provided in `qpcr.Assay` objects.
        `qpcr.Assay` objects now retain the `group_name` column in any case, but migrating names into `qpcr.Results` is still an additional step performed by `adopt_names()`.
        
        Returns 
        -------
        bool
            `True` if `group_name` column is present in the Results dataframe, else `False`.
        """
        return "group_name" in self._df.columns
    
    def names(self, as_set = False):
        """
        Returns 
        -------
        names : list or None
            The adopted `group_names` (only works if a `qpcr.Assay` have been linked using `adopt_names()`!)
        """
        if self._Assay is not None:
            return self._Assay.names(as_set)
        return None

    def get(self):
        """
        Returns 
        -------
        data : pd.DataFrame
            The results dataframe
        """
        return self._df

    def is_empty(self):
        """
        Checks if any results have been stored so far.

        Returns
        -------
        bool
            `True` if NO data is yet stored, else `False`.
        """
        return self._df is None

# ah for whatever funcking reason it tries to add the HNRNPL_rel28fuck twice
# dunno why, sucks big time... 

    def add(self, column:pd.Series):
        """
        Adds a new column of either DeltaCt 
        computed data or normalised DeltaCt data to the results dataframe.

        Parameters
        ----------
        column : pd.Series
            A named pandas Series or DataFrame that can be joined into the already
            stored dataframe.
        """
        # print(self._df)
        # print(column)
        self._df = self._df.join(column)
        
    def merge(self, *Results):
        """
        Merges any number of other qpcr.Results objects into this one.
        The source id of the results is added as column-name suffix. 

        Parameters
        ----------
        *Results
            An arbitrary number of qpcr.Results objects.

        """
        new_df = self._df
        for R in Results: 
            R_df = R.get()
            # we merge the dataframes based on their groups, and add the instance id as identifier
            new_df = pd.merge(new_df, R_df["dCt"], 
                                right_index = True, left_index = True, 
                                suffixes = [f"_{self.id()}", f"_{R.id()}"]
                            )
        self._df = new_df

    def drop_cols(self, *cols):
        """
        Drops all specified columns from the dataframes
        this is used for normaliser pre-processing.

        Parameters
        ----------
        *cols
            Any column names (as `str`) to be dropped.
            If no names are specified any/all `deltaCt` data-containing columns are dropped!
            If this is the case then the only columns retained are: `"group", "group_name", "Sample", "assay"`.
        """
        if cols == ():
            _to_drop = [c for c in self._df.columns if c not in ["group", "group_name", "Sample", "assay"]]
        else:
            _to_drop = [c for c in list(cols) if c in list(self._df.columns)]
        self._df = self._df.drop(columns = _to_drop)
        
    def rename_cols(self, cols:dict):
        """
        Renames all columns according to a dictionary as key -> value.

        Parameters
        ----------
        cols : dict
            A dictionary specifying old column names (keys) and new colums names (values).
        """
        self._df = self._df.rename(columns = cols)


    def stats(self, recompute = False) -> pd.DataFrame:
        """
        Computes summary statistis about the replicate groups: 
        `Mean`, `Median`, and `StDev` of all replicate groups, for all assays.
        
        Parameters
        ----------
        recompute : bool
            Statistics will only be once unless recompute is set to `True`.

        Returns
        -------
        stats_df : pd.DataFrame
            A new dataframe containing the computed statistics for each replicate group.

        """
        # if stats_df is already present, return but sorted according to samples, not groups (nicer for user to inspect)
        if self._stats_df is not None and not recompute:
            return self._stats_df.sort_values("assay")
        
        # get groups and samples 
        groups = aux.sorted_set(list(self._df["group"]))
        samples = [c for c in self._df.columns if c not in ["Sample", "group", "group_name", "assay"]]
     
        # compute stats for all samples per group
        for group in groups:
            group_subset = self._df.query(f"group == {group}")
            
            median = self._stat_var(group_subset, np.nanmedian)
            mean = self._stat_var(group_subset, np.nanmean)
            stdv = self._stat_var(group_subset, np.nanstd)
            self._add_stats(samples, group, median, mean, stdv)
            
        # add group names if present
        if self.is_named():
            self._add_stats_names(samples)

        self._stats_df = pd.DataFrame(self._stats_results)
        return self._stats_df.sort_values("assay")

    def save(self, path, df = True, stats = True):
        """
        Saves a csv file for each specified type of results.

        Parameters
        ----------
        path : str
            Path has to be a filepath if only one type of results shall be saved (i.e. either `df` or `stats`), 
            otherwise a path to the directory where both `df` and `stats` shall be saved.
        
        df : bool
            Save the results dataframe containing all replicate values (the full results).
            Default is `df = True`.
        
        stats : bool
            Save the results dataframe containing summary statistics for all replicate groups.
            Default is `stats = True`.
        
        """
        if df and stats and not os.path.isdir(path):
            aw.HardWarning("Results:save_need_dir")

        if df:
            self._save_single(path, self._df, "_df")
        if stats:
            if self._stats_df is None:
                self.stats()
            self._save_single(path, self._stats_df, "_stats")

    def drop_rel(self):
        """
        Crops the `X_rel_Y` column-names to just `X`.
        """
        colnames = self._df.columns
        to_change = {i : i.split("_rel_")[0] for i in colnames if "_rel_" in i }
        self.rename_cols(to_change)

    def split(self, reset_names = False, drop_rel = True):
        """
        Splits the stored results dataframe into separate qpcr.Results objects containing only a signle deltaCt column each.

        Parameters
        ----------
        reset_names : bool
            Resets the deltaCt column-name from `"X_rel_Y"` to just `"dCt"`.

        drop_rel : bool
            Crops `"X_rel_Y"` deltaCt column-names to just `"X"`. 

        Returns 
        -------
        objects : list
            A list of qpcr.Results objects containing only a single dCt column each (retaining group columns etc.)
        """
        shared_columns = [i for i in self._df.columns if i in ["group", "group_name", "Sample", "assay"]]
        dct_columns = [i for i in self._df.columns if i not in ["group", "group_name", "Sample", "assay"]]
        
        dfs = [self._df[shared_columns + [i]] for i in dct_columns]
        objects = [Results() for i in dfs]

        for o, df, dct_col in zip(objects, dfs, dct_columns): 
            o._df = df
            if reset_names:
                o.rename_cols({dct_col : "dCt"})
            if drop_rel: 
                o.drop_rel()

            o.id(dct_col)
        
        return objects

    def _save_single(self, path, src, suffix=""):
        """
        Saves either self._df or self._stats_df to a csv file based on a path
        (path can be either filename or directory)
        """
        filename = path if not os.path.isdir(path) else os.path.join(path, f"rel_{self.id()}{suffix}.csv")
        src.to_csv(filename)
        
    def _drop_setup_cols(self):
        """
        Removes unnnecessary columns from the Sample df during self._df setup with link()
        """
        self.drop_cols("Ct")


    def _add_stats_names(self, samples):
        """
        Adds a group_name column to self._stats_result with appropriate
        repetition of group_names for each sample...
        """
        self._stats_results["group_name"] = []
        group_names = aux.sorted_set(list(self._df["group_name"]))
        for group_name in group_names:
            self._stats_results["group_name"].extend([group_name] * len(samples))

    def _add_stats(self, samples, group, median, mean, stdv):
        """
        Adds new summary entries to self._stats_results
        """
        self._stats_results["group"].extend([group] * len(samples))
        self._stats_results["assay"].extend(samples)
        self._stats_results["median"].extend(median)
        self._stats_results["mean"].extend(mean)
        self._stats_results["stdev"].extend(stdv)


    def _stat_var(self, group_subset, func, **kwargs):
        """
        Performs a function (like mean or stdv) over all rows
        and returns the result as list with a float for each column in the df
        any function can be passed as long as it works with an iterable
        """
        # ignore group and group_name columns
        ignore = ["Sample", "group", "group_name", "assay"]
        all_cols = [g for g in group_subset.columns if g not in ignore]
        tmp = group_subset[all_cols]
        # compute stats based on func
        stats = [func(tmp[col], **kwargs) for col in tmp.columns]
        return stats
        

class Analyser(aux._ID):
    """
    Performs Single Delta CT (first normalisation within dataset) 
    Note
    ----
    Delta Delta CT (normalisation using second dataset), is handled by qpcr.Normaliser!

    Parameters
    ----------
    Assay : qpcr.Assay
        A `qpcr.Assay` object (optional) to compute DeltaCT on. 
        `qpcr.Assay` objects can be iteratively linked subsequently using `link()`.
    """
    def __init__(self, Assay:Assay = None):
        super().__init__()
        self._Assay = Assay
        self._Results = Results()

        # default settings
        self._anchor = "first"
        self._efficiency = 2
        self._deltaCt_function = self._get_deltaCt_function(exp = True)

        if self._Assay is not None: 
            self._Results.adopt_id(Assay)
            self._Results.adopt_names(self._Assay)
    
    def get(self):
        """
        Returns 
        -------
        Results
            A `qpcr.Results` object that contains the deltaCT results
        """
        return self._Results

    def has_results(self):
        """
        Returns
        -------
        bool
            `True` if any results were already computed, else `False`.
        """
        return not self._Results.is_empty()

    def link(self, Assay:Assay, force = False, silent = False):
        """
        Links a `qpcr.Assay` object to the Analyser.
        Note
        ----
        If there are any precomputed results, no new data will be linked, unless force=True is called. 
        The user is notified if results are already present and how to proceed. 

        Parameters
        ----------
        Assay : qpcr.Assay
            A `qpcr.Assay` object containing data.
        
        force : bool
            Any already linked `qpcr.Assay` objects (and their data and results) will be overwritten
            if `force = True` (default is `force = False`).
        
        silent : bool
            Warnings about overwriting data will be suppressed if `silent = True` (default is `silent = False`).
            This is only relevant if `force = True`. 
        """
        empty = self._Results.is_empty()
        dont_overwrite = not empty and not force
        if not dont_overwrite:
            self._Assay = Assay
            self.adopt_id(self._Assay)
            self._Results = Results()
            self._Results.adopt_names(self._Assay)
            self._Results.adopt_id(self._Assay)
            
        if not silent:
            # notify the user of changes to the Analyser data and results
            if not dont_overwrite and not empty:
                aw.SoftWarning("Analyser:newlinked")
            elif dont_overwrite and not empty:
                aw.SoftWarning("Analyser:not_newlinked")

    def pipe(self, Assay:Assay, **kwargs) -> Results:
        """
        A quick one-step implementation of link + DeltaCt.
        This is the suggested application of the `qpcr.Analyser` class!

        Note
        ----
        This will silently overwrite any previous results! 

        Parameters
        ----------
        Assay : qpcr.Assay
            A `qpcr.Assay` object to be linked to the Analyser for DeltaCt computation.
        
        **kwargs
            Any additional keyword arguments to be passed to the `DeltaCt()` method.

        Returns 
        -------
        results : qpcr.Results
            A `qpcr.Results` object. 
        """
        self.link(Assay, force=True, silent=True)
        self.DeltaCt(**kwargs)
        return self.get()

    def efficiency(self, e:float = None):
        """
        Sets an efficiency factor for externally calculated qPCR amplification efficiency.
        By default `efficiency = 2` is assumed.

        Parameters
        ----------
        e : float
            An amplification efficiency factor. Default is `e = 2`.

        """
        if isinstance(e, (int, float)):
            self._efficiency = float(e)
        elif e is None: 
            return self._efficiency

    def anchor(self, anchor):
        """
        Sets the anchor for DeltaCt for internal normalisation.

        Parameters
        ----------
        anchor : str or float
            The internal anchor for normalisation.
            This can be either `"first"` (default, the very first dataset entry),
            `"grouped"` (first entry for each replicate group), or 
            any specified numeric value (as `float`).
        """
        self._anchor = anchor

    def func(self, f:(str or function)):
        """
        Sets the function to be used for DeltaCt (optional)

        Parameters
        ----------
        f : str or function
            The function to be used for DeltaCt computation. Pre-defined functions are 
            either `"exponential"` (which uses  `efficiency^(-(s-r))`, default), or `"linear"` 
            (uses uses `s-r`), where `s` is any sample entry in the dataframe and `r` is the anchor.
            It is also possible to assign any defined function that accepts an anchor (1st!) and sample (2nd!) 
            numeric value each, alongside any kwargs (which will be forwarded from DeltaCt()...).
        """
        if f in ["exponential", "linear"]:
            f = True if f == "exponential" else False
            self._deltaCt_function = self._get_deltaCt_function(f)
        elif type(f) == type(aux.fileID):
            self._deltaCt_function = f
        else:
            aw.HardWarning("Analyser:cannot_set_func", func = f)

    def DeltaCt(self, **kwargs):
        """
        Calculates DeltaCt for all groups within the dataframe.
        Any specifics such as `anchor` or `func` must have already been 
        set using the respective methods prior to calling `DeltaCt()`!

        Parameters
        ----------
        **kwargs
            Any additional keyword arguments that a custom DeltaCt function may require.
        """
        if self._anchor == "first":
            self._DeltaCt_first_anchored(self._deltaCt_function, **kwargs)
        elif self._anchor == "grouped":
            self._DeltaCt_grouped_anchored(self._deltaCt_function, **kwargs)
        else: 
            self._DeltaCt_externally_anchored(self._anchor, self._deltaCt_function, **kwargs)


    def _DeltaCt_externally_anchored(self, anchor:float, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using a specified anchor
        """
        df = self._Assay.get()
        df["dCt"] = df["Ct"].apply(deltaCt_function, ref = anchor, **kwargs)
        self._Results.add(df["dCt"])


    def _DeltaCt_grouped_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the first entry of each group as anchor
        """
        # get set of sample groups and dataset
        groups = self._Assay.groups()
        df = self._Assay.get()

        dCt = pd.Series()
        for group in groups: 
            group_subset = df.query(f"group == {group}")
            anchor = group_subset["Ct"][0]
            delta_cts = group_subset["Ct"].apply(deltaCt_function, ref=anchor, **kwargs)
            dCt.append(delta_cts)
        self._Results.add(dCt)

    def _DeltaCt_first_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the very first entry of the dataset as anchor
        """
        df = self._Assay.get()
        anchor = df["Ct"][0]
        df["dCt"] = df["Ct"].apply(deltaCt_function, ref=anchor, **kwargs)
        self._Results.add(df["dCt"])

    def _exp_DCt(self, sample, ref, **kwargs):
        """
        Calculates deltaCt exponentially
        """
        factor = sample-ref 
        return self._efficiency **(-factor)

    def _simple_DCt(self, sample, ref, **kwargs):
        """
        Calculates deltaCt linearly
        """
        return sample-ref

    def _get_deltaCt_function(self, exp):
        """
        Returns the function to be used for DeltaCt based on 
        whether or not exponential shall be used.
        """
        if exp == True:
            dCt = self._exp_DCt
        else:
            dCt = self._simple_DCt
        return dCt


class Normaliser(aux._ID):
    """
    Handles the second normalisation normalisation of two (or more) datasets against one of them used as normaliser.
    This requires that all have been analysed in the same way before!
    """
    def __init__(self):
        super().__init__()
        self._Normalisers = []
        self._Assay = []
        self._Results = Results()
        self._normaliser = None
        self._prep_func = self._average
        self._norm_func = self._divide_by_normaliser

    def get(self, copy=False):
        """
        Parameters
        ----------
        copy : bool
            Will return a deepcopy of the Results object if `copy = True` (default is `copy = False`).
        
        Returns
        -------
        Results : qpcr.Results
            A `qpcr.Results` object containing the normalised dataframe
        """
        if copy: 
            return deepcopy(self._Results)
        return self._Results
    
    def link(self, samples:(list or tuple) = None, normalisers:(list or tuple) = None):
        """
        Links either normaliser- or sample-data-containing `qpcr.Results` objects coming from the same `qpcr.Analyser`.

        Parameters
        ----------
        samples : list or tuple
            A list of `qpcr.Results` objects coming from a `qpcr.Analyser` which shall be normalised against a normaliser.
        
        normalisers : list or tuple
            A list of `qpcr.Results` objects coming from a `qpcr.Analyser` which shall be used as normalisers. These will be
            combined into one single pseudo-normaliser which will then be used to normalise the samples. The method of 
            combining the normalisers can be specified using the `prep_func()` method.
        """
        self._link_normaliser(normalisers)
        self._link_samples(samples)
    
    def prep_func(self, f = None):
        """
        Sets any defined function for combined normaliser pre-processing.
        If no `f` is provided, it returns the current `prep_func`.

        Parameters
        ----------
        f : function
            The function may accept one list of qpcr.Results objects, and must return 
            one list (or iterable) of the same length as entries within the qpcr.Results dataframes.
            If the provided function does adhere to these criteria is NOT vetted by this method!
        """
        if type(f) == type(aux.fileID):
            self._prep_func = f
        elif f is None:
            return f
        else: 
            aw.HardWarning("Normaliser:cannot_set_prep_func", func = f)

    def norm_func(self, f = None):
        """
        Sets any defined function to perform normalisation of samples against normalisers.
        If no `f` is provided, it returns the current `norm_func`.

        Parameters
        ----------
        f : function
            The function may accept one numeric entry for a sample and a normaliser, and must return 
            a numeric value. By default `s/n` is used, where `s` is a column of sample deltaCt values, and `n` is the corresponding deltaCt column from the combined normaliser.
        """
        if type(f) == type(aux.fileID):
            self._norm_func = f
        elif f is None:
            return f
        else: 
            aw.HardWarning("Normaliser:cannot_set_norm_func", func = f)

    def normalise(self, **kwargs):
        """
        Normalises all linked samples against the combined pseudo-normaliser, and stores the results in a new Results object.

        Parameters
        ----------
        **kwargs
            Any additional keyword arguments that may be passed to a custom `norm_func`.
        """
        if self._normaliser is None: 
            self._preprocess_normalisers()

        if self._Assay == [] or self._normaliser is None:
            aw.SoftWarning("Normaliser:no_data_yet")

        # get normaliser dataframe
        normaliser = self._normaliser.get()

        # setup groups for _Results
        self._Results.adopt_names(self._Assay[0])
        self._Results.drop_cols()
        # print(self._Results.get())

        # combine normalised samples into unified dataframe
        for S in self._Assay:
            S_df = S.get()
            column_name = f"{S.id()}_rel_{self._normaliser.id()}"
            normalised = self._norm_func_wrapper(S_df, normaliser, **kwargs)
            normalised = normalised.rename(column_name)
            self._Results.add(normalised)

    def _norm_func_wrapper(self, sample_assay, normaliser, dCt_col="dCt", norm_col="dCt_combined"):
        """
        The wrapper that will apply the _norm_func to the sample and normaliser dataframes and return a normalised dataframe
        """
        # for double normalised we want the same columns as dct and norm...
        dCt_col, norm_col = self._prep_columns(sample_assay, dCt_col, norm_col)

        tmp_df = normaliser.join(sample_assay, lsuffix="_s")
        # tmp_df = sample_assay.join(normaliser, rsuffix = "_n")
        results = self._norm_func(tmp_df[[dCt_col, norm_col]])
        return results

    def _prep_columns(self, sample_assay, dCt_col, norm_col):
        """
        Returns the columns to use if named columns shall be used (named columns will be used for second-normalisation of entire runs)
        """
        if dCt_col == "named":
            dCt_col = [i for i in sample_assay.columns if i not in ["group", "group_name", "Sample", "assay"]]
            # assert len(dCt_col) == 1, f"length of dCt_col is: {len(dCt_col)}"
            dCt_col = dCt_col[0]

        if norm_col == "same": 
            norm_col = dCt_col + "_s"
        return dCt_col,norm_col

    def _divide_by_normaliser(self, df):
        """
        Performs normalisation of sample s against normaliser n
        s and n are specified as two pandas dataframe columns
        Note, that the dataframe must ONLY contain these two columns, first the dCt sample, then the normaliser!
        (default _norm_func)
        """
        dCt_col, norm_col = df.columns
        s, n = df[dCt_col], df[norm_col]
        return s / n

    def _link_samples(self, samples):
        """
        Links any provided samples and checks their datatype in the process...
        """
        if samples is not None:
            for sample in samples: 
                if isinstance(sample, Results):
                    self._Assay.append(sample)
                elif isinstance(sample, Analyser) and sample.has_results():
                    self._Assay.append(sample.get())
                elif isinstance(sample, Analyser) and not sample.has_results():
                    aw.SoftWarning("Normaliser:empty_data", s = sample)
                else: 
                    aw.SoftWarning("Normaliser:unknown_data", s = sample)
                
    def _link_normaliser(self, normalisers):
        """
        Checks if normaliser is provided and has proper datatype to be added...
        """
        if normalisers is not None:
            for normaliser in normalisers:
                if isinstance(normaliser, Results):
                    self._Normalisers.append(normaliser)
                elif isinstance(normaliser, Analyser) and normaliser.has_results():
                    self._Normalisers.append(normaliser.get())
                else: 
                    aw.SoftWarning("Normaliser:norm_unknown_data", s = normaliser)

    def _preprocess_normalisers(self):
        """
        Averages the provided normalisers row-wise for all normalisers into a 
        single combined normaliser, that will be stored as a Results instance.
        """
        combined = Results() # setup new dataframe for combined normalisers, intialise with first id
        combined.adopt_names(self._Normalisers[0])
        combined.adopt_id(self._Normalisers[0])
        combined.merge(*self._Normalisers[1:])

        tmp_df = self._prep_func(combined)
        tmp_df = tmp_df.rename("dCt_combined")
        combined.add(tmp_df)
        combined.drop_cols("group", "group_name", "Sample")

        self._normaliser = combined  
        if len(self._Normalisers) > 1:
            self._update_combined_id()
        
        # forward combined_id to self and _Results 
        self.adopt_id(self._normaliser)
        self._Results.adopt_id(self._normaliser)

    def _update_combined_id(self):
        """
        Generates a new id based on all normaliser ids,
        joining them as a+b+c,...
        """
        ids = [N.id() for N in self._Normalisers]
        ids = "+".join(ids)
        self._normaliser.id(ids)
        

    def _average(self, combined):
        """
        Averages row-wise all Normaliser entries and 
        generates a series of their per-row means
        (default preprocess_normalisers function)
        """
        tmp = combined.get()
        tmp_df = tmp.drop(columns = ["group"]) # drop group as it is a numeric column and would otherwise skew the average
        tmp_df = tmp_df.mean(axis = 1)
        return tmp_df


if __name__ == "__main__":
    
    files = ["Example Data/28S.csv", "Example Data/actin.csv", "Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]
    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    analysers = []

    reader = SampleReader()
    reader.replicates(6)
    reader.names(groupnames)

    analyser = Analyser()
    analyser.anchor("first")

    for file in files: 

        # reader = Reader(file)
        # reader.id(aux.fileID(file))

        # samples = Assay(reader)
        
        sample = reader.read(file)
        # sample.ignore((0,1,3,4))

        # analyser.link(sample, force=True, silent = False)
        # analyser.DeltaCt()
        # res = analyser.get()

        res = analyser.pipe(sample)
        # print(res)
        analysers.append(res)

    # for a in analysers: print(a.id(), "\n", a.get())

    normaliser = Normaliser()
    normaliser.link(normalisers = analysers[:2])
    normaliser.link(samples = analysers[2:])

    normaliser.normalise()
    
    result = normaliser.get()

    print(result.get())

    splitted = result.split(reset_names = False)

    i, j = splitted
    i.rename_cols({"HNRNPL_nmd": "HNRNPL"})
    j.rename_cols({"HNRNPL_prot": "HNRNPL"})

    # print(i, j)
    # print("-------")
    # print(i.get()["HNRNPL"] / j.get()["HNRNPL"])
    # print("-----")
    sn = Normaliser()

    sn.link(
        samples = [splitted[0]], 
        normalisers = [splitted[1]],
    )
    
    sn.normalise(dCt_col = "named", norm_col = "same")

    print(sn.get().get())
    print(sn.get().stats())

    # # result.save("..")
    
    # #result.add_names(samples)

    # print(result.stats())

    exit(0)