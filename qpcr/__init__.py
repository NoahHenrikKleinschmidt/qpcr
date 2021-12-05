"""
This module is designed to provide functions to analyse qPCR data. 
It is designed for maximal user-friendliness and streamlined data-visualisation.
"""
import statistics as stat 
import pandas as pd
import auxiliary as aux
from auxiliary import warnings as aw
import glob, os
import copy
import statistics as stats

# TODO: A class to read csv pqcr raw data << CHECK
# TODO: A class to perform data preprocessing (like grouping replicates etc...) << CHECK
#   - grouping replicates << CHECK
#   - assigning group names (optional) << CHECK
# TODO: A class to perform delta delta ct << CHECK
# TODO: A class to handle and store results << CHECK
# TODO: A class to visulalise results 
# TODO: A compilation of common pipelines (like the qpcr.Analysis module from earlier)

RAW_COL_NAMES = ["Sample", "Ct"]


class Reader(aux._ID):
    """
    This class reads qpcr raw data files in csv format. 
    It requires that two columns are present, one for sample names, one for Ct values
    it will load these into a pandas dataframe.
    """
    def __init__(self, filename:str) -> pd.DataFrame: 
        super().__init__()
        self._src = filename
        self._delimiter = ";" if self._is_csv2() else ","
        self.read()

    def get(self):
        """
        Returns the samples dataframe
        """
        return self._df

    def n(self):
        """
        Returns the number of samples
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


class Samples(aux._ID):
    """
    This class groups a set of samples into groups of replicates as specified by the user.
    It adds a "group" (numeric) column and "group_name" (string) to the Reader dataframe that specifies the replicate groups. 
    Optionally, users may re-name the groups manually (otherwise Group1,... will be used by default)
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
        return self._df

    def link(self, Reader:Reader):
        """
        Links a qpcr.Reader object to the samples
        """
        self._Reader = Reader
        self.adopt_id(Reader)

    def names(self, as_set = True):
        """
        Returns a set of sample group names (maintaing group order)
        or using as_set=False the full group_name column with replicate repeats.
        """
        if as_set:
            return aux.sorted_set(list(self._df["group_name"]))
        else: 
            return list(self._df["group_name"])
    
    def is_named(self): # not used so far...
        """
        Returns True if .rename() was performed and custom group names are provided
        """
        return self._renamed

    def groups(self):
        """
        Returns a set of sample groups (numeric)
        """
        return sorted(list(set(self._df["group"])))

    def replicates(self, replicates : (int or tuple) = None):
        """
        Either sets or gets the replicates to be used for grouping the samples
        Before they are assigned, replicates are vetted to ensure they cover all data entries.
        This may either be an int (equal group sizes), a tuple (uneven group sizes)
        """
        if replicates is None:
            return self._replicates
        else: 
            if self._vet_replicates(replicates):
                self._replicates = replicates
            else: 
                aw.HardWarning("Samples:reps_dont_cover", n_samples = self._Reader.n(), reps = replicates)

    def group(self):
        """
        Groups the samples according to replicates specified
        """
        df = self._Reader.get()
        
        # generate group and group_names columns
        if isinstance(self._replicates, int):
            samples = self._Reader.n()
            groups, group_names = self._make_equal_groups(samples)            
        elif isinstance(self._replicates, tuple):
            groups, group_names = self._make_unequal_groups()
        else:
            aw.HardWarning("Samples:no_reps_yet")

        df["group"], df["group_name"] = groups, group_names
        self._df = df

    def rename(self, names:(list or dict)):
        """
        Replaces the generic Group0,... in the "group_name" column,
        using a list (set) or dict specifying new group names. 
        Group names only need to be specified once, and are applied to all replicate entries.
        """
        # get new group names based on list (index) or dict (key)
        if isinstance(names, (list, tuple, set)):
            new_names = self._rename_per_index(names)       
        elif isinstance(names, dict):
            new_names = self._rename_per_key(names)
        else:
            aw.HardWarning("Samples:no_groupname_assignment", names = names)

        # update "group_name"
        self._df["group_name"] = new_names
        self._renamed = True

    def ignore(self, entries:tuple):
        """
        Removed lines based on index from the dataframe.
        This is useful when removing corrupted data entries.
        """
        self._df = self._df.drop(index = list(entries)).reset_index()

    def _rename_per_key(self, names):
        """
        Generates new name list based on current names in "group_name" and uses string.replace()
        to update groupnames, based on key (old name) : value (new name) indexing. 
        Before applying it checks if all groups are covered by new names
        """
        current_names = len(aux.sorted_set(self._df["group_name"]))
        all_groups_covered = len(names) == current_names
        if all_groups_covered:
            current_names = list(self._df["group_name"])
            new_names = "$".join(current_names)
            for old_name, new_name in names.items():
                new_names = new_names.replace(old_name, new_name)
            new_names = new_names.split("$")       
            return new_names
        else:
            aw.HardWarning("Samples:groupnames_dont_colver", current_groups = current_names)

    def _rename_per_index(self, names):
        """
        Generates new name list based on current names in "group_names" and uses string.replace()
        to update groupnames to new names based on index (using a the order 
        of groups as is currently present in "group_name"). 
        """
        current_names = len(aux.sorted_set(self._df["group_name"]))
        all_groups_covered = len(names) == current_names
        if all_groups_covered:
            current_names = list(self._df["group_name"])
            names = list(names)
            new_names = "$".join(current_names)

            current_names_set = aux.sorted_set(current_names)

            for old_name, new_name in zip(current_names_set, names):
                new_names = new_names.replace(old_name, new_name)
            new_names = new_names.split("$")
            return new_names
        else:
            aw.HardWarning("Samples:groupnames_dont_colver", current_groups = current_names)


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

class SampleReader(Samples):
    """
    This class reads in a sample file and handles the 
    stored raw data in a pandas dataframe
    """
    def __init__(self):
        super().__init__()
        self._replicates = None
        self._names = None
        self._Reader = None
        self._Samples = None

    def replicates(self, replicates:(int or tuple)):
        """
        Set the replicates to group samples.
        This can be either an int (all groups have same number of replicates)
        or a tuple that specifies the number of replicates for every group separately.
        """
        self._replicates = replicates

    def names(self, names:(list or dict)):
        """
        Set names for replicates groups. These may be either specified as a list (requires same length as number of groups, assignment per index), 
        or as a dictionary (also requires new_name for every group), assignment per key. 
        """
        self._names = names
        
    def read(self, filename):
        """
        Reads one raw datafile (csv) 
        """
        self._Reader = Reader(filename)
        self._Reader.id(aux.fileID(filename))

        self._Samples = Samples(self._Reader)
        self._Samples.adopt_id(self._Reader)

        if self._replicates is not None:
            self._Samples.replicates(self._replicates)
            self._Samples.group()
        
        if self._names is not None:
            self._Samples.rename(self._names)

        return self._Samples


class Results(aux._ID):
    """
    This class handles a pandas dataframe for the results from qpcr.Analyser
    """
    def __init__(self):
        super().__init__()
        self._results = {"group" : [], "dCt" : []}
        self._df = None
        self._Samples = None
        self._stats_results = {"group" : [], "sample" : [], "mean" : [], "stdev" : [], "median" : []}
        self._stats_df = None

    def link(self, Samples:Samples):
        """
        Links an instance of Samples to be used for group_names
        And adds a group_name column to the dataframe
        """
        self._Samples = Samples
        if self.is_empty():
            self._results = {"group" : [], "group_name" : [], "dCt" : []}
        else:
            aw.SoftWarning("Results:cannot_link")
    
    def is_named(self):
        """
        Returns True if group_name column is present
        """
        return "group_name" in self._results.keys()
    
    def names(self, as_set = False):
        """
        Returns the linked group_names (only works if Samples have been linked!)
        """
        if self._Samples is not None:
            return self._Samples.names(as_set)
        return None

    def get(self):
        """
        Returns the results dataframe
        """
        return self._df

    def is_empty(self):
        """
        Checks if any results have been stored so far (True if not)
        """
        return self._df is None

    def add_names(self, Sample):
        """
        Adds group_names to results dataframe using the group names
        specified by a qpcr.Sample object.
        """
        names = Sample.names(as_set = False)
        length = len(self._results["group"])
        if len(names) == length:
            self._results["group_name"] = names
            self._df = pd.DataFrame(self._results)
        else: 
            aw.SoftWarning("Results:cannot_add_names", length = length, length_names = len(names))
    
    def add(self, dCt:tuple, group:int):
        """
        Adds a new set of deltaCt results for a given group. 
        This is used during DeltaCt computation
        """
        entries = len(dCt)
        self._results["group"].extend([group] * entries)
        self._results["dCt"].extend(dCt)
        if self._Samples is not None:
            names = self._Samples.names()
            self._results["group_name"].extend([names[group]] * entries)
        
        self._df = pd.DataFrame(self._results)

    def add_column(self, column, name):
        """
        Adds an entire computed results column to the dataframe
        """
        length = len(self._results["group"])
        if len(column) == length:
            self._results[name] = column
            self._df = pd.DataFrame(self._results)
        else:
            aw.HardWarning("Results:cannot_add_column", length = length, length_column = len(column))
    
    def copy_groups(self, Results):
        """
        Copies the groups from another Results instance dataframe, 
        but drops the dCt column. This is used when doing normalisation of datasets,
        where no individual dCt columns are necessary...
        """
        self._results["group"] = list(Results.get()["group"])
        del self._results["dCt"]
        self._df = pd.DataFrame(self._results)

    def merge(self, *Results):
        """
        Merges any number of other Results Instances into this one
        """
        new_df = self._df
        for R in Results: 
            R_df = R.get()
            # we merge the dataframes based on their groups, and add the instance id as identifier
            new_df = pd.merge(new_df, R_df["dCt"], 
                                right_index = True, left_index = True, 
                                suffixes = [f"_{self.id()}", f"_{R.id()}"]
                            )
        self._df = new_df # update self._df and also self._results!
        self._results = self._df.to_dict(orient="list")

    def stats(self, recompute = False) -> pd.DataFrame:
        """
        Returns a new dataframe containing Mean, Median, and 
        StDev of all replicate groups, for all samples.
        """
        # if stats_df is already present, return but sorted according to samples, not groups (nicer for user to inspect)
        if self._stats_df is not None and not recompute:
            return self._stats_df.sort_values("sample")
        
        # get groups and samples 
        groups = aux.sorted_set(list(self._df["group"]))
        samples = [c for c in self._df.columns if c not in ["group", "group_name"]]
        
        # compute stats for all samples per group
        for group in groups:
            group_subset = self._df.query(f"group == {group}")
            
            median = self._stat_var(group_subset, stats.median)
            mean = self._stat_var(group_subset, stats.mean)
            stdv = self._stat_var(group_subset, stats.stdev)
            self._add_stats(samples, group, median, mean, stdv)
            
        # add group names if present
        if self.is_named():
            self._add_stats_names(samples)

        self._stats_df = pd.DataFrame(self._stats_results)
        return self._stats_df.sort_values("sample")

    def save(self, path, df = True, stats = True):
        """
        Saves a csv file for each specified type of results.
        Path has to be a directory, if both df and stats are True.
        """
        if df and stats and not os.path.isdir(path):
            aw.HardWarning("Results:save_need_dir")

        print(self._df)

        if df:
            self._save_single(path, self._df, "_df")
        if stats:
            if self._stats_df is None:
                self.stats()
            self._save_single(path, self._stats_df, "_stats")

    def _save_single(self, path, src, suffix=""):
        """
        Saves either self._df or self._stats_df to a csv file based on a path
        (path can be either filename or directory)
        """
        if not os.path.isdir(path):
            src.to_csv(path)
        else:
            src.to_csv(os.path.join(path, f"rel_{self.id()}{suffix}.csv"))
        

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
        self._stats_results["sample"].extend(samples)
        self._stats_results["median"].extend(median)
        self._stats_results["mean"].extend(mean)
        self._stats_results["stdev"].extend(stdv)


    def _stat_var(self, group_subset, func):
        """
        Performs a function (like mean or stdv) over all rows
        and returns the result as list with a float for each column in the df
        any function can be passed as long as it works with an iterable
        """
        # ignore group and group_name columns
        ignore = ["group", "group_name"]
        all_cols = [g for g in group_subset.columns if g not in ignore]
        tmp = group_subset[all_cols]
        # compute stats based on func
        stats = [func(tmp[col]) for col in tmp]
        return stats
        

class Analyser(aux._ID):
    """
    This class performs Single Delta CT (normalisation within dataset) 
    Delta Delta CT (normalisation using second dataset), is handled by Normaliser()!
    """
    def __init__(self, Samples:Samples = None):
        super().__init__()
        self._Samples = Samples
        self._Results = Results()

        # default settings
        self._anchor = "first"
        self._efficiency = 2
        self._deltaCt_function = self._deltaCt_function(exp = True)

        if self._Samples is not None: 
            self._Results.adopt_id(Samples)
            if self._Samples.is_named(): # link Sample to add group_name if groups are named...
                self._Results.link(self._Samples)
    
    def get(self):
        """
        Returns the Results object that contains the results
        """
        return self._Results

    def has_results(self):
        return not self._Results.is_empty()

    def link(self, Samples:Samples, force = False, silent = False):
        """
        Links a qpcr.Samples object to the Analyser
        Note: If there are any precomputed results, no new data will be linked, unless force=True is called. 
        The user is notified if results are present and how to proceed. 
        """
        empty = self._Results.is_empty()
        dont_overwrite = not empty and not force
        if not dont_overwrite:
            self._Samples = Samples
            self.adopt_id(self._Samples)
            self._Results = Results()
            self._Results.adopt_id(self._Samples)
            if self._Samples.is_named():
                self._Results.link(self._Samples)
            
        if not silent:
            # notify the user of changes to the Analyser data and results
            if not dont_overwrite and not empty:
                aw.SoftWarning("Analyser:newlinked")
            elif dont_overwrite and not empty:
                aw.SoftWarning("Analyser:not_newlinked")

    def pipe(self, Samples:Samples, **kwargs) -> Results:
        """
        A quick one-step implementation of link + DeltaCt (silently overwrite any previous results)! 
        Returns a Results instance. Pipe forwards any **kwargs to DeltaCt.
        """
        self.link(Samples, force=True, silent=True)
        self.DeltaCt(**kwargs)
        return self.get()

    def efficiency(self, e:float = None):
        """
        Sets an efficiency factor for externally calculated qPCR amplification efficiency.
        By default efficiency = 2 is assumed.
        """
        if isinstance(e, (int, float)):
            self._efficiency = float(e)
        elif e is None: 
            return self._efficiency

    def anchor(self, anchor):
        """
        Sets the anchor for DeltaCt
        This can be either 
        - "first" (default, very first dataset entry)
        - "grouped" (first entry for each replicate group)
        - a specified numeric value
        """
        self._anchor = anchor

    def func(self, f:(str or function)):
        """
        Sets the function to be used for DeltaCt (optional)
        Available inputs are 
        "exponential" --> uses efficiency^(-(s-r)) # default efficiency = 2
        "linear" --> uses s-r
        any defined function that accepts an anchor (1st!) and sample (2nd!) 
        numeric value each, alongside any kwargs (will be forwarded from .DeltaCt()...)
        """
        if f in ["exponential", "linear"]:
            f = True if f == "exponential" else False
            self._deltaCt_function = _deltaCt_function(f)
        elif type(f) == type(aux.fileID):
            self._deltaCt_function = f
        else:
            aw.HardWarning("Analyser:cannot_set_func", func = f)

    def DeltaCt(self, **kwargs):
        """
        Calculates DeltaCt for all groups within the samples.
        As anchor for normalisation either three options may be specified
        first   -  the very first row in the dataset
        grouped -  the first replicate within each group
        specified_value (some external numeric value that will be used directly).
        Any additional arguments that a custom deltaCt function will require may be passed
        via the kwargs.
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
        # get set of sample groups and dataset
        groups = self._Samples.groups()
        df = self._Samples.get()
        for group in groups: 
            group_subset = df.query(f"group == {group}").reset_index()
            delta_cts = [deltaCt_function(anchor, i, **kwargs) for i in group_subset["Ct"]]
            self._Results.add(delta_cts, group)

    def _DeltaCt_grouped_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the first entry of each group as anchor
        """
        # get set of sample groups and dataset
        groups = self._Samples.groups()
        df = self._Samples.get()
        for group in groups: 
            group_subset = df.query(f"group == {group}").reset_index()
            anchor = group_subset["Ct"][0]
            delta_cts = [deltaCt_function(anchor, i, **kwargs) for i in group_subset["Ct"]]
            self._Results.add(delta_cts, group)

    def _DeltaCt_first_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the very first entry of the dataset as anchor
        """
        # get set of sample groups and dataset
        groups = self._Samples.groups()
        df = self._Samples.get()
        anchor = df["Ct"][0]
        for group in groups:
            group_subset = df.query(f"group == {group}").reset_index()
            delta_cts = [deltaCt_function(anchor, i, **kwargs) for i in group_subset["Ct"]]
            self._Results.add(delta_cts, group)
        return anchor

    def _exp_DCt(self, ref, sample, **kwargs):
        """
        Calculates deltaCt exponentially
        """
        factor = sample-ref 
        return self._efficiency **(-factor)

    def _simple_DCt(self, ref, sample, **kwargs):
        """
        Calculates deltaCt linearly
        """
        return sample-ref

    def _deltaCt_function(self, exp):
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
    This class handles normalisation of two (or more) datasets against, 
    using one as normaliser.
    This requires that all have been Analysed in the same way before!
    """
    def __init__(self):
        super().__init__()
        self._Normalisers = []
        self._Samples = []
        self._Results = Results()
        self._normaliser = None
        self._prep_func = self._average
        self._norm_func = self._divide_by_normaliser

    def get(self):
        """
        Returns the normalised dataframe
        """
        return self._Results
    
    def link(self, samples:(list or tuple) = None, normalisers:(list or tuple) = None):
        """
        Links either normalisers or any number of samples
        """
        self._link_normaliser(normalisers)
        self._link_samples(samples)
    
    def prep_func(self, f = None):
        """
        Sets any defined function for combined normaliser preprocessing...
        The function may accept one list of qpcr.Results instances, and must return 
        one list (iterable) of the same length as entries within the qpcr.Results dataframes.
        """
        # isinstance(f, function) didn't work...
        if type(f) == type(aux.fileID):
            self._prep_func = f
        elif f is None:
            return f
        else: 
            aw.HardWarning("Normaliser:cannot_set_prep_func", func = f)

    def norm_func(self, f = None):
        """
        Sets any defined function to perform normalisation of samples against normalisers.
        The function may accept one numeric entry for a sample and a normaliser, and must return 
        a numeric value. 
        Be default s/n is used...
        """
        if type(f) == type(aux.fileID):
            self._norm_func = f
        elif f is None:
            return f
        else: 
            aw.HardWarning("Normaliser:cannot_set_norm_func", func = f)

    def normalise(self):
        """
        Normalises all linked samples against the normaliser set. 
        Stores the results in a new Results 
        """
        if self._normaliser is None: 
            self._preprocess_normalisers()

        if self._Samples == [] or self._normaliser is None:
            raise Warning("Normalisation cannot be performed, as either samples are missing or no normaliser has been specified yet, or could not be processed!")

        # get normaliser dataframe
        normaliser = self._normaliser.get()

        # setup groups for _Results
        self._Results.copy_groups(self._Samples[0])

        # combine normalised samples into unified dataframe
        for S in self._Samples:
            S_df = S.get()
            column_name = f"{S.id()}_rel_{self._normaliser.id()}"
            normalised = [self._norm_func(s, n) for s,n in zip(S_df["dCt"], normaliser["dCt_combined"])]
            self._Results.add_column(normalised, column_name)

        # forward group_name column if present 
        if self._Samples[0].is_named():
            self._Results.add_names(self._Samples[0])

    def _divide_by_normaliser(self, s, n):
        """
        Performs normalisation of sample s against normaliser n
        (default _norm_func)
        """
        return s / n

    def _link_samples(self, samples):
        """
        Links any provided samples and checks their datatype in the process...
        """
        if samples is not None:
            for sample in samples: 
                if isinstance(sample, Results):
                    self._Samples.append(sample)
                elif isinstance(sample, Analyser) and sample.has_results():
                    self._Samples.append(sample.get())
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
        try: 
            combined = Results() # setup new dataframe for combined normalisers, intialise with first id
            combined.adopt_id(self._Normalisers[0])
            combined.copy_groups(self._Normalisers[0])
            combined.merge(*self._Normalisers)
            tmp_df = self._prep_func(combined)
            combined.add_column(tmp_df, "dCt_combined")
            self._normaliser = combined  
            if len(self._Normalisers) > 1:
                self._update_combined_id()
            
            # forward combined_id to self and _Results 
            self.adopt_id(self._normaliser)
            self._Results.adopt_id(self._normaliser)

        except: 
            pass

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
        generates a list of their per-row means
        (default preprocess_normalisers function)
        """
        tmp = combined.get()
        tmp_df = tmp.drop(["group"], axis = 1)
        tmp_df = tmp_df.mean(axis = 1)
        return list(tmp_df)


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

        # samples = Samples(reader)
        
        sample = reader.read(file)
        # sample.ignore((0,1,3,4))

        # analyser.link(sample, force=True, silent = False)
        # analyser.DeltaCt()
        # res = analyser.get()

        res = analyser.pipe(sample)
        # print(res)
        analysers.append(res)

    normaliser = Normaliser()
    normaliser.link(normalisers = analysers[:2])
    normaliser.link(samples = analysers[2:])

    normaliser.normalise()
    
    result = normaliser.get()

    print(result.get())
    
    # result.save("..")
    
    #result.add_names(samples)

    print(result.stats())

    # a = Reader("Example Data/28S.csv")
    # x = Reader("Example Data/28S.csv")
    # a.id("28S")
    # x.id("also28S")
    # b = Samples(a)
    # y = Samples(x)
    # b.replicates(6)
    # y.replicates(6)
    # b.group()
    # y.group()
    # b.rename(["wt-", "wt+", "ko-", "ko+"])
    # y.rename(["wt-", "wt+", "ko-", "ko+"])

    # c = Analyser(b)
    # z = Analyser(y)
    # c.DeltaCt(anchor = "first")
    # z.DeltaCt(anchor = "first")
    # d = c.get()
    # g = z.get()

    # print(d.id())

    # e = Normaliser()
    # e.link(normaliser = d)
    # e.link(g, c)
    # e.normalise()
    # r = e.get()
    # print(r.get())

    exit(0)


def open_csv_file(filename, export="dict"):
    """
    This function opens a given input csv file. To function properly, the file is required to have 
    a first column of sample names and a second column with Ct values. 
    The first line will be read as Column Titles and be cropped during processing.
    """

    with open(filename, "r") as openfile:
        contents = openfile.read()

    open_args = dict(header = 0, names = ["Sample", "Mean Ct"])
    
    if ";" in contents:
        contents = pd.read_csv(filename, sep = ";", **open_args )
    else: 
        contents = pd.read_csv(filename, **open_args )
    
    if export == "dict":
        contents = contents.to_dict(orient = "list")

    return contents

def _define_replicates(sample_dict, column_name, replicates):
    """
    Auxiliary Function that creates groups of indices to be used by group_samples for replicate assignment.
    """
    reps, equal_reps = _prep_replicates(replicates)
    
    samples = sample_dict[column_name]
    intervals = _make_intervals(reps, equal_reps, samples)
    if intervals is None:
        return None

    replicate_samples = []
    if column_name == "Mean Ct": # the numeric column
        for i in intervals:
            tmp = [float(samples[j]) for j in i]
            replicate_samples.append(tmp)
    else:
         for i in intervals:
            tmp = [str(samples[j]) for j in i]
            replicate_samples.append(tmp)
    
    return replicate_samples

def _make_intervals(reps, equal_reps, samples):
    """
    Makes equal or unequal intervals for samples according to reps and 
    checks if intervals cover all provided sample data entries.
    Returns None if not all samples are covered...
    """
    if equal_reps == True:
        intervals = _equal_intervals(samples, reps)
    else:
        intervals = _unequal_intervals(samples, reps)
    return intervals

def _unequal_intervals(samples, reps):
    if sum(reps) != len(samples): 
        print(f"You are not providing replicate annotation for each sample!\n{len(samples)-sum(reps)} samples are missing intervals!")
        return None
    intervals = _create_unequal_intervals(reps)
    return intervals

def _equal_intervals(samples, reps):
    """
    Makes equal intervals according to reps specified, 
    and checks it they are appropriate for the given sample. 
    Will return none if not...
    """
    intervals = _create_equal_intervals(len(samples), reps)
    if _vet_equal_intervals(intervals, samples) == False:
        return None
    return intervals

def _vet_equal_intervals(intervals, samples):
    #test out if equal intervals include all samples
    verdict = True
    ints = []
    for i in intervals: ints.extend(i)
    if len(ints) != len(samples):
        print(f"Warning: With the current replicate settings you are loosing samples!\nYou specify {len(intervals)} intervals for {len(ints)} entries but data has {len(samples)} entries...")
        verdict = False
    return verdict 


def _prep_replicates(replicates):
    """Establishes what kind of replicates are passed"""
    try: 
        reps = int(replicates)
        equal_reps = True
    except TypeError:
        reps = list(replicates)
        equal_reps = False
    return reps,equal_reps

#interval auxiliaries for replicate finding

#creates equalually spaces index groups
def _create_equal_intervals(length, reps):
    indices = []
    i = 0
    while i + reps <= length:
        indices.append(list(range(i, i+reps)))
        i+= reps
    return indices

#creates index groups based on the entry list provided in reps
def _create_unequal_intervals(reps):
    indices = []
    i = 0
    for r in reps:
        tmp = list(range(i, i+r))
        i+=r
        indices.append(tmp)
    return indices


def group_samples(sample_dict, replicates):
    """
    Second step in manual qpcr analysis. Groups samples into replicate groups as defined by replicates.
    """
    samples = _define_replicates(sample_dict, "Sample", replicates)
    cts = _define_replicates(sample_dict, "Mean Ct", replicates)
    
    groupings = {}
    for i in range(0,len(samples)):
        key = "Group {}".format(i+1)
        values = [samples[i],cts[i]]
        tmp = {key : values}
        groupings.update(tmp)
    return groupings



def Delta_Ct(grouped_dict, exp=True, anchor="first"):
    """
    Calculates Delta_Ct within one sample_dict. Requires group_samples to have been performed before use.
    exp allows to adjust to either 2^(d1-d2) (True) or only (d1-d2) (False).
    anchor sets the reference (d1). Possible are None (group internal, default), "first" (very first entry), or any specified numeric value. 
    """
    keys = list(grouped_dict.keys())
    new_dict = {}
    dCt = _deltaCt_function(exp)
        
    if anchor is None:
        _grouped_anchored(grouped_dict, keys, new_dict, dCt)
    elif anchor == "first":
        _first_anchored(grouped_dict, keys, new_dict, dCt)
    else:
        _specified_anchored(grouped_dict, anchor, keys, new_dict, dCt)
    return new_dict

def _deltaCt_function(exp):
    """
    Returns the function to be used for DeltaCt based on 
    whether or not exponential shall be used.
    """
    if exp == True:
        dCt = _exp_DCt
    else:
        dCt = _simple_DCt
    return dCt

# functions to calculate delta ct based on the anchor specified 

def _specified_anchored(grouped_dict, anchor, keys, new_dict, dCt):
    first = anchor
    for k in keys:
        cts = grouped_dict[k][1]
        tmp = [dCt(first, i) for i in cts]
        new_dict.update({k : tmp})

def _first_anchored(grouped_dict, keys, new_dict, dCt):
    first = grouped_dict[keys[0]][1][0]
    for k in keys:
        cts = grouped_dict[k][1]
        tmp = [dCt(first, i) for i in cts]
        new_dict.update({k : tmp})

def _grouped_anchored(grouped_dict, keys, new_dict, dCt):
    for k in keys:
        cts = grouped_dict[k][1]
        first = cts[0]
        tmp = [dCt(first, i) for i in cts]
        new_dict.update({k : tmp})

def _exp_DCt(ref, sample):
    factor = sample-ref 
    return 2**(-factor)

def _simple_DCt(ref, sample):
    return sample-ref


def normalise(normaliser:dict, sample:dict, no_head=True):
    """
    This function normalises a sample_dict against a corresponding normaliser. 
    Requires that both have undergone Delta_Ct before. 
    This function corresponds to the second Delta in DeltaDelta CT analysis of qPCR data.
    no_head=False would add "Legend" : "DDCT" as very first entry to the returned dictionary.
    """
    if no_head == False:
        normalised = {"Legend" : "DDCT"}
    else: 
        normalised = {}
    keys_n = normaliser.keys()
    keys_s = sample.keys()
    if keys_n != keys_s:
        print("Normaliser and Samples do not share the same group names!")
        return None
    for k in keys_s:
        normals = normaliser[k]
        samples = sample[k]
        if len(normals) != len(samples):
            print("Normaliser and Samples do not have the same length!")
            return None
        
        temp = []
        for i in range(0, len(normals)):
            n = normals[i]
            s = samples[i]
            rel = s/n
            temp.append(rel)
        
        normalised.update({k : temp})
    return normalised

#if several normalisers should be used they first have to be pre-processed to get their Delta_Ct values before they can be combined using combine_normalisers
def preprocess_normalisers(normalisers:list, replicates, run_names=None, group_names=None, anchor=None, return_type="list"):
    """
    This function takes in a list of input csv files of normalisers and automatically performes grouping, Delta_Ct and renaming.
    Note that the normalisers have to be named and group_names must correspond to group_names later used for the actual samples!
    """
    normalisers_dict = {}
    ndx = 0
    if run_names is None: 
        names = _ddCt_generate_run_names(normalisers)
        run_names = names

    for norm in normalisers:
        name = run_names[ndx]
        tmp = open_csv_file(norm)
        tmp = group_samples(tmp, replicates=replicates)
        if group_names is not None:
            tmp = rename_groups(tmp, new_names=group_names)
        tmp = Delta_Ct(tmp, anchor=anchor)
        tmp = {name : tmp}
        normalisers_dict.update(tmp)
        ndx +=1
    if return_type == "list":
        temp = []
        for i in run_names:
            temp.append(normalisers_dict[i])
        normalisers_dict = list(temp)
    return normalisers_dict

#if several normalisers should be averaged, use combine_normalisers
def combine_normalisers(normalisers:list):
    combined_normaliser = {}
    keys = list(normalisers[0].keys()) #get all groupings
    for k in keys:
        temp = []
        for norm in normalisers:
            length = len(norm[k])
            assert isinstance(norm[k][0], float), "norm[k][0] is not a simple number but: {} ... \nAre you sure you already performed qpcr.Delta_Ct on your normaliser?".format(norm[k][0])
            for i in norm[k]:
                temp.append(i)
        temp = {k : [stat.mean(temp) for i in range(length)]} #we must conserve dimensionality...
        combined_normaliser.update(temp)
    return combined_normaliser


def rename_groups(sample_dict, new_names:list):
    """
    Sample names from the original csv file are not imported! 
    Samples will be only named according to their groups. Group names can be set using this function. If no names are provided, "Group 1", "Group 2" etc. will be used by default.
    rename_groups can be used on any dict at any stage in the process. 
    """
    keys = sample_dict.keys()
    new_dict = {}
    i = 0
    for k in keys:
        new_name = new_names[i]
        tmp = { new_name : sample_dict[k] }
        new_dict.update(tmp)
        i +=1
    return new_dict


def get_stats(deltaCT_dict, export = ["avg", "stdv"]):
    """
    By default Delta_Ct (and normalise afterward) will work with individual replicate values, and return a dict of the same dimensions.
    To concentrate group information from many replicates directly to average and stdev, use get_stats.
    This will return a dict where each group assigns a list of two values, the first being average, the second being stdev.
    """
    keys = deltaCT_dict.keys()
    new_dict = {}
    legend = []
    if "avg" in export: legend.append("Avg")
    if "stdv" in export: legend.append("Stdev")
    if "med" in export: legend.append("Med")
    new_dict.update({"Legend" : legend})
    for k in keys:
        dcts = deltaCT_dict[k]
        avg = stat.mean(dcts)
        std = stat.stdev(dcts)
        med = stat.median(dcts)
        
        tmp = []
        if "avg" in export:
            tmp.append(avg)
        if "stdv" in export:
            tmp1 = [i for i in tmp]
            tmp1.append(std)
            tmp = tmp1
        if "med" in export:
            tmp1 = [i for i in tmp]
            tmp1.append(med)
            tmp = tmp1

        tmp = {
            k : tmp
        }
        new_dict.update(tmp)
    return new_dict


def export_to_csv(data, filename, transpose=False):
    """
    This function writes a csv file with the results generated. It may be called at any stage in the analysis process to save intermediary data.
    transpose will transpose columns and rows.
    """
    export_data = pd.DataFrame(data)
    if transpose == True:
        export_data = export_data.transpose()
    export_data.to_csv(filename)
    #and now re-open to remove the first line that only contains 0 1 and empty commas
    contents = []
    with open(filename, "r") as openfile:
        for i in openfile:
            contents.append(i)
    contents = contents[1:]
    with open(filename, "w") as openfile:
        if "Legend" not in contents[0]: 
            openfile.write("Sample,Value\n")
        for i in contents: 
            openfile.write(i)

#export grouped raw CT values to plot if desired
def export_raw_data(filename, replicates, group_names=None, export_location=None):
    """
    To also group raw data, use this function. This function reads original csv files, groups them and exports them back to the same location, creating a new file with _raw.csv as ending.
    """
    contents_dict = open_csv_file(filename)
    contents_dict = group_samples(contents_dict, replicates=replicates)
    if group_names is not None:
        contents_dict = rename_groups(contents_dict, group_names)
    tmp = {}
    for k in list(contents_dict.keys()):
        tmp.update({k : [str(i) for i in contents_dict[k][1]]  })
    savefile_dict = tmp

    if export_location is None:
        new_file = "{}_raw.csv".format(filename.replace(".csv", ""))
    else:
        new_file = export_location
    
    with open(new_file, "w") as openfile:
        for k in list(savefile_dict.keys()):
            newstr = "{},{}\n".format(k, ",".join(savefile_dict[k]))
            openfile.write(newstr)


#if one wants to revisit already computed results
def load_results(filename, mode="individual", *args, **kwargs):
    """
    This function opens pre-computed results of the ops.biotools.qpcr Module from the generated csv files. 
    It supports two modes: 
    mode = "individual" (default) where filename specifies the file to be opened. It returns a dictionary containing the grouped computed values (replciates, or avg, stdev).
    mode = "pairs" allows users to load an entire set of samples and normalisers that were previously computed, to then be used by ops.biotools.qpcr.Analysis.normalise_pairs. It returns two separate dictionaries each containing the set of sample files it was given by kwargs parameters samples:list, normalisers:list. Optionally, names may be additionally assigned to samples and normalisers using sample_names:list and norm_names:list. 

    To facilitate working with replicate assays (i.e. same qPCR assay normalised against different normalisers separately), samples / normalisers support only partial naming and need no full filepath to function. Like this multiple analysis result with e.g. "HNRNPL NMD" in their name will all be loaded. 
    """
    if mode == "individual":

        export_dict = _load_individual_results(filename)

    elif mode == "pairs":
        #get necessary kwargs for _load_pair_results
        location, samples, normalisers, sample_names, norm_names = _check_loadpair_kwargs(filename=filename, **kwargs)

        #load samples into dict
        export_dict = _load_pair_results(location=location, 
                                            samples=samples, 
                                            normalisers=normalisers, 
                                            sample_names=sample_names, 
                                            norm_names=norm_names)
    
    elif mode == "multiple":
        assert isinstance(filename, list), "When using mode='multiple', please provide a list of valid filenames for the filename argument"
        
        keys = oo.from_kwargs("keys", None, **kwargs)
        export_dict = _bulk_load(filename, keys=keys)
           
    else:
        print("Please, use either mode='individual', mode='multiple', or mode = 'pairs' to load your results")
        return None

    return export_dict

def _load_individual_results(filename):
    contents = []
    with open(filename, "r") as openfile:
        contents = openfile.read().split("\n")

    contents = [i.split(",") for i in contents]
    conditions = [i[0] for i in contents]
    values = [i[1:] for i in contents]
    values = values[1:] #crop the first labelling
    values = [[float(j) for j in i] for i in values] #get the numeric values back
    
    #now pack everything back into a dict
    export_dict = {}
    idx = 0
    for i in conditions[1:]:
        tmp = {i : values[idx]}
        export_dict.update(tmp)
        idx +=1
    return export_dict

def _check_loadpair_kwargs(filename, samples, **kwargs):
    pair_args = oo.get_kwargs(_load_pair_results, **kwargs)
    location = filename
    
    normalisers = oo.from_kwargs(("norms", "normalisers"), None, pair_args)
    sample_names = oo.from_kwargs("sample_names", None, pair_args)
    norm_names = oo.from_kwargs("norm_names", None, pair_args)

    if normalisers is None:
        raise NameError("No normalisers were found in kwargs! Make sure to provide normalisers using either 'norms' or 'normalisers' arguments!")
    return location, samples, normalisers, sample_names, norm_names

def _load_pair_results(location, samples, normalisers, sample_names = None, norm_names = None):
    if isinstance(location, (tuple, list)) and len(location) == 2:
        sample_loc = location[0]
        norm_loc = location[1]
    else:
        sample_loc = location
        norm_loc = location
    
    
    sample_items = fc.item_finder(path=sample_loc, filter=[".csv"])
    sample_items.sort()
    sample_list = fc.filter_files(items_list=sample_items, target_names=samples)
    
    norm_items = fc.item_finder(path=norm_loc, filter=[".csv"])
    norm_items.sort()
    norm_list = fc.filter_files(items_list=norm_items, target_names=normalisers)
    
    sample_dict = _bulk_load(sample_list, keys=sample_names)
    norm_dict = _bulk_load(norm_list, keys=norm_names)

    return sample_dict, norm_dict

def _bulk_load(filenames, keys=None):
    return_dict = {}
    k = 0
    print("Loading Result files...")
    if keys is not None:
        assert len(keys) == len(filenames), "Keys and filenames lists are not of equal length!"
        for f in filenames:
            name = keys[k]
            print("""
            {}\t\t Src: {}
            """.format(name, f))
            tmp = load_results(f)
            return_dict.update({name : tmp})
            k+=1
    else:
        for f in filenames:
            name = "Assay {}".format(k)
            print("""
            {}\t\t Src: {}
            """.format(name, f))
            tmp = load_results(f)
            return_dict.update({name : tmp})
            k+=1
    return return_dict            

# generate run_names for qA.delta_deltaCt in case the user does not want to assign manually
def _ddCt_generate_run_names(data_files):
    target_names = [_remove_suffix(i) for i in data_files]
    return target_names

def _remove_suffix(file):
    if "/" in file:
        norm_name = file.split("/")
        norm_name = norm_name[-1]
    norm_name = norm_name.split(".")
    norm_name = norm_name[0]
    return norm_name


def help():
    from ops.biotools.qpcr import Analysis
    Analysis.help()
    
def Info():
    from ops.biotools.qpcr import Analysis
    Analysis.Info()

def Example():
    from ops.biotools.qpcr import Analysis
    Analysis.Example()



if __name__ == "__main__":
    test = "Example Data/28S.csv"

    a = open_csv_file(test)

    grouped = group_samples(a, 3)
    # print(grouped)

    loaded = load_results("./Example Data/HNRNPL_nmd_DeltaDelta_Ct.csv")
    print(loaded)
