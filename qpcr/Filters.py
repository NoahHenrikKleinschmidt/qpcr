"""
This submodule defines a number of filters that can be used to 
remove faulty reads prior to deltact analysis. 

Filters implemented include
- RangeFilter (filter out all raw Ct values that do not comply to a user-specified range, default +-1)
    - here the user will have the option of specifying an anchor (default = median of replicates)

- IQRFilter (filter out any outliers by IQR)
"""

import __init__ as qpcr
import pandas as pd
import numpy as np
import auxiliary.warnings as wa
import auxiliary as aux
import os 

class Filter(aux._ID):
    """
    The super filtering class that takes in a qpcr.Assay() object and updates its _df to a filtered version.
    """
    def __init__(self):
        super().__init__()
        self._Assay = None
        self._report_loc = None
        self._id = type(self).__name__
    

    def link(self, Assay:qpcr.Assay):
        """
        Link a qpcr.Assay to be filtered
        """
        self._Assay = Assay

    def pipe(self, Assay:qpcr.Assay, **kwargs):
        """
        A shortcut for link+filter
        It returns directly the filtered Assay object
        """
        self.link(Assay)
        return self.filter(**kwargs)

    def filter(self, **kwargs):
        """
        Applies the filter and returns an updates Assay object
        """
        if self._Assay is not None:
            return self._filter(**kwargs)
        else: 
            wa.HardWarning("Filter:no_assay")

    def report(self, filename):
        """
        Stores a report of any replicates that were filtered out
        """
        self._report_loc = filename

    def reset(self):
        """
        Resets the excluded indices
        """
        self._faulty_indices = []

    def _filter(self, **kwargs):
        """
        The actual filtering function that each FilterObject will define.
        """
        print("The actual filtering function that each FilterObject will define")
        # do stuff
        return self._Assay

    def _write_report(self, faulty_indices):
        """
        Generates a filtering report file
        """
        filename = "filter_" + self._Assay.id() + ".txt"
        filename = os.path.join(self._report_loc, filename)

        report_string = f"""
Filtering Report
Filter: {self._id}
Sample: {self._Assay.id()}
Found faulty Replicates: {len(faulty_indices)}
Found Indices: {faulty_indices}
        """
        report_string = report_string.strip()
        if os.path.exists(filename):
            with open(filename, "a") as f:
                f.write(report_string.replace("Filtering Report", "\n"))
        else:
            with open(filename, "w") as f:
                f.write(report_string)


class RangeFilter(Filter):
    """
    This class filters out any replicate that lie outside a user-specified range.
    Default are +- 1 or the replicate group median. 
    """
    def __init__(self):
        super().__init__()
        self._upper = 1
        self._lower = 1
        self._anchor = None
        

    def set_lim(self, lim = None, upper = None, lower = None):
        """
        Sets the range limits for inclusion range.
        lim will set symmetric upper and lower bounds, otherwise
        upper and lower can be directly specified.
        """
        if lim is not None: self._upper, self._lower = lim, lim
        if upper is not None: self._upper = upper
        if lower is not None: self._lower = lower

    def set_anchor(self, anchor):
        """
        Set the range anchor (center of inclusion range)
        Supported types for anchor are:
        - a numeric value (int or float)
        - an iterable of same length as groups
            - if a dict keys must be group indices (starting from 0)
        - a function that works with a pandas dataframe as stored by qpcr.Assay objects
            - function must return a single numeric value for anchor
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
        
        # self._faulty_indices.extend(faulty_indices)

        # exclude faulty entries
        if len(faulty_indices) > 0:
            self._Assay.ignore(faulty_indices)

        if self._report_loc is not None: 
            self._write_report(faulty_indices)

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


if __name__ == "__main__":
    
    files = ["Example Data/28S.csv", "Example Data/actin.csv", "Example Data/HNRNPL_nmd.csv", "Example Data/HNRNPL_prot.csv"]
    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    analysers = []

    reader = qpcr.SampleReader()
    reader.replicates(6)
    reader.names(groupnames)

    analyser = qpcr.Analyser()
    analyser.anchor("first")

    range_filter = RangeFilter()
    range_filter.report(".")

    for file in files: 
        
        sample = reader.read(file)
        
        # here comes in the filter...
        sample = range_filter.pipe(sample)
        # print(sample.get())

        res = analyser.pipe(sample)
        analysers.append(res)

    # for a in analysers: print(a.id(), "\n", a.get())

    normaliser = qpcr.Normaliser()
    normaliser.link(normalisers = analysers[:2])
    normaliser.link(samples = analysers[2:])

    normaliser.normalise()
    
    result = normaliser.get()

    print(result.get())
    
    # result.save("..")
    
    # result.add_names(samples)

    print(result.stats())

    exit(0)