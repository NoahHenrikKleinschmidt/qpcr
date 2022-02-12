"""
This module contains classes designed to work with irregularly structured datafiles.
It provides `Parsers` that are able to extract the replicate and Ct values as pandas DataFrames
from irregular `csv` and `excel` files.


## Working with "irregular files"
----

Any datafile that does not only consist of a replicate identifer and Ct column is called "irregular". 
In fact, most excel sheets or csv exports from qPCR machines are actually irregular as they often contain some
information about the run, and melting curve data, and so forth. Such data is not relevant or of interest to `qpcr`, however, 
so we have extract the data we are intersted in. This is the job of the `qpcr.Parsers`. They read in an irregular datafile and use 
a guided-parsing approach to find the relevant sections within the datafile. There are essentially two ways how they can do this. 

### "Finding" relevant datasets through `assay_patterns`
The Parsers are quipped with a method called `_CORE_Parser.find_assays` which locates assays (or just "datasets") within the datafile
using `regex`. Of course, in order to do that it needs to know the patterns it's supposed to look for. Some patterns are already pre-specified
in the `qpcr.Parsers.assay_patterns` dictionary and can simply be specified using their key. If your own pattern is not yet pre-defined, 
[post an issue on github](https://github.com/NoahHenrikKleinschmidt/qpcr/issues) and supply some samples of how your assays usually appear in your datafiles 
alongside with the name of the machine that produces your datafiles.
Of course, you can also manually specify your own regex pattern. The only constraint is that is _must_ have *one* capturing group for the assay name.

Once the assays in your datafile are identified, the data columns belonging to them are searched for. The constraint here is that they must start either
in the row exactly below the assay headers or have at most one single row in between them. Anything else is no good! The data columns _must_ be labelled 
(i.e. have a header). By default `Name` and `Ct` are assumed as data column labels / headers, but this can be changed. 

If your datafiles contain multiple datasets / assays, the `qpcr.Parsers` will be able to extract all of them and store their extracted data. 


## Working with `multi-assay` files
----

While working with datafiles that contain multiple assays you have essentially two options. First, you are only interested in one of the assays in your file,
or you want to work with _all_ assays in your file at the same time. There is two different methods to deal with either of these situations.

### Getting a single assay from a multi-assay file
Parsers are actually implemented in the core of the `qpcr.Reader` but if you use a multi-assay file using the `qpcr.Reader` 
you will have to specify which assay you are specifically interested in using the `assay` argument. So you can simply read 
a multi-assay file as you would a regular file using the `qpcr.SampleReader` (or `qpcr.Reader` manually) and add the `assay` argument as kwarg to the `read` method.

### Getting all assays from a multi-assay file
Here it depends on what you want to do. Do you wish to extract the individual assays to make individual files from them or do you want to use them directly for an analysis?
Both options are possible. 

#### Making indivdual assay files from a multi-assay file
This is the core-business of the `qpcr.Parsers`. So you can simply set up a Parser, set a saving location using the 
Parser's `save_to` method and then `pipe` your file throuhg. All done at this point. Of course, you can also work with the dataframes directly
using the Parser's `get` method.

#### Using a multi-assay file directly for my analysis
So, you want to just feed in your one datafile and expect to get a table with your Delta-Delta-Ct values for all assays against all normalisers?!
Sure, no problem! You can easily read a multi-assay file using the `qpcr.MultiReader` which allows you (after setup) to simply `pipe` through your datafile
and returns immediately a list of your assays-of-interest and normaliser-assays. How does it know which is which though? 
That's the super-important topic of the next paragraph.

#### Decorators
A `decorator` technically is a function that wraps another function when coding. Well, that's not quite the case for the `qpcr` decorators but the idea is the same. 
Instead of wrapping functions the `qpcr` decorators decorate assays in a multi-assay file. What does "decorate" mean? It means they provide meta-data about the assay
in question. What does that mean? There are two implemented `qpcr` decorators `@qpcr:assay` and `@qpcr:normaliser` – now clear what they do? They are placed in the 
cell _exactly above_ the cell where the assay header is located (seriously, anything else won't do!) and tell the `qpcr.MultiReader` 
if a specific assay is an "assay-of-interest" or a "normaliser-assay".


`qpcr.Parsers` can identify assays either through _de novo_ finding using `regex` patterns *or* through decorators. To tell a Parser to use a specific decorator 
for finding assays you can specify the `decorator` argument in `pipe` or `parse`. To specify a decorator like this you only write `qpcr:assay` `qpcr:normaliser` (no @).
A third option is available only to `qpcr.Parsers` (not to the `qpcr.MultiReader`) and that is `qpcr:all` which will tell the Parser to extract _all_ decorated assays.
The implemented decorators are specified in the `qpcr.Parsers.decorators` dictionary and can be accessed via keys (you have seen the three keys already).


> ##### Two things of Note: 
> - When specifying a decorator any non-decorated assay will be *ignored*!
> - If you are using `excel` you may have to add a single tick `'` in front of your decorators.
"""

import qpcr
import qpcr._auxiliary as aux
import qpcr._auxiliary.defaults as defaults
import qpcr._auxiliary.warnings as aw
import pandas as pd
import numpy as np
import re
from io import StringIO
from copy import deepcopy
import os

__pdoc__ = {
    "_CORE_Parser" : True
}

# this is the dictionary where we store pre-defined 
# patterns of headers above assays within the datafiles
# important here is that they must specify a capturing group for the assay name.

assay_patterns = {
                    "all"           : r"([A-Za-z0-9.:, ()_\-/]+)",
                    "Rotor-Gene"    : r"Quantitative analysis of .+(?<=\()([A-Za-z0-9.:, _\-/]+)",
                }

decorators = {
                    "qpcr:all"          : "(@qpcr:|'@qpcr:)",
                    "qpcr:assay"        : "(@qpcr:assay\s{0,}|'@qpcr:assay\s{0,})",
                    "qpcr:normaliser"   : "(@qpcr:normaliser\s{0,}|'@qpcr:normaliser\s{0,})",     
                    "qpcr:group"        : "(@qpcr:group\s{0,}|'@qpcr:group\s{0,})",
                    "qpcr:column"       : "(@qpcr|'@qpcr)",

            }

plain_decorators = {
                    "qpcr:all"          : "@qpcr:",
                    "qpcr:assay"        : "@qpcr:assay",
                    "qpcr:normaliser"   : "@qpcr:normaliser",
                    "qpcr:group"        : "@qpcr:group",
                    "qpcr:column"       : "@qpcr",
            }

# get the standard column headers to use for the 
# replicate id and Ct column of the finished dataframes
standard_id_header = defaults.raw_col_names[0]
standard_ct_header = defaults.raw_col_names[1]

default_group_name = defaults.default_group_name
default_dataset_header = defaults.default_dataset_header

# set a dummy default value for any np.nan values
# in the column storing the assay headers
# we do this in case the "all" assay_pattern is used without decorators
# in this case any cell would be selected as "nan" also matches the pattern
dummy_blank = "$"

class _CORE_Parser:
    """
    This is the functional core for the irregular multi-assay file-reader classes.
    It handles the regex searching and numpy indexing of relevant column subsets of the datafiles.
    """
    def __init__(self):
        self._src = None
        self._pattern = None
        self._data = None

        # the found assays, these will be arrays/lists that store the indices, 

        self._assay_indices = None                  # indices of the assay identifiers
        self._assay_names = None                    # names of the assays

        self._assay_names_start_indices = None      # indices of the rep. id headers
        self._assay_names_end_indices = None        # indices of the last entry of the rep. id columns

        self._assay_ct_start_indices = None         # indices of the ct headers
        self._assay_ct_end_indices = None           # indices of the last entry of the ct columns

        # a dictionary to store all assay dataframes
        self._dfs = {}

        # setup the labels for replicate ids and ct value column headers
        self.labels()

        # we must specify a maximum allowed length for the assay names before hand 
        # (since we're using numpy arrays for storing the names, which require enough open slots to store the characters)
        self._max_assay_name_length = 20
    
        # a folder into which the new assay-split datafiles should be stored
        self._save_loc = None

        # set transpose option in case datasets are stored not on separate row ranges but separate column ranges
        self._transpose = False


        # set up a BigTable data range
        self._bigtable_range = None

        # setup a warning for decorators without patterns warning that it will set to default
        # we do this here in case the MultiSheetReader calls read+parse here multiple times
        # but we don't need to know that it defaults every time.
        self._WARNING_decorators_but_no_patterns = aw.SoftWarning("Parser:decorators_but_no_pattern", once = True)

    def prune(self):
        """
        Completely resets the Parser, clearing all data and preset-specifics such as the assay_pattern.
        """
        self.__init__()

    def clear(self):
        """
        Clears all datasets that were extracted.
        """
        self._dfs = {}

        self._assay_indices = None                  # indices of the assay identifiers
        self._assay_names = None                    # names of the assays

        self._assay_names_start_indices = None      # indices of the rep. id headers
        self._assay_names_end_indices = None        # indices of the last entry of the rep. id columns

        self._assay_ct_start_indices = None         # indices of the ct headers
        self._assay_ct_end_indices = None           # indices of the last entry of the ct columns

    def transpose(self):
        """
        Inverts the `col` index used by `qpcr.Parsers._CORE_Parser.find_assays` and `qpcr.Parsers._CORE_Parser.find_by_decorator`.
        By default the `col` refers to a column. After using `transpose` it will be interpreted as a `row`.
        Note
        ----
        This is method is dynamic, so repeated calling of `transpose` will keep reverting the interpretation from row to col, back to row, etc.
        """
        self._transpose = not self._transpose

    def save_to(self, location : str = None):
        """
        Sets the location into which the individual assay datafiles should be saved.
        Parameters
        ----------
        location : str
            The path to a directory where the newly generated assay datafiles shall be saved.
            If this directory does not yet exist, it will be automatically made.
        """
        if location is not None: 
            self._save_loc = location
            if not os.path.exists(self._save_loc):
                os.mkdir(self._save_loc)
        return self._save_loc
    
    def get(self, assay : str = None):
        """
        Parameters
        ----------
        assay : str
            The name of an assay found in the datafile. Available assays can be assessed using the `self.assays` method.
    
        Returns
        -------
        data : pd.DataFrame or dict
            Either a specific pandas dataframe of one of the assays (if an `assay` name was specified)
            or the entire dictionary of all found dataframes from all assays.
        """
        if assay is not None:
            data = self._dfs[assay]
        else: 
            data = self._dfs
        return data

    def save(self):
        """
        Saves the individual assays as separate csv files.
        This requires that a saving directory has been set using `self.save_to`.
        The files will simply be named according to the assay name (i.e. `ActinB.csv` for instance).
        """
        if self._save_loc is None:
            aw.SoftWarning("Parser:no_save_loc")
        else:
            for assay, df in self._dfs.items():
                assay_path = os.path.join(self.save_to(), f"{assay}.csv")
                df.to_csv(assay_path, index = False)

    def labels(self, id_label : str = "Name", ct_label : str = "Ct"):
        """
        Sets the headers for the relevant data columns for each assay within the datafile.

        Parameters
        ----------
        id_label : str
            The header above the column containing replicate identifiers. 
        
        ct_label : str
            The header above the column containing the replicates' Ct values.
        """
        self._id_label = id_label
        self._ct_label = ct_label

    def assays(self):
        """
        Returns
        -------
        names : list
            The names of the found assays of the datafile
        """        
        return list(self._dfs.keys())

    def assay_pattern(self, pattern : str = None, *flags):
        """
        Sets up a regex pattern defining the assay declarations within the datafile.

        Parameters
        ----------
        pattern : str
            A string containing either the key to a predefined pattern from the `assay_patterns` dictionary, 
            or directly regex pattern. 
            If a regex pattern is directly provided, that pattern must contain a capturing group
            for the assay name that can be extracted.
        *flags 
            Any additional flags to pass to `re.compile()` for the regex pattern

        Returns
        -------
        pattern : re.Pattern
            The currently used regex pattern to identify assays within the datafile.
        """
        if pattern is not None: 
            # try to get the pattern from the predefined patterns via key
            _pattern = aux.from_kwargs(pattern, None, assay_patterns)
            _pattern = pattern if _pattern is None else _pattern
            self._pattern = re.compile(_pattern, *flags)
        return self._pattern

    def max_assay_name_length(self, length = 20):
        """
        Sets the maximum allowed name length (number of characters) assay names.
        
        Parameters
        ----------
        length : int
            The maximum number of characters to store for the assay name. 
            Default is `length = 20` characters.
        """
        self._max_assay_name_length = length
    
    def parse(self, **kwargs):
        """
        A wrapper for find_assays+find_columns+make_dataframes
        This is the functional core of the Parser's `pipe` method.

        Parameters
        -------
        **kwargs
            Any additional keyword argument that will be passed to any of the wrapped methods.
        """
        decorator = aux.from_kwargs("decorator", None, kwargs, rm = True)
        if decorator is not None:
            self.find_by_decorator(decorator = decorator, **kwargs)
        else: 
            self.find_assays(**kwargs)
        
        # ignore if no assays were found (default is false, unless we use multi-assay multi-sheet files)
        ignore_empty = aux.from_kwargs("ignore_empty", False, kwargs)
        if ignore_empty:
            try: 
                self.find_columns()
                self.make_dataframes(**kwargs)
            except: 
                pass
        else: 
            self.find_columns()
            self.make_dataframes(**kwargs)

    def find_by_decorator(self, decorator : str, col = 0, **kwargs):
        """
        Parses through a column of the datafile and finds all assays that are decorated with a specific decorator.
        Note that this requires that the decorator is in the cell above the assay header. Also, make sure to specify
        an `assay_pattern` to extract the assay name. If no `assay_pattern` is provided, it will simply take the entire cell content.
        
        Parameters
        -----------
        decorator : str
            One of the available `qpcr-decorator`'s for irregular multi-assay files. 
            Available decorators can be assessed via the `qpcr.Parsers.decorators` dictionary keys.
        col : int
            The column in which to look for assay identifiers. 
            By default the first column `col = 0`.
        
        Returns
        -------
        assay_indices : np.ndarray
            The indices (row, col) of all assays found.
        names : np.ndarray
            The extracted names of all assays found.
        """
        # ignore if no assays were found (default is false, unless we use multi-assay multi-sheet files)
        ignore_empty = aux.from_kwargs("ignore_empty", False, kwargs)
        # get the pattern required (or raise error if invalid decorators are provided)
        if decorator not in decorators.keys():
            aw.HardWarning("Parser:invalid_decorator", d = decorator, all_d = list(decorators.keys()))

        decorator_pattern = re.compile( decorators[decorator] )
        decorator_indices, decorator_names = self.find_assays(col = col, pattern = decorator_pattern, **kwargs )

        # check if decorators were identified
        found_indices = decorator_indices.size > 0
        if not found_indices:
            if not ignore_empty: # if none were found either raise error or ignore
                aw.HardWarning("Parser:no_decorators_found")
            else: 
                return

        # if no assay_pattern was specified then default to generic "all" to get full cell contents
        if self.assay_pattern() is None:
            self._WARNING_decorators_but_no_patterns.trigger()
            self.assay_pattern("all")

        assay_indices = decorator_indices
        # get assay indices as the cells IMMEDIATELY BELOW the decorators
        # we adjust either col or the row indices depending on the transposition
        if self._transpose:
            col = col + 1
        else:
            assay_indices = assay_indices + 1
      
        # get all assay header cells into an array to extract their names
        array = self._prep_header_array(col, assay_indices)

        # adjust avaliable length of stored assay names
        max_length = max(
                            list(  map(len, array)  )
                        )
        self.max_assay_name_length(max_length)

        names = np.array(["-"*self._max_assay_name_length for _ in range(len(array))]) # we need to pre-specify the max allowed length for the assay names by filling an array with some dummy placeholders ('-')
        idx = 0
        for entry in array:
            # try:
            match = self.assay_pattern().search(entry)
            if match is not None: 
                name = match.group(1)
                names[idx] = name
            # except: 
            #     continue
            idx += 1

        self._assay_indices = assay_indices
        self._assay_names = names
        
        return assay_indices, names

    def find_assays(self, col = 0, **kwargs):
        """
        Parses through a column of the datafile and identifies all indices of cells that match the provided `assay_pattern``.
        It stores these values internally and also returns the results as numpy arrays.

        Parameters
        -----------
        col : int
            The column in which to look for assay identifiers. 
            By default the first column `col = 0`.

        Returns
        -------
        indices : np.ndarray
            The indices (row, col) of all assays found.
        names : np.ndarray
            The extracted names of all assays found.
        """
        
        custom_pattern = aux.from_kwargs("pattern", None, kwargs)
        if self._pattern is None and custom_pattern is None: 
            aw.HardWarning("Parser:no_pattern_yet")

        pattern_to_use = self._pattern if custom_pattern is None else custom_pattern

        array = self._prep_header_array(col = col)

        indices = np.zeros(len(array))
        names = np.array(["-"*self._max_assay_name_length for _ in range(len(array))]) # we need to pre-specify the max allowed length for the assay names by filling an array with some dummy placeholders ('-')
        idx = 0
        for entry in array:
            # try: 
            match = pattern_to_use.search(entry)
            if match is not None: 
                name = match.group(1)
                names[idx] = name
                indices[idx] = 1
            # except: 
            #     continue
            idx += 1
        indices = np.argwhere(indices == 1)

        if indices.size == 0:
            # ignore if no assays were found (default is false, 
            # unless we use multi-assay multi-sheet files)
            ignore_empty = aux.from_kwargs("ignore_empty", False, kwargs)
            if not ignore_empty:
                aw.HardWarning("Parser:no_assays_found", traceback = False)

        names = names[indices]
        names = names.reshape(len(names))

        self._assay_indices = indices
        self._assay_names = names
        
        return indices, names
    
    def find_columns(self):
        """
        Identifies the relevant data column belonging to each assay within the datafile.
        """
        # search indices of the starts of id and ct columns
        # these are now the row, col coordinates of each name_column header
        name_col_starts = self._find_column_starts(
                                                    label = self._id_label, 
                                                    ref_indices = self._assay_indices
                                                )
        # these are now the row, col coordinates of each ct_column header
        ct_col_starts = self._find_column_starts(
                                                    label = self._ct_label, 
                                                    ref_indices = self._assay_indices
                                                )

        # now we need to generate know also the end indices of the datacolumns
        name_col_ends = self._find_column_ends(name_col_starts)
        
        # now that we know the end indices for the replicate id column we will adopt the end row indices
        # onto the ct column as well (we don't parse through the Ct column because it might have missing 
        # Ct values intersperced which would prematurely terminate the parsing...)

        # (1) we transpose to have all row indices easily accessible in the first line
        # (2) we adopt row indices from the transposed name col
        # (3) and transpose back to get our final ct end indices
        ct_col_ends = deepcopy( np.transpose(ct_col_starts) )
        name_col_ends_t = np.transpose(name_col_ends)
        ct_col_ends[0] = name_col_ends_t[0]
        ct_col_ends = np.transpose(ct_col_ends)

        # now store the data
        self._assay_names_start_indices = name_col_starts
        self._assay_names_end_indices = name_col_ends
        self._assay_ct_start_indices = ct_col_starts
        self._assay_ct_end_indices = ct_col_ends

    def make_dataframes(self, allow_nan_ct : bool = True, default_to : float = None, **kwargs):
        """
        Generates a set of `pandas DataFrame`s each containing two columns 
        (one for the replicate identifiers, one for the Ct values)
        for subsequent use with the main `qpcr` module.

        Parameters
        ------
        allow_nan_Ct : bool
            Allows Ct values to be NaN within the final dataframe (if `True`, default).
            If no NaN Ct values should be maintained a default value for NaN Ct values must be specified
            using `default_to`.
        default_to : float
            The default value to replace NaN Ct values with. 
            This is ignored if `allow_nan_ct = True`.
        """  

        adx = 0
        # print(self._assay_names, self._assay_indices, self._assay_names_start_indices, self._assay_names_end_indices)
        for assay in self._assay_names:
            
            # get the assay's indices of both replicate id and ct columns
            names_start = self._assay_names_start_indices[adx]
            names_end = self._assay_names_end_indices[adx]
            ct_start = self._assay_ct_start_indices[adx]
            ct_end = self._assay_ct_end_indices[adx]

            # generate the final index slices from the total array
            # of both replicate id (names) and ct columns
            names_range, names_col = self._make_index_range(names_start, names_end, crop_first = True)
            ct_range, ct_col = self._make_index_range(ct_start, ct_end, crop_first = True)

            # get the assay data
            assay_names = self._data[names_range, names_col]
            assay_cts = self._data[ct_range, ct_col]

            # and convert to numeric data
            try: 
                assay_cts = assay_cts.astype(float)
            except ValueError as e:
                assay_cts = np.genfromtxt(  np.array(assay_cts, dtype=str)  )
                bad_value = e.__str__().split(": ")[1]
                aw.SoftWarning("Parser:found_non_readable_cts", assay = assay, bad_value = bad_value)

            # assemble the assay dataframe 
            assay_df = pd.DataFrame(
                                    {
                                        standard_id_header : assay_names, 
                                        standard_ct_header : assay_cts,
                                    }
                                )

            if not allow_nan_ct:
                if not isinstance(default_to, (int, float)): 
                    aw.HardWarning("Parser:no_ct_nan_default", d = default_to)
                
                # apply defaulting lambda function
                assay_df[ standard_ct_header ] = assay_df[ standard_ct_header ].apply(
                                                                                        lambda x: x if x == x else default_to
                                                                                    )
            # and store dataframe
            self._dfs.update(
                                { assay : assay_df }
                            )

            adx += 1

    def _make_BigTable_range(self, **kwargs):
        """
        Generates a pandas DataFrame of a subsection of an irregular datafile
        containing a "big data table" with multiple assays specified in it. 

        It makes use of the `id_label` specified using `_CORE_Parser.labels` as the
        anchor. The resulting dataframe fill contain all rows from the cell where `id_label``
        as located until the data is empty. 

        If additionally `replicates` are specified in the `kwargs` 
        the starting positions of assay replicates are inferred based on `decorators`. 
        Note, this only works for `horizontal` Big Tables!
        """
        is_horizontal = aux.from_kwargs("is_horizontal", False, kwargs)
        
        # get the main data
        data = self._data.astype("str")
        ref_col_header = self._id_label   
        
        # find big table starting row
        idx = np.argwhere(data == ref_col_header)

        # vet that we actually found the big table
        if idx.size == 0:
            aw.HardWarning("Parser:no_bigtable_header", header = ref_col_header)
        idx = idx.reshape(idx.size)
        start, col = idx

        idx = 1
        while True:
            try: 
                entry = data[start+idx, col]
            except: 
                break
            if entry == "nan":
                break
            idx += 1

        end = start + idx

        # put a -1 offset on the rows if horziontal, as the decorators 
        # are in the row above the actual column headers.
        if is_horizontal:
            start -= 1

        # generate bigtable data range and store
        relevant_data = data[start : end, : ]
        self._bigtable_range = relevant_data  


    def _infer_BigTable_groups(self, **kwargs):
        """
        Gets the group ranges from the bigtable datarange.
        Note, this is only used in case of horizontal big tables.
        """
        # get the relevant data
        array = self._bigtable_range
        maxrows, allcols = array.shape

        ignore_empty = aux.from_kwargs("ignore_empty",False,kwargs)

        # get and vet replicates
        replicates = aux.from_kwargs("replicates", None, kwargs, rm = True)
        if replicates is None: 
            aw.HardWarning("Parser:bigtable_no_replicates", traceback = False)
        replicates, names = self._vet_replicates(ignore_empty, replicates, array, **kwargs)

        
        rdx = 0 # counter for the replicate groups
        
        data_array = None # this array will store the entire transposed data

        decorator = plain_decorators["qpcr:group"]

        # now get the assays in question
        # we already vetted if the file is properly 
        # decorated during _vet_replicates
            
        # find decorated starting columns
        # get only column indices            
        indices = np.argwhere(array == decorator)
        indices = indices[ : , 1 ]

        if indices.size == 1:
            indices = [indices]

        # get total slice of bigtable rows
        rows = slice( 1, maxrows )
        
        # iterate over each group
        for col in indices: 
            
            rep = replicates[rdx]
            name = names[rdx]

            # get data columns
            cols = slice( col, col + rep )
            
            # get data
            data = array[  rows, cols  ]
            
            # rename data cols uniformly
            data[ 0 , : ] = name
            
            # concatenate data into a single array
            if data_array is None: 
                data_array = data
            else:
                data_array = np.concatenate(
                                                ( data_array, data ),
                                                axis = 1,
                                            )
            rdx += 1

        # remove groups from the data array
        groups = data_array[ 0, : ]
        data_array = data_array[ 1:, : ]

        # reshape data into a single column
        data_array = np.concatenate(data_array, axis = 0)

        # now repeat the groups to match the stacked new data column
        groups_tiled = np.tile(groups, data_array.size // groups.size )

        # now get the dataset id column
        id_col = self._id_label
        id_start = np.where(array == id_col)
        id_row, id_col = id_start
        id_rows = slice( int(id_row + 1), maxrows )
        id_col = array[  id_rows, id_col  ]
        
        # repeat dataset ids to match stacked new data column
        ids_tiled = np.repeat(  id_col , groups.size  )

        # assemble all data
        total_data = [ids_tiled, groups_tiled, data_array]
        headers = [ default_dataset_header, standard_id_header, standard_ct_header ]  

        # check for qpcr column and if present, get and adjust shape
        self._BigTable_horizontal_qpcr_col(array, maxrows, groups, total_data, headers)


        # combine the three columns (dataset id, groups, and Ct (actual data_array))
        data_array = np.stack( total_data, axis = 1 )

        # add default names into the first row
        data_array = np.concatenate(
                                        ( 
                                            [headers],
                                            data_array
                                       ),   axis = 0
                                )
        
        # actually return the finished array
        return data_array

    def _BigTable_horizontal_qpcr_col(self, array, maxrows, groups, total_data, headers):
        """
        Checks if a "@qpcr" column is present in the data and if so, adjusts its shape and 
        adds it to the data to be assembled for the assays.
        """
        column_decorator = plain_decorators["qpcr:column"]
        qpcr_col = np.where(array == column_decorator)
        qpcr_row, qpcr_col = qpcr_col
        if len(qpcr_col) != 0:
            
            qpcr_rows = slice( int(qpcr_row + 1), maxrows )
            qpcr_col = array[  qpcr_rows, qpcr_col  ]

            qpcr_tiled = np.repeat(  qpcr_col , groups.size  )
            total_data.append(qpcr_tiled)
            headers.append("@qpcr")

    def _vet_replicates(self, ignore_empty, replicates, array, **kwargs):
        """
        Checks if provided replicates cover all data groups (annoated columns).
        And it also gets the names supposed to be used for the columns.
        """
        # get assays for each decorator
        groups = 0
        decorator = plain_decorators["qpcr:group"]
    
        # find decorated starting columns
        indices = np.argwhere(array == decorator)
        if indices.size == 0 and not ignore_empty:
            aw.HardWarning("Parser:no_decorators_found", traceback = False)

        # get only column indices            
        indices = indices[ : , 1 ]
        groups += indices.size

        # check if replicates are an integer, if so transform to 
        # tuple that cover all found assays
        if aux.same_type(replicates, 1): 
            replicates = np.tile([replicates], groups)
        
        # if replicates are already a tuple, make sure they cover all rows
        elif aux.same_type(replicates, ()):
            all_covered = groups == len(replicates)
            if not all_covered:
                aw.HardWarning("Assay:reps_dont_cover", n_samples = groups, reps = replicates, traceback = False)
            
        # get names for assays
        group_names = aux.from_kwargs("names", [ default_group_name.format(i) for i in range( len(replicates) )  ], kwargs)

        # vet that names cover
        if len(group_names) != len(replicates):
            aw.HardWarning("Assay:groupnames_dont_colver", traceback = False, current_groups = f"None, but needs to be {len(replicates)} names.", new_received = group_names)
        
        # return tranformed replicates
        return replicates, group_names

        

    def _prep_header_array(self, col = None, row = None):
        """
        Generates the array in which header entries should be searched for
        """
        
        if row is None and col is not None: 
            array = self._data[:, col] if not self._transpose else self._data[col, :]
        elif row is not None and col is None:
            array = self._data[row, :] if not self._transpose else self._data[:, row]
        elif row is not None and col is not None:
            array = self._data[row, col] if not self._transpose else self._data[col, row]
        else:
            aw.HardWarning("Parser:invalid_range")

        # re-format to str and reset "nan" to dummy_blank
        array = array.astype(str)
        array[ np.argwhere(array == "nan") ] = dummy_blank
        array = array.reshape(array.size)
        return array


    def _make_index_range(self, start_indices, end_indices, crop_first = True):
        """
        Generates an index range for a data column based on start and stop indices.
        This assumes that the column entry (i.e. entry[1]) is always the same and only the rows are different.
        
        Parameters
        ----------
        start_indices : np.ndarray
            Row, col indices of the header of the data column
        end_indices : np.ndarray
            Row, col indices of the last entry of the data column
        crop_first : bool
            If set to True it will offset the start row indices by +1 to exclude the header.
        """

        start = start_indices[0] + 1 if crop_first else start_indices[0]
        end = end_indices[0]

        row_range = slice(  start, end  )
        col = start_indices[1]
        return row_range, col

    def _find_column_starts(self, label, ref_indices):
        """
        This function uses the assays' found reference row indices to 
        now search for the coordinates of the labeled cell so we know where a data column starts
        """

        # get index to match to ref_indices
        idx_to_match = 0 if not self._transpose else 1

        data = self._data
        all_found = np.argwhere(data == label)
        row_indices = np.transpose(all_found)[ idx_to_match ]

    
        # adjust coordinates +1 as the headers would be in the row below the assay declaration
        # we only have to do this if we use the default setting of assays in the same column
        if not self._transpose:
            ref_indices = ref_indices + 1

        # ref_indices = ref_indices + 1 
        matching_rows = np.where(np.isin(row_indices, ref_indices))
                
        # if no matches were found, try incrementing the index offset once more 
        # (we'll allow for a single row between the header and the start of the data)
        no_matches = len(matching_rows) == 1 and matching_rows[0].size == 0
        if no_matches:
            ref_indices = ref_indices + 1
            matching_rows = np.where(np.isin(row_indices, ref_indices))

        # check again, and raise Error if still no matches are found
        no_matches = len(matching_rows) == 1 and matching_rows[0].size == 0
        if no_matches:
            aw.HardWarning("Parser:no_data_found", label = label)

        matching_rows = all_found[matching_rows]
        return matching_rows

    def _find_column_ends(self, indices):
        """
        Determines the end index of a column within the datafile based on the starting indices of
        its header label
        """
        data = self._data
        end_indices = np.zeros(indices.shape, dtype = int)
        adx = 0
        for i in indices:
            row, col = i
            idx = 0
            value = 0
            while True:
                try: value = data[row + idx, col]
                except: break
                if value != value: 
                    break
                idx += 1
            finals = np.array([row + idx, col], dtype = int)
            end_indices[adx] += finals
            adx += 1
        return end_indices


class CsvParser(_CORE_Parser):
    """
    Handles reading and parsing irregular `csv` files that contain multiple assays.
    It extracts datasets either through regex pattern matching or/and through provided
    decorators within the datafile.
    """
    def __init__(self):
        super().__init__()

    def pipe(self, filename :str, **kwargs):
        """
        A wrapper for read+parse

        Note 
        ----
        This is the suggested use of `CsvParser`. 
        If a directory has been specified into which the datafiles shall be saved, 
        then saving will automatically be done.

        Parameters
        -------
        filename : str
            A filepath to an input csv file.
        **kwargs
            Any additional keyword argument that will be passed to any of the wrapped methods.
        Returns
        -------
        assays : dict
            A dictionary of all the extracted assays from the datafile storing the data as pandas DataFrames.
            Individual assays can also be accessed using the `get` method.
        """
        try: 
            self.read(filename, **kwargs)
        except: 
            self.read(filename)
            aw.SoftWarning("Parser:incompatible_read_kwargs", func = "pandas.read_csv")
        
        self.parse(**kwargs)
        assays = self.get()
        
        if self._save_loc is not None: 
            self.save()
        
        return assays

    def read(self, filename : str, **kwargs):
        """
        Reads an input csv file. 

        Parameters
        -------
        filename : str
            A filepath to an input csv file.
        **kwargs 
            Any additional keyword arguments to be passed to pandas' `read_csv` function.
        """
        self._src = filename

        contents = self._prepare_commas()
        contents = StringIO(contents) # convert to StringIO for pandas to be able to read
        
        delimiter = ";" if self._is_csv2() else ","
        delimiter = aux.from_kwargs("sep", delimiter, kwargs, rm = True)

        # now read the data and convert to numpy array
        try: 
            df = pd.read_csv(contents, header = None, sep = delimiter, **kwargs)
        except: 
            aw.SoftWarning("Parser:incompatible_read_kwargs", func = "pandas.read_csv()")
            df = pd.read_csv(contents, header = None, sep = delimiter)

        df = df.dropna(axis = 0, how = "all").reset_index(drop=True)
        data = df.to_numpy()

        self._data = data

    def _is_csv2(self):
        """
        Tests if csv file is ; delimited (True) or common , (False)
        """
        with open(self._src, "r") as openfile: 
            content = openfile.read()
        if ";" in content: 
            return True
        return False

    def _prepare_commas(self):
        """
        This function reads the datafile and adjusts the number of commas 
        within each line to ensure equal commas in the entire file.

        Note
        -------
        Although the method uses the term "commas" it also works with semicolons for csv2

        Returns
        -------
        new_content : str
            A string containing the entire file contents with adjusted commas.
        """

        delimiter = ";" if self._is_csv2() else ","

        # check if quotes are in datafile and adjust comma-patterns to use
        empty_comma_filler = f'{delimiter}""' if self._has_quotes() else f"{delimiter}"
        comma_sep = f'"{delimiter}"' if self._has_quotes() else f"{delimiter}"
        comma_sep = re.compile(comma_sep)

        with open(self._src, "r") as f:
            content = f.read()
            lines = content.split("\n")
            comma_counts = [len(comma_sep.findall(i)) for i in lines]
            max_commas = max(comma_counts)
            lines = [i + (max_commas - j) * empty_comma_filler for i, j in zip(lines, comma_counts)]
        new_content = "\n".join(lines)
        return new_content

    def _has_quotes(self):
        """
        Checks if cells from the csv input file have quotes around them.
        Essentially it checks if there are any "," patterns in the file.
        """
        delimiter = ";" if self._is_csv2() else ","
        with open(self._src, "r") as f:
            content = f.read()
        has_quotes = f'"{delimiter}"' in content
        return has_quotes

class ExcelParser(_CORE_Parser):
    """
    Handles reading and parsing `excel` files that may contain multiple assays.
    It extracts datasets either through regex pattern matching or/and through provided
    decorators within the datafile.
    """
    def __init__(self):
        super().__init__()

    def read(self, filename : str, sheet_name : (str or int) = 0, **kwargs):
        """
        Reads an input excel file. 

        Parameters
        -------
        filename : str
            A filepath to an input excel file.
        sheet_name : int or str
            The name of a specific spreadsheet of the file to read.
            If none is provided by default the first sheet will be read.
            Only one single sheet can be read at a time. 
            If an `integer` is provided the sheets will be accessed by their order, otherwise by their name (if a `string` is provided).
        **kwargs
            Any additional keyword arguments to be passed to pandas `read_excel` function.
        """
        self._src = filename

        # read data and convert to numpy array
        try: 
            data = pd.read_excel(self._src, sheet_name = sheet_name, **kwargs)
        except: 
            data = pd.read_excel(self._src, sheet_name = sheet_name)
            aw.SoftWarning("Parser:incompatible_read_kwargs", func = "pandas.read_excel()")

        data = data.to_numpy()

        self._data = data

    def pipe(self, filename :str, **kwargs):
        """
        A wrapper for read+parse

        Note 
        ----
        This is the suggested use of `ExcelParser`. 
        If a directory has been specified into which the datafiles shall be saved, 
        then saving will automatically be done.

        Parameters
        -------
        filename : str
            A filepath to an input excel file.
        **kwargs
            Any additional keyword argument that will be passed to any of the wrapped methods.
        Returns
        -------
        assays : dict
            A dictionary of all the extracted assays from the datafile storing the data as pandas DataFrames.
            Individual assays can also be accessed using the `get` method.
        """
        try: 
            self.read(filename, **kwargs)
        except: 
            self.read(filename)
            aw.SoftWarning("Parser:incompatible_read_kwargs", func = "pandas.read_excel")

        self.parse(**kwargs)
        assays = self.get()

        if self._save_loc is not None: 
            self.save()

        return assays

if __name__ == "__main__":
    
    # parser = CsvParser()
    # parser.assay_pattern("Rotor-Gene")
    # parser.save_to("__csvparser")
    # mycsv = "./__parser_data/Brilliant III Ultra Fast SYBR Green 2019-01-07 (1).csv"
    # parser.pipe(mycsv)

    # print("""\n\n\n ========================= \n All good with CsvParser \n ========================= \n\n\n""")

    # parser2 = ExcelParser()
    # parser2.assay_pattern("Rotor-Gene")
    # parser2.save_to("./__excelparser")
    # myexcel = "./__parser_data/excel 3.9.19.xlsx"
    # parser2.pipe(myexcel, sheet_name = 1)

    # print("""\n\n\n ========================= \n All good with ExcelParser \n ========================= \n\n\n""")

    # parser3 = ExcelParser()
    # decorated_excel = "./__parser_data/excel 3.9.19_decorated.xlsx"
    # parser3.save_to("./__decorated_excelparser")
    # parser3.read(decorated_excel)
    # parser3.assay_pattern("Rotor-Gene")
    # parser3.find_by_decorator(decorator = "qpcr:all")
    # parser3.find_columns()
    # parser3.make_dataframes()
    # parser3.save()
    # # print(parser3.get())

    # print("""\n\n\n ========================= \n All good with decorated ExcelParser \n ========================= \n\n\n""")

    # parser4 = CsvParser()
    # decorated_csv = "./__parser_data/Brilliant III Ultra Fast SYBR Green 2019-01-07 (1)_decorated.csv"
    # parser4.save_to("./__decorated_csvparser")
    # parser4.read(decorated_csv)
    # parser4.assay_pattern("Rotor-Gene")
    # parser4.find_by_decorator(decorator = "qpcr:all")
    # parser4.find_columns()
    # parser4.make_dataframes()
    # parser4.save()
    # # print(parser4.get())

    # print("""\n\n\n ========================= \n All good with decorated CsvParser \n ========================= \n\n\n""")


    # parser4 = CsvParser()
    # decorated_csv = "./__parser_data/Brilliant III Ultra Fast SYBR Green 2019-01-07 (1)_decorated.csv"
    # parser4.save_to("./__decorated_csvparser_pipe")
    # parser4.assay_pattern("Rotor-Gene")
    # # parser4.pipe(decorated_csv, decorator = "qpcr:assay")
    # # print(parser4.get())

    # parser4.pipe("./__parser_data/manual_decorated.csv")

    # print("""\n\n\n ========================= \n All good with decorated CsvParser using pipe \n ========================= \n\n\n""")

    # parser3 = ExcelParser()
    # decorated_excel = "./__parser_data/excel 3.9.19_decorated.xlsx"
    # parser3.save_to("./__decorated_excelparser_pipe_nodec")
    # parser3.assay_pattern("Rotor-Gene")
    # parser3.pipe(decorated_excel)
    # # print(parser3.get())

    # print("""\n\n\n ========================= \n All good with decorated ExcelParser using pipe without dec\n ========================= \n\n\n""")


    # same_row_assays = "/Users/NoahHK/Downloads/qPCR cytokines upon treatment_decorated.xlsx"

    # parser5 = ExcelParser()
    # parser5.transpose()
    # parser5.read(same_row_assays, sheet_name = 1)
    # parser5.save_to("./__transposed_parser")
    # parser5.assay_pattern("all")
    # parser5.labels(  id_label = "Sample Name", ct_label = "CT"  )

    # print("""\n\n\n ========================= \n Transposed excel (FIND)\n ========================= \n\n\n""")

    # parser5.find_assays(col = 1)
    # parser5.find_columns()
    # parser5.make_dataframes()
    # r = parser5.get()
    # print(r)

    # # assert parser5._assay_indices is not None, "(find_assays) No assay_indices could be found!!!"
    # parser5.clear()

    # print("""\n\n\n ========================= \n Transposed excel (DECO)\n ========================= \n\n\n""")


    # parser5.find_by_decorator("qpcr:all")
    # # assert parser5._assay_indices is not None, "(find_by_decorator) No assay_indices could be found!!!"
    # parser5.find_columns()
    # parser5.make_dataframes()
    # r = parser5.get()
    
    # parser5.clear()

    # parser5.read(same_row_assays, sheet_name = 1)
    # parser5.parse(decorator = "qpcr:all")
    # r = parser5.get()
    # print(r)
    # parser5.save()

    bigtable_horiztonal = "/Users/NoahHK/Downloads/Local_cohort_Adenoma_qPCR_rawdata_decorated.xlsx"

    parser_bigtable = ExcelParser()
    parser_bigtable.read(bigtable_horiztonal)

    parser_bigtable.labels( id_label = "tissue_number" )
    parser_bigtable._make_BigTable_range(is_horizontal = True)
    r = parser_bigtable._infer_BigTable_groups( replicates = (3, 2) )
    print(r)