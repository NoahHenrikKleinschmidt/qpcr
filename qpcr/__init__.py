"""
This module is designed to provide functions to analyse qPCR data. 
It is designed for maximal user-friendliness and streamlined data-visualisation.

Terminology
-------

Let's talk some terminology first. There are a number of important terms that you will meet when using the `qpcr` module that you may feel have nothing to do with qPCR. 
Read on to learn what these terms mean and why they are important. 

#### File (or "datafile")
A `file` is simply a datafile. For the `qpcr` module this means either a `csv` or an `excel` file. Files are the at the very basis of our data pipeline. 
However, we do not strictly assume that these files correspond to qPCR assays per se. The default settings of the `qpcr` processing classes are designed to work with 
datafiles that contain Ct values of _one single qPCR assay_. However, if your experimental setup looks differently, there are ways to adapt `qpcr` to handle different setups.
However, if your data follows the above mentioned standard arrangements, you can think of a file as identical to a "qPCR assay".

#### An `qpcr.Assay`
The `qpcr.Assay` class is used to store the Ct values from a datafile. The class was named with the assumption in mind that each qPCR assay is stored as a separate datafile, but again,
if this does not match your setup then this class will still handle the data of your datafiles.

#### Replicates
Each row within your datafiles corresponds to one replicate. Hence, a _replicate_ is just a single pair of some identifier and a corresponding Ct value. 
That's it for what "replicates" _are_. Now to how we deal with them. The `qpcr.Assay` class has an attribute called `replicates` which must be set 
if you want to read your data. `replicates` can be either an `integer` or a `tuple`. Why's that? Well, normally we perform experiments as "triplicates", or "duplicates", or whatever multiplets. 
Hence, we have different individual replicates that actually belong together (they form a "group", see the next paragraph!). To tell the `qpcr` module that your experiment was done as 
triplicates for instance, you would specify `replicates = 3`. Of course, most qPCR machines already output averaged Ct values alongside the individual ones. If you wish to use these, 
then you might specify `replicates = 1`, for instance (there is no strict contraint to have multipltets). On the other hand, if you did biological triplicates and technical duplicates you can choose
if you would like keep the triplicates separately or merge them together. If you choose the latter then you would assign `replicates = 6` (2x3 = 6), otherwise you might keep to `replicates = 2`. 
And so forth, you get the idea... Why then is a tuple also a valid input for `replicates`? Well, if you manually removed outliers from your data, then your groups of replicates would no longer all have the same length. 
You can account for this by specifying the number of replicates in each group in a tuple. 
However, at this point it must be mentioned that `qpcr` comes with `Filters` that will handle this for you automatically, so there is no need to tamper with your datafiles manually!

#### Groups (of Replicates)
This is one of the most fundamental terms. In the previous paragraph we already started to talk about "groups of replicates". 
Groups of replicates are, as their name already implies, well, a number of replicates that somehow belong together.
For most cases, if your datafiles contain one assay each, then the groups of replicates are most likely your different qPCR samples / experimental conditions. 
We talk about _groups_ of replicates instead of samples or conditions, primarily, because there might be different data setups so that these terms might not be always appropriate.
If your data follows default arrangements, however, then a _group of replicates_ is just what you would think of as a "qPCR sample". 
Groups are assigned a numeric index starting from 0, which is how they are identified by the classes. 
However, they also come with a text label called the `group_name` (you can manually set and re-set the group names as you like). Some classes such as the `qpcr.SampleReader` class will actually just use the term `names` instead of the full `group_names`. 
Whenever you see anything "names"-related it is (super-duper most likely) a reference to the `group_names`.

#### "Sample"
We described in the previous section how we use the term "group of replicates". However, you may find that there is also a term "Sample" within `qpcr`'s vocabulary. 
As far as the `qpcr` module is concerned, the term "Sample" is intercangable with the term "Replicate" (see above) but has nothing to do with qPCR samples or experimental conditions as a wet-lab biochemist would understand them. 
You will find that your dataframes contain a "Sample" column that contains the original replicate identifiers from your datafiles. 
Why is this? Essentially because the `qpcr` module does not make strict assumptions of what data is actually contained within the files, it simply considers all replicate values as its "sample data" to work on.
So would another name have also worked here? Sure, it would have.

#### `Delta-Ct` vs `Delta-Delta-Ct` vs `normalisation` ???
Now it gets even more technical (sorry).
The default analysis workflow in $\Delta \Delta Ct$ analysis is to first calculate a $\Delta Ct$ using an intra-assay reference and then calculate the $\Delta \Delta Ct$ using a normaliser assay. 
The first $\Delta Ct$ step is performed by a class called `qpcr.Analyser` using its native method `DeltaCt()`. 
Why just `DeltaCt()` and not `DeltaDeltaCt()`? Well, we call this second `Delta`-step in `DeltaDeltaCt` differently. 
We name it `normalisation`, and it is handled by a class called `qpcr.Normaliser` using its native method `normalise()`. 
So, as far as the `qpcr` module is concerned there is only `DeltaCt()` which performs the first $\Delta Ct$, and `normalise()` which later handles the second "delta"-step to get to $\Delta \Delta Ct$. 
Of course, this means that the `qpcr.Normaliser` will need to have knowledge about which `qpcr.Assay` objects contain actual assays-of-interest and which ones contain normaliser-assays (specifying that is easy, though, so don't worry about that).

#### The `anchor` and the "reference group"
Next to the groups of replicates, this is probably one of the most important terms to come to grips with. The `anchor` is simply the intra-dataset reference used by the `qpcr.Analyser` to perform its first $\Delta Ct$. 
If your datafiles contain one assay each, and your groups of replicates are your qPCR samples, then you will likely have some "wildtype", "untreated", or "control" sample, right? Well, in `qpcr` terms that would be your _reference group_.
Usually your `anchor` is part of or generated from the Ct values of your _reference group_ (like their mean for instance).
By default it is assumed that your reference group is the _very first_ group of replicates. However, it's not a big problem if this is not the case, as you can specify different anchors easily.
So, again, the `anchor` is the dataset-internal reference value used for the first $\Delta Ct$.

#### "assays" vs "normalisers"
You will likely encounter methods and/or arguments that speak of "assays" and "normalisers", especially with the `qpcr.Normaliser`. 
For all intents and purposes, these terms simply refer to datasets stored as `qpcr.Assay` or `qpcr.Results` objects (ultimately they refer to your datafiles). 
"Assays" are the short notation for assays-of-interest. 
"Normalisers" on the other hand are your normaliser-assays (from housekeeping genes like ActinB for instance). 
You will also find that the term "assays" is used within the final results dataframe (when using the summary-statistics mode). 
In this setting "assays" refers to the assay-of-interst whose data was analysed according to the provided normaliser-assays. 
In fact, the entries within the "assay" column of the final results are a combined ID of all the assays (-of-interest AND -normaliser) that were involved in computing the results (check out what the final tables look like and it'll be immediately clear).
Please, do not confuse these with the `qpcr.Assay` or `qpcr.Normaliser` classes (these are always referenced in singular never in plural and are technical terms while the others are conceptual terms)!

The most crucial terms have now been discussed. There are some more important terms that are discussed below, but they should not pose problems if you refer to the API documentation.

#### "Results" = my final results?
Yes and no. Anything that is computed through any of the `qpcr` classes is called a "result" of some kind.
In practice as soon as you pass your data through an `qpcr.Analyser` or `qpcr.Normaliser` you generate some "results" which are actually stored in a separate class called `qpcr.Results` (that's not so important, though).
The `qpcr.Results` that comes out of your `qpcr.Normaliser` will probably be what you would think of as your "final results", yes. The one from the `qpcr.Analyser` will probably not be.

#### `get`ting your data
Too many classes and objects? Well, no worries, the underlying data is stored as `pandas DataFrames`. To get your data from the clutches of the `qpcr` processing classes you can always use the `get()` method. 
`get` is pretty universal in the `qpcr` module, so whenever you want to extract your data, there's a `get()` method to help you. 

#### `link` vs `add` vs `pipe`
Different classes have slightly different methods of adding data to them. Classes that only accept one single data input (such as a single `qpcr.Assay` object or a single filepath as `string`)
usually have a `link()` method that, well, links the data to them. After that the classes are ready to perform whatever actions they can perform (an `qpcr.Analyser` would perform `DeltaCt()` for instance). 
Some classes such as the `qpcr.Analyser` have a wrapper that will call both their `link()` as well as their actual core-functional method together in one go. This wrapper is called `pipe()`. 
So for the `qpcr.Analyser` you could either manually use `link()` and then `DeltaCt()`, or simply call `pipe()` which does both for you. It is noteworthy that `pipe` methods actually _return_ whatever
their output is, which is *not* normally the case otherwise (normally you'd use the `get()` method to extract your data, see above).
Alright, we know about `link()` and `pipe()` now, what about `add`?  Classes that accept multiple inputs have `add` methods, which tells the class where exactly to store the input data. 
`add`-methods are especially implemented within the pre-defined analysis pipelines of the `qpcr.Pipes` submodule. You will probably often use the methods `add_assays()` and `add_normalisers()` if you plan on using these predefined pipelines.
However, these classes usually still have a `link()` method somewhere that you can use as well. For instance, the `qpcr.Normaliser` uses a `link` method to add both assays-of-interest and normaliser-assays simultaneously.


"""

import pandas as pd
import qpcr._auxiliary as aux
from qpcr._auxiliary import warnings as aw
import qpcr.Parsers as Parsers
import os
import numpy as np 
from copy import deepcopy 
from io import StringIO

# default column names for raw Ct data files (don't change this!)
RAW_COL_NAMES = ["Sample", "Ct"]

class _CORE_Reader(aux._ID):
    """
    The class handling the core functions of the Reader class. 
    Both the standard qpcr.Reader as well as the qpcr._Qupid_Reader
    inherit from this. 
    """
    def __init__(self):
        super().__init__()
        self._src = None
        self._delimiter = None
        self._df = None

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


    def read(self, **kwargs):
        """
        Reads the given data file.

        If the data file is an Excel file replicates and their Ct values will be 
        extracted from the first excel sheet of the file. Note, this assumes by default
        that the replicates are headed by the label `"Name"` and the corresponding Ct values
        are headed by the label `"Ct"`. Both labels have to be on the same row. 

        If these labels do not match your excel file, you may
        specify `name_label` and `Ct_label` as additional arguments.
        """
        suffix = self._filesuffix()
        if suffix == "csv":
            try: 
                self._csv_read()
            except:
                parser = Parsers.CsvParser()
                assay_pattern = aux.from_kwargs("assay_pattern", "Rotor-Gene", kwargs)
                assay_of_interest = aux.from_kwargs("assay", None, kwargs, rm=True)
                parser.assay_pattern(assay_pattern)
                parser.pipe(self._src, **kwargs)

                if len(parser.assays()) > 1:
                    if assay_of_interest is None: 
                        aw.HardWarning("Reader:cannot_read_multifile", file = self._src, assays = parser.assays())
                    self._df = parser.get(assay_of_interest)
                    self.id(assay_of_interest)
                else:
                    assay_of_interest = parser.assays()[0]
                    self._df = parser.get(assay_of_interest)
                    self.id(assay_of_interest)

        elif suffix == "xlsx":
            parser = Parsers.ExcelParser()
            assay_pattern = aux.from_kwargs("assay_pattern", "Rotor-Gene", kwargs)
            assay_of_interest = aux.from_kwargs("assay", None, kwargs, rm=True)
            parser.assay_pattern(assay_pattern)
            parser.pipe(self._src, **kwargs)

            if len(parser.assays()) > 1:
                if assay_of_interest is None: 
                    aw.HardWarning("Reader:cannot_read_multifile", file = self._src, assays = parser.assays())
                self._df = parser.get(assay_of_interest)
                self.id(assay_of_interest)
            else:
                assay_of_interest = parser.assays()[0]
                self._df = parser.get(assay_of_interest)
                self.id(assay_of_interest)

            # self._excel_read(**kwargs)

    def _parse_read(self, data, name_label = "Name", Ct_label = "Ct", **kwargs):
        """
        Reads a given data file and locates the relevant columns for sample ids (names) and Ct values.
        This is the funtional core of the _excel_read() method, but also works as a backup read method
        in case a csv converted input file is provided that does not conform to the default two-column structure
        requirements. 

        Parameters
        ----------
        data : np.ndarray
            A numpy array of the raw data file.
        name_label : str
            The header above the replicate sample ids.
        Ct_label : str
            The header above the replicates' Ct values.
        """
        # locate name_label and Ct_label
        names_loc = np.argwhere(data == name_label)
        cts_loc = np.argwhere(data == Ct_label)

        # check if locations match, if not, get the one where both share the same row
        names_loc, cts_loc = self._adjust_loc_shapes(names_loc, cts_loc)

        # get the values from the excel sheet array
        sample_ids = self._get_values_from_range(data, names_loc)
        ct_values = self._get_values_from_range(data, cts_loc)

        # assemble the dataframe
        df = pd.DataFrame(
                            dict(
                                    Sample = sample_ids, 
                                    Ct = ct_values
                                )
                        )
        return df

    def _excel_read(self, name_label = "Name", Ct_label = "Ct", **kwargs):
        """
        Reads the given data file if it's an excel file

        Parameters
        ----------
        name_label : str
            The header above the replicate sample ids.
        Ct_label : str
            The header above the replicates' Ct values.
        """
        # read data and convert to numpy array
        data = pd.read_excel(self._src)
        data = data.to_numpy()

        # parse file and get data
        df = self._parse_read(data, name_label, Ct_label)
        self._df = df
        

    def _get_values_from_range(self, data, indices):
        """
        Gets values from an excel spreadhseet as numpy ndarray starting from the indices
        provided until the cells are filled with np.nan
        """
        row, col = indices
        idx = 0
        value = 0
        while True:
            idx += 1
            try: value = data[row + idx, col]
            except: break
            if value == np.nan: 
                break
        
        values = data[  # we start at row+1 to exclude the header
                        slice(row+1, row + idx+1), 
                        col
                    ]
        return values

    def _adjust_loc_shapes(self, names_loc, cts_loc):
        """
        Adjusts the location indices of the start of sample ids and ct values,
        in case the headers were identified multiple times in the spreadsheet.

        It returns the location indices for both names (ids) and ct values that are 
        on the same row.
        """
        trans_names = np.transpose(names_loc)
        trans_cts = np.transpose(cts_loc)

        # get common index for the row, first apply to names
        common_row_index = np.argwhere(trans_names[0] == trans_cts[0])
        names_loc = names_loc[common_row_index]
        # now also apply to cts
        common_row_index = np.argwhere(trans_cts[0] == names_loc[0][0][0])
        cts_loc = cts_loc[common_row_index]

        # now reduce shape back to 1D
        names_loc, cts_loc = names_loc.reshape(2), cts_loc.reshape(2)

        return names_loc, cts_loc
            
        
    def _csv_read(self, **kwargs):
        """
        Reads the given data file if it's a csv file
        """
        df = None
        try: 
            df = pd.read_csv(
                                self._src, 
                                sep = self._delimiter, 
                                header = self._has_header(), 
                                names = RAW_COL_NAMES
                            )
        except: 
            aw.HardWarning("Reader:cannot_read_csv", file = self._src)

        self._df = df

    def _filesuffix(self):
        """
        Returns the file-suffix of the provided file
        """
        suffix = self._src.split(".")[-1]
        return suffix

class Reader(_CORE_Reader):
    """
    Reads qpcr raw data files in csv or excel format. 

    Parameters
    ----------
    filename : str
        A filepath to a raw data file.
        If the file is a `csv` file, it has to have two named columns; one for replicate names, one for Ct values. 
        Both csv (`,` spearated) and csv2 (`;` separated) are accepted.
        If the file is an `excel` file it the relevant sections of the spreadsheet are identified automatically. 
        But they require identifying headers. By default `Name` and `Ct` are assumed but these can be changed using 
        the `name_label` and `Ct_label` arguments that can be passed as kwargs (they will be forwarded to the `.read()` method). 

    **kwargs
        Any additional keyword arguments that shall be passed to the `read()` method which is immediately called during init.
    """
    def __init__(self, filename:str, **kwargs) -> pd.DataFrame: 
        super().__init__()
        self._src = filename
        if self._filesuffix() == "csv":
            self._delimiter = ";" if self._is_csv2() else ","
        self.read(**kwargs)

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


class _Qupid_Reader(_CORE_Reader):
    """
    This Reader class works with streamlit's UploadedFile class.
    
    Note
    -------
    We have to use a little hack to make the UploadedFile readable.
    We convert to a string and then back to a StringIO which we can then pass to pandas...
    
    Parameters
    ----------
    file
        A streamlit UploadedFile object
    """
    def __init__(self, file) -> pd.DataFrame: 
        super().__init__()
        self._filename = file.name
        if self._filesuffix() == "csv":
            self._content = file.read().decode()
            self._src = StringIO(self._content)
            self._delimiter = ";" if self._is_csv2() else ","
        else: 
            self._src = file
        self.read()

    def _filesuffix(self):
        """
        Returns the file-suffix of the provided file
        """
        suffix = self._filename.split(".")[-1]
        return suffix

    def _is_csv2(self):
        """
        Tests if csv file is ; delimited (True) or common , (False)
        """
        if ";" in self._content: 
            return True
        return False

    def _has_header(self):
        """
        Checks if column headers are provided in the data file
        It does so by checking if the second element in the first row is numeric
        if it is numeric (returns None << False) no headers are presumed. Otherwise
        it returns 0 (as in first row has headers)...
        """
        content = self._content.split("\n")[0]
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

    def replicates(self, replicates : (int or tuple or str) = None):
        """
        Either sets or gets the replicates to be used for grouping the samples
        Before they are assigned, replicates are vetted to ensure they cover all data entries.

        Parameters
        ----------
        replicates : int or tuple or str
            Can be an `integer` (equal group sizes, e.g. `3` for triplicates), 
            or a `tuple` (uneven group sizes, e.g. `(3,2,3)` if the second group is only a duplicate). 
            Another method to achieve the same thing is to specify a "formula" as a string of how to create a replicate tuple.
            The allowed structure of such a formula is `n:m,` where `n` is the number of replicates in a group and `m` is the number of times
            this pattern is repeated (if no `:m` is specified `:1` is assumed). So, as an example, if there are 12 groups which are triplicates, but
            at the end there is one which only has a single replicate (like the commonly measured diluent qPCR sample), we could either specify the tuple
            individually as `replicates = (3,3,3,3,3,3,3,3,3,3,3,3,1)` or we use the formula to specify `replicates = "3:12,1"`. Of course, this works for
            any arbitrary setting such as `"3:5,2:5,10,3:12"` (which specifies five triplicates, followed by two duplicates, a single decaplicate, and twelve triplicates again – truly a dataset from another dimension)...
        """
        if replicates is None:
            return self._replicates
        else: 
            # convert a string formula to tuple if one was provided
            if isinstance(replicates, str): 
                replicates = self._reps_from_formula(replicates)
            # vet replicate coverage
            if self._vet_replicates(replicates):
                self._replicates = replicates
            else: 
                aw.HardWarning("Assay:reps_dont_cover", n_samples = self._Reader.n(), reps = replicates)

    def group(self, infer_names = True):
        """
        Groups the samples according to replicates specified.

        Parameters
        ----------
        infer_names : bool
            Try to infer names of replicate groups based on the individual replicate sample identifiers.
            Note that this only works if all replicates have an identical sample name!
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
        
        if infer_names:
            # infer group names
            self._infer_names()
            

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
    
    def _reps_from_formula(self, replicates):
        """
        Generates a replicate tuple from a string formula. 
        See the docstring of `replicates()` for more info on the formula.

        Example:
        "3:4,1:4,2:3,9" -> (3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 9)
        """

        # split the formula and adjust standard formatting
        replicates = replicates.split(",")
        replicates = [i + ":1" if ":" not in i else i for i in replicates]

        # convert to numeric values and extend
        replicates = [np.array(i.split(":"), dtype = int) for i in replicates]
        replicates = [np.tile(i[0], i[1]) for i in replicates]
        
        # generate replicate tuple
        replicates = np.concatenate(replicates)
        replicates = tuple(replicates)
        
        return replicates

    def _infer_names(self):
        """
        Infers replicate group names from the given sample identifier column
        """
        # get the first group of replicates
        group0 = self._df.query("group == 0")["Sample"]
        all_identical = all(group0 == group0[0])
        
        if all_identical:
            self._df["group_name"] = self._df["Sample"]
        else: 
            aw.SoftWarning("Assay:groupnames_not_inferred")

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
        verdict = None
        
        # for INT -> modulo will be 0 if all samples are covered
        # for TUPLE -> sum(replicates) should cover all samples...

        if isinstance(replicates, int):
            verdict = True if samples % replicates == 0 else False
        elif isinstance(replicates, tuple): 
            verdict = True if sum(replicates) == samples else False
        
        if verdict is None: 
            aw.HardWarning("Assay:reps_could_not_vet", reps = replicates)

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
        
    def read(self, filename, **kwargs):
        """
        Reads one raw datafile (csv or excel format).

        Parameters
        ----------
        filename : str
            A filepath to a raw data file.
            If the file is a `csv` file, it has to have two named columns; one for replicate names, one for Ct values. 
            Both csv (`,` spearated) and csv2 (`;` separated) are accepted.
            If the file is an `excel` file it the relevant sections of the spreadsheet are identified automatically. 
            But they require identifying headers. By default `Name` and `Ct` are assumed but these can be changed using 
            the `name_label` and `Ct_label` arguments that can be passed to `read()` as kwargs.

        Returns
        -------
        Assay : qpcr.Assay
            A `qpcr.Assay` object containing the grouped and renamed data.
        """
        self._Reader = Reader(filename, **kwargs)
        self._Reader.id(aux.fileID(filename))

        self._Assay = Assay(self._Reader)
        self._Assay.adopt_id(self._Reader)

        if self._replicates is not None:
            self._Assay.replicates(self._replicates)
            self._Assay.group()
        else: 
            aw.HardWarning("SampleReader:no_reps_yet")

        if self._names is not None:
            self._Assay.rename(self._names)

        return self._Assay

class _Qupid_SampleReader(SampleReader):
    """
    Sets up a Reader+Assay pipeline that reads in a sample file and handles the 
    stored raw data in a pandas dataframe. 
    Its `read()` method directly returns a `qpcr.Assay` object that can be piped to Analyser. 
    Note
    ----
    This is the Qupid applicable version of the SampleReader
    """
    def __init__(self):
        super().__init__()
        

    def read(self, file):
        """
        Reads one raw datafile (csv format).

        Parameters
        ----------
        file
            An UploadedFile object from streamlit.

        Returns
        -------
        Assay : qpcr.Assay
            A `qpcr.Assay` object containing the grouped and renamed data.
        """
        self._Reader = _Qupid_Reader(file)
        self._Reader.id(aux.fileID(file.name)) # use the .name method to get the filename

        self._Assay = Assay(self._Reader)
        self._Assay.adopt_id(self._Reader)

        if self._replicates is not None:
            self._Assay.replicates(self._replicates)
            self._Assay.group()
        else: 
            aw.HardWarning("SampleReader:no_reps_yet")

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
        stats = []
        for col in tmp.columns:
            try: 
                stat = func(tmp[col], **kwargs)
            except: 
                stat = np.nan
            stats.append(stat)
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
            `"grouped"` (first entry for each replicate group), 
            any specified numeric value (as `float`), or a `function` that will calculate the anchor and returns a single numeric value. 
            If you wish to use a function to compute the anchor, you can access the dataframe stored by the `qpcr.Assay` that is being analysed through the 
            `data` argument. `data` will be automatically forwarded to your custom anchor-function, unless you specify it directly. Please, make sure your function can handle
            `**kwargs` because any kwargs supplied during `DeltaCt()`-calling will be passed per default to both the anchor-function and DeltaCt-function. 
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
        elif isinstance(self._anchor, (float, int)): 
            self._DeltaCt_externally_anchored(self._anchor, self._deltaCt_function, **kwargs)
        elif type(self._anchor) == type(aux.fileID):
            self._DeltaCt_function_anchor(self._anchor, self._deltaCt_function, **kwargs)

    def _DeltaCt_function_anchor(self, anchor_function, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using a function as anchor
        """
        df = self._Assay.get()
        # update a "data" argument into the kwargs for the anchor_function
        if "data" not in kwargs.keys(): 
            kwargs.update(dict(data = df))
        anchor = anchor_function(**kwargs)

        df["dCt"] = df["Ct"].apply(deltaCt_function, ref = anchor, **kwargs)
        self._Results.add(df["dCt"])



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
            group_subset = df.query(f"group == {group}").reset_index(drop = True)
            anchor = group_subset["Ct"][0]
            delta_cts = group_subset["Ct"].apply(deltaCt_function, ref=anchor, **kwargs)
            dCt = dCt.append(delta_cts).reset_index(drop = True)
        dCt.name = "dCt"
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
    files = ["Example Data 2/28S.csv", "Example Data 2/actin.csv", "Example Data 2/HNRNPL_nmd.csv", "Example Data 2/HNRNPL_prot.csv"]

    groupnames = ["wt-", "wt+", "ko-", "ko+"]

    analysers = []

    reader = SampleReader()
    reader.replicates("6:4")
    # reader.names(groupnames)

    analyser = Analyser()

    def myanchor(data):
        """
        computes a custom anchor
        """
        df = data.query("group == 0")["Ct"].reset_index(drop=True)
        return df.mean()

    analyser.anchor(myanchor)

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