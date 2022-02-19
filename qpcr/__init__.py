"""
This module is designed to provide functions to analyse qPCR data. 
It is designed for maximal user-friendliness and streamlined data-visualisation.

## Terminology
---

Let's talk some terminology first. There are a number of important terms that you will meet when using the `qpcr` module that you may feel have nothing to do with qPCR. 
Read on to learn what these terms mean and why they are important. 

### `qpcr` vs qPCR
Throughout this documentation, we will refer to the python module as `qpcr` in all lowercase, while we will talk about the qPCR experiment itself 
by using `qPCR` with uppercase. 

### File (or "datafile")
Let's start simple. A `file` is simply one of your input datafiles (shocker). For the `qpcr` module this means either a `csv` or an `excel` file. 
Files are the at the very basis of our data pipeline. 
However, `qpcr` does not strictly require that these files correspond to anything in particular from "qPCR vocabulary" _per se_ 
(although the basic assumption underlying the default settings
is that a single "regular" datafile contains exactly one single qPCR assay). 

The main `qpcr` classes are designed to work with datafiles that contain Ct values and replicate identifiers (see below for "replicates").
In the very most basic case, a simple datafile contains extactly an `"id"` or `"Name"` or whatever column and a `"Ct"` column (can also be named differently),
that store values from exactly _one single qPCR assay_. 

However, if your experimental setup looks differently, there are ways to adapt `qpcr` to handle different setups.
But if your data follows the above mentioned standard arrangements, you can think of a file as identical to a "qPCR assay".

### "Regular" vs "irregular" datafiles
`Regular` datafiles only contain the `id` and `Ct` columns and nothing else, storing values of exactly one qPCR assay. 
`Irregular` datafiles also have to have these two columns, but they allow for more irregular stuff
around these columns. 
Irregular files may contain multiple datasets / assays (see below). Check out the documentation of the `qpcr.Parsers` and `qpcr.Readers` for 
details on how to work with irregular and multi-assay datafiles.

### An `qpcr.Assay` and a "dataset" / "assay"
Let's start by the term "dataset". Well, a "dataset" is simply a collection or replicate identifiers and corresponding Ct values belonging together. 
By default this usually corresponds to a single qPCR assay, hence you will often encounter the term "assay" used instead of "dataset" to be more 
intuitive. However, the two terms are really interchangeable. The `qpcr.Assay` class is used to store the Ct values and replicate identifiers 
from a single dataset extracted from a datafile. 

### Replicates
As far as the `qpcr` module is concerned, each row within your datasets corresponds to one replicate. 
Hence, a _replicate_ is just a single pair of some identifier and a corresponding Ct value. 
That means that `qpcr` does not terminologically distinguish between multiplets (i.e. "true replicates" as an experimentor would understand them) and unicates. 
If you feel more comfortable with a term like "measurement" then you can also think of the replicates like that. 

### Groups (of Replicates)
This is one of the most important terms. In the previous paragraphs we already started to talk about replicates "belonging together". 
Groups of replicates are, as their name already implies, well, a number of replicates that somehow do belong together.
For most cases, if your datafiles contain one assay each, then the groups of replicates are most likely your different qPCR samples / experimental conditions. 
We talk about _groups_ of replicates instead of samples or conditions, primarily, because there might be different data setups so that these terms might not be always appropriate.
If your data follows default arrangements, however, then a _group of replicates_ is just what you would think of as a "qPCR sample". 
Groups are assigned a numeric index starting from 0, which is how they are identified by the classes. 
However, they also come with a text label called the `group_name` (you can manually set and re-set the group names as you like). 
Many classes such as the `qpcr.DataReader` will actually just use the term `names` instead of the full `group_names`. 
Whenever you see anything "names"-related it is (super-duper most likely) a reference to the `group_names`.

### Specifying (Groups of) Replicates
Here's the best part: usually, we don't necessarily need to do anything because `qpcr.Assay` objects are able to infer the groups of replicates in your data 
automatically from the replicate identifiers (yeah!). However, you will be asked to manually provide replicate settings in case this fails. 
In case you want to / have to manually specify replicate settings, an `qpcr.Assay` accepts an input `replicates` which is where you can specify this information. 

`replicates` can be either an `integer`, a `tuple`, or a `string`. Why's that? Well, normally we perform experiments as "triplicates", or "duplicates", or whatever multiplets.
Hence, if we always have the same number of replicates in each group (say all triplicates) we can simply specify this number as `replicates = 3`. 
However, some samples might only be done in unicates (such as the diluent sample), while others are triplicates.

In these cases your dataset does not have uniformly sized groups of replicates and a single number will not do to describe the groups of replicates. 
For these cases you can specify the number of replicates in each group separately as a `tuple` such as `replicates = (3,3,3,3,1)` or as a `string` "formula"
which allows you to avoid repeating the same number of replicates many times like `replicates = "3:4,1"`, which will translate into the same tuple as we specified manually. 
Check out the documentation of `qpcr.Assay.replicates` for more information on this. 

### `Delta-Ct` vs `Delta-Delta-Ct` vs `normalisation`
The default analysis workflow in $\Delta \Delta Ct$ analysis is to first calculate a $\Delta Ct$ using an intra-assay reference and then calculate the $\Delta \Delta Ct$ using a normaliser assay. 
The first $\Delta Ct$ step is performed by a class called `qpcr.Analyser` using its native method `qpcr.Analyser.DeltaCt`. 
Why just `DeltaCt` and not `DeltaDeltaCt`? Well, we call this second `Delta`-step in `DeltaDeltaCt` differently. 
We name it `normalisation`, and it is handled by a class called `qpcr.Normaliser` using its native method `qpcr.Normaliser.normalise`. 
So, as far as the `qpcr` module is concerned there is only `qpcr.Analyser.DeltaCt` which performs the first $\Delta Ct$, and `qpcr.Normaliser.normalise` which later handles the second "delta"-step to get to $\Delta \Delta Ct$. 
Of course, this means that the `qpcr.Normaliser` will need to have knowledge about which `qpcr.Assay` objects contain actual assays-of-interest and which ones contain normaliser-assays (specifying that is easy, though, so don't worry about that).

> Please, note at this point that, as described in more detail in the documentation of the `qpcr.Analyser`, Delta-Ct values are directly computed as 
> exponential values $ \Delta Ct' = 2^{ - \Delta Ct}$, while normalisation later performs $ \mathrm{norm. } \Delta\Delta Ct = \\frac{  \Delta Ct'_s  }{  \Delta Ct'_n  }$, where $s$ is an assay of interest's Delta-Ct ($\Delta Ct'$) value of some replicate, 
> and $n$ is the corresponding value of the normaliser assay. This is based on the mathemathical equivalence of $n^{  a - b  } \equiv \\frac{  n^{ a } } {  n^{ b } }$. 
> Hence, while the documentation will continuously use the terms $\Delta Ct$ and $\Delta\Delta Ct$, they are in fact the exponential deriviative of the conventional values.

### The `anchor` and the "reference group"
Next to the "groups of replicates", this is probably one of the most important terms. The `anchor` is simply the intra-dataset reference used by the `qpcr.Analyser` to perform its first $\Delta Ct$. 
If your datafiles contain one assay each, and your groups of replicates are your qPCR samples, then you will likely have some "wildtype", "untreated", or "control" sample, right? 
Well, in `qpcr` terms that would be your _reference group_.
Usually your `anchor` is part of or generated from the Ct values of your _reference group_ (like their `mean` for instance).
By default it is assumed that your reference group is the _very first_ group of replicates. However, it's not a big problem if this is not the case, as you can specify different anchors easily.
So, again, the `anchor` is the dataset-internal reference value used for the first $\Delta Ct$.

### "assays" vs "normalisers"
You will likely encounter methods and/or arguments that speak of "assays" and "normalisers", especially with the `qpcr.Normaliser`. 
For all intents and purposes, an "assay" is simply one of your datasets (we know this already).
However, in practice "assays" are the short notation for specifically "assays-of-interest" 
(or more formally "datasets-of-interest"), while "normalisers" refer to your normaliser-assays (from housekeeping genes like ActinB for instance). 
But again, if your datafiles do not conform to standard data arrangements, do not be distracted from the terminology here.

You will also find that the term "assays" is used within the final results dataframe (when using the summary-statistics mode). 
In this setting "assays" refers to the assay-of-interst whose data was analysed according to the provided normaliser-assays. 
In fact, this is a new "hybrid" assay identifier taht includes the names of all the normaliser-assays used during computation (check out what the final results look like and it'll be immediately clear).

### "Samples"
You may find that there is also a term "sample" within `qpcr`'s vocabulary. 
As far as the `qpcr` module is concerned, the term "sample" is not very important in itself and usually appears in the context of "sample assays".
In this setting it is used interchangeably with "assays-of-interest". 
Actually, we try to phase out the term "sample" and it currently mainly appears in hidden auxiliary functions which have retained the term from earlier development versions.


## Some more basics 
----

### `pipeline`s 
A `pipeline` is essentially any workflow that starts from one or multiple input datafiles and ultimately pops out some results table you are happy with.
Pipelines can be manually created by assembling the main `qpcr` classes, usually starting with a Reader, passing to an Analyser, to an Normaliser, and you're good to go.
When manually assembling your workflow you can extract your data at any point and perform your own computations on it as you like. However, if you wish to "just do some good ol' Delta-Delta-Ct"
there are pre-defined pipelines that will handle writing the workflow and only require a very basic setup. You can find these in the `qpcr.Pipes` submodule.

### "Results" = my final results?
Yes and no. Anything that is computed through any of the `qpcr` classes is called a "result" of some kind.
In practice as soon as you pass your data through a `qpcr.Normaliser` you generate "results" which are actually stored in a separate class called `qpcr.Results` (that's not so important, though).
So, when using the term "result" we usually mean just anything that was explicitly computed by `qpcr`. 

### `get`ting your data
Too many classes and objects? Well, no worries, the underlying data is stored as `pandas DataFrames`. To get your data from the clutches of the `qpcr` classes you can always use the `get()` method. 
`get` is pretty universal in the `qpcr` module, so whenever you want to extract your data, there's a `get()` method to help you.

### `link` vs `add` vs `pipe`
Different classes have slightly different methods of adding data to them. Classes that only accept one single data input (such as a single `qpcr.Assay` object or a single filepath as `string`)
usually have a `link()` method that, well, links the data to them. After that the classes are ready to perform whatever actions they can perform (an `qpcr.Analyser` would perform `DeltaCt()` for instance). 
Some classes such as the `qpcr.Analyser` have a wrapper that will call both their `link()` as well as their actual core-functional method together in one go. This wrapper is called `pipe()`. 
So for the `qpcr.Analyser` you could either manually use `link()` and then `DeltaCt()`, or simply call `pipe()` which does both for you. It is noteworthy that `pipe` methods actually _return_ whatever
their output is, which is *not* normally the case otherwise (normally you'd use the `get()` method to extract your data, see above). Most `qpcr.Readers` and `qpcr.Parsers` are also equipped with `pipe` methods.
Alright, we know about `link` and `pipe` now, what about `add`?  Classes that accept multiple inputs have `add` methods, which tells the class where exactly to store the input data. 
`add`-methods are especially implemented within the pre-defined analysis pipelines of the `qpcr.Pipes` submodule. You will probably often use the methods `add_assays()` and `add_normalisers()` if you plan on using these predefined pipelines.
However, these classes usually still have a `link()` method somewhere that you can use as well. For instance, the `qpcr.Normaliser` uses a `link` method to add both assays-of-interest and normaliser-assays simultaneously.

## Getting started
----
You can find useful tutorials and applied examples as `jupyter notebooks` [on GitHub](https://github.com/NoahHenrikKleinschmidt/qpcr/tree/main/Examples).

"""

import pandas as pd
import qpcr._auxiliary as aux
from qpcr._auxiliary import warnings as aw
import qpcr._auxiliary.defaults as defaults
import qpcr.Parsers as Parsers
import qpcr.Readers as Readers
import os
import numpy as np 
from copy import deepcopy 
from io import StringIO
import re

__pdoc__ = {
    "_CORE_Reader" : True
}

# default column names for raw Ct data files (don't change this!)
raw_col_names = defaults.raw_col_names

supported_filetypes = defaults.supported_filetypes

default_group_name = defaults.default_group_name

# the default columns that are structurally part of an Assay object
ref_cols = [ raw_col_names[0], "group", "group_name", defaults.default_dataset_header ]


# At some future point we will remove the _CORE_Reader from the 
# __init__ as it's now part of the Readers submodule...
class _CORE_Reader(aux._ID):
    """
    The class handling the core functions of the Reader class. 
    The standard qpcr.Reader inherits from this. 

    Note
    ----
    This implementation is deprecated. The _CORE_Reader has been moved
    (and further updated since) to the `qpcr.Readers` submodule. 

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
            The dataframe from the datafile.
        """
        return self._df

    def n(self):
        """
        Returns
        -------
        n : int
            The number of replicates (entries) in the dataframe.
        """
        return len(self._df[raw_col_names[0]])


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
                # setup parser
                parser = Parsers.CsvParser()
                # check it file should be read transposed
                transpose = aux.from_kwargs("transpose", False, kwargs, rm = True)
                if transpose:
                    parser.transpose()
                
                # setup patterns and store assay-of-interest
                assay_pattern = aux.from_kwargs("assay_pattern", "Rotor-Gene", kwargs)
                assay_of_interest = aux.from_kwargs("assay", None, kwargs, rm=True)
                parser.assay_pattern(assay_pattern)
                
                # get data column labels
                id_label = aux.from_kwargs("id_label", "Name", kwargs, rm = True)
                ct_label = aux.from_kwargs("ct_label", "Ct", kwargs, rm = True)
                parser.labels(id_label,ct_label)
                
                # pipe the datafile through the parser
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
            # setup parser
            parser = Parsers.ExcelParser()
            # check it file should be read transposed
            transpose = aux.from_kwargs("transpose", False, kwargs, rm = True)
            if transpose:
                parser.transpose()
            
            # check for sheet_name
            sheet_name = aux.from_kwargs("sheet_name", 0, kwargs, rm = True)

            # setup patterns and store assay-of-interest
            assay_pattern = aux.from_kwargs("assay_pattern", "Rotor-Gene", kwargs)
            assay_of_interest = aux.from_kwargs("assay", None, kwargs, rm=True)
            parser.assay_pattern(assay_pattern)
            
            # get data column labels
            id_label = aux.from_kwargs("id_label", "Name", kwargs, rm = True)
            ct_label = aux.from_kwargs("ct_label", "Ct", kwargs, rm = True)
            parser.labels(id_label,ct_label)

            # pipe the datafile through the parser
            parser.read(self._src, sheet_name = sheet_name)
            parser.parse(**kwargs)

            if len(parser.assays()) > 1:
                if assay_of_interest is None: 
                    aw.HardWarning("Reader:cannot_read_multifile", file = self._src, assays = parser.assays(), traceback = False)
                self._df = parser.get(assay_of_interest)
                self.id(assay_of_interest)
            else:
                assay_of_interest = parser.assays()[0]
                self._df = parser.get(assay_of_interest)
                self.id(assay_of_interest)
        
    def _csv_read(self, **kwargs):
        """
        Reads the given data file if it's a csv file
        """
        df = None
        try: 
            df = pd.read_csv(
                                self._src, 
                                sep = self._delimiter, 
                                header = aux.from_kwargs("header", 0, kwargs), 
                                names = raw_col_names
                            )
        except: 
            aw.HardWarning("Reader:cannot_read_csv", file = self._src)

        # check if a valid Ct column was found
        Ct = raw_col_names[1]
        full_valid_Ct_col = len(  df[ df[Ct] == df[Ct] ]  ) == len(df)
        if not full_valid_Ct_col:
            aw.HardWarning("Reader:cannot_read_csv", file = self._src)

        self._df = df

    def _filesuffix(self):
        """
        Returns the file-suffix of the provided file
        """
        try: 
            suffix = self._src.split(".")[-1]
        except: 
            pass
        return suffix

class Reader(_CORE_Reader):
    """
    Reads qpcr raw data files in csv or excel format to get a single dataset. 

    Note
    -----
    This implementation of the default Reader is now deprecated and will be removed
    at some point in the future. The new impelentation is the `qpcr.Readers.SingleReader`.
    Please, use that one instead.

    Input Data Files
    ----------------
    Valid input files are either regular `csv` files, or  irregular `csv` or `excel` files, 
    that specify assays by one replicate identifier column and one Ct value column.

    Irregular input files may specify multiple assays as separate tables, 
    one assay has to be selected using the `assay` argument. 
    Separate assay tables may be either below one another (separated by blank lines!)
    or besides one another (requires `transpose = True`).

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


# class _Qupid_Reader(_CORE_Reader):
#     """
#     This Reader class works with streamlit's UploadedFile class.
    
#     Note
#     -------
#     We have to use a little hack to make the UploadedFile readable.
#     We convert to a string and then back to a StringIO which we can then pass to pandas...
    
#     Parameters
#     ----------
#     file
#         A streamlit UploadedFile object
#     """
#     def __init__(self, file) -> pd.DataFrame: 
#         super().__init__()
#         self._filename = file.name
#         if self._filesuffix() == "csv":
#             self._content = file.read().decode()
#             self._src = StringIO(self._content)
#             self._delimiter = ";" if self._is_csv2() else ","
#         else: 
#             self._src = file
#         self.read()

#     def _filesuffix(self):
#         """
#         Returns the file-suffix of the provided file
#         """
#         suffix = self._filename.split(".")[-1]
#         return suffix

#     def _is_csv2(self):
#         """
#         Tests if csv file is ; delimited (True) or common , (False)
#         """
#         if ";" in self._content: 
#             return True
#         return False

#     def _has_header(self):
#         """
#         Checks if column headers are provided in the data file
#         It does so by checking if the second element in the first row is numeric
#         if it is numeric (returns None << False) no headers are presumed. Otherwise
#         it returns 0 (as in first row has headers)...
#         """
#         content = self._content.split("\n")[0]
#         content = content.split(self._delimiter)
#         try: 
#             second_col = content[1]
#             second_col = float(second_col)
#         except ValueError:
#             return 0 # Headers in row 0
#         return None  # no headers

class Assay(aux._ID):
    """
    The central storing unit of single datasets that were read from datafiles.
    An `qpcr.Assay` stores the replicate identifiers and Ct values, and also 
    groups these according to the `replicates` information (which is automatically
    inferred by default). Groups of replicates can be arbitrarily renamed by the user.

    Note
    -------
    The new implementation of the `qpcr.Assay` works directly with a DataFrame
    that was generated by any one of the `qpcr.Readers` or `qpcr.Parsers` 
    instead of a (now depcrecated) `qpcr.Reader` object. 
    However, a `qpcr.Reader` can still be passed as `df` argument and will be read-in using `link`. 
    Support for this will be removed at some point in the future. 

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame produces by one of the `qpcr.Readers` 
        containing an `id` column for the replicate identifiers 
        and a `Ct` value column. 
    id : str
        The identifer of the assays (the Assay name, essentially). 

    replicates : int or tuple or str
            Can be an `integer` (equal group sizes, e.g. `3` for triplicates), 
            or a `tuple` (uneven group sizes, e.g. `(3,2,3)` if the second group is only a duplicate). 
            Another method to achieve the same thing is to specify a "formula" as a string of how to create a replicate tuple.
            The allowed structure of such a formula is `n:m,` where `n` is the number of replicates in a group and `m` is the number of times
            this pattern is repeated (if no `:m` is specified `:1` is assumed). See `qpcr.Assay.replicates` for an example. 

    group_names : list
        A list of names to use for the replicates groups. If replicates of the same group share the same identifier, then the 
        group will be inferred automatically. Otherwise, default group names will be set if no `group_names` are provided. 
    """
    def __init__(self, df : pd.DataFrame = None, id : str = None, replicates : (int or tuple or str) = None, group_names : list = None) -> dict:
        super().__init__()
        
        # FUTURE DROP HERE
        # check if a qpcr.Reader was supplied instead of a dataframe
        # drop this at some point
        if isinstance(df, Reader):
            self.link(df)
        else: 
            self._df = df
            if id is not None: self._id = id

        # setup length of the found data        
        self._length = None if self._df is None else len(self._df)

        # get replicates
        self._replicates = replicates

        # store names 
        self._names = group_names

        # if we got data, try to read it 
        if self._df is not None: 
            try: 
                self.replicates(self._replicates)
                self.group()
            except Exception as e:
                aw.SoftWarning("Assay:setup_not_grouped")
            
            # and try to change names, provided that we could group yet...
            if self._names is not None and self.groups() is not None: 
                self.rename(self._names)

    def save(self, filename : str):
        """
        Saves the data from the `Assay` to a `csv` file.
        Parameters
        ----------
        filename : str
            The filename into which the assay should be stored.
            If this is a `directory`, then the assay `id` will automatically
            be used as filename. 
        """
        if os.path.isdir(filename):
            filename = os.path.join(filename, f"{self.id()}.csv")
        self.to_csv(filename, index = False)

    def get(self):
        """
        Returns
        -------
        data : pandas.DataFrame
            The stored dataframe
        """
        return self._df

    def Ct(self):
        """
        Returns
        ------
        Ct : pandas.Series
            A pandas Series with the assay's Ct values. The column is renamed 
            from "Ct" to the assay's `id`.
        """
        Ct = self._df[ raw_col_names[1] ]
        Ct.name = self.id()
        return Ct


    def dCt(self):
        """
        Returns
        -------
        dCt : pandas.Series
            A pandas Series with the computed Delta-Ct values. The column is renamed 
            from "dCt" to the assay's `id`.
        """
        dCt = self._df["dCt"]
        dCt.name = self.id()
        return dCt

    def ddCt(self):
        """
        Returns
        -------
        ddCt : pandas.DataFrame
            A pandas DataFrame with all Delta-Delta-Ct values that the Assay has stored. 
            All `"rel_{}"` columns are renamed to include the assay `id` to `"{id}_rel_{}"`.
        """
        # get all ddCt columns
        ddCt = [ i for i in self._df.columns if "rel_" in i ]
        id = self._id
        # make new names and generate renaming dictionary 
        new_names = [ f"{id}_{i}" for i in ddCt ]
        new_names = {  old : new for new, old in zip(new_names, ddCt)  }

        # get the data and rename
        ddCt = self._df[ ddCt ]
        if not isinstance(ddCt, pd.DataFrame):
            ddCt = pd.DataFrame(ddCt)
        ddCt = ddCt.rename(columns = new_names)
        
        return ddCt

    # FUTURE FEATURE HERE
    # def fc(self):
        # some method to also return the fold change columns... 


    def rename_cols(self, cols:dict):
        """
        Renames columns according to a dictionary as key -> value.

        Parameters
        ----------
        cols : dict
            A dictionary specifying old column names (keys) and new colums names (values).
        """
        self._df = self._df.rename(columns = cols)


    def link(self, Reader:Reader):
        """
        Links a `qpcr.Reader` object to the Assay.

        Note
        ------
        This is deprecated since the `qpcr.Reader` has been 
        replaced by the `qpcr.Readers.SingleReader`. 
        This method will be removed in the future.  

        Parameters
        ----------
        Reader : qpcr.Reader
            A qpcr.Reader object.
        """
        self._Reader = Reader
        self.adopt_id(Reader)
        df = self._Reader.get()
        self._df = df
        self._length = self._Reader.n()
        
    def n(self):
        """
        Returns 
        ------

        int 
            The number of entries (individual replicates) within the Assay.
        """
        return self._length

    def add_dCt(self, dCt : pd.Series): 
        """
        Adds results from Delta-Ct (first Delta-Ct performed by a `qpcr.Analyser`).

        Parameters
        -----------
        dCt : pandas.Series
            A pandas Series of Delta-Ct values that will be stored in a column `"dCt"`.
            Note, that each `Assay` can, of course, only store one single Delta-Ct column. 
        """
        self._df["dCt"] = dCt
    
    def add_ddCt(self, normaliser_id : str, ddCt : pd.Series):
        """
        Adds results from Delta-Delta-Ct ("normalisation" performed by a `qpcr.Normaliser`).
        These will be stored in a column named `"rel_{normaliser_id}"`. Hence, an Assay can store
        an arbitrary number of Delta-Delta-Ct columns against an arbitrary number of different normalisers. 
        
        Parameters
        ----------
        normaliser_id : str
            The id of the normaliser Assay used to compute the Delta-Delta-Ct values.
        ddCt : pandas.Series
            A pandas Series of Delta-Delta-Ct values.
        """
        name = f"rel_{normaliser_id}"
        self._df[name] = ddCt

    # FUTURE FEATURE HERE
    # some method to add fc columns here...

    def adopt(self, df : pd.DataFrame, force = False):
        """
        Adopts an externally computed dataframe as its own.
        This is supposed to be used when setting up new `qpcr.Assay` objects that do not 
        inherit data from one of the `qpcr.Readers`. If you wish to alter an existing `qpcr.Assay` use `force = True`.
        When doing this, please, make sure to retain the proper data structure!

        Parameters
        ----------
        df : pd.DataFrame
            A pandas DataFrame.
        force : bool
            If a dataframe is already stored the new dataframe will only be stored if `force = True`.
        """
        if self._df is None: 
            self._df = df
        elif force:
            self._df = df
        else:
            aw.HardWarning("Assay:no_data_adopted")
        self._length = len(self._df)

    def names(self, as_set = True):
        """
        Returns a set of the replicate group names (maintaing group order).

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
        if "group_name" in self._df.columns: 
            if as_set:
                return aux.sorted_set(list(self._df["group_name"]))
            else: 
                return self._df["group_name"]
        else: 
            aw.SoftWarning("Assay:no_groupname_assignment")
            return None
    
    def groups(self, as_set = True):
        """
        Returns a set of sample groups (numeric).

        Parameters
        ----------
        as_set : bool
            If `as_set = True` (default) it returns a set (as list without duplicates) 
            of assigned group names for replicate groups.
            If `as_set = False` it returns the full group_name column (including all repeated entries).
        
        Returns
        -------
        groups : list
            The given numeric group identifiers of all replicate groups.
        """
        if "group" in self._df.columns:
            groups = sorted(list(set(self._df["group"]))) if as_set else self._df["group"]
            return groups
        else:
            aw.SoftWarning("Assay:setup_not_grouped")
            return None

    def replicates(self, replicates : (int or tuple or str) = None):
        """
        Either sets or gets the replicates settings to be used for grouping
        Before they are assigned, replicates are vetted to ensure they cover all data entries.

        Parameters
        ----------
        replicates : int or tuple or str
            Can be an `integer` (equal group sizes, e.g. `3` for triplicates), 
            or a `tuple` (uneven group sizes, e.g. `(3,2,3)` if the second group is only a duplicate). 
            Another method to achieve the same thing is to specify a "formula" as a string of how to create a replicate tuple.
            The allowed structure of such a formula is `n:m,` where `n` is the number of replicates in a group and `m` is the number of times
            this pattern is repeated (if no `:m` is specified `:1` is assumed). 
            
            So, as an example, if there are 12 groups which are triplicates, but
            at the end there is one which only has a single replicate (like the commonly measured diluent qPCR sample), we could either specify the tuple
            individually as `replicates = (3,3,3,3,3,3,3,3,3,3,3,3,1)` or we use the formula to specify `replicates = "3:12,1"`. Of course, this works for
            any arbitrary setting such as `"3:5,2:5,10,3:12"` (which specifies five triplicates, followed by two duplicates, a single decaplicate, and twelve triplicates again – truly a dataset from another dimension)...
        """
        if replicates is not None and self._df is not None: 
            # convert a string formula to tuple if one was provided
            if isinstance(replicates, str): 
                replicates = self._reps_from_formula(replicates)
            # vet replicate coverage
            if self._vet_replicates(replicates):
                self._replicates = replicates
            else: 
                aw.HardWarning("Assay:reps_dont_cover", n_samples = self._length, reps = replicates)
        return self._replicates

    def group(self, infer_names = True):
        """
        Groups the data according to replicates-settings specified.

        Parameters
        ----------
        infer_names : bool
            Try to infer names of replicate groups based on the individual replicate sample identifiers.
            Note that this only works if all replicates have an identical sample name!
        """
        
        # generate group and group_names columns
        if isinstance(self._replicates, int):
            groups, group_names = self._make_equal_groups()            
        elif isinstance(self._replicates, tuple):
            groups, group_names = self._make_unequal_groups()
        else:
            if self._identically_named():
                groups = self._infer_replicates()
                group_names = [default_group_name.format(i) for i in groups]
            else: 
                aw.HardWarning("Assay:no_reps_inferred", assay = self.id())
        
        # add numeric group identifiers
        self._df["group"] = groups
        self._df["group_name"] = group_names
        
        if infer_names: #and self._names is None:
            # infer group names
            self._infer_names()
            

    def rename(self, names:(list or dict)):
        """
        Replaces the current names of the replicate groups 
        (stored in the "group_name" column).

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
            aw.SoftWarning("Assay:no_groupname_assignment", names = names)

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

    def _infer_replicates(self):
        """
        Infers the replicate groups based on the replicate ids in case all replicates of the same group have the same name.
        """
        names = self._df[raw_col_names[0]]
        names_set = aux.sorted_set(names)
        groups = [i for i in range(len(names_set))]
        for name, group in zip(names_set, groups):
            names = names.replace(name, group)
        
        indices = np.array(names, dtype = int)
        return indices

    def _infer_names(self):
        """
        Infers replicate group names from the given replicate identifier column
        """
        if self._identically_named():
            self._df["group_name"] = self._df[raw_col_names[0]]
        elif self._names is None: 
            aw.SoftWarning("Assay:groupnames_not_inferred")

    def _identically_named(self):
        """
        Checks if all replicates in the same group have the same name / id
        It checks simply the first group, if that is identical then it's fine.
        """
        if "group" not in self._df.columns:
            names = self._df[raw_col_names[0]]
            names_set = aux.sorted_set(names)
            first_name = names_set[0]
            group0 = self._df.query(f"{raw_col_names[0]} == '{first_name}'")[raw_col_names[0]]
            entries = len(group0)
            all_identical = entries > 1                
        else: 
            group0 = self._df.query("group == 0")[raw_col_names[0]]
            all_identical = all(group0 == group0[0])
        return all_identical

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
        to cover all data entries.
        (this function works with a tuple for replicate group sizes)
        """
        groups = []
        group_names = []
        for rep, idx in zip(self._replicates, range(len(self._replicates))): 
            groups.extend([idx] * rep)
            group_names.extend([ default_group_name.format(idx) ] * rep)
        return groups, group_names

    def _make_equal_groups(self):
        """
        Returns two lists of [0,0,0,1,1,1] and 
        [Group0, Group0, Group0, Group1,...] 
        to cover all data entries.
        (this function works with an integer group size, 
        assuming all groups have the same size)
        """
        assays = self._length
        groups = []
        group_names = []
        slices = range(int(assays / self._replicates))
        for i in slices:
            groups.extend([i] * self._replicates)
            group_names.extend([ default_group_name.format(i) ] * self._replicates)
        return groups, group_names

    def _vet_replicates(self, replicates : (int or tuple)):
        """
        Checks if provided replicates will place all data entries into a group
        returns True if all replicates are covered, False if not...
        """
        current_entries = self._length
        verdict = None

        # for INT -> modulo will be 0 if all replicates are covered
        # for TUPLE -> sum(replicates) should cover all replicates...

        if isinstance(replicates, int):
            verdict = True if current_entries % replicates == 0 else False
        elif isinstance(replicates, tuple): 
            verdict = True if sum(replicates) == current_entries else False
        
        if verdict is None: 
            aw.HardWarning("Assay:reps_could_not_vet", reps = replicates)

        return verdict

"""
## Setting up a `qpcr.Assay`
---

Here is a manual example of creating a `qpcr.Assay` object. You can use either the `qpcr.DataReader` or any one of `qpcr.Readers` directly to 
read in your data and generate a pandas DataFrame. Note, the `qpcr.Readers` are already equipped with `make_Assay(s)` methods that will handle
setting up `qpcr.Assay` objects for you. 

However, setting up `qpcr.Assay`s manually can be as simple as:

```python
# get the dataframe from one of the qpcr.Readers
mydata = some_reader.get()

assay = Assay( df = mydata, id = "my_assay" )

```

If your replicate identifiers are the same for all replicates within each group then the groups are automatically inferred. And your assay is 
ready at this point already to be passed to an `qpcr.Analyser`. If not, you can specify the replicates manually like: 

```python
# manually specify triplicates during setup
assay = Assay( df = mydata, id = "my_assay", replicates = 3 )

# or you can change the replicates after initial setup like 
assay = Assay( df = mydata, id = "my_assay" )
assay.replicates(3)
assay.group()

```

"""

class SampleReader(Assay):
    """
    Sets up a Reader+Assay pipeline that reads in a single datafile and handles the 
    extracted dataset in a pandas dataframe. 
    Its `read()` method directly returns a `qpcr.Assay` object that can be piped to Analyser. 
    
    Note
    ----
    This is now deprecated and will be removed in a future version! Please, use the `qpcr.DataReader` instead.
    """
    def __init__(self):
        super().__init__()
        
        aw.SoftWarning("Versions:Deprecation", old = "SampleReader", new = "DataReader")

        self._replicates = None
        self._names = None
        self._Reader = None
        self._Assay = None

    # def __str__(self):
    #     header = f"qpcr.SampleReader ({self._id})"
    #     reps = f"Replicate settings:\t{self._replicates}"
    #     names = f"Name settings:\t\t{self._names}"

    #     header_line = max([len(i) for i in [header, reps, names]])
    #     header_line = "-" * header_line

    #     string = f"{header_line}\n{header}\n{header_line}\n\n{reps}\n{names}\n\n{header_line}"
    #     return string

    def replicates(self, replicates:(int or tuple)):
        """
        Set the replicates specifics to use for grouping.

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
        **kwargs
            Any additional keyword arguments that should be passed to the `qpcr.Reader`.

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
        else: 
            pass 
            # VITAL CHANGE HERE
            # since Assays can now infer replicates we 
            # don't raise an Error if no replicates are specified manually...
            # aw.HardWarning("SampleReader:no_reps_yet")

        if self._names is not None:
            self._Assay.group(infer_names = False)
            self._Assay.rename(self._names)
        else:
            self._Assay.group()
            
        return self._Assay

# class _Qupid_SampleReader(SampleReader):
#     """
#     Sets up a Reader+Assay pipeline that reads in a datafile and handles the 
#     stored raw data in a pandas dataframe. 
#     Its `read()` method directly returns a `qpcr.Assay` object that can be piped to Analyser. 
#     Note
#     ----
#     This is the Qupid applicable version of the SampleReader
#     """
#     def __init__(self):
#         super().__init__()
        

#     def read(self, file):
#         """
#         Reads one raw datafile (csv format).

#         Parameters
#         ----------
#         file
#             An UploadedFile object from streamlit.

#         Returns
#         -------
#         Assay : qpcr.Assay
#             A `qpcr.Assay` object containing the grouped and renamed data.
#         """
#         self._Reader = _Qupid_Reader(file)
#         self._Reader.id(aux.fileID(file.name)) # use the .name method to get the filename

#         self._Assay = Assay(self._Reader)
#         self._Assay.adopt_id(self._Reader)

#         if self._replicates is not None:
#             self._Assay.replicates(self._replicates)
#             self._Assay.group()
#         else: 
#             aw.HardWarning("SampleReader:no_reps_yet")

#         if self._names is not None:
#             self._Assay.rename(self._names)

#         return self._Assay

class DataReader(aux._ID):
    """
    Handles reading a single file containing input data
    for `qpcr`. 
    
    Note
    -----
    This is a top-level class that is designed
    as the central port through which data is read into `qpcr.Assay` objects
    from both regular and irregular, single- and multi-assay files.
    This is the suggested way to read your data for most users, which should
    work in most cases.
    
    However, due to the automated setup of the inferred Readers there may be cases where you 
    will either have a hard time or be unable to read your datafiles using the `DataReader`. 
    In such cases, don't try too long to make it work with the DataReader, 
    just use one of the `qpcr.Readers` or even `qpcr.Parsers`directly.  

    """
    def __init__(self):
        super().__init__()
        self._replicates = None
        self._names = None
        self._Reader = None             # the functional core will be either a Reader
        self._Data = {}                 # the _Data attribute will store any output from the Reader as a dictionary with filename : data structure.
        self._tmp_data = None           # this will not attempt to distinguish between assays / normalisers, or anything. It's just an archive of whatever data we got.
                                        # by default data is not stored by the DataReader, it's read only pipes through, but data can be stored using the .store() method.

    def Reader(self, Reader = None):
        """
        Gets or sets the functional core Reader    
        """
        if Reader is not None:
            self._Reader = Reader
        return self._Reader

    def clear(self):
        """
        Clears all data that was extracted
        """
        self._Data = {}
    
    def reset(self):
        """
        Resets the core Reader
        """
        self._Reader = None
    
    def prune(self):
        """
        Resets the DataReader completely
        """
        self.__init__()

    def store(self):
        """
        Will store the read data. 


        Note 
        -------
        `DataReader` does NOT have a specific data storage facility to distinguish between 
        assays / normalisers, data types, etc. It simply keeps a dictionary of `{filename : data}`
        that can be accessed. This is designed in case multiple files should be read using the same 
        `DataReader` to allow an easier access of the data in case the data outputs are of the same type.

        However, the main intended application of `DataReader` is to use the`read` method's returned data directly.
        """
        self._Data.update(self._tmp_data)
    
    def get(self):
        """
        Returns the stored data
        """
        return self._Data

    def read(self, filename : str, multi_assay : bool = False, big_table : bool = False, decorator : (bool or str) = None, reset = False, **kwargs):
        """
        Reads an input file and extracts available datasets using the
        specified `Reader` or by setting up an approproate `Reader`. 

        Parameters
        ----------
        filename : str
            A filepath to an input datafile.

        multi_assay : bool
            Set to `True` if the file contains multiple assays you wish to read.

        big_table : bool
            Set to `True` if the file is a "Big Table" file. 
            Check out the documentation of the `qpcr.Readers` for more 
            information on "Big Table" files.

        decorator : str or bool
            Set if the file is decorated. This can be set either to `True` for `multi_assay` and multi-sheet (excel) or `big_table` files,
            or it can be set to a valid `qpcr decorator` for single assay files or single-sheet files.
            Check out the documentation of the `qpcr.Parsers` for more information on decorators.

        reset : bool
            If multiple input files should be read but they do not all 
            adhere to the same filetype / datastructure, use `reset = True` 
            to set up a new Reader each time `read` is called.

        **kwargs
            Any additional keyword arguments to be passed to the core Reader.
            Note, while this tries to be utmost versatile there is a limitation
            to costumizibility through the kwargs. If you require streamlined datareading
            use dedicated `qpcr.Readers` and/or `qpcr.Parsers` directly.
        """
        self._src = filename
        # vet filesuffix
        suffix = self._filesuffix()
        if suffix not in supported_filetypes:
            aw.HardWarning("MultiReader:unknown_datafile", file = self._src)
        
        if reset or self._Reader is None: 
            self.reset()
            self._setup_Reader(
                                multi_assay = multi_assay, 
                                big_table = big_table,
                                decorator = decorator,
                                **kwargs
                            )
        
        # read file and return data
        data = self._Reader._DataReader( 
                                            filename = self._src, 
                                            decorator = decorator, 
                                            **kwargs 
                                    )

        self._tmp_data = {self._src : data}
        return data

    def _filesuffix(self):
        """
        Returns the filesuffix
        """
        return self._src.split(".")[-1]

    def _setup_Reader(self, **kwargs):
        """
        Sets up the core Reader 
        """
        suffix = self._filesuffix()

        use_multi = aux.from_kwargs("multi_assay", False, kwargs, rm = True)
        is_bigtable = aux.from_kwargs("big_table", False, kwargs, rm = True)

        if suffix == "csv":
            if not use_multi and not is_bigtable:
                reader = Readers.SingleReader()
            elif is_bigtable:
                reader = Readers.BigTableReader()
            else:
                reader = Readers.MultiReader()

        elif suffix == "xlsx":
            
            # check if not only a single sheet should be read
            multi_sheet = "sheet_name" not in kwargs

            if use_multi and multi_sheet:
                reader = Readers.MultiSheetReader()
            elif use_multi and not multi_sheet:
                reader = Readers.MultiReader()
            elif is_bigtable:
                reader = Readers.BigTableReader()
            else:
                reader = Readers.SingleReader()
        self._Reader = reader


class Results(aux._ID):
    """
    Handles a pandas dataframe for data and computed results from a `qpcr` class. 
    
    Note
    -----
    This is a central data collection that can inherit directly from `qpcr.Assay` objects and from 
    extrenally computed sources. Please, note that it will not perform extensive vetting on its data input, 
    so make sure to only provide proper data input when manually assembling your `qpcr.Results`!
    """
    def __init__(self):
        super().__init__()
        self._df = None
        self._Assay = None
        self._stats_results = {
                                "group" : [], 
                                defaults.default_dataset_header : [], 
                                "mean" : [], "stdev" : [], "median" : []
                            }
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
        if self.is_empty():
            self._df = Assay.get()
            self._drop_setup_cols()
        else:
            named_identically = all( self._df["group_name"] == Assay.get()["group_name"] )
            if not named_identically:
                aw.SoftWarning("Results:cannot_link")


    def add_Ct(self, assay : Assay):
        """
        Adds a `"Ct"` column with Delta-Ct values from an `qpcr.Assay`.
        It will store these as a new column using the Assay's `id` as header.

        Parameters
        -------
        assay : qpcr.Assay
            An `qpcr.Assay` object from which to import.
        """
        id = assay.id()
        Ct = assay.get()[ raw_col_names[1] ]
        Ct.name = id
        self.add(Ct)

    def add_dCt(self, assay : Assay):
        """
        Adds a `"dCt"` column with Delta-Ct values from an `qpcr.Assay`.
        It will store these as a new column using the Assay's `id` as header.

        Parameters
        -------
        assay : qpcr.Assay
            An `qpcr.Assay` object from which to import.
        """
        dCt = assay.dCt()
        self.add(dCt)

    def add_ddCt(self, assay : Assay):
        """
        Adds all `"rel_{}"` columns with Delta-Delta-Ct values from an `qpcr.Assay`.
        It will store these as new columns using the Assay's `id` + the `_rel_{}` composite id.

        Parameters
        -------
        assay : qpcr.Assay
            An `qpcr.Assay` object from which to import.
        """
        df = assay.ddCt()
        # add data
        self.add(df)

    # FUTURE FEATURE HERE
    # this feature will be for fold-change columns of second normalisation
    # currently this is not yet supported but is supposed to be an implemented feature of the 
    # Normaliser. The Normaliser is currently perfectly able to perform second normalisation but
    # does not yet have an integrated method of properly storing the results and there is also not
    # # yet a method for pairing assays together for second normalisation.
    # def add_fc(self, assay : Assay):
    #     """
    #     Adds all `"fc_{}"` columns with Delta-Delta-Ct values from an `qpcr.Assay`.
    #     It will store these as new columns using the Assay's `id` + the `_fc_{}` composite id.

    #     Parameters
    #     -------
    #     assay : qpcr.Assay
    #         An `qpcr.Assay` object from which to import.
    #     """
    #     id = assay.id()
    #     df = assay.get()
    #     # get the ddCt containing columns
    #     rel_cols = [ i for i in df.columns if "fc_" in i ]
    #     # generate new composite ids 
    #     new_names = [ f"{id}_{i}" for i in rel_cols ]
    #     # get ddCt columns 
    #     df = df[rel_cols]
    #     # rename the columns to include the assay id
    #     df = df.rename( columns = { new : old for new, old in zip(new_names, rel_cols) } )

    #     # add data
    #     self.add(df)

    def names(self, as_set = False):
        """
        Returns 
        -------
        names : list or None
            The adopted `group_names` 
            (only works if a `qpcr.Assay` has already been linked 
            using `adopt_names()`!)
        """
        if self._df is not None:
            names = self._df["group_name"]
            if as_set: names = aux.sorted_set(names)
            return names
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

    def drop_groups(self, groups : (list or str)):
        """
        Removes specific groups of replicates from the DataFrame.

        Parameters
        ----------
        groups : list
            Either the numeric group identifiers or the group names
            of the groups to be removed, or a `regex` pattern defining which groups
            should be dropped (this is useful for systematically removing RT- groups etc.)
        """
        # check for regex pattern
        # and get corresponding group names 
        if isinstance(groups, str):
            groups = [i for i in self._df["group_name"] if re.match(groups, i) is not None]
        
        # get the right reference column and query to use to be 
        # used (either group or group_name)
        ref_query = "group != {group}" if isinstance( groups[0], int ) else "group_name != '{group}'"
        
        # remove groups from dataset
        for group in groups: 
            self._df = self._df.query(ref_query.format(group = group))

            # also drop from stats df
            if self._stats_df is not None:
                self._stats_df = self._stats_df.query(ref_query.format(group = group))
    

    def add(self, column:pd.Series, replace : bool = False):
        """
        Adds some new column of data.

        Note
        ----
        The `column` argument has to be named for this to work. However, there are 
        already implemented methods dedicated to adding specifically Delta-Ct, Delta-Delta-Ct or just
        Ct values to the Results.

        Parameters
        ----------
        column : pd.Series
            A named pandas Series or DataFrame that can be joined into the already
            stored dataframe.
        replace : bool
            In case results from a computation with the same identifiers are already stored
            no new data can be stored under that id. Either the new data must be renamed or
            `replace = True` must be set to overwrite the presently stored data. 
        """
        if isinstance(column, pd.Series):
            if column.name in self._df.columns:
                if not replace:  
                    aw.SoftWarning("Results:name_overlap", name = column.name)
                else: 
                    self._df[column.name] = column
            else: 
                self._df = self._df.join(column)
        else: 
            for i in column.columns:
                if i in self._df.columns: 
                    if not replace: 
                        aw.SoftWarning("Results:name_overlap", name = i )
                col = column[i]
                self._df[i] = col

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

            # get only the delta-delta-Ct columns
            cols = [i for i in R_df.columns if i not in ref_cols]
            R_df = R_df[cols]


            # we merge the dataframes first without adding 
            # some new id suffix, only do so if this fails
            try: 
                new_df = pd.merge(new_df, R_df, 
                                    right_index = True, left_index = True, 
                                )
            except: 
                new_df = pd.merge(new_df, R_df, 
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
            If this is the case then the only columns retained are: `"group", "group_name", "id", "assay"`.
        """
        if cols == ():
            _to_drop = [c for c in self._df.columns if c not in [ "group", "group_name", raw_col_names[0], defaults.default_dataset_header ]]
        else:
            _to_drop = [c for c in cols if c in list(self._df.columns)]
        self._df = self._df.drop(columns = _to_drop)
        
    def rename_cols(self, cols:dict):
        """
        Renames columns according to a dictionary as key -> value.

        Parameters
        ----------
        cols : dict
            A dictionary specifying old column names (keys) and new colums names (values).
        """
        self._df = self._df.rename(columns = cols)


    def stats(self, recompute = False) -> pd.DataFrame:
        """
        Computes summary statistis about the replicate groups: 
        `Mean`, `Median`, and `StDev` of all replicate groups, for all datasets (assays).
        
        Parameters
        ----------
        recompute : bool
            Statistics will only be once unless recompute is set to `True`.

        Returns
        -------
        stats_df : pd.DataFrame
            A new dataframe containing the computed statistics for each replicate group.

        """
        default_dataset_header = defaults.default_dataset_header
        # if stats_df is already present, return but sorted according to assays, not groups (nicer for user to inspect)
        if self._stats_df is not None and not recompute:
            return self._stats_df.sort_values(default_dataset_header)
        elif recompute: 
            self._stats_results = {"group" : [], default_dataset_header : [], "mean" : [], "stdev" : [], "median" : []}
            self._stats_df = None

        # get groups and corresponding assay columns 
        groups = aux.sorted_set(list(self._df["group"]))
        assays = [c for c in self._df.columns if c not in ref_cols]
     
        # compute stats for all replicates per group
        for group in groups:
            group_subset = self._df.query(f"group == {group}")
            
            median = self._stat_var(group_subset, np.nanmedian)
            mean = self._stat_var(group_subset, np.nanmean)
            stdv = self._stat_var(group_subset, np.nanstd)
            self._add_stats(assays, group, median, mean, stdv)
            
        # add group names
        self._add_stats_names(assays)

        self._stats_df = pd.DataFrame(self._stats_results)
        return self._stats_df.sort_values(default_dataset_header)

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
            # in case of raw results export we don't need the "assay" column as all 
            # assays are stored as separate columns anyaway, so it doesn't store any useful data
            _df = self._df
            if "assay" in _df.columns: _df = self._df.drop( columns = ["assay"] )
            self._save_single(path, _df, "_df")
        if stats:
            # compute stats if none have been computed yet...
            if self._stats_df is None:
                self.stats()
            self._save_single(path, self._stats_df, "_stats")

    def drop_rel(self):
        """
        Crops the `X_rel_Y` column-names of Delta-Delta-Ct results to just `X`.
        I.e. reduces back to the assay-of-interest name only.
        """
        colnames = self._df.columns
        to_change = {i : i.split("_rel_")[0] for i in colnames if "_rel_" in i }
        self.rename_cols(to_change)

        # also recompute the stats df with new names...
        if self._stats_df is not None: 
            self.stats(recompute = True)

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
        shared_columns = [i for i in self._df.columns if i in ref_cols]
        dct_columns = [i for i in self._df.columns if i not in ref_cols]
        
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
        src.to_csv(filename, index = False)
        
    def _drop_setup_cols(self):
        """
        Removes unnnecessary columns from the df during self._df setup with link()
        """
        # drop the Ct column
        self.drop_cols(
                        raw_col_names[1], "dCt"
                    )


    def _add_stats_names(self, samples):
        """
        Adds a group_name column to self._stats_result with appropriate
        repetition of group_names for each group of replicates...
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
        ignore = [raw_col_names[0], "group", "group_name", "assay"]
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
    Performs Single Delta-Ct (first normalisation 
    within dataset against the `anchor`) 
    """
    def __init__(self):
        super().__init__()
        self._Assay = None

        # default settings
        self._anchor = "first"
        self._ref_group = 0
        self._ref_group_col = "group" # used in case of "mean" anchor where the ref_group must be located either from a numeric (group) or string (group_name) id
        
        self._efficiency = 1                # the formal effiency in percent
        self._eff = 2 * self._efficiency    # the actual doubplciation factor used for calculation
        self._deltaCt_function = self._get_deltaCt_function(exp = True)

    def get(self):
        """
        Returns 
        -------
        Assay : qpcr.Assay
            The analysed `qpcr.Assay` object that contains now deltaCT values.
        """
        return self._Assay

    def link(self, Assay:Assay):
        """
        Links a `qpcr.Assay` object to the Analyser.
       
        Parameters
        ----------
        Assay : qpcr.Assay
            A `qpcr.Assay` object containing data.
        
        """
        self._Assay = Assay

    def pipe(self, Assay:Assay, **kwargs) -> Results:
        """
        A quick one-step implementation of link + DeltaCt.

        Note
        ----
        This is the suggested application of the `qpcr.Analyser`.


        Parameters
        ----------
        Assay : qpcr.Assay
            A `qpcr.Assay` object to be linked to the Analyser for DeltaCt computation.
        
        **kwargs
            Any additional keyword arguments to be passed to the `DeltaCt()` method.

        Returns 
        -------
        assay : qpcr.Assay
            The same `qpcr.Assay` with computed Delta-Ct values. 

        """
        self.link(Assay)
        self.DeltaCt(**kwargs)
        assay = self.get()
        return assay

    def efficiency(self, e:float = None):
        """
        Sets an efficiency factor for externally calculated qPCR amplification efficiency.
        By default `efficiency = 1` (100%) is assumed.

        Parameters
        ----------
        e : float
            An amplification efficiency factor. Default is `e = 1`, 
            which is then treated as `eff = 2 * e`, so `e = 1` corresponds to true duplication
            each cycle.

        Returns
        -------
        efficiency : float
            The current efficiency used.

        """
        if isinstance(e, (int, float)):
            self._efficiency = float( e )
            self._eff = 2 * self._efficiency
        return self._efficiency

    def anchor(self, anchor : (str or float or function) = None, group : (int or str) = 0):
        """
        Sets the anchor for DeltaCt for internal normalisation.

        Parameters
        ----------
        anchor : str or float or function
            The internal anchor for normalisation.
            This can be either `"first"` (default, the very first dataset entry),
            `"mean"` (mean of the reference group), 
            `"grouped"` (first entry for each replicate group), 
            any specified numeric value (as `float`), 
            or a `function` that will calculate the anchor and returns a single numeric value. 
            If you wish to use a function to compute the anchor, you can access the dataframe stored by the `qpcr.Assay` that is being analysed through the 
            `data` argument. `data` will be automatically forwarded to your custom anchor-function, unless you specify it directly. Please, make sure your function can handle
            `**kwargs` because any kwargs supplied during `DeltaCt()`-calling will be passed per default to both the anchor-function and DeltaCt-function. 
        group : int or str
            The reference group identifier. This can be either the numeric identifier or the `group_name`. This is only used for `anchor = "mean"`. 
            By default the first group is assumed. 
        
        Returns
        -------
        anchor 
            The currently selected `anchor`.
        ref_group
            The current reference group.
        """
        if anchor is not None:
            self._anchor = anchor
        if group != self._ref_group:
            self._ref_group = group
        
        # update column where to search for ref_group identifier
        if isinstance(group, str):
            self._ref_group_col = "group_name"
        else: 
            self._ref_group_col = "group"
        
        return self._anchor, self._ref_group


    def func(self, f:(str or function)):
        """
        Sets the function to be used for DeltaCt (optional)

        Parameters
        ----------
        f : str or function
            The function to be used for DeltaCt computation. Pre-defined functions are 
            either `"exponential"` (which uses  `dCt = eff ** ( -(s-r) )`, default), or `"linear"` 
            (uses `dCt = s-r`), where `s` is any replicate entry in the dataframe and `r` is the anchor. `eff = 2 * efficiency` is the 
            numeric duplication factor (default assumed `efficiency = 1`).
            It is also possible to assign any defined function that accepts one `float` Ct value `s` (1st!) and anchor `r` value (2nd!), 
            alongside any `kwargs` (which will be forwarded from DeltaCt()...). It must return also a single `float`. 
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
        Calculates Delta-Ct for all groups within the dataframe.
        Any specifics such as `anchor` or `func` must have already been 
        set using the respective methods prior to calling `DeltaCt()`!

        Parameters
        ----------
        **kwargs
            Any additional keyword arguments that a custom DeltaCt function may require.
        """

        predefined = {
                        "first"     : self._DeltaCt_first_anchored,
                        "grouped"   : self._DeltaCt_grouped_anchored,
                        "mean"      : self._DeltaCt_mean_anchored,
                    }

        if self._anchor in predefined.keys():
        
            deltaCt_func = predefined[ self._anchor ]
            deltaCt_func(
                                            self._deltaCt_function, 
                                            **kwargs
                                    )
        
        elif isinstance(self._anchor, (float, int)): 
            self._DeltaCt_externally_anchored(
                                                self._anchor, 
                                                self._deltaCt_function, 
                                                **kwargs
                                            )
        elif type(self._anchor) == type(aux.fileID):
            self._DeltaCt_function_anchor(
                                            self._anchor, 
                                            self._deltaCt_function, 
                                            **kwargs
                                        )

    def _DeltaCt_function_anchor(self, anchor_function, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using a function as anchor
        """
        df = self._Assay.get()
        # update a "data" argument into the 
        # kwargs for the anchor_function
        if "data" not in kwargs.keys(): 
            kwargs.update(dict(data = df))
        anchor = anchor_function(**kwargs)

        # apply deltaCt_function
        Ct = raw_col_names[1]
        dCt = df[Ct].apply(deltaCt_function, r = anchor, **kwargs)
        dCt.name = "dCt"
        # store results
        self._Assay.add_dCt(dCt)

    def _DeltaCt_mean_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the mean of the reference group as anchor
        """
        df = self._Assay.get()

        # get  reference group
        ref_query = "{} == {}" if isinstance(self._ref_group, int) else "{} == '{}'" 
        ref = df.query(
                        ref_query.format(  self._ref_group_col, self._ref_group  )
                    )
        # get Ct values from ref group and make anchor
        ref = ref[raw_col_names[1]]
        anchor = ref.mean()
        
        # apply DeltaCt function
        Ct = raw_col_names[1]
        dCt = df[Ct].apply(deltaCt_function, r = anchor, **kwargs)
        dCt.name = "dCt"
        # store results
        self._Assay.add_dCt(dCt)


    def _DeltaCt_externally_anchored(self, anchor:float, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using a specified anchor
        """
        # get Ct column label 
        Ct = raw_col_names[1]
        df = self._Assay.get()

        # apply DeltaCt function
        dCt = df[Ct].apply(deltaCt_function, r = anchor, **kwargs)
        dCt.name = "dCt"
        # store results
        self._Assay.add_dCt(dCt)


    def _DeltaCt_grouped_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the first entry of each group as anchor
        """
        # get set of sample groups and dataset
        groups = self._Assay.groups()
        df = self._Assay.get()

        # get Ct column label
        Ct = raw_col_names[1]

        dCt = pd.Series()
        for group in groups: 
            group_subset = df.query(f"group == {group}").reset_index(drop = True)
            anchor = group_subset[Ct][0]
            delta_cts = group_subset[Ct].apply(deltaCt_function, r = anchor, **kwargs)
            dCt = dCt.append(delta_cts).reset_index(drop = True)
        
        dCt.name = "dCt"

        # store results
        self._Assay.add_dCt(dCt)

    def _DeltaCt_first_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the very first entry of the dataset as anchor
        """
        # get Ct column label
        Ct = raw_col_names[1]
        
        # get first available entry
        # we do this instead of just 0 
        # because the truly first entry 
        # might have been filtered out 
        df = self._Assay.get()
        first = list(df.index)[0]

        # get anchor
        anchor = df[Ct][ first ]

        # apply DeltaCt function
        dCt = df[Ct].apply(deltaCt_function, r = anchor, **kwargs)
        dCt.name = "dCt"
        # store results
        self._Assay.add_dCt(dCt)

    def _exp_DCt(self, s, r, **kwargs):
        """
        Calculates deltaCt exponentially
        """
        factor = s - r 
        return self._eff **(-factor)

    def _simple_DCt(self, s, r, **kwargs):
        """
        Calculates deltaCt linearly
        """
        return s - r

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
    Handles the second step in Delta-Delta-Ct (normalisation against normaliser assays).
    
    Note
    -----
    This requires that all have been analysed in the same way before!
    """
    def __init__(self):
        super().__init__()

        self._Normalisers = []
        self._Assays = []
        
        self._Results = Results()
        
        # the actually used normaliser 
        # which will be a pre-processed version
        # from all supplied normalsier assays
        # by default this will be the averaged version
        # of all normaliser assays (but this can be changed)  
        self._normaliser = None

        # setup defaults
        self._prep_func = self._preprocess_normalisers
        self._norm_func = self._divide_by_normaliser

    def clear(self):
        """
        Will clear the presently stored results
        """
        self._Results = Results()

    
    def prune(self, assays = True, normalisers = True, results = True):
        """
        Will clear assays, normalisers, and/or results
        assays : bool
            Will clear any sample assays if True (default).
        
        results : bool
            Will clear any computed results if True (default).
        
        normalisers : bool
            Will clear any normalisers if True (default).
        """
        if assays:
            self._Assays = []
        if normalisers: 
            self._Normalisers = []
        if results: 
            self.clear()


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
    
    def link(self, assays:(list or tuple or Analyser) = None, normalisers:(list or tuple or Analyser) = None):
        """
        Links either normalisers or assays-of-interest `qpcr.Assay` objects coming from the same `qpcr.Analyser`.

        Parameters
        ----------
        assays : list or tuple or qpcr.Analyser
            A list of `qpcr.Assay` objects coming from a `qpcr.Analyser` or the `qpcr.Analyser` itself. These assays will be normalised against a normaliser.
        
        normalisers : list or tuple or qpcr.Analyser
            A list of `qpcr.Assay` objects coming from a `qpcr.Analyser` or the `qpcr.Analyser` itself. These assays will be used as normalisers. These will be
            combined into one single pseudo-normaliser which will then be used to normalise the assays. The method of 
            combining the normalisers can be specified using the `qpcr.Normaliser.prep_func` method.
        """
        # convert to lists (since the _link methods really want lists)
        if assays is not None and not isinstance(assays, (list, tuple)): 
            assays = [assays]
        self._link_assays(assays)

        if normalisers is not None and not isinstance(normalisers, (list, tuple)): 
            normalisers = [normalisers]
        self._link_normaliser(normalisers)
    
    def prep_func(self, f = None):
        """
        Sets any defined function for combined normaliser pre-processing.
        If no `f` is provided, it returns the current `prep_func`.

        Parameters
        ----------
        f : function
            The function may accept one list of `qpcr.Assay` objects, and must return 
            either an `qpcr.Assay` object directly or a `pandas.Dataframe` (that will be migrated to an `qpcr.Assay`).
            The returned dataframe must contain a `"dCt"` column which stores the delta-Ct values ultimately used as 
            "normaliser assay".  
        """
        
        if aux.same_type(f, aux.fileID):
            self._prep_func = f
        elif f is None:
            return f
        else: 
            aw.HardWarning("Normaliser:cannot_set_prep_func", func = f)

    def norm_func(self, f = None):
        """
        Sets any defined function to perform normalisation of assays against normalisers.
        If no `f` is provided, it returns the current `norm_func`.

        Parameters
        ----------
        f : function
            The function may accept one `pandas.DataFrame` containing two numeric columns of delta-Ct values from a sample (named "s") and a normaliser assay (named "n"),
            as well as a group identifier column (named "group"). It must return a numeric `pandas.Series` of the same length. 
            
            By default `s/n` is used, where `s` is a column of sample-assay deltaCt values, 
            and `n` is the corresponding `"dCt"` column from the normaliser.

        """
        if aux.same_type(f, aux.fileID):
            self._norm_func = f
        elif f is None:
            return f
        else: 
            aw.HardWarning("Normaliser:cannot_set_norm_func", func = f)

    def normalise(self, **kwargs):
        """
        Normalises all linked assays against the combined pseudo-normaliser 
        (by default, unless a custom `prep_func` has been specified), 
        and stores the results in a new Results object.

        Parameters
        ----------
        **kwargs
            Any additional keyword arguments that may be passed to a custom 
            `norm_func` and `prep_func` (both will receive the kwargs!).
        """
        if self._normaliser is None: 
            self._normaliser = self._prep_func(self._Normalisers, **kwargs)
            self._vet_normaliser()

        if self._Assays == [] or self._normaliser is None:
            aw.SoftWarning("Normaliser:no_data_yet")

        # get combined dataframe
        # normaliser = self._normaliser.get()

        # setup _Results by passing in the common columns 
        # from the first Assays (since all Assays should have 
        # the same id, group, and group_name columns, if they
        # come from the same experiment...)
        self._Results.adopt_names(self._Assays[0])
        self._Results.drop_cols()

        # perform normalisation for each assay 
        for assay in self._Assays:

            # get data
            # assay_df = assay.get()

            # apply normalisation (delta-delta-Ct)
            normalised = self._norm_func_wrapper(
                                                    assay, 
                                                    self._normaliser, 
                                                    **kwargs
                                            )

            # # store results in _Results
            # self._store_to_Results(assay, normalised)

            # and store results also in the Assay itself
            assay.add_ddCt( self._normaliser.id(), normalised )

            # and store to results
            self._Results.add_ddCt(assay)

    def _vet_normaliser(self):
        """
        Checks if the normaliser is already a qpcr.Assay object, and if not
        convert it to one. 
        """
        if not isinstance(self._normaliser, Assay):
            tmp = Assay()
            tmp.adopt(self._normaliser)
            tmp.id("combined_normaliser")
            self._normaliser = tmp
            

    def _store_to_Results(self, assay, normalised):
        """
        Stores computed Delta-delta-ct (normalisation)
        into the _Results object.
        """
        column_name = f"{assay.id()}_rel_{self._normaliser.id()}"
        normalised = normalised.rename(column_name)
        self._Results.add(normalised)


    def _norm_func_wrapper(self, sample_assay, normaliser, **kwargs):
        """
        The wrapper that will apply the _norm_func to the sample and normaliser dataframes and return a pandas series of normalised values
        """
        # for double normalised we want the same columns as dct and norm...

        sample_dCt = sample_assay.dCt()
        groups = sample_assay.groups( as_set = False )
        norm_dCt = normaliser.dCt()

        tmp_df = pd.DataFrame( dict( group = groups, s = sample_dCt, n = norm_dCt )  )

        results = self._norm_func(tmp_df, **kwargs)

        # this is the old call from before factoring out to Assays 
        # dCt_col, norm_col = self._prep_columns(sample_assay, dCt_col, norm_col)

        # tmp_df = normaliser.join(sample_assay, lsuffix="_s")
        # # tmp_df = sample_assay.join(normaliser, rsuffix = "_n")
        # results = self._norm_func(tmp_df[[dCt_col, norm_col]], **kwargs)
        return results

    # NOT USED ANYMORE
    # Used to be used before factoring out data storage to the Assays 
    # def _prep_columns(self, sample_assay, dCt_col, norm_col):
    #     """
    #     Returns the columns to use if named columns shall be used 
    #     (named columns will be used for second-normalisation of entire runs)
    #     Note
    #     ----
    #     Currently, second normalisation is not yet really implemented, so this is 
    #     kinda not really used and would probably get overhauled when we 
    #     try to seriously implement that at some point. 
    #     """
    #     if dCt_col == "named":
    #         dCt_col = [i for i in sample_assay.columns if i not in ["group", "group_name", raw_col_names[0], "assay"]]
    #         # assert len(dCt_col) == 1, f"length of dCt_col is: {len(dCt_col)}"
    #         dCt_col = dCt_col[0]

    #     if norm_col == "same": 
    #         norm_col = dCt_col + "_s"
    #     return dCt_col,norm_col

    def _divide_by_normaliser(self, df, **kwargs):
        """
        Performs normalisation of sample s against normaliser n
        s and n are specified as two pandas dataframe columns
        Note, that the dataframe must ONLY contain these two columns, first the dCt sample, then the normaliser!
        (default _norm_func)
        """
        s, n = df["s"], df["n"]
        # this is the old call from before factoring out to Assays
        # dCt_col, norm_col = df.columns
        # s, n = df[dCt_col], df[norm_col]
        return s / n

    def _link_assays(self, assays):
        """
        Links any provided assays and checks their datatype in the process...
        """
        if assays is not None:
            for sample in assays: 
                if aux.same_type(sample, Assay()):
                    self._Assays.append(sample)

                elif aux.same_type(sample, Analyser()):
                    self._Assays.append(sample.get())
                
                else: 
                    aw.SoftWarning("Normaliser:unknown_data", s = sample)
                
    def _link_normaliser(self, normalisers):
        """
        Checks if normaliser is provided and has proper datatype to be added...
        """
        if normalisers is not None:
            for normaliser in normalisers:

                if aux.same_type(normaliser, Assay()):
                    self._Normalisers.append(normaliser)

                elif aux.same_type(normaliser, Analyser()):
                    self._Normalisers.append(normaliser.get())

                else: 
                    aw.SoftWarning("Normaliser:norm_unknown_data", s = normaliser)

    def _preprocess_normalisers(self, *args, **kwargs):
        """
        Averages the provided normalisers row-wise for all normalisers into a 
        single combined normaliser, that will be stored as a new Assay object.
        """

        # initialise new Results to store the dCt values form all normalisers
        combined = Results()
        
        # setup names using the first normaliser
        combined.adopt_names(self._Normalisers[0])
        combined.drop_cols("dCt")
        combined.adopt_id(self._Normalisers[0])

        # now add all dCt columns from all normalisers
        for norm in self._Normalisers:
            combined.add_dCt(norm)
        
        # remove the non-dCt columns as they would interfere with
        # pre-processing. But we keep the "group" column because some
        # custom prep_func may want to use the group references.
        # i.e. we just replace any "Ct" columns that may have smuggled in if 
        # a non-default process is used for assembly
        combined.drop_cols( raw_col_names[1] )

        # now generate the combined normaliser
        combined_normaliser = self._average(combined)
        combined_normaliser = combined_normaliser.rename("dCt")
        combined.add(combined_normaliser)
        
        # now assemble the normaliser into a qpcr.Assay
        combined = combined.get()
        normaliser = Assay()
        normaliser.adopt(combined)
        normaliser.adopt_id(self._Normalisers[0])


        self._normaliser = normaliser  
        if len(self._Normalisers) > 1:
            self._update_combined_id()
  
        # forward combined_id to self and _Results 
        self.adopt_id(self._normaliser)
        self._Results.adopt_id(self._normaliser)

        return self._normaliser

    def _update_combined_id(self):
        """
        Generates a new id based on all normaliser ids,
        joining them as a+b+c,...
        """
        ids = [N.id() for N in self._Normalisers]
        ids = "+".join(ids)
        self._normaliser.id_reset()
        self._normaliser.id(ids)
        

    def _average(self, combined):
        """
        Averages row-wise all Normaliser entries and 
        generates a series of their per-row means
        (default preprocess_normalisers function)
        """
        tmp_df = combined.get()

        # drop group as it is a numeric column 
        # and would otherwise skew the average
        if "group" in tmp_df.columns:
            tmp_df = tmp_df.drop(columns = ["group"]) 
        tmp_df = tmp_df.mean(axis = 1, numeric_only = True)

        return tmp_df


if __name__ == "__main__":
    
    files = ["./Examples/Example Data/28S.csv", "./Examples/Example Data/actin.csv", "./Examples/Example Data/HNRNPL_nmd.csv", "./Examples/Example Data/HNRNPL_prot.csv"]
    # files = ["Example Data 2/28S.csv", "Example Data 2/actin.csv", "Example Data 2/HNRNPL_nmd.csv", "Example Data 2/HNRNPL_prot.csv"]

    groupnames = ["wt-", "wt+", "ko-", "ko+"]


    reader = DataReader()

    analyser = Analyser()
    analyser.anchor("first")

    normaliser = Normaliser()

    assays = []
    for file in files: 

        assay = reader.read(file, replicates = "6,6,6,6", names = groupnames)
        assay = analyser.pipe(assay)
        assays.append(assay)

    # first to files are normalisers
    normaliser.link(normalisers = assays[:2])

    # last to files are assays
    normaliser.link(assays = assays[2:])

    # def some_prepfunc(normalisers):
    #     first = normalisers[0].get()
    #     first = first["dCt"]
    #     return pd.DataFrame(dict(dCt_combined = first))
    
    # normaliser.prep_func(some_prepfunc)

    normaliser.normalise()

    print(assays[2].id())
    print(assays[2].get())
    
    result = normaliser.get()
    print(result.stats())



    # alternatively we could link the analyser directly 
    # to get the Assay from there like

    # for file in files[2:]:

    #     assay = reader.read(file, replicates = 6)
    #     analyser.pipe(assay)
    #     normaliser1.link(assays = analyser)

    # for file in files[:2]:

    #     assay = reader.read(file, replicates = 6)
    #     analyser.pipe(assay)
    #     normaliser1.link(normalisers = analyser)

    # normaliser1.normalise()
    
    # result1 = normaliser1.get()

    # result1.drop_groups([3, 1])
    # print(result1.stats())

    # print(result.get() == result1.get())

# down here is some perliminary trial at making second
# normalisation... But this should be properly addressed at some point...

#     splitted = result.split(reset_names = False)

#     i, j = splitted
#     i.rename_cols({"HNRNPL_nmd": "HNRNPL"})
#     j.rename_cols({"HNRNPL_prot": "HNRNPL"})

# # print(i, j)
#     # print("-------")
#     # print(i.get()["HNRNPL"] / j.get()["HNRNPL"])
#     # print("-----")
    
    
    
#     sn = Normaliser()

#     sn.link(
#         assays = [splitted[0]], 
#         normalisers = [splitted[1]],
#     )
    
#     sn.normalise(dCt_col = "named", norm_col = "same")

#     print(sn.get().get())
#     print(sn.get().stats())

#     # # result.save("..")
    
#     # #result.add_names(samples)

#     # print(result.stats())
