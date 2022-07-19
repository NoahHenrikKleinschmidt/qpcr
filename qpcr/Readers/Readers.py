"""
This module provides different `Reader` classes that allow reading simple and complex datafiles
of various architectures.

## Available Data Readers
---

### SingleReader
The `SingleReader` is able to read both regular and irregular single-assay datafiles. 
It can also read multi-assay datafiles but requires an `assay` argument, specifying
which assay specifically to extract from it.

### MultiReader
The `MultiReader` can read irregular multi-assay datafiles and extract all assays from them
either using a specific `assay_pattern` to find them or using `decorators` (check out the documentation
of the `qpcr.Parsers` for more information).

### MultiSheetReader
The `MultiSheetReader` is able to read irregular multi-assay datafiles that contain assays in multiple 
datasheets.

### BigTableReader
The `BigTableReader` is able to read datafiles that store their assays in one single "big table". It
can extract all assays from that big table using either simple extraction methods or `decorators` depending
on the type of big table (check out the documentation of the `BigTableReader` for more information on the
types of "big tables").


> ### Kwarg incompatibility Warning
> When using the `qpcr.DataReader` or a `pipe` method you will regularly observe the following warning: 
> 
> ```
> Warning:
> It appears as if some provided kwargs were incompatible with pandas.read_excel()! Defaulting to standard settings for file-reading...
> If the kwargs you specified are actually important for file reading, try manually reading and parsing to avoid kwarg incompatibilities.
> ```
>
> This is because the `pipe` method (and the `qpcr.DataReader`, which usually calls the `pipe` method of a specific Reader) pass all kwargs to both `read` and `parse`. 
> However, `pandas`' `read_excel` and `read_csv` are rather picky with the arguments they accept. So, in case you observe this warning, just know that the kwargs were removed from the 
> `read` call.
"""


# Concept to link the Readers to the DataReader
# All Readers define a _DataReader method that 
# specifies which of its methods is supposed 
# to be used for (mostly pipe, sometimes read...)
# _DataReader methods *must* return the data they read!

import pandas as pd
import qpcr
import qpcr.defaults as defaults
import qpcr._auxiliary as aux
from qpcr._auxiliary import warnings as aw
import qpcr.Parsers as Parsers
import os
import numpy as np 
from copy import deepcopy 
import re 


__pdoc__ = {
    "_CORE_Reader" : True
}

# migrate default settings from __init__
raw_col_names = defaults.raw_col_names
supported_filetypes = defaults.supported_filetypes
default_dataset_header = defaults.dataset_header
default_id_header = defaults.id_header
default_ct_header = defaults.ct_header
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
        self._header = 0
        self._df = None
        self._replicates = None
        self._names = None


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

    def make_Assay(self):
        """
        Converts the extracted dataset into an `qpcr.Assay`.
        Returns
        --------
        Assay : qpcr.Assay
            The `qpcr.Assay` from the extracted dataset.
        """
        assay = self._make_new_Assay(self.id(), self._df)
        return assay

    def read(self, **kwargs):
        """
        Reads the given data file.

        If the data file is an Excel file replicates and their Ct values will be 
        extracted from the first excel sheet of the file. Note, this assumes by default
        that the replicates are headed by the label `"Name"` and the corresponding Ct values
        are headed by the label `"Ct"`. Both labels have to be on the same row. 

        If these labels do not match your excel file, you may
        specify `id_label` and `ct_label` as additional arguments.
        """
        suffix = self._filesuffix()
        # check for a valid input file
        if suffix not in supported_filetypes:
            aw.HardWarning("MultiReader:empty_data", file = self._src)

        if suffix == "csv":
            
            # first try simple read of "regular files"
            try: 
                self._csv_read(**kwargs)
            except Exception as e:

                # users can force-regular reading mode
                is_regular = aux.from_kwargs("is_regular", False, kwargs, rm= True)
                if is_regular:
                    # print out warning
                    print(e)
                    return

                # setup parser
                parser = Parsers.CsvParser()
                self._prep_Parser(kwargs, parser)
                
                assay_of_interest = aux.from_kwargs("assay", None, kwargs, rm=True)
                
                # pipe the datafile through the parser
                parser.pipe(self._src, **kwargs)

                # get the data
                self._get_single_assay(parser, assay_of_interest)

        elif suffix == "xlsx":

            try: 
                self._excel_read(**kwargs)
            except Exception as e:

                # users can force-regular reading mode
                is_regular = aux.from_kwargs("is_regular", False, kwargs, rm= True)
                if is_regular:
                    # print out warning
                    print(e)
                    return

                # setup parser
                parser = Parsers.ExcelParser()
                self._prep_Parser(kwargs, parser)
                
                # check for sheet_name
                sheet_name = aux.from_kwargs("sheet_name", 0, kwargs, rm = True)

                # store assay-of-interest
                assay_of_interest = aux.from_kwargs("assay", None, kwargs, rm=True)
                
                # pipe the datafile through the parser
                parser.read(self._src, sheet_name = sheet_name)
                parser.parse(**kwargs)

                # get the data
                self._get_single_assay(parser, assay_of_interest)

    def names(self, names:(list or dict)):
        """
        Set names for replicates groups.

        Parameters
        ----------
        names : list or dict
            Either a `list` (new names without repetitions) or `dict` (key = old name, value = new name) specifying new group names. 
            Group names only need to be specified once, and are applied to all replicate entries.
        """
        if names is not None: 
            self._names = names
        return self._names

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
            this pattern is repeated (if no `:m` is specified `:1` is assumed). So, as an example, if there are 12 groups which are triplicates, but
            at the end there is one which only has a single replicate (like the commonly measured diluent qPCR sample), we could either specify the tuple
            individually as `replicates = (3,3,3,3,3,3,3,3,3,3,3,3,1)` or we use the formula to specify `replicates = "3:12,1"`. Of course, this works for
            any arbitrary setting such as `"3:5,2:5,10,3:12"` (which specifies five triplicates, followed by two duplicates, a single decaplicate, and twelve triplicates again – truly a dataset from another dimension)...
        """
        if replicates is not None:
            self._replicates = replicates
        return self._replicates

    def _get_single_assay(self, parser, assay_of_interest):
        """
        Gets a single dataset from the Parser
        """

        # check if there are multiple datasets
        # and if so, check if we got a specified assay_of_interest
        if len(parser.assays()) > 1:
            
            if assay_of_interest is None: 
                aw.HardWarning("Reader:cannot_read_multifile", file = self._src, assays = parser.assays(), traceback = False)
            self._df = parser.get(assay_of_interest)
            self.id_reset()
            self.id(assay_of_interest)

        # if only one assay is present anyway, get that one
        else:

            assay_of_interest = parser.assays()[0]
            self._df = parser.get(assay_of_interest)
            self.id_reset()
            self.id(assay_of_interest)

    def _prep_Parser(self, kwargs, parser):
        transpose = aux.from_kwargs("transpose", False, kwargs, rm = True)
        if transpose:
            parser.transpose()
                
        # setup patterns and store assay-of-interest
        assay_pattern = aux.from_kwargs("assay_pattern", "Rotor-Gene", kwargs)
        parser.assay_pattern(assay_pattern)
                
        # get data column labels
        id_label = aux.from_kwargs("id_label", "Name", kwargs, rm = True)
        ct_label = aux.from_kwargs("ct_label", "Ct", kwargs, rm = True)
        parser.labels(id_label,ct_label)
  
    def _csv_read(self, **kwargs):
        """
        Reads the given data file if it's a csv file

        This is the basic default reading method for 
        regular csv files.
        """
        df = None
        # header = aux.from_kwargs("header", 0, kwargs, rm  = True)
        try: 
            df = pd.read_csv(
                                self._src, 
                                sep = self._delimiter, 
                                header = self._header, 
                                # names = raw_col_names
                            )
        except: 
            aw.HardWarning("Reader:cannot_read_csv", file = self._src)


        # try to get a FileID, from the kwargs
        # if that fails, try to get the one from fileID
        id = aux.from_kwargs("id", None, kwargs, rm = True)
        if id is not None:
            self.id_reset()
            self.id(id)
        elif isinstance(self._src, str):
            self.id_reset()
            self.id(aux.fileID(self._src))
        
        # vet and crop the dataframe where necessary
        df = self._vet_single_assay_df(kwargs, df)

        self._df = df

    def _excel_read(self, **kwargs):
        """
        Reads the given data file if it's an excel file

        This is the basic default reading method for 
        regular excel files.
        """
        df = None
        sheet_name = aux.from_kwargs("sheet_name", 0, kwargs, rm = True)
        # header = aux.from_kwargs("header", 0, kwargs, rm  = True)
        try: 
            df = pd.read_excel(
                                self._src, 
                                sheet_name = sheet_name, 
                                header = self._header, 
                                # names = raw_col_names
                            )
        except: 
            aw.HardWarning("Reader:cannot_read_csv", file = self._src)


        # try to get a FileID, from the kwargs
        # if that fails, try to get the one from fileID
        id = aux.from_kwargs("id", None, kwargs, rm = True)
        if id is not None:
            self.id_reset()
            self.id(id)
        elif isinstance(self._src, str):
            self.id_reset()
            self.id(aux.fileID(self._src))

        # vet and crop the dataframe where necessary
        df = self._vet_single_assay_df(kwargs, df)

        self._df = df

    def _vet_single_assay_df(self, kwargs, df):
        """
        Vets that both Id and Ct columns are present in the data 
        and if so crops the df to the relevant columns, or checks if 
        only two columns are present anyway and then assumes Id+Ct as these two.
        """
        
        # check if we got exactly two columns only
        if len( df.columns ) == 2:

            # just get the current column names for later renaming
            Id, Ct = df.columns

        else: 

            # check if a valid Ct column was found
            Ct = aux.from_kwargs( "ct_label", default_ct_header, kwargs )
            Id = aux.from_kwargs( "id_label", default_id_header, kwargs )
            
            valid_data = Ct in df.columns and Id in df.columns
            if not valid_data:
                aw.HardWarning("Reader:cannot_find_datacols", id_label = Id, ct_label = Ct)
            else:
                # get only the relevant data columns 
                df = df[[Id, Ct]]

        # make sure to convert Ct values to float
        tmp_parser = Parsers.CsvParser()
        Ct_col = df[Ct].to_numpy()
        df[Ct] = tmp_parser._convert_to_numeric(self.id(), Ct_col)

        # rename to qpcr default headers (id + Ct)
        df = df.rename( columns = { Id : raw_col_names[0] , Ct : raw_col_names[1] }  )
        return df

    def _filesuffix(self):
        """
        Returns the file-suffix of the provided file
        """
        try: 
            suffix = self._src.split(".")[-1]
            return suffix
        except: 
            pass

    def _make_new_Assay(self, name, df):
        """
        Makes a new Assay object and performs group() already...
        """
        new_assay = qpcr.Assay(
                                df = df, 
                                id = name, 
                                replicates = self._replicates,
                                group_names = self._names
                            )
        return new_assay

class SingleReader(_CORE_Reader):
    """
    Reads qpcr raw data files in csv or excel format to get a single dataset. 

    Input Data Files
    ----------------
    Valid input files are either regular `csv` or `excel` files, or  irregular `csv` or `excel` files, 
    that specify assays by one replicate identifier column and one Ct value column.

    Irregular input files may specify multiple assays as separate tables, 
    one assay has to be selected using the `assay` argument. 
    Separate assay tables may be either below one another (separated by blank lines!)
    or besides one another (requires `transpose = True`).

    #### Example of a "regular" single-assay datafile
    |id|Ct| other data |
    |---|---| --- |
    | ctrl1| 5.67 | ... |
    | ctrl2| 5.79 | ... |
    | ctrl3 | 5.86 | ... |
    | condA1 | 5.34 | ... |
    | ... | ... | ... |


    #### Example of an "irregular" single-assay datafile
    |                     |                    |            |      |      |
    | ------------------- | ------------------ | ---------- | ---- | ---- |
    | Some meta-data here | maybe today's date |            |      |      |
    |                     |                    |            |      |      |
    | Assay 1             |                    |            |      |      |
    | id                  | Ct                 | other_data |      |      |
    | ctrl1               | 5.67               | ...        |      |      |
    | ctrl2               | 5.79               | ...        |      |      |
    | ...                 | ...                |            |      |      |


    Note
    ----
    This is the successor of the original `qpcr.Reader` (not the `qpcr.SampleReader`!).
    Hence, the `SingleReader` will return a pandas DataFrame of the dataset 
    directly using `get` but not an `qpcr.Assay`. 
    An `qpcr.Assay` object will be returned after calling `make_Assay`, however. 
    Furthermore, if the provided file cannot be read as a "regular" file the Reader will automatically
    switch to parsing. However, if your file _is_ a regular input file, you can force regular reading 
    by passing the argument `is_regular = True` to the `read` method, which will prevent parsing and allow 
    you to figure out why regular reading may have failed instead (the Reader will not 
    provide further insight into why regular reading failed if it switches to parsing).

    Parameters
    ----------
    filename : str
        A filepath to a raw data file.
        If the file is a `csv` file, it has to have two named columns; one for replicate names, one for Ct values. 
        Both csv (`,` spearated) and csv2 (`;` separated) are accepted.
        If the file is an `excel` file it the relevant sections of the spreadsheet are identified automatically. 
        But they require identifying headers. By default it is assumed that replicate identifiers and Ct values are
        stored in columns named `Name` and `Ct` but these can be changed using 
        the `id_label` and `ct_label` arguments that can be passed as kwargs. 
        Also the assay's `id` can be set as a kwarg. 

    **kwargs
        Any additional keyword arguments that shall be passed to the `read()` method which is immediately called during init.
    """
    def __init__(self, filename:str = None, **kwargs) -> pd.DataFrame: 
        super().__init__()
        self._src = filename
        self._delimiter = None
        self._header = 0
        if self._src is not None:
            self.read(**kwargs)

    def read(self, filename : str, **kwargs):
        """
        Reads the given data file.

        Note
        -----
        If the data file is an Excel file replicates and their Ct values will be 
        extracted from the first excel sheet of the file by default. 
        A separate sheet can be specified using `sheet_name`.

        Parameters
        ----------
        filename : str
            A filepath to a raw data file.
            If the file is a `csv` file, it has to have two named columns; one for replicate names, one for Ct values. 
            Both csv (`,` spearated) and csv2 (`;` separated) are accepted.
            If the file is an `excel` file it the relevant sections of the spreadsheet are identified automatically. 
            But they require identifying headers. 
            By default it is assumed that replicate identifiers and Ct values are
            stored in columns named `Name` and `Ct` but these can be changed using 
            the `id_label` and `ct_label` arguments that can be passed as kwargs. 
            Note, if only two columns are present anyway, they are assumed to be Id (1st) and Ct (2nd) column, 
            and inputs for `id_label` and `ct_label` are being ignored!
            The assay's `id` can be set as a kwarg. By default the filename is adopted as id.  
        """
        self._src = filename

        self._replicates = aux.from_kwargs("replicates", None, kwargs)
        self._names = aux.from_kwargs("names", None, kwargs)
        self._header = aux.from_kwargs("header", 0, kwargs, rm = True)

        if self._filesuffix() == "csv":
            self._delimiter = ";" if self._is_csv2() else ","
        super().read(**kwargs)

    def pipe(self, filename : str, **kwargs):
        """
        A wrapper for read+parse+make_Assay

        Returns
        -------
        assay : qpcr.Assay
            An `qpcr.Assay` object of the extracted data
        """
        self.read(filename = filename, **kwargs)
        assay = self.get()
        assay = self.make_Assay()
        return assay

    def _DataReader(self, **kwargs):
        """
        The DataReader interacting method
        """
        data = self.pipe(**kwargs)
        return data

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



# removed qpcr.Assay from inheritance here...
class MultiReader(SingleReader, aux._ID):
    """
    Reads a single multi-assay datafile and reads assays-of-interest and normaliser-assays based on decorators.
    
    Input Data Files
    ----------------
    Valid input files are multi-assay irregular `csv` or `excel` files, 
    that specify assays by one replicate identifier column and one Ct value column.

    Separate assay tables may be either below one another (separated by blank lines!)
    or besides one another (requires `transpose = True`), but ALL in the SAME sheet!

    Assays of interest and normaliser assays *must* be marked using `decorators`.

    #### Example of an "irregular" multi-assay datafile
    |                     |                    |            |      |      |
    | ------------------- | ------------------ | ---------- | ---- | ---- |
    | Some meta-data here | maybe today's date |            |      |      |
    |                     |                    |            |      |      |
    | Assay 1             |                    |            |      |      |
    | id                  | Ct                 | other_data |      |      |
    | ctrl1               | 5.67               | ...        |      |      |
    | ctrl2               | 5.79               | ...        |      |      |
    | ...                 | ...                |            |      |      |
    |                     |                    |            |   <- blank line here!   |      |
    | Assay 2             |                    |            |      |      |
    | id                  | Ct                 | other_data |      |      |
    | ctrl1               | 10.23              | ...        |      |      |
    | ctrl2               | 10.54              | ...        |      |      |
    | ...                 | ...                |            |      |      |

    Note
    ------
    `MultiReader` can transform the extracted datasets directly into `qpcr.Assay` objects using `MultiReader.make_Assays()`.
    It will perform grouping of assays if possible but will return raw-assays if not! `get` will either return a dictionary
    of the raw dataframes or a list of `qpcr.Assay`s.

    Parameters
    ----------
    filename : str
        A filepath to a raw data file, containing multiple assays that were decorated. 
        Check out the documentation of the `qpcr.Parsers`'s to learn more about decorators.
    **kwargs
            Any additional keyword arguments that should be passed to the `read` method which is immediately called during init if a filename is provided.
    """
    def __init__(self, filename : str = None, **kwargs):
        super(aux._ID, self).__init__()
        self._src = filename
        self._save_loc = None
        self._replicates = None
        self._names = None
        self._Parser = None
        self._assay_pattern = None
        self._assays = {}
        self._normalisers = {}
        if self._src is not None: 
            self._Parser = Parsers.CsvParser() if self._filesuffix() == "csv" else Parsers.ExcelParser()
            self.read(filename = self._src, **kwargs)

    def clear(self):
        """
        Clears all the extracted data from the Reader
        """
        self._assays = {}
        self._normalisers = {}

    def assays(self, which : str = None):
        """
        Parameters
        ----
        which : str
            If specified it only returns the data for the specified assay.
            Otherwise (default) it returns all assays.

        Returns
        -------
        data : dict or list
            Returns either the raw dictionary of dataframes returned by the Parser 
            (if `make_Assays` has not been run yet)
            or a list of `qpcr.Assay` objects.
        names : list
            A list of the names of all extracted assays.
        """
        return self._get_from_which(self._assays, which)

    def normalisers(self, which : str = None):
        """
        Parameters
        ----
        which : str
            If specified it only returns the data for the specified normaliser.
            Otherwise (default) it returns all normalisers.

        Returns
        -------
        data : dict or list
            Returns either the raw dictionary of dataframes returned by the Parser 
            (if `make_Assays` has not been run yet)
            or a list of `qpcr.Assay` objects.
        names : list
            A list of the names of all extracted normalisers.
        """
        return self._get_from_which(self._normalisers, which)

    def get(self, which : str):
        """
        Returns the stored assays or normalisers.

        Parameters
        ----------
        which : str
            Can be either `"assays"` or `"normalisers"` or any specific assay identifier.

        Returns
        -------
        data : dict or list
            Returns either the raw dictionary of dataframes returned by the Parser 
            (if `make_Assays` has not been run yet)
            or a list of `qpcr.Assay` objects.
        """
        data = None
        if which == "assays":
            data = self._assays
        elif which == "normalisers":
            data = self._normalisers
        else: 
            try:
                data = self._get_from_which(self._assays, which)
            except: 
                data = self._get_from_which(self._normalisers, which)
        return data

    def read(self, filename : str, **kwargs):
        """
        Reads a multi-assay datafile with decorated assays. 
        Any non-decorated assays are ignored!

        Parameters
        ----------
        filename : str
            A filepath to a raw data file, containing multiple assays that were decorated. 
            Check out the documentation of the `qpcr.Parsers`'s to learn more about decorators.
        **kwargs
                Any additional keyword arguments that should be passed to the `qpcr.Parsers`''s `read` method that extracts the datasets.
        """
        self._src = filename

        # check for a valid input file
        if self._filesuffix() not in supported_filetypes:
            aw.HardWarning("MultiReader:empty_data", file = self._src)

        self._Parser = Parsers.CsvParser() if self._filesuffix() == "csv" else Parsers.ExcelParser()

        # check if file should be read transposed
        transpose = aux.from_kwargs("transpose", False, kwargs, rm = True)
        if transpose:
            self._Parser.transpose()
        
        # setup a saving location if it was provided
        if self.save_to() is not None: 
            self._Parser.save_to(self.save_to())

        self._Parser.read(self._src, **kwargs)

    def parse(self, **kwargs):
        """
        Extracts the datasets (assays) from the read datafile.
        
        Parameters
        ----------
        **kwargs
            Any additional keyword arguments that should be passed to the `qpcr.Parsers`''s `parse` method that extracts the datasets.
        """
        # check for the two required inputs (either decorator must be specified or assay_pattern)
        # if neither is specified we default to using decorators!  
        decorator = aux.from_kwargs("decorator", None, kwargs)
        assay_pattern = aux.from_kwargs("assay_pattern", None, kwargs)

        # pass kwargs to Parser for setup
        self._prep_parser(kwargs)

        # set up parsing
        if decorator is not None or assay_pattern is None:
            self._parse_by_decorators(**kwargs)
        elif assay_pattern is not None: 
            self._parse_by_pattern(kwargs, assay_pattern)
        else: 
            # ERROR HERE
            aw.SoftWarning("MultiReader:no_decorator_or_pattern")
             
    def make_Assays(self):
        """
        Convert all found assays and normalisers into `qpcr.Assay` objects.
        """
        # convert assays to qpcr.Assay and overwrite current dict by new list
        new_assays = []
        for name, df in self._assays.items():
            new_assay = self._make_new_Assay(name, df)
            new_assays.append(new_assay)
        self._assays = new_assays

        # do the same for normalisers
        new_normalisers = []
        for name, df in self._normalisers.items():
            new_assay = self._make_new_Assay(name, df)
            new_normalisers.append(new_assay)
        self._normalisers = new_normalisers

    def pipe(self, filename :str, **kwargs):
        """
        A wrapper for read+parse+make_Assays

        Note 
        ----
        This is the suggested use of `MultiReader`. 
        If a directory has been specified into which the datafiles shall be saved, 
        then saving will automatically be done.

        Parameters
        -------
        filename : str
            A filepath to an input datafile.
        **kwargs
            Any additional keyword argument that will be passed to any of the wrapped methods.
        Returns
        -------
        data : tuple
            A tuple of the found assays-of-interst (first element) and normaliser-assays (second element).
        """
        
        # clear previously read data 
        self.clear()

        # read new data
        try: 
            self.read(filename, **kwargs)
        except: 
            self.read(filename)
            aw.SoftWarning("Parser:incompatible_read_kwargs", func = f"{type(self._Parser).__name__}'s read method")
        
        # parse and make assays
        self.parse(**kwargs)
        self.make_Assays()

        # return new data
        assays = self.get( which = "assays" )
        normalisers = self.get( which = "normalisers" )
        return assays, normalisers

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

    def _DataReader(self, **kwargs):
        """
        The DataReader interacting method
        """
        replicates = aux.from_kwargs("replicates", None, kwargs, rm = True)
        self.replicates(replicates)
        names = aux.from_kwargs("names", None, kwargs, rm = True)
        self.names(names)
        data = self.pipe(**kwargs)
        return data


    def _get_from_which(self, dataset, which):
        """
        The core for assayS() and normalisers() to get either all or a specific one
        """
        if which is not None: 
            if aux.same_type(dataset, {}):
                assay = dataset[which]
            else: 
                assay = [i for i in dataset if i.id() == which][0]
            return assay
        else:
            if aux.same_type(dataset, {}):
                names = dataset.keys()
                assays = dataset.values()
            else:
                names = [i.id() for i in dataset]
                assays = dataset
            return assays, names

    def _parse_by_pattern(self, kwargs, assay_pattern):
        """
        Parses the file only based on assay_pattern.
        Note this will also work if a decorator has been specified additionally.
        """
        self._Parser.parse( assay_pattern = assay_pattern, **kwargs )
        assays = self._Parser.get()       
        self._assays = assays

    def _parse_by_decorators(self, **kwargs):
        """
        Parses the file and idenifies assays and normalisers
        based on decorators
        """
        aux.from_kwargs("decorator", None, kwargs, rm = True)
        
        # get assays-of-interest
        self._Parser.parse( decorator = "qpcr:assay", **kwargs )
        assays = self._Parser.get()       
        self._assays = assays

        # save extracted files if so desired...
        if self.save_to() is not None: self._Parser.save()

        # clear results and run again for normalisers
        self._Parser.clear()

        # get normaliser-assays
        self._Parser.parse( decorator = "qpcr:normaliser", **kwargs )
        normalisers = self._Parser.get()
        self._normalisers = normalisers
        if self.save_to() is not None: self._Parser.save()

    def _prep_parser(self, kwargs):
        """
        Passes kwargs to Parser and performs additional setup
        """
        # setup assay_patterns if they were provided
        assay_pattern = aux.from_kwargs("assay_pattern", None, kwargs, rm = True)
        self._Parser.assay_pattern(assay_pattern)

        # check if file should be read transposed
        transpose = aux.from_kwargs("transpose", False, kwargs, rm = True)
        if transpose:
            self._Parser.transpose()

        # get data column labels
        id_label = aux.from_kwargs("id_label", "Name", kwargs, rm = True)
        ct_label = aux.from_kwargs("ct_label", "Ct", kwargs, rm = True)
        self._Parser.labels(id_label,ct_label)

class MultiSheetReader(MultiReader):
    """
    Reads a single multi-assay datafile and reads assays-of-interest and normaliser-assays based on decorators.
    
    Input Data Files
    ----------------
    Valid input files are multi-assay irregular `excel` files, 
    that specify assays by one replicate identifier column and one Ct value column.

    Separate assay tables may be either below one another (separated by blank lines!)
    or besides one another (requires `transpose = True`), but may be in DIFFERENT sheets.
    All assays from all sheets will be read!

    Assays of interest and normaliser assays *must* be marked using `decorators`.


    Parameters
    ----------
    filename : str
        A filepath to a raw data file, containing multiple assays that were decorated. 
        Check out the documentation of the `qpcr.Parsers`'s to learn more about decorators.
    **kwargs
    """
    def __init__(self):
        super().__init__()

    def read(self, *args, **kwargs):
        """
        The `MultiSheetReader` **only** offers a `pipe` method!
        Hence, neither `read` nor `parse` will work directly!
        """
        print("Sorry, the MultiSheetReader can currently only be used, through it's pipe() method!")

    def parse(self, *args, **kwargs):
        """
        The `MultiSheetReader` **only** offers a `pipe` method!
        Hence, neither `read` nor `parse` will work directly!
        """
        print("Sorry, the MultiSheetReader can currently only be used, through it's pipe() method!")



    def pipe(self, filename : str, **kwargs):
        """
        Reads a multi-assay and multi-sheet datafile with decorated assays. 
        Any non-decorated assays are ignored!

        Parameters
        ----------
        filename : str
            A filepath to a raw data file, containing multiple assays that were decorated. 
            Check out the documentation of the `qpcr.Parsers`'s to learn more about decorators.
        **kwargs
                Any additional keyword arguments that should be passed to the `qpcr.Parsers`''s `read` method that extracts the datasets.
        
        Returns
        -------
        assays : dict or list
            Returns either the raw dictionary of dataframes returned by the Parser (if no qpcr.Assays could be made automatically)
            or a list of `qpcr.Assay` objects.
        normalisers : dict or list
            Returns either the raw dictionary of dataframes returned by the Parser 
            (if no qpcr.Assays could be made automatically)
            or a list of `qpcr.Assay` objects.
        """
        self._src = filename

        # read file to get all sheets
        sheets = pd.read_excel(filename, sheet_name = None)

        all_assays = {}
        all_normalisers = {}
        
        # now repetitively read all sheets and extract data
        reader = MultiReader()
        for sheet in sheets.keys():
            try: 
                # read file and parse data
                kws = deepcopy(kwargs)
                reader.read(filename, sheet_name = sheet)
                reader.parse(ignore_empty = True, **kws)

                # get assays
                assays, normalisers = reader.get("assays"), reader.get("normalisers")
                all_assays.update(assays)
                all_normalisers.update(normalisers)

            except Exception as e: 
                # ERROR HERE
                aw.SoftWarning("MultiSheetReader:sheet_unreadable", sheet = sheets, e = e)

        # store data
        self._assays = all_assays
        self._normalisers = all_normalisers

        # try making Assays directly, return dictioanries if not possible...
        try: 
            self.make_Assays()
        except Exception as e:
            print(e)

        assays, normalisers = self._assays, self._normalisers
        return assays, normalisers
        
    def _DataReader(self, **kwargs):
        """
        The DataReader interacting method
        """
        replicates = aux.from_kwargs("replicates", None, kwargs, rm = True)
        self.replicates(replicates)
        names = aux.from_kwargs("names", None, kwargs, rm = True)
        self.names(names)
        data = self.pipe(**kwargs)
        return data

class BigTableReader(MultiReader):
    """
    Reads a single multi-assay datafile and reads assays-of-interest and normaliser-assays based on decorators.
    
    ### Input Data Files
    ----------------
    Valid input files are multi-assay irregular `csv` or `excel` files, 
    that specify assays as one big table containing all information together.
    Note that this implies that the entire data is stored in a single sheet (if using `excel` files).

    Two possible data architectures are allowed:
    
    #### `Vertical` Big Tables
    Big Tables of this kind require three columns (any additional columns are disregarded): 
    one specifying the assay, one specifying the replicate identifiers, and one specifying the Ct values. 
    An additional fourth column (`@qpcr`) may be filled with decorators but this is not necessary in this setup.

    Example:

    | assay | id   | Ct    | @qpcr |
    | ----- | ---- | ----- | ---- |
    | assay 1   | group0 | 7.65  | normaliser | 
    | assay 1   | group0 | 7.74  | normaliser | 
    | assay 1   | group0 | 7.54  | normaliser | 
    | assay 1   | group1   | 7.86  | normaliser | 
    | assay 1   | group1   | 7.57  | normaliser | 
    | assay 1   | group1   | 7.67  | normaliser | 
    | assay 2 | group0 | 16.67 | assay | 
    | assay 2 | group0 | 16.54  | assay | 
    | ...   | ...  | ...   | ... | 


    #### `Horizontal` Big Tables
    Big Tables of this kind store replicates from assays in side-by-side columns.
    The replicates may be labelled numerically or all have the same column header. 
    A second column is required specifying the replicate identifier. 

    Note, this kind of setup *requires* decorators above the first replicate of each assay,
    as well as user-defined `replicates`!

    Example:

    |      | @qpcr:group |      |      | @qpcr:group |      |      |    |
    | ---- | ---------------- | ---- | ---- | ---------------- | ---- | ---- | ---- |
    | assay   | group0_1  | group0_2  | group0_3  | group1_1 | group1_2 | ...  | @qpcr   |
    | assay1 | 7.74 | 7.65 | 7.54 | 11.54 | 11.67 | ...  |  normaliser  |
    | assay 2   | 16.67 | 16.54 | 16.97 |  16.43 |  16.56 | ...  | assay   |
    | ...  | ...  | ...  | ...  | …     | ...   | ...  |  ...  |

    > Note
    >
    > The column headers have to be **unique** to the table!
    >
    > Also, a word of warning with regard to replicate _assays_. The entries in the `assay` defining column *must* be unique! If you have multiple assays from the same gene which therefore also have the same id they will be interpreted as belonging together and will be assembled into the same `qpcr.Assay` object. However, this will result in differently sized `Assays` which will cause problems downstream when you (or a `qpcr.Normaliser`) try to assemble a `qpcr.Results` object!

    #### `Hybrid` Big Tables
    Big Tables of this kind store Ct values of different assays in separate side-by-side columns, 
    but they store the replicate identifiers as a separate column. Hence, they combine aspects of vertical and horizontal Big Tables.
    

    Example: 

    |      | @qpcr:assay| @qpcr:normaliser |  |
    | ------ | ------- | ------- | ----- |
    | id     | assay 1 | assay 2 | other_data |
    | group0 | 7.65    | 11.78   |     ...  |
    | group0 | 7.87    | 11.56   |  ...     |
    | group0 | 7.89    | 11.76   |   ...    |
    | group1 | 7.56    | 11.98   |  ...     |
    | group1 | 7.34    | 11.56   |   ...    |
    | ...    | ...     | ...     |  ...     |


    > Note
    >
    > Two options exist to read this kind of setup. 
    > - A `list` of `ct_col` values can be passed which contains the column header of each assay.
    > - The table can be `decorated`, in which case only decorated assays (columns) are extracted.
    >
    > Please, note that the two methods of reading this table are mutually exclusive! So,
    > if you decorate your table you cannot pass specific assay headers to the `ct_col` argument anymore.
    """
    def __init__(self):
        super().__init__()

        self._assays = {}
        self._normalisers = {}

        self._data = None

        self._Parser = None
        self._kind = None  # horizontal or vertical 
        self._id_col = None
        self._ct_col = None
        self._assay_col = None
        self._is_regular = False    # store if the datafile was regular and does not 
                                    # have to be converted to a dataframe based on a numpy array...
        self._hybrid_decorated = False  # because hybrid bigtables require decorator input separately for both read and parse, we store the info so it only needs to be passed during read...
    
    def pipe(self, filename : str, kind : str, id_col : str, **kwargs):
        """
        A wrapper for read+parse+make_Assays

        Note 
        -------
        This is the suggested use of the `BigTableReader`.

        Parameters
        ----------
        filename : str
            A filepath to a raw data file, containing multiple assays that were decorated. 
            Check out the documentation of the `qpcr.Parsers`'s to learn more about decorators.
        kind : str
            Specifies the kind of Big Table from the file. 
            This may either be `"horizontal"`, `"vertical"`, or `"hybrid"`.
        id_col : str
            The column header specifying the replicate identifiers 
            (or "assays" in case of `horizontal` big tables).
        **kwargs
            Any additional columns or keyword arguments.
        Returns
        -------
        assays : dict or list
            Returns either the raw dictionary of dataframes returned by the Parser (if no qpcr.Assays could be made automatically)
            or a list of `qpcr.Assay` objects.
        normalisers : dict or list
            Returns either the raw dictionary of dataframes returned by the Parser 
            (if no qpcr.Assays could be made automatically)
            or a list of `qpcr.Assay` objects.
        """
        replicates = aux.from_kwargs("replicates", None, kwargs)
        names = aux.from_kwargs("names", None, kwargs, rm = True)

        self.read(
                    filename = filename, 
                    kind = kind, 
                    id_col = id_col,
                    **kwargs
                )
        self.parse(**kwargs)

        self.replicates(replicates)
        self.names(names)
        self.make_Assays()

        assays, normalisers = self.get("assays"), self.get("normalisers")
        return assays, normalisers


    def read(self, filename : str, kind : str, id_col : str, **kwargs):
        """
        Reads a regular or irregular `csv` or `excel` datafile that contains data stored 
        in a single big table. Files are first tried to be read regularly, if this fails, 
        the Reader resorts to parsing to identify the relevant sections of the data. 

        Parameters
        ----------
        filename : str
            A filepath to a raw data file, containing multiple assays that were decorated. 
            Check out the documentation of the `qpcr.Parsers`'s to learn more about decorators.
        kind : str
            Specifies the kind of Big Table from the file. 
            This may either be `"horizontal"`, `"vertical"`, or `"hybrid"`.
        id_col : str
            The column header specifying the replicate identifiers 
            (or "assays" in case of `horizontal` big tables).
        **kwargs
            Any additional columns or keyword arguments.
        """
        self._src = filename
        self._kind = kind
        is_horizontal = self._kind == "horizontal"
        self._id_col = id_col
        self._ct_col = aux.from_kwargs("ct_col", None, kwargs, rm = True)
        self._assay_col = aux.from_kwargs("assay_col", None, kwargs, rm = True)

        if not is_horizontal:
            # first try default pd.read_csv or read_excel
            self._try_simple_read(**kwargs)

            # check if we got data, and abort if so
            if self._data is not None: 
                self._df = self._data
                self._is_regular = True
                return

        # if haven't got data, then we go to parsing...

        # set is_horizontal to True for hybrid 
        # tables that are decorated
        if self._kind == "hybrid" and aux.from_kwargs("decorator", False, kwargs):
            is_horizontal = True
            self._hybrid_decorated = True

        # setup Parser to get data
        self._Parser = Parsers.CsvParser() if self._filesuffix() == "csv" else Parsers.ExcelParser()

        self._Parser.read(self._src, **kwargs)
        self._Parser.labels(id_label = self._id_col)
        self._Parser._make_BigTable_range(is_horizontal = is_horizontal)
        self._data = self._Parser._bigtable_range

    def parse(self, **kwargs):
        """
        Parses the big table and extracts the individual assays.
        """
        
        if self._kind == "vertical":
            self._parse_vertical(**kwargs)
        elif self._kind == "horizontal": 
            self._parse_horizontal(**kwargs)
        elif self._kind == "hybrid":
            decorator = aux.from_kwargs("decorator", self._hybrid_decorated, kwargs, rm = True)
            self._parse_hybrid(decorator = decorator, **kwargs)
    
    def _parse_hybrid(self, **kwargs):
        """
        Extracts assay datasets for hybrid big tables,
        it gets first the id_col and then all ct_cols, 
        and assembles new dfs of them into a dict
        """
        # first check the kind of data we got because it 
        # can either be a pandas dataframe or a numpy ndarray
        # depending on whether or not the file is regular and/or decorated

        # it's a dataframe if it's a "regular" big table
        if isinstance(self._data, pd.DataFrame):

            # check if we got ct cols to extract
            if self._ct_col is None: 
                aw.HardWarning("BigTableReader:no_ct_cols", traceback = False)

            # we got a nice dataframe with columns in to extract directly
            data = self._extract_from_hybrid_dataframe(  data = self._data, to_extract = self._ct_col  )
            
            # since in this setting the data was not decorated, 
            # we store all datasets into the assays
            self._assays = data

        # it's a numpy ndarray if it's an "irregular" big table
        else: 
            # check if the file is supposed to be decorated
            decorator = aux.from_kwargs("decorator", False, kwargs)
            if decorator: 

                # in case it is decorated we extract 
                # data based on the decorators
                self._hybrid_bigtable_extract_by_decorator()
            
            # if we dont have decorated data, we must have received a list for 
            # ct_col values to use. So we can just transform the np.ndarray to a datafram
            # and use the same approach as for "regular" big tables...
            else: 
                
                # check if we got ct cols to extract
                if self._ct_col is None: 
                    aw.HardWarning("BigTableReader:no_ct_cols", traceback = False)

                # transform the numpy array from the 
                # hybrid bigtable to a pandas dataframe
                data = self._convert_hybrid_nparray_to_dataframe()

                # we got a nice dataframe with columns in to extract directly
                data = self._extract_from_hybrid_dataframe(  data = data, to_extract = self._ct_col  )
                
                # since in this setting the data was not decorated, 
                # we store all datasets into the assays
                self._assays = data

    def _convert_hybrid_nparray_to_dataframe(self):
        """
        Converts a numpy array of a hybrid bigtable 
        generated by one of the Parsers into a 
        pandas DataFrame. It returns the DataFrame
        """
        # first get the row in which the id_col header is located
        # this should be the second row because if it was the first, 
        # then we would have gotten a "regular" table anyway, so there
        # must actually be some decorators present, but the user decided
        # to ignore them. But to be save, we specifically search again. 
        starting_index = np.argwhere( self._data == self._id_col )

        # and get the row index. We assume there is only one hit for the 
        # id_col since if it were otherwise, the Parsers would have noticed
        starting_index = starting_index.reshape(starting_index.size)
        starting_index = starting_index[0]

        # and now we can assemble the data for the dataframe
        names = self._data[ starting_index, : ]
        data = self._data[ (starting_index + 1):, : ]
        data = pd.DataFrame( data, columns = names )

        return data


    def _generate_subset_dataframe(self, data, decorator):
        """
        Generates a pandas DataFrame of a subset of columns 
        from a decorated hybrid big table. 

        Note
        ----
        This is a downsized version of the `qpcr.Parsers` `find_assays` core.

        Parameters
        -------
        data : np.ndarray
            A numpy array to search in. The **first** row will be searched. 
        """
        
        # compile decorator to search for 
        decorator = Parsers.decorators[ decorator ]
        decorator = re.compile( decorator )

        # get first row of the data
        array = data[ 0, : ]

        # set up an index array for the columns that match
        indices = np.zeros(len(array))

        # iterate over all entries in the first row
        idx = 0 
        for entry in array:
        # try: 
            match = decorator.search(entry)
            if match is not None: 
                indices[idx] = 1
        # except: 
        #     continue
            idx += 1

        # get matching indices and reduce dimensionality
        indices = np.argwhere(indices == 1)
        indices = indices.reshape(indices.size)
            
        # now also add the id_col column to the set of relevant indices
        # for this we search in the second row, but with exact matching.
        array = data[ 1, : ]
        id_col = np.argwhere( array == self._id_col )
        id_col = id_col.reshape(id_col.size)

        # and merge the id_col to the found indices
        indices = np.concatenate( (id_col, indices) )
        

        # now get the datframe relevant subset and convert into a DataFrame
        # we get all rows except the first (since there are the decorators)
        # and only the columns with matching indices.
        names = data[ 1, indices ]
        df = data[ 2:, indices ]

        # convert to numeric
        # we use the _convert_to_numeric method from the Parsers to do that
        # we iterate over each column (except the first where the identifiers are stored)
        # and convert all entries to numeric...
        tmp_parser = Parsers.CsvParser()
        for i in range( 1, df.shape[1] ):
            df[ :, i ] = tmp_parser._convert_to_numeric( "None (from BigTableReader)", df[ :, i ] )

        # convert to dataframe
        df = pd.DataFrame( df, columns = names )
        return df

    def _extract_from_hybrid_dataframe(self, data, to_extract):
        """
        Extracts datasets from a hybrid dataframe 

        Parameters
        ----------
        data : pd.DataFrame
            The dataframe to extract from. This needs to include both the id_col and all ct_cols.
        to_extract : str or list
            The column names to of Ct columns extract (the id_col is referenced from self._id_col )

        Returns 
        -------
        dfs : dict
            A dictionary of all extracted assays (assay id as key, df as value).
        """

        # check which columns we should extract
        # and convert to list if we only have a single one
        # just so we can use a loop uniformly
        if not isinstance(to_extract, (list,tuple)):
            to_extract = [to_extract]


        # get the id_column
        id_col = data[ self._id_col ]
        # and convert to string (just to be sure)
        id_col = id_col.astype(str)

         # setup a dict for the assay dataframes
        dfs = {}

        # iterate over all assays to extract
        for col in to_extract:
            # get ct values                
            ct_col = data[ col ]
            
            # convert to float (i.e. introduce nan were necessary)
            # to that end we use the same approach as the Parsers .make_dataframes()
            try: 
                ct_col = ct_col.astype(float)
            except: 

                ct_col = np.array( ct_col, dtype=str )
                try: 
                    ct_col = np.genfromtxt(  ct_col  )
                except: 
                    faulties = np.argwhere(    [  Parsers.float_pattern.match(i) is None for i in ct_col ]  ) 
                    ct_col[ faulties ] = "nan"
                    ct_col = np.genfromtxt(  ct_col  )


            # and assemble new dataframe
            tmp = pd.DataFrame(
                                    {   # using qpcr default column headers
                                        raw_col_names[0] : id_col, 
                                        raw_col_names[1] : ct_col
                                    }
                                )
            # and save dataframe
            dfs[ col ] = tmp

        # and store to data
        return dfs

    def _parse_horizontal(self, **kwargs):
        """
        Extracts assay datasets for a horizontal big table
        by first transforming it to a vertical one and then using 
        the vertical parse
        """
        # transform data into vertical
        self._data = self._Parser._infer_BigTable_groups(**kwargs)
        
        # transform array into df 
        # and set _is_regular to True so _parse_vertical 
        # wont try to also convert to df
        self._make_vertical_range_df()
        self._is_regular = True 

        # and parse_vertical to get assays
        ct_col = raw_col_names[1]
        assay_col = default_dataset_header
        self._id_col = raw_col_names[0]
        self._parse_vertical(ct_col = ct_col, assay_col = assay_col, **kwargs)



    def _parse_vertical(self, **kwargs):
        """
        Extracts assay datasets for vertical big tables
        """
        
        self._ct_col = aux.from_kwargs("ct_col", None, kwargs, rm = True)
        self._assay_col = aux.from_kwargs("assay_col", None, kwargs, rm = True)

        # in case of vertical tables, check if we got ct_col and assay_col
        got_no_cols = (self._ct_col is None or self._assay_col is None)
        if self._kind == "vertical" and got_no_cols:
            aw.HardWarning("BigTableReader:no_cols", ct_col = self._ct_col, assay_col = self._assay_col)
            
        # convert to pandas dataframe in case 
        # the data had to be parsed
        if not self._is_regular:
            self._make_vertical_range_df()
        
        # and test if the ones we have are good
        if not self._test_cols_are_good():
            aw.HardWarning("BigTableReader:cols_no_good", ct_col = self._ct_col, assay_col = self._assay_col)


        df = self._data
        assay_col_header = self._assay_col
        ct_col_header = self._ct_col
        id_col_header = self._id_col

        # get columns to include 
        cols_to_use = [id_col_header, ct_col_header]

        # now read the separate assays and store in assays and normalisers
        if self._vertical_decorated(): 
            # cols_to_use.append("@qpcr")
            self._get_vertical_assays_decorated(
                                                    df, 
                                                    assay_col_header, 
                                                    ct_col_header, 
                                                    cols_to_use
                                                )
        else:
            self._get_vertical_assays_not_decorated(
                                                        df, 
                                                        assay_col_header, 
                                                        ct_col_header, 
                                                        cols_to_use,
                                                        self._assays
                                                    )


    def _get_vertical_assays_decorated(self, df, assay_col_header, ct_col_header, cols_to_use):
        """
        Gets assays based on decorators, storing them in _assays and _normalisers
        """

        # get assays 
        tmp = df.query("`@qpcr` == 'assay'")
        self._get_vertical_assays_not_decorated(
                                                        tmp, 
                                                        assay_col_header, 
                                                        ct_col_header, 
                                                        cols_to_use,
                                                        self._assays
                                                    )

        # get normalisers
        tmp = df.query("`@qpcr` == 'normaliser'")
        self._get_vertical_assays_not_decorated(
                                                        tmp, 
                                                        assay_col_header, 
                                                        ct_col_header, 
                                                        cols_to_use,
                                                        self._normalisers
                                                    )

    def _get_vertical_assays_not_decorated(self, df, assay_col_header, ct_col_header, cols_to_use, store_in):
        """
        Gets assays without a decorator, storing them all in self._assays...
        """

        # get default names for the id and ct columns
        to_defaults = { _from : _to for _from, _to in 

                            zip(
                                    [self._id_col, self._ct_col], 
                                    raw_col_names
                                )
                    }

        # iterate over all assays
        assays = set(df[assay_col_header])
        for assay in assays:
            # get the assay subset
            subset = df.query(f"`{assay_col_header}` == '{assay}'")
            subset = subset[cols_to_use]

            # make Cts numeric
            cts = subset[ct_col_header].to_numpy()
            cts = np.genfromtxt(  np.array(cts, dtype=str)  )
            subset[ct_col_header] = cts
            
            subset = subset.reset_index(drop = True)
            # rename to defaults
            subset = subset.rename(columns = to_defaults)

            # save assay
            store_in.update({ assay : subset })

    def _vertical_decorated(self):
        """
        Checks if a vertical bigtable is decorated
        """
        return "@qpcr" in self._data.columns

    def _make_vertical_range_df(self):
        """
        Converts the numpy array from the Parser 
        to a pandas dataframe in case of irregular vertical big tables.
        """
        data = self._data
        headers = data[0, :]
        data = data[ 1: ,: ]
        df = pd.DataFrame(data, columns = headers)
        self._data = df

    def _test_cols_are_good(self):
        """
        Tests if specified ct and assay cols are valid
        """
        cols = self._data.columns
        all_good = self._ct_col in cols and self._assay_col in cols
        return all_good

    def _try_simple_read(self, **kwargs):
        """
        Try default readings without parsing in case the file is a regular
        csv or excel file (only works in case of vertical big tables).
        """
        if self._filesuffix() == "csv":
            
            delimiter = ";" if self._is_csv2() else ","
            data = pd.read_csv(self._src, delimiter = delimiter)
        
        else:
            
            sheet_name = aux.from_kwargs("sheet_name", 0, kwargs, rm = True)
            data = pd.read_excel(self._src, sheet_name = sheet_name)
        
        # check if we got the data we looked for...
        got_regular_data = self._id_col in data.columns
        if got_regular_data:
            self._data = data
        
    def _DataReader(self, **kwargs):
        """
        The DataReader interacting method
        """
        # replicates = aux.from_kwargs("replicates", None, kwargs)
        # self.replicates(replicates)
        # names = aux.from_kwargs("names", None, kwargs, rm = True)
        # self.names(names)
        data = self.pipe(**kwargs)
        return data

    def _hybrid_bigtable_extract_by_decorator(self):
        """
        Extracts assays and normalisers based 
        on decorators from a decorated hybrid big table.
        """
        # find all assays first
        data = self._generate_subset_dataframe( 
                                                    data = self._data, 
                                                    decorator = "qpcr:assay" 
                                            )
        # now that we have a nice dataframe 
        # we can select all columns that are not the id_col as our ct_cols of interest
        # and use the same _extract_from_hybrid_dataframe as for the undecorated hybrid tables.
        ct_cols = [ i for i in data.columns if i != self._id_col ]
        data = self._extract_from_hybrid_dataframe(  data = data, to_extract = ct_cols  )

        # and save to the assays
        self._assays = data

        # now repeat the same for the normalisers
        data = self._generate_subset_dataframe( 
                                                    data = self._data, 
                                                    decorator = "qpcr:normaliser" 
                                            )
        ct_cols = [ i for i in data.columns if i != self._id_col ]
        data = self._extract_from_hybrid_dataframe(  data = data, to_extract = ct_cols  )

        # and save to the normalisers
        self._normalisers = data


if __name__ == "__main__":

    # multisheet_file = "/Users/NoahHK/Downloads/Corti IPSCs July 2019_decorated.xlsx"
    # decorated_excel = "./__parser_data/excel 3.9.19_decorated.xlsx"

    # reader = MultiReader()
    # reader.read(decorated_excel, sheet_name = 1)
    # reader.parse(decorator = True, ignore_empty = True, assay_pattern = "Rotor-Gene")
    # reader.make_Assays()
    # r = reader.get("assays")
    # print(r[0].get())
    # assert r is not None, "MultiReader failed somewhere..."

    # reader = MultiSheetReader()
    # reader.pipe(
    #             multisheet_file, 
    #             # decorator = True, 
    #             assay_pattern = "Rotor-Gene"
    #         )
    # print(reader.get("Actin"))

    # bigtable_horiztonal = "/Users/NoahHK/Downloads/Local_cohort_Adenoma_qPCR_rawdata_decorated.xlsx"
    # bigtable_vertical = "/Users/NoahHK/Downloads/qPCR all plates.xlsx"

    # reader = BigTableReader()


    # reader.read(bigtable_vertical, kind = "vertical", id_col = "Individual")
    # reader.parse(ct_col = "Ct", assay_col = "Gene")
    # reader.make_Assays()
    # r = reader.get("assays")
    # print(r[0].get(), r[0].id())
    # reader.clear()


    # assays, normalisers = reader._DataReader(
    #                                     filename = bigtable_horiztonal, 
    #                                     kind = "horizontal", 
    #                                     id_col = "tissue_number",
    #                                     replicates = (3,4), 
    #                                     names = ["Gapdh", "Sord1"]
    #                                 )
    # r = normalisers
    # print(r[0].get())

    # reader = qpcr.DataReader()
    # r = reader.read( "./Examples/Example Data/actin.xlsx", header = 0, replicates = None, id = "myActin")
    # # r = reader.make_Assay()
    # print(r.get(), r.id())


    # reader.read( "./Examples/Example Data/actin_nan.csv", replicates = 6, id_label = "Hii", id = "myActin", is_regular = True )
    # r = reader.make_Assay()
    # print(r.get(), r.id())

    hybrid_bigtable = "./Examples/Example Data/Big Table Files/hybrid_bigtable.xlsx"

    # sheet 0 is not decorated 
    # sheet 1 is decorated 

    reader = BigTableReader()
    reader.read(
                    hybrid_bigtable, 
                    kind = "hybrid", 
                    id_col = "group", 
                    # ct_col = ["TLR1", "TLR4", "GAPDH"],
                    decorator = True, 
                    sheet_name = 1
            )
    reader.parse()
    reader.make_Assays()
    r = reader.assays()
    print(r)
    print(r[0][0].get())

    print(r[0][0].get().dtypes)
    r = reader.normalisers()
    print(r)