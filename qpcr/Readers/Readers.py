"""
This module provides different `Reader` classes that allow reading simple and complex datafiles
of various architectures.

## Available Data Readers
---


"""


# Concept to link the Readers to the DataReader
# All Readers define a _DataReader method that 
# specifies which of its methods is supposed 
# to be used for (mostly pipe, sometimes read...)
# _DataReader methods *must* return the data they read!

from attr import asdict
import pandas as pd
import qpcr
import qpcr._auxiliary as aux
from qpcr._auxiliary import warnings as aw
import qpcr._auxiliary.defaults as defaults
import qpcr.Parsers as Parsers
import os
import numpy as np 
from copy import deepcopy 
from io import StringIO


__pdoc__ = {
    "_CORE_Reader" : True
}

# migrate default settings from __init__
raw_col_names = defaults.raw_col_names
supported_filetypes = defaults.supported_filetypes

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
        specify `name_label` and `Ct_label` as additional arguments.
        """
        suffix = self._filesuffix()
        
        # check for a valid input file
        if suffix not in supported_filetypes:
            aw.HardWarning("MultiReader:empty_data", file = self._src)

        if suffix == "csv":
            try: 
                self._csv_read()
            except:
                # setup parser
                parser = Parsers.CsvParser()
                self._prep_Parser(kwargs, parser)
                
                assay_of_interest = aux.from_kwargs("assay", None, kwargs, rm=True)
                
                # pipe the datafile through the parser
                parser.pipe(self._src, **kwargs)

                # get the data
                self._get_single_assay(parser, assay_of_interest)

        elif suffix == "xlsx":
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
        self._names = names

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
        if len(parser.assays()) > 1:
            if assay_of_interest is None: 
                aw.HardWarning("Reader:cannot_read_multifile", file = self._src, assays = parser.assays(), traceback = False)
            self._df = parser.get(assay_of_interest)
            self.id(assay_of_interest)
        else:
            assay_of_interest = parser.assays()[0]
            self._df = parser.get(assay_of_interest)
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
        try: 
            df = pd.read_csv(
                                self._src, 
                                sep = self._delimiter, 
                                header = self._has_header(), 
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
        suffix = self._src.split(".")[-1]
        return suffix

    def _make_new_Assay(self, name, df):
        """
        Makes a new Assay object and performs group() already...
        """
        new_assay = qpcr.Assay()
        new_assay.adopt(df)
        new_assay.id(name)
        
        new_assay.replicates(self._replicates)
        new_assay.group()

        if self._names is not None:
            new_assay.rename(self._names)
        return new_assay

class SingleReader(_CORE_Reader):
    """
    Reads qpcr raw data files in csv or excel format to get a single dataset. 

    Input Data Files
    ----------------
    Valid input files are either regular `csv` files, or  irregular `csv` or `excel` files, 
    that specify assays by one replicate identifier column and one Ct value column.

    Irregular input files may specify multiple assays as separate tables, 
    one assay has to be selected using the `assay` argument. 
    Separate assay tables may be either below one another (separated by blank lines!)
    or besides one another (requires `transpose = True`).

    Note
    ----
    This is the successor of the original `qpcr.Reader` (not the `qpcr.SampleReader`!).
    Hence, the `SingleReader` will return a pandas DataFrame of the dataset 
    directly using `get` but not an `qpcr.Assay`.

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
    def __init__(self, filename:str = None, **kwargs) -> pd.DataFrame: 
        super().__init__()
        self._src = filename
        self._delimiter = None
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
        
        Note, this assumes by default
        that the replicates are headed by the label `"Name"` and the corresponding Ct values
        are headed by the label `"Ct"`. Both labels have to be on the same row. 
        If these labels do not match your excel file, you may
        specify `name_label` and `Ct_label` as additional arguments.

        Parameters
        ----------
        filename : str
            A filepath to an input datafile.
        """
        self._src = filename
        if self._filesuffix() == "csv":
            self._delimiter = ";" if self._is_csv2() else ","
        super().read(**kwargs)

    def _DataReader(self, **kwargs):
        """
        The DataReader interacting method
        """
        self.read(**kwargs)
        data = self.get()
        data = self.make_Assay()
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
            self._parse_by_decorators(kwargs)
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
                data = assays, names
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

    def _parse_by_decorators(self, kwargs):
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
    | 28S   | ctrl | 7.65  | normaliser | 
    | 28S   | ctrl | 7.74  | normaliser | 
    | 28S   | ctrl | 7.54  | normaliser | 
    | 28S   | kd   | 7.86  | normaliser | 
    | 28S   | kd   | 7.57  | normaliser | 
    | 28S   | kd   | 7.67  | normaliser | 
    | Actin | ctrl | 11.67 | normaliser | 
    | Actin | ctrl | 11.54  | normaliser | 
    | ...   | ...  | ...   | ... | 


    #### `Horizontal` Big Tables
    Big Tables of this kind store replicates from assays in side-by-side columns.
    The replicates may be labelled numerically or all have the same column header. 
    A second column is required specifying the replicate identifier. 

    Note, this kind of setup *requires* decorators above the first replicate of each assay,
    as well as user-defined `replicates`!

    Example:

    |      | @qpcr:normaliser |      |      | @qpcr:normaliser |      |      |
    | ---- | ---------------- | ---- | ---- | ---------------- | ---- | ---- |
    | id   | 28S1  | 28S2  | 28S3  | Actin1 | Actin2 | ...  |
    | ctrl | 7.74 | 7.65 | 7.54 | 11.54 | 11.67 | ...  |
    | kd   | 7.86 | 7.57 | 7.67 | 11.43 | 11.56 | ...  |
    | ...  | ...  | ...  | ...  | …     | ...   | ...  |

    """
    def __init__(self):
        super().__init__()
        self._Parser = None
        self._kind = None  # horizontal or vertical 
        self._id_col = None
        self._ct_col = None
        self._assay_col = None

    
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
            This may either be `"horizontal"` or `"vertical"`.
        id_col : str
            The column header specifying the replicate identifiers 
            (or "group names" in case of `horizontal` big tables).
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
            self._try_simple_read(kwargs)

            # check if we got data, and abort if so
            if self._data is not None: 
                print(self._data)
                return

        # if haven't got data, then we go to parsing...

        # setup Parser to get data
        self._Parser = Parsers.CsvParser() if self._filesuffix() == "csv" else Parsers.ExcelParser()

        self._Parser.read(self._src)
        self._Parser.labels(id_label = self._id_col)
        print(self._Parser._data)
        self._Parser._make_BigTable_range(is_horizontal = is_horizontal)
        print(self._Parser._bigtable_range)

    def _try_simple_read(self, kwargs):
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
        


if __name__ == "__main__":

    multisheet_file = "/Users/NoahHK/Downloads/Corti IPSCs July 2019_decorated.xlsx"
    decorated_excel = "./__parser_data/excel 3.9.19_decorated.xlsx"

    # reader = MultiReader()
    # reader.read(multisheet_file, sheet_name = 1)
    # reader.parse(decorator = "qpcr:assay")
    # r = reader.get("assays")
    # print(r)
    # exit(1) 

    reader = MultiSheetReader()
    reader.pipe(
                multisheet_file, 
                # decorator = True, 
                assay_pattern = "Rotor-Gene"
            )
    print(reader.get("Actin"))

    # not decorated yet!
    bigtable_horiztonal = "/Users/NoahHK/Downloads/Local_cohort_Adenoma_qPCR_rawdata_decorated.xlsx"
    bigtable_vertical = "/Users/NoahHK/Downloads/qPCR all plates.xlsx"

    reader = BigTableReader()
    reader.read(bigtable_vertical, kind = "vertical", id_col = "Individual")
    reader.read(bigtable_horiztonal, kind = "horizontal", id_col = "tissue_number")
