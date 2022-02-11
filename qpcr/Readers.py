"""
This module provides different `Reader` classes that allow reading simple and complex datafiles
of various architectures.

## Available Data Readers
---


"""

import pandas as pd
import __init__ as qpcr
import qpcr._auxiliary as aux
from qpcr._auxiliary import warnings as aw
import qpcr.Parsers as Parsers
import os
import numpy as np 
from copy import deepcopy 
from io import StringIO


__pdoc__ = {
    "_CORE_Reader" : True
}

# migrate default settings from __init__
raw_col_names = qpcr.raw_col_names
supported_filetypes = qpcr.supported_filetypes



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



class MultiReader(qpcr.Assay, SingleReader, aux._ID):
    """
    Reads a single multi-assay datafile and reads assays-of-interest and normaliser-assays based on decorators.
    
    Input Data Files
    ----------------
    Valid input files are multi-assay irregular `csv` or `excel` files, 
    that specify assays by one replicate identifier column and one Ct value column.

    Separate assay tables may be either below one another (separated by blank lines!)
    or besides one another (requires `transpose = True`), but ALL in the SAME sheet!

    Assays of interest and normaliser assays *must* be marked using `decorators`.


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

    def get(self, which : str):
        """
        Returns the stored assays or normalisers.

        Parameters
        ----------
        which : str
            Can be either `"assays"` or `"normalisers"`.

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

        # setup assay_patterns if they were provided
        assay_pattern = aux.from_kwargs("assay_pattern", None, kwargs, rm = True)
        self._Parser.assay_pattern(assay_pattern)

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
        Extracts the decorated datasets (assays) from the read datafile.
        
        Parameters
        ----------
        **kwargs
            Any additional keyword arguments that should be passed to the `qpcr.Parsers`''s `parse` method that extracts the datasets.
        """
        # remove any decorator argument that the user may have tried to pass...  
        aux.from_kwargs("decorator", None, kwargs, rm = True)

        # check if file should be read transposed
        transpose = aux.from_kwargs("transpose", False, kwargs, rm = True)
        if transpose:
            self._Parser.transpose()

        # get data column labels
        id_label = aux.from_kwargs("id_label", "Name", kwargs, rm = True)
        ct_label = aux.from_kwargs("ct_label", "Ct", kwargs, rm = True)
        self._Parser.labels(id_label,ct_label)

        # setup assay_patterns if they were provided
        assay_pattern = aux.from_kwargs("assay_pattern", None, kwargs, rm = True)
        self._Parser.assay_pattern(assay_pattern)

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

        try: 
            self.read(filename, **kwargs)
        except: 
            self.read(filename)
            aw.SoftWarning("Parser:incompatible_read_kwargs", func = f"{type(self._Parser).__name__}'s read method")
        
        self.parse(**kwargs)
        self.make_Assays()

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
        self._all_assays = []
        self._all_normalisers = []

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
        """
        self._src = filename

        # read file to get all sheets
        sheets = pd.read_excel(self._src, sheet_name = None)
        sheets = sheets.keys()

        # now repetitively read all sheets and extract data
        for sheet in sheets:
            try: 
                super().read(self._src, sheet_name = sheet)
                self.parse(**kwargs)
                self.make_Assays()
                
                assays = self.get( which = "assays" )
                normalisers = self.get( which = "normalisers" )
                
                self._all_assays.extend(assays)
                self._all_normalisers.extend(normalisers)

            except: 
                # ERROR HERE
                print("sheet : ", sheet, " could not be read")

if __name__ == "__main__":

    multisheetreader = MultiSheetReader()

    multisheet_file = "/Users/NoahHK/Downloads/Corti IPSCs July 2019_decorated.xlsx"
    multisheetreader.read(
                            multisheet_file,
                            assay_pattern = "Rotor-Gene",
                            
                        )
    
    r = multisheetreader._all_normalisers
    print(r)