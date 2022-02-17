"""
This module contains classes that help link the main qpcr classes 
(especially the Readers) to the streamlit interface of Qupid.

Note
----

Most of these methods are pretty hacky as the Readers and Parsers are very much designed to work with actual files and not the `streamlit.UploadedFile` data objects.
Hence, A lot of the data has to be extracted and then manually fed into the Readers and Parsers to then be usible for parsing and Assay extraction. 

The `QupidReader` is the central hub for allowing a Reader to work with streamlit uploaded data, 
however, it lacks a clear defined workflow like the Readers. Or uniform dataflow for different input datatypes.
Check out the actual `core.py` and `main.py` modules of the `Qupid` repository to learn more about actual usage.

### Main methods of the `QupidReader`: 

Note the "main" methods all return lists of `qpcr.Assay` objects
for assays and normalisers! However, they differ in the type of 
input they require, so make sure to check out their definitions
to learn more. 

-----------------------------------------------------------------------------------------
Method                          Use
-----------------------------------------------------------------------------------------
read_regular_csv_files          Reads a list of regular csv files.

read_irregular_csv              Reads a single irregular csv file.

parse_one_excel_sheet           Extracts assays from a single excel sheet.

parse_csv_file                  Extracts assays from an irregular csv file.

parse_excel_BigTable            Extracts assays from a bigtable (vertical or horizontal)

parse_csv_BigTable              Extracts assays from a bigtable (Vertical or horizontal)
-----------------------------------------------------------------------------------------
"""

import streamlit as st
import qpcr 
from io import StringIO, BytesIO
import pandas as pd 
from copy import deepcopy
import re


# we also define session here, as we would 
# otherwise have a circular import with controls
def session(key, value = None, reset = False):
    """
    Adds a variable to the st.session state or gets it
    """
    if value is not None :
        st.session_state[key] = value
    elif key in st.session_state     and value is None and not reset: 
        return st.session_state[key]
    elif key in st.session_state     and reset:
        st.session_state[key] = None
    elif key not in st.session_state:
        return None



class QupidReader:
    """
    This class will handle the reading of a streamlit `UploadedFile` object into something that can be read by 
    a `qpcr Reader or Parser`
    """
    def __init__(self):
        self._src = None        # the UploadedFile object
        self._name = None       # the name attribute of the object (filename)
        self._data = None       # the data that will be passed on to a Reader
        self._raw = None        # the raw data not yet ready for a Reader...

        self._is_read = False   # set to True after the _src has been read... 

        # parameters for an excel file
        self._sheetnames = None
    
    def get(self):
        """
        Returns readable data
        """
        return self._data
    
    def name(self):
        """
        Returns the filename
        """
        return self._name
    
    def setup(self, UploadedFile):
        """
        Sets up the QupidReader with an UploadedFile
        """
        # get the filename
        self._name = UploadedFile.name

        # store the object (deepcopy to be save...)
        self._src = deepcopy(  UploadedFile  )
    
    def filesuffix(self):
        """
        Returns the filesuffix
        """
        suffix = self._name.split('.')[-1]
        return suffix

    def is_csv(self):
        """
        Returns True if the file is a CSV file
        """
        return self.filesuffix() == "csv"
    
    def is_excel(self):
        """
        Returns True if the file is an Excel file
        """
        return self.filesuffix() == "xlsx"
    
    def read_excel(self):
        """
        preliminarliy reads an excel file
        """
        data = pd.read_excel(self._src, sheet_name = None)

        self._sheetnames = list(  data.keys()  )
        self._raw = data
        self._is_read = True

    def sheets(self):
        """
        Returns a list of sheet names of an excel file...
        """
        return self._sheetnames
    
    def is_multisheet(self):
        """
        Returns True if the file is a multi-sheet excel file
        """
        return len(self._sheetnames) > 1


    def read_irregular_csv(self):
        """
        Prepares the content of an irregular csv file for a MultiReader's CsvParser
        The returned data can be directly set to the Parser's `_data` attribute.

        Returns
        -----
        data : np.ndarray
            A numpy array that can directly be stored to a Parser
        """

        # prep file and adjust commas
        self._prep_csv(prepare_commas = True)
        
        delimiter = session("delimiter")

        contents = self.get()
        # now read the data and convert to numpy array
        df = pd.read_csv(contents, header = None, sep = delimiter)
        df = df.dropna(axis = 0, how = "all").reset_index(drop=True)
        data = df.to_numpy()

        self._data = data
        return data


    def read_regular_csv_files(self, session_source):
        """
        Reads a list of regular csv files from a session saved list
        and returns a list of qpcr.Assay objects.

        Parameters
        --------------------------------
        session_source : str
            The key of the session state variable to source from

        Returns
        -----
        assays : list
            A list of qpcr.Assay objects.
        """

        # set up a SingeReader        
        reader = qpcr.Readers.SingleReader()

        files = session(session_source)
        assays = []
        for file in files:
            
            # link the file to and return the prepped csv
            self.setup(file)
            self._prep_csv()
            prepped = self.get()

            # extract the assay name from the filename
            name = self.name()
            name = "".join(name.split(".")[0])
            # and adopt name as id
            # so that the assay later on knows what to use for the id
            reader._id = name

            # hackishly add to the reader
            reader._delimiter = session("delimiter")
            reader._header = 0
            reader._src = prepped
            reader._csv_read() 

            # pass on replicate information
            replicates_from_session(reader)

            # get new Assay and store to tmp list
            new = reader.make_Assay()
            assays.append(new)
        
        return assays



    def parse_one_excel_sheet(self, file, reader, col, sheet):
        """
        Parses one excel file sheet for data and returns two lists
        for assays, and normalisers of qpcr.Assay objects.
        Note this requires that a `reader` is passed that already has a set up `Parser`!

        Parameters
        --------
        file : st.UploadedFile
            The uploaded input file
        reader : qpcr.Readers.MultiReader
            A MultiReader with a set up `Parser`.
        col : int
            The column in which to search for assays.
        sheet : str
            The sheet name to read.

        Returns
        -----
        assays : list
            A list of qpcr.Assay objects.
        normalisers : list
            A list of qpcr.Assay objects.
        """
        # link the data to the Parser
        reader._Parser.read(file, sheet_name = sheet)
                    
        # now parse by decorators
        reader._parse_by_decorators(ignore_empty = True, col = col)
                    
        reader.make_Assays()
        assays, normalisers = reader.get("assays"), reader.get("normalisers")

        # clear Parser memory after each sheet...
        reader._Parser.clear()
        return assays,normalisers

    def parse_csv_file(self, data, reader, col):
        """
        Parses a csv file. Contrary to the parse_one_sheet here 
        both decorators have to be present since there will be no
        other sheets available.

        Parameters
        --------
        data : np.ndarray
            A numpy array that can be directly set to a Parser's `_data`
            attribute. 
        reader : qpcr.Readers.MultiReader
            A MultiReader with a set up `Parser`.
        col : int
            The column in which to search for assays.

        Returns
        -----
        assays : list
            A list of qpcr.Assay objects.
        normalisers : list
            A list of qpcr.Assay objects.
        """
        # add the data
        reader._Parser._data = data

        # parse by decorators
        reader._parse_by_decorators( col = col )
                    
        reader.make_Assays()
        assays, normalisers = reader.get("assays"), reader.get("normalisers")

        # clear Parser memory after each sheet...
        reader._Parser.clear()
        return assays, normalisers


    def parse_excel_BigTable(self, file, reader, sheet_name, is_horizontal):
        """
        Extracts datsets from a big table from an excel file.

        Parameters
        -----
        file : st.UploadedFile
            A file object than can be read directly by a Parser
        reader : qpcr.Readers.BigTableReader
            A BigTableReader with an already set up Parser
        sheet_name : str
            The sheet name to read.
        is_horizontal : bool
            True if the big table is a "horizontal" big table. 
        Returns
        -------------------
        assays : list
            List of qpcr.Assay objects
        normalisers : list 
            List of qpcr.Assay objects
        """
        # set the id_label to the assay_col in case of horizontal bigtables
        if is_horizontal:
            reader._Parser.labels( id_label = session("assay_col") )

        # link the data to the Parser
        reader._Parser.read(file, sheet_name = sheet_name)

                # generate the bigtable range
        reader._Parser._make_BigTable_range(is_horizontal = is_horizontal)
        reader._data = reader._Parser._bigtable_range

                # get the assays
        if is_horizontal:
            reader.parse( 
                                    replicates = reader._replicates, 
                                    names = reader._names,    
                                )
        else:
            reader._data = _make_vertical_range_df(reader._data)
                    # set data cols
            reader._assay_col = session("assay_col")
            reader._ct_col = session("ct_col")
            reader._id_col = session("id_col")
                    
            # parse vertical 
            _parse_vertical(reader)

        # convert to assays
        reader.make_Assays()
        assays, normalisers = reader.get("assays"), reader.get("normalisers")
        return assays,normalisers

    def parse_csv_BigTable(self, reader, is_horizontal):
        """
        Extracts datsets from a big table from a csv file.

        Parameters
        -----
        file : st.UploadedFile
            A file object than can be read directly by a Parser
        reader : qpcr.Readers.BigTableReader
            A BigTableReader with an already set up Parser
        is_horizontal : bool
            True if the big table is a "horizontal" big table. 
        Returns
        -------------------
        assays : list
            List of qpcr.Assay objects
        normalisers : list 
            List of qpcr.Assay objects
        """
        # prepare the csv file to be read 
        data = self.read_irregular_csv()


        # set the id_label to the assay_col in case of horizontal bigtables
        if is_horizontal:
            reader._Parser.labels( id_label = session("assay_col") )

                # link the data to the Parser
        reader._Parser._data = data

                # generate the bigtable range
        reader._Parser._make_BigTable_range(is_horizontal = is_horizontal)
        reader._data = reader._Parser._bigtable_range

                        # get the assays
        if is_horizontal:
            reader.parse( 
                                            replicates = reader._replicates, 
                                            names = reader._names,    
                                        )
        else:
            reader._data = _make_vertical_range_df(reader._data)
                            # set data cols
            reader._assay_col = session("assay_col")
            reader._ct_col = session("ct_col")
            reader._id_col = session("id_col")
                            
            # parse vertical 
            _parse_vertical(reader)

        # convert to assays
        reader.make_Assays()
        assays, normalisers = reader.get("assays"), reader.get("normalisers")
        return assays, normalisers

    def _prep_csv(self, prepare_commas = False):
        """
        Prepares a csv file to be read by a Reader
        """
        data = self._src.read().decode('utf-8')
        if prepare_commas:
            data = self._prepare_commas(data)
        data = StringIO(data)
        self._data = data


    def _prepare_commas(self, data):
        """
        Performs the `qpcr.Parsers.CsvParser._prepare_commas()` method to 
        make the commas equal within the entire csv file. We need this because
        otherwise pandas would read nonsense and the numpy array would also be crap.
        """
        
        content = data
        delimiter = session("delimiter")

        has_quotes = f'"{delimiter}"' in content

        # check if quotes are in datafile and adjust comma-patterns to use
        empty_comma_filler = f'{delimiter}""' if has_quotes else f"{delimiter}"
        comma_sep = f'"{delimiter}"' if has_quotes else f"{delimiter}"
        comma_sep = re.compile(comma_sep)

        # update commas to make equal...
        lines = content.split("\n")
        comma_counts = [len(comma_sep.findall(i)) for i in lines]
        max_commas = max(comma_counts)
        lines = [i + (max_commas - j) * empty_comma_filler for i, j in zip(lines, comma_counts)]
        
        new_content = "\n".join(lines)
        return new_content



def _make_vertical_range_df(data):
    """
    Replaces the BigTableReader _make_vertical_range_df
    ---
    Converts the numpy array from the Parser 
    to a pandas dataframe in case of irregular vertical big tables.

    Returns
    -----
    df : pd.DataFrame
        Of the numpy array
    """
    rows, cols = data.shape
    headers = data[0, slice(0, cols)]
    data = data[ slice(1, rows) , slice(0, cols) ]
    df = pd.DataFrame(data, columns = headers)
    return df


def _parse_vertical(reader):
        """
        A cropped version of the original _parse_vertical of the BigTableReader without the vetting
        """
        df = reader._data
        assay_col_header = reader._assay_col
        ct_col_header = reader._ct_col
        id_col_header = reader._id_col

        # get columns to include 
        cols_to_use = [id_col_header, ct_col_header]

        # now read the separate assays and store in assays and normalisers
        if reader._vertical_decorated(): 
            # cols_to_use.append("@qpcr")
            reader._get_vertical_assays_decorated(
                                                    df, 
                                                    assay_col_header, 
                                                    ct_col_header, 
                                                    cols_to_use
                                                )
        else:
            reader._get_vertical_assays_not_decorated(
                                                        df, 
                                                        assay_col_header, 
                                                        ct_col_header, 
                                                        cols_to_use,
                                                        reader._assays
                                                    )


def replicates_from_session(reader):
    """
    Sets up replicates and groupnames to a Reader from session variables
    """
    reader._replicates = session("replicates")
    reader._names = session("names")


def setup_parser_from_session(reader):
    """
    Sets up parameters of a Reader's Parser from session_state variables
    
    Parameters
    --------
    reader
        A qpcr.Readers object that has a ._Parser attribute.
    """
    # setup assay_pattern
    reader._Parser.assay_pattern(  session("assay_pattern")  )
            
    # set transposed
    if session("transpose"): 
        reader._Parser.transpose()

    # add the data column labels
    reader._Parser.labels(
                            id_label = session("id_col"), 
                            ct_label = session("ct_col")
                        )

def sheet_name_from_session(Qreader):
    """
    Gets the sheet_name to be read from a multi-sheet excel file
    where only a single sheet should be read using either the session
    specified name or the one found by the Qreader.
    """
    # get the sheet_name to read and parse
    # if no sheet_name was specified because there is only a single 
    # sheet anyway, then just use the first / only one found by Qreader
    sheet_name = session("sheet_name")
    sheet_name = Qreader.sheets()[0] if sheet_name is None else sheet_name
    return sheet_name


