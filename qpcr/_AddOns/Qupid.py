"""
This module is designed to serve as the intermediary between the streamlit framework
for handling and storing data and the qpcr module. 

At its core it defines the QupidReader which handles file traffic from streamlit
to the qpcr.Readers and qpcr.Parsers. We require this intermediary because the Readers and Parsers
are only able to work with true files stored in  a system whereas the streamlit produced 
`UploadedFile` objects are stored in memory. 

To do it's job, QupidReader reads the UploadedFile objects and processes the data
into pandas DataFrames or numpy arrays and feeds these directly to the Readers and 
Parsers. To help with that, we also defined the ArrayParser that directly accepts a 
numpy array as its default data input and will not try to vet any data through reading
filenames or accessing the files using `open()` (like the Csv and Excel Parsers try to do). 

QupidReader defines a main reading method for each type of data input from each reader, 
and will set up the corresonding Reader accordingly.

The main reading methods of QupidReader are:
------------------------------------------------------------------------

SingleReader_read_regular           To read regular single-assay files
MultiReader_read                    To read multi-assay files
MultiSheetReader_read               To read multi-sheet multi-assay files
BigTableReader_read                 To read bigtable files
------------------------------------------------------------------------
->  Note, all main reading methods **return** lists of qpcr.Assay objects for both
    assays and normalisers, except for the SingleReader_read_regular which only returns
    a single qpcr.Assay object.
"""
import streamlit as st
import qpcr 
import qpcr._auxiliary as aux
from qpcr.Readers import * 
from qpcr.Parsers import ArrayParser

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
    Handles the Inferface between the streamlit app environment and the qpcr.Readers and qpcr.Parsers.
    It essentially takes the role of the qpcr.DataReader.

    It performs preliminary filereading to get the actual data from the UploadedFile Objects but then sets 
    defines methods for reading each specific type of input datafile through setting up the dedicated parser.
    However, contrary to the qpcr.DataReader there is no automatic inference of the correct Reader.

    The main reading methods all **return** qpcr.Assay objects. Either a list for assays and normalisers each, or just a single one in case of the SingleReader method.
    But they vary in their required inputs, naturally. 
    
    The QupidReader stores the actual raw data of a read file in a _data attribute, but non of the actual Assays!

    """
    def __init__(self):
        self._src = None
        self._filename = None
        self._data = None
        
        # parameters for an excel file
        self._sheetnames = None

    
    def get(self):
        """
        Returns 
        -------
        data
            The raw read data
        """
        return self._data
    
    def name(self):
        """
        Returns 
        -------
        name
            The filename
        """
        return self._name

    
    def sheets(self):
        """
        Returns 
        -------
        sheets
            A list of sheet names of an excel file.
        """
        if self._sheetnames is None:
            sheets = [0]
        else: 
            sheets = self._sheetnames
        return sheets

    def link(self, UploadedFile):
        """
        Sets up the QupidReader with an UploadedFile
        Parameters
        -----------
        UploadedFile : st.UploadedFile
            An UploadedFile object.
        """
        # get the filename
        self._name = UploadedFile.name

        # store the object (deepcopy to be save...)
        self._src = deepcopy(  UploadedFile  )

        # and reset the data
        self._data = None
        self._sheetnames = None

    def reset(self):
        """
        Resets the Reader
        """
        self.__init__()

    def filesuffix(self):
        """
        Returns
        -------
        suffix : str
            The filesuffix
        """
        suffix = self._name.split('.')[-1]
        return suffix

    def is_csv(self):
        """
        Returns
        -------
        bool
            True if the file is a CSV file
        """
        return self.filesuffix() == "csv"
    
    def is_excel(self):
        """
        Returns
        -------
        bool
            True if the file is an EXCEL file
        """
        return self.filesuffix() == "xlsx"


    def is_multisheet(self):
        """
        Returns 
        -------
        bool
            True if the file is a multi-sheet excel file
        """
        if self._sheetnames is not None: 
            return len(self._sheetnames) > 1
        return False


    def ncols(self):
        """
        Returns 
        -------
        cols : int 
            The number of columns within the stored raw data.
        """
        if isinstance(self._data, pd.DataFrame):
            cols = len( self._data.columns )
        else:
            cols = self._data.shape[1]
        return cols

    # ----------------------------------------------------------------
    #   Main Prelim reading methods
    # ----------------------------------------------------------------

    def read_excel(self, header = 0, to_numpy = False, drop_nan = True, **kwargs):
        """
        Preliminarliy reads an excel file.
        Parameters
        -------
        header : int or None
            Specifies the header row to use while reading (default is first row).

        to_numpy : bool
            If True the dataframe is converted to a numpy 
            array (if only a single sheet is present!)
            Otherwise, the dataframe will be stored directly.
            However, in either case, the dictionary that is normally
            returned by pd.read_excel is dropped.
        drop_nan : bool
            Will remove any all-nan lines from the data if True.
        """
        data = pd.read_excel(self._src, header = header, sheet_name = None)

        self._sheetnames = list(  data.keys()  )

        if drop_nan:
        # drop all-nan lines
            for sheet in self._sheetnames:
                data[ sheet ] = data[ sheet ].dropna(axis = 0, how = "all").reset_index(drop=True)

        # # and store data
        # self._data = data

        # check if we got only a single sheet anyway
        if len(self._sheetnames) == 1:
            # get the data from that only sheet
            data = data[ self.sheets()[0] ]
            if to_numpy:
                data = data.to_numpy()

        # or check if we are supposed to look only for a specific sheet
        elif self.is_multisheet() and not session("multi_sheet") and session( "sheet_name" ) is not None:
            data = data[  session( "sheet_name" )  ]
            if to_numpy:
                data = data.to_numpy()

        # and store new data
        self._data = data
        
    
    def read_csv(self, header = 0, to_numpy = False, drop_nan = True, **kwargs):
        """
        Preliminarliy reads a csv file.
        Parameters
        -------
        header : int or None
            Specifies the header row to use while reading (default is first row).
            
        to_numpy : bool
            If True the dataframe is converted to a numpy 
            array.
        drop_nan : bool
            Will remove any all-nan lines from the data if True.
        """
        # prep file and adjust commas
        self._prep_csv(prepare_commas = True)
        
        # get the appropriate delimiter from the session
        delimiter = session("delimiter")

        contents = self.get()

        # now read the data and drop all-nan lines
        df = pd.read_csv(contents, header = header, sep = delimiter)
        if drop_nan: 
            df = df.dropna(axis = 0, how = "all").reset_index(drop=True)
    
        if to_numpy:
            df = df.to_numpy()
        
        # and store data
        self._data = df



    # ---------------------------------------------------------------
    #   Main Reading Methods
    # ---------------------------------------------------------------

    # The nomenclature is always {CoreReader}_read[_{filetype}]

    def SingleReader_read_regular(self, file, **kwargs):
        """
        Reads a regular datafile in either excel or csv format.

        Parameters
        ----------
        file : st.UploadedFile
            An UploadedFile object to read.
        Returns
        -------
        assay : qpcr.Assay
            A qpcr.Assay object of the file's assay data.
        """

        # first setup file
        self.link(file)
        if self.is_csv():
            self.read_csv(**kwargs)
        else:
            self.read_excel(**kwargs)

        # setup SingleReader
        reader = SingleReader()

        # link setup data from session
        replicates_from_session(reader)

        # vet the dataframe
        self._data = reader._vet_single_assay_df(kwargs, self._data)

        # link the data
        reader._df = self._data
        reader._id = self.name().split(".")[0]

        # convert to assay
        assay = reader.make_Assay()

        return assay

    def MultiReader_read(self, file, **kwargs):
        """
        Reads an irregular decorated multi-assay datafile 
        in either excel or csv format.

        Parameters
        ----------
        file : st.UploadedFile
            An UploadedFile object to read.
        Returns
        -------
        assays : list
            A list of qpcr.Assay object of the file's assay data.
        normalisers : list
            A list of qpcr.Assay object of the file's normaliser-assay data.
        """

        # first setup file
        self.link(file)
        
        # if we got vertical arrangement (i.e. not transposed)
        # we drop all-nan columns because the decorator will 
        # inevitably pose a non-all-nan but still "blank" line 
        # within the data, so we can reduce some memory usage here.
        # however, in case of transposed datasets, we will need to have
        # the blank lines in there since the decorators are no longer valid
        # demarkations ... 

        drop_nan = not session("transpose")
        if self.is_csv():
            self.read_csv(to_numpy = True, drop_nan = drop_nan, **kwargs)
        else:
            self.read_excel(to_numpy = True, drop_nan = drop_nan, **kwargs)

        # setup a MultiReader
        reader = MultiReader()


        # link setup data from session
        replicates_from_session(reader)

        data = self._data #.astype(str)

        # setup an ArrayParser
        reader._Parser = ArrayParser()
        setup_parser_from_session(reader)
        reader._Parser.read( data )
        
        id_label, ct_label = session("id_col"), session("ct_col")

        reader.parse( 
                        decorator = True, 
                        id_label = id_label, 
                        ct_label = ct_label, 
                        # replicates = reader._replicates,
                        # names = reader._names,
                        **kwargs 
                )
        
        
        reader.make_Assays()

        assays, _ = reader.assays()
        normalisers, _ = reader.normalisers()

        return assays, normalisers

    def MultiSheetReader_read(self, file, **kwargs):
        """
        Reads an irregular decorated multi-assay and multi-sheet excel file.

        Parameters
        ----------
        file : st.UploadedFile
            An UploadedFile object to read.
        Returns
        -------
        assays : list
            A list of qpcr.Assay object of the file's assay data.
        normalisers : list
            A list of qpcr.Assay object of the file's normaliser-assay data.
        """

        # first setup file
        self.link(file)
        if session("transpose"):
            self.read_excel(drop_nan = False, **kwargs)
        else:
            self.read_excel(**kwargs)

        # setup MultiReader (we manually iterate over the sheets)
        reader = MultiReader()

        # link setup data from session
        replicates_from_session(reader)

        # setup an ArrayParser
        reader._Parser = ArrayParser()
        setup_parser_from_session(reader)


        # setup assays and normalsiers lists
        assays = []
        normalisers = []

        # iterate over all sheets
        for sheet in self.sheets():
            
            # get data sheet and 
            # convert to numpy array
            data = self._data[ sheet ]
            data = data.to_numpy()

            # feed data to Parser and parse for assays
            # since we now have multiple sheets, it's okey not to find
            # any assays or normalisers on a given sheet
            reader._Parser.read( data )
            id_label, ct_label = session("id_col"), session("ct_col")
            reader.parse( 
                            decorator = True, 
                            ignore_empty = True, 
                            id_label = id_label, 
                            ct_label = ct_label,
                            **kwargs 
                    )
            reader.make_Assays()

            # get assays and normalisers
            a, _ = reader.assays()
            n, _ = reader.normalisers()

            assays.extend( a )
            normalisers.extend ( n )

            # clear memory after each sheet
            reader.clear()
            reader._Parser.clear()

        return assays, normalisers

    def BigTableReader_read(self, file, **kwargs):
        """
        Reads a decorated big table file in excel or csv format.

        Parameters
        ----------
        file : st.UploadedFile
            An UploadedFile object to read.
        Returns
        -------
        assays : list
            A list of qpcr.Assay object of the file's assay data.
        normalisers : list
            A list of qpcr.Assay object of the file's normaliser-assay data.
        """

        self.link( file )
        if self.is_csv():
            self.read_csv()
        else:
            self.read_excel()

        # setup BigTableReader
        reader = BigTableReader()
        replicates_from_session(reader)

        # setup kind
        reader._kind = session( "kind" )

        # setup columns
        reader._id_col = session( "id_col" )
        reader._ct_col = session( "ct_col" )
        reader._assay_col = session( "assay_col" )

        # check if we have a regular vertical table, because in this case
        # we don't have to set up a parser and stuff. We check if we
        # got a pandas DataFrame from read_{} that already contains the Id col.

        is_regular = False
        if reader._kind == "vertical":
            is_regular = reader._id_col in self._data.columns
        
        # if we got a regular vertical big table, just use the df we got
        if is_regular: 

            reader._data = self._data
            reader._df = self._data
            reader._is_regular = True

        # if not, then we gotta start parsing over...
        else: 
            
            # since the file was not regular, 
            # we need to re-read again but without headers this time
            self.link( file )
            if self.is_csv():
                self.read_csv(header = None, to_numpy = True)
            else:
                self.read_excel(header = None, to_numpy = True)

            # setup Parser
            reader._Parser = ArrayParser()
            setup_parser_from_session(reader)

            # link the data
            data = self._data.astype(str)
            reader._Parser.read(  data  )

            is_horizontal = not reader._kind == "vertical"

            # set up "hybrid_horizontal" reading memory
            if reader._kind == "hybrid" and is_horizontal:
                reader._hybrid_decorated = True
            
            # switch id and assay cols in case of horizontal
            if reader._kind == "horizontal":
                reader._Parser.labels( id_label = reader._assay_col )

            # generate bigtable_range
            reader._Parser._make_BigTable_range( is_horizontal = is_horizontal )
            reader._data = reader._Parser._bigtable_range


        # update kwargs to specific parameters required by the different
        # kinds of big tables. This may be slightly redundant, from what 
        # was specified before but is actually required. 

        # for vertical big tables we need to add data columns 
        # again to kwargs manually...
        if reader._kind == "vertical":

            kwargs = dict(
                            id_col = reader._id_col,
                            ct_col = reader._ct_col, 
                            assay_col = reader._assay_col,
                            replicates = reader._replicates,
                            names = reader._names, 
                            **kwargs
                            
                        )

        # for horizontal big tables we need to 
        # specify replicates manually...
        elif reader._kind == "horizontal":
            
            kwargs = dict(
                            replicates = reader._replicates,
                            names = reader._names, 
                            **kwargs
                            
                        )
        
        # no additional setup required for hybrid big tables
        elif reader._kind ==  "hybrid":

            pass
       
        # extract the data through parsing
        reader.parse( decorator = True, **kwargs )
        reader.make_Assays()

        assays, _ = reader.assays()
        normalisers, _ = reader.normalisers()
        
        return assays, normalisers


    
        

    # ---------------------------------------------------------------
    #   Auxiliary methods
    # ---------------------------------------------------------------

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


def replicates_from_session(reader):
    """
    Sets up replicates and groupnames to a Reader from session variables
    Parameters
    -------
    reader 
        A qpcr.Readers objects
    """
    reader._replicates = session("replicates")
    reader._names = session("names")

def sheet_name_from_session(Qreader):
    """
    Gets the sheet_name to be read from a multi-sheet excel file
    where only a single sheet should be read using either the session
    specified name or the one found by the Qreader.

    Parameters
    -------
    Qreader : QupidReader
        A QupidReader
    Returns
    --------
    sheet_name : str
        The sheet_name to read,
    """
    # get the sheet_name to read and parse
    # if no sheet_name was specified because there is only a single 
    # sheet anyway, then just use the first / only one found by Qreader
    sheet_name = session("sheet_name")
    sheet_name = Qreader.sheets()[0] if sheet_name is None else sheet_name
    return sheet_name

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
