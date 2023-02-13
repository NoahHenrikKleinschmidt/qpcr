"""
This is the `qpcr.DataReader` that serves as a general Hub for the `qpcr.Readers`
and allows versatile file reading. 

Reading Data Files
==================

Setting up the DataReader is really easy and works just as setting up any other of the classes in ``qpcr``.

.. code-block:: python

    reader = qpcr.DataReader()

    # now we can read the file using
    assay = reader.read( some_file )


If the file contains multiple assays we likely want to specify

.. code-block:: python

    assays = reader.read( some_file, multi_assay = True )

In case there are both "assays-of-interest" and "normaliser" assays in our file, we will have to decorate the datafile (or manually sort the read assays in our script).
To learn more about file pre-processing so that qpcr can automatically read your setup, check out the `Decorator tutorial <https://qpcr.readthedocs.io/en/latest/tutorials/8_decorating_datafiles.html>`_ .

.. code-block:: python

    # if we have a decorated file, then we can simply do
    assays, normalisers = reader.read( some_file, multi_assay = True, decorator = True )


Reading multiple files
----------------------

The DataReader is able to read multiple files successively when passed a ``list``. 
However, the DataReader functions as a *wrapper* around the :ref:`qpcr.Readers <Readers>` and to save computations it sets up a suitable Reader
and then re-uses that same reader to read all successive files. However, if your files are differently formatted, you may supply ``reset = True``
to force the DataReader to set up a new Reader for each file it reads.

.. code-block:: python

    many_files = [ ... ]

    assays = reader.read( many_files, reset = True )

Using dedicated ``qpcr.Readers``
--------------------------------

Not all datafiles will be (easily) readable by the ``qpcr.DataReader``. This is not necessarily because the files are bad, simply because the DataReader makes use of mostly default Readers.
Hence, it may be that your file will not be readable by the DataReader but will be readable by a dedicated Reader such as a ``BigTableReader`` for instance. Check out the :ref:`qpcr.Readers <Readers>` for more details.
"""

import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import qpcr._auxiliary.warnings as aw
from qpcr import Readers


logger = aux.default_logger()


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
    just use one of the `qpcr.Readers` or even `qpcr.Parsers` directly.
    """

    __slots__ = ["_src", "_Reader", "_Data", "_tmp_data"]

    def __init__(self):
        super().__init__()
        self._Reader = None  # the functional core will be either a Reader
        self._Data = {}  # the _Data attribute will store any output from the Reader as a dictionary with filename : data structure.
        self._tmp_data = None  # this will not attempt to distinguish between assays / normalisers, or anything. It's just an archive of whatever data we got.
        # by default data is not stored by the DataReader, it's read only pipes through, but data can be stored using the .store() method.

    def Reader(self, Reader=None):
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

    def read_multi_assay(self, filename: str, decorator: (bool or str) = True, reset: bool = False, **kwargs):
        """
        Reads a single irregular *multi assay* datafile.

        Parameters
        ----------
        filename : str
            A filepath to an input datafile.

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

        Returns
        -------
        assays
            Two lists will be returned, one for assays and one for normalisers.
            In case of a non-decorated file, the second (normaliser) list is empty.
        """

        return self.read(filename=filename, multi_assay=True, decorator=decorator, reset=reset, **kwargs)

    def read_bigtable(self, filename: str, kind: str, decorator: (bool or str) = True, assay_col: str = None, id_col: str = None, ct_col: str = None, reset: bool = False, **kwargs):
        """
        Reads a single BigTable datafile.

        Parameters
        ----------
        filename : str
            A filepath to an input datafile.

        kind : str
            Specifies the kind of Big Table from the file.
            This may either be `"horizontal"`, `"vertical"`, or `"hybrid"`.

        decorator : str or bool
            Set if the file is decorated. This can be set either to `True` for `multi_assay` and multi-sheet (excel) or `big_table` files,
            or it can be set to a valid `qpcr decorator` for single assay files or single-sheet files.
            Check out the documentation of the `qpcr.Parsers` for more information on decorators.

        assay_col : str
            The column header specifying the assay identifiers.

        id_col : str
            The column header specifying the replicate identifiers
            (or "assays" in case of `horizontal` big tables).

        ct_col : str
            The column header specifying the Ct values.

        reset : bool
            If multiple input files should be read but they do not all
            adhere to the same filetype / datastructure, use `reset = True`
            to set up a new Reader each time `read` is called.

        **kwargs
            Any additional keyword arguments to be passed to the core Reader.
            Note, while this tries to be utmost versatile there is a limitation
            to costumizibility through the kwargs. If you require streamlined datareading
            use dedicated `qpcr.Readers` and/or `qpcr.Parsers` directly.

        Returns
        -------
        assays
            Two lists will be returned, one for assays and one for normalisers.
            In case of a non-decorated file, the second (normaliser) list is empty.
        """
        if assay_col is None:
            assay_col = defaults.dataset_header
        if id_col is None:
            id_col = defaults.id_header
        if ct_col is None:
            ct_col = defaults.ct_header
        return self.read(filename=filename, big_table=True, kind=kind, id_col=id_col, assay_col=assay_col, ct_col=ct_col, decorator=decorator, reset=reset, **kwargs)

    def read(self, filename: str, multi_assay: bool = False, big_table: bool = False, decorator: (bool or str) = None, reset=False, **kwargs):
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
            use dedicated `qpcr.Readers` and/or `qpcr.Parsers` directly.$

        Returns
        -------
        assays
            Either a single `qpcr.Assay` object or a list thereof.
            In case of a decorated file, two lists will be returned, one for assays and one for normalisers.
        """
        if isinstance(filename, list):
            return [self.read(i, multi_assay=multi_assay, big_table=big_table, decorator=decorator, reset=reset, **kwargs) for i in filename]

        self._src = filename
        # vet filesuffix
        suffix = self._filesuffix()
        if suffix not in defaults.supported_filetypes:
            e = aw.MultiReaderError("unknown_datafile", file=self._src)
            logger.critical(e)
            raise e

        if reset or self._Reader is None:
            self.reset()
            self._setup_Reader(multi_assay=multi_assay, big_table=big_table, decorator=decorator, **kwargs)

        # read file and return data
        data = self._Reader.__dreader__(filename=self._src, decorator=decorator, **kwargs)

        self._tmp_data = {self._src: data}
        return data

    def __str__(self):
        s = f"""
{self.__class__.__name__}:\t{self._id}
Current Reader:\t{type(self._Reader).__name__}
Self-stored data:\t{self._Data}
        """.strip()
        return s

    def __repr__(self):
        base = type(self._Reader).__name__
        data = self._Data
        return f"{self.__class__.__name__}({base=}, {data=})"

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

        use_multi = kwargs.pop("multi_assay", False)
        is_bigtable = kwargs.pop("big_table", False)

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


__default_DataReader__ = DataReader()
"""The default DataReader"""
