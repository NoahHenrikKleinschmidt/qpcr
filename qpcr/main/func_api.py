"""
These are the stand-along functions to perform analyses using the ``qpcr`` classes. 
These functions make use of default instances of the qpcr classes and offer quicker workflow.
"""

import numpy as np
import qpcr.main.Assay as Assay
import qpcr.main.Analyser as Analyser
import qpcr.main.Normaliser as Normaliser
import qpcr.main.DataReader as DataReader
import qpcr.main.Calibrator as Calibrator


def read( filename : str, multi_assay : bool = False, big_table : bool = False, decorator : (bool or str) = None, reset = False, **kwargs):
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
            Set to `True` if you wish to reset the default DataReader before reading the file. This is required if
            already another file was read before which was of a different format than the current file.

        **kwargs
            Any additional keyword arguments to be passed to the core Reader.
            Note, while this tries to be utmost versatile there is a limitation
            to costumizibility through the kwargs. If you require streamlined datareading
            use dedicated `qpcr.Readers` and/or `qpcr.Parsers` directly.

        Returns
        -------
        assays
            Either a single `qpcr.Assay` object or a list thereof. 
            In case of a decorated file, two lists will be returned, one for assays and one for normalisers.
        """
        DataReader.__default_DataReader__.reset()
        return DataReader.__default_DataReader__.read( filename = filename, multi_assay = multi_assay, big_table = big_table, decorator = decorator, reset = reset, **kwargs )


def read_multi_assay( filename : str, decorator : (bool or str) = True, reset = False, **kwargs):
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
            Set to `True` if you wish to reset the default DataReader before reading the file. This is required if
            already another file was read before which was of a different format than the current file.

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
        DataReader.__default_DataReader__.reset()
        return DataReader.__default_DataReader__.read_multi_assay( filename = filename, decorator = decorator, reset = reset, **kwargs )

def read_bigtable( filename : str, kind : str, decorator : (bool or str) = True, assay_col : str = None, id_col : str = None, ct_col : str = None, reset : bool = False, **kwargs ):
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
            Set to `True` if you wish to reset the default DataReader before reading the file. This is required if
            already another file was read before which was of a different format than the current file.
            
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
        DataReader.__default_DataReader__.reset()
        return DataReader.__default_DataReader__.read_bigtable( filename = filename, kind = kind, id_col = id_col, assay_col = assay_col, ct_col = ct_col, decorator = decorator, reset = reset, **kwargs )



def delta_ct( assay : (Assay or list) ):
    """
    Computes DeltaCt values for a single or multiple `qpcr.Assay` objects.
    This is a synonym to the `analyse` function. 

    Note
    ----
    This is essentially a Wrapper for a default `Analyser` and it's `pipe` method.
    If you require non-default behaviour, please, set up your own `Analyser`. 
    
    Parameters
    -------
    assay : qpcr.Assay or list
        A single `qpcr.Assay` object or a list thereof.

    Returns
    -------
    qpcr.Assay or list
        The same as input but with computed DeltaCt values.
    """
    return Analyser.__default_Analyser__.pipe( assay )


def analyse( assay : (Assay or list) ):
    """
    Computes DeltaCt values for a single or multiple `qpcr.Assay` objects.
    This is a synonym to the `delta_ct` function. 

    Note
    ----
    This is essentially a Wrapper for a default `Analyser` and it's `pipe` method.
    If you require non-default behaviour, please, set up your own `Analyser`. 

    Parameters
    -------
    assay : qpcr.Assay or list
        A single `qpcr.Assay` object or a list thereof.

    Returns
    -------
    qpcr.Assay or list
        The same as input but with computed DeltaCt values.
    """
    return delta_ct( assay )


def normalise( assays : list, normalisers : list, mode : str = "pair-wise" ):
    """
    Computes DeltaDeltaCt values (Normalised fold changes) for a number of assays-of-interest and normalisers.

    Note
    ----
    This is essentially a Wrapper for a default `Normaliser` and it's `pipe` method.
    If you require non-default behaviour, please, set up your own `Normaliser`. 

    Parameters
    ----------
    assays : list
        A list `qpcr.Assay` objects that are "assays of interest".

    normalisers : list
        A list of `qpcr.Assay` objects that are "normalisers".

    mode : str
            The normalisation mode to use. This can be either `pair-wise` (default), 
            or `combinatoric`, or `permutative`.
            `pair-wise` will normalise replicates only by their partner (i.e. first against first, 
            second by second, etc.). `combinatoric` will normalise all possible combinations of a replicate 
            with all partner replicates of the same group from a normaliser (i.e. first against first, then second, then third, etc.).
            This will generate `n^2` normalised Delta-Delta-Ct values, where `n` is the number of replicates in a group.
            `permutative` will scramble the normaliser replicates randomly and then normalise pair-wise. This mode supports 
            a parameter `k` which specifies the times this process should be repeated, thus generating `n * k` normalised Delta-Delta-Ct values.
            Also, through setting `replace = True` replacement may be allowed during normaliser scrambling.
            Note, this setting will be ignored if a custom `norm_func` is provided.

    Returns
    -------
    results : Results
        A `qpcr.Results` object containing the computed results.
    """
    Normaliser.__default_Normaliser__.prune()
    return Normaliser.__default_Normaliser__.pipe( assays, normalisers )


def calibrate( assay : (Assay or list), dilution : (float or np.ndarray or tuple) = None, remove_calibrators : bool = True ):
    """
    Computes an efficiency from an `qpcr.Assay` object.
    
    This method will try to compute a new efficiency. To do this, it will check autonomously if
    `calibrator : {}` replicates are present and use these for computation. If none are 
    found it will assume the entire assay is to be used as calibrator.

    Note
    ----
    This will use a blank default Calibrator. If you have efficiencies already computed
    and wish to assign them, set up a `Calibrator` manually and load the data. 

    Parameters
    ----------
    assay : qpcr.Assay or list
        A `qpcr.Assay` object or a list thereof.

    dilution : float or np.ndarray
        The dilution step to be used. This must be a `float` fraction
        e.g. `0.5` for a `1 : 2` dilution series or `0.1` for a `1 : 10` series etc.
        If there are multiple steps because there is a gap in the dilution series. It is 
        necessary to supply a step for each group individually e.g. `[1,0.5,0.25,0.0625,0.03125]`.
        if there are 5 dilution steps (originally six but 0.125 was discarded). Note, both of the above also work with the inverse dilutions e.g. `2` or `[1,2,4,16,32]`.
        If the calibrators specify a dilution step already in their `group names` then the dilutions can be inferred automatically.
        More information about this can be found in the documentation of the `qpcr.Calibrator.dilution` method.

    remove_calibrators : bool
        If calibrators are present in the assay alongside other groups, 
        remove the calibrator replicates after efficiency calculation. 

    Returns
    -------
    assay
        The same as input but with updated efficiency. 
    """
    Calibrator.__default_Calibrator__.dilution( dilution )
    return Calibrator.__default_Calibrator__.calibrate( assay, remove_calibrators )
