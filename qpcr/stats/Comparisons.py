"""
This is the ``Comparisons`` module that stores the ``Comparison`` classes that are responsible for handling the results of ``qpcr.stats`` computations.
Multiple ``Comparison`` objects are stored together in ``MultipleComparisons`` objects from where they are easily accessible.
"""

from itertools import permutations
import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import numpy as np
import pandas as pd
from statsmodels.stats import multitest

logger = aux.default_logger()

class Comparison(aux._ID):
    """
    The base class for all single ``Comparison`` classes.

    Parameters
    -------
    id : str
        The comparison id (such as the name of the assay).
    pvalues : np.ndarray
        The p-values for the comparison. These may be corrected or not (this class can correct them).
    """
    __slots__ = ["_pvalues", "labels", "subset_groups"]
    def __init__(self, pvalues : np.ndarray, id : str = None, labels : list = None, subset : list = None ):
        super().__init__()
        self._pvalues = pvalues
        self.id(id)

        self.labels = self._set_labels(pvalues, labels)
        self.subset_groups = subset if subset is not None else labels


    def subset( self, subset ):
        """
        Set a subset of interest for the comparison.

        Parameters
        ----------
        subset : list
            A list of all groups that were tested that are of interest. This can be any subset of the labels.
        """
        self.subset_groups = self._set_labels(self._pvalues, subset)
        return self

    
    @property
    def pvalues(self):
        """
        Returns
        -------
        np.ndarray
            The p-values for the comparison (corrected if correction was performed, else the originally provided ones).
        """
        return self._pvalues

    @property
    def pvalues_subset(self):
        """
        Returns
        -------
        pd.DataFrame
            The p-values for the comparison (corrected if correction was performed, else the originally provided ones) With index and columns labeled by the compared groups.
            This will be cropped to the groups "of interest" specified.
        """
        if self._pvalues is None: 
            return None
        p = self.to_df()
        p = p.loc[ self.subset_groups[0], self.subset_groups[1] ]
        return p 

    def _to_df(self, data):
        """
        This is the core function for to_df()
        
        Converts a data array into a pandas DataFrame with labeled index and columns.
        This will retain the 2D grid arrangement of the data.

        Parameters
        -------
        which : str
            Either `"raw"` or `"adjusted"` (for the respective pvalues) or `"effects"` (for the effect sizes).
            If correction was performed, then `"adjusted"` is the default else `"raw"`.

        Returns
        -------
        pd.DataFrame
            The p-values for the comparison (corrected if correction was performed, else the originally provided ones) With index and columns labeled by the compared groups.
            This will include all groups (labels) present in the data. Use ``pvalues_tested`` to get a dataframe cropped to "groups of interest".
        """
        if data is None:
            return None
        p = pd.DataFrame( data, columns = self.labels[0], index = self.labels[1] )
        return p

    @staticmethod
    def _set_labels(pvalues, labels):
        """
        Set labels for columns and rows for the pvalues and effect size dataframes
        """
        if labels is not None:
            if isinstance( labels[0], (tuple, list, np.ndarray) ):
                labels = labels
            elif isinstance( labels[0], (int, str) ):
                labels = [labels, labels]
            else:
                raise TypeError( f"labels must be a list of ints, strings, or tuples. Got '{type(labels[0]).__name__}' instead" )
        else:
            labels = [ np.arange(pvalues.shape[0]), np.arange(pvalues.shape[1]) ]
        return labels

    def __collection_export__( self ):
        return self._to_df( self._pvalues )

    def __getitem__( self, row, col ):
        return self._pvalues[row, col]
    
    def __gt__( self, other ):
        if isinstance( other, self.__class__ ):
            return self._pvalues > other._pvalues
        else:
            return self._pvalues > other

    def __lt__( self, other ):
        if isinstance( other, self.__class__ ):
            return self._pvalues < other._pvalues
        else:
            return self._pvalues < other
    
    def __ge__( self, other ):
        if isinstance( other, self.__class__ ):
            return self._pvalues >= other._pvalues
        else:
            return self._pvalues >= other

    def __le__( self, other ):
        if isinstance( other, self.__class__ ):
            return self._pvalues <= other._pvalues
        else:
            return self._pvalues <= other
    
    def __eq__( self, other ):
        if isinstance( other, self.__class__ ):
            return self._pvalues == other._pvalues
        else:
            return self._pvalues == other


class MultiTestComparison(Comparison):
    """
    Stores the results of a multiple testing comparison.
    """
    __slots__ = ['_orig_pvalues', '_corrected_pvalues', '_p_are_adjusted', '_asymmetric_pvalues']
    def __init__(self, pvalues : np.ndarray, id : str = None, labels : list = None, subset : list = None ):
        super().__init__( pvalues = pvalues, id = id, labels = labels, subset = subset )
        self._orig_pvalues = pvalues.copy()
        self._p_are_adjusted = False
        self._corrected_pvalues = None
        self._asymmetric_pvalues = None

    def adjust_pvalues(self, make_symmetric : bool = False, **kwargs):
        """
        Adjusts the p-values for the comparison by benjamini-hochberg.

        Parameters
        ----------
        make_symmetric : bool
            Convert the assymmetric pvalues array symmetric along the diagonal.
            This will also convert the effect size array symmetric along the diagonal.

        **kwargs
            Any additional keyword arguments to be passed to ``multitest.fdrcorrection``.

        Returns
        -------
        pvalues_adjusted : np.ndarray
            The adjusted p-values.
        """
        if self._pvalues is None:
            return None
        pval_mask = np.isfinite( self._pvalues )
        adjusted = multitest.fdrcorrection( self._pvalues[pval_mask], **kwargs )[1]
        logger.debug( f"adjusted values are:\n{adjusted}" )

        self._pvalues[pval_mask] = adjusted
        self._p_are_adjusted = True

        if make_symmetric:
            self.make_symmetric()

        return self._pvalues

    def make_symmetric(self):
        """
        Fills up an assymettric 2D array of p-values to a symmetric 2D array around the diagonal.
        """

        self._asymmetric_pvalues = self._pvalues.copy()
        self._pvalues = self._make_symmetric( self._pvalues )
        return self 

    def stack(self):
        """
        Stacks the 2D data arrays of pvalues, adjusted pvalues and effect sizes 
        into multi-column column dataframe format. 
        Where `a` and `b` denote the two partners in the comparison, 
        and `pval` is the unadjusted p-value for the comparison, 
        `pval_adj` is the adjusted p-value for the comparison (if performed), and finally
        and `effect_size` is the effect size for the comparison (if stored).

        Returns
        -------
        pd.DataFrame
            The stacked data.
        """
        if self._pvalues is None:
            return None
        
        pvals = self._to_df( self._orig_pvalues )
        pvals.name = "pval"
        final = self._melt( pvals )
        if self._p_are_adjusted:
            pvals_adj = self._to_df( self._pvalues )
            pvals_adj.name = "pval_adj"
            pvals_adj = self._melt( pvals_adj )
            final[ "pval_adj" ] = pvals_adj[ "pval_adj" ]
        return final 

    @property
    def pvalues_adjusted(self):
        """
        Returns
        -------
        np.ndarray
            The p-values for the comparison, corrected for multiple comparisons.
        """
        if not self._p_are_adjusted:
            return None
        return self._pvalues
    
    @property
    def pvalues_raw(self):
        """
        Returns
        -------
        np.ndarray
            The p-values for the comparison, originally provided.
        """
        return self._orig_pvalues

    @staticmethod
    def _make_symmetric( data ):
        """
        The core function of make_symmetric
        """
        rows, cols = data.shape
        if rows != cols: 
            raise IndexError( "The p-values array is not square." )
        
        perms = tuple( permutations( np.arange(rows), r = 2 ) ) 
        for i,j in perms:
            if data[j,i] == data[j,i]:
                data[i,j] = data[j,i]
        return data

    def _melt( self, data ):
        """
        Melts a 2D dataframe into a three column dataframe.
        """
        # now melt the data
        res = data.melt( ignore_index = False )
        res.rename( columns = { "variable": "a", "value": data.name }, inplace = True )
        res["b"] = res.index
        res.reset_index( inplace = True, drop = True )
        res = res[ ["a", "b", data.name] ]
        return res

    def __collection_export__( self ):
        return self.stack()

class PairwiseComparison(MultiTestComparison):
    """
    Stores results from a pair-wise comparison.

    Parameters
    -------
    id : str
        The comparison id (such as the name of the assay).
    pvalues : np.ndarray
        The p-values for the comparison. These may be corrected or not (this class can correct them).
    effect_size : np.ndarray
        The actual effect sizes for each comparison. This needs to be of the same shape as the pvalues array.
    statistic : np.ndarray
        The t-statsistics from the comparison. This needs to be of the same shape as the pvalues array. 
    labels : list
        A list of the group labels that were compared. These need to be in the same order as the pvalues array.
    subset : list
        A list of all groups that were tested that are of interest. This can be any subset of the labels.
    """
    __slots__ = ["_orig_pvalues", "_effect_size", "_tstats", "_p_are_adjusted"]
    def __init__( self, pvalues : np.ndarray, effect_size : np.ndarray = None, statistic : np.ndarray = None, id : str = None, labels : list = None, subset : list = None ):
        super().__init__( pvalues = pvalues, id = id, labels = labels, subset = subset )

        self._tstats = statistic
        self._effect_size = effect_size

    def get( self, which : str = "pvalues" ):
        """
        Parameters
        -------
        which : str
            Either the `pvalues` or the `effect sizes` as panda DataFrame.
        """
        data = {
                    "pvalues": self.pvalues,
                    "effect_size": self.effect_size,
                    "tstats": self.tstats
                }
        if which in data:
            return data[ which ]
        else:
            raise ValueError( f"which must be either 'pvalues', 'tstats', 'effect_size'. Got '{which}' instead" )

    def make_symmetric(self):
        """
        Fills up an assymettric 2D array of p-values, and t-statistics (if provided), and effect sizess (if provided) to a symmetric 2D array around the diagonal.
        """

        super().make_symmetric()
        if self._effect_size is not None:
            self._effect_size = self._make_symmetric( self._effect_size )
        if self._tstats is not None:
            self._tstats = self._make_symmetric( self._tstats )
        return self 


    def stack(self):
        """
        Stacks the 2D data arrays of pvalues, adjusted pvalues and effect sizes 
        into multi-column column dataframe format. 
        Where `a` and `b` denote the two partners in the comparison, 
        and `pval` is the unadjusted p-value for the comparison, 
        `pval_adj` is the adjusted p-value for the comparison (if performed), and finally
        and `effect_size` is the effect size for the comparison (if stored).

        Returns
        -------
        pd.DataFrame
            The stacked data.
        """
        if self._pvalues is None:
            return None
        final = super().stack()
        if self._tstats is not None:
            tstats = self._to_df( self._tstats )
            tstats.name = "t_stat"
            tstats = self._melt( tstats )
            final[ "t_stat" ] = tstats[ "t_stat" ]
        if self._effect_size is not None:
            effects = self._to_df( self._effect_size )
            effects.name = "effect_size"
            effects = self._melt( effects )
            final[ "effect_size" ] = effects[ "effect_size" ]
        return final 

    


    def to_df(self, which : str = None):
        """
        Converts a data array into a pandas DataFrame with labeled index and columns.
        This will retain the 2D grid arrangement of the data. 
        To export a column-arranged dataframe, use the ``stack`` method.

        Parameters
        -------
        which : str
            Either `"raw"` or `"adjusted"` (for the respective pvalues) or `"tstats"` (for the t-statistics) or `"effects"` (for the effect sizes).
            If correction was performed, then `"adjusted"` is the default else `"raw"`.

        Returns
        -------
        pd.DataFrame
            The p-values for the comparison (corrected if correction was performed, else the originally provided ones) With index and columns labeled by the compared groups.
            This will include all groups (labels) present in the data. Use ``pvalues_tested`` to get a dataframe cropped to "groups of interest".
        """
        data = {
                    "adjusted" : self._pvalues,
                    "raw" : self._orig_pvalues,
                    "effects" : self._effect_size,
                    "tstats" : self._tstats,
                    None : self._pvalues
                }
        data = data[ which ]
        return self._to_df( data )

    def save( self, filename : str ):
        """
        Saves the comparison to a file.

        Parameters
        ----------
        filename : str
            The filename to save the comparison to.
        """
        df = self.stack()
        df.to_csv( filename, index = False )

    
    @property
    def effect_size(self):
        """
        Returns
        -------
        np.ndarray
            The effect sizes for the comparison.
        """
        return self._effect_size
    
    @property
    def effects_subset( self ):
        """
        Returns
        -------
        pd.DataFrame
            The effect sizes for the comparison. With index and columns labeled by the compared groups.
            This will be cropped to the groups "of interest" specified.
        """
        if self._effect_size is None:
            return None
        p = self.to_df( which = "effects" )
        p = p.loc[ self.subset_groups[0], self.subset_groups[1] ]
        return p
    
    @property
    def statistic(self):
        """
        Returns
        -------
        np.ndarray
            The t-statistics for the comparison.
        """
        return self._tstats
    
    @property
    def statistic_subset( self ):
        """
        Returns
        -------
        pd.DataFrame
            The t-statistics for the comparison. With index and columns labeled by the compared groups.
            This will be cropped to the groups "of interest" specified.
        """
        if self._tstats is None:
            return None
        p = self.to_df( which = "effects" )
        p = p.loc[ self.subset_groups[0], self.subset_groups[1] ]
        return p

    def __str__(self):
        length = max( [ len( str(  self.to_df( i )  ).split("\n")[0] ) for i in ("effects", None) ] )
        adjusted = " (adjusted)" if self._p_are_adjusted else ""
        s = f"""
{"-" * length}
Pairwise Comparison: {self._id}
{"-" * length}
Pvalues{adjusted}:
{"-" * length}
{self.to_df()}
{"-" * length}
t-statistics:
{"-" * length}
{self.to_df("tstats")}
{"-" * length}
Effect Sizes:
{"-" * length}
{self.to_df("effects")}
{"-" * length}
        """.strip()
        return s
    
    def __repr__(self):
        s = f"""PairwiseComparison(id={self._id}, labels={self.labels}, subset={self.subset_groups})"""
        return s


class ComparisonsCollection:
    """
    A collection of multiple PairWiseComparison objects.
    This is being returned by the Evaluator when calling a pair-wise comparison test.
    """
    __slots__ = ["_dict"]
    def __init__( self, comparisons : dict ):
        self._dict = comparisons

    def get( self ):
        """
        Returns
        -------
        list
            A list of the stored Comparison objects.
        """
        return self.comparisons

    def assemble_df( self ):
        """
        Assembles all stored comparisons into a single stacked dataframe.

        Returns
        -------
        pd.DataFrame
            A stacked dataframe with the p-values and effect sizes of all comparisons.
        """
        stacked = [c.__collection_export__() for c in self.comparisons]
        for obj,df in zip( self.comparisons, stacked ):
            df[ defaults.raw_col_names[0] ] = obj.id()
        stacked = pd.concat( stacked, axis = 0 )
        return stacked

    def save( self, directory : str ):
        """
        Saves the comparisons to files.

        Parameters
        ----------
        directory : str
            The directory to save the files to.
        """
        fname = "{directory}/{id}.csv"
        for i in self:
            i.save( filename = fname.format( directory = directory, id = i.id() ) )

    @property
    def comparisons(self):
        """
        Returns
        -------
        A list of all stored Comparison objects.
        """
        return list(self._dict.values())

    
    @property
    def ids(self):
        """
        Returns
        -------
        A list of the ids of all stored Comparison objects.
        """
        return list(self._dict.keys())

    def __getitem__( self, id ):
        if id in self.ids:
            idx = self.ids.index(id)
        elif isinstance( id, ( int, list, tuple ) ):
            idx = id
        else:
            raise ValueError( f"id must be one of the ids in the comparison (or a valid index between 0-{len(self)}). Got '{id}' instead" )
        return self.comparisons[idx]

    def __add__( self, other ):
        if not isinstance( other, self.__class__ ):
            raise TypeError( f"other must be a MultipleComparisons object. Got '{type(other).__name__}' instead" )
        return ComparisonsCollection( self._dict + other._dict )

    def __iter__( self ):
        return iter(self.comparisons)
    
    def __len__( self ):
        return len(self.comparisons)
    
    def __str__(self):
        s = f"""Stored Comparisons"""
        names = [str(c) for c in self.ids]
        length = max( [len(n) for n in names] + [len(s)] )
        s = f"""{'-' * length}\n{s}\n{'-' * length}\n"""
        s += "\n".join([str(c) for c in self.ids])
        s += f"\n{'-' * length}"
        return s
    
    def __repr__(self):
        s = f"""ComparisonsCollection(comparisons={self.ids})"""
        return s


