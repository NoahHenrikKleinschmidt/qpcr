"""
This is the ``PairwiseComparison`` class, which handles data from multiple pairwise-t-tests, conducted by the ``Evaluator``.
"""

from itertools import permutations
import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import numpy as np
import pandas as pd
from statsmodels.stats import multitest

logger = aux.default_logger()

class PairwiseComparison(aux._ID):
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
    labels : list
        A list of the group labels that were compared. These need to be in the same order as the pvalues array.
    subset : list
        A list of all groups that were tested that are of interest. This can be any subset of the labels.
    """
    __slots__ = ["_pvalues", "_orig_pvalues", "labels", "_effect_size", "_p_are_adjusted", "subset_groups"]
    def __init__( self, pvalues : np.ndarray, effect_size : np.ndarray = None, id : str = None, labels : list = None, subset : list = None ):
        super().__init__()
        self._id = id
        self._pvalues = pvalues
        self._effect_size = effect_size

        self._orig_pvalues = pvalues.copy()
        self._asymmetric_pvalues = None
        self._p_are_adjusted = False 

        self.labels = self._set_labels(pvalues, labels)
        self.subset_groups = subset if subset is not None else labels

    def get( self, which : str = "pvalues" ):
        """
        Parameters
        -------
        which : str
            Either the `pvalues` or the `effect sizes` as panda DataFrame.
        """
        if which == "pvalues":
            return self.pvalues
        elif which == "effect_size":
            return self.effect_size
        else:
            raise ValueError( f"which must be either 'pvalues' or 'effect_size'. Got '{which}' instead" )

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
        Fills up an assymettric 2D array of p-values / effect sizes (if provided) to a symmetric 2D array around the diagonal.
        """

        self._asymmetric_pvalues = self._pvalues.copy()

        rows, cols = self._pvalues.shape
        if rows != cols: 
            raise IndexError( "The p-values array is not square." )
        
        perms = tuple( permutations( np.arange(rows), r = 2 ) ) 
        for i,j in perms:
            if self._pvalues[j,i] == self._pvalues[j,i]:
                self._pvalues[i,j] = self._pvalues[j,i]

        if self._effect_size is not None:
            for i,j in perms:
                if self._effect_size[j,i] == self._effect_size[j,i]:
                    self._effect_size[i,j] = self._effect_size[j,i]
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
        else:
            pvals = self._melt( "pval" )
            final = pvals
        if self._p_are_adjusted:
            pvals_adj = self._melt( "pval_adj" )
            final[ "pval_adj" ] = pvals_adj[ "pval_adj" ]
        if self._effect_size is not None:
            effects = self._melt( "effect_size" )
            final[ "effect_size" ] = effects[ "effect_size" ]
        return final 


    def to_df(self, which : str = None):
        """
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
        if which == "adjusted":
            data = self._pvalues
        elif which == "raw":
            data = self._orig_pvalues
        elif which == "effects":
            data = self._effect_size
        if which is None: 
            if self._pvalues is None:
                return None
            data = self._pvalues
        p = pd.DataFrame( data, columns = self.labels[0], index = self.labels[1] )
        return p

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
    

    def _melt( self, which : str ):
        """
        Melts a 2D dataframe into a three column dataframe.
        """
        if which == "pval":
            res = self.to_df( "raw" )
        elif which == "pval_adj":
            res = self.to_df( "adjusted" )
        elif which == "effect_size":
            res = self.to_df( "effects" )
        
        # now melt the data
        res = res.melt( ignore_index = False )
        res.rename( columns = { "variable": "a", "value": which }, inplace = True )
        res["b"] = res.index
        res.reset_index( inplace = True, drop = True )
        res = res[ ["a", "b", which] ]
        return res

    def _set_labels(self, pvalues, labels):
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
Effect Sizes:
{"-" * length}
{self.to_df("effects")}
{"-" * length}
        """.strip()
        return s
    
    def __repr__(self):
        s = f"""PairwiseComparison(id={self._id}, labels={self.labels}, subset={self.subset_groups})"""
        return s


class MultipleComparisons:
    """
    A collection of multiple PairWiseComparison objects.
    This is being returned by the Evaluator when calling a pair-wise comparison test.
    """
    __slots__ = ["comparisons", "ids"]
    def __init__( self, comparisons : dict ):
        self.comparisons = list(comparisons.values())
        self.ids = list(comparisons.keys())
    
    def get( self ):
        """
        Returns
        -------
        list
            A list of PairWiseComparison objects.
        """
        return self.comparisons

    def stack( self ):
        """
        Stacks the stored comparisons into a single stacked dataframe.

        Returns
        -------
        pd.DataFrame
            A stacked dataframe with the p-values and effect sizes of all comparisons.
        """
        stacked = [c.stack() for c in self.comparisons]
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

    def __getitem__( self, id ):
        if id in self.ids:
            idx = self.ids.index(id)
        elif isinstance( id, ( int, list, tuple ) ):
            idx = id
        else:
            raise ValueError( f"id must be one of the ids in the comparison (or a valid index between 0-{len(self)}). Got '{id}' instead" )
        return self.comparisons[idx]

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
        s = f"""MultipleComparisons(comparisons={self.ids})"""
        return s