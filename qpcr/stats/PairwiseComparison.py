"""
This is the ``PairwiseComparison`` class, which handles data from multiple pairwise-t-tests, conducted by the ``Evaluator``.
"""

import qpcr._auxiliary as aux
import numpy as np
import pandas as pd
from statsmodels.stats import multitest

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

        self._orig_pvalues = pvalues
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
        
    def adjust_pvalues(self, **kwargs):
        """
        Adjusts the p-values for the comparison by benjamini-hochberg.

        Parameters
        ----------
        **kwargs
            Any additional keyword arguments to be passed to ``multitest.fdrcorrection``.

        Returns
        -------
        pvalues_adjusted : pd.DataFrame
            The adjusted p-values.
        """
        if self._pvalues is None:
            return None
        pval_mask = np.isfinite( self._pvalues )
        self._pvalues[pval_mask] = multitest.fdrcorrection( self._pvalues[pval_mask], **kwargs )[1]
        self._p_are_adjusted = True
        
        return self.pvalues

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
        pd.DataFrame
            The p-values for the comparison (corrected if correction was performed, else the originally provided ones) With index and columns labeled by the compared groups.
            This will include all groups (labels) present in the data. Use ``pvalues_tested`` to get a dataframe cropped to "groups of interest".
        """
        if self._pvalues is None: 
            return None
        p = pd.DataFrame( self._pvalues, columns = self.labels[0], index = self.labels[1] )
        return p

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
        p = self.pvalues
        p = p.loc[ self.subset_groups[0], self.subset_groups[1] ]
        return p 

    @property
    def effect_size(self):
        """
        Returns
        -------
        pd.DataFrame
            The effect sizes for the comparison. With index and columns labeled by the compared groups.
        """
        if self._effect_size is None: 
            return None
        p = pd.DataFrame( self._effect_size, columns = self.labels[0], index = self.labels[1] )
        p = p.loc[ self.subset_groups[0], self.subset_groups[1] ]
        return p
    
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
        p = self.effect_size
        p = p.loc[ self.subset_groups[0], self.subset_groups[1] ]
        return p
    
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
        length = len( str(self.pvalues).split("\n")[0] ) 
        adjusted = " (adjusted)" if self._p_are_adjusted else ""
        s = f"""
{"-" * length}
Pairwise Comparison: {self._id}
{"-" * length}
Pvalues{adjusted}:
{"-" * length}
{self.pvalues}
{"-" * length}
Effect Sizes:
{"-" * length}
{self.effect_size}
{"-" * length}
        """.strip()
        return s
    
    def __repr__(self):
        s = f"""PairwiseComparison(id={self._id}, labels={self.labels}, subset={self.subset_groups})"""
        return s


class PairwiseComparisons:
    """
    A collection of multiple PairWiseComparison objects.
    This is being returned by the Evaluator when calling a pair-wise comparison test.
    """
    __slots__ = ["comparisons", "ids"]
    def __init__( self, comparisons : dict ):
        self.comparisons = list(comparisons.values())
        self.ids = list(comparisons.keys())
    
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
        names = [str(c) for c in self.ids]
        length = max( [len(n) for n in names] )
        s = f"""{'-' * length}\nPairwise Comparisons\n{'-' * length}\n"""
        s += "\n".join([str(c) for c in self.ids])
        s += f"\n{'-' * length}"
        return s
    
    def __repr__(self):
        s = f"""PairwiseComparisons(comparisons={self.ids})"""
        return s