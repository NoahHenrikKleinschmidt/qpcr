"""
This is the ``qpcr.Evaluator`` responsible for statistical evaluation of the Results from an analysis.
"""

from qpcr import _auxiliary as aux
import qpcr.main as main
import qpcr.stats as qstats
import qpcr.stats.PairwiseComparison as PairwiseComparison

from itertools import permutations
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

logger = aux.default_logger()

class Evaluator(aux._ID):
    """
    Performs statistical evaluations of Results.
    """
    def __init__(self, id : str = None):
        super().__init__()
        self._id = id
        self._obj = None
        self.groupwise_results = None
        self._results = None
        self._effect_size_func = self._default_effect_size_func
    
    def link(self, obj : main.Results ):
        """
        Links a new object to evaluate.
        """
        self._obj = obj
    
    def get( self ):
        """
        Returns
        ------
        results 
            The results of the last performed test
        """
        return self._results

    def set_effect_size( self, f ):
        """
        Set the function to apply to compute the effect size.
        By default the absolute mean difference is used. 

        Parameters
        ----------
        f : function
            The function to apply to compute the effect size.
            The function must accept two arguments ``a`` and ``b``, which are numpy ndarrays.
            Alongside with any other keyword arguments. The function must return a single number.
        """
        self._effect_size_func = f
             
    def groupwise_ttests( self, obj : main.Results = None, groups : (list or dict) = None, columns : list = None, **kwargs ):
        """
        Perform multiple pairwise t-tests comparing the different `groups` within each `assay separately`.
        Hence, this method will compare for instance `ctrl-HNRNPL` against `KO-HNRNPL` but not `ctrl-SRSF11`. 
        Note, by default all non-setup cols are interpreted as data columns to perform tests on. You can restrict
        to only valid Delta-Delta-Ct columns (i.e. `{}_rel_{}` columns using ``restrict_ddCt = True``). 
        
        Parameters
        ----------
        obj : qpcr.Results or list
            A Results object to use for the comparison (if none is already linked). 
            Also a list of qpcr.Results can be passed.

        groups : list or dict
            The groups to pair-wise compare. If this is a ``list``
            then all listed groups will be compared pair-wise. If this is a ``dict``
            then all key-value pairs will be compared. The group declaration can be either
            through their `group_names` or their numeric `group identifiers`. 
            By default all groups that are present will be compared pair-wise.

        columns : list
            The columns of the Results dataframe to use as input data. By default this will all non-setup columns.
            You can pass a list of any subset of non-setup-cols here.

        Returns
        -------
        results : PairwiseComparisons
            A collection of ``PairwiseComparison`` objects for each assay in the `Results` dataframe.
        """

        if isinstance( obj, list ):
            return [ self.groupwise_ttests( i, groups, columns, **kwargs ) for i in obj ]

        if obj is not None: 
            self.link(obj)

        results = self._obj
        df = results._df

        # check if we should restrict to only conventional ddCt cols or all non-setup cols
        if columns is None:
            columns = results.data_cols
        elif kwargs.pop( "restrict_ddCt", False ):
            columns = results.ddCt_cols
            

        # get the groups to compare
        # the labels are for rows / columns annotations later
        # for the dataframes in PairwiseComparison
        groups, labels = self._prepare_pairwise_groups(groups, results)

        logger.debug( f"{labels=}" )
        # check the ref column to use
        ref_col = self._get_ref_col(groups)

        # generate subsets for each data column to pair-wise evaluate the groups
        subsets = ( df[ [ref_col, i] ] for i in columns )
        
        # setup a dictionary for the overall results
        self.groupwise_results = { i : None for i in columns }
        
        for subset in subsets:
            
            # transpose the dataframe
            name = subset.columns[1]
            subset = self._squash_groups(subset, ref_col)
            
            # compute pairwise t-tests and effect size
            pvalues = self._pairwise_ttest( df = subset, groups = groups, **kwargs )
            effect_sizes = self._pairwise_effect_size( df = subset, groups = groups, **kwargs )
            
            logger.debug( f"{pvalues=}" )
            logger.debug( f"{effect_sizes=}" )
             
            # assemble results and store
            r = PairwiseComparison.PairwiseComparison( id = name, pvalues = pvalues, effect_size = effect_sizes, labels = subset.columns, subset = labels )
            r.adjust_pvalues()
            self.groupwise_results[name] = r

        self.groupwise_results = PairwiseComparison.PairwiseComparisons( self.groupwise_results )           
        self._results = self.groupwise_results
        return self.groupwise_results

    @staticmethod
    def _get_ref_col(groups):
        """
        Checks if we have numeric or string groups and sets the reference column for subsetting accordingly.
        This happens AFTER the groups have been permuted into tuples.
        """
        logger.debug( f"{groups=}" )
        logger.debug( f"{type(groups[0][0])=}" )
        
        if isinstance(groups[0][0], (int, np.int64)):
            ref_col = "group"
        elif isinstance( groups[0][0], str ):
            ref_col = "group_name"
        else:
            raise ValueError( "groups must be a ints or strings" )
        return ref_col

    @staticmethod
    def _squash_groups(subset, ref_col):
        """
        Transposes the dataframe of a single ddCt_col to turn groups into columns.
        """
        rows = subset[ ref_col ].value_counts().max()
        rows = np.arange( rows )

        _prepped = pd.DataFrame( {"__blank" : rows } )
        for group, d in subset.groupby( ref_col ):
            d = d[ d.columns[-1] ]
            d.name = group
            d.reset_index( inplace = True, drop = True )
            _prepped = _prepped.join( d )
            logger.debug( d )
        
        del _prepped["__blank"]
        logger.debug( _prepped )
        return _prepped

    @staticmethod
    def _prepare_pairwise_vars(df):
        """
        Prepares an output array for outputs, and the index method
        for the pairwise ttest and effect size comparison.
        """
        # get pair-wise group permutations
        length = len(df.columns)
        cols = list(df.columns)

        # setup an empty array for the pvalues later
        out_array = np.full( (length, length), fill_value = np.nan ) 

        # and setup an index function to assign the pvalues to the right position
        index = lambda i, j: (cols.index(i), cols.index(j))
        return index, out_array

    def _prepare_pairwise_groups(self, groups, results):
        """
        Prepares the groups to be compared for pairwise comparison.
        """
        if groups is None:
            labels = [ results.groups(), results.groups() ]
            groups = list( permutations( results.groups(), r = 2 ) ) 
        elif isinstance( groups, (list,tuple) ):
            labels = [ groups, groups ]
            groups = list( permutations( groups, r = 2 ) ) 
        elif isinstance( groups, dict ):
            labels = [ i for i in groups.items() ]
            _labels = []
            for i in labels: _labels += [ j for j in i if j not in _labels ]
            labels = _labels
            labels = list( sorted( labels ) )
            labels = [ labels, labels ]
            groups = list( groups.items() )
        return groups,labels

    def _default_effect_size_func( self, a, b, **kwargs ):
        return np.abs( np.nanmean(a) - np.nanmean(b) )

    def _pairwise_effect_size(self, df, groups, **kwargs):
        """
        Computes the effect size of a pairwise comparison.
        """
        
        index, effect_sizes = self._prepare_pairwise_vars( df )

        for group in groups:
            j, i = index( *group )
            logger.debug( f"i={i}, j={j}" )

            # if we already have computed this permutation in reverse
            # we will skip this step (no need to compute it twice)
            if effect_sizes[ j,i ] == effect_sizes[ j,i ]:
                continue
            
            a, b = group
            effect_sizes[i,j] = self._effect_size_func( df[a], df[b], **kwargs )
        
        return effect_sizes

    def _pairwise_ttest(self, df, groups, **kwargs):
        """
        Performs a pair-wise t-test comparison between a given set of groups.
        """
        
        index, pvalues = self._prepare_pairwise_vars( df )

        # now we can loop through the permutations
        for group in groups:
            j, i = index( *group )
            logger.debug( f"i={i}, j={j}" )

            # if we already have computed this permutation in reverse
            # we will skip this step (no need to compute it twice)
            if pvalues[ j,i ] == pvalues[ j,i ]:
                continue

            a, b = group
            pvalues[i,j] = ttest_ind( df[a].dropna(), df[b].dropna(), **kwargs ).pvalue
            logger.debug( f"pvalues[{a},{b}]={pvalues[i,j]}" )
            
        logger.debug( pvalues )
        return pvalues 


