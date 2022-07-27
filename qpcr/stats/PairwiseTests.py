"""
This defines the TTests class which performs multiple pairwise t-tests on the groups or assays of a qpcr class.
Two modes are supported: ``groupwise`` and ``assaywise``. In `assaywise` mode, the t-tests are performed to compare the groups within each data column (assay) with each other.
In `groupwise` mode, the t-tests are performed to compare the data columns within each group of the dataframe overall.
P-values are corrected for multiple testing using the Benjamini-Hochberg procedure.


"""

import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import qpcr.main as main
import qpcr.stats.Comparisons as Comparisons
import qpcr.stats.StatsTest as StatsTest

from itertools import permutations
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

logger = aux.default_logger()


class PairwiseTests(StatsTest.StatsTest):
    """
    Performs statistical evaluations of Results using pairwise t-tests.
    Two modes are supported: ``groupwise`` and ``assaywise``.

    In `assaywise` mode, the t-tests are performed to compare the groups within each data column (assay) with each other.
    In `groupwise` mode, the t-tests are performed to compare the data columns within each group of the dataframe overall.
    """
    def __init__(self, id : str = None):
        super().__init__( id )
        self._effect_size_func = self._default_effect_size_func
    
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

                  
    def assaywise_ttests( self, obj : (main.Results or main.Assay) = None, groups : (list or dict) = None, columns : list = None, **kwargs ):
        """
        Perform multiple pairwise t-tests comparing the different `groups` within each `assay` within the Results dataframe separately`.
        Hence, this method will compare for instance `ctrl-HNRNPL` against `KO-HNRNPL` but not `ctrl-SRSF11`. 
        
        Note
        ----
        This will compute any combination `a,b` only once as the t-test of `b,a` yields the same. By default any skipped
        inverse combination is left blank. The blank fields can be filled with the corresponding values using the ``PairwiseComparison.make_symmetric`` method. 

        Parameters
        ----------
        obj : qpcr.Results or qpcr.Assay or list
            A Results or Assay object to use for the comparison (if none is already linked). 
            Also a list can be passed.

        groups : list or dict
            The groups to pair-wise compare. If this is a ``list``
            then all listed groups will be compared pair-wise. If this is a ``dict``
            then all key-value pairs will be compared. The group declaration can be either
            through their `group_names` or their numeric `group identifiers`. 
            By default all groups that are present will be compared pair-wise.

        columns : list
            The columns of the dataframe to use as input data. By default this will all non-setup columns.
            You can pass a list of any subset of non-setup-cols here. As a shortcut you can restrict to only 
            valid Delta-Delta-Ct columns (i.e. `{}_rel_{}` columns using the `kwarg` ``restrict_ddCt = True``). 

        Returns
        -------
        results : ComparisonsCollection
            A collection of ``PairwiseComparison`` objects for each assay in the `Results` object's dataframe.
        """
        if isinstance( obj, list ):
            return [ self.assaywise_ttests( i, groups, columns, **kwargs ) for i in obj ]

        if obj is not None: 
            self.link(obj)

        results = self._obj
        df = results._df.copy()

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
        ref_col = self._assaywise_get_ref_col(groups)

        # generate subsets for each data column to pair-wise evaluate the groups
        subsets = ( df[ [ref_col, i] ] for i in columns )
        
        # setup a dictionary for the overall results
        self.assaywise_results = {}
        
        for subset in subsets:
            
            # transpose the dataframe
            name = subset.columns[1]
            subset = self._squash_groups(subset, ref_col)
            
            # compute pairwise t-tests and effect size
            pvalues, tstats = self._pairwise_ttest( subset, groups, **kwargs )
            effect_sizes = self._pairwise_effect_size( subset, groups, **kwargs )
            
            logger.debug( f"{pvalues=}" )
            logger.debug( f"{effect_sizes=}" )
             
            # assemble results and store
            r = Comparisons.PairwiseComparison( id = name, pvalues = pvalues, effect_size = effect_sizes, statistic = tstats, labels = subset.columns, subset = labels )
            r.adjust_pvalues()
            self.assaywise_results[name] = r

        self.assaywise_results = Comparisons.ComparisonsCollection( self.assaywise_results )           
        self._results = self.assaywise_results
        return self.assaywise_results

    
    def groupwise_ttests( self, obj : (main.Results or main.Assay) = None, groups : (list) = None, columns : (list or dict) = None, **kwargs ):
        """
        Perform multiple pairwise t-tests comparing the different `assays` within each `group separately`.
        Hence, this method will compare for instance `ctrl-HNRNPL` against `ctrl-SRSF11` but not `KO-HNRNPL`. 
        
        Note
        ----
        This will compute any combination `a,b` only once as the t-test of `b,a` yields the same. By default any skipped
        inverse combination is left blank. The blank fields can be filled with the corresponding values using the ``PairwiseComparison.make_symmetric`` method. 

        Parameters
        ----------
        obj : qpcr.Results or qpcr.Assay or list
            A Results or Assay object to use for the comparison (if none is already linked). 
            Also a list can be passed.

        groups : list
            The groups to include in the comparison. This can 
            This can be a list of any valid subset of the `group_names` 
            or numeric `group identifiers` of the `Results` object.
            
        columns : list or dict
            The columns (assays) to pair-wise compare. By default this will all non-setup columns.
            You can pass a list of any subset of non-setup-cols here. As a shortcut you can restrict to only 
            valid Delta-Delta-Ct columns (i.e. `{}_rel_{}` columns using the `kwarg` ``restrict_ddCt = True``). 
            If a ``list`` is passed then all listed columns will be compared pair-wise. In case of a ``dict``
            then all key-value pairs will be compared. 

        Returns
        -------
        results : ComparisonsCollection
            A collection of ``PairwiseComparison`` objects for each group in the `Results` object's dataframe.
        """
        if isinstance( obj, list ):
            return [ self.groupwise_ttests( i, groups, columns, **kwargs ) for i in obj ]

        if obj is not None: 
            self.link(obj)

        results = self._obj
        df = results._df.copy()

        # check if we should restrict to only a subset of groups
        if groups is not None:
            if isinstance( groups[0], (int, np.int64) ):
                df = df[ df.group.isin( groups ) ]
                ref_col = "group"
            elif isinstance( groups[0], str ):
                df = df[ df.group_name.isin( groups ) ]
                ref_col = "group_name"
            else:
                raise ValueError( f"Invalid group type. Groups must be either integers or strings. Got: {type(groups[0])}" )
        else:
            ref_col = "group_name"

        # check if we should restrict to only conventional ddCt cols or all non-setup cols
        columns, comparisons, labels = self._prepare_pairwise_assays( columns, results, **kwargs )

        if len(columns) == 1:
            raise IndexError( "You must pass at least two columns to compare." )


        # now make the subsets over which to iterate
        subsets = df.groupby( ref_col )


        logger.debug( f"{df=}" )

        # setup the results dictionary
        self.groupwise_results = {}

        # now iterate over all groups
        for name,subset in subsets:

            # now drop the setup cols    
            subset = subset.drop( defaults.setup_cols, axis = 1 )
        
            # compute pairwise t-tests and effect size
            pvalues, tstats = self._pairwise_ttest(subset, comparisons, **kwargs)
            effect_sizes = self._pairwise_effect_size(subset, comparisons, **kwargs)
            
            # assemble results and store
            r = Comparisons.PairwiseComparison( id = name, pvalues = pvalues, effect_size = effect_sizes, statistic = tstats, labels = subset.columns, subset = labels )
            r.adjust_pvalues()
            self.groupwise_results[name] = r
        
        self.groupwise_results = Comparisons.ComparisonsCollection( self.groupwise_results )           
        self._results = self.groupwise_results
        return self.groupwise_results


    @staticmethod
    def _assaywise_get_ref_col(groups):
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
        for the pairwise ttest and effect size comparison, for the assaywise ttests.
        """

        # get pair-wise group permutations
        length = len(df.columns)
        cols = list(df.columns)

        # setup an empty array for the pvalues later
        out_array = np.full( (length, length), fill_value = np.nan ) 

        # and setup an index function to assign the pvalues to the right position
        index = lambda i, j: (cols.index(i), cols.index(j))
        return index, out_array

    @staticmethod
    def _prepare_pairwise_groups(groups, results):
        """
        Prepares the groups to be compared for pairwise comparison.
        """
        if groups is None:
            # labels = [ results.groups(), results.groups() ]
            # groups = list( permutations( results.groups(), r = 2 ) ) 
            labels = [ results.names(), results.names() ]
            groups = list( permutations( results.names(), r = 2 ) ) 
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
    
    @staticmethod
    def _prepare_pairwise_assays(assays, results, **kwargs):
        """
        Prepares the assays to be compared for pairwise comparison.
        """
        if assays is None:
            columns = results.data_cols
            labels = [ columns, columns ]
            combinations = list( permutations( columns, r = 2 ) ) 
        elif isinstance( assays, (list,tuple) ):
            columns = list(assays)
            labels = [ assays, assays ]
            combinations = list( permutations( assays, r = 2 ) ) 
        elif isinstance( assays, dict ):
            labels = [ i for i in assays.items() ]
            _labels = []
            for i in labels: _labels += [ j for j in i if j not in _labels ]
            labels = _labels
            labels = list( sorted( labels ) )
            columns = list(labels)
            labels = [ labels, labels ]
            combinations = list( assays.items() )
        elif kwargs.pop( "restrict_ddCt", False ):
            columns = results.ddCt_cols
            labels = [ columns, columns ]
            combinations = list( permutations( columns, r = 2 ) )
        return columns, combinations, labels
        

    def _default_effect_size_func( self, a, b, **kwargs ):
        return np.abs( np.nanmean(a) - np.nanmean(b) )

    def _pairwise_ttest( self, df, combinations, **kwargs ):
        """
        Performs a pair-wise t-test comparison between a given set of combinations.
        These combinations can be either "groups" or "assays".
        """

        # first prepare the pvalues and index function 
        # which will be different for groups or assays.
        index, pvalues = self._prepare_pairwise_vars( df )
        tstats = pvalues.copy()

        # now we can loop through the permutations
        for comb in combinations:
            j, i = index( *comb )

            # if we already have computed this permutation in reverse
            # we will skip this step (no need to compute it twice)
            if pvalues[ j,i ] == pvalues[ j,i ]:
                continue

            a, b = comb
            r = ttest_ind( df[a].dropna(), df[b].dropna(), **kwargs )
            pvalues[i,j] = r.pvalue
            tstats[i,j] = r.statistic

        logger.debug( pvalues )
        return pvalues, tstats 

    def _pairwise_effect_size(self, df, combinations, **kwargs):
        """
        Performs a pair-wise effect size between a given set of combinations.
        These combinations can be either "groups" or "assays".
        """
        
        index, effect_sizes = self._prepare_pairwise_vars( df )

        for comp in combinations:
            j, i = index( *comp )

            # if we already have computed this permutation in reverse
            # we will skip this step (no need to compute it twice)
            if effect_sizes[ j,i ] == effect_sizes[ j,i ]:
                continue
            
            a, b = comp
            effect_sizes[i,j] = self._effect_size_func( df[a], df[b], **kwargs )
        
        return effect_sizes

    def __str__(self):
        s = "Pairwise T-Tests"
        if self._results is not None:
            length = len( self._results.split("\n")[0] )
        else:
            length = len( s )
        s = f"""{'-' * length}\n{s}\n{'-' * length}\n{self._results}\n{'-' * length}"""

    def __repr__(self):
        return f"{self.__class__.__name__}(results={self._results})"


__default_PairwiseTests__ = PairwiseTests()
"""The default PairwiseTests"""