"""
This is the ``qpcr.Evaluator`` responsible for statistical evaluation of the Results from an analysis.
"""

import qpcr._auxiliary as aux
import qpcr.main as main
import qpcr.stats.PairwiseTests as PairwiseTests

logger = aux.default_logger()

class Evaluator(aux._ID):
    """
    Performs statistical evaluations of Results.
    """
    def __init__(self, id : str = None):
        super().__init__()
        self.id(id)
        self._results = None
    
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

                  
    def groupwise_ttests( self, obj : main.Results = None, groups : (list or dict) = None, columns : list = None, **kwargs ):
        """
        Perform multiple pairwise t-tests comparing the different `groups` within each `assay` within the Results dataframe separately`.
        Hence, this method will compare for instance `ctrl-HNRNPL` against `KO-HNRNPL` but not `ctrl-SRSF11`. 
        
        Note
        ----
        This will compute any combination `a,b` only once as the t-test of `b,a` yields the same. By default any skipped
        inverse combination is left blank. The blank fields can be filled with the corresponding values using the ``PairwiseComparison.make_symmetric`` method. 

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
            You can pass a list of any subset of non-setup-cols here. As a shortcut you can restrict to only 
            valid Delta-Delta-Ct columns (i.e. `{}_rel_{}` columns using the `kwarg` ``restrict_ddCt = True``). 

        Returns
        -------
        results : MultipleComparisons
            A collection of ``PairwiseComparison`` objects for each assay in the `Results` object's dataframe.
        """
        if isinstance( obj, list ):
            return [ self.groupwise_ttests( i, groups, columns, **kwargs ) for i in obj ]

        if obj is not None: 
            self.link(obj)

        results = PairwiseTests.__default_PairwiseTests__.groupwise_ttests( self._obj, groups, columns, **kwargs )
        self._results = results
        return results

    
    def assaywise_ttests( self, obj : main.Results = None, groups : (list) = None, columns : (list or dict) = None, **kwargs ):
        """
        Perform multiple pairwise t-tests comparing the different `assays` within each `group separately`.
        Hence, this method will compare for instance `ctrl-HNRNPL` against `ctrl-SRSF11` but not `KO-HNRNPL`. 
        
        Note
        ----
        This will compute any combination `a,b` only once as the t-test of `b,a` yields the same. By default any skipped
        inverse combination is left blank. The blank fields can be filled with the corresponding values using the ``PairwiseComparison.make_symmetric`` method. 

        Parameters
        ----------
        obj : qpcr.Results or list
            A Results object to use for the comparison (if none is already linked). 
            Also a list of qpcr.Results can be passed.

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
        results : MultipleComparisons
            A collection of ``PairwiseComparison`` objects for each group in the `Results` object's dataframe.
        """
        if isinstance( obj, list ):
            return [ self.assaywise_ttests( i, groups, columns, **kwargs ) for i in obj ]

        if obj is not None: 
            self.link(obj)

        results = PairwiseTests.__default_PairwiseTests__.assaywise_ttests( self._obj, groups, columns, **kwargs )
        self._results = results
        return results