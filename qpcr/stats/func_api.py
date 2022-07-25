"""
This module defines stand-alone functions for performing multiple comparisons.
"""

import qpcr.stats as stats

def groupwise_ttests( results, groups : (list or dict) = None, columns : list = None, **kwargs ):
    """
    This function will perform multiple pairwise t-tests comparing the different `groups` within each `assay` within the Results dataframe separately`.
    Hence, this function will compare for instance `ctrl-HNRNPL` against `KO-HNRNPL` but not `ctrl-SRSF11`. 
        
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
    results : PairwiseComparisons
        A collection of group-wise comparisons for each assay within the results object.
    """
    results = stats.__default_Evaluator__.groupwise_ttests( obj = results, groups = groups, columns = columns, **kwargs )
    return results 

def assaywise_ttests( results, groups : list = None, columns : (list or dict) = None, **kwargs ):
    """
    This function will perform multiple pairwise t-tests comparing the different `assays` within each `group separately`.
        Hence, this method will compare for instance `ctrl-HNRNPL` against `ctrl-SRSF11` but not `KO-HNRNPL`. 
        
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
        results : PairwiseComparisons
            A collection of ``PairwiseComparison`` objects for each group in the `Results` object's dataframe.
    """
    results = stats.__default_Evaluator__.assaywise_ttests( obj = results, groups = groups, columns = columns, **kwargs )
    return results