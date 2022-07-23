"""
This module defines stand-alone functions for performing multiple comparisons.
"""

import qpcr.stats as stats

def groupwise_ttests( results, groups : (list or dict) = None, columns : list = None, **kwargs ):
    """
    This function will perform multiple pairwise t-tests.
        
    Parameters
    -------
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
        A collection of group-wise comparisons for each assay within the results object.
    """
    results = stats.__default_Evaluator__.groupwise_ttests( obj = results, groups = groups, columns = columns, **kwargs )
    return results 