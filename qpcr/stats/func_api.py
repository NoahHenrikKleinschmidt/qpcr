"""
This module defines stand-alone functions for performing statistical analysis of results. 
"""

import qpcr.stats as stats

def assaywise_ttests( results, groups : list = None, columns : list = None, **kwargs ):
    """
    Perform multiple pairwise t-tests comparing the different `groups` within each `assay` within the Results dataframe `separately`.
    Hence, this function will compare for instance `ctrl-HNRNPL` against `KO-HNRNPL` but not `ctrl-SRSF11`. 
    
    Parameters
    ----------
    obj : qpcr.Results or qpcr.Assay or list
        A Results or Assay object to use for the comparison (if none is already linked). 
        Also a list can be passed.

    groups : list
        The groups to pair-wise compare. If this is a simple ``list``
        then all listed groups will be compared pair-wise. If this is a ``list of lists (or tuples)``
        then all provided pairs will be compared. By default all possible group pairings are compared.

    columns : list
        The columns of the dataframe to use as input data. By default this will all non-setup columns.
        You can pass a list of any subset of non-setup-cols here. As a shortcut you can restrict to only 
        valid Delta-Delta-Ct columns (i.e. `{}_rel_{}` columns using the `kwarg` ``restrict_ddCt = True``). 

    Returns
    -------
    results : ComparisonsCollection
        A collection of ``PairwiseComparison`` objects for each assay in the `Results` object's dataframe.
    """
    results = stats.__default_Evaluator__.assaywise_ttests( obj = results, groups = groups, columns = columns, **kwargs )
    return results 

def groupwise_ttests( results, groups : list = None, columns : list = None, **kwargs ):
    """
    Perform multiple pairwise t-tests comparing the different `assays` within each `group separately`.
    Hence, this method will compare for instance `ctrl-HNRNPL` against `ctrl-SRSF11` but not `KO-HNRNPL`. 
    
    Parameters
    ----------
    obj : qpcr.Results or qpcr.Assay or list
        A Results or Assay object to use for the comparison (if none is already linked). 
        Also a list can be passed.

    groups : list
        The groups to include in the comparison. This can 
        This can be a list of any valid subset of the `group_names` 
        or numeric `group identifiers` of the `Results` object.
        
    columns : list
        The columns (assays) to pair-wise compare. Any subset of non-setup-cols can be passed here. As a shortcut one can restrict to only 
        valid Delta-Delta-Ct columns (i.e. `{}_rel_{}` columns using the `kwarg` ``restrict_ddCt = True``), naturally this will not work if `drop_rel` has been called before. 
        If a simple ``list`` is passed then all listed columns will be compared pair-wise. In case of a ``list of lists (or tuples)``
        then all provided pairs will be compared. By default all possible column pairings are compared.

    Returns
    -------
    results : ComparisonsCollection
        A collection of ``PairwiseComparison`` objects for each group in the `Results` object's dataframe.
    """
    results = stats.__default_Evaluator__.groupwise_ttests( obj = results, groups = groups, columns = columns, **kwargs )
    return results

def assaywise_anova( results, equal_var : bool = True, groups : list = None, columns : list = None, **kwargs ):
    """
    Compare the different `groups` within each `assay` within the object's dataframe separately.
    Hence, this method will test for variance within `HNRNPL` and within `SRSF11` separately using all groups. 
    
    Parameters
    ----------
    results : qpcr.Results or qpcr.Assay or list
        A Results or Assay object to use for the comparison (if none is already linked). 
        Also a list can be passed.

    equal_var : bool
        Assume equal variance among all compared groups. 
        If this is `True` then a `oneway ANOVA` will be performed.
        Otherwise a `Kruskal-Wallis H-test` is performed. 

    groups : list
        The groups to include. The group declaration can be either
        through their `group_names` or their numeric `group identifiers`. 
        By default all groups that are present are included.

    columns : list
        The columns of the dataframe to use as input data. By default this will all non-setup columns.
        You can pass a list of any subset of non-setup-cols here. As a shortcut you can restrict to only 
        valid Delta-Delta-Ct columns (i.e. `{}_rel_{}` columns using the `kwarg` ``restrict_ddCt = True``). 

    Returns
    -------
    results : MultipleComparisons
        A collection of ``AnovaComparison`` objects for each assay in the `Results` object's dataframe.
    """
    results = stats.__default_Evaluator__.assaywise_anova( results, equal_var, groups, columns, **kwargs )
    return results
    
def groupwise_anova( results, equal_var : bool = True, groups : list = None, columns : list = None, **kwargs ):
    """
    Compare the different `assays` within each `group` within the object's dataframe separately.
    Hence, this method will test for variance within `ctrl` and within `knockout` separately using all data columns. 
    
    Parameters
    ----------
    results : qpcr.Results or qpcr.Assay or list
        A Results or Assay object to use for the comparison (if none is already linked). 
        Also a list can be passed.

    equal_var : bool
        Assume equal variance among all compared groups. 
        If this is `True` then a `oneway ANOVA` will be performed.
        Otherwise a `Kruskal-Wallis H-test` is performed. 

    groups : list
        The groups to include. The group declaration can be either
        through their `group_names` or their numeric `group identifiers`. 
        By default all groups that are present are included.

    columns : list
        The columns of the dataframe to use as input data. By default this will all non-setup columns.
        You can pass a list of any subset of non-setup-cols here. As a shortcut you can restrict to only 
        valid Delta-Delta-Ct columns (i.e. `{}_rel_{}` columns using the `kwarg` ``restrict_ddCt = True``). 

    Returns
    -------
    results : MultipleComparisons
        A collection of ``AnovaComparison`` objects for each assay in the `Results` object's dataframe.
    """
    results = stats.__default_Evaluator__.groupwise_anova( results, equal_var, groups, columns, **kwargs )
    return results