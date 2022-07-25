"""
This is the ``Anova`` class which is able to perform a oneway ANOVA analysis on provided data.
Two modes are available: ``groupwise`` and ``assaywise``. In `assaywise` mode, the ANOVA compares the groups within each data column (assay).
In `groupwise` mode, the ANOVA compares the data columns within each group of the dataframe overall.

If equal variance cannot be assumed for the data to be compared, instead of an ANOVA, a Kruskal-Wallis Test is performed.
"""

import numpy as np
import scipy.stats as scistats

import qpcr._auxiliary as aux
import qpcr.main as main
import qpcr.stats.Comparisons as Comparisons

logger = aux.default_logger()

class Anova(aux._ID):
    """
    Performs statistical evaluations of Results using an ANOVA model.
    Two modes are supported: ``groupwise`` and ``assaywise``. 
    Using either an ANOVA or Kruskal-Wallis Test.

    In `assaywise` mode, the ANOVA compares the groups within each data column (assay).
    In `groupwise` mode, the ANOVA compares the data columns within each group of the dataframe overall.
    """
    def __init__(self, id : str = None):
        super().__init__()
        self.id(id)
        self._obj = None
        self.assaywise_results = None
        self.groupwise_results = None
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
    
    def assaywise_anova( self, obj : (main.Results or main.Assay) = None, equal_var : bool = True, groups : list = None, columns : list = None, **kwargs ):
        """
        Compare the different `groups` within each `assay` within the object's dataframe separately.
        Hence, this method will test for variance within `HNRNPL` and within `SRSF11` separately using all groups. 
        
        Parameters
        ----------
        obj : qpcr.Results or qpcr.Assay or list
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
            A collection of ``AnovaComparison`` objects for each assay in the object's dataframe.
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

        # now restrict to only groups to be included
        ref_col = "group"
        if groups is not None:
            if not isinstance( groups, (list, tuple) ):
                raise ValueError( "groups must be a list or tuple" )
            if isinstance( groups[0], str ):
                ref_col = "group_name"
            elif isinstance( groups[0], int ):
                ref_col = "group"
            else: 
                raise ValueError( "groups must be a list of group names or group identifiers" )

            df = df[ df[ ref_col ].isin( groups ) ]
        

        # now set the method to employ, either ANOVA or kruskal
        method = self._oneway_anova if equal_var else self._kruskal

        self.assaywise_results = {}
        
        subsets = [ (i, df[ [ref_col, i] ])  for i in columns ]
        logger.debug( subsets )
        for name,subset in subsets:
            
            data = [ i for _,i in subset.groupby( ref_col ) ]
            data = [ np.squeeze( i.drop( ref_col, axis = 1 ).to_numpy().T  )      for i in data ]
            logger.debug( f"{data=}" )
            result = method( data, **kwargs )

            logger.debug( f"{result.pvalue}" )
            logger.debug( f"{result.statistic}" )
            
            result = Comparisons.AnovaComparison( id = name, pvalue = result.pvalue, statistic = result.statistic ) # , labels = columns )
            logger.debug(result)
            self.assaywise_results[name] = result

        self.assaywise_results = Comparisons.ComparisonsCollection( self.assaywise_results )
        self._results = self.assaywise_results
        return self.assaywise_results    
    
    def groupwise_anova( self, obj : (main.Results or main.Assay) = None, equal_var : bool = True, groups : list = None, columns : list = None, **kwargs ):
        """
        Compare the different `assays` within each `group` within the object's dataframe separately.
        Hence, this method will test for variance within `ctrl` and within `knockout` separately using all data columns. 
        
        Parameters
        ----------
        obj : qpcr.Results or qpcr.Assay or list
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
            A collection of ``AnovaComparison`` objects for each assay in the  object's dataframe.
        """
        if isinstance( obj, list ):
            return [ self.groupwise_ttests( i, groups, columns, **kwargs ) for i in obj ]

        if obj is not None: 
            self.link(obj)

        results = self._obj
        df = results._df.copy()

        # check if we should restrict to only conventional ddCt cols or all non-setup cols
        if columns is None:
            columns = results.data_cols
        elif kwargs.pop( "restrict_ddCt", False ):
            columns = results.ddCt_cols

        # now restrict to only groups to be included
        ref_col = "group"
        if groups is not None:
            if not isinstance( groups, (list, tuple) ):
                raise ValueError( "groups must be a list or tuple" )
            if isinstance( groups[0], str ):
                ref_col = "group_name"
            elif isinstance( groups[0], int ):
                ref_col = "group"
            else: 
                raise ValueError( "groups must be a list of group names or group identifiers" )

            df = df[ df[ ref_col ].isin( groups ) ]
        
        # now restrict to only columns to be included
        df = df[ [ ref_col ] + columns ]

        # now set the method to employ, either ANOVA or kruskal
        method = self._oneway_anova if equal_var else self._kruskal

        self.groupwise_results = {}
        
        subsets = [ (i[0], i[1][ columns ])  for i in df.groupby( ref_col ) ]
        logger.debug( subsets )
        for name,subset in subsets:
            
            data = [ subset[i] for i in subset.columns ]
            result = method( data, **kwargs )

            logger.debug( f"{result.pvalue}" )
            logger.debug( f"{result.statistic}" )

            result = Comparisons.AnovaComparison( id = name, pvalue = result.pvalue, statistic = result.statistic ) # , labels = groups )
            logger.debug(result)
            self.groupwise_results[name] = result

        self.groupwise_results = Comparisons.ComparisonsCollection( self.groupwise_results )
        self._results = self.groupwise_results
        return self.groupwise_results    

    @staticmethod
    def _oneway_anova( data , axis : int = 0, **kwargs ):
        """
        The core ANOVA method
        """
        logger.debug( data )
        return scistats.f_oneway( *data, axis = axis )
    
    @staticmethod
    def _kruskal( data, nan_policy : str = "omit", **kwargs ):
        """
        The core Kruskal method for unequal variance testing.
        """
        logger.debug( data )
        return scistats.kruskal( *data, nan_policy = nan_policy )


__default_Anova__ = Anova()