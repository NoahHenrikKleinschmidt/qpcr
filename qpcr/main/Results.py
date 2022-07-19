"""
This is the `qpcr.Results` class whose function is to accumulate results from various
`qpcr.Assay` objects and summarize them.
"""

import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import qpcr._auxiliary.warnings as aw
import qpcr.main as main

import logging
import re
import pandas as pd
import numpy as np
from scipy.stats import sem, t
import os

class Results(aux._ID):
    """
    Handles a pandas dataframe for data and computed results from a `qpcr` class. 
    
    Note
    -----
    This is a central data collection that can inherit directly from `qpcr.Assay` objects and from 
    externally computed sources. Please, note that it will not perform extensive vetting on its data input, 
    so make sure to only provide proper data input when manually assembling your `qpcr.Results`!
    """
    def __init__(self, id:str = None):
        super().__init__()
        if id is not None:
            self.id( id )

        self._df = pd.DataFrame()
        self._stats_df = pd.DataFrame()
    
    def __str__(self):
        _length = len( str(self._df).split("\n")[0] ) 
        s = f"""
{"-" * _length}
{self._df}
{"-" * _length}
        """.strip()
        if self.id_was_set():
            s = f"{'-' * _length}\nResults: {self._id}\n{s}"
        return s 

    def __repr__( self ):
        return self.__str__()

    def __len__(self):
        return len(self._df)

    def get(self):
        """
        Returns 
        -------
        data : pd.DataFrame
            The Results dataframe
        """
        return self._df


    def add_Ct(self, assay : main.Assay ):
        """
        Adds a `"Ct"` column with Delta-Ct values from an `qpcr.Assay`.
        It will store these as a new column using the Assay's `id` as header.

        Parameters
        -------
        assay : qpcr.Assay
            An `qpcr.Assay` object from which to import.
        """
        Ct = assay[ defaults.raw_col_names[1] ]
        Ct.name = assay.id()
        self.add( Ct )

    def add_dCt(self, assay : main.Assay ):
        """
        Adds a `"dCt"` column with Delta-Ct values from an `qpcr.Assay`.
        It will store these as a new column using the Assay's `id` as header.

        Parameters
        -------
        assay : qpcr.Assay
            An `qpcr.Assay` object from which to import.
        """
        self.add( assay.dCt )

    def add_ddCt(self, assay : main.Assay ):
        """
        Adds all `"rel_{}"` columns with Delta-Delta-Ct values from an `qpcr.Assay`.
        It will store these as new columns using the Assay's `id` + the `_rel_{}` composite id.

        Parameters
        -------
        assay : qpcr.Assay
            An `qpcr.Assay` object from which to import.
        """
        self.add( assay.ddCt )

    def add(self, data:(pd.Series or pd.DataFrame), replace : bool = False):
        """
        Adds some new datacolumn.

        Note
        ----
        The `column` argument has to be named for this to work. However, there are 
        already implemented methods dedicated to adding specifically Delta-Ct, Delta-Delta-Ct or just
        Ct values to the Results. In order to add a generic column from a numpy array or some other iterable
        just use default item setting (e.g. `results["new column"] = [1,2,3,4]`).

        Parameters
        ----------
        data : pd.Series or pd.DataFrame
            A named pandas Series or DataFrame that can be joined into the already
            stored dataframe. Note, a DataFrame may contain multiple columns.
        replace : bool
            In case results from a computation with the same identifiers are already stored
            no new data can be stored under that id. Either the new data must be renamed or
            `replace = True` must be set to overwrite the presently stored data. 
        """
        if isinstance(data, pd.Series):
            if data.name in self._df.columns:
                if not replace:  
                    e = aw.ResultsError( "name_overlap", name = data.name )
                    logging.error( e )
                    return
            self._df[data.name] = data
            # else: 
            #     self._df = self._df.join(data)
        
        elif isinstance(data, pd.DataFrame): 
            
            new = data.columns.unique()
            current = self._df.columns.unique()
            to_add = new
            if not replace:
                to_add = set(new) ^ set(current)
                
                # this line preserves the original order which is lost by the set()
                to_add = [ i for i in new if i in to_add ]
                
                if len(to_add) != len(new):
                    logging.info( f"Excluding {tuple(new.intersection(current))} due to name overlap. Use replace=True to force replacement." )

            to_add = list( to_add )
            self._df[ to_add ] = data[ to_add ]

    def __setitem__( self, key, value ):
        self._df[ key ] = value 

    def __getitem__( self, key ):
        if key in self._df.columns:
            return self._df[ key ]
        if key in self._stats_df.columns:
            return self._stats_df[ key ]
        
    def __delitem__( self, key ):
        if key in self._df.columns:
            del self._df[ key ]
        if key in self._stats_df.columns:
            del self._stats_df[ key ]

    def merge( self, *Results, all_cols : bool = False ):
        """
        Merge any number of `qpcr.Results` objects into this one.
        The same can be achieved using the + operator. 

        Note
        -----
        This operation will merge the columns of the Results' dataframes!

        Parameters
        ----------
        *Results
            An arbitrary number of `qpcr.Results` objects.
        all_cols : bool
            Set to `True` to merge not only the Delta-Delta-Ct columns (_rel_ columns)
            but also any additional columns. 
        """
        new_df = self._df.copy()
        for result in Results:
            df = result.get()

            if not all_cols:
                df = df[ result.ddCt_cols ]
            # we merge the dataframes first without adding 
            # some new id suffix, only do so if this fails
            try:

                # check if we have an overlap of column names
                intersect = set(df.columns).intersection( set(new_df.columns) )
                if intersect != {}:
                    raise IndexError( f"Duplicate column names were found: {intersect}" )
                new_df = pd.merge(new_df, df, 
                                    right_index = True, left_index = True, 
                                )
            except Exception as e:
                new_df = pd.merge(new_df, df, 
                                right_index = True, left_index = True, 
                                suffixes = [f"_{self.id()}", f"_{result.id()}"]
                            )
                logging.critical( e )
        self._df = new_df
        return self

    def __add__(self, other):
        self.merge( other )
        return self

    def rename( self, cols : dict ):
        """
        Renames columns according to a dictionary as key -> value.
        This is the same as calling `Results.rename_cols`. 

        Parameters
        ----------
        cols : dict
            A dictionary specifying old column names (keys) and new colums names (values).
        """
        self.rename_cols( cols )
        return self 

    def rename_cols(self, cols:dict):
        """
        Renames columns according to a dictionary as key -> value.
        This is the same as calling `Results.rename`.

        Parameters
        ----------
        cols : dict
            A dictionary specifying old column names (keys) and new colums names (values).
        """
        self._df = self._df.rename(columns = cols)
              
    def drop( self, *cols):
        """
        Drops all specified columns from the dataframe.
        This is used for normaliser pre-processing. 
        This is the same as calling `Results.drop_cols`.

        Parameters
        ----------
        *cols
            Any column names (as `str`) to be dropped.
        """
        self.drop_cols( *cols )
        return self

    def drop_cols(self, *cols):
        """
        Drops all specified columns from the dataframe.
        This is used for normaliser pre-processing.
        This is the same as calling `Results.drop`.

        Parameters
        ----------
        *cols
            Any column names (as `str`) to be dropped.
        """
        for c in cols: 
            del self[ c ] 

    def setup_cols( self, obj : (main.Assay or pd.DataFrame)):
        """
        Adopts the setup columns: `id, group, group_name` from another object.
       

        Parameters
        -------
        obj qpcr.Assay or qpcr.Results or pd.DataFrame
            Either a `qpcr.Assay` or a `qpcr.Results` or a pandas DataFrame 
            that has the given columns. 
        """
        self[ "id" ] = obj[ "id" ]
        self[ "group" ] = obj[ "group" ]
        self[ "group_name" ] = obj[ "group_name" ]

    def names(self, as_set = True):
        """
        Parameters
        ----------
        as_set : bool
            If `as_set = True` (default) it returns a set (as list without duplicates) 
            of assigned group names for replicate groups.
            If `as_set = False` it returns the full group_name column (including all repeated entries).
        
        Returns 
        -------
        names : list or None
            The adopted `group_names` 
            (only works if a `qpcr.Assay` has already been linked 
            using `adopt_names()`!)
        """
        if as_set: names = list( self._df["group_name"].unique() ) 
        return names

    def groups(self, as_set = True):
        """
        Parameters
        ----------
        as_set : bool
            If `as_set = True` (default) it returns a set (as list without duplicates) 
            of assigned group names for replicate groups.
            If `as_set = False` it returns the full group column (including all repeated entries).
        
        Returns
        -------
        groups : list
            The given numeric group identifiers of all replicate groups.
        """
        groups = list( self._df["group"].unique() ) if as_set else self._df["group"]
        return groups

    def drop_groups(self, groups : (list or str or int) ):
        """
        Removes specific groups of replicates from the DataFrame.

        Parameters
        ----------
        groups : list
            Either the numeric group identifiers or the group name, or an iterable thereof,
            of the groups to be removed, or a `regex` pattern defining which groups
            should be dropped (this is useful for systematically removing RT- groups etc.)
            A `regex pattern` can be supplied as well to match multiple group names.
        """
        # check for regex pattern
        # and get corresponding group names 
        if isinstance(groups, str):
            groups = [i for i in self._df["group_name"] if re.match(groups, i) is not None]
        elif isinstance(groups, int):
            groups = [groups]

        # get the right reference column and query to use to be 
        # used (either group or group_name)
        ref_query = "group != {group}" if isinstance( groups[0], int ) else "group_name != '{group}'"
        
        # remove groups from dataset
        for group in groups: 
            self._df = self._df.query(ref_query.format(group = group))

            # also drop from stats df
            if not len(self._stats_df) == 0:
                self._stats_df = self._stats_df.query(ref_query.format(group = group))
    
    def drop_rel(self):
        """
        Crops the `X_rel_Y` column-names of Delta-Delta-Ct results to just `X`.
        I.e. reduces back to the assay-of-interest name only.
        """
        to_change = {i : i.split("_rel_")[0] for i in self._df.columns if "_rel_" in i }
        self.rename_cols(to_change)

        # also recompute the stats df with new names...
        if not len(self._stats_df) == 0: 
            self.stats(recompute = True)

    def stats( self, recompute = False, iqr_limits : tuple = None, ci_level : float = 0.95 ):
        """
        Computes summary statistis about the replicate groups: 
        - `N (count)`
        - `Mean`
        - `Median`
        - `StDev` 
        - `IQR` 
        - `CI`

        of all replicate groups, for all datasets (assays).
        
        Parameters
        ----------
        recompute : bool
            Statistics will only be once unless recompute is set to `True`.
            The same dataframe can be directly accessed via this method once is has been computed.
        iqr_limits : tuple
            The lower and upper quantiles for the IQR computation. By default `(0.25, 0.75)`
        
        Returns
        -------
        stats_df : pd.DataFrame
            A new dataframe containing the computed statistics for each replicate group.
        """
        iqr_limits = (0.25, 0.75) if iqr_limits is None else iqr_limits
        _stats = {
                            "mean" : lambda x: np.nanmean( x, axis = 0), 
                            "stdev" : lambda x: np.nanstd( x, axis = 0), 
                            "median" : lambda x: np.nanmedian(x, axis = 0),
                            f"IQR_{iqr_limits}" : lambda x: np.nanquantile(x, iqr_limits[1], axis = 0) - np.nanquantile(x, iqr_limits[0], axis = 0),  
                            f"CI_{ci_level}" : lambda x: [ i for i in np.array( t.interval( ci_level, len(x)-1, loc = np.nanmean(x, axis = 0), scale = sem(x, nan_policy = "omit" ) ) ).transpose() ],
                        }
        
        # if stats_df is already present, return but sorted according to assays, not groups (nicer for user to inspect)
        if not len(self._stats_df) == 0 and not recompute:
            return self._stats_df

        self._stats_df = pd.DataFrame()

        for group, name in zip( self.groups(), self.names() ):
            subset = self._df.query( f"group == {group}" )
            _subset = subset.drop( columns = defaults.setup_cols, errors = "ignore" )
            logging.debug( _subset )
            
            # setup a stats dataframe with the right columns (in the right order)
            _stat = pd.DataFrame( columns = ["group", "group_name", "assay", "n"] + list(_stats.keys()) )
            
            # compute all statistics
            for label, func in _stats.items():
                s = func(_subset)
                logging.debug( f"{label}: {s}" )
                _stat[ label ] = s
            
            # fill in with groups and group names and assay identifiers in the right length
            _stat[ "group" ] = group
            _stat[ "group_name" ] = name
            _stat[ defaults.dataset_header] = _subset.columns
            _stat[ "n" ] = len(_subset)

            # and add to the stats dataframe
            self._stats_df = pd.concat( (self._stats_df, _stat), ignore_index = True ) 
            
        self._stats_df = self._stats_df.sort_values( defaults.dataset_header )
        return self._stats_df

    def save(self, path, df = True, stats = True):
        """
        Saves a csv file for each specified type of results.

        Parameters
        ----------
        path : str
            Path has to be a filepath if only one type of results shall be saved (i.e. either `df` or `stats`), 
            otherwise a path to the directory where both `df` and `stats` shall be saved.
        
        df : bool
            Save the results dataframe containing all replicate values (the full results).
            Default is `df = True`.
        
        stats : bool
            Save the results dataframe containing summary statistics for all replicate groups.
            Default is `stats = True`.
        
        """
        if df and stats and not os.path.isdir(path):
            e = aw.ResultsError( "save_need_dir" )
            logging.error( e )
            raise e 

        if df:
            # in case of raw results export we don't need the "assay" column as all 
            # assays are stored as separate columns anyaway, so it doesn't store any useful data
            _df = self._df
            if "assay" in _df.columns: _df = self._df.drop( columns = ["assay"] )
            self._save_single(path, _df, "_df")
        if stats:
            # compute stats if none have been computed yet...
            if len(self._stats_df) == 0:
                self.stats()
            self._save_single(path, self._stats_df, "_stats")


    def preview( self, kind : str = "AssayBars", mode : str = "static", **kwargs ):
        """
        A shortcut to call on a `qpcr.Plotters.PreviewResults` wrapper to visualise 
        the results.

        Parameters
        ----------
        kind : str
            The kind of Plotter to call. This can be any of the four wrapped 
            Plotters, e.g. `kind = "GroupBars"`.
        mode : str
            The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).

        Returns
        -------
        fig : plt.figure or plotly.figure
            The figure generated by `PreviewResults`.
        """
        import qpcr.Plotters as Plotters
        preview_results = Plotters.PreviewResults( mode = mode, kind = kind )
        preview_results.params( **kwargs )
        preview_results.link( self )
        fig = preview_results.plot()
        return fig 

    def __qplot__( self, **kwargs ):
        return self.preview

    @property
    def ddCt_cols( self ):
        """
        Returns
        -------
        cols
            A list of all {}_rel_{} columns within the Results's dataframe.
        """
        return [i for i in self._df.columns if "_rel_" in i]

    @property
    def is_empty(self):
        """
        Checks if any results have been stored so far.

        Returns
        -------
        bool
            `True` if NO data is yet stored, else `False`.
        """
        return len(self) == 0

    def _save_single(self, path, src, suffix=""):
        """
        Saves either self._df or self._stats_df to a csv file based on a path
        (path can be either filename or directory)
        """
        filename = path if not os.path.isdir(path) else os.path.join(path, f"rel_{self.id()}{suffix}.csv")
        src.to_csv(filename, index = False)
        