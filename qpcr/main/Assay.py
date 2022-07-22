"""
This is the ``qpcr.Assay`` class whose job is to store qPCR datasets. It is a central data-handling class in ``qpcr``. 

Setting up a ``qpcr.Assay``
===========================

Here is a manual example of creating a ``qpcr.Assay`` object. You can use either the ``qpcr.DataReader`` or any one of :ref:`qpcr.Readers <qpcr.Readers>` directly to 
read in your data and generate a pandas DataFrame. Note, the ``qpcr.Readers`` are already equipped with ``make_Assay(s)`` methods that will handle
setting up ``qpcr.Assay`` objects for you. 

However, setting up a ``qpcr.Assay`` manually can be as simple as:

.. code-block:: python

    # get the dataframe from one of the qpcr.Readers
    mydata = some_reader.get()

    assay = Assay( df = mydata, id = "my_assay" )

If your replicate identifiers are the same for all replicates within each group then the groups are automatically inferred. And your assay is 
ready at this point already to be passed to an `qpcr.Analyser`. If not, you can specify the replicates manually like this: 

.. code-block:: python

    # manually specify triplicates during setup
    assay = Assay( df = mydata, id = "my_assay", replicates = 3 )

    # or you can change the replicates after initial setup like 
    assay = Assay( df = mydata, id = "my_assay" )
    assay.group( replicates = 3 )

We can now actually interact with the `qpcr.Assay`. Assays support direct item setting, getting, and deleting on their dataframes.

.. code-block:: python

    # we could for instance fill a new column with only ones
    assay[ "my_new_column" ] = 1 

    # or get the id column from the assay
    ids = assay[ "id" ]

Specifying (Groups of) Replicates
=================================

The groups are essential to analysing our data, so ``qpcr`` needs to know about how the data is grouped. By the way, if you are unfamiliar with "groups" check out this 
Here's the best part: usually, we don't necessarily need to do anything here because ``qpcr.Assay`` are able to infer the groups of replicates in your data 
automatically from the replicate identifiers (yeah!). However, you will be asked to manually provide replicate settings in case this fails. 
In case you want to / have to manually specify replicate settings, a ``qpcr.Assay`` accepts an input ``replicates`` which is where you can specify this information. 

This input can be either an ``integer``, a ``tuple``, or a ``string``. Why's that? 
Well, normally we perform experiments as "triplicates", or "duplicates", or whatever multiplets.
Hence, if we always have the same number of replicates in each group (say all triplicates) we can simply specify this number as ``replicates = 3``. 
However, some samples might only be done in unicates (such as the diluent sample), while others are triplicates.

In these cases your dataset does not have uniformly sized groups of replicates and a single number will not do to describe the groups of replicates. 
For these cases you can specify the number of replicates in each group separately as a ``tuple`` such as ``replicates = (3,3,3,3,1)`` or as a ``string`` "formula"
which allows you to avoid repeating the same number of replicates many times like ``replicates = "3:4,1"``, which will translate into the same tuple as we specified manually. 

"""

import pandas as pd
import numpy as np
import logging
import os 
import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import qpcr._auxiliary.warnings as aw

from copy import deepcopy

logger = aux.default_logger()

raw_col_names = defaults.raw_col_names
class Assay(aux._ID):
    """
    The central storing unit of single datasets that were read from datafiles.
    An `qpcr.Assay` stores the replicate identifiers and Ct values, and also 
    groups these according to the `replicates` information (which is automatically
    inferred by default). Groups of replicates can be arbitrarily renamed by the user.

    Note
    -------
    The new implementation of the `qpcr.Assay` works directly with a DataFrame
    that was generated by any one of the `qpcr.Readers` or `qpcr.Parsers`.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame produces by one of the `qpcr.Readers` 
        containing an `id` column for the replicate identifiers 
        and a `Ct` value column. 
    id : str
        The identifer of the assays (the Assay name, essentially). 

    replicates : int or tuple or str
            Can be an `integer` (equal group sizes, e.g. `3` for triplicates), 
            or a `tuple` (uneven group sizes, e.g. `(3,2,3)` if the second group is only a duplicate). 
            Another method to achieve the same thing is to specify a "formula" as a string of how to create a replicate tuple.
            The allowed structure of such a formula is `n:m,` where `n` is the number of replicates in a group and `m` is the number of times
            this pattern is repeated (if no `:m` is specified `:1` is assumed). See `qpcr.Assay.replicates` for an example. 

    group_names : list
        A list of names to use for the replicates groups. If replicates of the same group share the same identifier, then the 
        group will be inferred automatically. Otherwise, default group names will be set if no `group_names` are provided. 
    """
    __slots__ = ["_df", "_id", "_efficiency", "_efficiency", "_replicates", "_group_names", "_groups"]

    def __init__(self, df : pd.DataFrame, id : str = None, replicates : (int or tuple or str) = None, group_names : list = None):
        super().__init__()
       
        if isinstance( df, pd.DataFrame ):
            self._df = df
        else:
            raise TypeError( f"df argument must be a pandas DataFrame (got {type(df).__name__})" )

        if id is not None: self._id = id

        # get replicates
        self._replicates = replicates

        # store names 
        self._names = group_names

        # setup the amplification efficiency
        self._efficiency = 1.0 
        self._eff = 2 * self._efficiency

        # now try to group the data    
        try: 
            self.replicates(self._replicates)
            self.group()
        except Exception as e:
            raise aw.AssayError( "setup_not_grouped" ) 
        
        # and try to change names, provided that we could group yet...
        if self._names is not None and self.groups() is not None: 
            self.rename(self._names)

    def __str__(self):

        _length = len( str(self._df).split("\n")[0] ) 
        s = f"""
{"-" * _length}
Assay: {self._id}
Amplif. Eff.: {self._efficiency}
{"-" * _length}
{self._df}
{"-" * _length}
        """.strip()
        return s 

    def __repr__( self ):
        id = self._id
        eff = self._efficiency
        n = len(self)
        return f"Assay({id=}, {eff=}, {n=})"

    def efficiency( self, eff : float = None ):
        """
        Gets or sets the amplification efficiency of the Assay.

        Parameters
        -------
        eff : float
            A new efficiency to assign to the assay.

        Returns
        -------
        float 
            The currently assigned efficiency.
        """
        logger.debug( f"{eff=}" )
        if isinstance( eff, float ):
            self._efficiency = eff 
            self._eff = 2 * self._efficiency
            logger.info( f"New efficiency set to {self._efficiency} (computes as binary factor {self._efficiency} * 2 = {self._eff})" )
        elif eff is None:
            return self._efficiency
        else:
            raise TypeError( f"Expected a float efficiency but got type {type(eff).__name__}" )

    def save(self, filename : str):
        """
        Saves the data from the `Assay` to a `csv` file.
        Parameters
        ----------
        filename : str
            The filename into which the assay should be stored.
            If this is a `directory`, then the assay `id` will automatically
            be used as filename. 
        """
        if os.path.isdir(filename):
            filename = os.path.join(filename, f"{self.id()}.csv")
        self.to_csv(filename, index = False)

    def get(self, copy : bool = False ):
        """
        Parameters
        -------
        copy : bool
            If `True` returns a deepcopy of the stored dataframe.

        Returns
        -------
        data : pandas.DataFrame
            The stored dataframe
        """
        if copy: 
            data = deepcopy( self._df )
        else: 
            data = self._df 
        return data

    def boxplot( self, mode : str = None, **kwargs ):
        """
        A shortcut to call a `qpcr.Plotters.ReplicateBoxPlot` plotter
        to visualise the loaded replicates.

        Parameters
        -------
        mode : str
            The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).
        **kwargs
            Any additional keyword arguments to be passed to the plotter.

        Returns
        -------
        fig : plt.figure or plotly.figure
            The figure generated by `ReplicateBoxPlot`.
        """
        import qpcr.Plotters as Plotters
        plotter = Plotters.ReplicateBoxPlot( mode = mode )
        plotter.link( self )
        fig = plotter.plot( **kwargs )
        return fig 

    def __qplot__( self, **kwargs):
        return self.boxplot

    def tile(self, n : int = 1):
        """
        Expands the dataframe to the square number of entries for each group.
        This is useful for combinatoric normalisation wherein each replicate is normalised
        against each replicate group-wise from the normaliser, instead of only its supposed partner value.
        
        Parameters
        -------
        n : int
            The number of tiles to produce. By default `1 tile` will effectively *square* the number of entries within the dataframe.
        """
        df = self._df
        groups = self.groups()

        new = None

        for group in groups: 
            subset = df.query(f"group == {group}")
            length = len(subset) * n
            subset = pd.concat( [subset for i in range(length) ], ignore_index = True )
            if new is None:
                new = subset
            else:
                new = pd.concat( [new, subset], ignore_index = True )

        self.adopt( new )

    def stack(self, n : int = 2):
        """
        Expands the dataframe entry-wise `n` times. 

        Parameters
        -------
        n : int
            The number of stacks to produce. `1 stack` will introduce one more copy of each replicate.
            Note, `n == 1` will keep the current entries!
        """
        df = self.get()
        groups = self.groups()

        n = int(n)

        new = None
        if n > 1:
            for group in groups: 
                subset = df.query(f"group == {group}")
                length = n
                subset = pd.concat( [subset for i in range(length) ], ignore_index = True )
                if new is None:
                    new = subset
                else:
                    new = pd.concat( [new, subset], ignore_index = True )

            self.adopt( new )

    @property
    def Ct(self):
        """
        Returns
        ------
        Ct : pandas.Series
            A pandas Series with the assay's Ct values. The column is renamed 
            from "Ct" to the assay's `id`.
        """
        Ct = self._df[ raw_col_names[1] ]
        Ct.name = self.id()
        return Ct

    @property         
    def dCt(self):
        """
        Returns
        -------
        dCt : pandas.Series
            A pandas Series with the computed Delta-Ct values. The column is renamed 
            from "dCt" to the assay's `id`.
        """
        dCt = self._df["dCt"]
        dCt.name = f"{self.id()}_dCt"
        return dCt

    @property
    def ddCt(self):
        """
        Returns
        -------
        ddCt : pandas.DataFrame
            A pandas DataFrame with all Delta-Delta-Ct values that the Assay has stored. 
            All `"rel_{}"` columns are renamed to include the assay `id` to `"{id}_rel_{}"`.
        """
        # get all ddCt columns
        ddCt = [ i for i in self._df.columns if "rel_" in i ]
        id = self._id
        # make new names and generate renaming dictionary 
        new_names = [ f"{id}_{i}" for i in ddCt ]
        new_names = {  old : new for new, old in zip(new_names, ddCt)  }

        # get the data and rename
        ddCt = self._df[ ddCt ]
        if not isinstance(ddCt, pd.DataFrame):
            ddCt = pd.DataFrame(ddCt)
        ddCt = ddCt.rename(columns = new_names)
        
        return ddCt

    @property
    def ddCt_cols( self ):
        """
        Returns
        -------
        cols
            A list of all rel_{} columns within the Assays's dataframe.
        """
        return [i for i in self._df.columns if "rel_" in i]
        # for i in self._df.columns: 
        #     if "rel_" in i:
        #         yield i

    # FUTURE FEATURE HERE
    # def fc(self):
        # some method to also return the fold change columns... 

    @property
    def columns(self):
        return self._df.columns

    def rename_cols(self, cols:dict):
        """
        Renames columns according to a dictionary as key -> value.

        Parameters
        ----------
        cols : dict
            A dictionary specifying old column names (keys) and new colums names (values).
        """
        self._df = self._df.rename(columns = cols)
    
    def n(self):
        """
        Returns 
        ------

        int 
            The number of entries (individual replicates) within the Assay.
        """
        return len( self._df )  # self._length

    def __len__(self):
        return len( self._df )

    def __setitem__( self, key, value ):
        self._df[ key ] = value
        

    def __getitem__( self, key ):
        return self._df[ key ]
    
    def __delitem__( self, key ):
        del self._df[ key ]

    def add_dCt(self, dCt : pd.Series): 
        """
        Adds results from Delta-Ct (first Delta-Ct performed by a `qpcr.Analyser`).

        Parameters
        -----------
        dCt : pandas.Series
            A pandas Series of Delta-Ct values that will be stored in a column `"dCt"`.
            Note, that each `Assay` can, of course, only store one single Delta-Ct column. 
        """
        self._df["dCt"] = dCt
    
    def add_ddCt(self, normaliser_id : str, ddCt : pd.Series):
        """
        Adds results from Delta-Delta-Ct ("normalisation" performed by a `qpcr.Normaliser`).
        These will be stored in a column named `"rel_{normaliser_id}"`. Hence, an Assay can store
        an arbitrary number of Delta-Delta-Ct columns against an arbitrary number of different normalisers. 
        
        Parameters
        ----------
        normaliser_id : str
            The id of the normaliser Assay used to compute the Delta-Delta-Ct values.
        ddCt : pandas.Series
            A pandas Series of Delta-Delta-Ct values.
        """
        name = f"rel_{normaliser_id}"
        self._df[name] = ddCt

    # FUTURE FEATURE HERE
    # some method to add fc columns here...

    def adopt(self, df : pd.DataFrame):
        """
        Adopts an externally computed dataframe as its own.
        This is supposed to be used when setting up new `qpcr.Assay` objects that do not 
        inherit data from one of the `qpcr.Readers`. If you wish to alter an existing `qpcr.Assay` use `force = True`.
        When doing this, please, make sure to retain the proper data structure!

        Parameters
        ----------
        df : pd.DataFrame
            A pandas DataFrame.
        """
        self._df = df

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
        names : list or pd.Series
            The given group names of all replicate groups.
        """
        if "group_name" in self._df.columns: 
            if as_set:
                return list( self._df["group_name"].unique() ) 
            else: 
                return self._df["group_name"]
        else: 
            logger.warning( aw.AssayError( "no_groupnames" ) )
            return None
    
    def groups(self, as_set = True):
        """
        Parameters
        ----------
        as_set : bool
            If `as_set = True` (default) it returns a set (as list without duplicates) 
            of assigned group names for replicate groups.
            If `as_set = False` it returns the full group_name column (including all repeated entries).
        
        Returns
        -------
        groups : list
            The given numeric group identifiers of all replicate groups.
        """
        if "group" in self._df.columns:
            groups = list( self._df["group"].unique() ) if as_set else self._df["group"]
            return groups
        else:
            logger.warning( aw.AssayError("setup_not_grouped") ) 
            return None

    def replicates(self, replicates : (int or tuple or str) = None):
        """
        Either sets or gets the replicates settings to be used for grouping
        Before they are assigned, replicates are vetted to ensure they cover all data entries.

        Parameters
        ----------
        replicates : int or tuple or str
            Can be an `integer` (equal group sizes, e.g. `3` for triplicates), 
            or a `tuple` (uneven group sizes, e.g. `(3,2,3)` if the second group is only a duplicate). 
            Another method to achieve the same thing is to specify a "formula" as a string of how to create a replicate tuple.
            The allowed structure of such a formula is `n:m,` where `n` is the number of replicates in a group and `m` is the number of times
            this pattern is repeated (if no `:m` is specified `:1` is assumed). 
            
            So, as an example, if there are 12 groups which are triplicates, but
            at the end there is one which only has a single replicate (like the commonly measured diluent qPCR sample), we could either specify the tuple
            individually as `replicates = (3,3,3,3,3,3,3,3,3,3,3,3,1)` or we use the formula to specify `replicates = "3:12,1"`. Of course, this works for
            any arbitrary setting such as `"3:5,2:5,10,3:12"` (which specifies five triplicates, followed by two duplicates, a single decaplicate, and twelve triplicates again – truly a dataset from another dimension)...
        """
        if replicates is not None and self._df is not None: 
            # convert a string formula to tuple if one was provided
            if isinstance(replicates, str): 
                replicates = self._reps_from_formula(replicates)
            # vet replicate coverage
            if self._vet_replicates(replicates):
                self._replicates = replicates
            else: 
                logger.critical( aw.AssayError( "reps_dont_cover", n_samples = self.n(), reps = replicates ) )
                raise aw.AssayError( "reps_dont_cover", n_samples = self.n(), reps = replicates )
        return self._replicates

    def group(self, replicates : (int or tuple or str) = None, infer_names = True):
        """
        Groups the data according to replicates-settings specified.

        Parameters
        ----------
        replicates : int or tuple or str
            The replicate settings after which to group the `Assay`. This will just
            get forwarded to the `replicates` method, so there is no need to specify replicates
            here if the replicates method has already been called. 
            See the documentation of the `Assay.replicates` method for more details.

        infer_names : bool
            Try to infer names of replicate groups based on the individual replicate sample identifiers.
            Note that this only works if all replicates have an identical sample name!
        """

        if replicates is not None: 
            self.replicates( replicates )
        
        # generate group and group_names columns
        if isinstance(self._replicates, int):
            groups, group_names = self._make_equal_groups()            
        elif isinstance(self._replicates, tuple):
            groups, group_names = self._make_unequal_groups()
        else:
            if self._identically_named():
                groups = self._infer_replicates()
                group_names = [defaults.group_name.format(i) for i in groups]
            else: 
                e = aw.AssayError("no_reps_inferred", assay = self.id())
                logger.critical( e )
                raise e
        
        # add numeric group identifiers
        self._df["group"] = groups
        self._df["group_name"] = group_names
        
        if infer_names: #and self._names is None:
            # infer group names
            self._infer_names()
            

    def rename(self, names:(list or dict)):
        """
        Replaces the current names of the replicate groups 
        (stored in the "group_name" column).

        Parameters
        ----------
        names : list or dict
            Either a `list` (new names without repetitions) or `dict` (key = old name, value = new name) specifying new group names. 
            Group names only need to be specified once, and are applied to all replicate entries.
        """
        # get new group names based on list (index) or dict (key)
        if isinstance(names, (list, tuple, set)):
            new_names = self._rename_per_index(names)       
        elif isinstance(names, dict):
            new_names = self._rename_per_key(names)
        else:
            logger.error( aw.AssayError("no_groupname_assignment", names = names ) )
            return

        # update "group_name"
        self._df["group_name"] = new_names
        self._renamed = True

    def ignore(self, entries:tuple, drop = False):
        """
        Remove lines based on index from the dataframe.
        This is useful when removing corrupted data entries.

        Parameters
        ----------
        entries : tuple
            Tuple of row indices from the dataframe to drop.
        drop : bool 
            If True the provided entries will be entirely removed from the 
            dataset. If False, ignore entries will be set to NaN. 
        """
        if drop:
            self._df = self._df.drop(index = list(entries))
        else: 
            Cts = np.array( self.Ct )
            Cts[ entries ] = np.nan
            self._df["Ct"] = Cts
    
    def _reps_from_formula(self, replicates):
        """
        Generates a replicate tuple from a string formula. 
        See the docstring of `replicates()` for more info on the formula.

        Example:
        "3:4,1:4,2:3,9" -> (3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 9)
        """

        # split the formula and adjust standard formatting
        replicates = replicates.split(",")
        replicates = [i + ":1" if ":" not in i else i for i in replicates]

        # convert to numeric values and extend
        replicates = [np.array(i.split(":"), dtype = int) for i in replicates]
        replicates = [np.tile(i[0], i[1]) for i in replicates]
        
        # generate replicate tuple
        replicates = np.concatenate(replicates)
        replicates = tuple(replicates)
        
        return replicates

    def _infer_replicates(self):
        """
        Infers the replicate groups based on the replicate ids in case all replicates of the same group have the same name.
        """
        names = self._df[raw_col_names[0]]
        names_set = names.unique()
        groups = [i for i in range(len(names_set))]
        for name, group in zip(names_set, groups):
            names = names.replace(name, group)
        
        indices = np.array(names, dtype = int)
        return indices

    def _infer_names(self):
        """
        Infers replicate group names from the given replicate identifier column
        """
        if self._identically_named():
            self._df["group_name"] = self._df[raw_col_names[0]]
        elif self._names is None: 
            logger.warning( aw.AssayError("groupnames_not_inferred") )

    def _identically_named(self):
        """
        Checks if all replicates in the same group have the same name / id
        It checks simply the first group, if that is identical then it's fine.
        """
        if "group" not in self._df.columns:
            names = self._df[raw_col_names[0]]
            names_set = names.unique()
            # names_set = aux.sorted_set(names)
            first_name = names_set[0]
            group0 = self._df.query(f"{raw_col_names[0]} == '{first_name}'")[raw_col_names[0]]
            entries = len(group0)
            all_identical = entries > 1                
        else: 
            group0 = self._df.query("group == 0")[raw_col_names[0]]
            all_identical = all(group0 == group0[0])
        return all_identical

    def _rename_per_key(self, names):
        """
        Generates new name list based on current names in "group_name" and uses string.replace()
        to update groupnames, based on key (old name) : value (new name) indexing. 
        Before applying it checks if all groups are covered by new names
        """
        current_names = self.names()
        # current_names = aux.sorted_set(self._df["group_name"])
        all_groups_covered = len(names) == len(current_names)
        if all_groups_covered:
            current_names = list(self._df["group_name"])
            new_names = "$".join(current_names)
            for old_name, new_name in names.items():
                new_names = new_names.replace(old_name, new_name)
            new_names = new_names.split("$")       
            return new_names
        else:
            e = aw.AssayError("groupnames_dont_colver", current_groups = current_names, new_received = names)
            logger.critical( e )
            raise e 

    def _rename_per_index(self, names):
        """
        Generates new name list based on current names in "group_names" and uses string.replace()
        to update groupnames to new names based on index (using a the order 
        of groups as is currently present in "group_name"). 
        """
        current_names_set = self.names()
        # current_names_set = aux.sorted_set(self._df["group_name"])
        all_groups_covered = len(names) == len(current_names_set)
        if all_groups_covered:
            current_names = list(self._df["group_name"])
            new_names = "$".join(current_names)
            names = list(names)
            for old_name, new_name in zip(current_names_set, names):
                new_names = new_names.replace(old_name, new_name)
            new_names = new_names.split("$")
            return new_names
        else:
            e = aw.AssayError("groupnames_dont_colver", current_groups = current_names, new_received = names)
            logger.critical( e )
            raise e 

    def _make_unequal_groups(self):
        """
        Returns two lists of [0,0,0,1,1,1] and 
        [Group0, Group0, Group0, Group1,...] 
        to cover all data entries.
        (this function works with a tuple for replicate group sizes)
        """
        groups = []
        group_names = []
        for rep, idx in zip(self._replicates, range(len(self._replicates))): 
            groups.extend([idx] * rep)
            group_names.extend([ defaults.group_name.format(idx) ] * rep)
        return groups, group_names

    def _make_equal_groups(self):
        """
        Returns two lists of [0,0,0,1,1,1] and 
        [Group0, Group0, Group0, Group1,...] 
        to cover all data entries.
        (this function works with an integer group size, 
        assuming all groups have the same size)
        """
        assays = self.n()
        groups = []
        group_names = []
        slices = range(int(assays / self._replicates))
        for i in slices:
            groups.extend([i] * self._replicates)
            group_names.extend([ defaults.group_name.format(i) ] * self._replicates)
        return groups, group_names

    def _vet_replicates(self, replicates : (int or tuple)):
        """
        Checks if provided replicates will place all data entries into a group
        returns True if all replicates are covered, False if not...
        """
        current_entries = self.n()
        verdict = None

        # for INT -> modulo will be 0 if all replicates are covered
        # for TUPLE -> sum(replicates) should cover all replicates...

        if isinstance(replicates, int):
            verdict = True if current_entries % replicates == 0 else False
        elif isinstance(replicates, tuple): 
            verdict = True if sum(replicates) == current_entries else False
        
        if verdict is None:
            e = aw.AssayError( "reps_could_not_vet", reps = replicates ) 
            logger.error( e )
            raise e 

        return verdict