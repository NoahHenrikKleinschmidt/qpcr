"""
This is the ``qpcr.Normaliser`` class whose function is to compute Fold-Changes of Delta-Ct values
from assays and normalisers using ``qpcr.Assay`` objects.


Modes of *normalisation*
========================

The ``qpcr.Normaliser`` supports custom functions for normalisation. However, it also comes with three built-in methods to normalise sample-assays against normalisers.
These are accessible via the `mode` argument of the `qpcr.Normaliser.normalise` method, which can be set to ``"pair-wise"`` (default), ``"combinatoric"``, or ``"permutative"``. 

The default option `"pair-wise"` is computationally the fastest and will rigidly normalise replicates against their corresponding partner from the normaliser. I.e. first against first, 
second against second, etc. This mode is appropriate for multiplex qPCR experiments. 

For qPCR reactions that were pipetted individually, there is not reason to strictly only pair first
with first, second with second etc. For these cases there are other two options `"combinatoric"` and ``"permutative"``. ``"combinatoric"`` normalisation will calculate all possible group-wise 
combinations of a sample-assay replicate with all available normaliser replicates of the same group. I.e. first against first, and against second, etc. This will generate :math:`n^2` values where 
:math:`n` is the number of replicates within a group. This mode is appropriate for small-scale datasets but will substantially increase computation times for larger datasets.   

``"permutative"`` on the other hand will reflect the equivalence of replicates within a group through random permutations wihtin the normaliser replicates. Hence, 
first may be normalised against first, or second, etc. This normalisation method may by used iteratively to increase the accuracy. Replacement during permutations are allowed (although disabled by default).
If replacement is desired by the user, the probability of each replicate to be chosen will be weighted based on a fitted normal distribution. This method is appropriate for larger datasets for which combinatoric normalisation
is not desired. 


Normalising data
================

Setting up a ``qpcr.Normaliser`` is really easy and we see it in virtually every GitHub tutorial.

.. code-block::
    
    normaliser = qpcr.Normaliser()

    # and now directly pipe the data through
    results = normaliser.pipe( some_assays, some_normalisers )

Alternatively we can also directly use the ``qpcr.normalise`` function that will call on a Normaliser for us.

.. code-block::
    
    results = qpcr.normalise( some_assays, some_normalisers )


Preprocessing "normalisers"
---------------------------

The ``qpcr.Normaliser`` can work with multiple ``qpcr.Assay`` as "normaliser-assays". However, it requires for computation a single set of numbers to compute fold changes.
Therefore, the Normaliser first performs some *pre-processing* of all normaliser assays it receives. By default it will **compute their mean** and use this to compute fold changes.
However, the ``qpcr.Normaliser`` is equipped with a method ``prep_func`` which allows you to pass a custom function for preprocessing.

Normalisation
-------------

By default the ``qpcr.Normaliser`` will compute normalised fold changes by dividing the assays-of-interst by the pre-processed normaliser. 
As a side note here, the ``qpcr.Analyser`` already stores the Delta-Ct values it computes as :math:`efficiency^{-\Delta Ct}`. 
Hence, by default the Delta-Ct values stored by ``qpcr.Assay``s after having been "analysed" are exponentials.
However, using the method ``norm_func`` also a custom normalisation function can be specified.

"""

import qpcr._auxiliary.warnings as aw
import qpcr.defaults as defaults
import qpcr._auxiliary as aux

from qpcr.main.Assay import Assay
from qpcr.main.Results import Results

import pandas as pd
import numpy as np
import scipy.stats as scistats


from copy import deepcopy
import logging 

logger = aux.default_logger()

class Normaliser(aux._ID):
    """
    Handles the second step in Delta-Delta-Ct (normalisation against normaliser assays).
    
    Note
    -----
    This requires that all have been analysed in the same way before!
    """
    __slots__ = ["_Normalisers", "_Assays", "_Results", "_normaliser", "_prep_func", "_norm_func", "_norm_func_is_set"]

    def __init__(self):
        super().__init__()

        self._Normalisers = []
        self._Assays = []
        
        self._Results = Results()
        
        # the actually used normaliser 
        # which will be a pre-processed version
        # from all supplied normalsier assays
        # by default this will be the averaged version
        # of all normaliser assays (but this can be changed)  
        self._normaliser = None

        # setup defaults
        self._prep_func = self._preprocess_normalisers
        self._norm_func = self._divide_by_normaliser

        # store the state if a custom norm func was provided
        self._norm_func_is_set = False

    def __str__(self):
        s = f"""
Normaliser: \t{self.id()}
Prep.Function: \t{self._prep_func.__name__}
Prep.Function: \t{self._norm_func.__name__}
        """.strip()
        if not self._Results.is_empty:
            _length = len( str(self._Results._df).split("\n")[0] ) 
            s = f"{'-' * _length}\n{s}\n{'-' * _length}\n{self._Results._df}\n{'-' * _length}"
        return s

    def __repr__(self):
        prep = self.prep_func.__name__
        norm = self.norm_func.__name__
        return f"Normaliser({prep=}, {norm=})"

    def clear(self):
        """
        Will clear the presently stored results
        """
        self._Results = Results()

    
    def prune(self, assays = True, normalisers = True, results = True):
        """
        Will clear assays, normalisers, and/or results
        assays : bool
            Will clear any sample assays if True (default).
        
        results : bool
            Will clear any computed results if True (default).
        
        normalisers : bool
            Will clear any normalisers if True (default).
        """
        if assays:
            self._Assays = []
        if normalisers: 
            self._Normalisers = []
        if results: 
            self.clear()

    def get(self, copy=False):
        """
        Parameters
        ----------
        copy : bool
            Will return a deepcopy of the Results object if `copy = True` (default is `copy = False`).
        
        Returns
        -------
        Results : qpcr.Results
            A `qpcr.Results` object containing the normalised dataframe
        """
        if copy: 
            return deepcopy(self._Results)
        return self._Results
    
    def link(self, assays:(list or tuple) = None, normalisers:(list or tuple) = None):
        """
        Links either normalisers or assays-of-interest `qpcr.Assay` objects coming from the same `qpcr.Analyser`.

        Parameters
        ----------
        assays : list or tuple or qpcr.Analyser
            A list of `qpcr.Assay` objects coming from a `qpcr.Analyser`. These assays will be normalised against a normaliser.
        
        normalisers : list or tuple
            A list of `qpcr.Assay` objects coming from a `qpcr.Analyser`. These assays will be used as normalisers. These will be
            combined into one single pseudo-normaliser which will then be used to normalise the assays. The method of 
            combining the normalisers can be specified using the `qpcr.Normaliser.prep_func` method.
        """
        # convert to lists (since the _link methods really want lists)
        if assays is not None and not isinstance(assays, (list, tuple)): 
            assays = [assays]
        self._link_assays(assays)

        if normalisers is not None and not isinstance(normalisers, (list, tuple)): 
            normalisers = [normalisers]
        self._link_normaliser(normalisers)
    
    def prep_func(self, f = None):
        """
        Sets any defined function for combined normaliser pre-processing.
        If no `f` is provided, it returns the current `prep_func`.

        Parameters
        ----------
        f : function
            The function may accept one list of `qpcr.Assay` objects, and must return 
            either an `qpcr.Assay` object directly or a `pandas.Dataframe` (that will be migrated to an `qpcr.Assay`).
            The returned dataframe must contain a `"dCt"` column which stores the delta-Ct values ultimately used as 
            "normaliser assay".  
        """
        
        if aux.same_type(f, aux.fileID):
            self._prep_func = f
        elif f is None:
            return f
        else: 
            aw.HardWarning("Normaliser:cannot_set_prep_func", func = f)

    def norm_func(self, f = None):
        """
        Sets any defined function to perform normalisation of assays against normalisers.
        If no `f` is provided, it returns the current `norm_func`.

        Parameters
        ----------
        f : function
            The function may accept two `qpcr.Assay` objects (named `assay` and `normaliser` which will be forwarded from the `qpcr.Normaliser`). 
            The function may also accept one `pandas.DataFrame` (named `df`) containing two numeric columns of delta-Ct values from a sample-assay (named "s") and a normaliser-assay (named "n"),
            as well as a group identifier column (named "group"). 
            Whatever inputs it works with, it must return a named numeric `pandas.Series` of the same length as entries in the Assays' dataframes. 
            
            ##### Note
            Support for the dataframe direct usage will be dropped at some point in the future. 

            By default `s/n` is used, where `s` is a column of sample-assay deltaCt values, 
            and `n` is the corresponding `"dCt"` column from the normaliser.

        """
        if aux.same_type(f, aux.fileID):
            self._norm_func = f
            self._norm_func_is_set = True
        elif f is None:
            return f
        else: 
            aw.HardWarning("Normaliser:cannot_set_norm_func", func = f)

    def normalise(self, mode = "pair-wise", **kwargs):
        """
        Normalises all linked assays against the combined pseudo-normaliser 
        (by default, unless a custom `prep_func` has been specified), 
        and stores the results in a new Results object.

        Parameters
        ----------
        mode : str
            The normalisation mode to use. This can be either `pair-wise` (default), 
            or `combinatoric`, or `permutative`.
            `pair-wise` will normalise replicates only by their partner (i.e. first against first, 
            second by second, etc.). `combinatoric` will normalise all possible combinations of a replicate 
            with all partner replicates of the same group from a normaliser (i.e. first against first, then second, then third, etc.).
            This will generate `n^2` normalised Delta-Delta-Ct values, where `n` is the number of replicates in a group.
            `permutative` will scramble the normaliser replicates randomly and then normalise pair-wise. This mode supports 
            a parameter `k` which specifies the times this process should be repeated, thus generating `n * k` normalised Delta-Delta-Ct values.
            Also, through setting `replace = True` replacement may be allowed during normaliser scrambling.
            Note, this setting will be ignored if a custom `norm_func` is provided.

        **kwargs
            Any additional keyword arguments that may be passed to a custom 
            `norm_func` and `prep_func` (both will receive the kwargs!).
        """
        if self._normaliser is None: 
            self._normaliser = self._prep_func(self._Normalisers, **kwargs)
            self._vet_normaliser()

        if self._Assays == [] or self._normaliser is None:
            e = aw.NormaliserError( "no_data_yet" )
            logger.error( e )
            
        # check which kind of norm_func we should use
        if not self._norm_func_is_set:
            # pair-wise is already default so we don't check...
            if mode == "combinatoric":
                self._norm_func = self._tile_normalise
                tiled = deepcopy(self._Assays[0])
                tiled.tile()
                self._Results.setup_cols( tiled )
                del tiled 
            elif mode == "permutative":
                self._norm_func = self._permutate_normalise
                n = aux.from_kwargs("k", 1, kwargs)
                tiled = deepcopy(self._Assays[0])
                tiled.stack(n)
                self._Results.setup_cols( tiled )
                del tiled 
            elif mode == "pair-wise":
                self._Results.setup_cols( self._Assays[0] )
        else:
            self._Results.setup_cols( self._Assays[0] )

        # perform normalisation for each assay 
        for assay in self._Assays:

            # get data
            # assay_df = assay.get()

            # apply normalisation (delta-delta-Ct)
            normalised = self._norm_func_wrapper(
                                                    assay = assay, 
                                                    normaliser = self._normaliser, 
                                                    **kwargs
                                            )

            # # store results in _Results
            # self._store_to_Results(assay, normalised)

            # and store results also in the Assay itself
            assay.add_ddCt( self._normaliser.id(), normalised )

            # and store to results
            self._Results.add_ddCt(assay)

    def pipe( self, assays : list, normalisers : list, mode = "pair-wise", **kwargs ):
        """
        A wrapper for `Normaliser.link` and `Normaliser.normalise`

        Parameters
        ----------
        assays : list 
            A list of  `qpcr.Assay` objects.
        normalisers : list
            A list of `qpcr.Assay` objects.
        mode : str
            The normalisation mode to use. This can be either `pair-wise` (default), 
            or `combinatoric`, or `permutative`.
            `pair-wise` will normalise replicates only by their partner (i.e. first against first, 
            second by second, etc.). `combinatoric` will normalise all possible combinations of a replicate 
            with all partner replicates of the same group from a normaliser (i.e. first against first, then second, then third, etc.).
            This will generate `n^2` normalised Delta-Delta-Ct values, where `n` is the number of replicates in a group.
            `permutative` will scramble the normaliser replicates randomly and then normalise pair-wise. This mode supports 
            a parameter `k` which specifies the times this process should be repeated, thus generating `n * k` normalised Delta-Delta-Ct values.
            Also, through setting `replace = True` replacement may be allowed during normaliser scrambling.
            Note, this setting will be ignored if a custom `norm_func` is provided.

        **kwargs
            Any additional keyword arguments that may be passed to a custom 
            `norm_func` and `prep_func` (both will receive the kwargs!).
        Returns
        -------
        results : qpcr.Results
            A `qpcr.Results` object of the assembled results.
        """
        self.link( assays = assays, normalisers = normalisers )
        self.normalise( mode = mode )
        return self.get()


    def _vet_normaliser(self):
        """
        Checks if the normaliser is already a qpcr.Assay object, and if not
        convert it to one. 
        """
        if not isinstance(self._normaliser, Assay):
            tmp = Assay( df = self._normaliser )          
            tmp.id("combined_normaliser")
            self._normaliser = tmp
            

    def _store_to_Results(self, assay, normalised):
        """
        Stores computed Delta-delta-ct (normalisation)
        into the _Results object.
        """
        column_name = f"{assay.id()}_rel_{self._normaliser.id()}"
        normalised = normalised.rename(column_name)
        self._Results.add(normalised)


    def _norm_func_wrapper(self, assay, normaliser, **kwargs):
        """
        The wrapper that will apply the _norm_func to the sample and normaliser dataframes and return a pandas series of normalised values
        """
        # for double normalised we want the same columns as dct and norm...

        # FUTURE DROP HERE
        # In the future we will not be creating the tmp_df 
        # directly anymore, but intead will only pass assay and normaliser. 

        sample_dCt = assay.dCt
        groups = assay.groups( as_set = False )
        norm_dCt = normaliser.dCt

        tmp_df = pd.DataFrame( dict( group = groups, s = sample_dCt, n = norm_dCt )  )

        results = self._norm_func(df = tmp_df, assay = assay, normaliser = normaliser, **kwargs)

        # this is the old call from before factoring out to Assays 
        # dCt_col, norm_col = self._prep_columns(sample_assay, dCt_col, norm_col)

        # tmp_df = normaliser.join(sample_assay, lsuffix="_s")
        # # tmp_df = sample_assay.join(normaliser, rsuffix = "_n")
        # results = self._norm_func(tmp_df[[dCt_col, norm_col]], **kwargs)
        return results

    def _tile_normalise(self, assay, normaliser, **kwargs):
        """
        Normalises assays and normalisers group wise, iteratively normalising
        each individual replicate against all replicates from the normaliser.
        This generates `n**2` normalised Delta-Delta-Ct values where `n` is the
        group size. 
        """

        # get the untiled data
        adf = assay.get()
        groups = assay.groups()
        ndf = normaliser.get()
        
        # get the column to draw data from
        col = aux.from_kwargs("col", "dCt", kwargs)

        # tile the assay
        assay.tile() 

        # generate results array
        ddCts = np.zeros( len( assay.get() ) )
       
        # now compute ddCt
        idx = 0
        for group in groups: 

            a_dCt = adf.query( f"group == {group}" )[ col ].to_numpy()
            n_dCt = ndf.query( f"group == {group}" )[ col ].to_numpy()
            
            for a in a_dCt:
                for n in n_dCt:

                    r = a / n
                    
                    try: 
                        ddCts[idx] = r 
                    except: 
                        break 

                    idx += 1

        ddCts = pd.Series( ddCts )       
        ddCts.name = "ddCt"
        return ddCts
        
    
    def _permutate_normalise(self, assay, normaliser, **kwargs):
        """
        Scrambles randomly the normaliser's replicate values group-wise 
        """
        # get the untiled data
        adf = assay.get()
        groups = assay.groups()
        ndf = normaliser.get()
        
        # get the column to draw data from
        col = aux.from_kwargs("col", "dCt", kwargs)

        # stack the assay
        # get the number of permutations to perform
        n = aux.from_kwargs("k", 1, kwargs)  
        assay.stack(n) 

        # get replace argument for random choice
        replace = aux.from_kwargs("replace", False, kwargs)

        # generate results array
        ddCts = np.zeros( len( assay.get() ) )
       
        # now compute ddCt
        idx = 0
        for group in groups: 
            for i in range(n):
                a_dCt = adf.query( f"group == {group}" )[ col ].to_numpy()
                n_dCt = ndf.query( f"group == {group}" )[ col ].to_numpy()
                
                # randomly scramble the normaliser replicates
                np.random.seed( defaults.default_seed )
                if replace:
                    # in case of replacement, we generate a normal 
                    # distribution from the replicate values and weigh
                    # the random.choice with the given probabilities of the values 
                    # being chosen. Since the probabilities do not themselves correspond to 1
                    # we normalise against the limited available subset to generate a probabilities
                    # array that sums up to 1.
                    mu, sd = scistats.norm.fit(n_dCt)
                    probs = scistats.norm.pdf(n_dCt, loc = mu, scale = sd)
                    probs = probs / np.sum(probs)
                    n_dCt = np.random.choice( n_dCt, size = n_dCt.size, replace = replace, p = probs )
                else:
                    n_dCt = np.random.choice( n_dCt, size = n_dCt.size, replace = replace )
                length = a_dCt.size
                r = a_dCt / n_dCt
                        
                try: 
                    ddCts[idx : idx + length] = r 
                except: 
                    break 

                idx += length

        ddCts = pd.Series( ddCts )       
        ddCts.name = "ddCt"
        return ddCts

    def _divide_by_normaliser(self, df, **kwargs):
        """
        Performs normalisation of sample s against normaliser n
        s and n are specified as two pandas dataframe columns
        Note, that the dataframe must ONLY contain these two columns, first the dCt sample, then the normaliser!
        (default _norm_func)
        """
        s, n = df["s"], df["n"]
        # this is the old call from before factoring out to Assays
        # dCt_col, norm_col = df.columns
        # s, n = df[dCt_col], df[norm_col]
        return s / n

    def _link_assays(self, assays):
        """
        Links any provided assays and checks their datatype in the process...
        """
        if assays is not None:
            for assay in assays: 
                # if aux.same_type(assay, Assay()):
                # for some reason the isinstance is not working again...
                # UPDATE: that's probably due to the way we import it without
                # the whole module stuff in front. Anyway, this one works and
                # it should be fine for our purposes although it's a bit hacky...
                if type(assay).__name__ == "Assay" :
                    self._Assays.append(assay)
                else: 
                    e = aw.NormaliserError( "unknown_data", s = assay )
                    logger.error( e )
                
    def _link_normaliser(self, normalisers):
        """
        Checks if normaliser is provided and has proper datatype to be added...
        """
        if normalisers is not None:
            for normaliser in normalisers:
                # if aux.same_type(normaliser, Assay()):
                if type(normaliser).__name__ == "Assay" :
                    self._Normalisers.append(normaliser)

                else: 
                    e = aw.NormaliserError( "norm_unknown_data", s = normaliser )
                    logger.error( e )

    def _preprocess_normalisers(self, *args, **kwargs):
        """
        Averages the provided normalisers row-wise for all normalisers into a 
        single combined normaliser, that will be stored as a new Assay object.
        """

        # initialise new Results to store the dCt values form all normalisers
        combined = Results()
        
        # setup names using the first normaliser
        combined.setup_cols( self._Normalisers[0] )
        combined.adopt_id(self._Normalisers[0])

        # now add all dCt columns from all normalisers
        for norm in self._Normalisers:
            combined.add_dCt(norm)

        # now generate the combined normaliser
        combined_normaliser = self._average(combined)
        combined_normaliser = combined_normaliser.rename("dCt")
        combined.add(combined_normaliser)

        logger.debug( combined )
        
        # now assemble the normaliser into a qpcr.Assay
        combined = combined.get()
        normaliser = Assay( df = combined, replicates = self._Normalisers[0].replicates(), group_names = self._Normalisers[0].names() )

        self._normaliser = normaliser  
        self._update_combined_id()
        logger.debug( f"combined normaliser id: {self._normaliser.id()}" )

        # forward combined_id to self and _Results 
        self.adopt_id(self._normaliser)
        self._Results.adopt_id(self._normaliser)

        return self._normaliser

    def _update_combined_id(self):
        """
        Generates a new id based on all normaliser ids,
        joining them as a+b+c,...
        """
        ids = [N.id() for N in self._Normalisers]
        ids = "+".join(ids)
        self._normaliser.id_reset()
        self._normaliser.id(ids)
        

    def _average(self, combined):
        """
        Averages row-wise all Normaliser entries and 
        generates a series of their per-row means
        (default preprocess_normalisers function)
        """
        tmp_df = combined.get()

        # drop group as it is a numeric column 
        # and would otherwise skew the average
        if "group" in tmp_df.columns:
            tmp_df = tmp_df.drop(columns = ["group"]) 
        tmp_df = tmp_df.mean(axis = 1, numeric_only = True)

        return tmp_df


__generic__Normaliser__ = Normaliser()
"""The default Normaliser"""

def normalise( assays : list, normalisers : list, mode : str = "pair-wise" ):
    """
    Computes DeltaDeltaCt values (Normalised fold changes) for a number of assays-of-interest and normalisers.

    Note
    ----
    This is essentially a Wrapper for a default `Normaliser` and it's `pipe` method.
    If you require non-default behaviour, please, set up your own `Normaliser`. Also,
    since this function will setup a new Normaliser each time it is called, it is
    less efficient than the conventional way of directly using the `pipe` method.

    Parameters
    ----------
    assays : list
        A list `qpcr.Assay` objects that are "assays of interest".

    normalisers : list
        A list of `qpcr.Assay` objects that are "normalisers".

    mode : str
            The normalisation mode to use. This can be either `pair-wise` (default), 
            or `combinatoric`, or `permutative`.
            `pair-wise` will normalise replicates only by their partner (i.e. first against first, 
            second by second, etc.). `combinatoric` will normalise all possible combinations of a replicate 
            with all partner replicates of the same group from a normaliser (i.e. first against first, then second, then third, etc.).
            This will generate `n^2` normalised Delta-Delta-Ct values, where `n` is the number of replicates in a group.
            `permutative` will scramble the normaliser replicates randomly and then normalise pair-wise. This mode supports 
            a parameter `k` which specifies the times this process should be repeated, thus generating `n * k` normalised Delta-Delta-Ct values.
            Also, through setting `replace = True` replacement may be allowed during normaliser scrambling.
            Note, this setting will be ignored if a custom `norm_func` is provided.

    Returns
    -------
    results : Results
        A `qpcr.Results` object containing the computed results.
    """
    return __generic__Normaliser__.pipe( assays, normalisers )
