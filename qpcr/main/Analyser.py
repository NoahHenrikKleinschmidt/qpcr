"""
This is the `qpcr.Analyser` whose function is to perform dataset-internal 
normalisation to compute the first-step Delta-Ct Values within an `qpcr.Assay` object.



Computing Delta-Ct values
=========================

Setting up a ``qpcr.Analyser`` is really easy and we see it in virtually every GitHub tutorial.

.. code-block:: python
    
    analyser = qpcr.Analyser()

    # and now directly pipe the data through
    some_assays = analyser.pipe( some_assays )

Alternatively we can also directly use the ``qpcr.delta_ct`` function that will call on a Analyser for us.

.. code-block:: python
    
    some_assays = qpcr.delta_ct( some_assays )


Delta-Ct values
---------------
The computed values are stored in the respective ``qpcr.Assay`` dataframe into a column called ``"dCt"``.
By default, the ``qpcr.Analyser`` will compute the Delta-Ct values already in exponential form. I.e. as :math:`efficiency^{-\Delta Ct}`.
This behaviour can be changed by changing the applied function using the provided ``func`` method. 
Check out `this tutorial about custom anchors <https://github.com/NoahHenrikKleinschmidt/qpcr/blob/main/Examples/4_custom_anchor.ipynb>`_  which will give you an idea of the flavour of editing the ``qpcr.Analyser``. 

"""

import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import qpcr._auxiliary.warnings as aw
import pandas as pd
from qpcr.main.Assay import Assay

logger = aux.default_logger()

raw_col_names = defaults.raw_col_names
class Analyser(aux._ID):
    """
    Performs Single Delta-Ct (first normalisation 
    within dataset against the `anchor`) 
    """
    __slots__ = ["_anchor", "_ref_group", "_ref_group_col", "_deltaCt_function", "_Assay"]

    def __init__(self):
        super().__init__()
        self._Assay = None

        # default settings
        self._anchor = "first"
        self._ref_group = 0
        self._ref_group_col = "group"       # used in case of "mean" anchor where the ref_group must be located either from a numeric (group) or string (group_name) id
        
        self._deltaCt_function = self._get_deltaCt_function(exp = "exponential")

        # Efficiency setting in the Analyser is disabled now...
        # self._eff_src = self._Assay         # By default use the efficencies stored in the Assay directly 
        #                                     # this is set to self if self.efficiency() is called...
        # self._efficiency = 1                # the formal effiency in percent
        # self._eff = 2 * self._efficiency    # the actual doubplciation factor used for calculation


    def get(self):
        """
        Returns 
        -------
        Assay : qpcr.Assay
            The analysed `qpcr.Assay` object that contains now deltaCT values.
        """
        return self._Assay

    def link(self, Assay:Assay):
        """
        Links a `qpcr.Assay` object to the Analyser.
       
        Parameters
        ----------
        Assay : qpcr.Assay
            A `qpcr.Assay` object containing data.
        
        """
        self._Assay = Assay

        # Efficiency setting in the Analyser is disabled now...
        # # check if no efficiency was specifically set for the Analyser
        # # and if not so, use the Assay's efficiency...
        # if self._eff_src is not self: 
        #     self._eff_src = self._Assay

    def pipe(self, Assay:Assay, **kwargs):
        """
        A quick one-step implementation of link + DeltaCt.

        Note
        ----
        This is the suggested application of the `qpcr.Analyser`.


        Parameters
        ----------
        Assay : qpcr.Assay
            A `qpcr.Assay` object to be linked to the Analyser for DeltaCt computation.
        
        **kwargs
            Any additional keyword arguments to be passed to the `DeltaCt()` method.

        Returns 
        -------
        assay : qpcr.Assay
            The same `qpcr.Assay` with computed Delta-Ct values. 

        """
        if isinstance( Assay, list ):
            return [ self.pipe( A ) for A in Assay ]
        else:
            self.link(Assay)
            self.DeltaCt(**kwargs)
            assay = self._Assay
            return assay

    # Efficiency setting in the Analyser is disabled now...
    # def efficiency(self, e:float = None):
    #     """
    #     Sets an efficiency factor for externally calculated qPCR amplification efficiency.
    #     By default `efficiency = 1` (100%) is assumed.

    #     Note
    #     ----
    #     By default the efficiency is now (`qpcr 3.2.0`) handled by the `qpcr.Assay` objects directly.
    #     The Analyser will directly read the efficiency from the `qpcr.Assay`s. However,
    #     it is still possible to set the efficiency via the `qpcr.Analyser` in this fashion.

    #     Parameters
    #     ----------
    #     e : float
    #         An amplification efficiency factor. Default is `e = 1`, 
    #         which is then treated as `eff = 2 * e`, so `e = 1` corresponds to true duplication
    #         each cycle.

    #     Returns
    #     -------
    #     efficiency : float
    #         The current efficiency used.

    #     """
    #     # deprecation warning
    #     aw.SoftWarning( "blank", msg = "The use of efficiency() from the Analyser is deprecated and will be dropped in a future version! Please, set efficiencies directly in the Assay." )
    #     if isinstance(e, (int, float)):
    #         self._efficiency = float( e )
    #         self._eff = 2 * self._efficiency
    #         # and update the efficiency source 
    #         # to self instead of the assay...
    #         self._eff_src = self 
    #     return self._efficiency

    def anchor(self, anchor : (str or float or function) = None, group : (int or str) = 0):
        """
        Sets the anchor for DeltaCt for internal normalisation.

        Parameters
        ----------
        anchor : str or float or function
            The internal anchor for normalisation.
            This can be either `"first"` (default, the very first dataset entry),
            `"mean"` (mean of the reference group), 
            `"grouped"` (first entry for each replicate group), 
            any specified numeric value (as `float`), 
            or a `function` that will calculate the anchor and returns a single numeric value. 
            If you wish to use a function to compute the anchor, you can access the dataframe stored by the `qpcr.Assay` that is being analysed through the 
            `data` argument. `data` will be automatically forwarded to your custom anchor-function, unless you specify it directly. Please, make sure your function can handle
            `**kwargs` because any kwargs supplied during `DeltaCt()`-calling will be passed per default to both the anchor-function and DeltaCt-function. 
        group : int or str
            The reference group identifier. This can be either the numeric identifier or the `group_name`. This is only used for `anchor = "mean"`. 
            By default the first group is assumed. 
        
        Returns
        -------
        anchor 
            The currently selected `anchor`.
        ref_group
            The current reference group.
        """
        if anchor is not None:
            self._anchor = anchor
        if group != self._ref_group:
            self._ref_group = group
        
        # update column where to search for ref_group identifier
        if isinstance(group, str):
            self._ref_group_col = "group_name"
        else: 
            self._ref_group_col = "group"
        
        return self._anchor, self._ref_group


    def func(self, f:(str or function)):
        """
        Sets the function to be used for DeltaCt (optional)

        Parameters
        ----------
        f : str or function
            The function to be used for DeltaCt computation. Pre-defined functions are 
            either `"exponential"` (which uses  `dCt = eff ** ( -(s-r) )`, default), or `"linear"` 
            (uses `dCt = s-r`), where `s` is any replicate entry in the dataframe and `r` is the anchor. `eff = 2 * efficiency` is the 
            numeric duplication factor (default assumed `efficiency = 1`).
            It is also possible to assign any defined function that accepts one `float` Ct value `s` (1st!) and anchor `r` value (2nd!), 
            alongside any `kwargs` (which will be forwarded from DeltaCt()...). It must return also a single `float`. 
        """
        if f in ["exponential", "linear"]:
            # f = True if f == "exponential" else False
            self._deltaCt_function = self._get_deltaCt_function(f)
        elif type(f) == type(aux.fileID):
            self._deltaCt_function = f
        else:
            e = aw.AnalyserError( "cannot_set_func", func = f )
            logger.critical( e )
            raise e

    def DeltaCt(self, **kwargs):
        """
        Calculates Delta-Ct for all groups within the dataframe.
        Any specifics such as `anchor` or `func` must have already been 
        set using the respective methods prior to calling `DeltaCt()`!

        Parameters
        ----------
        **kwargs
            Any additional keyword arguments that a custom DeltaCt function may require.
        """

        predefined = {
                        "first"     : self._DeltaCt_first_anchored,
                        "grouped"   : self._DeltaCt_grouped_anchored,
                        "mean"      : self._DeltaCt_mean_anchored,
                    }

        if self._anchor in predefined.keys():
        
            deltaCt_func = predefined.get( self._anchor )
            deltaCt_func(
                                            self._deltaCt_function, 
                                            **kwargs
                                    )
        
        elif isinstance(self._anchor, (float, int)): 
            self._DeltaCt_externally_anchored(
                                                self._anchor, 
                                                self._deltaCt_function, 
                                                **kwargs
                                            )
        elif type(self._anchor) == type(aux.fileID):
            self._DeltaCt_function_anchor(
                                            self._anchor, 
                                            self._deltaCt_function, 
                                            **kwargs
                                        )

    def _DeltaCt_function_anchor(self, anchor_function, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using a function as anchor
        """
        df = self._Assay.get()
        # update a "data" argument into the 
        # kwargs for the anchor_function
        if "data" not in kwargs.keys(): 
            kwargs.update(dict(data = df))
        anchor = anchor_function(**kwargs)

        # apply deltaCt_function
        Ct = raw_col_names[1]
        dCt = df[Ct].apply(deltaCt_function, r = anchor, **kwargs)
        dCt.name = "dCt"
        # store results
        self._Assay.add_dCt(dCt)

    def _DeltaCt_mean_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the mean of the reference group as anchor
        """
        df = self._Assay.get()

        # get  reference group
        ref_query = "{} == {}" if isinstance(self._ref_group, int) else "{} == '{}'" 
        ref = df.query(
                        ref_query.format(  self._ref_group_col, self._ref_group  )
                    )
        # get Ct values from ref group and make anchor
        ref = ref[raw_col_names[1]]
        anchor = ref.mean()
        
        # apply DeltaCt function
        Ct = raw_col_names[1]
        dCt = df[Ct].apply(deltaCt_function, r = anchor, **kwargs)
        dCt.name = "dCt"
        # store results
        self._Assay.add_dCt(dCt)


    def _DeltaCt_externally_anchored(self, anchor:float, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using a specified anchor
        """
        # get Ct column label 
        Ct = raw_col_names[1]
        df = self._Assay.get()

        # apply DeltaCt function
        dCt = df[Ct].apply(deltaCt_function, r = anchor, **kwargs)
        dCt.name = "dCt"
        # store results
        self._Assay.add_dCt(dCt)


    def _DeltaCt_grouped_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the first entry of each group as anchor
        """
        # get set of sample groups and dataset
        groups = self._Assay.groups()
        df = self._Assay.get()

        # get Ct column label
        Ct = raw_col_names[1]

        dCt = pd.Series()
        for group in groups: 
            group_subset = df.query(f"group == {group}").reset_index(drop = True)
            anchor = group_subset[Ct][0]
            delta_cts = group_subset[Ct].apply(deltaCt_function, r = anchor, **kwargs)
            dCt = dCt.append(delta_cts).reset_index(drop = True)
        
        dCt.name = "dCt"

        # store results
        self._Assay.add_dCt(dCt)

    def _DeltaCt_first_anchored(self, deltaCt_function, **kwargs):
        """
        Performs DeltaCt using the very first entry of the dataset as anchor
        """
        # get Ct column label
        Ct = raw_col_names[1]
        
        # get first available entry
        # we do this instead of just 0 
        # because the truly first entry 
        # might have been filtered out 
        df = self._Assay.get()
        first = list(df.index)[0]

        # get anchor
        anchor = df[Ct][ first ]

        # apply DeltaCt function
        dCt = df[Ct].apply(deltaCt_function, r = anchor, **kwargs)
        dCt.name = "dCt"
        # store results
        self._Assay.add_dCt(dCt)

    def _exp_DCt(self, s, r, **kwargs):
        """
        Calculates deltaCt exponentially
        """
        factor = s - r 
        return self._Assay._eff **(-factor)

    def _simple_DCt(self, s, r, **kwargs):
        """
        Calculates deltaCt linearly
        """
        return s - r

    def _get_deltaCt_function(self, exp):
        """
        Returns the function to be used for DeltaCt based on 
        whether or not exponential shall be used.
        """
        funcs = {
            "exponential" : self._exp_DCt,
            "linear" : self._simple_DCt
        }
        return funcs.get( exp )
        # if exp == True:
        #     dCt = self._exp_DCt
        # else:
        #     dCt = self._simple_DCt
        # return dCt

    def __str__(self):
        s = f"""
{self.__class__.__name__}: {self._id}\n

Anchor: {self._anchor}\n
Ref.Group: {self._ref_group}\n

DeltaCt: {self._deltaCt_function}
        """.strip()
        return s
    
    def __repr__(self):
        anchor = self._anchor
        ref = self._ref_group
        return f"{self.__class__.__name__}({anchor=}, {ref=})"

__default_Analyser__ = Analyser()
"""The default Analyser"""