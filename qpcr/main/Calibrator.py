"""
This is the `qpcr.Calibrator` class that is able to compute qPCR amplification efficiencies from `qpcr.Assay` objects.


Amplification Efficiencies in ``qpcr``
------------------------

By default ``qpcr`` sets the amplification efficiency of a new ``qpcr.Assay`` to ``1`` (100%).  However, they can be set to any percentage (also > 1) using the ``efficiency`` method of the ``qpcr.Assay``.
``qpcr`` stores the efficiency as percentage but actually calculates with :math:`2 \ \cdot \ efficiency` when computing Delta-Ct values.

Assigning an assay's efficiency
------

The ``qpcr.Calibrator`` is dedicated to either computing the amplification efficiency of an assay or assigning existing effiencies that have been calculated elsewhere.
In order to compute new efficiencies the Calibrator requires a set of *decorated replicates* that come from a **dilution series**. You can check out `this tutorial <https://github.com/NoahHenrikKleinschmidt/qpcr/blob/main/Examples/9_custom_efficiencies.ipynb>`_ to learn more.

If appropriate data is available we can use the ``qpcr.Calibrator`` to compute and assign new efficiencies by:

.. code-block:: python

    calibrator = qpcr.Calibrator()

    assay = calibrator.calibrate( assay )

This will read the data, perform a linear regression to determine the efficiency, and assign the computed efficiency to the assay. Also, the efficiency is now stored by the Calibrator.
The stored efficiencies can be easily saved to a file using:

.. code-block:: python

    calibrator.save( "my_efficiencies.csv" )
    
We usually do not have a dilution series in each of our datasets, however. So most often you will wish to assign efficiencies that have been already computed from other qPCR runs.
For these cases, the Calibrator can read a reference "database file" and assign existing efficiencies to new Assays as long as their ``id`` is present among the reference efficiencies.

.. code-block:: python

    # read already existing efficiencies from a file
    calibrator.load( "my_efficiencies.csv" )

    # and now simple assign an existing efficiency to the assay
    assays = calibrator.assign( assay )

If we have both assays with existing efficiencies and such with new dilution series data, we can actually just use the Calibrator's ``pipe`` method to process all of them at once. 

.. code-block:: python

    many_assays = [ ... ]

    calibrator.load( "my_efficiencies.csv" )

    # pipe all assays, which will assign where possible, and calibrate anew where necessary (and data is available)
    many_assays = calibrator.pipe( many_assays )

    # finally, save the (now updated) database of efficiencies 
    # this will by default update the file "my_efficiencies.csv" which was already loaded.
    calibrator.save()

"""

from qpcr import defaults
from qpcr import _auxiliary as aux
from qpcr._auxiliary import warnings as aw

import pandas as pd
import numpy as np
import scipy.stats as scistats


from qpcr.Curves import EfficiencyCurve
from qpcr.main import Assay

logger = aux.default_logger()


class Calibrator(aux._ID):
    """
    Calculates qPCR primer efficiency based on a dilution series.
    The dilution series may either be represented as an entire assay
    or as a subset of groups within an assay denoted as `calibrator : {some_name}`.
    In this mode, calibrator replicates will be removed after calibration is done.

    It is possible to specify the dilution steps directly in the groupnames as:
    `calibrator: {some_name}: dil` where `dil` is the inverse dilution step, e.g.
    `calibrator: my_sample: 2` for a `1 : 2` dilution or `calibrator: my_sample: 100`
    for a `1 : 100`. Note, this will have to be present in **each** groupname!

    Alternatively, if no dilution is specified in the groupnames or they cannot be inferred
    for some other reason, it is possible to supply a dilution step via
    the `qpcr.Calibrator.dilution` method.
    """

    __slots__ = ["_id", "_loaded_file", "_eff_dict", "_computed_values", "_orig_dilution", "_manual_dilution_set"]

    def __init__(self):
        super().__init__()
        self._eff_dict = {}  # stores the efficiencies as assay_id : efficiency
        self._computed_values = {}  # stores newly computed efficiency data as:
        # assay_id :  EfficiencyCurve(...)
        #             The EfficiencyCurve object stores:
        #             - dilutions
        #             - Ct values
        #             - the Linreg model
        self._dilution_step = None  # the dilution step(s) used
        self._manual_dilution_set = False
        self._orig_dilution = None

        self._loaded_file = None

    def __str__(self):
        effs = str(pd.DataFrame(self._eff_dict))
        _length = len(effs.split("\n")[0])
        s = f"""
{"-" * _length}
{self.__class__.__name__}:\t{self._id}
Loaded File:\t{self._loaded_file}
        """.strip()
        if self._manual_dilution_set:
            s = f"{s}\nDilution:\t{self._orig_dilution}"
        s = f""""{s}\n{"-" * _length}\n{effs}\n{"-" * _length}"""
        return s

    def __repr__(self):
        file = self._loaded_file
        effs = self._eff_dict
        return f"{self.__class__.__name__}({file=}, {effs=})"

    def save(self, filename: str = None, mode: str = "write"):
        """
        Saves the calculated efficiencies to a `csv` file.

        Parameters
        -------
        filename : str
            The filepath in which to store the efficiencies. If a file
            was already loaded then by default the same file will be
            used to save values again.

        mode : str
            Can be either `"write"` to fully overwrite an existing file,
            with the newly computed data, or `"append"` to only add newly
            computed efficiencies.
        """

        if filename is None and self._loaded_file is not None:
            filename = self._loaded_file

        if mode == "write":
            self._save(filename, self._eff_dict)

        elif mode == "append":
            current = self._load(filename)
            new = {**current, **self._eff_dict}
            self._save(filename, new)

        else:
            e = aw.CalibratorError("unknown_savemode", mode=mode)
            logger.info(e)

    def load(self, filename, merge: bool = True, supersede: bool = False):
        """
        Loads a `csv` file of previously computed efficiencies.

        Parameters
        -------
        filename : str
            The filepath to load efficiencies from.
        merge : bool
            In case efficiencies are already loaded, merge the
            new and existing ones. If `False` the current ones will
            be replaced completely.
        supersede : bool
            In case efficiencies of the same assay are already loaded
            they will be overwritten by the newly incoming ones if `supersede = True`.
        """
        try:
            current = self._load(filename)
            if merge:
                current = {**self._eff_dict, **current} if supersede else {**current, **self._eff_dict}
            self.adopt(current)
            self._loaded_file = filename
            return current
        except Exception as e:
            logger.error(e)
            e = aw.CalibratorError("unknown_filetype", filename=filename)
            logger.critical(e)
            raise e

    def get(self, which="efficiencies"):
        """
        Returns
        -------
        dict
            Either the stored efficiencies (if
            `which = "efficiencies"`) or
            the computed values of newly computed
            efficiencies (if `which = "values"`).
        """
        if which == "efficiencies":
            return self.efficiencies
        elif which == "values":
            return self.computed_values

    def merge(self, *filenames, outfile=None, adopt=True):
        """
        Merges multiple efficiency files together into a single one.

        Parameters
        -------
        filenames : iterable
            Filepaths to load data from which should be merged together.

        outfile : str
            The filepath in which to store the merged efficiencies.
            Not saved if set to `None`.

        adopt : bool
            Will adopt the merged dictionary as its own if `True` (default).

        Returns
        -------
        all_effiencies : dict
            The merged dictionary of all efficiencies from all files.
        """
        all_efficiencies = {}
        for filename in filenames:
            new = self._load(filename)
            all_efficiencies = {**all_efficiencies, **new}

        if outfile is not None:
            self._save(outfile, all_efficiencies)

        if adopt:
            self.adopt(all_efficiencies)

        return all_efficiencies

    def reset(self):
        """
        Resets the Calibrator to initial settings. This will
        clear all stored efficiency values and computed data!
        """
        self.__init__()
        return self

    def clear(self):
        """
        Will clear all stored efficiency values and computed data.
        """
        self._eff_dict = {}
        self._computed_values = {}
        return self

    def adopt(self, effs: dict):
        """
        Adopts an externally generated dictionary of `assay : efficiency`
        structure as its own.

        Parameters
        -------
        effs : dict
            A dictionary where keys are Assay Ids (`str`)
            and values are `float` efficiencies.
        """
        if aux.same_type(effs, {}):
            self._eff_dict = effs
        else:
            aw.SoftWarning("Calibrator:cannot_adopt", effs=effs, eff_type=type(effs).__name__)
        return self

    def dilution(self, step: float or np.ndarray or tuple = None):
        """
        Gets or sets the dilution steps used. This must be a `float` fraction
        e.g. `0.5` for a `1 : 2` dilution series or `0.1` for a `1 : 10` series etc.
        If there are multiple steps because there is a gap in the dilution series. It is
        necessary to supply a step for each group individually e.g. `[1,0.5,0.25,0.0625,0.03125]`.
        if there are 5 dilution steps (originally six but 0.125 was discarded).

        Note, both of the above also work with the inverse dilutions e.g. `2` or `[1,2,4,16,32]`.

        By default the `qpcr.Calibrator` tries to infer the dilutions automatically.
        This only works, however, if the calibrator groupnames specify `calibrator: {some name} : dil` where
        `dil` is the inverse dilution step (e.g. `calibrator: my_sample: 2` for a `1 : 2` dilution). Note,
        it is important that the dilution step is given as the inverse (i.e. *not* as `1:2 or 1/2` or something else! )

        Parameters
        ----------
        step : float or np.ndarray
            The dilution step used.

        Returns
        -------
        dilution : float or np.ndarray
            The currently used dilution step.
        """

        dilution = self._dilution(step)
        self._orig_dilution = self._dilution_step

        # if the dilution is now set to a valid
        # value we set the _manual_dilution_set check to True
        if dilution is not None:
            self._manual_dilution_set = True
        return dilution

    def pipe(self, assay: Assay, remove_calibrators: bool = True, ignore_uncalibrated: bool = False):
        """
        A wrapper for calibrate / assign.

        This method will first try to assign pre-computed efficiencies
        and if no matching ones are found it will try to calculate a new efficiency
        from the assay.

        Parameters
        ----------
        assay : qpcr.Assay
            A `qpcr.Assay` object.
        remove_calibrators : bool
            If calibrators are present in the assay alongside other groups,
            remove the calibrator replicates after assignment or efficiency calculation.
        ignore_uncalibrated : bool
            If `True` assays that could neither be newly calibrated nor be assigned an existing
            efficiency will be ignored. Otherwise, and error will be raised.

        Returns
        -------
        assay : qpcr.Assay
            The now calibrated `qpcr.Assay`.
        """
        if isinstance(assay, list):
            return [self.pipe(A, remove_calibrators=remove_calibrators, ignore_uncalibrated=ignore_uncalibrated) for A in assay]
        if self._eff_dict != {}:
            # first try to assign (will leave the assay unchanged if nothing is found)
            eff = self._get_efficiency(assay)
            if eff is not None:
                assay = self.assign(assay, remove_calibrators=remove_calibrators)
            else:
                try:
                    assay = self.calibrate(assay, remove_calibrators=remove_calibrators)
                except Exception as e:
                    logger.info(e)
                    e = aw.CalibratorError("cannot_process_assay", id=assay.id())
                    if not ignore_uncalibrated:
                        logger.critical(e)
                        raise e
                    else:
                        logger.info(e)
        else:
            try:
                assay = self.calibrate(assay, remove_calibrators=remove_calibrators)
            except Exception as e:
                logger.info(e)
                e = aw.CalibratorError("cannot_process_assay", id=assay.id())
                if not ignore_uncalibrated:
                    logger.critical(e)
                    raise e
                else:
                    logger.info(e)
        return assay

    def calibrate(self, assay: Assay, remove_calibrators: bool = True):
        """
        Computes an efficiency from an `qpcr.Assay` object.

        This method will try to compute a new efficiency. To do this, it will check autonomously if
        `calibrator : {}` replicates are present and use these for computation. If none are
        found it will assume the entire assay is to be used as calibrator.

        Note
        ----
        Calibrators are searched for through the group `names` not the replicate ids!

        Parameters
        ----------
        assay : qpcr.Assay
            A `qpcr.Assay` object.

        remove_calibrators : bool
            If calibrators are present in the assay alongside other groups,
            remove the calibrator replicates after efficiency calculation.
        """

        if isinstance(assay, list):
            return [self.calibrate(assay=i, remove_calibrators=remove_calibrators) for i in assay]

        # get the assay's groupnames and check for the calibrator prefix.
        names = assay.names(as_set=False).unique()

        # check if any groups are declared as calibrators
        calibrators = np.array([self._has_calibrator_prefix(i) for i in names])
        has_calibrators = any(calibrators)

        # now get the relevant dataframe for the computation
        # this will either be the entire df (if no calibrator groups are present)
        # or just the subset of calibrators
        df = assay.get()

        if has_calibrators:
            df = self._subset_calibrators(names, calibrators, df)

        # get Ct column name
        ct_name = defaults.raw_col_names[1]
        # drop NaN cols as they are incompatible with linregress anyway...
        df = df[df[ct_name] == df[ct_name]]

        # now sort the dataframe by Ct values as they need to strictly
        # increase for dilution series.
        df = df.sort_values(ct_name).reset_index()
        df = df.rename(columns={"index": "orig_index"})

        # now generate dilution steps ( i.e. "concentrations" )
        # to do that we first need to check if dilutions have not
        # been supplied, and then try to infer them based on the groupnames
        # if we got an input for dilution() we use  that to generate a
        # dilution steps array...

        # NOTE: The non-log-scaled dilutions are now stored in self._dilution_steps
        #       while the log-scaled versions are returned. Hence, the dilutions
        #       variable below is the log-scaled version!
        if not self._manual_dilution_set:
            dilutions = self._infer_dilution_steps(df)
        else:
            dilutions = self._generate_dilution_steps(df)

        # now interpolate a line through the log dilutions and the ct values
        cts = df[ct_name].to_numpy()

        regression_line = scistats.linregress(x=dilutions, y=cts)
        # and now compute the efficiency from the regression line
        efficiency = self._compute_efficiency(regression_line)

        # and now assign the efficiency to the assay
        assay.efficiency(efficiency)

        # save the efficiency in self._eff_dict
        # self._eff_dict.update( { assay.id() : assay.efficiency() } )
        self._eff_dict[assay.id()] = assay.efficiency()

        # and, finally, save the computed values and the efficiency
        self._save_computation(assay, dilutions, cts, regression_line)

        # now remove the calibrators from the assay
        # but only do so in case there were calibrators found!
        # If the entire assay was used, we do not delete the entries...
        if has_calibrators and remove_calibrators:
            self._remove_calibrators(assay, df)

        return assay

    def assign(self, assay: Assay, remove_calibrators: bool = True):
        """
        Assigns an efficiency to an `qpcr.Assay` based on its Id.
        This requires that an efficiency corresponding to the Assay's Id
        is present in the currently loaded / computed effiencies.

        Parameters
        ----------
        assay : qpcr.Assay
            A `qpcr.Assay` object.

        remove_calibrators : bool
            If calibrators are present in the assay alongside other groups,
            remove the calibrator replicates.
        """
        eff = self._get_efficiency(assay)
        if eff is not None:
            # set assay's efficiency
            assay.efficiency(eff)

            # check if any groups are declared as calibrators
            # that should be removed from the df
            names = assay.names(as_set=False).unique()
            calibrators = np.array([self._has_calibrator_prefix(i) for i in names])
            has_calibrators = any(calibrators)
            if has_calibrators and remove_calibrators:
                df = assay.get()
                df = self._subset_calibrators(names, calibrators, df)
                df = df.sort_values(defaults.raw_col_names[1]).reset_index()
                df = df.rename(columns={"index": "orig_index"})
                self._remove_calibrators(assay, df)
        else:
            aw.SoftWarning("Calibrator:could_not_assign", id=assay.id())
        return assay

    def plot(self, mode: str = None, **kwargs):
        """
        A shortcut to call a `qpcr.Plotters.EfficiencyLines` plotter
        to visualise the regression lines from de novo efficiency computations.

        Parameters
        -------
        mode : str
            The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).
        **kwargs
            Any additional keyword arguments to be passed to the plotter.

        Returns
        -------
        fig : plt.figure or plotly.figure
            The figure generated by `EfficiencyLines`.
        """
        import qpcr.Plotters as Plotters

        plotter = Plotters.EfficiencyLines(mode=mode)
        plotter.link(self)
        fig = plotter.plot(**kwargs)
        return fig

    def __qplot__(self, **kwargs):
        return self.plot

    @property
    def efficiencies(self):
        """
        Returns
        ------
        dict
            The currently stored efficienies.
        """
        return self._eff_dict

    @property
    def computed_values(self):
        """
        Returns
        ------
        dict
            The currently stored values from newly
            computed efficiencies.
        """
        return self._computed_values

    def _dilution(self, step=None):
        """
        The functional core of self.dilution() the only difference is
        that self.dilution also sets a boolean attribute self._manual_dilution_set to True...
        Which signals that the manually supplied dilutions should be used rather than that
        they should be inferred from the dataset...
        """
        if step is not None:
            unknown_datatype = not isinstance(step, (float, int, np.ndarray, pd.Series, tuple, list))
            if unknown_datatype:
                aw.HardWarning("Calibration:cannot_interpret_dilution", step=step, step_type=type(step).__name__)

            # check if we need to convert to numpy array
            # because we have an iterable without math operations support.
            if isinstance(step, (tuple, list)):
                step = np.array(step)

            # check for an ndarray and make sure to invert if
            # the dilution steps are given as 2 4 instead of
            # 0.5 0.25 etc., also do the same for a single number...

            is_inverse_array = isinstance(step, (np.ndarray, pd.Series)) and any(step > 1)
            is_inverse_number = isinstance(step, (float, int)) and step > 1
            need_inverse = is_inverse_array or is_inverse_number

            if need_inverse:
                step = 1 / step

            # and store new steps
            self._dilution_step = step

        return self._dilution_step

    def _remove_calibrators(self, assay, df):
        """
        Drops all calibrator replicates from the dataframe of the Assay.

        Note
        -------
        This leaves the index unchanged! Possible we might wish to also reset the
        index during this step...
        """
        to_remove = df["orig_index"].to_numpy()
        index = np.zeros(assay.n())
        for i in to_remove:
            index[i] = 1
        index = np.argwhere(index == 1)
        index = np.squeeze(index)
        assay.ignore(index, drop=True)

    def _save_computation(self, assay, dilutions, ct_values, linreg):
        """
        Creates a new entry in self._computed_values for the newly computed
        efficiency.
        """
        self._computed_values[assay.id()] = EfficiencyCurve(dilutions=dilutions, ct_values=ct_values, model=linreg, efficiency=assay.efficiency())
        # also add the Id of the Assay to the _EfficiencyCurve
        self._computed_values[assay.id()].id(assay.id())

    def _compute_efficiency(self, regression_line):
        """
        Calculates the efficiency from the regression line slope.
        """
        slope = regression_line.slope
        efficiency = -1 / slope
        efficiency = np.exp(efficiency)
        efficiency -= 1
        efficiency = round(efficiency, 4)
        return efficiency

    def _subset_calibrators(self, names, calibrators, df):
        """
        Generates a dataframe subset containing only calibrator replicates.
        """
        # generate a total query formula for all found calibrators
        q = names[calibrators]
        q = "' or group_name == '".join(q)
        q = "group_name == '" + q + "'"
        df = df.query(q)
        return df

    def _infer_dilution_steps(self, df):
        """
        Infers the dilution steps from the group names if they are specified
        as `calibrator: some_name: dilution` e.g. `calibrator: mysample: 5`.
        """
        try:
            # get dilution steps from groupnames in format calibrator: name : dil
            steps = df["group_name"].apply(lambda x: float(x.split(":")[2]))
            steps = steps.to_numpy()

            # preprocess to get proper format
            self._dilution(steps)

            # get and transform to log-scale
            dilutions = self._dilution()
            dilutions = np.log(dilutions)
            self._reset_dilution()
            return dilutions
        except Exception as e:
            logger.error(e)
            e = aw.CalibratorError("could_not_infer_dilution")
            logger.critical(e)
            raise e

    def _generate_dilution_steps(self, df):
        """
        Generates a numpy ndarray of log-scaled dilution steps
        for the calibrators.
        """
        # generate steps range
        # we shall use the concept of (dilution)^m to generate
        # the dilution steps. To get m we use a re-anchored df["group"]
        # column

        self._reset_groups(df)
        steps = df["group"]
        counts = df["group"].value_counts(sort=False)
        repeats = steps.size if isinstance(self._dilution(), float) else counts

        # repeat the dilution steps to match the group replicate numbers
        dilutions = np.repeat(self._dilution(), repeats)

        # if we only had a single value for the dilution we
        # also need to now scale the dilutions to generate the
        # actual dilution scale. If we already got a dilution series
        # as input, we must not do this as we otherwise doubly scale...
        if isinstance(self._dilution(), float):
            dilutions = dilutions**steps

        # save dilutions
        self._dilution(dilutions)

        # and transform to log scale
        dilutions = np.log(dilutions)

        # and reset the diltion back to the what was
        # originally set (or None from init)
        self._reset_dilution()
        return dilutions

    def _reset_dilution(self):
        """
        Resets the to the default dilution step to ensure
        the same starting conditions are met for each assay
        as it is passed through the Calibrator.
        """
        self._dilution_step = self._orig_dilution

    def _reset_groups(self, df):
        """
        Resets the numeric group identifiers to start continuously from 0.
        This method sets the initial group to 0 and then successively resets any
        gaps to match e.g. a 0,1,3 -> 0,1,2...
        """

        # get counts of each group_name in the df
        counts = df["group_name"].value_counts(sort=False)

        # generate new numeric identifiers for each group
        new_groups = np.arange(len(counts))

        # transform to match the right repeats
        new_groups = np.repeat(new_groups, counts)

        # and set new groups
        df["group"] = new_groups

    def _has_calibrator_prefix(self, string):
        """
        Checks if the calibrator prefix is the start of a string
        The string is a replicate group name in this case...
        """
        calibrator_prefix = defaults.calibrator_prefix
        return string.startswith(calibrator_prefix)

    def _get_efficiency(self, assay: Assay):
        """
        Returns the efficiency from the currently loaded effiencies
        that corresponds to the assay's Id. Returns None if no match is present.
        """
        id = assay.id()
        effs = self.get()

        if id not in effs.keys():
            return None
        else:
            return effs[id]

    def _load(self, filename):
        """
        Loads a csv file but does not adobt the data as its own yet.
        Returns a dictionary.
        """
        # current = json.load( open(filename, "r") )
        current = pd.read_csv(filename)
        current = current.to_dict("split")["data"]
        current = {id: eff for id, eff in current}
        return current

    def _save(self, filename, dict_to_save):
        """
        Saves a dictionary to a csv file.
        """
        # json.dump( dict_to_save , open(filename, "w") )
        df = pd.DataFrame(dict_to_save, index=["eff"])
        df = df.transpose().reset_index()
        df.to_csv(filename, index=False)


__default_Calibrator__ = Calibrator()
"""The default Calibrator"""
