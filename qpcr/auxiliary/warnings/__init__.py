"""
This submodule contains the warning strings for the various classes in qpcr __init__.py
It also defines two classes, a SoftWarning (which only prints the warning string) and a HardWarning which raises an Exception.
"""

# this dictionary stores all the warnings for qpcr
WARNINGS = {

"Analyser:newlinked" : "You have linked a new dataset to your Analyser, your previous results have been cleared!\nIf this was not intended, use .get() to get your results before using .link()",
"Analyser:not_newlinked" : "You have precomputed results in your Analyser, no new dataset was linked!\n\If you want to link new data anyway (and clear your current results) use force=True.",
"Analyser:cannot_set_func" : "Only 'exponential'/'linear' or a defined function may be passed as .func()\nReceived func = {func}", 

"Samples:reps_dont_cover" : "The replicates you provide do not cover all data entries!\nExpected: either reps % len(samples) == 0 (if reps is int) or sum(reps) == {n_samples} (if reps is tuple), received reps = {reps}",
"Samples:no_reps_yet" : "No replicates have been defined so far!\nAdd replicates using .replicates()",
"Samples:no_groupname_assignment" : "Could not establish new group names!\nGroup names can only be established based on a list (per index) or dictionary (per key)! New group names received = {names}",
"Samples:groupnames_dont_colver" : "The new group names specified do not cover all current groups of replicates!\nRequired are new_names for all groups in {current_groups} either as dictionary old_name : new_name or as list, replaced per index.", 

"Results:cannot_link" : "Names could not be linked!\nIf names shall be added after results have been computed, use .add_names() instead! .link() is only allowed on empty Results!",
"Results:cannot_add_names" : "Names could not be added! Dimensionality did not match...\nExpected length: {length}, but got {length_names}",
"Results:cannot_add_column" : "The column you try to add to Results does not match the current dimensions! Have: {length}, got: {length_column}",

"Normaliser:cannot_set_prep_func" : "Unknown function supplied for prep_func!\n Received func = {func}", 
"Normaliser:cannot_set_norm_func" : "Unknown function supplied for norm_func!\n Received func = {func}", 
"Normaliser:empty_data" : "Sample {s} was not added because it did not contain any results data!", 
"Normaliser:unknown_data" : "Sample {s} was not added because it could not be read!\nOnly qpcr.Results or qpcr.Analysis objects are allowed!"
}

class SoftWarning(Exception):
    """
    Only prints out the warning string. 
    It takes the warning identifier alongside with any formatted input that the warning may display
    """
    def __init__(self, warning, **kwargs):
        self._message = WARNINGS[warning]
        self._message = self._message.format(**kwargs)
        print("Warning:", self._message, sep = "\n")

class HardWarning(Exception):
    """
    Raises an Exception with the warning string
    """
    def __init__(self, warning, **kwargs):
        self._message = WARNINGS[warning]
        self._message = self._message.format(**kwargs)
        raise Warning(f"\n{self._message}")
