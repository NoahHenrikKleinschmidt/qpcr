"""
This submodule contains the warning strings for the various classes in qpcr __init__.py
It also defines two classes, a SoftWarning (which only prints the warning string) and a HardWarning which raises an Exception.
"""

# this dictionary stores all the warnings for qpcr
WARNINGS = {

"Reader:cannot_read_csv" : "The file \"{file}\" could not be read!\nMake sure the csv file adheres to the two column structure required of csv input files.",
"Reader:cannot_read_multifile" : "The file \"{file}\" appears to specify multiple assays!\nIf you wish to use this file for your input, please specify an assay name through the argument 'assay'.\n\nAvailable assays from this file are: {assays}",

"Analyser:newlinked" : "You have linked a new dataset to your Analyser, your previous results have been cleared!\nIf this was not intended, use .get() to get your results before using .link()",
"Analyser:not_newlinked" : "You have precomputed results in your Analyser, no new dataset was linked!\n\If you want to link new data anyway (and clear your current results) use force=True.",
"Analyser:cannot_set_func" : "Only 'exponential'/'linear' or a defined function may be passed as .func()\nReceived func = {func}", 

"Assay:reps_dont_cover" : "The replicates you provide do not cover all data entries!\nExpected: either reps % len(samples) == 0 (if reps is int) or sum(reps) == {n_samples} (if reps is tuple), received reps = {reps}",
"Assay:reps_could_not_vet" : "The replicates could not be vetted! Replicates must either be integers or tuples of integers. Currently received replicates = {reps}",
"Assay:no_reps_yet" : "No replicates have been defined so far!\nAdd replicates using .replicates()",
"Assay:no_reps_inferred" : "Replicates could be inferred for assay '{assay}'!\nSpecify replicates manually using replicates()",
"Assay:no_groupname_assignment" : "Could not establish new group names!\nGroup names can only be established based on a list (per index) or dictionary (per key)! New group names received = {names}",
"Assay:groupnames_dont_colver" : "The new group names specified do not cover all current groups of replicates!\nRequired are new_names for all groups either as dictionary old_name : new_name or as list, replaced per index.\nCurrent group names \t = {current_groups}\nReceived new group names = {new_received}", 
"Assay:groupnames_not_inferred" : "Could not infer group names based on the provided replicate sample names.\nPlease, specify group names directly.",
"Assay:no_data_adopted" : "No data was adopted by the Assay!\nData is already stored by this Assay. If you wish to overwrite it use the force=True option.",

"Results:cannot_link" : "Names could not be linked!\nIf names shall be added after results have been computed, use .add_names() instead! .link() is only allowed on empty Results!",
"Results:cannot_add_names" : "Names could not be added! Dimensionality did not match...\nExpected length: {length}, but got {length_names}",
"Results:cannot_add_column" : "The column you try to add to Results does not match the current dimensions! Have: {length}, got: {length_column}",
"Results:save_need_dir" : "When saving both df and stats, a directory must be specified to store files in!",

"Normaliser:cannot_set_prep_func" : "Unknown function supplied for prep_func!\n Received func = {func}", 
"Normaliser:cannot_set_norm_func" : "Unknown function supplied for norm_func!\n Received func = {func}", 
"Normaliser:empty_data" : "Assay {s} was not added because it did not contain any results data!", 
"Normaliser:unknown_data" : "Assay {s} was not added because it could not be read!\nOnly qpcr.Results or qpcr.Analysis objects are allowed!",
"Normaliser:norm_unknown_data" : "Normaliser {s} was not added because it could not be read!\nOnly qpcr.Results or qpcr.Analysis objects are allowed!",
"Normaliser:no_data_yet" : "Normalisation cannot be performed, as either samples are missing or no normaliser has been specified yet, or could not be processed!",

"SampleReader:no_reps_yet" : "Could not read data as no replicates have been specified yet!\nPlease, make sure to provide replicate information using the .replicate() method.",

"Plotter:unknown_data" : "Unknown data linkage!\nOnly qpcr.Results or pd.DataFrame objects are allowed!\nReceived: {obj}",
"Plotter:already_compiled" : "Figure is laready compiled! Make sure to use force=True if you wish to re-compile the figure!",
"Plotter:no_fig_yet" : "No figure has been generated yet, nothing to save...",

"Pipeline:no_reps" : "Pipeline requires at least one set of replicates!, Make sure to provide these using .replicates()!",
"Pipeline:no_data" : "Pipeline requires at least one Normaliser and/or Sample Assay!, Make sure to provide these using .add_() or link_() functions!",
"Pipeline:faulty_index" : "The index file does not appear to be properly constructed!\nMake sure to use a standard csv format with the following labeled columns: 'id', 'condition', 'normaliser', and 'path' (plus any additional columns of your own, which will be ignored)",

"Filter:no_assay" : "Filtering requires a sample Assay! First link an Assay object before applying filter()",

"Parser:no_save_loc" : "Cannot save data files because no target directory has been specified yet!\nMake sure to specify a directory using save_to()",
"Parser:no_pattern_yet" : "Assays cannot be identified as no pattern has been specified yet!\nMake sure to specify a pattern using assay_pattern(). Either specify your own pattern, use a pre-defined pattern from the qpcr.Parser.assay_patterns dictionary, or add decorators to your datafile.",
"Parser:no_assays_found" : "No assays could be be identified with the provided pattern!\nMake sure to specify an appropriate pattern for your datafile using assay_pattern(). Either specify your own pattern or use a pre-defined one from the qpcr.Parser.assay_patterns dictionary (patterns from there can be specified by just providing their key to assay_pattern() ).",
"Parser:no_ct_nan_default" : "If allow_nan_ct = False is set then a numeric value must be specified for default_to! Currently default_to = {d}",
"Parser:incompatible_read_kwargs" : "It appears as if some provided kwargs were incompatible with {func}! Defaulting to standard settings for file-reading...\nIf the kwargs you specified are actually important for file reading, try manually reading and parsing to avoid kwarg incompatibilities.",
"Parser:no_data_found" : "No data could be found for at least one assay!\nThis could either be because different headers than the current ones are used above the data column (current header = '{label}'), or because there are too many rows between the header and the data.\nMax allowed rows between assay identifier and data are 2!",
"Parser:invalid_decorator" : "Invalid decorator provided! The decorator '{d}' could not be understood. Available decorators are: {all_d}",
"Parser:no_decorators_found" : "No decorators could be be identified with the provided pattern!\nMake sure to specify the right decorators at the appropriate cells within your datafile.",
"Parser:decorators_but_no_pattern" : "No assay_pattern has been specified yet!\nWill default to just extracting the entire cell content below the decorators. To deal more properly with your assays, please, specify an assay_pattern.",

"MultiReader:empty_data" : "No data is currently stored by the MultiReader!\nIf you already read a file then this could either be because the file did not contain valid decorators, or because it used different headers than the current ones above the data column, or because there are too many rows between the header and the data. Max allowed rows between assay identifier and data are 2! There must not be any rows between decorators and assay identifiers!",
"MultiReader:unknown_datafile" : "Could not read file '{file}'!\nCurrently, only 'csv' and 'excel' files are supported. Make sure to provide either of those formats!",

}

class SoftWarning:
    """
    Only prints out the warning string. 
    It takes the warning identifier alongside with any formatted input that the warning may display
    """
    def __init__(self, warning, **kwargs):
        self._message = WARNINGS[warning]
        self._message = self._message.format(**kwargs)
        print("Warning:", self._message, sep = "\n")

class HardWarning:
    """
    Raises an Exception with the warning string
    """
    def __init__(self, warning, **kwargs):
        self._message = WARNINGS[warning]
        self._message = self._message.format(**kwargs)
        raise Warning(f"\n{self._message}")
