"""
This submodule contains the warning strings for the various classes in qpcr __init__.py
It also defines two classes, a SoftWarning (which only prints the warning string) and a HardWarning which raises an Exception.
"""

# this dictionary stores all the warnings for qpcr
WARNINGS = {

"Reader:cannot_read_csv" : "The file \"{file}\" could not be read!\nMake sure the csv file adheres to the two column structure required of csv input files.",
"Reader:cannot_read_multifile" : "The file \"{file}\" appears to specify multiple assays!\nIf you wish to use this file for your input, please specify an assay name through the argument 'assay'.\n\nAvailable assays from this file are: {assays}\n\nIf you wish to read all assays in this file use the qpcr.MultiReader instead!",
"Reader:cannot_find_datacols" : "No data could be identified using the given column names!\nCurrently provided are id_label = '{id_label}' and ct_label = '{ct_label}'. Make sure to provide correct data column headers.",

"Analyser:newlinked" : "You have linked a new dataset to your Analyser, your previous results have been cleared!\nIf this was not intended, use .get() to get your results before using .link()",
"Analyser:not_newlinked" : "You have precomputed results in your Analyser, no new dataset was linked!\n\If you want to link new data anyway (and clear your current results) use force=True.",
"Analyser:cannot_set_func" : "Only 'exponential'/'linear' or a defined function may be passed as .func()\nReceived func = {func}", 

"Assay:reps_dont_cover" : "The replicates you provide do not cover all data entries!\nExpected: either reps % len(samples) == 0 (if reps is int) or sum(reps) == {n_samples} (if reps is tuple), received reps = {reps}",
"Assay:reps_could_not_vet" : "The replicates could not be vetted! Replicates must either be integers or tuples of integers. Currently received replicates = {reps}",
"Assay:no_reps_yet" : "No replicates have been defined so far!\nAdd replicates using .replicates()",
"Assay:no_reps_inferred" : "Replicates could not be inferred for assay '{assay}'!\nSpecify replicates manually using replicates()",
"Assay:no_groupname_assignment" : "Could not establish new group names!\nGroup names can only be established based on a list (per index) or dictionary (per key)! New group names received = {names}",
"Assay:groupnames_dont_colver" : "The new group names specified do not cover all current groups of replicates!\nRequired are new_names for all groups either as dictionary old_name : new_name or as list, replaced per index.\nCurrent group names \t = {current_groups}\nReceived new group names = {new_received}", 
"Assay:groupnames_not_inferred" : "Could not infer group names based on the provided replicate sample names.\nPlease, specify group names directly.",
"Assay:no_data_adopted" : "No data was adopted by the Assay!\nData is already stored by this Assay. If you wish to overwrite it use the force=True option.",
"Assay:setup_not_grouped" : "Could not group the dataset! This is probably because the replicates could either not be inferred or the provided replicates do not cover all the data. Make sure to properly specify replicates.",

"Results:cannot_link" : "Names could not be linked!\nIf names shall be added after results have been computed, use .add_names() instead! .link() is only allowed on empty Results!",
"Results:cannot_add_names" : "Names could not be added! Dimensionality did not match...\nExpected length: {length}, but got {length_names}",
"Results:cannot_add_column" : "The column you try to add to Results does not match the current dimensions! Have: {length}, got: {length_column}",
"Results:save_need_dir" : "When saving both df and stats, a directory must be specified to store files in!",
"Results:name_overlap" : "Could not add results from computation '{name}' because it appears results from this assay against these normalisers are already stored!\n",

"Normaliser:cannot_set_prep_func" : "Unknown function supplied for prep_func!\n Received func = {func}", 
"Normaliser:cannot_set_norm_func" : "Unknown function supplied for norm_func!\n Received func = {func}", 
"Normaliser:empty_data" : "Assay {s} was not added because it did not contain any results data!", 
"Normaliser:unknown_data" : "Assay {s} was not added because it could not be read!\nOnly qpcr.Assay or qpcr.Analyser objects are allowed!",
"Normaliser:norm_unknown_data" : "Normaliser {s} was not added because it could not be read!\nOnly qpcr.Assay or qpcr.Analyser objects are allowed!",
"Normaliser:no_data_yet" : "Normalisation cannot be performed, as either samples are missing or no normaliser has been specified yet, or could not be processed!",

"SampleReader:no_reps_yet" : "Could not read data as no replicates have been specified yet!\nPlease, make sure to provide replicate information using the .replicate() method.",

"Plotter:unknown_data" : "Unknown data linkage!\nOnly qpcr.Results or pd.DataFrame objects are allowed!\nReceived: {obj}",
"Plotter:already_compiled" : "Figure is laready compiled! Make sure to use force=True if you wish to re-compile the figure!",
"Plotter:no_fig_yet" : "No figure has been generated yet, nothing to save...",

"Pipeline:no_reps" : "Pipeline requires at least one set of replicates!, Make sure to provide these using .replicates()!",
"Pipeline:no_data" : "Pipeline requires at least one Normaliser and Sample Assay!, Make sure to provide these using .add_() or link() methods!",
"Pipeline:faulty_index" : "The index file does not appear to be properly constructed!\nMake sure to use a standard csv format with the following labeled columns: 'id', 'condition', 'normaliser', and 'path' (plus any additional columns of your own, which will be ignored)",
"Pipeline:no_data_input" : "Could not interpret data input: '{file}'!\nIt appears the given data input is neither an existing file nor a directory containing files. Make sure to provide valid data inputs.",


"Filter:no_assay" : "Filtering requires a sample Assay! First link an Assay object before applying filter()",

"Parser:no_save_loc" : "Cannot save data files because no target directory has been specified yet!\nMake sure to specify a directory using save_to()",
"Parser:no_pattern_yet" : "Assays cannot be identified as no pattern has been specified yet!\nMake sure to specify a pattern using assay_pattern(). Either specify your own pattern, use a pre-defined pattern from the qpcr.Parser.assay_patterns dictionary, or add decorators to your datafile.",
"Parser:no_assays_found" : "No assays could be be identified with the provided pattern!\nMake sure to specify an appropriate pattern for your datafile using assay_pattern(). Either specify your own pattern or use a pre-defined one from the qpcr.Parser.assay_patterns dictionary (patterns from there can be specified by just providing their key to assay_pattern()).\nIf you are sure you provided the correct pattern or used decorators, maybe check if you provide the correct input 'col' argument to search in (by default the first column / row is used).\nIn case of multi-sheet excel files also make sure to provide the correct 'sheet_name' (if you are not using the MultiSheetReader).",
"Parser:no_ct_nan_default" : "If allow_nan_ct = False is set then a numeric value must be specified for default_to! Currently default_to = {d}",
"Parser:incompatible_read_kwargs" : "It appears as if some provided kwargs were incompatible with {func}! Defaulting to standard settings for file-reading...\nIf the kwargs you specified are actually important for file reading, try manually reading and parsing to avoid kwarg incompatibilities.",
"Parser:no_data_found" : "No data could be found for at least one assay!\nThis could either be because different headers than the current ones are used above the data column (current header = '{label}'), or because there are too many rows between the header and the data.\nMax allowed rows between assay identifier and data are 2!",
"Parser:invalid_decorator" : "Invalid decorator provided! The decorator '{d}' could not be understood. Available decorators are: {all_d}",
"Parser:no_decorators_found" : "No decorators could be be identified with the provided pattern!\nMake sure to specify the right decorators at the appropriate cells within your datafile.",
"Parser:decorators_but_no_pattern" : "No assay_pattern has been specified yet!\nWill default to just extracting the entire cell content below the decorators. To deal more properly with your assays, please, specify an assay_pattern.",
"Parser:invalid_range" : "No data range could be generated from provided inputs!",
"Parser:found_non_readable_cts" : "Assay: '{assay}'\nAt least one Ct value of this assay could not be read and was set to NaN!\nThe value responsible for this warning was {bad_value}",
"Parser:bigtable_no_replicates" : "No replicates have been specified!\nCannot infer a horizontal Big Table without replicate information!",
"Parser:no_bigtable_header" : "The provided Big Table anchor: '{header}' could not be identified within the data!\nMake sure to provide a valid column header as id_label using labels().",

"MultiReader:empty_data" : "No data is currently stored by the MultiReader!\nIf you already read a file then this could either be because the file did not contain valid decorators, or because it used different headers than the current ones above the data column, or because there are too many rows between the header and the data. Max allowed rows between assay identifier and data are 2! There must not be any rows between decorators and assay identifiers!",
"MultiReader:unknown_datafile" : "Could not read file '{file}'!\nCurrently, only 'csv' and 'excel' files are supported. Make sure to provide either of those formats!",
"MultiReader:no_decorator_or_pattern" : "No assays could be identified!\nMake sure to either specify a valid assay_pattern or decorate your assays in your file.",

"MultiSheetReader:sheet_unreadable" : "Sheet: {sheet}\n{e}",

"BigTableReader:no_cols" : "In order to parse a vertical Big Table both a 'ct_col' and 'assay_col' have to be specified.\nMake sure to specify valid columns to your data for these inputs! Currently ct_col = {ct_col} and assay_col = {assay_col}",
"BigTableReader:cols_no_good" : "In order to parse a vertical Big Table both a 'ct_col' and 'assay_col' have to be specified.\nHowever, the currently specified columns cannot be found! Currently ct_col = {ct_col} and assay_col = {assay_col}",

"Versions:Deprecation" : "Class {old} is deprecated and will be dropped in a future release! Please, use {new} instead.",

}

class SoftWarning:
    """
    Only prints out the warning string. 
    It takes the warning identifier alongside with any formatted input that the warning may display
    """
    def __init__(self, warning, **kwargs):
        self._message = WARNINGS[warning]
        self._was_raised = False
        self._once_only = kwargs.pop("once", False)
        if not self._once_only:
            self.trigger(**kwargs)
            
    def trigger(self, **kwargs):
        should_not_raise = self._once_only and self._was_raised
        if not should_not_raise:
            self._message = self._message.format(**kwargs)
            self._was_raised = True
            print("\nWarning:", self._message, sep = "\n")

class HardWarning:
    """
    Raises an Exception with the warning string
    """
    def __init__(self, warning, **kwargs):
        self._message = WARNINGS[warning]
        traceback = kwargs.pop("traceback", True)
        self._message = self._message.format(**kwargs)
        if traceback:
            raise Exception(f"\n{self._message}")
        else:
            print("\nError:", self._message, sep = "\n")
            exit(1)
