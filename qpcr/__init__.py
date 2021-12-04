"""
This module is designed to provide functions to analyse qPCR data. 
It is designed for maximal user-friendliness and streamlined data-visualisation.
"""
import statistics as stat 
import pandas as pd
import qpcr.auxiliary as aux
import glob, os


# TODO: A class to read csv pqcr raw data
# TODO: A class to perform delta delta ct 
# TODO: A class to handle and store results
# TODO: A class to visulalise results
# TODO: A compilation of common pipelines (like the qpcr.Analysis module from earlier)

RAW_COL_NAMES = ["Sample", "Ct"]

class Reader:
    """
    This class reads qpcr raw data files in csv format. 
    It requires that two columns are present, one for sample names, one for Ct values
    it will load these into a pandas dataframe.
    """
    def __init__(self, filename:str) -> pd.DataFrame: 
        self._src = filename
        self._delimiter = ";" if self._is_csv2() else ","
        self.read()

    def get(self):
        """
        Returns the self._df dataframe
        """
        return self._df

    def read(self):
        """
        Reads the data file
        """
        self._df = pd.read_csv(
                                self._src, 
                                sep = self._delimiter, 
                                header = self._has_header(), 
                                names = RAW_COL_NAMES
                            )

    def _is_csv2(self):
        """
        Tests if csv file is ; delimited (True) or common , (False)
        """
        with open(self._src, "r") as openfile: 
            content = openfile.read()
        if ";" in content: 
            return True
        return False

    def _has_header(self):
        """
        Checks if column headers are provided in the data file
        It does so by checking if the second element in the first row is numeric
        if it is numeric (returns False) no headers are presumed. 
        """
        with open(self._src, "r") as openfile: 
            content = openfile.read().split("\n")[0]
        try: 
            second_col = content[1]
            float(second_col)
            return False
        except TypeError:
            return True


if __name__ == "__main__":
    a = Reader("Example Data/28S.csv")
    print(a.get())

    exit(0)

    
def open_csv_file(filename, export="dict"):
    """
    This function opens a given input csv file. To function properly, the file is required to have 
    a first column of sample names and a second column with Ct values. 
    The first line will be read as Column Titles and be cropped during processing.
    """

    with open(filename, "r") as openfile:
        contents = openfile.read()

    open_args = dict(header = 0, names = ["Sample", "Mean Ct"])
    
    if ";" in contents:
        contents = pd.read_csv(filename, sep = ";", **open_args )
    else: 
        contents = pd.read_csv(filename, **open_args )
    
    if export == "dict":
        contents = contents.to_dict(orient = "list")

    return contents

def _define_replicates(sample_dict, column_name, replicates):
    """
    Auxiliary Function that creates groups of indices to be used by group_samples for replicate assignment.
    """
    reps, equal_reps = _prep_replicates(replicates)
    
    samples = sample_dict[column_name]
    intervals = _make_intervals(reps, equal_reps, samples)
    if intervals is None:
        return None

    replicate_samples = []
    if column_name == "Mean Ct": # the numeric column
        for i in intervals:
            tmp = [float(samples[j]) for j in i]
            replicate_samples.append(tmp)
    else:
         for i in intervals:
            tmp = [str(samples[j]) for j in i]
            replicate_samples.append(tmp)
    
    return replicate_samples

def _make_intervals(reps, equal_reps, samples):
    """
    Makes equal or unequal intervals for samples according to reps and 
    checks if intervals cover all provided sample data entries.
    Returns None if not all samples are covered...
    """
    if equal_reps == True:
        intervals = _equal_intervals(samples, reps)
    else:
        intervals = _unequal_intervals(samples, reps)
    return intervals

def _unequal_intervals(samples, reps):
    if sum(reps) != len(samples): 
        print(f"You are not providing replicate annotation for each sample!\n{len(samples)-sum(reps)} samples are missing intervals!")
        return None
    intervals = _create_unequal_intervals(reps)
    return intervals

def _equal_intervals(samples, reps):
    """
    Makes equal intervals according to reps specified, 
    and checks it they are appropriate for the given sample. 
    Will return none if not...
    """
    intervals = _create_equal_intervals(len(samples), reps)
    if _vet_equal_intervals(intervals, samples) == False:
        return None
    return intervals

def _vet_equal_intervals(intervals, samples):
    #test out if equal intervals include all samples
    verdict = True
    ints = []
    for i in intervals: ints.extend(i)
    if len(ints) != len(samples):
        print(f"Warning: With the current replicate settings you are loosing samples!\nYou specify {len(intervals)} intervals for {len(ints)} entries but data has {len(samples)} entries...")
        verdict = False
    return verdict 


def _prep_replicates(replicates):
    """Establishes what kind of replicates are passed"""
    try: 
        reps = int(replicates)
        equal_reps = True
    except TypeError:
        reps = list(replicates)
        equal_reps = False
    return reps,equal_reps

#interval auxiliaries for replicate finding

#creates equalually spaces index groups
def _create_equal_intervals(length, reps):
    indices = []
    i = 0
    while i + reps <= length:
        indices.append(list(range(i, i+reps)))
        i+= reps
    return indices

#creates index groups based on the entry list provided in reps
def _create_unequal_intervals(reps):
    indices = []
    i = 0
    for r in reps:
        tmp = list(range(i, i+r))
        i+=r
        indices.append(tmp)
    return indices


def group_samples(sample_dict, replicates):
    """
    Second step in manual qpcr analysis. Groups samples into replicate groups as defined by replicates.
    """
    samples = _define_replicates(sample_dict, "Sample", replicates)
    cts = _define_replicates(sample_dict, "Mean Ct", replicates)
    
    groupings = {}
    for i in range(0,len(samples)):
        key = "Group {}".format(i+1)
        values = [samples[i],cts[i]]
        tmp = {key : values}
        groupings.update(tmp)
    return groupings



def Delta_Ct(grouped_dict, exp=True, anchor="first"):
    """
    Calculates Delta_Ct within one sample_dict. Requires group_samples to have been performed before use.
    exp allows to adjust to either 2^(d1-d2) (True) or only (d1-d2) (False).
    anchor sets the reference (d1). Possible are None (group internal, default), "first" (very first entry), or any specified numeric value. 
    """
    keys = list(grouped_dict.keys())
    new_dict = {}
    if exp == True:
        dCt = _exp_DCt
    else:
        dCt = _simple_DCt
        
    if anchor is None:
        _grouped_anchored(grouped_dict, keys, new_dict, dCt)
    elif anchor == "first":
        _first_anchored(grouped_dict, keys, new_dict, dCt)
    else:
        _specified_anchored(grouped_dict, anchor, keys, new_dict, dCt)
    return new_dict

# functions to calculate delta ct based on the anchor specified 

def _specified_anchored(grouped_dict, anchor, keys, new_dict, dCt):
    first = anchor
    for k in keys:
        cts = grouped_dict[k][1]
        tmp = [dCt(first, i) for i in cts]
        new_dict.update({k : tmp})

def _first_anchored(grouped_dict, keys, new_dict, dCt):
    first = grouped_dict[keys[0]][1][0]
    for k in keys:
        cts = grouped_dict[k][1]
        tmp = [dCt(first, i) for i in cts]
        new_dict.update({k : tmp})

def _grouped_anchored(grouped_dict, keys, new_dict, dCt):
    for k in keys:
        cts = grouped_dict[k][1]
        first = cts[0]
        tmp = [dCt(first, i) for i in cts]
        new_dict.update({k : tmp})

def _exp_DCt(ref, sample):
    factor = sample-ref 
    return 2**(-factor)

def _simple_DCt(ref, sample):
    return sample-ref


def normalise(normaliser:dict, sample:dict, no_head=True):
    """
    This function normalises a sample_dict against a corresponding normaliser. 
    Requires that both have undergone Delta_Ct before. 
    This function corresponds to the second Delta in DeltaDelta CT analysis of qPCR data.
    no_head=False would add "Legend" : "DDCT" as very first entry to the returned dictionary.
    """
    if no_head == False:
        normalised = {"Legend" : "DDCT"}
    else: 
        normalised = {}
    keys_n = normaliser.keys()
    keys_s = sample.keys()
    if keys_n != keys_s:
        print("Normaliser and Samples do not share the same group names!")
        return None
    for k in keys_s:
        normals = normaliser[k]
        samples = sample[k]
        if len(normals) != len(samples):
            print("Normaliser and Samples do not have the same length!")
            return None
        
        temp = []
        for i in range(0, len(normals)):
            n = normals[i]
            s = samples[i]
            rel = s/n
            temp.append(rel)
        
        normalised.update({k : temp})
    return normalised

#if several normalisers should be used they first have to be pre-processed to get their Delta_Ct values before they can be combined using combine_normalisers
def preprocess_normalisers(normalisers:list, replicates, run_names=None, group_names=None, anchor=None, return_type="list"):
    """
    This function takes in a list of input csv files of normalisers and automatically performes grouping, Delta_Ct and renaming.
    Note that the normalisers have to be named and group_names must correspond to group_names later used for the actual samples!
    """
    normalisers_dict = {}
    ndx = 0
    if run_names is None: 
        names = _ddCt_generate_run_names(normalisers)
        run_names = names

    for norm in normalisers:
        name = run_names[ndx]
        tmp = open_csv_file(norm)
        tmp = group_samples(tmp, replicates=replicates)
        if group_names is not None:
            tmp = rename_groups(tmp, new_names=group_names)
        tmp = Delta_Ct(tmp, anchor=anchor)
        tmp = {name : tmp}
        normalisers_dict.update(tmp)
        ndx +=1
    if return_type == "list":
        temp = []
        for i in run_names:
            temp.append(normalisers_dict[i])
        normalisers_dict = list(temp)
    return normalisers_dict

#if several normalisers should be averaged, use combine_normalisers
def combine_normalisers(normalisers:list):
    combined_normaliser = {}
    keys = list(normalisers[0].keys()) #get all groupings
    for k in keys:
        temp = []
        for norm in normalisers:
            length = len(norm[k])
            assert isinstance(norm[k][0], float), "norm[k][0] is not a simple number but: {} ... \nAre you sure you already performed qpcr.Delta_Ct on your normaliser?".format(norm[k][0])
            for i in norm[k]:
                temp.append(i)
        temp = {k : [stat.mean(temp) for i in range(length)]} #we must conserve dimensionality...
        combined_normaliser.update(temp)
    return combined_normaliser


def rename_groups(sample_dict, new_names:list):
    """
    Sample names from the original csv file are not imported! 
    Samples will be only named according to their groups. Group names can be set using this function. If no names are provided, "Group 1", "Group 2" etc. will be used by default.
    rename_groups can be used on any dict at any stage in the process. 
    """
    keys = sample_dict.keys()
    new_dict = {}
    i = 0
    for k in keys:
        new_name = new_names[i]
        tmp = { new_name : sample_dict[k] }
        new_dict.update(tmp)
        i +=1
    return new_dict


def get_stats(deltaCT_dict, export = ["avg", "stdv"]):
    """
    By default Delta_Ct (and normalise afterward) will work with individual replicate values, and return a dict of the same dimensions.
    To concentrate group information from many replicates directly to average and stdev, use get_stats.
    This will return a dict where each group assigns a list of two values, the first being average, the second being stdev.
    """
    keys = deltaCT_dict.keys()
    new_dict = {}
    legend = []
    if "avg" in export: legend.append("Avg")
    if "stdv" in export: legend.append("Stdev")
    if "med" in export: legend.append("Med")
    new_dict.update({"Legend" : legend})
    for k in keys:
        dcts = deltaCT_dict[k]
        avg = stat.mean(dcts)
        std = stat.stdev(dcts)
        med = stat.median(dcts)
        
        tmp = []
        if "avg" in export:
            tmp.append(avg)
        if "stdv" in export:
            tmp1 = [i for i in tmp]
            tmp1.append(std)
            tmp = tmp1
        if "med" in export:
            tmp1 = [i for i in tmp]
            tmp1.append(med)
            tmp = tmp1

        tmp = {
            k : tmp
        }
        new_dict.update(tmp)
    return new_dict


def export_to_csv(data, filename, transpose=False):
    """
    This function writes a csv file with the results generated. It may be called at any stage in the analysis process to save intermediary data.
    transpose will transpose columns and rows.
    """
    export_data = pd.DataFrame(data)
    if transpose == True:
        export_data = export_data.transpose()
    export_data.to_csv(filename)
    #and now re-open to remove the first line that only contains 0 1 and empty commas
    contents = []
    with open(filename, "r") as openfile:
        for i in openfile:
            contents.append(i)
    contents = contents[1:]
    with open(filename, "w") as openfile:
        if "Legend" not in contents[0]: 
            openfile.write("Sample,Value\n")
        for i in contents: 
            openfile.write(i)

#export grouped raw CT values to plot if desired
def export_raw_data(filename, replicates, group_names=None, export_location=None):
    """
    To also group raw data, use this function. This function reads original csv files, groups them and exports them back to the same location, creating a new file with _raw.csv as ending.
    """
    contents_dict = open_csv_file(filename)
    contents_dict = group_samples(contents_dict, replicates=replicates)
    if group_names is not None:
        contents_dict = rename_groups(contents_dict, group_names)
    tmp = {}
    for k in list(contents_dict.keys()):
        tmp.update({k : [str(i) for i in contents_dict[k][1]]  })
    savefile_dict = tmp

    if export_location is None:
        new_file = "{}_raw.csv".format(filename.replace(".csv", ""))
    else:
        new_file = export_location
    
    with open(new_file, "w") as openfile:
        for k in list(savefile_dict.keys()):
            newstr = "{},{}\n".format(k, ",".join(savefile_dict[k]))
            openfile.write(newstr)


#if one wants to revisit already computed results
def load_results(filename, mode="individual", *args, **kwargs):
    """
    This function opens pre-computed results of the ops.biotools.qpcr Module from the generated csv files. 
    It supports two modes: 
    mode = "individual" (default) where filename specifies the file to be opened. It returns a dictionary containing the grouped computed values (replciates, or avg, stdev).
    mode = "pairs" allows users to load an entire set of samples and normalisers that were previously computed, to then be used by ops.biotools.qpcr.Analysis.normalise_pairs. It returns two separate dictionaries each containing the set of sample files it was given by kwargs parameters samples:list, normalisers:list. Optionally, names may be additionally assigned to samples and normalisers using sample_names:list and norm_names:list. 

    To facilitate working with replicate assays (i.e. same qPCR assay normalised against different normalisers separately), samples / normalisers support only partial naming and need no full filepath to function. Like this multiple analysis result with e.g. "HNRNPL NMD" in their name will all be loaded. 
    """
    if mode == "individual":

        export_dict = _load_individual_results(filename)

    elif mode == "pairs":
        #get necessary kwargs for _load_pair_results
        location, samples, normalisers, sample_names, norm_names = _check_loadpair_kwargs(filename=filename, **kwargs)

        #load samples into dict
        export_dict = _load_pair_results(location=location, 
                                            samples=samples, 
                                            normalisers=normalisers, 
                                            sample_names=sample_names, 
                                            norm_names=norm_names)
    
    elif mode == "multiple":
        assert isinstance(filename, list), "When using mode='multiple', please provide a list of valid filenames for the filename argument"
        
        keys = oo.from_kwargs("keys", None, **kwargs)
        export_dict = _bulk_load(filename, keys=keys)
           
    else:
        print("Please, use either mode='individual', mode='multiple', or mode = 'pairs' to load your results")
        return None

    return export_dict

def _load_individual_results(filename):
    contents = []
    with open(filename, "r") as openfile:
        contents = openfile.read().split("\n")

    contents = [i.split(",") for i in contents]
    conditions = [i[0] for i in contents]
    values = [i[1:] for i in contents]
    values = values[1:] #crop the first labelling
    values = [[float(j) for j in i] for i in values] #get the numeric values back
    
    #now pack everything back into a dict
    export_dict = {}
    idx = 0
    for i in conditions[1:]:
        tmp = {i : values[idx]}
        export_dict.update(tmp)
        idx +=1
    return export_dict

def _check_loadpair_kwargs(filename, samples, **kwargs):
    pair_args = oo.get_kwargs(_load_pair_results, **kwargs)
    location = filename
    
    normalisers = oo.from_kwargs(("norms", "normalisers"), None, pair_args)
    sample_names = oo.from_kwargs("sample_names", None, pair_args)
    norm_names = oo.from_kwargs("norm_names", None, pair_args)

    if normalisers is None:
        raise NameError("No normalisers were found in kwargs! Make sure to provide normalisers using either 'norms' or 'normalisers' arguments!")
    return location, samples, normalisers, sample_names, norm_names

def _load_pair_results(location, samples, normalisers, sample_names = None, norm_names = None):
    if isinstance(location, (tuple, list)) and len(location) == 2:
        sample_loc = location[0]
        norm_loc = location[1]
    else:
        sample_loc = location
        norm_loc = location
    
    
    sample_items = fc.item_finder(path=sample_loc, filter=[".csv"])
    sample_items.sort()
    sample_list = fc.filter_files(items_list=sample_items, target_names=samples)
    
    norm_items = fc.item_finder(path=norm_loc, filter=[".csv"])
    norm_items.sort()
    norm_list = fc.filter_files(items_list=norm_items, target_names=normalisers)
    
    sample_dict = _bulk_load(sample_list, keys=sample_names)
    norm_dict = _bulk_load(norm_list, keys=norm_names)

    return sample_dict, norm_dict

def _bulk_load(filenames, keys=None):
    return_dict = {}
    k = 0
    print("Loading Result files...")
    if keys is not None:
        assert len(keys) == len(filenames), "Keys and filenames lists are not of equal length!"
        for f in filenames:
            name = keys[k]
            print("""
            {}\t\t Src: {}
            """.format(name, f))
            tmp = load_results(f)
            return_dict.update({name : tmp})
            k+=1
    else:
        for f in filenames:
            name = "Assay {}".format(k)
            print("""
            {}\t\t Src: {}
            """.format(name, f))
            tmp = load_results(f)
            return_dict.update({name : tmp})
            k+=1
    return return_dict            

# generate run_names for qA.delta_deltaCt in case the user does not want to assign manually
def _ddCt_generate_run_names(data_files):
    target_names = [_remove_suffix(i) for i in data_files]
    return target_names

def _remove_suffix(file):
    if "/" in file:
        norm_name = file.split("/")
        norm_name = norm_name[-1]
    norm_name = norm_name.split(".")
    norm_name = norm_name[0]
    return norm_name


def help():
    from ops.biotools.qpcr import Analysis
    Analysis.help()
    
def Info():
    from ops.biotools.qpcr import Analysis
    Analysis.Info()

def Example():
    from ops.biotools.qpcr import Analysis
    Analysis.Example()



if __name__ == "__main__":
    test = "Example Data/28S.csv"

    a = open_csv_file(test)

    grouped = group_samples(a, 3)
    # print(grouped)

    loaded = load_results("./Example Data/HNRNPL_nmd_DeltaDelta_Ct.csv")
    print(loaded)
