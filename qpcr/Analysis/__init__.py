import qpcr
import matplotlib.pyplot as plt
import statistics as stat 
import qpcr.aux.graphical.auxiliaries as gaux
import pandas as pd 
import difflib

#now let us try to warp this all up into a more user friendly package
def single_deltaCt(data_file, replicates, mode = 'replicate', transpose=True, export=True, group_names=None, stats = ['avg', 'stdv'], anchor = None, dCt_exp = True, exportname_addon=None):
    contents = qpcr.open_csv_file(data_file)
    grouped_dict = qpcr.group_samples(sample_dict=contents, replicates=replicates)
    deltaCt_dict = qpcr.Delta_Ct(grouped_dict=grouped_dict, anchor=anchor, exp=dCt_exp)
    if group_names is not None:
        if group_names == "auto":
            deltaCt_dict = qpcr.rename_groups(sample_dict=deltaCt_dict, new_names=contents['Sample'])
        else:
            deltaCt_dict = qpcr.rename_groups(sample_dict=deltaCt_dict, new_names=group_names) 
    if mode == 'replicate':
        export_dict = deltaCt_dict
    elif mode == 'stats':
        export_dict = qpcr.get_stats(deltaCt_dict, export=stats)
    
    if export == True:
        #new_file = '{}_SingleCT.csv'.format(data_file.replace('.csv', ''))
        if exportname_addon is None: 
            addon = ""
        else: 
            addon = "{}_".format(exportname_addon)
        new_file = "{}_{}SingleDelta_Ct.csv".format(data_file.replace('.csv', ''), addon)
            
        qpcr.export_to_csv(data=export_dict, filename=new_file, transpose=transpose)
        print('Exported SingleDelta_Ct analysis to:\n{}'.format(new_file))
    else:
        return export_dict


def delta_deltaCt(data_files:list, replicates, normaliser, run_names="auto", mode = 'replicate',  group_names=None, transpose=True, stats = ['avg', 'stdv'], anchor = None, dCt_exp = True, export=True, exportname_addon=None, export_location=None):
    samples_dict = {}
    i = 0
    for f in data_files:
        key = 'Group {}'.format(i)
        delta_Cts = single_deltaCt(f, replicates=replicates, mode='replicate',
                                    transpose=transpose, export=False,
                                    group_names=group_names, stats=stats, anchor=anchor, dCt_exp=dCt_exp)
        tmp = { key : delta_Cts }
        samples_dict.update(tmp)
        i+=1
    
    #rename dict to the gene names specified
    if run_names == "auto":
        names = qpcr._ddCt_generate_run_names(data_files)
        run_names = names
    samples_dict = qpcr.rename_groups(samples_dict, new_names=run_names)
    
    keys = list(samples_dict.keys())

    if isinstance(normaliser, dict):
        normaliser_dict = normaliser
    else: 
        if normaliser not in keys:
            print('No match for the normaliser could be found! Make sure to have the normaliser name also specified identically in the gene_names.')
            return None
        normaliser_dict = samples_dict[normaliser]
        keys.remove(normaliser) #remove normaliser from list of samples
    
    #now normalise relative to normaliser
    normalised_dict = {}
    for k in keys:
        temp = qpcr.normalise(normaliser=normaliser_dict, sample=samples_dict[k])
        temp = {k : temp}
        normalised_dict.update(temp)
    
    if mode == "stats":
        stats_dict = {}
        for k in normalised_dict.keys():
            normdict = normalised_dict[k]
            tmp = qpcr.get_stats(normdict, export=stats)
            stats_dict.update({k : tmp})
        normalised_dict = stats_dict

    if export==True:
        #now save all normalised Delta_Delta Ct Values to new files
        file_index = 0
        for k in normalised_dict.keys():
            
            if export_location is None:
                file_location = data_files[file_index].split("/")
                file_location = "/".join(file_location[0:len(file_location)-1])
            else: 
                file_location = export_location
            
            if exportname_addon is None: 
                addon = ""
            else: 
                addon = "{}_".format(exportname_addon)
            new_file = "{}/{}_{}DeltaDelta_Ct.csv".format(file_location, k, addon)
            
            qpcr.export_to_csv(data=normalised_dict[k],filename=new_file, transpose=transpose)
            print('Exported DeltaDelta_Ct analysis to:\n{}'.format(new_file))
        
    return normalised_dict

 

def normalise_pairs(samples:dict, normalisers:dict, pair_names=None, export=True, export_location=None, transpose=True):
    """
    This function allows users to normalise entire sets of sample dicts against a corresponding set of normalisers. This is 
    especially useful when genes of different primer pairs (Assays) were measured separately and are now supposed to be normalized
    against each other. 
    """
    sample_keys = list(samples.keys())
    norm_keys = list(normalisers.keys())
    assert len(sample_keys) == len(norm_keys), "The samples and normalisers do not share the same length!"

    if pair_names is None:
        names = []
    elif len(pair_names) == len(sample_keys):
        names = pair_names
    else:
        raise IndexError("The pair_names provided do not cover every pair!")

    export_dict = {}
    for i in range(len(sample_keys)):
        skey = sample_keys[i]
        nkey = norm_keys[i]
        
        if pair_names is None:
            name = "{sample}_rel_{norm}".format(sample=skey, norm=nkey)
            names.append(name)
        else:
            name = names[i]

        string = """
        Acting on:\t{}
        Sample:\t\t{}
        Normaliser:\t{}
        """.format(name, skey, nkey)
        print(i, string)        
       
        sample = samples[skey]
        norm = normalisers[nkey]
        tmp = qpcr.normalise(sample=sample, normaliser=norm)
        export_dict.update({name : tmp})

    if export == True:
        if export_location is None:
            raise Exception("Missing Export Location!")
        for name in names:
            data = export_dict[name]
            filename = "{}/{}.csv".format(export_location, name)
            if "//" in filename: 
                filename = filename.replace("//", "/")
            qpcr.export_to_csv(data=data, filename=filename, transpose=transpose)
    
    return export_dict


def preview_results(results_dict, transpose=False, figsize = (8,8)):
    """
    This function plots the results of all Delta-Delta Ct analyses as individual bar charts into one figure.
    It returns a matplotlib.pyplot.figure object.
    """
    results_number = len(list(results_dict.keys()))
    rows, cols = gaux.adjust_layout(graph_number=results_number)

    if transpose == True:
        row_count = cols
        col_count = rows
        fig, axs = plt.subplots(cols, rows, squeeze=False, sharex=True, figsize=figsize)
    else:
        row_count = rows
        col_count = cols
        fig, axs = plt.subplots(rows, cols, squeeze=False, figsize=figsize)
    
    coordinates = []
    for r in range(0, row_count):
        for c in range(0, col_count):
            coordinates.append([r,c])
    #print('Creating a {}x{} Layout'.format(row_count, col_count))

    cdx = 0
    for k in list(results_dict.keys()):
        title = k
        sample_dict = results_dict[k]
        sample_keys = list(sample_dict.keys())
        if 'Legend' == sample_keys[0]:
            values = [sample_dict[i][0] for i in sample_keys[1:]]
            yerrs = [sample_dict[i][1] for i in sample_keys[1:]]
            labels =sample_keys[1:]
        else:
            values = [sample_dict[i] for i in sample_keys]
            labels = sample_keys
            yerrs = [stat.stdev(i) for i in values]
            values = [stat.mean(i) for i in values]

        r = coordinates[cdx][0]
        c = coordinates[cdx][1]
        #print(r, c)
        axs[r, c].bar(x=labels, height=values, color='lightgray', edgecolor = "dimgrey", linewidth=1)
        axs[r, c].errorbar(x = labels, y = values, yerr = yerrs, fmt = ".", markersize = 0, capsize = 3, ecolor="black", )
        axs[r, c].set(title=title)
        cdx +=1

    fig.tight_layout()
    plt.show()
    return fig

def make_aggregate_plot(result, figsize=(9,5), colormap = "GnBu_r", title=""):
    """
    This function creates a grouped bar chart of all Delta-Delta Ct analyses into one figure. 
    It returns a figure object.
    """
    result_df = _convert_to_stats(result)
    r = pd.DataFrame() # dataframe for the Ct values
    stv = pd.DataFrame() # dataframe for the StDevs

    for k in result_df.keys():
        tmp = pd.DataFrame(result_df[k])
        del tmp["Legend"] # remove the "Legend" column
        r[k] = tmp.iloc[0]
        stv[k] = tmp.iloc[1]

    fig, ax = plt.subplots(figsize=figsize)
    ax.set(ylabel="$\Delta\Delta C_T$", title=title)
    r.plot.bar(yerr = stv, ax=ax, colormap = colormap,rot=0, edgecolor = "black", linewidth = 1)
    plt.close()
    return fig

# if the result_dict is not yet in stats mode it will be converted into stats
def _convert_to_stats(result):
    try: 
        keys = list(result.keys())
        tmp = pd.DataFrame(result[keys[0]])
        tmp["Legend"]
        result_df = result
    except:
        conv_result = {}
        for k in result:
            tmp = qpcr.get_stats(result[k])
            conv_result.update({k : tmp})
        result_df = conv_result
    return result_df

def make_grouped_plots(result, subplots=True, figsize=(8,8), colormap = "GnBu"):
    """
    This function generates grouped barplots for related samples such as "GeneX_ctrl" and "GeneX_treatment". 
    The function (rather the two _ functions) tries to assess which samples might belong together by reading their run_names (assay names). 
    If subplot=True is set, it will produce one single figure with subplots (just like preview_results). If subplot=False it will produce one figure for each sample-grouping and return a list of figures.
    """
    if subplots == True:
        # fig will be only one figure
        fig = _make_grouped_plots_subplots(result, figsize, colormap)
    elif subplots == False:
        # fig will be a list of figures !
        fig = _make_grouped_plots_individuals(result, figsize, colormap)
    return fig


def _make_grouped_plots_individuals(result, figsize, colormap):
    """
    This function creates a grouped bar graph for each sample grouping (stored in the list all_matches). 
    It will then generate a figure object for each grouping and store and return these as a list.
    """
    result_df = _convert_to_stats(result)
    # extract which samples might belong together
    all_matches = find_matches(result_df)

    figs = []
    for m in all_matches:
        r = pd.DataFrame()
        stv = pd.DataFrame()
        for k in m:
            tmp = pd.DataFrame(result_df[k])
            del tmp["Legend"]
            r[k] = tmp.iloc[0]
            stv[k] = tmp.iloc[1]
        fig, ax = plt.subplots(figsize = figsize)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set(ylabel="$\Delta\Delta C_T$")
        r.plot.bar(yerr = stv, ax = ax, colormap = colormap,rot=0, edgecolor = "black", linewidth = 1)
        figs.append(fig)
        plt.close()
    return figs


def _make_grouped_plots_subplots(result, figsize, colormap):
    """
    This function will generate a grouped bargraph for each sample grouping (stored in the list all_matches) and plot these as subplots into one single figure which is returned.
    """
    result_df = _convert_to_stats(result)
    # extract which samples might belong together
    all_matches = find_matches(result_df)

    # generate a subplots figure
    rows, cols = gaux.adjust_layout(len(all_matches))
    fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=figsize)

    # assemble coordinates for the subplot axes
    coordinates = []
    for r in range(rows): 
        for c in range(cols):
            coordinates.append([r, c])
    
    # fill each ax with a grouped bar plot
    cdx = 0    
    for m in all_matches:

        r = pd.DataFrame() # dataframe for Ct values (heights)
        stv = pd.DataFrame() # dataframe for stdevs
        for k in m:
            tmp = pd.DataFrame(result_df[k])
            del tmp["Legend"] # remove the "legend" column
            r[k] = tmp.iloc[0] # store Cts
            stv[k] = tmp.iloc[1] # store stdevs
        
        # assign coordinate of subplot ax
        print(coordinates[cdx])
        _r, _c = coordinates[cdx]

        #ax.spines["right"].set_visible(False)
        #ax.spines["top"].set_visible(False)
        ax[_r, _c].set(ylabel="$\Delta\Delta C_T$")
        r.plot.bar(yerr = stv, ax = ax[_r, _c], colormap=colormap,rot=0, edgecolor = "black", linewidth = 1)
        fig.tight_layout()
        plt.close()
        cdx +=1

    return fig

def find_matches(result_df):
    """
    This function tries to estimate which samples (runs) belong together by estimating the similarity of their run_names.
    """
    all_matches = []
    for k in result_df.keys():
        matches = difflib.get_close_matches(k, result_df.keys(), cutoff=0.63)
        matches.sort()
        if matches not in all_matches: # store groups of related samples as list
            all_matches.append(matches)
    return all_matches

def help():
    helpstring = """
This package is designed to help analyse qPCR data generated by Qiagen RotorGene® 
and works with the Excel Spreadhseet exported data from this device. 

1. Prepare Data Sheet as csv
    -> copy only the lanes with Sample Names, Rep. Ct, and Rep. Ct Std. into a new Excel Sheet
    -> save this new sheet as .csv -- this is your source_file

2. Either perform a predefined Analysis or analyise manually
    -> Predefined Analysis: import ops.biotools.qpcr.Analysis

3. Manual Analysis:
    -> read your source file using qpcr.open_csv_file
    -> group your replicate samples using qpcr.group_samples
        >> make sure to enter either an integer if replicates are equal (always triplicate)
           or provide a tuple where the replicate number is specified (2, 3) = a douplicate then a triplicate
        >> if your replicate samples need to be renamed use qpcr.rename_groups
    -> now calculate Delta_Cts using qpcr.Delta_Ct
        >> for additional information about the Delta_Ct method used use qpcr.Info
    -> to normalise relative to a normaliser use qpcr.normalise
        >> note qpcr.normalise takes in two Delta_Ct dicts whose groups have the same names!
    -> if you wish to obtain Average/Stdv or Median of your calculations use qpcr.get_stats
    -> to export your analysis to a new csv file use qpcr.export_to_csv
        >> use transpose=False to have each Sample Group appear as one column
    """
    print(helpstring)

def Info():
    infostring = """
The Delta Ct method used to compute Delta Ct and Delta Delta Ct values is as follows:

1. Delta Ct: 
Each sample condition (i.e. WT+, WT– etc.) is treated as a 'Group' and each group is being related to its first respective entry. 
Hence, if WT+ is a triplicate measurement, in a first instance each will be normalised to the first WT+ Ct.
The formula is 𝛥Ct = 2^(nth_replicate - first_replicate), generating list of relative replicate Cts for each sample condition (always starting with 1.0)


2. Delta Delta Ct: 
Each gene of interest will be first analysed by Single Delta Ct as described above. Subsequently each sample condition (Group) is taken relative to the same Group  within the reference gene (normaliser). 
The formula applied is 𝛥𝛥Ct = 𝛥Ct(gene of interest) / 𝛥Ct(reference gene), generating a list of replicate 𝛥𝛥Cts relative within each group and to the reference gene. 


If different anchors for relation should be used, please use the functions defined in the ops.biotools.qpcr module to create your own workflow.
-> qpcr.Delta_Ct has the option exp=False to only compute Ct1-Ct0 instead of 2^(Ct1-Ct0)
-> qpcr.Delta_Ct also takes in an optional anchor argument. There any numeric value may be specified as anchor for Delta Ct computation, alternatively using anchor = "first" one may simply use the very first Delta Ct of the dataframe as reference, thereby ignoring Groupings.

The delta_deltaCt and single_deltaCt functionS accept both exp (as dCt_exp) and anchor arguments but their applicability may be limited within the fixed workflow of these functions.

    """
    print(infostring)

def Example():
    examplestring = """
Example of Manual usage of ops.biotools.qpcr module functions to perform Delta_Delta_Ct analysis

======================================
import ops.biotools.qpcr as qpcr

#define source_files
actin = "/Users/NoahHK/OneDrive - Universitaet Bern/Bachelor Project/qPCR/test_Actin.csv"
genex = "/Users/NoahHK/OneDrive - Universitaet Bern/Bachelor Project/qPCR/test_GeneX.csv"

#open files and group files into triplicates
actin_dict = qpcr.open_csv_file(actin)
genex_dict = qpcr.open_csv_file(genex)

actin_dict = qpcr.group_samples(actin_dict, 3)
genex_dict = qpcr.group_samples(genex_dict, 3)

#Compute Delta_Cts for each gene (relative within each triplicate condition)
actin_dict = qpcr.Delta_Ct(actin_dict)
genex_dict = qpcr.Delta_Ct(genex_dict)

#normalise GeneX relative to ActinB
normalised_dict = qpcr.normalise(normaliser=actin_dict, sample=genex_dict)

#get Average and StDev (= predefined settings of qpcr.get_stats)
normalised_dict = qpcr.get_stats(normalised_dict)

#inspect data and export to a new csv file
print(normalised_dict)
qpcr.export_to_csv(data = normalised_dict, filename= "/Users/NoahHK/OneDrive - Universitaet Bern/Bachelor Project/qPCR/test_analysis.csv")
======================================

Example of using a predefined Analysis routine to perform Delta_Delta_Ct analysis

======================================
import ops.biotools.qpcr.Analysis as qA

#define source_files
actin = "/Users/NoahHK/OneDrive - Universitaet Bern/Bachelor Project/qPCR/test_Actin.csv"
genex = "/Users/NoahHK/OneDrive - Universitaet Bern/Bachelor Project/qPCR/test_GeneX.csv"

#run predefined Delta_Delta_Ct procedure
qA.delta_deltaCt(data_files=[actin, genex], replicates=3, 
                normaliser='ActinB', gene_names=['ActinB', 'GeneX'], 
                mode = "stats")
   
    """
    print(examplestring)
