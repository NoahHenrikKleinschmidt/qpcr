"""
Stores basic default settings
"""

#  =================================================================
#                       Default FileIO settings
#  =================================================================


supported_filetypes = ["csv", "xlsx"]



#  =================================================================
#                       Default Inference settings
#  =================================================================

# default column names for raw Ct data files (don't change this!)
raw_col_names = ["id", "Ct"]

default_dataset_header = "assay"
default_group_name = "group{}"

# default assumed column names for Id and Ct values
default_id_header = "Name"
default_ct_header = "Ct"


#  =================================================================
#                       Default core settings
#  =================================================================

# default seed for randomised processes
default_seed = 11299114


#  =================================================================
#                       Default Calibrator settings
#  =================================================================

# the default prefix for calibrator replicates that must be specified 
# within the group names.
calibrator_prefix = "calibrator"


#  =================================================================
#                       Default Figure settings
#  =================================================================

# AssayDots / GroupDots
# -------------------------------

static_PreviewDots = dict(
                                    title = "Preview of Results",
                                    ylabel = "norm. Fold Change",
                                    linewidth = 0.8,
                                    rot = 30,
                                    frame = False,
                                    palette = "Blues"
                    )
                
interactive_PreviewDots = dict(
                                    title = "Preview of Results",
                                    ylabel = "norm. Fold Change",
                                    template = "plotly_white",
                            )

# Assay / GroupBars
# -------------------------------

static_PreviewBars = dict(
                                    title = "Preview of Results",
                                    ylabel = "norm. Fold Change",
                                    linewidth = 0.8,
                                    style = "dark",
                                    frame = False,
                                    legend = False,
                                    rot = 30,
                                    color = "Blues"
                    )
                
interactive_PreviewBars = dict(
                                    title = "Preview of Results",
                                    ylabel = "norm. Fold Change",
                                    template = "plotly_white",
                            )

# ReplicateBoxPlot
# --------------------------------
static_ReplicateBoxPlot = dict(
                                    title = "Summary of Replicates",
                                    linewidth = 0.8,
                                    frame = False,
                                    style = "dark",
                                    palette = "Blues"
                                )

interactive_ReplicateBoxPlot = dict(
                                        template = "plotly_white",
                                        title = "Summary of Replicates"
                                    )

# FilterSummary
# --------------------------------
static_FilterSummary = dict(
                                    title = "Filtering Summary",
                                    linewidth = 0.8,
                                    frame = False,
                                    style = "dark", 
                                    rot = 30,
                                    palette = "coolwarm"
                                )

interactive_FilterSummary = dict(
                                        template = "plotly_white",
                                        title = "Filtering Summary"
                                    )

# EfficiencyLines
# --------------------------------
static_EfficiencyLines = dict(
                                    title = "Efficiency Computations",
                                    xlabel = "Log Dilution",
                                    ylabel = "Ct",
                                    palette = "Blues",
                                    style = "dark",   
                            )

interactive_EfficiencyLines = dict(
                                    title = "Efficiency Computations",
                                    xlabel = "Log Dilution",
                                    ylabel = "Ct",
                            )