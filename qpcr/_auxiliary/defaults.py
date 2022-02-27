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
#                       Default Figure settings
#  =================================================================

# PreviewResults
# --------------------------------
static_PreviewResults = dict(
                                color = "white",
                                edgecolor = "xkcd:pure blue",
                                ecolor = "black",
                                edgewidth = 1,
                                rot = 0,
                                legend = False, 
                                frame = False,
                                title = "Preview of Results",
                                ylabel  = "norm. Fold Change",
                            )

interactive_PreviewResults = dict(
                                    template = "plotly_white",
                                    title = "Preview of Results",
                                    ylabel = "norm. Fold Change"
                                )

# ReplicateBoxPlot
# --------------------------------
static_ReplicateBoxPlot = dict(
                                    title = "Summary of Replicates",
                                    linewidth = 0.8,
                                    frame = False,
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
                                    frame = True,
                                    palette = "coolwarm"
                                )

interactive_FilterSummary = dict(
                                        template = "plotly_white",
                                        title = "Filtering Summary"
                                    )

# PreviewDots
# -------------------------------

static_PreviewDots = dict(
                                    title = "Preview of Results",
                                    ylabel = "norm. Fold Change",
                                    linewidth = 0.8,
                                    frame = False,
                                    palette = "Blues"
                    )
                
interactive_PreviewDots = dict(
                                    title = "Preview of Results",
                                    ylabel = "norm. Fold Change",
                                    template = "plotly_white",
                            )