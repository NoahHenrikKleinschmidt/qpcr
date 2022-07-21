"""
Stores basic default settings
"""

import logging


#  =================================================================
#                       Default FileIO settings
#  =================================================================


supported_filetypes = ["csv", "xlsx"]



#  =================================================================
#                       Default Inference settings
#  =================================================================

raw_col_names = ["id", "Ct"]
"""default column names for raw Ct data files (don't change this!)"""

dataset_header = "assay"
"""Default assumed column name for dataset/assay entries in BigTables"""

group_name = "group{}"
"""Default name to use when group names cannot be inferred"""

# default assumed column names for Id and Ct values
id_header = "Name"
"""Default assumed column name for replicate id entries"""
ct_header = "Ct"
"""Default assumed column name for replicate Ct entries"""

setup_cols = ["id", "group", "group_name"]
"""Default Assay dataset columns that supply metadata, and not specific values like Ct values."""

#  =================================================================
#                       Default core settings
#  =================================================================

seed = 11299114
"""default seed for randomised processes"""

strict_id = False
"""Set to True if object ids may strictly only be set once!"""

log_level = logging.WARNING
"""The default logging level"""

log_format = "%(asctime)s  |  %(levelname)s  |  %(module)s.%(funcName)s (%(lineno)d)  |  %(message)s"
"""The default logging format"""

init_log_loc = "stdout"
"""The default logging location."""

init_log_format = "%(levelname)s  |  %(module)s.%(funcName)s  |  %(message)s"
"""The default initial logging format that is used if aux.log() has not been called yet"""


#  =================================================================
#                       Default Calibrator settings
#  =================================================================

calibrator_prefix = "calibrator"
"""the default prefix for calibrator replicates that must be specified within the group names."""


#  =================================================================
#                       Default Figure settings
#  =================================================================

plotmode = "static"
"""The default setting for interactive or static figures"""

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
                                    color = "cornflowerblue",
                                    style = "dark",   
                            )

interactive_EfficiencyLines = dict(
                                    title = "Efficiency Computations",
                                    xlabel = "Log Dilution",
                                    ylabel = "Ct",
                            )