"""
Stores basic default settings
"""

#  =================================================================
#                       Default FileIO settings
#  =================================================================

# default column names for raw Ct data files (don't change this!)
raw_col_names = ["id", "Ct"]

supported_filetypes = ["csv", "xlsx"]



#  =================================================================
#                       Default Figure settings
#  =================================================================

# PreviewResults
# --------------------------------
static_PreviewResults = dict(
                                color = "white",
                                edgecolor = "navy",
                                edgewidth = 1,
                                rot = 0,
                                legend = False, 
                                frame = False,
                                title = "Preview of Results"
                            )

interactive_PreviewResults = dict(
                                    template = "plotly_white",
                                    title = "Preview of Results"
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