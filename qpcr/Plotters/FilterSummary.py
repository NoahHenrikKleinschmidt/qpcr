"""
Defines the `FilterSummary` class that is responsible for showing the before-after Cts of filtered assays.
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import plotly.graph_objs as go


import qpcr.defaults as defaults
from qpcr import _auxiliary as aux
import qpcr._auxiliary.graphical as gx
import qpcr.Plotters._base as base
import qpcr.main.Results as Results
import qpcr.main.Assay as Assay


class FilterSummary(base.AssayPlotter):
    """
    Generates a summary figure of the replicate
    Ct values for each assay before and after filtering.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).



    `"static"` Kwargs

            Static FilterSummary figures accept the following kwargs:

        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | Argument               | Description                                                                      | Example                        |
        +========================+==================================================================================+================================+
        | show : `bool`          | Whether or not to show the figure                                                | `show = True` (default)        |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | figsize : `tuple`      | The figure size                                                                  | `figsize = (10, 4)`            |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | subplots : `tuple`     | A tuple specifying the number of colums and rows (in that order) for the figure  | `subplots = (2, 3)`            |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | title : `str`          | The overall figure title                                                         | `title = "Today's run"         |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | xlabel : `str`         | The x-axis label of each subplot                                                 | `xlabel = "Conditions"`        |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | ylabel : `str`         | The y-axis label of each subplot                                                 | `ylabel = "My Ct values"`      |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | rot : `float`          | The rotation of x-axis labels                                                    | `rot = 0.3`                    |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | frame   : `bool`       | Show left and top spines of subplots (if True)                                   | `frame = False` (default)      |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | style : `str`          | A `seaborn` style to set. Check out available styles [1].                        | `style = "darkgrid"`           |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | color : `str or list`  | The color for the boxes.                                                         | `color = ["yellow", "green"]`  |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+
        | \*\*kwargs             | Any additional kwargs that can be passed to the `seaborn`'s `boxplot()`.         |                                |
        +------------------------+----------------------------------------------------------------------------------+--------------------------------+



    `"interactive"` Kwargs

            Interactive FilterSummary figures accept the following kwargs:

        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                      | Example                                       |
        +=============================+==================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                                | `show = True` (default)                       |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                         | `title = "Today's run"`                       |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                                 | `ylabel = "Raw Ct values"`                    |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                             | `height = 50`                                 |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                              | `width = 50`                                  |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                         | `padding = 0.2`                               |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The hoverinfo to display.                                                                                                                                        | `hoverinfo = "name"`                          |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | colors : `list`             | List of two colors for before and after filtering boxes.                                                                                                         | `colors = ["red", "green"]`                   |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                                 | `template = "plotly_dark"`                    |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                          | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default `"y+x+name"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                  | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Box()`.                                                                                       |                                               |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+

    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """

    def __init__(self, mode=None):
        self._setup_default_params(static=defaults.static_FilterSummary, interactive=defaults.interactive_FilterSummary)
        super().__init__(mode=mode)
        self._before = Results.Results()
        self._after = Results.Results()

    def clear(self):
        """
        Clears the pre- and post-filtering Ct value records.
        """
        self._before = Results.Results()
        self._after = Results.Results()

    def add_before(self, assay: Assay.Assay):
        """
        Add a pre-filtered set of Ct values
        from an `qpcr.Assay` object.

        Parameters
        ----------
        assay : qpcr.Assay
            An `qpcr.Assay` object
        """
        if isinstance(assay, list):
            return [self.add_before(i) for i in assay]
        # setup groups and stuff information
        if self._before.is_empty:
            self._before.setup_cols(assay)
            self._after.setup_cols(assay)

        self._before.add_Ct(assay)

    def add_after(self, assay: Assay.Assay):
        """
        Add a post-filtered set of Ct values
        from an `qpcr.Assay` object.

        Parameters
        ----------
        assay : qpcr.Assay
            An `qpcr.Assay` object
        """
        if isinstance(assay, list):
            return [self.add_after(i) for i in assay]
        self._after.add_Ct(assay)

    def _interactive_plot(self, **kwargs):
        """
        Generates an interactive Filter Summary Figure
        """
        kwargs = self.update_params(kwargs)

        # get the data
        before = self._before.get().copy()
        after = self._after.get().copy()

        title, show, assays, headers, ncols, nrows, xlabel, ylabel = self._prep_shared_kwargs(kwargs)
        hpad, vpad = self._get_padding(kwargs)
        height, width = self._get_height_and_width(kwargs)

        hoverinfo = kwargs.pop("hoverinfo", "y+x+name")
        template = kwargs.pop("template", None)

        # setup color scheme
        colors = kwargs.pop("colors", ["blue", "crimson"])
        before_color, after_color = colors

        fig, Coords = self._setup_interactive_figure(ncols, nrows, xlabel, ylabel, title, headers, hpad, vpad, height, width, template)

        fig.update_layout(
            boxmode="group",
        )

        show_legend = True
        # now iterate over each assay and make a BoxPlot
        for assay in assays:

            # get before and after filter datasets for each assay
            pre_Cts = before[["group", "group_name", assay]]
            post_Cts = after[["group", "group_name", assay]]

            row, col = Coords.get()

            # pre-filter boxes
            fig.add_trace(
                go.Box(
                    x=pre_Cts["group_name"],
                    y=pre_Cts[assay],
                    name="Before",
                    hoverinfo=hoverinfo,
                    marker_color=before_color,
                    legendgroup='group1',
                    showlegend=show_legend,
                    **kwargs,
                ),
                row,
                col,
            )
            # post-filter boxes
            fig.add_trace(
                go.Box(
                    x=post_Cts["group_name"],
                    y=post_Cts[assay],
                    name="After",
                    hoverinfo=hoverinfo,
                    marker_color=after_color,
                    legendgroup='group2',
                    showlegend=show_legend,
                    **kwargs,
                ),
                row,
                col,
            )

            # set show_legend to False after the first assay,
            # to avoid repetitive legends...
            show_legend = False

        fig.update_layout(showlegend=True)
        if show:
            fig.show()

        return fig

    def _static_plot(self, **kwargs):
        """
        Generates a static Filter Summary Fig
        """
        kwargs = self.update_params(kwargs)

        # get the data
        before = self._before.get().copy()
        after = self._after.get().copy()
        title, show, assays, headers, ncols, nrows, xlabel, ylabel = self._prep_shared_kwargs(kwargs)

        rot = kwargs.pop("rot", None)
        show_spines = kwargs.pop("frame", True)

        palette, is_palette = gx.generate_palette(kwargs, True)
        if palette:
            if is_palette:
                kwargs["palette"] = palette
            else:
                kwargs["color"] = palette
        else:
            kwargs["palette"] = defaults.default_palette

        style = kwargs.pop("style", defaults.default_style)
        sns.set_style(style)

        fig, Coords = self._setup_static_figure(ncols, nrows, title, kwargs)

        # set a kind column for grouping in the boxplots
        before["kind"] = "before"
        after["kind"] = "after"

        show_legend = True
        for assay, header in zip(assays, headers):

            # get before and after filter datasets for each assay
            pre_Cts = before[["group", "group_name", "kind", assay]]
            post_Cts = after[["group", "group_name", "kind", assay]]

            # assemble data to one dataframe
            df = pd.concat((pre_Cts, post_Cts))

            ax = Coords.subplot()

            # main box plot
            sns.boxplot(data=df, x="group_name", y=assay, hue="kind", ax=ax, **kwargs)

            # format the subplot
            ax.set(title=header, ylabel=ylabel, xlabel=xlabel)
            if rot is not None and rot != 0:
                plt.setp(ax.xaxis.get_majorticklabels(), rotation=-rot, ha="left", rotation_mode="anchor")
            ax.legend().remove()

        if show_legend:
            ax.legend(bbox_to_anchor=(1, 1), loc=None, frameon=False)

        if not show_spines:
            sns.despine()

        plt.tight_layout()
        if show:
            plt.show()

        return fig

    def _prep_shared_kwargs(self, kwargs):
        title, show = self._get_title_and_show(kwargs)
        assays = self._before.data_cols  # [ i for i in before.columns if i not in defaults.setup_cols ]
        headers = kwargs.pop("headers", assays)
        ncols, nrows = gx.make_layout_from_list(assays)
        xlabel, ylabel = self._axeslabels(kwargs)

        return title, show, assays, headers, ncols, nrows, xlabel, ylabel


if __name__ == "__main__":

    import qpcr

    # we set up the paths to 28S+actin as our normalisers
    normaliser_files = ["../Examples/Example Data/28S.csv", "../Examples/Example Data/actin.csv"]

    # we also set up the paths to the HNRNPL and SRSF11 transcripts
    assay_files = [
        "../Examples/Example Data/HNRNPL_nmd.csv",
        "../Examples/Example Data/HNRNPL_prot.csv",
        "../Examples/Example Data/SRSF11_nmd.csv",
        "../Examples/Example Data/SRSF11_prot.csv",
    ]

    a = qpcr.read(assay_files, replicates=(6, 5, 7, 6))

    # works :-)
    # p = FilterSummary( mode = "interactive" )
    # p.add_before(a)
    # p.add_after(a)
    # fig = p.plot()
    # fig.show()

    # works :-)
    p = FilterSummary()
    p.add_before(a)
    p.add_after(a)
    fig = p.plot(palette="viridis")
    plt.show()
