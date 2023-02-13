"""
Defines the `ReplicateBoxPlot` responsible for visualising the of raw Ct values from `qpcr.Assay` objects.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import plotly.graph_objs as go

import qpcr.defaults as defaults
from qpcr import _auxiliary as aux
import qpcr._auxiliary.graphical as gx

import qpcr.Plotters._base as base
import qpcr.main as main


class ReplicateBoxPlot(base.AssayPlotter):
    """
    Generate a boxplot figure summary for the input sample replicates.


    Note
    ------
    This used to be the core of making filter summary figures. However, this has now
    been replaced by a dedicated `FilterSummary` figure class. Hence, support for linking
    `qpcr.Filters` directly to this figure class will be dropped in a future release!

    Parameters
    ----------
    Filter : qpcr.Filters.Filter
        A qpcr.Filters.Filter object. The Filter can also be set using the `filter()` method.

    mode : str
        The plotting mode (either `"interactive"` or `"static"`).



    Plotting Kwargs
    ===============

    `"static"` Kwargs

        Static ReplicateBoxPlot figures accept the following kwargs:

        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | Argument               | Description                                                                      | Example                  |
        +========================+==================================================================================+==========================+
        | show : `bool`          | Whether or not to show the figure                                                | `show = True` (default)  |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | figsize : `tuple`      | The figure size                                                                  | `figsize = (10, 4)`      |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | subplots : `tuple`     | A tuple specifying the number of colums and rows (in that order) for the figure  | `subplots = (2, 3)`      |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | title : `str`          | The overall figure title                                                         | `title = "My assays"     |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | style : `str`          | A `seaborn` style to set. Check out available styles [1].                        | `style = "darkgrid"`     |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | color : `str or list`  | The fillcolor for the boxes.                                                     | `color = "yellow"`       |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | hue : `str`            | The data column by which to colo.                                                | `hue = "Ct"`             |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+
        | \*\*kwargs             | Any additional kwargs that can be passed to the `seaborn`'s `boxplot()`.         |                          |
        +------------------------+----------------------------------------------------------------------------------+--------------------------+


    `"interactive"` Kwargs

        Interactive ReplicateBoxPlot figures accept the following kwargs:

        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | Argument          | Description                                                                 | Example                    |
        +===================+=============================================================================+============================+
        | show : `bool`     | Whether or not to show the figure                                           | `show = True` (default)    |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | title : `str`     | The overall figure title                                                    | `title = "My run"`         |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | ylabel : `str`    | The y-axis title                                                            | `ylabel = "Raw Ct value"`  |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | height : `int`    | Height of the figure                                                        | `height = 50`              |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | width : `int`     | Width of the figure                                                         | `width = 50`               |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | template : `str`  | The `plotly` template to use. Check out available templates [2].            | `template = "plotly_dark"` |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+
        | \*\*kwargs        | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Box()`.  |                            |
        +-------------------+-----------------------------------------------------------------------------+----------------------------+


    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    """

    def __init__(self, mode: str = None):
        self._setup_default_params(static=defaults.static_ReplicateBoxPlot, interactive=defaults.interactive_ReplicateBoxPlot)
        # __init__ is going to require default_params to be already set!
        super().__init__(mode=mode)
        self._data = None

    def link(self, obj: main.Assay):
        """
        Links an Assay object to the BoxPlotter.
        This will simply add the Ct column to the current overall data!
        Hence, repeated linking of Assay objects will add data and NOT replace any existing.

        Parameters
        ----------
        obj : qpcr.Assay
            A qpcr.Assay object.
        """
        if self._data is None:
            super().link(obj)
            self._data[defaults.dataset_header] = self._obj.id()
            return
        data = self._data.copy()
        super().link(obj)

        # add itentifier column
        self._data[defaults.dataset_header] = self._obj.id()

        # and concat together
        self._data = pd.concat([data, self._data], ignore_index=True)

    def _interactive_plot(self, **kwargs):
        """
        Generates an interactive Boxplot summary of the input Ct values
        """
        kwargs = self.update_params(kwargs)
        data = self._data

        title, show = self._get_title_and_show(kwargs)
        xlabel, ylabel = self._axeslabels(kwargs)
        groups, group_names = self._get_groups_and_names(data)
        hpad, vpad, height, width, template, hoverinfo = self._get_interactive_setup_kwargs(kwargs)

        fig = go.Figure()
        fig.update_layout(
            title=title,
            boxmode="group",
            height=height,
            width=width,
            template=template,
        )

        # add default ylabel
        fig.update_yaxes(title_text=ylabel)
        fig.update_xaxes(showgrid=True, title_text=xlabel)

        for group, name in zip(groups, group_names):
            tmp_df = data.query(f"group == {group}")

            fig.add_trace(
                go.Box(
                    x=tmp_df["assay"],
                    y=tmp_df["Ct"],
                    name=name,
                    hoverinfo="y+name",
                    **kwargs,
                ),
            )

        if show:
            fig.show()

        return fig

    def _static_plot(self, **kwargs):
        """
        Generates a static Boxplot summary of the input Ct values
        """
        kwargs = self.update_params(kwargs)
        data = self._data

        title, show = self._get_title_and_show(kwargs)
        xlabel, ylabel = self._axeslabels(kwargs)

        style = kwargs.pop("style", defaults.default_style)
        sns.set_style(style)

        palette, is_palette = gx.generate_palette(kwargs, True)
        if is_palette:
            kwargs["palette"] = palette
        else:
            kwargs["color"] = palette

        hue = kwargs.pop("hue", "group_name")
        show_spines = kwargs.pop("frame", True)
        ncols, nrows = kwargs.pop("subplots", gx.make_layout(data, "assay"))

        fig, Coords = self._setup_static_figure(ncols, nrows, title, kwargs)

        # for assay in aux.sorted_set(data["assay"]):
        for assay, tmp in data.groupby(defaults.dataset_header):
            ax = Coords.subplot()

            # tmp = data.query(f"assay == '{assay}'")

            sns.boxplot(
                data=tmp,
                x=defaults.dataset_header,
                y=ylabel,
                hue=hue,
                # palette = palette,
                ax=ax,
                **kwargs,
            )

            ax.legend(bbox_to_anchor=(1, 1), loc=None).remove()
            ax.set(
                title=assay,
                xlabel=xlabel,
                xticklabels=[],
            )

            if not show_spines:
                sns.despine()

        # add one single legend to the last plot
        ax.legend(bbox_to_anchor=(1, 1), loc=None)

        plt.tight_layout()

        if show:
            fig.show()

        return fig


if __name__ == "__main__":

    import qpcr

    # we set up the paths to 28S+actin as our normalisers
    f = "../Examples/Example Data/28S.csv"

    a = qpcr.read(f)

    # works :-)
    p = ReplicateBoxPlot(mode="interactive")
    p.link(a)
    fig = p.plot()
    fig.show()

    # works :-)
    p = ReplicateBoxPlot(mode="static")
    p.link(a)
    fig = p.plot()
    plt.show()
