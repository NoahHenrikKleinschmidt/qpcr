"""
Defines the `AssayBars` and `AssayDots` class, which is used to preview the results of a qpcr.Results object, visualising the different assays in different subplots.
"""

from qpcr import defaults
from qpcr import _auxiliary as aux
import qpcr._auxiliary.graphical as gx
import qpcr.Plotters._base as base

import matplotlib.pyplot as plt
import seaborn as sns

import plotly.graph_objs as go

import numpy as np

logger = aux.default_logger()


class AssaySubplotsResults(base.ResultsPlotter):
    """
    A super class for Figureclasses that use a Results object and operate assay-wise
    (different assays in different subplots).
    """

    def __init__(self, mode: str):
        super().__init__(mode)


class AssayBars(AssaySubplotsResults):
    """
    Generate a Preview of all results from all Assays in subplots as Bar Charts.
    Each assay will be shown in a separate subplot.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).


    Plotting Kwargs
    ===============

    `"static"` Kwargs

        Static AssayBars figures accept the following kwargs:

        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                  | Description                                                                                                            | Example                                       |
        +===========================+========================================================================================================================+===============================================+
        | show: `bool`              | Whether or not to show the figure                                                                                      | `show = True` (default)                       |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | figsize: `tuple`          | The figure size                                                                                                        | `figsize = (10, 4)`                           |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title: `str`              | The overall figure title                                                                                               | `title = "Today's results"                    |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel: `str`             | The x-axis label of each subplot                                                                                       | `xlabel = "Conditions"`                       |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel: `str`             | The y-axis label of each subplot                                                                                       | `ylabel = "Mean Fold Change"`                 |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | rot: `float`              | The rotation of x-axis labels                                                                                          | `rot = 0.3`                                   |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers: `list`           | A list of titles for each subplot in the preview figure                                                                | `headers = ["transcript A", "transcript B"]`  |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | label_subplots: `bool`    | Label each subplot with A, B, C ... (if True, default)                                                                 | `label_subplots = True` (default)             |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | labeltype: `str`          | The starting character for subplot labelling. By default an `"A"`.                                                     | `labeltype = "a"`                             |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | annotate_pvals: `bool`    | Add annotations of p-value significance to each subplot ( True by default if statistics are available)                 | `annotate_pvals = True` (default)             |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | pval_kws: `dict`          | Keywords for pvalue formatting using `encode_pvalues`                                                                  | `pval_kws = dict(style = "*")`                |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | frame: `bool`             | Show left and top spines of subplots (if True)                                                                         | `frame = False` (default)                     |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | color: `str or list`      | The fillcolor for the individual bars                                                                                  | `color = "yellow"`                            |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | style: `str`              | A `seaborn` style to set. Check out available styles [1].                                                              | `style = "darkgrid"`                          |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgecolor: `str or list`  | The edgecolor for the individual bars                                                                                  | `edgecolor = "black"`                         |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgewidth: `float`        | The width of the edge of individual bars                                                                               | `edgewidth = 0.5`                             |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ecolor: `str`             | The color of errorbars                                                                                                 | `ecolor = "orange"`                           |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                | Any additional kwargs that can be passed to the `matplotlib`-backend pandas `.plot.bar()` API.                         |                                               |
        +---------------------------+------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+


    `"interactive"` Kwargs

        Interactive AssayBars figures accept the following kwargs:

        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                    | Example                                       |
        +=============================+================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                              | `show = True` (default)                       |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                       | `title = "Today's results"`                   |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`              | The x axis label                                                                                                                                               | `xlabel = "My super qPCR samples"`            |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                               | `ylabel = "Mean of ddCt"`                     |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                           | `height = 50`                                 |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                            | `width = 50`                                  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                       | `padding = 0.2`                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                               | `template = "plotly_dark"`                    |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                        | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | legend_title : `str`        | The title to be displayed above the legend                                                                                                                     | `legend_title = "my assays"`                  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default just `"y"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                  | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Bar()`.                                                                                     |                                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+


    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """

    def __init__(self, mode: str = None):
        self._setup_default_params(static=defaults.static_PreviewBars, interactive=defaults.interactive_PreviewBars)
        # __init__ is going to require default_params to be already set!
        super().__init__(mode)

    def _static_plot(self, **kwargs):
        """
        The static Preview Results Figure
        """
        kwargs = self.update_params(kwargs)
        data = self.get()

        ref_col, ncols, nrows, headings, x, y, sterr, query, xlabel, ylabel = self._prep_properties(kwargs)
        label_subplots, start_character, show_spines, rot = self._get_labels_and_spines_and_rot(kwargs)
        title, show = self._get_title_and_show(kwargs)

        headers = kwargs.pop("headers", None)
        legend = kwargs.pop("legend", True)

        style = kwargs.pop("style", defaults.default_style)
        sns.set_style(style)
        palette, is_palette = gx.generate_palette(kwargs, True)
        if palette:
            if is_palette:
                kwargs["palette"] = palette
            else:
                kwargs["color"] = palette
        else:
            kwargs["palette"] = defaults.default_palette

        edgecolor = kwargs.pop("edgecolor", "white")
        edgewidth = kwargs.pop("edgewidth", None)
        ecolor = kwargs.get("ecolor", "black")

        # get kwargs for pvalue annotations
        annotate_pvals = kwargs.pop("annotate_pvals", self._obj.comparisons is not None)
        pval_kws = kwargs.pop("pval_kws", {})

        if annotate_pvals:
            if self._obj.comparisons is None:
                raise AttributeError("Cannot annotate pvalues if no comparisons have been made! Please first perform a statistical assaywise comparison.")

        if edgewidth is None:
            edgewidth = kwargs.pop("linewidth", 1)
        else:
            kwargs.pop("linewidth", None)

        fig, Coords = self._setup_static_figure(ncols, nrows, title, kwargs)

        idx = 0
        for assay, tmp_df in data.groupby(ref_col):

            # now plot a new bar chart
            subplot = Coords.subplot()

            sns.barplot(data=tmp_df, x=x, y=y, ax=subplot, edgecolor=edgecolor, linewidth=edgewidth, **kwargs)

            if legend:
                subplot.legend(title=kwargs.get("legend_title", None))

            subplot.errorbar(
                x=tmp_df[x],
                y=tmp_df[y],
                yerr=tmp_df[sterr],
                fmt=".",
                markersize=0,
                capsize=3,
                ecolor=ecolor,
            )

            subplot.set(
                title=assay if headers is None else headers[idx],
                xlabel=xlabel,
                ylabel=ylabel,
            )

            if rot is not None:
                self._set_xtick_rotation(subplot, rot)

            if not show_spines:
                sns.despine()

            # add ABCD... label to subplot
            if label_subplots:
                self._add_subplot_label(idx, subplot, start_character)

            # add stats annotations
            if annotate_pvals:
                if assay not in self._obj.comparisons:
                    logger.info(f"Could not find a comparison for {assay}. Perhaps you did not use an assaywise-comparison but a groupwise-comparison?")
                else:
                    comparison = self._obj.comparisons[assay]
                    self._annotate_pvalues(x, y, sterr, pval_kws, comparison, tmp_df, subplot)

            idx += 1

        plt.tight_layout()

        if show:
            plt.show()

        # print("<<< Static Check >>>")

        return fig

    def _interactive_plot(self, **kwargs):
        """
        The interactive Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self.get()

        # setup figure framework variables that have to be removed from
        ref_col, ncols, nrows, headings, x, y, sterr, query, xlabel, ylabel = self._prep_properties(kwargs)
        title, show = self._get_title_and_show(kwargs)
        headers = kwargs.pop("headers", headings)

        hpad, vpad, height, width, template, hoverinfo = self._get_interactive_setup_kwargs(kwargs)

        fig, Coords = self._setup_interactive_figure(ncols, nrows, xlabel, ylabel, title, headers, hpad, vpad, height, width, template)

        fig.update_layout(
            legend={
                "title": kwargs.pop("legend_title", ref_col),
            },
        )

        idx = 0
        # for assay in headings:
        for assay, tmp_df in data.groupby(ref_col):
            row, col = Coords.get()
            # tmp_df = data.query(query.format(ref_col = ref_col, q = assay))
            # tmp_df = tmp_df.sort_values("group")

            # now plot a new bar chart
            fig.add_trace(go.Bar(name=assay, y=tmp_df[y], x=tmp_df[x], error_y=dict(type='data', array=tmp_df[sterr]), hoverinfo=hoverinfo, **kwargs), row, col)

            idx += 1

        if show:
            fig.show()

        # print("<<< Interactive Check >>>")
        return fig


class AssayDots(AssaySubplotsResults):
    """
    Generate a Preview of all results from all Assays in subplots, plotting individual values
    in a Dot Plot rather than Bar Plot. Each assay will be shown in a separate subplot.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).


    Plotting Kwargs
    ================

    `"static"` Kwargs

            Static AssayDots figures accept the following kwargs:

        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        |         Argument         |                                              Description                                               |                   Example                    |
        +==========================+========================================================================================================+==============================================+
        | show : `bool`            | Whether or not to show the figure                                                                      | `show = True` (default)                      |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | figsize : `tuple`        | The figure size                                                                                        | `figsize = (10, 4)`                          |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | title : `str`            | The overall figure title                                                                               | `title = "Today's results"                   |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | xlabel : `str`           | The x-axis label of each subplot                                                                       | `xlabel = "Conditions"`                      |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | ylabel : `str`           | The y-axis label of each subplot                                                                       | `ylabel = "Mean Fold Change"`                |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | rot : `float`            | The rotation of x-axis labels                                                                          | `rot = 0.3`                                  |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | headers : `list`         | A list of titles for each subplot in the preview figure                                                | `headers = ["transcript A", "transcript B"]` |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | label_subplots  : `bool` | Add each subplot with A, B, C ... (if True, default)                                                   | `label_subplots = True` (default)            |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | labeltype : `str`        | The starting character for subplot labelling. By default an `"A"`.                                     | `labeltype = "a"`                            |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | annotate_pvals: `bool`   | Add annotations of p-value significance to each subplot ( True by default if statistics are available) | `annotate_pvals = True` (default)            |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | pval_kws: `dict`         | Keywords for pvalue formatting using `encode_pvalues`                                                  | `pval_kws = dict(style = "*")`               |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | frame   : `bool`         | Show left and top spines of subplots (if True)                                                         | `frame = False` (default)                    |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | violin   : `bool`        | Show symmetric kde of the dots of each group (if True).                                                | `violin = False`                             |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | color : `str or list`    | The fillcolor for the individual dots from replicate groups                                            | `color = "yellow"`                           |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | style : `str`            | A `seaborn` style to set. Check out available styles [1].                                              | `style = "darkgrid"`                         |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | \*\*kwargs               | Any additional kwargs that can be passed to the `seaborn`'s `stripplot`.                               |                                              |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+



    `"interactive"` Kwargs

            Interactive AssayDots figures accept the following kwargs:

        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                    | Example                                       |
        +=============================+================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                              | `show = True` (default)                       |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                       | `title = "Today's results"`                   |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`              | The x axis label                                                                                                                                               | `xlabel = "My super qPCR samples"`            |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                               | `ylabel = "Mean of ddCt"`                     |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                           | `height = 50`                                 |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                            | `width = 50`                                  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                       | `padding = 0.2`                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                               | `template = "plotly_dark"`                    |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                        | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | violin   : `bool`           | Show symmetric kde of the dots of each group (if True).                                                                                                        | `violin = False`                              |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default just `"y"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\*kwargs                  | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Violin()`.                                                                                  |                                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+

    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """

    def __init__(self, mode: str = None):
        self._setup_default_params(static=defaults.static_PreviewDots, interactive=defaults.interactive_PreviewDots)
        # __init__ is going to require default_params to be already set!
        super().__init__(mode)

    def _static_plot(self, **kwargs):
        """
        The static Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self._rep_data

        # get assays to plot
        assays, ncols, nrows = self._make_assays_and_subolot_layout(data)

        # get kwargs incompatible with the main plotting method
        title, show = self._get_title_and_show(kwargs)
        label_subplots, start_character, show_spines, rot = self._get_labels_and_spines_and_rot(kwargs)
        xlabel, ylabel = self._axeslabels(kwargs)

        headers = kwargs.pop("headers", assays)
        show_violins = kwargs.pop("violin", True)
        alpha = kwargs.pop("alpha", 1)

        # generate a custom color palette in case color kwarg is provided
        palette, is_palette = gx.generate_palette(kwargs, True)
        if palette:
            if is_palette:
                palette = dict(palette=palette)
            else:
                palette = dict(color=palette)
        else:
            palette = dict(palette=defaults.default_palette)

        # set a seaborn style
        style = kwargs.pop("style", defaults.default_style)
        sns.set_style(style)

        # get kwargs for pvalue annotations
        annotate_pvals = kwargs.pop("annotate_pvals", self._obj.comparisons is not None)
        pval_kws = kwargs.pop("pval_kws", {})

        if annotate_pvals:
            if self._obj.comparisons is None:
                raise AttributeError("Cannot annotate pvalues if no comparisons have been made!, First perform a statistical assaywise comparison.")

        # make figure
        fig, Coords = self._setup_static_figure(ncols, nrows, title, kwargs)

        idx = 0
        for assay in assays:

            tmp_df = data[["group", "group_name", assay]]
            tmp_df = tmp_df.sort_values("group")

            # now plot a new violin chart
            subplot = Coords.subplot()

            if show_violins:
                sns.violinplot(
                    x=tmp_df["group_name"],
                    y=tmp_df[assay],
                    color=None,
                    inner=None,
                    ax=subplot,
                    **palette,
                )
                for i in subplot.collections:
                    i.set_alpha(alpha * 0.3)

            sns.stripplot(x=tmp_df["group_name"], y=tmp_df[assay], alpha=alpha, ax=subplot, **palette, **kwargs)

            subplot.set(
                title=assay if headers is None else headers[idx],
                xlabel=xlabel,
                ylabel=ylabel,
            )

            # adjust xtick rotation
            if rot is not None:
                self._set_xtick_rotation(subplot, rot)

            if not show_spines:
                sns.despine()

            # add ABCD... label to subplot
            if label_subplots:
                self._add_subplot_label(idx, subplot, start_character)

            # add stats annotations
            if annotate_pvals:
                if assay not in self._obj.comparisons:
                    logger.warning(f"Could not find a comparison for {assay}. Perhaps you did not use an assaywise-comparison but a groupwise-comparison?")
                else:
                    comparison = self._obj.comparisons[assay]
                    self._annotate_pvalues("group_name", assay, None, pval_kws, comparison, tmp_df, subplot)

            idx += 1

        plt.tight_layout()

        if show:
            plt.show()

        # print("<<< Static Check >>>")

        return fig

    def _interactive_plot(self, **kwargs):
        """
        The interactive Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self._rep_data

        assays, ncols, nrows = self._make_assays_and_subolot_layout(data)

        # get the info of groups present in the data
        # for later use as xtick labels
        group_names = data[["group", "group_name"]]
        group_names = group_names.sort_values("group")
        group_names = aux.sorted_set(group_names["group_name"])
        ticks = np.arange(len(group_names))

        # get incompaltible kwargs
        headers = kwargs.pop("headers", assays)
        show = kwargs.pop("show", True)
        show_violins = kwargs.pop("violin", True)

        title, show = self._get_title_and_show(kwargs)
        hpad, vpad, height, width, template, hoverinfo = self._get_interactive_setup_kwargs(kwargs)
        xlabel, ylabel = self._axeslabels(kwargs)

        # make figure
        fig, Coords = self._setup_interactive_figure(ncols, nrows, xlabel, ylabel, title, headers, hpad, vpad, height, width, template)

        idx = 0
        for assay, header in zip(assays, headers):

            row, col = Coords.get()

            tmp_df = data[["group", "group_name", assay]]
            tmp_df = tmp_df.sort_values("group")

            # now plot a new violin chart
            fig.add_trace(go.Violin(name=header, y=tmp_df[assay], x=tmp_df["group"], points="all", pointpos=0, hoverinfo=kwargs.pop("hoverinfo", "y"), **kwargs), row, col)

            # remove violins if not desired (leaving only dot plot)
            if not show_violins:
                fig.update_traces(
                    fillcolor="rgba(0,0,0,0)",
                    line_width=0,
                    selector=dict(type='violin'),
                )

            # update x axis to categorical group names
            fig.update_layout(
                {
                    f"xaxis{idx+1}": dict(
                        tickmode='array',
                        tickvals=ticks,
                        ticktext=group_names,
                    )
                }
            )
            idx += 1

        if show:
            fig.show()

        # print("<<< Interactive Check >>>")
        return fig

    def _make_assays_and_subolot_layout(self, data):
        assays = [i for i in data.columns if i not in defaults.setup_cols]
        # generate subplot layout
        ncols, nrows = gx.make_layout_from_list(assays)
        return assays, ncols, nrows


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
    n = qpcr.read(normaliser_files, replicates=(6, 5, 7, 6))

    a, n = qpcr.delta_ct([a, n])
    r = qpcr.normalise(a, n)
    r.drop_rel()

    # works :-)
    p = AssayBars(mode="interactive")
    p.link(r)
    fig = p.plot()
    fig.show()

    # works :-)
    p = AssayDots(mode="static")
    p.link(r)
    fig = p.plot()
    plt.show()
