"""
This is the `EfficiencyLines` FigureClass that visualises the regression lines performed by a ``qpcr.Calibrator`` to compute new amplification efficiencies.
"""

from matplotlib.lines import Line2D
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go

import qpcr.defaults as defaults
from qpcr import _auxiliary as aux
import qpcr._auxiliary.graphical as gx
import qpcr.Plotters._base as base
import qpcr.main as main


class EfficiencyLines(base.Plotter):
    """
    Generates a Figure for the linear regressions used for Assay efficiency
    calculations. This FigureClass specifically works with the `qpcr.Calibrator`
    class.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).



    Plotting Kwargs
    ================

    `"static"` Kwargs

            Static EfficiencyLines figures accept the following kwargs:

        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                  | Description                                                                                    | Example                                       |
        +===========================+================================================================================================+===============================================+
        | show : `bool`             | Whether or not to show the figure                                                              | `show = True` (default)                       |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | figsize : `tuple`         | The figure size                                                                                | `figsize = (10, 4)`                           |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`             | The overall figure title                                                                       | `title = "Today's efficiencies"               |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`            | The x-axis label of each subplot                                                               | `xlabel = "Log Dilution"`                     |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`            | The y-axis label of each subplot                                                               | `ylabel = "My Ct Values"`                     |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | rot : `float`             | The rotation of x-axis labels                                                                  | `rot = 0.3`                                   |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`          | A list of titles for each subplot in the preview figure                                        | `headers = ["transcript A", "transcript B"]`  |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | label_subplots  : `bool`  | Add each subplot with A, B, C ... (if True, default)                                           | `label_subplots = True` (default)             |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | labeltype : `str`         | The starting character for subplot labelling. By default an `"A"`.                             | `labeltype = "a"`                             |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | frame   : `bool`          | Show left and top spines of subplots (if True)                                                 | `frame = False` (default)                     |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | color : `str or list`     | The fillcolor for the individual dots                                                          | `color = "yellow"`                            |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | style : `str`             | A `seaborn` style to set. Check out available styles [1].                                      | `style = "darkgrid"`                          |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgecolor : `str or list` | The edgecolor for the individual dots for the datapoints.                                      | `edgecolor = "black"`                         |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | edgewidth : `float`       | The width of the edge of individual dots.                                                      | `edgewidth = 0.5`                             |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | linecolor : `str or list` | The color for regression line.                                                                 | `linecolor = "crimson"`                       |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | linewidth : `float`       | The width of the regression line.                                                              | `edgewidth = 0.5`                             |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\kwargs                 | Any additional kwargs that can be passed to `seaborn`'s `scatterplot` and `lineplot` (both!).  |                                               |
        +---------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------+



    `"interactive"` Kwargs

            Interactive EfficiencyLines figures accept the following kwargs:

        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | Argument                    | Description                                                                                                                                                      | Example                                       |
        +=============================+==================================================================================================================================================================+===============================================+
        | show : `bool`               | Whether or not to show the figure                                                                                                                                | `show = True` (default)                       |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | title : `str`               | The overall figure title                                                                                                                                         | `title = "Today's efficiencies"`              |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | xlabel : `str`              | The x axis label                                                                                                                                                 | `xlabel = "Log Dilution"`                     |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | ylabel : `str`              | The y axis label                                                                                                                                                 | `ylabel = My super Ct values"`                |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | height : `int`              | Height of the figure                                                                                                                                             | `height = 50`                                 |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | width : `int`               | Width of the figure                                                                                                                                              | `width = 50`                                  |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | padding : `float or tuple`  | Padding between subplots. This can be a single float (interpreted as horizontal padding), or a tuple of (horizontal, vertical) paddings.                         | `padding = 0.2`                               |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | template : `str`            | The `plotly` template to use. Check out available templates [2].                                                                                                 | `template = "plotly_dark"`                    |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | headers : `list`            | A list of titles for each subplot in the preview figure                                                                                                          | `headers = ["transcript A", "transcript B"]`  |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | hoverinfo : `str`           | The type of hoverinfo to display. By default just `"y+x"`. Learn more about plotly hoverinfo [3]. Please, note that `hovertemplate` is not currently supported.  | `hoverinfo = "name+y"`                        |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
        | \*\kwargs                   | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Scatter()`.                                                                                   |                                               |
        +-----------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+



    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """

    def __init__(self, mode: str = None):
        self._setup_default_params(static=defaults.static_EfficiencyLines, interactive=defaults.interactive_EfficiencyLines)
        super().__init__(mode=mode)
        self._Calibrator = None

    def link(self, calibrator: main.Calibrator):
        """
        Links a `qpcr.Calibrator` object to source data from.

        Note
        -----
        Only a `qpcr.Calibrator` that has actually *de novo*
        computed efficiencies will have data to plot, and only for
        the newly computed effiencies!

        Parameters
        ----------
        calibrator : qpcr.Calibrator
            A `qpcr.Calibrator` object that has computed new efficiencies.
        """
        self._Calibrator = calibrator

    def _prep_shared_kwargs(self, kwargs):
        data = self._Calibrator._computed_values
        headers = list(data.keys())
        headers = kwargs.pop("headers", headers)
        data = {h: i for h, i in zip(headers, list(data.values()))}

        title, show = self._get_title_and_show(kwargs)
        xlabel, ylabel = self._axeslabels(kwargs)

        # NOTE: interactive used to have this line transposed (i.e. nrows, ncols = ...)
        ncols, nrows = gx.make_layout_from_list(headers)

        return data, headers, title, show, xlabel, ylabel, ncols, nrows

    def _static_plot(self, **kwargs):
        """
        Generates a static Effiency Regression Figure
        """
        kwargs = self.update_params(kwargs)

        data, headers, title, show, xlabel, ylabel, ncols, nrows = self._prep_shared_kwargs(kwargs)
        label_subplots, start_character, show_spines, rot = self._get_labels_and_spines_and_rot(kwargs)

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

        edgecolor = kwargs.pop("edgecolor", "white")
        edgewidth = kwargs.pop("edgewidth", None)

        linecolor = kwargs.pop("linecolor", "black")
        linewidth = kwargs.pop("linewidth", 1)

        fig, Coords = self._setup_static_figure(ncols, nrows, title, kwargs)

        idx = 0
        for id, obj in data.items():

            ax = Coords.subplot()

            # get data to plot
            dilutions, cts = obj.values()
            eff = obj.efficiency()
            model = obj.model()
            Rsquare = model.rvalue**2

            # plot the regression line
            yvals = model.slope * dilutions + model.intercept
            sns.lineplot(x=dilutions, y=yvals, color=linecolor, linewidth=linewidth, ax=ax, **kwargs)

            # plot the original data
            sns.scatterplot(x=dilutions, y=cts, hue=dilutions if is_palette else None, edgecolor=edgecolor, linewidth=edgewidth, ax=ax, **palette, **kwargs)

            # add additional info to a legend
            ax.legend(
                handles=[
                    # R^2 Value
                    Line2D([0], [0], color="black", visible=False, label=f"$R^2$ \t = {Rsquare:.4f}"),
                    # Efficiency Value
                    Line2D([0], [0], color="black", visible=False, label=f"$eff.$\t= {eff:.4f}"),
                ],
                loc="upper right",
                frameon=False,
            )

            # some formatting...
            ax.set(
                title=id,
                xlabel=xlabel,
                ylabel=ylabel,
            )

            if rot is not None:
                self._set_xtick_rotation(ax, rot)

            # add ABCD... label to subplot
            if label_subplots:
                self._add_subplot_label(idx, ax, start_character)

            if not show_spines:
                sns.despine()

            idx += 1

        plt.tight_layout()
        if show:
            plt.show()

        return fig

    def _interactive_plot(self, **kwargs):
        """
        Generates an interactive EfficiencyLines figure
        """
        kwargs = self.update_params(kwargs)
        data, headers, title, show, xlabel, ylabel, ncols, nrows = self._prep_shared_kwargs(kwargs)
        hpad, vpad = self._get_padding(kwargs)
        height, width = self._get_height_and_width(kwargs)
        template = kwargs.pop("template", "plotly_white")
        hoverinfo = kwargs.pop("hoverinfo", "y+x")

        fig, Coords = self._setup_interactive_figure(ncols, nrows, xlabel, ylabel, title, headers, hpad, vpad, height, width, template)

        # iterate over the data...
        idx = 0
        for id, obj in data.items():

            row, col = Coords.get()

            # get data to plot
            dilutions, cts = obj.values()
            eff = obj.efficiency()
            model = obj.model()
            Rsquare = model.rvalue**2

            # plot the original data
            fig.add_trace(go.Scatter(name=id, x=dilutions, y=cts, mode="markers", showlegend=False, hoverinfo=hoverinfo, **kwargs), row, col)

            # plot the regression line
            yvals = model.slope * dilutions + model.intercept

            fig.add_trace(go.Scatter(name=id, x=dilutions, y=yvals, mode="lines", showlegend=False, **kwargs), row, col)

            # add infos as annotations (because custom legends don't work here)
            # NOTE: The xref x{idx} is a pretty nice lifehack to essentially
            #       emulate the ax behaviour from matplotib.subplots...
            xref = "x" if idx == 0 else f"x{idx+1}"
            yref = "y" if idx == 0 else f"y{idx+1}"
            info_text = f"R^2 = {Rsquare:.4f}<br>eff. = {eff:.4f}"

            fig.add_annotation(dict(x=0, y=np.max(cts), xref=xref, yref=yref, text=info_text, showarrow=False))
            idx += 1

        if show:
            fig.show()

        return fig


if __name__ == '__main__':

    import qpcr

    file = "../Examples/Example Data/Calibration/calibrators_only.csv"

    a = qpcr.read(file)

    c = qpcr.Calibrator()
    c.dilution(2)
    a = c.calibrate(a, remove_calibrators=False)
    for i in range(5):
        a.id(f"{i}")
        a = c.calibrate(a)

    # works :-)
    p = EfficiencyLines(mode="static")
    p.link(c)
    fig = p.plot()
    plt.show()

    # works :-)
    p = EfficiencyLines(mode="interactive")
    p.link(c)
    fig = p.plot()
    fig.show()
