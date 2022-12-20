"""
Defines the `PairwisePlot` class, which is used to preview the results of a pairwise-comparison.
"""

import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import qpcr._auxiliary.graphical as gx
import qpcr.Plotters._base as base
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D


logger = aux.default_logger()


class PairwisePlot(base.Plotter):
    """
    Generate a summary figure for a ``PairwiseComparison`` object.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).

    """

    def __init__(self, mode: str = None):
        self._setup_default_params(static=defaults.static_Pairwise, interactive=defaults.interactive_Pairwise)
        super().__init__(mode)

    def link(self, obj: ("ComparisonsCollection" or "PairwiseComparison")):
        """
        Link the PairwisePlot to a ``PairwiseComparison`` object or a ``ComparisonsCollection`` storing multiple ``PairwiseComparison``s.
        """
        self.obj = obj
        return self

    def get(self):
        """
        Get the data for the plot.
        """
        data = self.obj.get()
        return data

    def _static_plot(self, **kwargs):
        """
        The static Pairwise Plot
        """
        style = kwargs.get("style", None)
        sns.set_style(style)
        palette, is_palette = gx.generate_palette(kwargs, True)
        if palette:
            if is_palette:
                kwargs["palette"] = palette
            else:
                kwargs["color"] = palette
        else:
            kwargs["palette"] = defaults.default_palette

        if aux.pseudo_isinstance(self.obj, "ComparisonsCollection"):
            ncols, nrows = gx.make_layout_from_list(self.obj)
            fig = plt.figure(figsize=kwargs.pop("figsize", (ncols * 5, nrows * 5)))
            fig.set_tight_layout(False)
            subfigs = fig.subfigures(
                nrows=nrows,
                ncols=ncols,
                wspace=0.5,
            )

            for data, f in zip(self.obj, subfigs.flatten()):
                ax = f.add_subplot(111)
                self._static_single(data, ax, **kwargs)

            return fig

        elif aux.pseudo_isinstance(self.obj, "PairwiseComparison"):
            ax = kwargs.pop("ax", None)
            if ax is None:
                fig, ax = plt.subplots(figsize=kwargs.pop("figsize", None))
                fig.set_tight_layout(False)

            self._static_single(self.obj, ax, **kwargs)
            return ax.get_figure()
        else:
            raise TypeError("Cannot plot object of type '{}'".format(type(self.obj)))

    def _static_single(self, data, ax, **kwargs):
        """
        A single static Pairwise Plot
        """
        kwargs = self.update_params(kwargs)

        fig = ax.get_figure()
        title, data = data.id(), data.stack()

        style = kwargs.pop("style", "")
        grid = kwargs.pop("grid", "grid" in style)
        ax.grid(grid)
        if grid:
            ax.set_axisbelow(True)

        xlabel, ylabel = kwargs.pop("xlabel", ""), kwargs.pop("ylabel", "")
        title = kwargs.pop("title", title)

        x, y, hue, size = kwargs.pop("x", "b"), kwargs.pop("y", "a"), kwargs.pop("hue", "pval_adj"), kwargs.pop("size", "effect_size")
        sizes = kwargs.pop("sizes", (20, 200))
        cbar_title = kwargs.pop("colorbar_title", hue)
        size_title = kwargs.pop("sizelegend_title", size)
        cbar_loc = kwargs.pop("colorbar_loc", [0.99, 0.1, 0.018, 0.3])

        sns.scatterplot(data=data, x=x, y=y, hue=hue, size=size, sizes=sizes, legend=False, ax=ax, **kwargs)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        ax.xaxis.set_tick_params(rotation=90)

        # add the colorbar legend
        cmap = mpl.cm.get_cmap(kwargs.get("palette", defaults.default_palette))
        vmin = data[hue].min()
        vmax = data[hue].max()
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        cax = fig.add_axes(cbar_loc)
        cax.grid(False)
        cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation="vertical", label=cbar_title)

        # add the size legend
        summary = data[size].describe()

        l0 = [0]
        handles = [
            # Line2D( [0], [0], marker = "o", color = "w", label = f"{summary['min']:.2f}", markerfacecolor = "k", markersize = summary["min"] ),
            Line2D(l0, l0, marker="o", color="w", label=f"{summary['25%']:.2f}", markerfacecolor="k", markersize=summary["25%"]),
            Line2D(l0, l0, marker="o", color="w", label=f"{summary['50%']:.2f}", markerfacecolor="k", markersize=summary["50%"]),
            Line2D(l0, l0, marker="o", color="w", label=f"{summary['75%']:.2f}", markerfacecolor="k", markersize=summary["75%"]),
            # Line2D( [0], [0], marker = "o", color = "w", label = f"{summary['max']:.2f}", markerfacecolor = "k", markersize = summary["max"] )
        ]

        # for vertical
        cax.legend(handles=handles, title=size_title, loc="lower center", bbox_to_anchor=(1.1, 1.1), frameon=False)

        return fig

    def _interactive_plot(self, **kwargs):
        raise NotImplementedError("So far the PairwisePlot is only available through matplotlib")


if __name__ == "__main__":

    import qpcr

    assays = [
        '/Users/noahhk/GIT/qpcr/Examples/Example Data/HNRNPL_nmd.csv',
        '/Users/noahhk/GIT/qpcr/Examples/Example Data/HNRNPL_prot.csv',
        '/Users/noahhk/GIT/qpcr/Examples/Example Data/SRSF11_nmd.csv',
        '/Users/noahhk/GIT/qpcr/Examples/Example Data/SRSF11_prot.csv',
    ]
    normalisers = ['/Users/noahhk/GIT/qpcr/Examples/Example Data/28S.csv', '/Users/noahhk/GIT/qpcr/Examples/Example Data/actin.csv']

    assays = qpcr.read(assays, replicates=3)
    normalisers = qpcr.read(normalisers, replicates=6)

    assays = qpcr.delta_ct(assays)
    normalisers = qpcr.delta_ct(normalisers)

    results = qpcr.normalise(assays, normalisers)
    results.drop_rel()

    comparisons = qpcr.stats.groupwise_ttests(results)

    plotter = PairwisePlot()
    plotter.link(comparisons)
    fig = plotter.plot(figsize=(12, 12), s=500, palette="viridis")
    # plt.tight_layout(w_pad=80)
    plt.savefig("test.png", bbox_inches="tight")
    # plt.show()
