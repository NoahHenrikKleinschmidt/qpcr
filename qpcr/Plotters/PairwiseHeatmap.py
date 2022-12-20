"""
Defines the `PairwiseHeatmap` class, which is used to preview the results of a pairwise-comparison.
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


class PairwiseHeatmap(base.Plotter):
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
        if "cmap" not in kwargs:
            palette, is_palette = gx.generate_palette(kwargs, True)
            if palette:
                if is_palette:
                    kwargs["cmap"] = palette
            else:
                kwargs["cmap"] = defaults.default_palette

        if aux.pseudo_isinstance(self.obj, "ComparisonsCollection"):
            ncols, nrows = gx.make_layout_from_list(self.obj)
            fig, axs = plt.subplots(nrows, ncols, figsize=kwargs.pop("figsize", (ncols * 5, nrows * 5)))
            for comparison, ax in zip(self.obj, axs.flat):
                data = comparison.get()
                ax.set(
                    title=kwargs.pop("title", comparison.id()),
                    xlabel=kwargs.pop("xlabel", ""),
                    ylabel=kwargs.pop("ylabel", ""),
                )
                sns.heatmap(data, ax=ax, **kwargs)
                ax.set_title(comparison.id())
            return fig

        elif aux.pseudo_isinstance(self.obj, "PairwiseComparison"):
            ax = kwargs.pop("ax", None)
            if ax is None:
                fig, ax = plt.subplots(figsize=kwargs.pop("figsize", None))

            kwargs.setdefault("cbar_kws", {})
            kwargs["cbar_kws"].setdefault("label", "p-value")

            data = self.obj.get()
            ax.set(
                title=kwargs.pop("title", self.obj.id()),
                xlabel=kwargs.pop("xlabel", ""),
                ylabel=kwargs.pop("ylabel", ""),
            )
            sns.heatmap(data, ax=ax, **kwargs)

            return ax.figure
        else:
            raise TypeError("The object must be a PairwiseComparison or a ComparisonsCollection.")

    def _interactive_plot(self, **kwargs):
        """
        The interactive Pairwise Plot
        """
        raise NotImplementedError("Interactive Pairwise Plots are not yet implemented.")


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

    comparisons = qpcr.stats.assaywise_ttests(results)
    comparisons.make_symmetric()

    # plotter = PairwiseHeatmap()
    fig = qpcr.plot(comparisons[0])
    # fig = plotter.plot(figsize=(12, 12), cmap="mako_r")
    # plt.tight_layout(w_pad=80)
    plt.savefig("test.png", bbox_inches="tight")
    # plt.show()
