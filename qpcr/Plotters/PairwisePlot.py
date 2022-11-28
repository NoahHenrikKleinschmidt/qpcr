"""
Defines the `PairwisePlot` class, which is used to preview the results of a pairwise-comparison.
"""

import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import qpcr._auxiliary.graphical as gx
import qpcr.Plotters._base as base
import seaborn as sns
import matplotlib.pyplot as plt

logger = aux.default_logger()

class PairwisePlot(base.Plotter):
    """
    Generate a summary figure for a ``PairwiseComparison`` object.
    
    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).

    """
    def __init__(self, mode : str = None):
        self._setup_default_params(
                                    static = defaults.static_Pairwise, 
                                    interactive = defaults.interactive_Pairwise
                                )
        super().__init__(mode)

    def link(self, obj : ("ComparisonsCollection" or "PairwiseComparison")):
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
        sns.set_style(kwargs.pop("style", None))
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
            fig, axes = plt.subplots( nrows = nrows, ncols = ncols, figsize = kwargs.pop("figsize", (ncols * 5, nrows * 5)) )
            for ax, data in zip(axes.flatten(), self.obj):
                self._static_single(data, ax, **kwargs)
        elif aux.pseudo_isinstance(self.obj, "PairwiseComparison"):
            ax = kwargs.pop("ax", None)
            if ax is None:
                fig, ax = plt.subplots( figsize=kwargs.pop("figsize", None) )
            self._static_single(self.obj, ax, **kwargs)
        
        else:
            raise TypeError("Cannot plot object of type '{}'".format(type(self.obj)))
        
    def _static_single(self, data, ax, **kwargs):
        """
        A single static Pairwise Plot
        """
        kwargs = self.update_params(kwargs)

        title, data = data.id(), data.stack()
        
        grid = kwargs.pop("grid", True)
        ax.grid(grid)
        if grid:
            ax.set_axisbelow(True)

        palette = gx.generate_palette(kwargs)

        xlabel, ylabel = kwargs.pop("xlabel", ""), kwargs.pop("ylabel", "")
        title = kwargs.pop("title", title)
        x, y, hue, size = kwargs.pop("x", "b"), kwargs.pop("y", "a"), kwargs.pop("hue", "pval_adj"), kwargs.pop("size", "effect_size")
        sizes = kwargs.pop("sizes", (20, 200))
        
        sns.scatterplot( data = data, x = x, y = y, hue = hue, size = size, 
                         sizes = sizes, palette = palette, legend = False, ax = ax, **kwargs )
        ax.set(xlabel = xlabel, ylabel = ylabel, title = title)
        ax.xaxis.set_tick_params( rotation = 90 )

if __name__ == "__main__":

    import qpcr

    assays = ['/Users/noahhk/GIT/qpcr/Examples/Example Data/HNRNPL_nmd.csv',
 '/Users/noahhk/GIT/qpcr/Examples/Example Data/HNRNPL_prot.csv',
 '/Users/noahhk/GIT/qpcr/Examples/Example Data/SRSF11_nmd.csv',
 '/Users/noahhk/GIT/qpcr/Examples/Example Data/SRSF11_prot.csv']
    normalisers = ['/Users/noahhk/GIT/qpcr/Examples/Example Data/28S.csv',
 '/Users/noahhk/GIT/qpcr/Examples/Example Data/actin.csv']
    
    assays = qpcr.read( assays, replicates = 3 )
    normalisers = qpcr.read( normalisers, replicates = 6 )


    assays = qpcr.delta_ct( assays )
    normalisers = qpcr.delta_ct( normalisers )

    results = qpcr.normalise( assays, normalisers )
    results.drop_rel()

    comparisons = qpcr.stats.assaywise_ttests( results )

    for i in comparisons:
        i.make_symmetric()

    plotter = PairwisePlot()
    plotter.link( comparisons[0] )
    plotter.plot( s = 500, style = "darkgrid", palette = "viridis_r" )
    plt.show()