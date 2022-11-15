"""
Defines the `GroupBars` and `GroupDots` class, which is used to preview the results of a qpcr.Results object, visualising the different groups in different subplots.
"""

import pandas as pd
import qpcr.defaults as defaults
import qpcr._auxiliary as aux
import qpcr._auxiliary.graphical as gx
import qpcr.Plotters._base as base

import matplotlib.pyplot as plt
import seaborn as sns 

import plotly.graph_objs as go

import numpy as np 

logger = aux.default_logger()

class GroupSubplotsResults(base.ResultsPlotter):
    """
    A superclass for FigureClasses that work with Results objects and operate group-wise
    (different groups in different subplots).
    """
    def __init__(self, mode : str ):
        super().__init__(mode)


    def _setup_groups_names_and_layout(self, data):
        """
        Generates unique groups and names from a dataframe, and sets up the layou (nrows,ncols)
        """
        groups = data["group"].unique()
        names = data["group_name"].unique()
        nrows, ncols = gx.make_layout_from_list(groups)
        return groups, names, nrows, ncols



class GroupBars(GroupSubplotsResults):
    """
    Generates a Bar plot figure with a separate subplot for each group. 

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).



    Plotting Kwargs
    ================
    
    `"static"` Kwargs

	    Static GroupBars figures accept the following kwargs:
    
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        |         Argument          |                                              Description                                               |                   Example                    |
        +===========================+========================================================================================================+==============================================+
        | show : `bool`             | Whether or not to show the figure                                                                      | `show = True` (default)                      |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | figsize : `tuple`         | The figure size                                                                                        | `figsize = (10, 4)`                          |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | title : `str`             | The overall figure title                                                                               | `title = "Today's results"                   |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | xlabel : `str`            | The x-axis label of each subplot                                                                       | `xlabel = "Conditions"`                      |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | ylabel : `str`            | The y-axis label of each subplot                                                                       | `ylabel = "Mean Mean Fold Change"`           |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | rot : `float`             | The rotation of x-axis labels                                                                          | `rot = 0.3`                                  |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | headers : `list`          | A list of titles for each subplot in the preview figure                                                | `headers = ["transcript A", "transcript B"]` |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | label_subplots  : `bool`  | Add each subplot with A, B, C ... (if True, default)                                                   | `label_subplots = True` (default)            |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | labeltype : `str`         | The starting character for subplot labelling. By default an `"A"`.                                     | `labeltype = "a"`                            |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | annotate_pvals: `bool`    | Add annotations of p-value significance to each subplot ( True by default if statistics are available) | `annotate_pvals = True` (default)            |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | pval_kws: `dict`          | Keywords for pvalue formatting using `encode_pvalues`                                                  | `pval_kws = dict(style = "*")`               |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | frame   : `bool`          | Show left and top spines of subplots (if True)                                                         | `frame = False` (default)                    |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | color : `str or list`     | The fillcolor for the individual bars                                                                  | `color = "yellow"`                           |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | style : `str`             | A `seaborn` style to set. Check out available styles [1].                                              | `style = "darkgrid"`                         |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | edgecolor : `str or list` | The edgecolor for the individual bars                                                                  | `edgecolor = "black"`                        |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | edgewidth : `float`       | The width of the edge of individual bars                                                               | `edgewidth = 0.5`                            |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | ecolor : `str`            | The color of errorbars                                                                                 | `ecolor = "orange"`                          |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+
        | **kwargs                  | Any additional kwargs that can be passed to the `matplotlib`-backend pandas `.plot.bar()` API.         |                                              |
        +---------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+


    
    `"interactive"` Kwargs

	    Interactive GroupBars figures accept the following kwargs:

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
    def __init__(self, mode : str = None ):
        self._setup_default_params(
                                    static = defaults.static_PreviewBars, 
                                    interactive = defaults.interactive_PreviewBars
                                )
        super().__init__( mode = mode )
    
    
    def _static_plot(self, **kwargs):
        """
        Generate a static figure
        """

        kwargs = self.update_params(kwargs)
        data = self.get()

        groups, names, nrows, ncols, xlabel, ylabel, x, y, sterr, headers, title, show = self._prep_shared_kwargs(kwargs, data)
        label_subplots, start_character, show_spines, rot = self._get_labels_and_spines_and_rot(kwargs)

        palette = gx.generate_palette(kwargs)

        # set a seaborn style
        style = kwargs.pop("style", "ticks")
        sns.set_style( style )

        edgecolor = kwargs.pop("edgecolor", "white")
        edgewidth = kwargs.pop("linewidth", 1)
        edgewidth = kwargs.pop("edgewidth", edgewidth)
        ecolor = kwargs.get("ecolor", "black")

        # get kwargs for pvalue annotations
        annotate_pvals = kwargs.pop("annotate_pvals", self._obj.comparisons is not None )
        pval_kws = kwargs.pop("pval_kws", {})

        if annotate_pvals:
            if self._obj.comparisons is None:
                raise AttributeError("Cannot annotate pvalues if no comparisons have been made!, please first perform a statistical groupwise comparison.")
     
        fig, coords = self._setup_static_figure(ncols, nrows, title, kwargs)

        idx = 0 
        for group, name in zip(groups, headers): 

            ax = coords.subplot()
            
            tmp_df = data.query( f"group == {group}" )
            tmp_df = tmp_df.sort_values( x )

            tmp_df.plot.bar(
                        x = x, 
                        y = y,
                        # yerr = "stdev",
                        color = palette,
                        edgecolor = edgecolor,
                        linewidth = edgewidth,
                        ax = ax,
                        **kwargs
                    )

            ax.errorbar(
                            x = tmp_df[x], y = tmp_df[y], 
                            yerr = tmp_df[sterr], 
                            fmt = ".", markersize = 0, capsize = 3, 
                            ecolor = ecolor,
                    )
            
            ax.set(
                            title = name,
                            xlabel = xlabel,
                            ylabel = ylabel,        
                )

            if rot is not None: 
                self._set_xtick_rotation( ax, rot )

            if not show_spines:
                sns.despine()

            # add ABCD... label to subplot
            if label_subplots:
                self._add_subplot_label(idx, ax, start_character)         

            # add stats annotations
            if annotate_pvals:
                logger.debug(f"{group=}")
                logger.debug(f"{names[group]=}")
                name = names[ group ]
                if name not in self._obj.comparisons:
                    logger.warning( f"Could not find a comparison for {name}. Perhaps you did not use a groupwise-comparison but an assaywise-comparison?" )
                else:
                    comparison = self._obj.comparisons[ name ]
                    self._annotate_pvalues(x, y, sterr, pval_kws, comparison, tmp_df, ax)


            idx += 1

        plt.tight_layout()

        if show:
            plt.show()
        
        return fig 


    def _interactive_plot(self, **kwargs):
        """
        The interactive Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self.get()

        # Get all the kwargs that are common to static and interactive
        groups, names, nrows, ncols, xlabel, ylabel, x, y, sterr, headers, title, show = self._prep_shared_kwargs(kwargs, data)
        hpad, vpad, height, width, template, hoverinfo = self._get_interactive_setup_kwargs(kwargs)

        fig, Coords = self._setup_interactive_figure(ncols, nrows, xlabel, ylabel, title, headers, hpad, vpad, height, width, template)


        fig.update_layout(
                            legend = {"title" : kwargs.pop("legend_title", "group name"),
                            },
                        )

        idx = 0
        for group, name in zip( groups, names ):
            row, col = Coords.get()
            tmp_df = data.query(f"group == {group}")
            tmp_df = tmp_df.sort_values( x )

            # now plot a new bar chart 
            fig.add_trace(

                go.Bar(
                    name = name,
                    y = tmp_df[y], x = tmp_df[x], 
                    error_y=dict(type='data', array = tmp_df[sterr]),
                    hoverinfo = kwargs.pop("hoverinfo", "y"), 
                    **kwargs
                ), 

                row, col
            )

            idx += 1

        if show:
            fig.show()

        # print("<<< Interactive Check >>>")
        return fig 

    def _prep_shared_kwargs(self, kwargs, data):
        """
        Pops all the shared kwargs between the static and interactive figure.
        """
        groups, names, nrows, ncols = self._setup_groups_names_and_layout(data)
        xlabel, ylabel = self._axeslabels(kwargs)
        x, y, sterr = self._get_xysterr(kwargs)
        x = "assay"
        headers = kwargs.pop("headers", names)
        title, show = self._get_title_and_show(kwargs)
        return groups,names,nrows,ncols,xlabel,ylabel,x,y,sterr,headers,title,show


class GroupDots(GroupSubplotsResults):
    """
    Generate a Preview of all results from all Assays in subplots, plotting individual values
    in a Dot Plot rather than Bar Plot.

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).


    Plotting Kwargs
    ================
    
    `"static"` Kwargs

	    Static GroupDots figures accept the following kwargs:

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
        | \*\*kwargs               | Any additional kwargs that can be passed to the `seaborn`'s `stripplot()`.                             |                                              |
        +--------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------+

    
    `"interactive"` Kwargs

	    Interactive GroupDots figures accept the following kwargs:

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
        | \*\kwargs                   | Any additional kwargs that can be passed to `plotly`'s`graphs_objs.Violin()`.                                                                                  |                                               |
        +-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+

    
    [1] `seaborn styles <https://www.python-graph-gallery.com/104-seaborn-themes>`_
    [2] `Plotly templates <https://plotly.com/python/templates/>`_
    [3] `Plotly hoverinfo <https://plotly.com/python/hover-text-and-formatting/>`_
    """
    def __init__(self, mode:str = None):
        self._setup_default_params(
                                    static = defaults.static_PreviewDots, 
                                    interactive = defaults.interactive_PreviewDots
                                )
        # __init__ is going to require default_params to be already set!
        super().__init__(mode)

    def _static_plot(self, **kwargs):
        """
        The static Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self._rep_data

        # get assays to plot
        assays = self._obj.data_cols #[ i for i in data.columns if i not in [ "group", "group_name", defaults.raw_col_names[0] ] ]
        assays = sorted( assays )

        # get kwargs incompatible with the main plotting method
        groups, names, nrows, ncols, xlabel, ylabel, headers, title, show, show_violins = self._prep_shared_kwargs(kwargs, data)       
        label_subplots, start_character, show_spines, rot = self._get_labels_and_spines_and_rot(kwargs)
        alpha = kwargs.pop("alpha", 1)

        # get kwargs for pvalue annotations
        annotate_pvals = kwargs.pop("annotate_pvals", self._obj.comparisons is not None )
        pval_kws = kwargs.pop("pval_kws", {})

        if annotate_pvals:
            if self._obj.comparisons is None:
                raise AttributeError("Cannot annotate pvalues if no comparisons have been made!, please first perform a statistical groupwise comparison.")
        
        # generate a custom color palette in case color kwarg is provided
        palette = gx.generate_palette(kwargs)

        # set a seaborn style
        style = kwargs.pop("style", "dark")
        sns.set_style( style )

        # make figure
        fig, Coords = self._setup_static_figure(nrows, ncols, title, kwargs)

        idx = 0
        for group, name in zip( groups, headers ): 
                
            # prepare a vertical bigtable format dataframe 
            tmp_df = self._prepare_df(data, assays, group)
            
            # now plot a new violin chart 
            subplot = Coords.subplot()

            if show_violins:
                sns.violinplot(
                                x = tmp_df[ "assay" ],
                                y = tmp_df[ "value" ],
                                color = None,
                                inner = None, 
                                palette = palette,
                                ax = subplot,

                            )
                for i in subplot.collections:
                    i.set_alpha( alpha * 0.3 )

            sns.stripplot(
                            x = tmp_df[ "assay" ],
                            y = tmp_df[ "value" ],
                            palette = palette,
                            alpha = alpha,
                            ax = subplot,
                            **kwargs
                        )
        
            subplot.set(
                        title = name,
                        xlabel = xlabel,
                        ylabel = ylabel,        
                    )

            # adjust xtick rotation   
            if rot is not None: 
                self._set_xtick_rotation( subplot, rot )

            if not show_spines:
                sns.despine()

            # add ABCD... label to subplot
            if label_subplots:
                self._add_subplot_label(idx, subplot, start_character)         

            # add stats annotations
            if annotate_pvals:
                name = names[ group ]
                if name not in self._obj.comparisons:
                    logger.warning( f"Could not find a comparison for {name}. Perhaps you did not use a groupwise-comparison but an assaywise-comparison?" )
                else:
                    comparison = self._obj.comparisons[ name ]
                    self._annotate_pvalues("assay", "value", None, pval_kws, comparison, tmp_df, subplot )

            idx += 1
            
        plt.tight_layout()

        if show:
            plt.show()

        # print("<<< Static Check >>>")

        return fig 

    def _prepare_df(self, data, assays, group):
        """
        Concatenates the different assay colums into a single 
        set of two columns, one for all "value"s and one for 
        the "assay" identifiers. 
        """
        # get a group subset 
        group_df = data.query( f"group == {group}" )

        # concatenate all separate assay columns into a set of 
        # 'value' and 'assay' columns.

        # setup the _assays concatenated dataframe   
        _assays = group_df[ [ "group", assays[0] ] ]
        _assays = _assays.rename( columns = { assays[0] : "value" } )
        _assays["assay"] = assays[0] # [ assays[0] for i in _assays["group"] ]

        # now iteratively add all remaining assays
        for assay in assays[1:]:
            tmp_df = group_df[ [ "group", assay ] ]
            tmp_df = tmp_df.rename( columns = { assay : "value" } )
            tmp_df["assay"] = assay # [ assay for i in tmp_df["group"] ]
            _assays = pd.concat( [_assays, tmp_df], ignore_index = True)
                
        # remove the group col, and sort
        tmp_df = _assays.drop( columns = ["group"] )
        tmp_df = tmp_df.sort_values( "assay" )

        return tmp_df

    def _interactive_plot(self, **kwargs):
        """
        The interactive Preview Results Figure
        """

        kwargs = self.update_params(kwargs)
        data = self._rep_data

        # get assays to plot
        assays = self._obj.data_cols # [ i for i in data.columns if i not in [ "group", "group_name", defaults.raw_col_names[0] ] ]
        assays = sorted( assays )
        ticks = np.arange( len(assays) )

        # get the groups
        # and incompatible kwargs
        groups, names, nrows, ncols, xlabel, ylabel, headers, title, show, show_violins = self._prep_shared_kwargs(kwargs, data)
        hpad, vpad, height, width, template, hoverinfo = self._get_interactive_setup_kwargs(kwargs)
        
        fig, Coords = self._setup_interactive_figure(ncols, nrows, xlabel, ylabel, title, headers, hpad, vpad, height, width, template)
        # setup figure
        fig.update_layout( showlegend = False )

        idx = 0
        for group, name in zip( groups, headers ):
            
            # prepare a vertical bigtable format dataframe 
            tmp_df = self._prepare_df(data, assays, group)

            row, col = Coords.get()

            # now plot a new violin chart 
            fig.add_trace(

                go.Violin(
                                name = name,
                                x = tmp_df[ "assay" ], 
                                y = tmp_df[ "value" ], 
                                points = "all",
                                pointpos = 0,
                                hoverinfo = hoverinfo, 
                                **kwargs
                        ), 
                                row, col
                    )

            # remove violins if not desired (leaving only dot plot)
            if not show_violins: 
                fig.update_traces(
                                    fillcolor="rgba(0,0,0,0)", 
                                    line_width = 0, 
                                    selector=dict(type='violin'),
                                )


            # update x axis to categorical group names
            fig.update_layout(
                                { 
                                    f"xaxis{idx+1}" : dict(
                                                            tickmode = 'array',
                                                            tickvals = ticks,
                                                            ticktext = assays,
                                                    )         
                                }
                            )

            idx += 1

        if show:
            fig.show()

        # print("<<< Interactive Check >>>")
        return fig 


    def _prep_shared_kwargs(self, kwargs, data):
        """
        Pops all the shared kwargs between the static and interactive figure.
        """
        groups, names, nrows, ncols = self._setup_groups_names_and_layout(data)
        xlabel, ylabel = self._axeslabels(kwargs)
        headers = kwargs.pop("headers", names)
        title, show = self._get_title_and_show(kwargs)
        show_violins = kwargs.pop("violin", True)
        return groups, names, nrows, ncols, xlabel, ylabel, headers, title, show, show_violins


if __name__ == "__main__":

    import qpcr
    # we set up the paths to 28S+actin as our normalisers
    normaliser_files = [
                            "../Examples/Example Data/28S.csv",
                            "../Examples/Example Data/actin.csv"
                    ]

    # we also set up the paths to the HNRNPL and SRSF11 transcripts
    assay_files = [
                        "../Examples/Example Data/HNRNPL_nmd.csv",
                        "../Examples/Example Data/HNRNPL_prot.csv",
                        "../Examples/Example Data/SRSF11_nmd.csv",
                        "../Examples/Example Data/SRSF11_prot.csv",
                ]

    a = qpcr.read( assay_files, replicates = (6,5,7,6) )
    n = qpcr.read( normaliser_files, replicates = (6,5,7,6) )

    a,n = qpcr.delta_ct( [a,n] )
    r = qpcr.normalise(a,n)
    r.drop_rel()

    # works :-)
    p = GroupBars( mode = "interactive" )
    p.link(r)
    fig = p.plot()
    fig.show()

    # works :-)
    p = GroupDots( mode = "static" )
    p.link(r)
    fig = p.plot()
    plt.show()

    
