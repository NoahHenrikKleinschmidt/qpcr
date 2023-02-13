"""
These are the stand-alone methods for interacting with the qpcr plotting environment.
"""

import qpcr.defaults as defaults


def plot(obj, mode: str = None, **kwargs):
    """
    A generic plotting shortcut to visualise the data from a `qpcr` class object, if it supports visualisation.

    Parameters
    ----------
    obj
        The object to visualise from.
    mode : str
        The plotting mode to use. This can be either "static" (matplotlib) or "interactive" (plotly)
    **kwargs
        Any additional keyword arguments.

    Returns
    -------
    fig
        The figure produced.
    """
    kwargs["mode"] = mode
    func = obj.__qplot__(**kwargs)
    fig = func(**kwargs)
    return fig


def interactive():
    """Set the default plotting mode to ``interactive``"""
    defaults.plotmode = "interactive"
    return defaults.plotmode


def static():
    """Set the default plotting mode to ``static``"""
    defaults.plotmode = "static"
    return defaults.plotmode
