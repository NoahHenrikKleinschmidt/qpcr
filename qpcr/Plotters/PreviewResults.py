"""
Defines the `PreviewResults` class, which is used to preview the results of a qpcr.Results object.

This is a wrapper for the Plotters:

- `AssayBars` (which was previously called `PreviewResults`, and is now the default setting)
- `GroupBars`
- `AssayDots`
- `GroupDots
"""

import qpcr.defaults as defaults
from qpcr import _auxiliary as aux
import qpcr.Plotters._base as base

import qpcr.Plotters.AssaySubplotResults as AssaySubplotResults
import qpcr.Plotters.GroupSubplotResults as GroupSubplotResults

# The plotters that are wrapped by PreviewResults
_plotters = {
    "AssayBars": AssaySubplotResults.AssayBars,
    "GroupBars": GroupSubplotResults.GroupBars,
    "AssayDots": AssaySubplotResults.AssayDots,
    "GroupDots": GroupSubplotResults.GroupDots,
}


class PreviewResults(base.Wrapper):
    """
    Generate a Preview of all results from all Assays in subplots.

    This is a wrapper for the Plotters:

    - `AssayBars` (which was previously called `PreviewResults`, and is now the default setting)
    - `GroupBars`
    - `AssayDots`
    - `GroupDots`

    Parameters
    ----------
    mode : str
        The plotting mode. May be either "static" (matplotlib) or "interactive" (plotly).

    kind : str
        The kind of Plotter to call. This can be any of the four wrapped
        Plotters, e.g. `kind = "GroupBars"`.
    """

    def __init__(self, mode: str = None, kind: str = None):

        super().__init__(kind=kind, mode=mode)

        if kind is None:
            kind = defaults.default_preview
        self._Plotter = _plotters[kind]
        self._Plotter = self._Plotter(mode=mode)
