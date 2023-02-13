"""
This is the superclass `StatsTest` that stores the basic attributes and methods that are common to the other Test classes.
"""

from qpcr import _auxiliary as aux
import qpcr.main as main


class StatsTest(aux._ID):
    """
    The superclass of other statistical Test classes.
    This class is able to link to a Results object (or Assay),
    and assign its results back, get the results to the user etc.
    """

    __slots__ = ['_obj', '_results', 'assaywise_results', 'groupwise_results']

    def __init__(self, id: str = None):
        super().__init__()
        self._id = id
        self._obj = None
        self.assaywise_results = None
        self.groupwise_results = None
        self._results = None

    def link(self, obj: main.Results):
        """
        Links a new object to evaluate.
        """
        self._obj = obj

    def get(self):
        """
        Returns
        ------
        results
            The results of the last performed test
        """
        return self._results

    def assign_comparisons(self):
        """
        Assigns the current results back to the source object.
        They can then be accessed via the ``comparisons`` attribute from there as well.
        """
        self._obj.add_comparisons(self._results)
