import qpcr._auxiliary as aux

class EfficiencyCurve(aux._ID):
    """
    A helper class that will handle dilutions, ct values and the linreg model
    when newly computing efficiencies from assays.
    """
    __slots__ = ['_dilutions', '_ct_values', '_model', '_efficiency']
    
    def __init__(self, dilutions, ct_values, model, efficiency):
        super().__init__()
        self._dilutions = dilutions
        self._ct_values = ct_values
        self._model = model
        self._efficiency = efficiency
    
    def values(self):
        """
        Returns
        ------
        dilutions : np.ndarray
            The dilutions (x-values) used for efficiency calculation.
        ct_value : np.ndarray
            The underlying Ct values (y-values) used for efficiency calculation.
        """
        return self._dilutions, self._ct_values
    
    def model(self):
        """
        Returns
        -------
        model : stats.LinregressResult
            The linear regression model used for efficiency calculation.
        """
        return self._model
    
    def efficiency(self):
        """
        Returns
        -------
        eff : float
            The efficiency calculated from the stored data.
        """
        return self._efficiency