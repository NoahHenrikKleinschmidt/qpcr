
import pandas as pd
import numpy as np 
from scipy import interpolate

_upsample_smoothness = 0.0
_upsample_spline = 1

# interpolate more datapoints into the limited scale dataset
def upsample_coords(coords, n = 500):
    """
    Interpolates more datapoints into the limited-scale dataset

    Parameters
    ----------
    coords : iterable
        Measured absorption values as an iterable (tuple, list, np.ndarray)
    n : int
        The number of data points to interpolate (default is 500).
    Returns
    -------
    upsampled_coords : np.ndarray
        Upsampled absorption values
    """
    _xs = np.linspace(0, 1, n)

    tck, u = interpolate.splprep(
                                    coords, 
                                    k = _upsample_spline, 
                                    s = _upsample_smoothness
                                )
    upsampled_coords = interpolate.splev(_xs, tck)
    return upsampled_coords


def get_intersect(data, tline, use_nan = False):
    """
    Gets the index of the intersect between the absorption data and the threshold line.

    Parameters
    ----------
    data : np.ndarray
        Upsampled absorption values
    tline : np.ndarray
        A threshold line of the same length as `data`.
    use_nan : bool
        Sets the returned `idx` to `np.nan` if no intersect can be found if `use_nan = True`.
        Default is `use_nan = False` in which case the max possible index is 
        arbitrarily assigned.

    Returns
    -------
    idx : int
        The index of intersection between `data` and `tline`.
    """
    idx = np.sign(data - tline)
    idx = np.diff(idx)
    idx = np.argwhere(idx != 0).reshape(-1)

    # if no intersect can be found, just set default
    if idx.size == 0: 
        idx = np.nan if use_nan else len(data)-1
    else: 
        idx = idx[0]
    
    return idx
