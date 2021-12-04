"""
This module contains auxiliary functions to the qpcr module, 
that are not directly linked to qpcr Analysis per se. 
"""

def from_kwargs(key, default, kwargs, rm = False):
    """
    This function will try to extract key from the kwargs, 
    and if it fails it will return default
    """
    try: 
        if rm == False: 
            r = kwargs[key]
        else: 
            r = kwargs.pop(key)
    except: 
        r = default
    return r 
