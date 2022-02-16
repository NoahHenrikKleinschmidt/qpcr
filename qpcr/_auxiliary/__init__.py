"""
This module contains auxiliary functions to the qpcr module, 
that are not directly linked to qpcr Analysis per se. 
"""
import uuid
import os 


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


def sorted_set(some_list):
    """
    Generates a sorted set of unique entries in a list.
    Importantly, sorted means it keeps the order of entries.
    """
    list_set = []
    for i in some_list:
        if i not in list_set: 
            list_set.append(i)
    return list_set

def same_type(obj1, obj2):
    """
    Compares the datatypes from two data objects and returns True if they are the same
    -> we use this to compare the types of the qpcr main classes since isinstance will not work for some reason...
    """
    type1 = type(obj1).__name__
    type2 = type(obj2).__name__
    same = type1 == type2
    return same

class _ID:
    """
    A meta_superclass that simply adds an ID getter-setter to itself
    """
    def __init__(self):
        self._id = uuid.uuid1() 
        self._id_was_set = False
    
    def id_was_set(self):
        return self._id_was_set

    def id(self, id:str = None):
        """
        Adds a string identifier or returns the given id
        """
        if id is not None and not self.id_was_set():
            self._id = id
            self._id_was_set = True
        else: 
            return self._id
    
    def adopt_id(self, obj):
        """
        Adopts the id of another objects that has an id() getter
        """
        self._id = obj.id()
        self._id_was_set = False
    
    def _id_reset(self):
        """
        Resets the memory if the id was already changed
        """
        self._id_was_set = False


def fileID(filename):
    """
    returns the basename of a filename for use as id
    """
    basename = os.path.basename(filename)
    basename = basename.split(".")[0]
    if not isinstance(basename, str):
        basename = ".".join(basename)
    return basename