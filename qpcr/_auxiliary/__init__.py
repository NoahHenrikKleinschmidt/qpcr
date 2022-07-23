"""
This module contains auxiliary functions to the qpcr module, 
that are not directly linked to qpcr Analysis per se. 
"""
import uuid
import os 
import re 
import qpcr.defaults as defaults
import logging


def log( filename = None, level = None, format = None, name = "qpcr" ):
    """
    Store a log file

    Parameters
    -------
    filename : str
        The file in which to save the log. By default this will be a file `qpcr.log`
        (logger name + .log) in the same directory as the current main script.
        The log can be directed to stdout using `filename = "stdout"`.
    level 
        The logging level. By default WARNING is used.
    format  
        The logging format.
    name
        The logger name to use. By default `qpcr` is used.
    
    """

    logger = logging.getLogger( name = name )

    if logger.hasHandlers():
        if any( [ i.level == level for i in logger.handlers ] ):
            return logger

    # setup the dedicated qpcr logger
    if level is None:
        level = defaults.log_level
    if not logger.hasHandlers():
        logger.setLevel( level )

    if filename == "stdout":
        handler = logging.StreamHandler()
    else:
        if filename is None:

            # try to use a fileHandler but if we have a jupyter notebook or terminal
            # just use StreamHandler instead...
            try:
                import __main__
                filename = f"{os.path.dirname( __main__.__file__ )}/{name}.log"
            except:
                return log( filename = "stdout", level = level, format = format, name = name )
        
        elif not filename.endswith(".log"): 
            filename += ".log"
        handler = logging.FileHandler( filename = filename )
   
    if format is None: 
        format = defaults.log_format
    
    f = logging.Formatter( fmt = format )
    handler.setFormatter( f )
    handler.setLevel( level )
    logger.addHandler( handler )
    return logger

def default_logger():
    """The default logger of qpcr"""
    logger = logging.getLogger( name = "qpcr" )
    if not logger.hasHandlers():
        return log( filename = defaults.init_log_loc )
    return logger

def extensive_logger():
    """
    The extensive logger of qpcr.
    """
    logger = logging.getLogger( name = "qpcr" )
    if not logger.hasHandlers():
        return log( level = logging.DEBUG )
    logger.setLevel( logging.DEBUG )
    return logger

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
    __slots__ = ["_id", "_id_was_set", "_id_label", "_id_label_was_set", "_id_func",  "__dict__"]

    def __init__(self):
        self._id = uuid.uuid1() 
        self._id_was_set = False
        self._id_label_was_set = False
        self._id_label = None

        self._id_func = self._strict_id if defaults.strict_id else self._relaxed_id
    
    def id( self, id : str = None ):
        """
        Adds a string identifier or returns the given id
        """
        return self._id_func( id )

    def id_was_set(self):
        """
        Returns True if an Id was set
        """
        return self._id_was_set
    
    def _strict_id(self, id:str = None):
        """
        Adds a string identifier or returns the given id
        """
        if id is not None and not self.id_was_set():
            self._id = id
            self._id_was_set = True
        return self._id 

    def _relaxed_id( self, id:str = None ):
        """
        Adds a string identifier or returns the given id
        """
        if id is not None:
            self._id = id
            self._id_was_set = True
        return self._id

    def id_label(self, label:str = None):
        """
        Adds a label and gets the current label.
        The label serves as a secondary id that helps distinguish
        separate objects from each other that still belong together, 
        like like qpcr.Assays that store data from transcript isoforms for instance.
        """
        if label is not None and not self.id_label_was_set():
            self._id_label = label
            self._id_label_was_set = True
        return self._id_label

    def split_id(self, by : str, regex = False):
        """
        Splits an id into id + label based on `by`.
        `by` may either be a string by which to apply `split` to the id
        or a `regex pattern` which includes exactly _two_ capturing groups (one for id, one for label).
        Note, this requires `regex = True` and may render the original id irrestorable through merge_id !
        """
        if regex:
            by = re.compile(by)
            split = by.search(self._id).groups()
        else: 
            split = self._id.split(by)
        
        if len(split) > 1:
            id, label = split
        else: 
            id, label = split[0], None
        
        # set new id and label
        self._id = id
        self._id_label = label

    def merge_id(self, by : str = "_"):
        """
        Merges the id and id_label together into just id, and places `by` 
        between them (default just an underscore). It resets the id_label.
        If no label is present, nothing will happen.
        """
        new_id = self.get_merged_id(by = by)
        self._id = new_id

        # reset the label settings
        self._id_label = None
        self._id_label_was_set = False

    def get_merged_id(self, by : str = "_"):
        """
        Returns a merged id but does not adopt it!
        """
        label = by + self._id_label if self._id_label is not None else ""
        new_id = self._id + label
        return new_id

    def id_to_label(self):
        """
        Switches id and id_label
        """
        self._id, self._id_label = self._id_label, self._id

    def id_label_was_set(self):
        """
        Returns True if an Id label was set
        """
        return self._id_label_was_set

    def adopt_id(self, obj):
        """
        Adopts the id of another objects that has an id() getter
        """
        self._id = obj.id()
        self._id_was_set = False
    
    def id_reset(self):
        """
        Resets the memory if the id and id_label were already changed
        """
        self._id_label_was_set = False
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



if __name__ == "__main__":

    obj1 = _ID()
    obj1.id("Hnrnp l_nmd")
    obj1.split_id("_")
    # obj1.id_to_label()
    # obj1.id_to_label()
    # obj1.merge_id()
    r = obj1.id()
    print(r)

    print("......")
    obj2 = _ID()
    obj2.id("Hnrnp l_prot")
    obj2.split_id("_")
    # obj2.id_to_label()
    # obj2.id_to_label()
    # obj2.merge_id()
    r = obj2.get_merged_id()
    print(r)
    r = obj2.id()
    print(r)