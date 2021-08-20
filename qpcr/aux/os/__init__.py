from os.path import dirname, abspath


def get_parent():
    """
    Will return the path of the current file's parent folder
    """
    return dirname(dirname(abspath(__file__)))
