import os
import os.path
import shutil
from qpcr.aux.os.auxiliaries import *
from qpcr.aux.os.folders import *
from qpcr.aux.ops import *


# Main Functions ==============================================

def path_join( path, string):
    new_path = '{}/{}'.format(path, string)
    return new_path

def ext_extract( file_path):
    from os.path import splitext
    filename, file_ext = splitext(file_path)
    return file_ext
    
def itemcheck (item_path, filter = []):
    if filter == []:
        if os.path.isfile(item_path) == True:
            item = True
        else:
            item = False
    else:
        for filt in filter:
            if item_path.endswith(filt) == True:
                item = True
            else:
                item = False
    return item

def item_count(path, targets = None):
    count = 0
    if targets is not None:
        for i in os.listdir(path):
            if os.path.isfile(i):
                if os.path.basename(i) in targets:
                    count += 1
    else:
        for i in os.listdir(path):
            if os.path.isfile(i):
                count += 1
    return count 
