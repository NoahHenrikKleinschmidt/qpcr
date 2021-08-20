
import os
import os.path
import shutil
from  control import *
import sys

def cfsc_core(source, dest, item):
    try:
        copy_folder_select_core(source, dest, item)
    except PermissionError:
        destpath = path_join(dest, item)
        os.remove(destpath)
        copy_folder_select_core(source, dest, item)
    except FileNotFoundError:
        destpath = path_join(dest, item)
        shutil.rmtree(destpath)
        copy_folder_select_core(source, dest, item)
    except:
        sys.exit()

def copy_folder_select_core(source, dest, item):
    itempath = path_join(source, item)
    destpath = path_join(dest, item)
    #print('the itempath is: ', itempath)
    #print('the destiation is: ', destpath)
    if os.path.isfile(itempath):
        #print('copying: ', item)
        shutil.copy(itempath, destpath)
    elif os.path.isdir(itempath):
        shutil.copytree(itempath, destpath)

def subfolder_finder_core(path, export=True, path_export = True, msg = True):
    import os
    import os.path as pth
    subfolders = []
    for i in os.listdir(path):
        ipath = path_join(path, i)
        if pth.isdir(ipath) == True:
            if path_export == True:
                subfolders.append(ipath)
            else:
                subfolders.append(str(i))
    if msg == True:
        print('I have identified the following subfolders: {}'.format(subfolders))
    if export == True:
        return subfolders

def clear_folder_core(path):
    try:
        os.remove(path)
    except Exception as e:
        shutil.rmtree(path)
    except:
        print('I am unable to remove the file in question')
                
