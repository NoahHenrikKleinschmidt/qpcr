import os
import os.path
import shutil
from qpcr.aux.os.auxiliaries import *
from qpcr.aux.os.control import *



# Main Functions ==============================================

def subfolder_finder(path, mode = 'std', export=True, path_export = True, targets = None, excluded = None, sloppy = False,  msg = True):
    subfolders = []
    if mode == 'std':
        subfolders = subfolder_finder_core(path, export = export, path_export = path_export, msg = msg)
    elif mode == 'dyn':
        subfolders = search_subfolders_dynamic(path, export, export_paths = path_export)
    elif mode =='target':
        if sloppy == False:
            subfolders = search_subfolders_targeted(path = path, targets = targets, excluded = excluded)
        elif sloppy == True:
            subfolders = search_targeted(base_path = path, targeted=targets, excluded = excluded)
    if export == True:
        return subfolders

def search_subfolders_targeted(path, export = True, export_paths = True, targets=None, excluded=None):
    if export == True:
        message = False
    else:
        message = True
    subfolders = []
    subfolders_done = []
    new_path = path
    j = 0
    start = True
    while True:
        if new_path not in subfolders_done:
            subfolder_new = subfolder_finder(new_path, path_export = export_paths, msg = message)
            if subfolder_new != []:
                for i in subfolder_new:
                    if excluded is not None and targets is None:
                        for exc in excluded:
                            if exc not in i and i not in subfolders:
                                subfolders.append(i)
                    if excluded is None and targets is not None:
                        for tgt in targets:
                            if tgt in i and i not in subfolders:
                                subfolders.append(i)
                    elif excluded is not None and targets is not None:
                       for exc in excluded:
                           if exc not in i:
                               for tgt in targets:
                                   if tgt in i:
                                       subfolders.append(i)
                if start == True:
                    j = -1
                    start = False
            subfolders_done.append(new_path)
            j = j+1
            if j >= len(subfolders):
                break
            new_path = str(subfolders[j])
        else:
            j +=1
    if export == True:
        return subfolders

def search_subfolders_dynamic(path, export = True, export_paths = True):
    if export == True:
        message = False
    else:
        message = True
    subfolders = []
    subfolders_done = []
    new_path = path
    j = 0
    start = True
    while True:
        if new_path not in subfolders_done:
            subfolder_new = subfolder_finder(new_path, path_export = export_paths, msg = message)
            if subfolder_new != []:
                for i in subfolder_new:
                    subfolders.append(i)
                if start == True:
                    j = -1
                    start = False
            subfolders_done.append(new_path)
            j = j+1
            if j >= len(subfolders):
                break
            new_path = str(subfolders[j])
        else:
            j +=1
    if export == True:
        return subfolders



def search_targeted(base_path, mode = 'dyn', targeted = None, excluded = None):
    subfolders = []
    rawfolders = subfolder_finder(base_path, mode = mode)
    #print(base_path, '================')
    for folder in rawfolders:
        foldername = os.path.basename(folder)
        if targeted is not None and excluded is not None:
                if foldername in targeted:
                    exclude = False
                    for excl in excluded:
                        if excl in str(folder):
                            exclude = True
                        if exclude == False:
                            subfolders.append(folder)
        elif targeted is not None and excluded is None:
            if foldername in targeted:
                subfolders.append(folder)
        elif excluded is not None and targeted is None:
            for excl in excluded:
                if excl not in str(folder):
                    subfolders.append(folder)
                else:
                    continue
    #print('My final subfolders report: ', subfolders)
    return subfolders

def item_finder(path, filter = [], export = True, path_export = True):
    import os 
    import os.path as pth
    items = []
    for i in os.listdir(path):
        itempath = path_join(path, i)
        if itemcheck(itempath, filter) == True:
            if export == True:
                if path_export == True:
                    to_append = itempath
                else:
                    to_append = i
                items.append(to_append)
            else:
                print('I have found: {}'.format(i))   
    if export == True:
        return items

def filter_files(items_list, target_names):
    return_list = []
    for i in target_names:
        for j in items_list:
            if i in j: return_list.append(j)
    return return_list

def item_rename ( path, name, filter, numbering, i):
    for item in os.listdir(path):
        item_path= path_join(path, item)
        if itemcheck(item_path, filter) == True:
            ext = ext_extract(item_path)
            name_final = '{}{}{}'.format(name, i, ext)
            name_final = path_join(path, name_final)
            if numbering == True:
                i+=1
            os.rename(item_path, name_final)

#Housemaid functions ============================================

def clear_folder(path, filter = [], filtermode = 'exclude', mode = 'all'):
    if mode == 'files':
        items = item_finder(path, filter)
    elif mode == 'folders':
        items = subfolder_finder(path, msg = False)
    elif mode == 'all':
        items = os.listdir(path)
    for item in items:
        if filtermode == 'exclude':
            if os.path.basename(item) not in filter:
                path_del = path_join(path, item)
                clear_folder_core(path_del)
        elif filtermode == 'include':
            if os.path.basename(item) in filter:
                path_del = path_join(path, item)
                clear_folder_core(path_del)

def copy_folder(source, dest, mode = 'replace'):
    if mode == 'replace':
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.copytree(source, dest)

def copy_folder_select(source, dest, excluded = None, targets = None, mode = 'all'):
    if mode == 'all':
        for item in os.listdir(source):
            #print(os.listdir(source))
            if excluded is not None:
                if item not in excluded:
                    cfsc_core(source, dest, item)
                else: 
                    continue
            elif targets is not None:
                if item in targets:
                    cfsc_core(source, dest, item)
                else:                
                    continue
            else:
                cfsc_core(source, dest, item)
    elif mode == 'files':
        for item in os.listdir(source):
            itempath = path_join(source, item)
            destpath = path_join(dest, item)
            if itemcheck(itempath):
                if item not in excluded or item in targets:
                    shutil.copy(itempath, destpath)
    elif mode == 'folders':
        for item in os.listdir(source):
            itempath = path_join(source, item)
            destpath = path_join(dest, item)
            if os.path.isdir(itempath):
                if item not in excluded or item in targets:
                    if os.path.exists(destpath):
                        shutil.rmtree(destpath)
                    shutil.copytree(itempath, destpath)

