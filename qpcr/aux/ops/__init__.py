import numpy as np 
import inspect 
import sys

def get_kwargs(function, check_kwargs = True, **kwargs):
    args = [k for k,v in inspect.signature(function).parameters.items()]
    if check_kwargs == True:
        temp = {k: kwargs.pop(k) for k in dict(kwargs) if k in args}
        args = temp
    return args

def in_kwargs(key, arg = None, **kwargs):
    if str(key) in kwargs:
        bool = True
    else:
        bool = False
    if arg is not None:
        if kwargs.get(str(key)) == arg:
            bool = True
        else: bool = False
    return bool

def ProgressBar(current, total, barLength = 20):
    percent = float(current) * 100 / total
    arrow   = '=' * int(percent/100 * barLength - 1) + '>'
    spaces  = ' ' * (barLength - len(arrow))
    print('Progress: [%s%s] %d %%' % (arrow, spaces, percent), end='\r')
    if percent == 100.0:
        print("")