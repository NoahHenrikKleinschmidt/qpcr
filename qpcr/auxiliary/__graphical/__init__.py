import numpy as np

def make_layout(df, ref_column:str):
    """
    Generates a tuple for col, rows for subplots 
    based on the distinct sets of values within ref_column.
    """
    ref_col = df[ref_column]
    ref_col = aux.sorted_set(ref_col)
    ref_length = ref_col.shape[0]
    if ref_length % 2 == 0:
        nrows = 2
        if ref_length % 4 == 0 and ref_length > 4:
            nrows = 4
    elif ref_length % 6 == 0:
            nrows = 6
    else: 
        nrows = 3
    ncols = int(np.ceil(ref_length / nrows))
    return ncols, nrows



def make_speclist(maxrows, maxcols, spectype):
    """
    This function sets up the plotly speclist for subplots 
    according to the number of plots and type of figure desired
    """
    spec_list = []
    for r in range(maxrows):
        tmp = []
        for c in range(maxcols):
            tmp.append({'type':spectype})
        spec_list.append(tmp)
    return spec_list