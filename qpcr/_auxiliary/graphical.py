import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns

def generate_palette(kwargs, return_check = False ):
    """
    Generates a color pallete for seaborn 
    plots to support a custom color argument.

    Optionally it can also return a boolean that checks if the entry was a color palette or a single color.
    """
    palette = kwargs.pop("palette", None )
    color = kwargs.pop("color", None )
    is_palette = True
    if color is not None: 
        try: 
            palette = sns.color_palette( color )
        except:
            palette = color
            is_palette = False

    if return_check:
        return palette, is_palette
    return palette

def make_layout(df, ref_column:str):
    """
    Generates a tuple for col, rows for subplots 
    based on the distinct sets of values within ref_column.
    """
    ref_col = df[ref_column]
    ref_col = list(set(ref_col))
    ref_length = len(ref_col)

    if ref_length == 1: 
        return (1, 1)

    if ref_length % 2 == 0:
        nrows = 2
        if ref_length % 4 == 0 and ref_length > 4:
            nrows = 4
    elif ref_length % 6 == 0:
            nrows = 6
    else: 
        nrows = 3

    ncols = int(np.ceil(ref_length / nrows))
    if ncols * nrows >= ref_length:
        return ncols, nrows
    else:
        ncols = ncols+1
        return ncols, nrows

def make_layout_from_list( ref_list ):
    """
    Generates a subplot layout based on a list instead of dataframe column
    """
    ref_col = list(set(ref_list))
    ref_length = len(ref_col)

    if ref_length == 1: 
        return (1, 1)

    if ref_length % 2 == 0:
        nrows = 2
        if ref_length % 4 == 0 and ref_length > 4:
            nrows = 4
    elif ref_length % 6 == 0:
            nrows = 6
    else: 
        nrows = 3

    ncols = int(np.ceil(ref_length / nrows))
    if ncols * nrows >= ref_length:
        return ncols, nrows
    else:
        ncols = ncols+1
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


class AxesCoords:
    """
    This class handles the coordinates for an Axes instance generated by subplots. 
    """
    def __init__(self, fig, axs, subplots:tuple):
        self._fig = fig
        self._axs = axs
        self._ncols, self._nrows = subplots
        self.col, self.row = None, None
        self._get_plot_type()
        self._init_row_and_col()
        self._coords = None

        # self._transpose = False
        
        self._first = self.col if self._type == "plt" else self.row
        self._second = self.row if self._type == "plt" else self.col
        self._limit_first = self._ncols if self._type == "plt" else self._nrows
        self._limit_second = self._nrows if self._type == "plt" else self._ncols
        
        self._check_2D()

        self._idx = 0

        self._2D_setup_ranges()

    def get(self):
        """
        Returns the current subplot coordinates
        """
        if self._is2D:
            coords = self._coords[self._idx]
        else:
            coords =  self._coords[self._idx][0]
        return coords

    def increment(self):
        """
        Switches to next subplot coordinates
        """
        self._idx += 1

    # def transpose(self, t:bool):
    #     """
    #     Switch col and row dimensions and reference attributes
    #     """    
    #     if t != self._transpose:
    #         self._transpose = t
    #         self._first = self.col if self._first == self.row else self.row
    #         self._second = self.col if self._second == self.row else self.row

    #         self._limit_first = self._ncols if self._limit_first == self._nrows else self._nrows
    #         self._limit_second = self._nrows if self._limit_second == self._ncols else self._ncols
    #         self._2D_setup_ranges()

    def subplot(self):
        """
        Returns the current subplot to be plotted (only for matplotlib figures)
        """
        coord = self.get()
        current_subplot = self._axs[coord]
        return current_subplot

    def _2D_setup_ranges(self):
        """
        Setup the possible coordinates for 2D grid
        """
        coords = []
        first_range = range(self._base, self._limit_first+self._base)
        second_range = range(self._base, self._limit_second+self._base)

        for first in first_range:
            for second in second_range:
                coords.append((first, second))
        self._coords = coords

    def _check_2D(self):
        """
        Checks if 2D coords will be required or not...
        """
        if self._type == "plotly": 
            self._is2D = True
        else:
            self._is2D = True if self._limit_second > 1 else False
        
    def _get_plot_type(self):
        """
        Sets if plotly or matplotlib figure are present, 
        and also checks if axs contain only one single subplot, 
        and makes it subscriptable if so...
        """
        fig, axs = plt.subplots()
        if type(self._fig) == type(fig):
            self._type = "plt"
        else:
            self._type = "plotly"
        
        if type(self._axs) == type(axs):
            self._axs = [self._axs]
        plt.close()

    def _setup_ranges(self):
        """
        Setup ranges of possible values
        """
        self._first_range = range(self._base, self._limit_first + self._base)
        self._second_range = range(self._base, self._limit_second + self._base)

    def _init_row_and_col(self):
        """
        Will set col and row to 1 if plotly figure is used, 
        otherwise 0 in case of matplotlib
        """
        if self._type == "plt":
            self._base = 0
        else: 
            self._base = 1
        self.col, self.row = (self._base, self._base)