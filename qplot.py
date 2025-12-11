"""
This module defines the class TextPlot, which represents a simple ASCII scatter plot.

Exported classes:
    TextPlot: Class to create a simple ASCII scatter plot.

Exported functions:
    qplot: Function to print a simple ASCII plot to stdout.

Exported exceptions:
    None
"""

# stadard imports
from array import array


# TODO: (1) Add axis labels that print as a title line of the form: 'x (units) vs y (units)'
# TODO: (2) Add a plot legend, with one string per curve, showing the symbol and a description of the curve, to be
# printed as a title line.
class TextPlot():
    """
    Class to create a simple ASCII plot.
    Usage:
        (1) Create an instance of the TextPlot class with the parameters that define the plot.
        (2) Use the str() function to get the ASCII plot as a string.
        (3) Print the string or write it to a file as desired.
    """
    def __init__(self, ncur=1, npts=10, x=None, y=None, symbol=None, titl1='', titl2=''):
        """
        Provide all of the data needed to create the ASCII plot.
        :parameter ncur: Number of curves to plot, int
        :parmaeter npts: Number of points per curve, int
        :parameter x: X-coordinates of points, list of ncur arrays of npts floats, [array('f') ... array('f')]
        :parameter y: Y-coordinates of points, list of ncur arrays of npts floats, [array('f') ... array('f')]
        :parameter symbol: Plotting symbol for each curve, array of ncur characters, array('u')
        :parameter titl1: Title line 1, string
        :parameter titl2: Title line 2, string
        """
        assert(int(ncur)>0)
        self._ncur=int(ncur)
        assert(int(npts)>0)
        self._npts=int(npts)
        assert(type(x)==list)
        # TODO: Strictly speaking, if the code in _qplot() is changed, not all curves would have to have the same number of points.
        # Alhough the number of x-points and y-points for each curve would still have to match.
        for curve in x:
            assert(type(curve)==array)
            assert(len(curve)==int(self._npts))
        self._x=x
        assert(type(y)==list)
        for curve in y:
            assert(type(curve)==array)
            assert(len(curve)==int(self._npts))
        self._y=y
        assert(type(symbol)==array)
        assert(len(symbol)==int(self._ncur))
        self._symbol=symbol
        self._titl1=str(titl1)
        self._titl2=str(titl2)

    def _qplot(self):
        """
        Utility function called by __str__() to build the list of strings that make up the plot.
        :return: list of strings that make up the plot
        """
        # Characters for line of the plot
        cpl = 60
        # Lines (not curves, but ASCII lines) per plot
        lpp = 20

        # Set up the plot array
        line=[]
        # First, blank out the entire plot area with ' ' including the borders
        for i in range(lpp+1): # 0..lpp inclusive
            a_line = array('u')
            for j in range(cpl+1): # 0..cpl inclusive
                a_line.append(' ')
            line.append(a_line)
        # Now create top and bottome axes
        for j in range(cpl+1):
            line[0][j]='-'
            line[lpp][j]='-'
        # Now create x tick marks
        line[0][round(cpl/4)]='|'
        line[lpp][round(cpl/4)]='|'
        line[0][round(cpl/2)]='|'
        line[lpp][round(cpl/2)]='|'
        line[0][round(3*cpl/4)]='|'
        line[lpp][round(3*cpl/4)]='|'
        # Now create left and right axes
        for i in range(lpp+1):
            line[i][0]='|'
            line[i][cpl]='|'
        # Now create y tick marks
        line[round(lpp/4)][0]='+'
        line[round(lpp/4)][cpl]='+'
        line[round(lpp/2)][0]='+'
        line[round(lpp/2)][cpl]='+'
        line[round(3*lpp/4)][0]='+'
        line[round(3*lpp/4)][cpl]='+'

        # Find maximum and minimum values of x and y accross all curves, for scaling
        # TODO: May be a more Pythonic way to do this
        xmax=self._x[0][0]
        xmin=self._x[0][0]
        ymax=self._y[0][0]
        ymin=self._y[0][0]
        for j in range(self._ncur):
            for i in range(self._npts):
                if self._x[j][i]>xmax:
                    xmax=self._x[j][i]
                if self._x[j][i]<xmin:
                    xmin=self._x[j][i]
                if self._y[j][i]>ymax:
                    ymax=self._y[j][i]
                if self._y[j][i]<ymin:
                    ymin=self._y[j][i]

        # Fill in the points
        for j in range(self._ncur):
            for i in range(self._npts):
                yi=round((lpp)*(self._y[j][i]-ymin)/(ymax-ymin))
                xi=round((cpl)*(self._x[j][i]-xmin)/(xmax-xmin))
                line[yi][xi]=self._symbol[j]

        # Now build the list of strings that make up the plot

        plot_strings=[]
                
        # Add the titles
        plot_strings.append(self._titl1)
        plot_strings.append(self._titl2)
        plot_strings.append(' ') # Blank line

        # Add the top-axis tick labels
        xtick1=xmin
        xtick2=xmin+.25*(xmax-xmin)
        xtick3=xmin+.5*(xmax-xmin)
        xtick4=xmin+.75*(xmax-xmin)
        xtick5=xmax
        plot_strings.append(f"{xtick1:.2e}.......{xtick2:.2e}.......{xtick3:.2e}.......{xtick4:.2e}.......{xtick5:.2e}")

        # Add most of the graph
        for i in range(lpp,-1,-1):
            ytick=None
            if i==lpp or i==round(lpp*3/4) or i==round(lpp/2) or i==round(lpp/4) or i==0:
                ytick=ymin+(i*(ymax-ymin)/lpp)
            prnt_str=''
            for j in range(cpl+1):
                prnt_str += line[i][j]
            if ytick is not None:
                # Add the ytick label to the end of the line
                prnt_str += f"  {ytick:.2f}"
            plot_strings.append(prnt_str)

        # Add the bottom-axis tick labels
        plot_strings.append(f"{xtick1:.2e}.......{xtick2:.2e}.......{xtick3:.2e}.......{xtick4:.2e}.......{xtick5:.2e}")
        
        return plot_strings
    
    def __str__(self):
        """
        Return the ASCII plot as a string.
        :return: string representation of the ASCII plot
        """
        plot_lines = self._qplot()
        return '\n'.join(plot_lines)+'\n'


def qplot(ncur=1,npts=10,x=None,y=None,symbol=None,titl1='',titl2=''):
    """
    Print a simple ASCII plot to stdout.
    :parameter ncur: Number of curves to plot, int
    :parmaeter npts: Number of points per curve, int
    :parameter x: X-coordinates of points, list of ncur arrays of npts floats, [array('f') ... array('f')]
    :parameter y: Y-coordinates of points, list of ncur arrays of npts floats, [array('f') ... array('f')]
    :parameter symbol: Plotting symbol for each curve, array of ncur characters, array('u')
    :parameter titl1: Title line 1, string
    :parameter titl2: Title line 2, string
    :return: None
    """
    plot=TextPlot(ncur,npts,x,y,symbol,titl1,titl2)
    print(plot)
    return None
