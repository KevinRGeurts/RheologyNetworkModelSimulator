# FeneTrout Developer Documentation

## qplot Module

Working through the Fortran code to determine if I have the right indexing for the python implementation:

1. 1. Vertically speaking, or "y" speaking, in Fortran, the plot area is indexed from 0 to lpp, with the 0 and lpp lines being the borders.
So that the total number of lines of the plot area is lpp+1, and the total number of lines exclusive of the borders is lpp-1
	a. In Python, since lists index starting at 0, range(lpp+1) should give the correct number of lines including borders.
2. Horizontally speaking, or "x" speaking, in Fortran, the plot area is indexed from 0 to cpl, with the 0 and cpl columns being the borders.
So that the total number of columns of the plot area is cpl+1, and the total number of columns exclusive of the borders is cpl-1
	a. In Python, since lists index starting at 0, range(cpl+1) should give the correct number of columns including borders.

Example:
0 1 2 3 4 5 (cpl=5)
6 "cells" per line
5 increments from 0 to 5
so if the plot x is a real range fom 0 to 100, then the real increments are
0, 20, 40, 60, 80, 100
or 20 units per increment, or 100/cpl
