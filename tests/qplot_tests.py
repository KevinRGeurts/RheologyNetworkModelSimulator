"""
This module contains unit tests for the TextPlot class and the qplot function.
"""


# standard imports
import unittest
from array import array
import sys
from io import StringIO

# local imports
from RheologyNetworkModelSimulator.qplot import TextPlot, qplot


class Test_qplot_tests(unittest.TestCase):
    def test_TextPlot_qplot(self):
        self.maxDiff = None
        ncur=2
        npts=10
        x=[array('f'),array('f')]
        y=[array('f'),array('f')]
        for i in range(npts):
            x[0].append(float(i))
            x[1].append(float(i))
            y[0].append(float(i*i))
            y[1].append(float(i*i*i))
        symbol=array('u')
        symbol.append('*')
        symbol.append('o')
        titl1='Test ASCII Plot'
        titl2='*: y = x^2; o: y = x^3'
        plot=TextPlot(ncur,npts,x,y,symbol,titl1,titl2)
        act_val = plot._qplot()
        exp_val=[]
        exp_val.append('Test ASCII Plot')
        exp_val.append('*: y = x^2; o: y = x^3')
        exp_val.append(' ')
        exp_val.append('0.00e+00.......2.25e+00.......4.50e+00.......6.75e+00.......9.00e+00')
        exp_val.append('|--------------|--------------|--------------|--------------o  729.00')
        exp_val.append('|                                                           |')
        exp_val.append('|                                                           |')
        exp_val.append('|                                                           |')
        exp_val.append('|                                                           |')
        exp_val.append('+                                                           +  546.75')
        exp_val.append('|                                                    o      |')
        exp_val.append('|                                                           |')
        exp_val.append('|                                                           |')
        exp_val.append('|                                                           |')
        exp_val.append('+                                                           +  364.50')
        exp_val.append('|                                              o            |')
        exp_val.append('|                                                           |')
        exp_val.append('|                                                           |')
        exp_val.append('|                                       o                   |')
        exp_val.append('+                                                           +  182.25')
        exp_val.append('|                                                           |')
        exp_val.append('|                                o                          |')
        exp_val.append('|                          o                         *      *')
        exp_val.append('|                   o            *      *      *            |')
        exp_val.append('o------o-----o-|----*------*--|--------------|--------------|  0.00')
        exp_val.append('0.00e+00.......2.25e+00.......4.50e+00.......6.75e+00.......9.00e+00')
        self.assertListEqual(exp_val, act_val)

    def test_TextPlot_dunder_str_dunder(self):
        self.maxDiff = None
        # Prepare test data
        ncur=2
        npts=10
        x=[array('f'),array('f')]
        y=[array('f'),array('f')]
        for i in range(npts):
            x[0].append(float(i))
            x[1].append(float(i))
            y[0].append(float(i*i))
            y[1].append(float(i*i*i))
        symbol=array('u')
        symbol.append('*')
        symbol.append('o')
        titl1='Test ASCII Plot'
        titl2='*: y = x^2; o: y = x^3'
        plot=TextPlot(ncur,npts,x,y,symbol,titl1,titl2)
        act_val = str(plot)
        print(f"Actual Value:")
        print(f"{act_val}")
        exp_val=''
        exp_val+='Test ASCII Plot\n'
        exp_val+='*: y = x^2; o: y = x^3\n'
        exp_val+=' \n'
        exp_val+='0.00e+00.......2.25e+00.......4.50e+00.......6.75e+00.......9.00e+00\n'
        exp_val+='|--------------|--------------|--------------|--------------o  729.00\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='+                                                           +  546.75\n'
        exp_val+='|                                                    o      |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='+                                                           +  364.50\n'
        exp_val+='|                                              o            |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                       o                   |\n'
        exp_val+='+                                                           +  182.25\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                o                          |\n'
        exp_val+='|                          o                         *      *\n'
        exp_val+='|                   o            *      *      *            |\n'
        exp_val+='o------o-----o-|----*------*--|--------------|--------------|  0.00\n'
        exp_val+='0.00e+00.......2.25e+00.......4.50e+00.......6.75e+00.......9.00e+00\n'
        print(f"Expected Value:")
        print(f"{exp_val}")
        self.assertEqual(exp_val, act_val)

    def test_qplot_function(self):
        self.maxDiff = None
        # Redirect standard output to a buffer
        captured_output = StringIO()
        sys.stdout = captured_output
        # Prepare test data
        ncur=2
        npts=10
        x=[array('f'),array('f')]
        y=[array('f'),array('f')]
        for i in range(npts):
            x[0].append(float(i))
            x[1].append(float(i))
            y[0].append(float(i*i))
            y[1].append(float(i*i*i))
        symbol=array('u')
        symbol.append('*')
        symbol.append('o')
        titl1='Test ASCII Plot'
        titl2='*: y = x^2; o: y = x^3'
        # Prepare expected output
        exp_val=''
        exp_val+='Test ASCII Plot\n'
        exp_val+='*: y = x^2; o: y = x^3\n'
        exp_val+=' \n'
        exp_val+='0.00e+00.......2.25e+00.......4.50e+00.......6.75e+00.......9.00e+00\n'
        exp_val+='|--------------|--------------|--------------|--------------o  729.00\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='+                                                           +  546.75\n'
        exp_val+='|                                                    o      |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='+                                                           +  364.50\n'
        exp_val+='|                                              o            |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                       o                   |\n'
        exp_val+='+                                                           +  182.25\n'
        exp_val+='|                                                           |\n'
        exp_val+='|                                o                          |\n'
        exp_val+='|                          o                         *      *\n'
        exp_val+='|                   o            *      *      *            |\n'
        exp_val+='o------o-----o-|----*------*--|--------------|--------------|  0.00\n'
        exp_val+='0.00e+00.......2.25e+00.......4.50e+00.......6.75e+00.......9.00e+00\n'
        # Call the qplot function, which prints to stdout
        qplot(ncur,npts,x,y,symbol,titl1,titl2)
        # Get the captured output
        act_val = captured_output.getvalue().strip()+'\n'
        # Reset the standard output
        sys.stdout = sys.__stdout__
        # Compare expected and actual output
        self.assertEqual(exp_val, act_val)


if __name__ == '__main__':
    unittest.main()
