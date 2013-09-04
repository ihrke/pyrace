import unittest
import numpy as np
import pandas as pd
import sys
from pyrace import Parameters
import pyrace

class testParameters(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        class mypar(Parameters):
            parnames=['ich', 'du', 'er']

        pars=mypar(1,2,3, er=5)
        print pars
        assert pars[0]==1
        assert pars[1]==2
        assert pars[2]==5
        assert len(pars)==3
        assert pars.ich==1
        assert pars.du==pars[1]
        assert pars.er==5

        assert tuple(list(pars))==(1,2,5), str(tuple(list(pars)))

    def test_random(self):
        class whatever(Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]

        pars=whatever()            
        for i in range(10):
            pars.random()
            assert pars.a>100 and pars.a<200
            assert pars.b>1 and pars.b<2
            assert pars.c>2 and pars.c<50

    def test_sub(self):
        class whatever(Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]

        pars=whatever()
        pars.random()
        assert np.all(np.abs( (pars-pars) ))<1e-10


    def test_in_range(self):
        class whatever(Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]
        pars=whatever().random()
        assert pars.in_range()
        out=whatever(a=50, b=1.5, c=5)
        assert not out.in_range()
        assert not pars.in_range(out)
        assert pars.in_range(pars)
    
    def test_abs(self):
        class whatever(Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]

        pars=whatever().random()
        pars2=whatever().random()
        pars3=abs(pars-pars2)
        assert np.all( pars3 )>=0
        
import pylab as pl
class testPlotting(unittest.TestCase):
    def test_plot_bar_parameters(self):
        class whatever(Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]
        pl.figure()
        a=whatever().random()
        b=whatever().random()
        pyrace.plot_bar_parameters(a,b,whatever().random())

    def test_plot_bar_parameters(self):
        class whatever2(Parameters):
            parnames=['a', 'b', 'c', 'd', 'ich', 'du']
            lower   =[100,   1,   2,   1,      1,   1]
            upper   =[200,   2,  50,   2,       2,   2]
        npar=10
        pars=[whatever2().random() for i in range(npar)]

        pl.figure()
        pyrace.plot_bar_parameters(*pars, title='%i tests'%npar)
        
        
    def tearDown(self):
        pl.show()
        
if __name__ == '__main__':
    unittest.main()
