import unittest
import numpy as np
import pandas as pd
import scipy
import sys
from pyrace import pnormP, dnormP, LBAAccumulator
import pyrace.crace 

class testCRace(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_pnormP(self):
        x=np.linspace(-10,10,100)
        y=pnormP(x)
        y2=np.array([pyrace.crace.pnormP(xx,0,1) for xx in x], dtype=np.double)
        assert np.all( np.abs(y-y2)<1e-9)

        y=pnormP(x, mean=2, sd=3)
        y2=np.array([pyrace.crace.pnormP(xx,2,3) for xx in x], dtype=np.double)
        assert np.all( np.abs(y-y2)<1e-9)

    def test_dnormP(self):
        x=np.linspace(-10,10,100)
        y=dnormP(x)
        y2=np.array([pyrace.crace.dnormP(xx,0,1) for xx in x], dtype=np.double)
        assert np.all( np.abs(y-y2)<1e-9)
        
        y=dnormP(x, mean=2, sd=3)
        y2=np.array([pyrace.crace.dnormP(xx,2,3) for xx in x], dtype=np.double)
        assert np.all( np.abs(y-y2)<1e-9)

    def test_lba_pdf(self):
        acc=LBAAccumulator(.2, .5, 1.0, 1.0, 1.0)
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        y=acc.pdf(x)
        y2=np.array([pyrace.crace.lba_pdf(xx, acc.ter, acc.A, acc.v, acc.sv, acc.b) for xx in x], dtype=np.double)
        assert np.all( np.abs(y-y2)<1e-9)

    def test_lba_cdf(self):
        acc=LBAAccumulator(.2, .5, 1.0, 1.0, 1.0)
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        y=acc.cdf(x)
        y2=np.array([pyrace.crace.lba_cdf(xx, acc.ter, acc.A, acc.v, acc.sv, acc.b) for xx in x], dtype=np.double)
        assert np.all( np.abs(y-y2)<1e-9)
    

if __name__ == '__main__':
    unittest.main()

