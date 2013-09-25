import unittest
import numpy as np
import scipy
import pandas as pd
import sys
import pyrace as pr
import pylab as pl
import os
import inspect

class PlottingEnabledTestCase(unittest.TestCase):
    def setUp(self):
        self.global_output_dir='test_output'
        cname=(self.__class__.__name__).split(".")[-1]        
        self.output_dir=os.path.join(self.global_output_dir, cname)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.saved=[]


    def numplots(self):
        return len(self.saved)
        
    def savefig(self, fname=None):
        if fname==None:
            curframe = inspect.currentframe()
            calname = inspect.getouterframes(curframe, 2)[1][3]
            fname="%03i_%s.png"%(self.numplots(),calname)
        fname=os.path.join(self.output_dir, fname)
        print "> saving ", fname
        self.saved.append(fname)
        pl.savefig(fname)
        pl.clf()
        pl.close()

    def tearDown(self):
        print "> saved ", self.saved
        
