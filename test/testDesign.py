import unittest
import numpy as np
import pandas as pd
import sys
from pyrace import Design, StopTaskDataSet, pnormP, dnormP, list_is_eq

class testDesign(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        factors=[{'sleepdep':['normal','deprived']},
                 {'stimulus':['left', 'right']}]
        responses=['left','right']
        design=Design(factors,responses, 'stimulus')
        assert design.nresponses()==len(responses)
        assert design.nconditions()==4
        assert len(design.factors_to_int)==len(design.factors_from_int)
        k=['deprived','left']
        assert design.condidx(  k )==2


class testDataSet(unittest.TestCase):
    def setUp(self):
        self.factors=[{'sleepdep':['normal','deprived']},
                      {'stimulus':['left', 'right']}]
        self.responses=['left','right']
        self.design=Design(self.factors,self.responses, 'stimulus')
        
    def test_init(self):
        factors=[{'sleepdep':['normal','deprived']},
                 {'stimulus':['left', 'right']}]
        responses=['left','right']
        design=Design(factors,responses, 'stimulus')
        dat=pd.read_csv('./data/sleep_stop_onesubj_test.csv')
        assert dat.shape[0]==800
        dat.columns=['sleepdep','stimulus','SSD','response','correct', 'RT']
        ds=StopTaskDataSet(design,dat)
        assert len(ds.RT)==len(ds.condition)
        assert len(ds.SSD)==len(ds.RT)
        assert len(ds.response)==len(ds.RT)

    def test_factorval(self):
        assert self.design.factorval(0,'sleepdep')=='normal'
        assert self.design.factorval(0,'stimulus')=='left'

        assert self.design.factorval(1,'sleepdep')=='normal'
        assert self.design.factorval(1,'stimulus')=='right'

        assert self.design.factorval(2,'sleepdep')=='deprived'
        assert self.design.factorval(2,'stimulus')=='left'

        assert self.design.factorval(3,'sleepdep')=='deprived'
        assert self.design.factorval(3,'stimulus')=='right'


class testTools(unittest.TestCase):
    def test_dnormP(self):
        x=np.linspace(-10,10,100)
        y=dnormP(x)
        assert y[0]==0
        assert y[-1]==0

        assert dnormP(-8.1, mean=-1)==0
        assert dnormP(100, mean=108)==0
        assert dnormP(116, mean=108)==0
        assert dnormP(110, mean=108)>0
        assert np.all(dnormP(x)>=0)
        assert dnormP(10, mean=10, sd=10)>0
        assert dnormP(0, mean=10, sd=10)>0
        assert dnormP(-100, mean=10, sd=10)==0
        

    def test_pnormP(self):
        x=np.linspace(-10,10,100)
        y=pnormP(x)
        assert y[0]==0
        assert y[-1]==1

        assert pnormP(-8.1, mean=-1)==0
        assert pnormP(6.1, mean=-1)==1   
        assert pnormP(100, mean=108)==0
        assert pnormP(116, mean=108)==1
        assert pnormP(110, mean=108)>0
        assert np.all(pnormP(x)>=0)
        assert np.all(pnormP(x)<=1)
        assert pnormP(10, mean=10, sd=10)>0
        assert pnormP(0, mean=10, sd=10)>0
        assert pnormP(-100, mean=10, sd=10)==0
        
    def test_list_is_eq(self):
        l1=[1,2,3]
        l2=[2,1,3]
        assert list_is_eq(l1,l2)
        assert list_is_eq(l1,l1)
        assert not list_is_eq( ['DU', 'hier'], [1,2])
        assert list_is_eq( ['du','hier'], ['hier','du'])
        assert not list_is_eq( ['du','hier'], ['hier'])
        
        
if __name__ == '__main__':
    unittest.main()
