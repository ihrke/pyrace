import unittest
import numpy as np
import pandas as pd
import scipy
import sys
from pyrace import Design, StopTaskDataSet
from pyrace import LBAAccumulator, pSSLBA, pSSLBA_modelA

class testLBA(unittest.TestCase):
    def setUp(self):
        pass

    def test_lba_acc(self):
        acc=LBAAccumulator(.2, .2, 1.0, 1.0, 1.0)
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        self.assertAlmostEqual(scipy.integrate.quad(acc.pdf, 0, np.infty)[0], 1)
        self.assertAlmostEqual(0, acc.cdf(0.01))
        assert abs( acc.cdf(10000)-1)<1e-4

    def test_sample(self):
        acc=LBAAccumulator(.2, .2, 2.0, 1.0, 1.0)
        nsamples=100000
        x=np.linspace(0,10, nsamples)
        
        import pylab as pl
        samp=acc.sample(nsamples)
        #dens=scipy.stats.gaussian_kde(samp[samp<10])
        #pl.hist(acc.sample(nsamples),200, normed=True)
        h,hx=np.histogram(samp, density=True, bins=1000)
        hx=hx[:-1]+(hx[1]-hx[0])/2.
        assert np.any(np.abs(h-acc.pdf(hx))<0.5)
        #pl.bar(hx, np.abs(h-acc.pdf(hx)))
        #pl.xlim(-1,1001)
        
        if False:
            #pl.subplot(2,1,1)
            pl.hist(samp[samp<10],300, normed=True, alpha=.3)
            pl.xlim(0,10)

            #pl.subplot(2,1,2)                
            pl.plot(x,acc.pdf(x), color='red', label='analytical')
            #pl.plot(x,dens(x),    color='green', label='kde')
            pl.legend()
            pl.show()
        


class testpSSLBA(unittest.TestCase):
    def setUp(self):
        self.factors=[{'sleepdep':['normal','deprived']},
                 {'stimulus':['left', 'right']}]
        self.responses=['left','right']
        self.design=Design(self.factors,self.responses, 'stimulus')
        self.dat=pd.read_csv('./data/sleep_stop_onesubj_test.csv')
        assert self.dat.shape[0]==800
        self.dat.columns=['sleepdep','stimulus','SSD','response','correct', 'RT']
        self.ds=StopTaskDataSet(self.design,self.dat)
        self.mod=pSSLBA_modelA(self.design, .2, .15, .2, 1.0, 1.0, 2, 1, 0.5)

    def test_simulate(self):
        print self.mod.simulate(100).as_dataframe(form='wide')
        
    def test_init(self):
        factors=[{'sleepdep':['normal','deprived']},
                 {'stimulus':['left', 'right']}]
        responses=['left','right']
        design=Design(factors,responses, 'stimulus')
        dat=pd.read_csv('./data/sleep_stop_onesubj_test.csv')
        assert dat.shape[0]==800
        dat.columns=['sleepdep','stimulus','SSD','response','correct', 'RT']
        ds=StopTaskDataSet(design,dat)
        mod=pSSLBA_modelA(design, .2, .15, .2, 1.0, 1.0, 2, 1, 0.5)

        print mod.parstring(full=True)
        print mod
        print mod.go_accumulators[0]
        nsamples=100
        x=np.linspace(0.01,10, nsamples)
        condition=0
        dens=mod.dens_acc_go(x,condition,1)

        cscore=0
        for i in range(design.nresponses()):
            score=scipy.integrate.quad(mod.dens_acc_go, 0.01, np.infty, args=(condition,i))[0]
            print "Chance of winning Acc %i (condition=%s): %f"%(i,str(design.condidx(condition)),score)
            cscore+=score
        assert abs(cscore-1)<1e-6
        print cscore

        L=mod.loglikelihood(ds)
        assert np.isfinite(L)
        print "L=",L


if __name__ == '__main__':
    unittest.main()
