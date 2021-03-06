import unittest
import numpy as np
import pandas as pd
import scipy
import sys
from pyrace import Design, StopTaskDataSet
from pyrace import LBAAccumulator, pSSLBA
from pyrace.models.psslba import pSSLBA_modelA

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
        pars=pSSLBA_modelA.paramspec(.2, .15, .2, 1.0, 1.0, 2, 1, 0.5)
        self.mod=pSSLBA_modelA(self.design, pars)

    def test_simulate(self):
        self.mod.simulate(100).as_dataframe(form='wide')
        
    def test_init(self):
        factors=[{'sleepdep':['normal','deprived']},
                 {'stimulus':['left', 'right']}]
        responses=['left','right']
        design=Design(factors,responses, 'stimulus')
        dat=pd.read_csv('./data/sleep_stop_onesubj_test.csv')
        assert dat.shape[0]==800
        dat.columns=['sleepdep','stimulus','SSD','response','correct', 'RT']
        ds=StopTaskDataSet(design,dat)
        pars=pSSLBA_modelA.paramspec(.2, .15, .2, 1.0, 1.0, 2, 1, 0.5)        
        mod=pSSLBA_modelA(design, pars)

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

    def test_set_params_c(self):
        #rpar=np.random.rand(8)
        rpar=pSSLBA_modelA.paramspec(*list(np.random.rand(8)))
        self.mod.set_params_c(rpar)
        cp1=self.mod.cpars
        self.mod.set_params(rpar)
        cp2=self.mod.get_cpars()
        parnames=['go_v','go_ter','go_A','go_b','go_sv', 'stop_v','stop_ter','stop_A','stop_b','stop_sv', 'pgf','ptf']
        ix=0
        for c1,c2 in zip(cp1,cp2):
            assert np.all(np.abs(c1-c2)<1e-10), "%i (=parameter %s): %s,%s"%(ix, parnames[ix], str(c1),str(c2))
            ix+=1

    def test_deviance_precalc(self):
        rpar=pSSLBA_modelA.paramspec(*list(np.random.rand(8)))
        self.mod.set_params_c(rpar)
        L1=self.mod.deviance_precalc(self.ds)

        self.mod.set_params(rpar)
        L2=self.mod.deviance(self.ds)

        assert abs(L1-L2)<1e-10, "L1,L2=%f, %f"%(L1,L2)
        
        
        
        
if __name__ == '__main__':
    unittest.main()
