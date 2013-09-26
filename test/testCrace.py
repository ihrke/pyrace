import unittest
import numpy as np
import pandas as pd
import scipy
import sys
from pyrace import pnormP, dnormP, LBAAccumulator, Design, \
     StopTaskRaceModel, StopTaskDataSet
from pyrace.models.psslba import pSSLBA_modelA
import pyrace.crace
import pyrace as pr


class ToyWald_paramspec(pr.Parameters):
    parnames=['ter',  'b',  'V', 'Vs']
    lower   =[ 1e-5, 1e-5,  -5,   -5]
    upper   =[    1,    5,   5,    5]


class ModToyWald(pr.SSWald):
    paramspec=ToyWald_paramspec;
    
    def __init__(self, design, pars=None):
        self.design=design
        self.sv=1.0
        if pars!=None:
            self.set_params(pars)
        else:
            self.set_params(self.__class__.paramspec().random())
        self.set_mixing_probabilities(0,0)

    def set_params(self, pars):
        self.params=pars
        go_acc=[]
        for cond in range(self.design.nconditions()):
            go_acc.append([ pr.ShiftedWaldAccumulator( pars.b, pars.V, pars.ter, 
                                                      name='go-'+":".join(self.design.condidx(cond)) ) ])
        
        stop_acc=[]
        for cond in range(self.design.nconditions()):
            stop_acc.append( pr.ShiftedWaldAccumulator( pars.b, pars.Vs, pars.ter, 
                                                       name='stop-'+":".join(self.design.condidx(cond)) ) )
        
        self.set_accumulators(go_acc, stop_acc)

class ToyLN_paramspec(pr.Parameters):
    parnames=['ter',  'mu',  'sigma', 'sigmas']
    lower   =[ 1e-5,    -5,        0,        0]
    upper   =[    1,     5,       10,       10]


class ModToyLN(pr.SSLogNormal):
    paramspec=ToyLN_paramspec;

    def __init__(self, design, pars=None):
        self.design=design
        self.sv=1.0
        if pars!=None:
            self.set_params(pars)
        else:
            self.set_params(self.__class__.paramspec().random())
        self.set_mixing_probabilities(0,0)

    def set_params(self, pars):
        self.params=pars
        go_acc=[]
        for cond in range(self.design.nconditions()):
            go_acc.append([ pr.ShiftedLogNormalAccumulator( pars.ter, pars.mu, pars.sigma,
                                                      name='go-'+":".join(self.design.condidx(cond)) ) ])

        stop_acc=[]
        for cond in range(self.design.nconditions()):
            stop_acc.append( pr.ShiftedLogNormalAccumulator( pars.ter, pars.mu, pars.sigmas,
                                                       name='stop-'+":".join(self.design.condidx(cond)) ) )

        self.set_accumulators(go_acc, stop_acc)

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
    
    def test_wald_pdf(self):
        acc=pr.ShiftedWaldAccumulator(.2, .5, 1.0)
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        y=acc.pdf(x)
        y2=np.array([pyrace.crace.wald_pdf(xx, acc.alpha, acc.gamma, acc.theta) for xx in x], dtype=np.double)
        assert np.all( np.abs(y-y2)<1e-9)

    def test_wald_cdf(self):
        acc=pr.ShiftedWaldAccumulator(.2, .5, 1.0)        
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        y=acc.cdf(x)
        y2=np.array([pyrace.crace.wald_cdf(xx, acc.alpha, acc.gamma, acc.theta) for xx in x], dtype=np.double)        
        assert np.all( np.abs(y-y2)<1e-9)

    def test_lognorm_pdf(self):
        acc=pr.ShiftedLogNormalAccumulator(.2, .5, 1.0)
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        y=acc.pdf(x)
        y2=np.array([pyrace.crace.slognorm_pdf(xx, acc.ter, acc.mu, acc.sigma) for xx in x], dtype=np.double)
        assert np.all( np.abs(y-y2)<1e-9)

    def test_lognorm_cdf(self):
        acc=pr.ShiftedLogNormalAccumulator(.2, .5, 1.0)
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        y=acc.cdf(x)
        y2=np.array([pyrace.crace.slognorm_cdf(xx, acc.ter, acc.mu, acc.sigma) for xx in x], dtype=np.double)
        assert np.all( np.abs(y-y2)<1e-9)

    def test_sslba_likelihood_trials(self):
        factors=[{'sleepdep':['normal','deprived']},
                 {'stimulus':['left', 'right']}]
        responses=['left','right']
        design=Design(factors,responses, 'stimulus')
        dat=pd.read_csv('./data/sleep_stop_onesubj_test.csv')
        dat.columns=['sleepdep','stimulus','SSD','response','correct', 'RT']
        ds=StopTaskDataSet(design,dat)

        pars=pSSLBA_modelA.paramspec(.1, .2, .2, .5, .8, 2, 1, 0)
        mod=pSSLBA_modelA(design, pars)
        print mod

        Ls=StopTaskRaceModel.likelihood_trials(mod, ds)
        L2=mod.likelihood_trials(ds)

        goodix=(np.abs(Ls-L2)<1e-5)
        badix=np.logical_not(goodix)
        assert np.sum(badix)==0, "num bad idx=%i"%np.sum(badix)

    def test_sswald_likelihood_trials(self):
        factors=[{'deprivation':['normal', 'sleep']},
                 {'stimulus':['go']}]
        responses=['go']
        design=pr.Design(factors, responses, 'stimulus', name='singlego_wald')

        pars=ModToyWald.paramspec(ter=.3, b=1.0, V=1.0, Vs=2.0);
        mod=ModToyWald(design, pars)
        ds=mod.simulate(300)

        print mod

        Ls=StopTaskRaceModel.likelihood_trials(mod, ds)
        L2=mod.likelihood_trials(ds)

        goodix=(np.abs(Ls-L2)<1e-5)
        badix=np.logical_not(goodix)
        assert np.sum(badix)==0, "num bad idx=%i"%np.sum(badix)

    def test_sslognorm_likelihood_trials(self):
        factors=[{'deprivation':['normal', 'sleep']},
                 {'stimulus':['go']}]
        responses=['go']
        design=pr.Design(factors, responses, 'stimulus', name='singlego_ln')

        pars=ModToyLN.paramspec(ter=.3, mu=0, sigma=.1, sigmas=.2);
        mod=ModToyLN(design, pars)
        ds=mod.simulate(300)

        print mod

        Ls=StopTaskRaceModel.likelihood_trials(mod, ds)
        L2=mod.likelihood_trials(ds)

        goodix=(np.abs(Ls-L2)<1e-5)
        badix=np.logical_not(goodix)
        assert np.sum(badix)==0, "num bad idx=%i"%np.sum(badix)

    def test_loglikelihood(self):
        factors=[{'sleepdep':['normal','deprived']},
                 {'stimulus':['left', 'right']}]
        responses=['left','right']
        design=Design(factors,responses, 'stimulus')
        dat=pd.read_csv('./data/sleep_stop_onesubj_test.csv')
        dat.columns=['sleepdep','stimulus','SSD','response','correct', 'RT']
        ds=StopTaskDataSet(design,dat)

        pars=pSSLBA_modelA.paramspec(.1, .2, .2, .5, .8, 2, 1, 0)        
        mod=pSSLBA_modelA(design, pars)
        print mod

        Ls=StopTaskRaceModel.deviance(mod, ds)
        L2=mod.deviance(ds)
        assert np.abs(Ls-L2)<1e-8
        
if __name__ == '__main__':
    unittest.main()

