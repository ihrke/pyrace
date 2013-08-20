import numpy as np
import pandas as pd
import scipy
from collections import namedtuple

from tools import *
from racemodel import *
from data import *
from lba import *
import crace

class pSSLBA(StopTaskRaceModel):
    """Need this layer for reimplementing fast likelihood"""
    def get_cpars(self):
        """
        return many numpy arrays containing the parameters for convenient passing to C
        
        go_v : c1r1, c2r1, ..., cnr1, c1r2, c2r2, ..., cnr2
        """
        go_v=np.zeros( self.design.nconditions()*self.design.nresponses(), dtype=np.float)
        go_ter=np.zeros_like(go_v)
        go_A=np.zeros_like(go_v)
        go_b=np.zeros_like(go_v)
        go_sv=np.zeros_like(go_v)
        idx=0
        for resp in range(self.design.nresponses()):
            for cond in range(self.design.nconditions()):
                go_v[idx]=self.go_accumulators[cond][resp].v
                go_ter[idx]=self.go_accumulators[cond][resp].ter
                go_A[idx]=self.go_accumulators[cond][resp].A
                go_b[idx]=self.go_accumulators[cond][resp].b
                go_sv[idx]=self.go_accumulators[cond][resp].sv
                idx+=1
                
        stop_v=np.zeros( self.design.nconditions(), dtype=np.float)
        stop_ter=np.zeros_like(stop_v)
        stop_A=np.zeros_like(stop_v)
        stop_b=np.zeros_like(stop_v)
        stop_sv=np.zeros_like(stop_v)
        for cond in range(self.design.nconditions()):
            stop_v[cond]  =self.stop_accumulators[cond].v
            stop_ter[cond]=self.stop_accumulators[cond].ter
            stop_A[cond]  =self.stop_accumulators[cond].A
            stop_b[cond]  =self.stop_accumulators[cond].b
            stop_sv[cond] =self.stop_accumulators[cond].sv
        
        pgf=np.array( self.prob_go_fail, dtype=np.float)
        ptf=np.array( self.prob_trigger_fail, dtype=np.float)
        return go_v,go_ter,go_A,go_b,go_sv, stop_v,stop_ter,stop_A,stop_b,stop_sv, pgf,ptf

    def likelihood_trials(self,dat):
        L=np.zeros(dat.ntrials, dtype=np.double)
        pars=self.get_cpars()+(L,)
        resp=np.array(dat.response, dtype=np.int32)
        crace.sslba_loglikelihood( self.design.nconditions(), self.design.nresponses(), dat.ntrials,
                                     dat.condition.astype(np.int32), resp, dat.RT, dat.SSD, *pars)
#        print "L(pSSLBA)=",L
#        LL=np.sum(np.log(np.maximum(L,1e-10)))
        return L
        
        
class pSSLBA_modelA(pSSLBA):
    paramspec=namedtuple('modelpars_pSSLBA_modelA', ['ster', 'ter', 'A', 'Bs', 'B', 'Vs', 'V', 'v'])
    
    def __init__(self, design, pars=None):
        """
        ster - non-decistion time for all stop-accumulators
        ter - non-decision time for all go-accumulators
        A   - starting point var for all go+stop accumulators
        Bs  - distance b-A for all stop-accs
        B   - distance b-A for all go-accs
        Vs  - mean drift for all stop-accs
        V   - mean drift for all correct go-accs
        v   - mean drift for all wrong go-accs
        """
        self.design=design
        self.sv=1.0
        self.set_params(pars)

        self.set_mixing_probabilities(0,0)

    def trans(self, pars):
        """pars is a paramspec
        
        Returns an array of ['ster', 'ter', 'A', 'Bs', 'B', 'Vs', 'V', 'v']
        in transformed space for the optimizer.
        """
        x=np.zeros(len(self.params), dtype=np.double)*np.nan
        x[0]=np.log(pars.ster)
        x[1]=np.log(pars.ter)
        x[2]=np.log(pars.A)
        x[3]=np.log(pars.Bs)
        x[4]=np.log(pars.B)
        x[5]=pars.Vs
        x[6]=pars.V
        x[7]=pars.v
        return x
        
    def untrans(self, x):
        """reverse the process in trans, i.e., return a dictionary form vector"""
        pars=pSSLBA_modelA.paramspec(
            ster=np.exp(x[0]),
            ter=np.exp(x[1]),
            A=np.exp(x[2]),
            Bs=np.exp(x[3]),
            B=np.exp(x[4]),
            Vs=x[5],
            V=x[6],
            v=x[7])
        return pars
        
    def set_params(self, pars):
        self.params=pars
        go_acc=[]
        for cond in range(self.design.nconditions()):
            correct=self.design.correct_response(cond)
            go_acc.append([ LBAAccumulator(self.params.ter, self.params.A, 
                                           self.params.V if resp==correct else self.params.v, 
                                           self.sv, self.params.A+self.params.B,
                                           name='go-'+":".join(self.design.condidx(cond))
                                             +'-%s'%('correct' if resp==correct else 'incorrect')) 
                           for resp in self.design.get_responses() ])
        
        stop_acc=[]
        for cond in range(self.design.nconditions()):
            stop_acc.append( LBAAccumulator( self.params.ster, self.params.A, self.params.Vs, 
                                             self.sv, self.params.A+self.params.Bs ) )
        
        self.set_accumulators(go_acc, stop_acc)

if __name__=="__main__":
    import pylab as pl


    factors=[{'sleepdep':['normal','deprived']},
             {'stimulus':['left', 'right']}]
    responses=['left','right']
    design=Design(factors,responses, 'stimulus')
    dat=pd.read_csv('./data/sleep_stop_onesubj_test.csv')
    print dat.shape[0]
    dat.columns=['sleepdep','stimulus','SSD','response','correct', 'RT']
    dat2=dat#.head(20)
    ds=StopTaskDataSet(design,dat2)
    
#    mod=pSSLBA_modelA(design, .2, .15, .2, 1.0, 1.0, 2, 1, 0.5)
    #start=c(ster=.1,ter=.2,A=.2,Bs=.5,B=.8,Vs=2,V=1,v=0)
    mod=pSSLBA_modelA(design, .1, .2, .2, .5, .8, 2, 1, 0)
    
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
    print cscore


    rts=np.arange(0.01,5,0.01)
    y=mod.dens_acc_go(rts,0,0)
    y2=mod.dens_acc_go(rts,0,1)
    pl.plot(rts,y)
    pl.plot(rts,y2)
    #pl.show()
    
    #L,Ls=mod.loglikelihood(ds)
    L2=mod.loglikelihood(ds)
    print L2
    #goodix=(np.abs(Ls-L2)<1e-5) & (L2>0)
    #badix=np.logical_not(goodix) & (L2>0)
    

    
#    print L2
#    print Ls
#    print L2
#    pl.show()

#    print -2*L
