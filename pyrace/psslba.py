import numpy as np
import pandas as pd
import scipy

from tools import *
from racemodel import *
from data import *
from lba import *


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
        
        
class pSSLBA_modelA(pSSLBA):
    def __init__(self, design, ster, ter, A, Bs, B, Vs, V, v, sv=1):
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
        self.params={'ster':ster, 'ter':ter, 'A':A, 'Bs':Bs, 'B':B, 'Vs':Vs, 'V':V, 'v':v}
        self.sv=sv
        self.set_params(**self.params)
        self.set_mixing_probabilities(1e-6,1e-6)
    
    def set_params(self, **kwargs):
        self.params.update(kwargs)
        go_acc=[]
        for cond in range(self.design.nconditions()):
            correct=self.design.correct_response(cond)
            print correct
            go_acc.append([ LBAAccumulator(self.params['ter'], self.params['A'], 
                                           self.params['V' if resp==correct else 'v'], 
                                           self.sv, self.params['A']+self.params['B'],
                                           name='go-'+":".join(self.design.condidx(cond))
                                             +'-%s'%('correct' if resp==correct else 'incorrect')) 
                           for resp in self.design.get_responses() ])
        
        stop_acc=[]
        for cond in range(design.nconditions()):
            stop_acc.append( LBAAccumulator( self.params['ster'], self.params['A'], self.params['Vs'], 
                                            self.sv, self.params['A']+self.params['Bs'] ) )
        
        self.set_accumulators(go_acc, stop_acc)


if __name__=="__main__":
    import pylab as pl


    factors=[{'sleepdep':['normal','deprived']},
             {'stimulus':['left', 'right']}]
    responses=['left','right']
    design=Design(factors,responses, 'stimulus')
    dat=pd.read_csv('../data/sleep_stop_onesubj_test.csv')
    print dat.shape[0]
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
    pl.plot(x, dens)

    cscore=0
    for i in range(design.nresponses()):
        score=scipy.integrate.quad(mod.dens_acc_go, 0.01, np.infty, args=(condition,i))[0]
        print "Chance of winning Acc %i (condition=%s): %f"%(i,str(design.condidx(condition)),score)
        cscore+=score
    print cscore

    L=mod.loglikelihood(ds)
    print L
    pl.show()
