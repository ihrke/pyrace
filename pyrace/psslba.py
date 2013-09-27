import numpy as np
from .racemodel import *
from .lba import *
import crace

class pSSLBA(StopTaskRaceModel):
    """Need this layer for reimplementing fast likelihood (in C)"""
    accumulator_type=LBAAccumulator

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
        self.cpars=(go_v,go_ter,go_A,go_b,go_sv, stop_v,stop_ter,stop_A,stop_b,stop_sv, pgf,ptf)
        return self.cpars

    def likelihood_trials(self,dat):
        L=np.zeros(dat.ntrials, dtype=np.double)
        pars=self.get_cpars()+(L,)
        resp=np.array(dat.response, dtype=np.int32)
        crace.sslba_loglikelihood( self.design.nconditions(), self.design.nresponses(), dat.ntrials,
                                     dat.condition.astype(np.int32), resp, dat.RT, dat.SSD, *pars)
        return L

    def deviance_precalc(self,dat):
        return -2*np.sum(np.log(np.maximum(self.likelihood_trials_precalc(dat),1e-10)))


    def copy(self):
        m=self.__class__(self.design, self.params)
        return m

    
    def likelihood_trials_precalc(self,dat):
        L=np.zeros(dat.ntrials, dtype=np.double)
        pars=self.cpars+(L,)
        resp=np.array(dat.response, dtype=np.int32)
        crace.sslba_loglikelihood( self.design.nconditions(), self.design.nresponses(), dat.ntrials,
                                     dat.condition.astype(np.int32), resp, dat.RT, dat.SSD, *pars)
        return L


"""
sv=1.0
V ~ deprivation * stop
v ~ deprivation
tf=0
gf=0

# hybrid II    
{'ter':None,
 'A':None,
 'B':['stop'], # -> B (all GO), Bs
 'V':['deprivation'], # -> Vc, Vd
 'v':['deprivation'],
 'sv':1.0,
 'tf':0,
 'gf':0}

# hybrid I
{'ter':None,
 'A':None,
 'B':['stop', 'deprivation'], # -> Bc, Bd, Bcs, Bds (all GO), Bs
 'V':['deprivation', 'stop'], # -> Vc, Vd
 'v':['deprivation'],
 'sv':1.0,
 'tf':0,
 'gf':0}



def create_model(**args):
    class paramspec(Parameters):
        parnames=[]
        lower=[]
        upper=[]
    class model(pSSLBA):
        def __init__(self, design, pars=None):
            self.design=design
            self.sv=1.0
            if pars!=None:
                self.set_params(pars)
            else:
                self.set_params(self.__class__.paramspec().random())
            self.set_mixing_probabilities(0,0)
            
    

model( B=['Bc', 'Bd', 'Bcs', 'Bds'],
       Bc={'deprivation':'control'},
       Bcs={'deprivation':'control', 'stop':True },
       Bd={'deprivation':'deprived'},
       Bds{'deprivation':'deprived'},
       bounds={'B':[]})
{'B': [ ('Bc', {'deprivation':'control'}
                
"""
