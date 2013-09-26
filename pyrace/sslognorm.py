import numpy as np
from .racemodel import *
import crace

class SSLogNormal(StopTaskRaceModel):
    """
    Stop-Signal LogNormal Race Model (only LN-accumulators).
    
    Need this layer for reimplementing fast likelihood (in C)
    """
    def get_cpars(self):
        """
        return many numpy arrays containing the parameters for convenient passing to C
        
        go_v : c1r1, c2r1, ..., cnr1, c1r2, c2r2, ..., cnr2
        """
        go_ter=np.zeros( self.design.nconditions()*self.design.nresponses(), dtype=np.float)
        go_mu=np.zeros_like(go_ter)
        go_sigma=np.zeros_like(go_ter)
        
        idx=0
        for resp in range(self.design.nresponses()):
            for cond in range(self.design.nconditions()):
                go_ter[idx]=self.go_accumulators[cond][resp].ter
                go_mu[idx]=self.go_accumulators[cond][resp].mu
                go_sigma[idx]=self.go_accumulators[cond][resp].sigma
                idx+=1
            
        stop_ter=np.zeros( self.design.nconditions(), dtype=np.float)
        stop_mu=np.zeros_like(stop_ter)
        stop_sigma=np.zeros_like(stop_ter)
        for cond in range(self.design.nconditions()):
            stop_ter[cond]  =self.stop_accumulators[cond].ter
            stop_mu[cond]=self.stop_accumulators[cond].mu
            stop_sigma[cond]  =self.stop_accumulators[cond].sigma
        
        pgf=np.array( self.prob_go_fail, dtype=np.float)
        ptf=np.array( self.prob_trigger_fail, dtype=np.float)
        self.cpars=(go_ter,go_mu,go_sigma,  stop_ter,stop_mu,stop_sigma, pgf,ptf)
        return self.cpars

    def likelihood_trials(self,dat):
        L=np.zeros(dat.ntrials, dtype=np.double)
        pars=self.get_cpars()+(L,)
        resp=np.array(dat.response, dtype=np.int32)
        crace.sslognorm_loglikelihood( self.design.nconditions(), self.design.nresponses(), dat.ntrials,
                                       dat.condition.astype(np.int32), resp, dat.RT, dat.SSD, *pars)
        return L

    def copy(self):
        m=self.__class__(self.design, self.params)
        return m

