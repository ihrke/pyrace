import numpy as np
from .racemodel import *
from .swald import *
import crace

class SSWald(StopTaskRaceModel):
    """
    Stop-Signal Wald Race Model (only Wald-accumulators).
    
    Need this layer for reimplementing fast likelihood (in C)
    """
    accumulator_type=ShiftedWaldAccumulator # define which accumulators are used by this model (exclusively)

    def get_cpars(self):
        """
        return many numpy arrays containing the parameters for convenient passing to C
        
        go_v : c1r1, c2r1, ..., cnr1, c1r2, c2r2, ..., cnr2
        """
        go_gamma=np.zeros( self.design.nconditions()*self.design.nresponses(), dtype=np.float)
        go_theta=np.zeros_like(go_gamma)
        go_alpha=np.zeros_like(go_gamma)
        
        idx=0
        for resp in range(self.design.nresponses()):
            for cond in range(self.design.nconditions()):
                go_gamma[idx]=self.go_accumulators[cond][resp].gamma
                go_theta[idx]=self.go_accumulators[cond][resp].theta
                go_alpha[idx]=self.go_accumulators[cond][resp].alpha
                idx+=1
                
        stop_gamma=np.zeros( self.design.nconditions(), dtype=np.float)
        stop_theta=np.zeros_like(stop_gamma)
        stop_alpha=np.zeros_like(stop_gamma)
        for cond in range(self.design.nconditions()):
            stop_gamma[cond]  =self.stop_accumulators[cond].gamma
            stop_theta[cond]=self.stop_accumulators[cond].theta
            stop_alpha[cond]  =self.stop_accumulators[cond].alpha
        
        pgf=np.array( self.prob_go_fail, dtype=np.float)
        ptf=np.array( self.prob_trigger_fail, dtype=np.float)
        self.cpars=(go_gamma,go_theta,go_alpha,  stop_gamma,stop_theta,stop_alpha, pgf,ptf)
        return self.cpars

    def likelihood_trials(self,dat):
        L=np.zeros(dat.ntrials, dtype=np.double)
        pars=self.get_cpars()+(L,)
        resp=np.array(dat.response, dtype=np.int32)
        crace.sswald_loglikelihood( self.design.nconditions(), self.design.nresponses(), dat.ntrials,
                                    dat.condition.astype(np.int32), resp, dat.RT, dat.SSD, *pars)
        return L

    def deviance_precalc(self,dat):
        return -2*np.sum(np.log(np.maximum(self.likelihood_trials_precalc(dat),1e-10)))


    def copy(self):
        m=self.__class__(self.design, self.params)
        return m

