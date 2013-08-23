import numpy as np

from .. import Parameters, pSSLBA, LBAAccumulator
#import ..crace


class pSSLBA_modelB_paramspec(Parameters):
    """
    define the parameter-specification for model B
    """
    parnames=['ster', 'ter',  'A', 'Bs',  'B', 'Vs', 'Vtut', 'Vont', 'v']
    lower   =[  1e-5,  1e-5, 1e-5, 1e-5, 1e-5,   -3,     -3,     -3,  -3]
    upper   =[     1,     1,    5,    5,    5,    3,      3,      3,   3]
        
    
class pSSLBA_modelB(pSSLBA):
    """
    Different v for correct/incorrect GO and
    all other parameters shared among GO-accumulators.

    Let V differ between TUT and on-task.
    
    Separate B, ter and V for STOP.
    """
    
    paramspec=pSSLBA_modelB_paramspec;
        
    def __init__(self, design, pars=None):
        """
        ster - non-decistion time for all stop-accumulators
        ter - non-decision time for all go-accumulators
        A   - starting point var for all go+stop accumulators
        Bs  - distance b-A for all stop-accs
        B   - distance b-A for all go-accs
        Vs  - mean drift for all stop-accs
        Vtut- mean drift for all correct go-accs in TUT condition
        Vont- mean drift for all correct go-accs in on-task condition        
        v   - mean drift for all wrong go-accs
        """
        self.init_cmodule()
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
        x[6]=pars.Vtut
        x[7]=pars.Vont
        x[8]=pars.v
        return x
        
    def untrans(self, x):
        """reverse the process in trans, i.e., return a dictionary form vector"""
        pars=pSSLBA_modelB.paramspec(
            ster=np.exp(x[0]),
            ter=np.exp(x[1]),
            A=np.exp(x[2]),
            Bs=np.exp(x[3]),
            B=np.exp(x[4]),
            Vs=x[5],
            Vtut=x[6],
            Vont=x[7],
            v=x[8])
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
