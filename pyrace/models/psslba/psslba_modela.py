import numpy as np
__all__=['pSSLBA_modelA', 'pSSLBA_modelA_paramspec']

from ... import Parameters, pSSLBA, LBAAccumulator

class pSSLBA_modelA_paramspec(Parameters):
    """
    define the parameter-specification for model A
    """
    parnames=['ster', 'ter',  'A', 'Bs',  'B', 'Vs', 'V', 'v']
    lower   =[  1e-5,  1e-5, 1e-5, 1e-5, 1e-5,   -3,  -3,  -3]
    upper   =[     1,     1,    5,    5,    5,    3,   3,   3]
        
    
class pSSLBA_modelA(pSSLBA):
    """
    Example model. Different v for correct/incorrect GO and
    all other parameters shared among GO-accumulators.

    Separate B, ter and V for STOP.
    """
    
    paramspec=pSSLBA_modelA_paramspec;
        
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
        self.init_cmodule()
        self.design=design
        self.sv=1.0
        if pars!=None:
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

    def set_params_c(self,pars):
        """
        Optional for additional speedup by reducing
        the overhead of setting the accumulators every time.
        """
        self.params=pars
        nc=self.design.nconditions()
        nr=self.design.nresponses()
        
        go_v=np.zeros( nc*nr, dtype=np.float)
        go_ter=np.zeros_like(go_v)+(pars.ter)
        go_A=np.zeros_like(go_v)+(pars.A)
        go_b=np.zeros_like(go_v)+(pars.B+pars.A)
        go_sv=np.zeros_like(go_v)+(self.sv)
        
        idx=0

        # this is the index of the correct response per condition
        corridx=np.array([self.design.responses.index(self.design.correct_response(ci)) for ci in range(nc)])
        for resp in range(nr):
            go_v[(resp*nc)+np.where(corridx==resp)[0]]=pars.V
            go_v[(resp*nc)+np.where(corridx!=resp)[0]]=pars.v
                
        stop_v=np.zeros( nc, dtype=np.float)+(pars.Vs)
        stop_ter=np.zeros_like(stop_v)+(pars.ster)
        stop_A=np.zeros_like(stop_v)+(pars.A)
        stop_b=np.zeros_like(stop_v)+(pars.Bs+pars.A)
        stop_sv=np.zeros_like(stop_v)+(self.sv)
        
        pgf=np.zeros(nc, dtype=np.float)
        ptf=np.zeros(nc, dtype=np.float)        
        self.cpars=(go_v,go_ter,go_A,go_b,go_sv, stop_v,stop_ter,stop_A,stop_b,stop_sv, pgf,ptf)

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
