import numpy as np
import pylab as pl
import scipy

def opt_func_deviance( x, mod, data, trace ):
    xp=mod.untrans(x)
    mod.set_params(xp)
    score=mod.deviance(data)
    if trace==100:
        print score, xp
    if np.isnan(score):
        print "opt_func_deviance: found nan at ",x, "-->", xp
        score=np.infty
    return score

def opt_func_deviance_precalc( x, mod, data, trace ):
    """
    works only for models which also implements set_params_c()
    """
    xp=mod.untrans(x)
    mod.set_params_c(xp)
    score=mod.deviance_precalc(data)
    if trace==100:
        print score, xp
    return score


class Optimizer:
    def set_general_opts( self, noptimizations=1, trace=10, pool=None, ncpu=None, **kwargs ):
        """
        trace : int; on a scale from 0-100 verbosity (no output to max output)
        
        For optimizers that support multicore within dataset
        (for cross-dataset parallelization, call optimize_multi):
        
        pool : multiprocessing.Pool instance or None
        ncpu : int or None
        """    
        self.trace=trace
        self.noptimizations=noptimizations
        self.opts={}
        self.opts.update(kwargs)
        self.results=[]
        self.result=None
        if pool!=None:
            self.pool=pool
            self.ncpu=pool._processes
        else:
            if ncpu!=None and ncpu>1:
                self.pool=mp.Pool(ncpu)
                self.ncpu=ncpu
            else:
                self.pool=None
                self.ncpu=None
        
    def optimize(self):
        raise NotImplementedError


    def plot_convergence_nopt(self):
        # objective function over noptimizations in repeated optimization
        pl.subplot(2,1,1)
        pl.plot( 0, self.optfunc( self.x0, self.model, self.data, self.trace, *(self.optfunc_pars)), 'x',
                 markersize=10, markeredgewidth=2, color='red')
        pl.plot( [v['nopt'] for v in self.results], [v['fopt'] for v in self.results], '-o')
        pl.xlabel('number of optimization run')
        pl.ylabel('score')
        pl.xlim(-.1, len(self.results))

        # num iterations
        pl.subplot(2,1,2)
        pl.plot( [v['nopt'] for v in self.results], [v['iter'] for v in self.results], '-.', label='niter')
        pl.plot( [v['nopt'] for v in self.results], [v['funccalls'] for v in self.results], '-o', label='funccalls')
        pl.legend()
        pl.xlabel('number of optimization run')
        pl.ylabel('num')
        pl.xlim(-.1, len(self.results))
                
    
    def get_best(self):
        return self.model.untrans(self.result['xopt'])

    def get_best_score(self):
        return self.result['fopt']
