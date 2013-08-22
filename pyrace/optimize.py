import numpy as np
import scipy
import pylab as pl
import multiprocessing as mp
import copy

from psslba import pSSLBA_modelA#, pSSLBA_modelA_paramspec
from design import Design



def opt_func_deviance( x, mod, data, trace ):
    xp=mod.untrans(x)
    mod.set_params(xp)
    score=mod.deviance(data)
    if trace=='full':
        print score, xp
    return score

def opt_func_deviance_precalc( x, mod, data, trace ):
    """
    works only for models which also implements set_params_c()
    """
    xp=mod.untrans(x)
    mod.set_params_c(xp)
    score=mod.deviance_precalc(data)
    if trace=='full':
        print score, xp
    return score

class Optimizer:
    def __init__(self, model, data, opttype='simplex', optfunc=opt_func_deviance, optfunc_pars=(), noptimizations=1, trace='some',  **kwargs):
        """
        trace : one of None, 'some', 'full'
        """
        x0pars=model.params
        self.x0=model.trans(x0pars)
        self.initial=self.x0.copy()
        self.model=model
        self.data=data
        self.opttype=opttype
        self.opts={}
        self.opts.update(kwargs)
        self.trace=trace
        self.noptimizations=noptimizations
        self.results=[]
        self.result=None
        self.optfunc=optfunc
        self.optfunc_pars=optfunc_pars

    def optimize(self):
        if self.opttype=='simplex':
            method='Nelder-Mead'
            self.opts['full_output']=1 ## setting this, because we rely on the results to be present
        else:
            raise ValueError("don't know optimization method %s"%self.opttype)

        if self.trace=='some':
            print "> Optimize with %s and options %s"%(method, str(self.opts))
        if self.opttype=='simplex':
            for i in range(self.noptimizations):
                r=scipy.optimize.fmin(self.optfunc, self.x0, (self.model, self.data, self.trace)+self.optfunc_pars,
                                      **(self.opts))
                rr={'opttype':self.opttype, 'nopt':(i+1), 'start':self.x0.copy(), 'xopt':r[0],'fopt':r[1],'iter':r[2],'funccalls':r[3]}
                self.results.append(rr)
                self.result=rr
                self.x0=rr['xopt'].copy()
                if self.trace=='some':
                    print "> Optimization run %i: Bestscore: "%(i+1),self.get_best_score()
                    print "> Optimization run %i: Best: "%(i+1),self.get_best()
                    
        return self.get_best(), self.get_best_score()



    def plot_convergence(self):
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


def _run_optim( optimizer ):
    optimizer.optimize()
    return optimizer

def optimize_multi(model, data, ncpu=2, start_points=None, optimizer_pars={}):
    """
    data : list of datasets which should be separately fit with the same model
    start_points : parameter set used as starting point for each dataset or None (than
              model.params is used for all datasets

    returns a list of Optimizer's
    """
    pool=mp.Pool(ncpu)

    optimizers=[Optimizer(model.copy(), dat, **optimizer_pars) for dat in data]
    if start_points!=None and len(start_points)!=len(optimizers):
       raise ValueError
    elif start_points!=None:
       for io,opt in enumerate(optimizers):
           opt.model.set_params(start_points[io])

    opts=pool.map( _run_optim, optimizers )

    return opts

    
if __name__=="__main__":
    factors=[{'sleepdep':['normal','deprived']},
             {'stimulus':['left', 'right']}]
    responses=['left','right']
    design=Design(factors,responses, 'stimulus')

    truepar =pSSLBA_modelA.paramspec(.1, .2, .2,   .5,  .8, 2.0, 1.0, 0.0)
    print truepar
    startpar=pSSLBA_modelA.paramspec(.5, .1,  1., 2.0, 1.5, 1.0, 0.5, 0.9)
    print startpar
    mod=pSSLBA_modelA(design, truepar)
    dat=mod.simulate(100, upper_limit=5)

    mod.set_params( startpar)
    opt=Optimizer(mod, dat, noptimizations=5, disp=1, xtol=1.0, ftol=1.0, full_output=1)
    best,bestscore=opt.optimize()

    if False:
        import pylab as pl
        pl.figure(figsize=(20,8))
        mod.set_params(best)    
        mod.plot_fit_go(dat)

        pl.figure(figsize=(20,8))
        mod.set_params(truepar)
        mod.plot_fit_go(dat)

        pl.show()

    

