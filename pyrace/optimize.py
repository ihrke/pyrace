import numpy as np
import scipy
from psslba import pSSLBA_modelA
from design import Design



def opt_func_deviance( x, mod, data, trace ):
    xp=mod.untrans(x)
    mod.set_params(xp)
    score=mod.deviance(data)
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
        else:
            raise ValueError("don't know optimization method %s"%self.opttype)

        if self.trace=='some':
            print "> Optimize with %s and options %s"%(method, str(self.opts))
        if self.opttype=='simplex':
            for i in range(self.noptimizations):
                r=scipy.optimize.fmin(self.optfunc, self.x0, (self.model, self.data, self.trace)+self.optfunc_pars,
                                      **(self.opts))
                rr={'opttype':self.opttype, 'nopt':i, 'start':self.x0.copy(), 'xopt':r[0],'fopt':r[1],'iter':r[2],'funccalls':r[3]}
                self.results.append(rr)
                self.result=rr
                self.x0=rr['xopt'].copy()
                if self.trace=='some':
                    print "> Optimization run %i: Bestscore: "%i,self.get_best_score()
                    print "> Optimization run %i: Best: "%i,self.get_best()
                    
        return self.get_best(), self.get_best_score()

    def get_best(self):
        return self.model.untrans(self.result['xopt'])

    def get_best_score(self):
        return self.result['fopt']
    

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

    

