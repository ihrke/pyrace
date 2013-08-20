import numpy as np
import scipy
from psslba import pSSLBA_modelA
from design import Design



def opt_func_deviance( x, mod, data ):
    xp=mod.untrans(x)
    mod.set_params(xp)
    score=mod.deviance(data)
#    print score, xp
    return score

class Optimizer:
    def __init__(self, model, data, opttype='simplex',  **kwargs):
        x0pars=model.params
        self.x0=model.trans(x0pars)
        self.model=model
        self.data=data
        self.opttype=opttype
        self.opts={}
        self.opts.update(kwargs)

    def optimize(self):
        if self.opttype=='simplex':
            method='Nelder-Mead'
        else:
            raise ValueError("don't know optimization method %s"%self.opttype)

        print "Optimize with %s and options %s"%(method, str(self.opts))
        if self.opttype=='simplex':
            r=scipy.optimize.fmin(opt_func_deviance, self.x0, (self.model, self.data),
                                  **(self.opts))

            self.xopt=r[0]
            self.fopt=r[1]
            self.iter=r[2]
            self.funcalls=r[3]
        return r

    def get_best(self):
        return self.model.untrans(self.xopt)

    def get_best_score(self):
        return self.fopt
    
#        return scipy.optimize.minimize( opt_func_deviance, self.x0, (self.model, self.data),
#                                        method=method, options=self.opts)
        

from collections import namedtuple


if __name__=="__main__":
    factors=[{'sleepdep':['normal','deprived']},
             {'stimulus':['left', 'right']}]
    responses=['left','right']
    design=Design(factors,responses, 'stimulus')

    
    #modelpars=namedtuple('modelpars_pSSLBA_mA', ['ster', 'ter', 'A', 'Bs', 'B', 'Vs', 'V', 'v'])
    
    truepar =pSSLBA_modelA.paramspec(.1, .2, .2,   .5,  .8, 2.0, 1.0, 0.0)
    print truepar
    startpar=pSSLBA_modelA.paramspec(.5, .1,  1., 2.0, 1.5, 1.0, 0.5, 0.9)
    print startpar
    mod=pSSLBA_modelA(design, truepar)
    dat=mod.simulate(1000, upper_limit=5)


    nopt=2
    for iopt in range(nopt):
        mod.set_params( startpar)
        opt=Optimizer(mod, dat, disp=1, xtol=1.0, ftol=1.0, full_output=1)
        res=opt.optimize()
        best=opt.get_best()
        print "Opt %i: Bestscore: "%iopt,opt.get_best_score()
        print "Opt %i: Best: "%iopt,best
        print "Opt %i: True: "%iopt,truepar
        print "Opt %i: Diff: "%iopt,np.abs(np.array(best)-np.array(truepar))
        startpar=best

    import pylab as pl
    pl.figure(figsize=(20,8))
    mod.set_params(best)    
    mod.plot_fit_go(dat)

    pl.figure(figsize=(20,8))
    mod.set_params(truepar)
    mod.plot_fit_go(dat)
    
    pl.show()

    

