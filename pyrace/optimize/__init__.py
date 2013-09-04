import numpy as np
import scipy
import pylab as pl
import multiprocessing as mp
import copy

#from .objective import *
from .simplex import *
from .diffev import *
    
# --------------------------------------------------------------------------------
# MULTI-core (simultaneously for different datasets)
# 
#
#
# --------------------------------------------------------------------------------
    
def _run_optim( optimizer ):
    optimizer.optimize()
    return optimizer

def optimize_multi(model, data, pool=None, ncpu=2, start_points=None, optimizer_pars={}):
    """
    data : list of datasets which should be separately fit with the same model
    start_points : parameter set used as starting point for each dataset or None (than
              model.params is used for all datasets

    pool : multiprocessing.Pool instance; if not given, construct own pool with ncpu cpus
    
    returns a list of Optimizer's
    """
    if pool==None:
        pool=mp.Pool(ncpu)

    optimizers=[SimplexOptimizer(model.copy(), dat, **optimizer_pars) for dat in data]
    if start_points!=None and len(start_points)!=len(optimizers):
       raise ValueError
    elif start_points!=None:
       for io,opt in enumerate(optimizers):
           optimizers[io].set_startpoint(start_points[io])

    opts=pool.map( _run_optim, optimizers )

    return opts

    
if __name__=="__main__":
    from .models.psslba_modela import pSSLBA_modelA
    from .design import Design

    
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
    opt=SimplexOptimizer(mod, dat, noptimizations=5, disp=1, xtol=1.0, ftol=1.0, full_output=1)
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

    

