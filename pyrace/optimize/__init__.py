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

    

