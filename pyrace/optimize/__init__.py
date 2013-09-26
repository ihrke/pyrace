import numpy as np
import scipy
import pylab as pl
import multiprocessing as mp
import time
import copy

#from .objective import *
from .simplex import *
from .diffev import *
from ..tools import ProgressBar
    
# --------------------------------------------------------------------------------
# MULTI-core (simultaneously for different datasets)
# 
#
#
# --------------------------------------------------------------------------------
    
def _run_optim( optimizer ):
    optimizer.optimize()
    return optimizer

def optimize_multi(model, data, pool=None, ncpu=2, start_points=None, progressbar=False, optimizer_pars={}):
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


    if not progressbar:
        opts=pool.map( _run_optim, optimizers )
    else:
        jobs=[pool.apply_async( _run_optim, (opt,) ) for opt in optimizers]
        prog=ProgressBar(len(jobs))
        waiting=True
        while waiting:
            time.sleep(.1)
            ndone=np.array([job.ready() for job in jobs]).sum()

            prog.animate(ndone)
            if ndone==len(jobs):
                waiting=False
        opts=[job.get() for job in jobs]

    return opts

    

