import numpy as np
import scipy
import pylab as pl
import multiprocessing as mp
import copy



def opt_func_deviance( x, mod, data, trace ):
    xp=mod.untrans(x)
    mod.set_params(xp)
    score=mod.deviance(data)
    if trace==100:
        print score, xp
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


def _cont_chunks(l, n):
    """
    Divide list l in n chunks (preserve order)
    """
    ll=[]
    c=int(np.ceil(len(l)/float(n)))
    for i in range(0, len(l), c):
        ll.append(l[i:i+c])
    return ll
    

def _flatten(l):
    return [item for sublist in l for item in sublist]


def de_rand_1_bin( pop_ini, objective, objargs=(), F=.5, CR=0, gen_max=1000, save_stats=None ):
    """
    Differential Evolution DE/rand/1/bin according to
    the Storn and Price paper (1997, Journal of Global Optimization).

    Direct implementation of algorithm in Figure 3.
    
    pop_ini : list of initial solutions (must have len and be in
              transformed space - np.array)

    save_stats : int or None; at which iterations, should stats be saved


    Returns: (final pop, final score, stats)
    """
    stopcrit=False
    npop=len(pop_ini)
    D=len(pop_ini[0])

    # initial evaluation
    score=np.zeros(npop)
    for i,cand in enumerate(pop_ini):
        score[i]=objective(cand, *objargs)

    cpop=copy.deepcopy(pop_ini)
    niter=0
    stats={'generation':[],
           'min':[],
           'max':[],
           'mean':[],
           'median':[]}
    while not stopcrit:
        # build candidate pop        
        pop_new=[]
        for i in range(npop): # population loop
            a,b,c=tuple(np.random.permutation(np.delete(np.arange(npop), i))[0:3])
            cand=np.zeros(D, dtype=np.float)
            j=np.random.randint(D)
            for k in range(D):
                if np.random.uniform()<CR or k==D-1:
                    cand[j]=cpop[c][j]+F*( cpop[a][j] - cpop[b][j])
                else:
                    cand[j] = cpop[i][j]
                j=(j+1) % D
            pop_new.append(cand)

        # update scores
        for i,cand in enumerate(pop_new):
            cscore=objective( cand, *objargs )
            if cscore<=score[i]:
                score[i]=cscore
                cpop[i]=cand

        # save stats
        if save_stats!=None and niter % save_stats == 0:
            stats['generation'].append(niter)
            stats['min'].append(np.min(score))
            stats['max'].append(np.max(score))
            stats['mean'].append(np.mean(score))
            stats['median'].append(np.median(score))

        if trace_stats!=None and niter % trace_stats == 0:
            print "DE(%i): mean=%.2f, median=%.2f, range=(%.2f, %.2f)"%(niter,stats['mean'][-1],
                                                                        stats['median'][-1],
                                                                        stats['min'][-1],stats['max'][-1])
            
        niter+=1
        if niter>=gen_max:
            stopcrit=True

    return cpop, score, stats if save_stats!=None else None


def _eval_pop( candidates, objective, objargs ):
    res=[objective(cand, *objargs) for cand in candidates]
    return res

def de_rand_1_bin_mp( pop_ini, objective, objargs=(), F=.5, CR=0, gen_max=1000, trace_stats=None, save_stats=None, pool=None, ncpu=4 ):
    """
    same as de_rand_1_bin except that each population is evaluated in parallel.

    pool : multiprocessing.Pool instance or None (own pool is created)
    ncpu : number of cores (not used if pool!=None)
    """
    if pool==None:
        pool=mp.Pool(ncpu)
    else:
        ncpu=pool._processes
    
    stopcrit=False
    npop=len(pop_ini)
    D=len(pop_ini[0])

    # initial evaluation
    candchunks=_cont_chunks(pop_ini, ncpu)
    jobs=[pool.apply_async( _eval_pop, (x, objective,objargs) ) for x in candchunks]
    score=np.array(_flatten([job.get() for job in jobs]))
    
    cpop=copy.deepcopy(pop_ini)
    niter=0
    stats={'generation':[],
           'min':[],
           'max':[],
           'mean':[],
           'median':[]}
    while not stopcrit:
        # build candidate pop        
        pop_new=[]
        for i in range(npop): # population loop
            a,b,c=tuple(np.random.permutation(np.delete(np.arange(npop), i))[0:3])
            cand=np.zeros(D, dtype=np.float)
            j=np.random.randint(D)
            for k in range(D):
                if np.random.uniform()<CR or k==D-1:
                    cand[j]=cpop[c][j]+F*( cpop[a][j] - cpop[b][j])
                else:
                    cand[j] = cpop[i][j]
                j=(j+1) % D
            pop_new.append(cand)

        # update score
        candchunks=_cont_chunks(pop_new, ncpu)
        jobs=[pool.apply_async( _eval_pop, (x, objective,objargs) ) for x in candchunks]
        new_score=np.array(_flatten([job.get() for job in jobs]), dtype=np.float)

        score=np.where( new_score<=score, new_score, score)
        cpop=[pop_new[i] if new_score[i]<=score[i] else cpop[i] for i in range(npop)]

        # save stats
        if save_stats!=None and niter % save_stats == 0:
            stats['generation'].append(niter)
            stats['min'].append(np.min(score))
            stats['max'].append(np.max(score))
            stats['mean'].append(np.mean(score))
            stats['median'].append(np.median(score))

        if trace_stats!=None and niter % trace_stats == 0:
            print "DE(%i): mean=%.2f, median=%.2f, range=(%.2f, %.2f)"%(niter,stats['mean'][-1],
                                                                        stats['median'][-1],
                                                                        stats['min'][-1],stats['max'][-1])
        niter+=1
        if niter>=gen_max:
            stopcrit=True

    return cpop, score, stats if save_stats!=None else None
    

        

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


class SimplexOptimizer(Optimizer):
    """Nelder-Mead Simplex"""
    def __init__(self, model, data, optfunc=opt_func_deviance, optfunc_pars=(), **kwargs):
        self.set_general_opts(**kwargs)
        self.opttype=self.__class__.__name__
        
        self.model=model
        self.data=data

        self.optfunc=optfunc
        self.optfunc_pars=optfunc_pars
        self.set_startpoint(model.params)

    def set_startpoint(self, pars):
        """pars is a model's paramspec"""
        self.x0pars=pars
        self.x0=self.model.trans(self.x0pars)
        self.initial=self.x0.copy()

    def optimize(self):
        if self.trace>=10:
            print "> Optimize with %s and options %s"%(self.opttype, str(self.opts))

        for i in range(self.noptimizations):
            r=scipy.optimize.fmin(self.optfunc, self.x0, (self.model, self.data, self.trace)+self.optfunc_pars,
                                  **(self.opts))
            rr={'opttype':self.opttype, 'nopt':(i+1), 'start':self.x0.copy(), 'xopt':r[0],'fopt':r[1],'iter':r[2],'funccalls':r[3]}
            self.results.append(rr)
            self.result=rr
            self.x0=rr['xopt'].copy()
            if self.trace>=10:
                print "> Optimization run %i: Bestscore: "%(i+1),self.get_best_score()
                print "> Optimization run %i: Best: "%(i+1),self.get_best()
        
        return self.get_best(), self.get_best_score()

class DEOptimizer(Optimizer):
    """Differential Evolution: DE/rand/1/bin"""
    def __init__(self, model, data, optfunc=opt_func_deviance, optfunc_pars=(), save_stats=None, pop_size=25, F=.5, CR=0, trace_stats=None, gen_max=1000, **kwargs):
        self.set_general_opts(**kwargs)
        self.opttype=self.__class__.__name__
        
        self.model=model
        self.data=data

        self.optfunc=optfunc
        self.optfunc_pars=optfunc_pars
        self.pop_size=pop_size
        self.opts.update({"F":F, 'CR':CR, "gen_max":gen_max, 'save_stats':save_stats, 'trace_stats':trace_stats})

    def optimize(self):
        if self.trace>=10:
            print "> Optimize with %s and options %s"%(self.opttype, str(self.opts))

        for i in range(self.noptimizations):
            pop_ini=[self.model.trans(self.model.paramspec().random()) for _ in range(self.pop_size)]
            if self.pool!=None:
                final_pop, score, stats=de_rand_1_bin_mp(pop_ini, self.optfunc, (self.model, self.data, self.trace)+self.optfunc_pars,
                                        pool=self.pool, **(self.opts))
            else:
                final_pop, score, stats=de_rand_1_bin(pop_ini, self.optfunc, (self.model, self.data, self.trace)+self.optfunc_pars,**(self.opts))
            idx=np.argsort(score)
            rr={'opttype':self.opttype, 'nopt':(i+1), 'xopt':final_pop[idx[0]],'fopt':score[idx[0]],
                'iter':self.opts['gen_max'],'funccalls':self.opts['gen_max']*self.pop_size, 'stats':stats}
            self.results.append(rr)
            self.result=rr
            self.x0=rr['xopt'].copy()
            if self.trace>=10:
                print "> Optimization run %i: Bestscore: "%(i+1),self.get_best_score()
                print "> Optimization run %i: Best: "%(i+1),self.get_best()
        
        return self.get_best(), self.get_best_score()

    def plot_stats(self):
        stats=self.result['stats']
        pl.plot( stats['generation'], stats['min'], label='min')
        pl.plot( stats['generation'], stats['max'], label='max')
        pl.plot( stats['generation'], stats['mean'], label='mean')
        pl.plot( stats['generation'], stats['median'], label='median')
        pl.legend()
        
    
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

    

