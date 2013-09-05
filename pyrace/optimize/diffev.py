from .objective import *
import copy

__all__=['DEOptimizer']

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

def de_rand_1_bin( pop_ini, objective, objargs=(), F=.5, CR=0, gen_min=10, gen_max=1000,
                      nxtol=10, xtol=None, ftol=None, verbosity=0,
                      trace_stats=None, save_stats=None):
    """
    Differential Evolution DE/rand/1/bin according to
    the Storn and Price paper (1997, Journal of Global Optimization).

    Direct implementation of algorithm in Figure 3.
    
    pop_ini : list of initial solutions (must have len and be in
              transformed space - np.array)

    save_stats : int or None; at which iterations, should stats be saved

    ftol : computed from one population to the next with respect to the mean
           score of the population

    xtol : computed for best member which should change no more than xtol over nxtol generations
    nxtol : int

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

    best=[None for _ in range(nxtol)] # keep best individuals for nxtol generations
    prevmean=np.nan
    nxtoli=0
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
        for i,cand in enumerate(pop_new):
            cscore=objective( cand, *objargs )
            if cscore<=score[i]:
                score[i]=cscore
                cpop[i]=cand

        bestix=np.argmin(score)
        best[nxtoli]=cpop[bestix]
        nxtoli = (nxtoli+1) % nxtol
        
        # save stats
        if save_stats!=None and niter % save_stats == 0:
            stats['generation'].append(niter)
            stats['min'].append(np.min(score))
            stats['max'].append(np.max(score))
            stats['mean'].append(np.mean(score))
            stats['median'].append(np.median(score))

        if (trace_stats!=None and niter % trace_stats == 0) or (verbosity==100):
            print "DE(%i): mean=%.2f, median=%.2f, range=(%.2f, %.2f)"%(niter,stats['mean'][-1],
                                                                        stats['median'][-1],
                                                                        stats['min'][-1],stats['max'][-1])
        niter+=1

        ## convergence criteria
        if niter>gen_min:
            if niter>=gen_max:
                if verbosity>=10:
                    print "DE(%i) reached gen_max=%i"%(niter,gen_max)
                stopcrit=True
            if ftol!=None and np.abs(prevmean-np.mean(score))<ftol:
                if verbosity>=10:
                    print "DE(%i) reached ftol=%f"%(niter,ftol)
                stopcrit=True
            if xtol!=None and np.all(np.array([np.linalg.norm( best[0]-cur) for cur in best[1:] if cur!=None])<=xtol):
                if verbosity>=10:
                    print "DE(%i) reached xtol=%f for nxtol=%i"%(niter,xtol,nxtol)
                stopcrit=True

        prevmean=np.mean(score)
                
    return cpop, score, stats if save_stats!=None else None, niter


def _eval_pop( candidates, objective, objargs ):
    res=[objective(cand, *objargs) for cand in candidates]
    return res

def de_rand_1_bin_mp( pop_ini, objective, objargs=(), F=.5, CR=0, gen_min=10, gen_max=1000,
                      nxtol=10, xtol=None, ftol=None, verbosity=0,
                      trace_stats=None, save_stats=None, pool=None, ncpu=4 ):
    """
    same as de_rand_1_bin except that each population is evaluated in parallel.

    pool : multiprocessing.Pool instance or None (own pool is created)
    ncpu : number of cores (not used if pool!=None)

    ftol : computed from one population to the next with respect to the mean
           score of the population

    xtol : computed for best member which should change no more than xtol over nxtol generations
    nxtol : int    
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

    best=[None for _ in range(nxtol)] # keep best individuals for nxtol generations
    prevmean=np.nan
    nxtoli=0
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

        bestix=np.argmin(score)
        best[nxtoli]=cpop[bestix]
        nxtoli = (nxtoli+1) % nxtol
        
        # save stats
        if save_stats!=None and niter % save_stats == 0:
            stats['generation'].append(niter)
            stats['min'].append(np.min(score))
            stats['max'].append(np.max(score))
            stats['mean'].append(np.mean(score))
            stats['median'].append(np.median(score))

        if (trace_stats!=None and niter % trace_stats == 0) or (verbosity==100):
            print "DE(%i): mean=%.2f, median=%.2f, range=(%.2f, %.2f)"%(niter,stats['mean'][-1],
                                                                        stats['median'][-1],
                                                                        stats['min'][-1],stats['max'][-1])
        niter+=1

        ## convergence criteria
        if niter>gen_min:
            if niter>=gen_max:
                if verbosity>=10:
                    print "DE(%i) reached gen_max=%i"%(niter,gen_max)
                stopcrit=True
            if ftol!=None and np.abs(prevmean-np.mean(score))<ftol:
                if verbosity>=10:
                    print "DE(%i) reached ftol=%f"%(niter,ftol)
                stopcrit=True
            if xtol!=None and np.all(np.array([np.linalg.norm( best[0]-cur) for cur in best[1:] if cur!=None])<=xtol):
                if verbosity>=10:
                    print "DE(%i) reached xtol=%f for nxtol=%i"%(niter,xtol,nxtol)
                stopcrit=True

        prevmean=np.mean(score)
                
    return cpop, score, stats if save_stats!=None else None, niter
    

class DEOptimizer(Optimizer):
    """Differential Evolution: DE/rand/1/bin"""
    def __init__(self, model, data, optfunc=opt_func_deviance, optfunc_pars=(), save_stats=None,
                 xtol=None, ftol=None, nxtol=10, pop_ini=None,
                 pop_size=25, F=.5, CR=0, trace_stats=None, gen_max=1000, **kwargs):
        """
        xtol, ftol, nxtol : convergence parameters
        pop_size : int; size of population
        F, CR : DE parameters (see paper)
        gen_max : maximum number of generations (only when no convergence)
        """
        self.set_general_opts(**kwargs)
        self.opttype=self.__class__.__name__

        self.pop_ini=pop_ini
        self.model=model
        self.data=data

        self.optfunc=optfunc
        self.optfunc_pars=optfunc_pars
        self.pop_size=pop_size
        self.opts.update({"F":F, 'CR':CR, "gen_max":gen_max, 'save_stats':save_stats, 'trace_stats':trace_stats,
                          'xtol':xtol, 'ftol':ftol, 'nxtol':nxtol})

    def optimize(self, pop_ini=None):
        if pop_ini!=None:
            self.pop_ini=pop_ini
        if self.trace>=10:
            print "> Optimize with %s and options %s"%(self.opttype, str(self.opts))

        for i in range(self.noptimizations):
            if self.pop_ini==None:
                self.pop_ini=[self.model.trans(self.model.paramspec().random()) for _ in range(self.pop_size)]
                
            if self.pool!=None:
                final_pop, score, stats, ngen=de_rand_1_bin_mp(self.pop_ini, self.optfunc, (self.model, self.data, self.trace)+self.optfunc_pars,
                                        pool=self.pool, verbosity=self.trace, **(self.opts))
            else:
                final_pop, score, stats, ngen=de_rand_1_bin(self.pop_ini, self.optfunc, (self.model, self.data, self.trace)+self.optfunc_pars,
                                                      verbosity=self.trace, **(self.opts))
            idx=np.argsort(score)
            rr={'opttype':self.opttype, 'nopt':(i+1), 'xopt':final_pop[idx[0]],'fopt':score[idx[0]],
                'iter':ngen,'funccalls':ngen*self.pop_size, 'stats':stats}
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
