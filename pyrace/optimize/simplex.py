from .objective import *

__all__=['SimplexOptimizer']

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
        self.opts.update({'full_output':True})

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
