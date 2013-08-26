import numpy as np

class Parameters(object):
    parnames=[]
    lower   =[]
    upper   =[]

    def __init__(self, *args, **kwargs):
        self.pars={}        
        if len(args)>len(self):
            raise ValueError('too many input args')
        for i,arg in enumerate(args):
            self.pars[self.__class__.parnames[i]]=arg

        for k,v in kwargs.items():
            self.__setattr__(k,v)

    def __repr__(self):
        r="{cname}({parlist})".format(cname=self.__class__.__name__,
                                      parlist=",".join(["%s=%s"%(k,self.pars[k]) for k in self.__class__.parnames]))
        return r

    def __sub__(self, b):
        return self.__class__(*list(np.array([self.pars[k] for k in self.__class__.parnames]) \
                                       - np.array([b.pars[k] for k in b.__class__.parnames])))

    def __abs__(self):
        return self.__class__(*list(np.abs(np.array([self.pars[k] for k in self.__class__.parnames]))))
    
    def __getattr__(self, k):
        if k=='pars':
            return self.pars
        elif k in self.__class__.parnames:
            return self.pars[k]
        else:
            raise ValueError("don't have parameter '%s'"%k)

    def __setstate__(self,state):
        self.__dict__=state
        
    def __setattr__(self, k, v):
        if k=='__dict__':
            super(Parameters,self).__setattr__(k,v)
        elif k=='pars':
            self.__dict__[k]=v
        elif k in self.__class__.parnames:
            self.pars[k]=v
        else:
#            self.__dict__[k]=v
            raise ValueError("don't have parameter '%s'"%k)

    def __iter__(self):
        for k in self.__class__.parnames:
            yield self.pars[k]
    
    def __getitem__(self, i):
        if i<len(self.__class__.parnames) and i>=0:
            return self.pars[self.__class__.parnames[i]]
        else:
            raise ValueError("index out of bound: %i"%i)

    def __setitem__(self, i, v):
        if i<len(self.__class__.parnames) and i>=0:
            self.pars[self.__class__.parnames[i]]=v
        else:
            raise ValueError("index out of bound: %i"%i)

    def __len__(self):
        return len(self.__class__.parnames)


    def random(self):
        """
        Set the values to random elements within the borders
        set by lower and upper.
        """
        for i,(l,u) in enumerate(zip( self.__class__.lower, self.__class__.upper)):
            self.__setitem__(i, np.random.uniform( l, u ))
        return self
