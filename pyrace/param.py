import numpy as np
import collections


class Parameters(object):
    parnames=[]
    lower   =[]
    upper   =[]

    def __init__(self, *args, **kwargs):
        self.pars={k:np.nan for k in self.__class__.parnames}
        if len(args)>len(self):
            raise ValueError('too many input args')
        if len(args)==1 and isinstance(args[0], collections.Iterable):
            argl=np.array(args[0])
            args=argl

        for i,arg in enumerate(args):
            self.pars[self.__class__.parnames[i]]=arg

        for k,v in kwargs.items():
            self.__setattr__(k,v)

    def __repr__(self):
        r="{cname}({parlist})".format(cname=self.__class__.__name__,
                                      parlist=",".join(["%s=%.2f"%(k,self.pars[k]) for k in self.__class__.parnames]))
        return r

    def _repr_html_(self):
        """used by ipython notebook"""
        r="<table>\n"
        r+="<caption>%s</caption>\n"%self.__class__.__name__
        r+="<tr>"
        for nam in self.__class__.parnames:
            r+="<th>%s</th>"%nam
        r+="</tr>\n"
        r+="<tr>"
        for nam in self.__class__.parnames:
            r+="<td>%.2f</td>"%self.pars[nam]
        r+="</tr>\n</table>"
        return r
    
    def __sub__(self, b):
        return self.__class__( np.array(self) - np.array(b) )

    def __rsub__(self, b):
        return self.__class__( np.array(self) - np.array(b) )

    def __add__(self, b):
        return self.__class__( np.array(self) + np.array(b) )

    def __radd__(self, b):
        return self.__class__( np.array(self) + np.array(b) )

    def __mul__(self, b):
        return self.__class__( np.array(self) * np.array(b) )

    def __rmul__(self, b):
        return self.__class__( np.array(self) * np.array(b) )
    
    def __div__(self, b):
        return self.__class__( np.array(self) / np.array(b) )

    def __rdiv__(self, b):
        return self.__class__( np.array(self) / np.array(b) )

    def __abs__(self):
        return self.__class__( np.abs(np.array(self))) 
    
    def __getattr__(self, k):
        if k=='pars':
            return self.pars
        elif k in self.__class__.parnames:
            return self.pars[k]
        else:
            raise ValueError("don't have parameter '%s'"%k)

    def __setstate__(self,state):
        self.__dict__=state

    def __getstate__(self):
        return self.__dict__
        
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


    def bounds(self, par):
        """return bounds of parameter par (can be index or name)"""
        if isinstance(par,int): # index
            return (self.__class__.lower[par], self.__class__.upper[par])
        elif isinstance(par,str): # name
            ix=self.__class__.parnames.index(par)
            return (self.__class__.lower[ix], self.__class__.upper[ix])

    def bound_lower(self, par):
        return self.bounds(par)[0]
    def bound_upper(self,par):
        return self.bounds(par)[1]
    
    def in_range(self, pars=None):
        """
        return whether or not a specific setting of the parameters are
        in range or not. if pars==None, use current setting.
        """
        if pars==None:
            pars=self
        for i,(l,u) in enumerate(zip( self.__class__.lower, self.__class__.upper)):
           if pars[i]<l or pars[i]>u:
               return False
        return True 
        
    def random(self):
        """
        Set the values to random elements within the borders
        set by lower and upper.
        """
        for i,(l,u) in enumerate(zip( self.__class__.lower, self.__class__.upper)):
            self.__setitem__(i, np.random.uniform( np.maximum(l,-1e20), np.minimum(u,1e20) ))
        return self
