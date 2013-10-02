import numpy as np
import scipy
import scipy.stats as stats

from .tools import *
from .racemodel import *

class LBAAccumulator(Accumulator):
    parnames=[   'ter',        'A',       'v',     'sv',      'b']
    lower   =[       0,          0, -np.infty,        0,        0]  # these bounds are as large as possible
    upper   =[np.infty,   np.infty,  np.infty, np.infty, np.infty]


    def __init__(self, ter=np.nan, A=np.nan, v=np.nan, sv=np.nan, b=np.nan, name='unknown'):
        self.ter=float(ter)
        self.A=float(A)
        self.v=float(v)
        self.sv=float(sv)
        self.b=float(b)
        self.name=name

    def pdf(self, t):
        """LBA PDF for a single accumulator"""
        t=np.maximum(t-self.ter, 1e-5) # absorbed into pdf 
        if self.A<1e-10: # LATER solution
            return np.maximum( 0, (self.b/(t**2)*dnormP(self.b/t, mean=self.v,sd=self.sv))
                              /pnormP(self.v/self.sv) )
        zs=t*self.sv
        zu=t*self.v
        bminuszu=self.b-zu
        bzu=bminuszu/zs
        bzumax=(bminuszu-self.A)/zs
        return np.maximum(0, ((self.v*(pnormP(bzu)-pnormP(bzumax)) +
                        self.sv*(dnormP(bzumax)-dnormP(bzu)))/self.A)/pnormP(self.v/self.sv))
    
    def cdf(self,t):
        """LBA CDF for a single accumulator"""
        t=np.maximum(t-self.ter, 1e-5) # absorbed into cdf         
        if self.A<1e-10: # LATER solution
            return np.minimum(1, np.maximum(0, (pnormP(self.b/t,mean=self.v,sd=self.sv))
                                            /pnormP(self.v/self.sv) ))
        zs=t*self.sv
        zu=t*self.v
        bminuszu=self.b-zu
        xx=bminuszu-self.A
        bzu=bminuszu/zs
        bzumax=xx/zs
        tmp1=zs*(dnormP(bzumax)-dnormP(bzu))
        tmp2=xx*pnormP(bzumax)-bminuszu*pnormP(bzu)
        return np.minimum(np.maximum(0,(1+(tmp1+tmp2)/self.A)/pnormP(self.v/self.sv)),1)

    def sample(self, n, upper_limit=np.infty):
        """draw n random samples from this accumulator's distribution"""
        nnot=n
        rts=[]
        while nnot>0:
            vs=stats.truncnorm.rvs((0-self.v)/float(self.sv), np.infty, loc=self.v, scale=self.sv, size=nnot)
            zs=stats.uniform.rvs(0,self.A,size=nnot)
            crts=self.ter+((self.b - zs)/vs)
            nnot=np.sum(crts>=upper_limit)
            rts+=list(crts[crts<upper_limit])
            
        return np.array(rts, dtype=np.double)
    
if __name__=="__main__":
    import pylab as pl
    # test LBA pdf
    acc=LBAAccumulator(.2, .2, 1.0, 1.0, 1.0)
    print acc
    nsamples=10000
    x=np.linspace(0,10, nsamples)
    pl.plot( x, acc.pdf(x))
    pl.plot( x, acc.cdf(x))
    
    print "PDF Sum: ", (np.sum(acc.pdf(x))/nsamples*10)
    print "PDF Integral: ",scipy.integrate.quad(acc.pdf, 0, np.infty)[0]
    print "CDF(0,1000) :", acc.cdf(np.array([0,1000]))
    pl.show()
