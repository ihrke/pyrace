import numpy as np
from tools import *
from racemodel import *
import scipy

class LBAAccumulator(Accumulator):
    def __init__(self, ter, A, v, sv, b, name='unknown'):
        self.parnames=['ter','A','v','sv','b']
        self.ter=ter
        self.A=A
        self.v=v
        self.sv=sv
        self.b=b
        self.name=name

    def pdf(self, t):
        """LBA PDF for a single accumulator"""
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
