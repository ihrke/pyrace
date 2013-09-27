import numpy as np
import scipy
import scipy.stats as stats
from scipy.special import erf

from .tools import *
from .racemodel import *

class ShiftedLogNormalAccumulator(Accumulator):
    """
    Shifted LogNormal-Distribution.
    """
    parnames=['ter', 'mu', 'sigma']

    def __init__(self, ter=np.nan, mu=np.nan, sigma=np.nan, name='unknown'):
        """
        Parametrization:

        ter - shift  (must be >0)
        mu  - mean (location) of the underlying normal (in [-infty,infty])
        sigma - sigma of underlying normal (>0)
        """
        self.ter=float(ter)
        self.mu=float(mu)
        self.sigma=float(sigma)
        self.name=name

    def pdf(self, t):
        t=np.maximum(t-self.ter, 1e-5) # absorbed into pdf
        r=1./(t*np.sqrt(2*np.pi)*self.sigma)*np.exp(-(np.log(t)-self.mu)**2/(2*self.sigma**2))
        return np.maximum(0.0, r)
        
    def cdf(self,t):
        t=np.maximum(t-self.ter, 1e-5) # absorbed into cdf
        r=.5+.5*erf((np.log(t)-self.mu)/(np.sqrt(2)*self.sigma))
        return np.minimum( np.maximum( 0., r ), 1.)

    def sample(self, n, upper_limit=np.infty):
        """draw n random samples from this accumulator's distribution"""
        nnot=n
        rts=[]
        while nnot>0:
            crts=self.ter+np.random.lognormal(self.mu, self.sigma, size=nnot)
            nnot=np.sum(crts>=upper_limit)
            rts+=list(crts[crts<upper_limit])
        return np.array(rts, dtype=np.double)
    
if __name__=="__main__":
    """Note: for tesing, run as

        python -m pyrace.slognorm
    """
    import pylab as pl
    # test LBA pdf
    acc=ShiftedLogNormalAccumulator(0, -2, 0.25)
    print acc
    nsamples=10000
    x=np.linspace(0,3, nsamples)
    pl.plot( x, acc.pdf(x))
    pl.plot( x, acc.cdf(x))

    y=acc.sample(10000, upper_limit=3)
    a=pl.hist(y, 100, normed=True)
    print np.sum(a[1])
    
    print "PDF Sum: ", (np.sum(acc.pdf(x))/nsamples*10)
    print "PDF Integral: ",scipy.integrate.quad(acc.pdf, 0, np.infty)[0]
    print "CDF(0,1000) :", acc.cdf(np.array([0,1000]))
    pl.show()
