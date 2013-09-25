import numpy as np
import scipy
import scipy.stats as stats

from .tools import *
from .racemodel import *

class ShiftedWaldAccumulator(Accumulator):
    """
    Shifted Wald distribution (see Watzke and Wagenmakers, 2009).
    
    Note: forget scipy.stats.wald or scipy.stats.invgauss
    """
    def __init__(self, alpha, theta, gamma, name='unknown'):
        """
        Shifted Wald distribution (see Watzke and Wagenmakers, 2009).

        Parametrization according to their paper:

        alpha - boundary
        gamma - drift-rate
        theta - non-decision-time

        the Wikipedia-parameters correspond to those as follows:

        mu = alpha/gamma
        lambda=alpha^2
        theta is just a shift (replace x with x-theta)
        """
        self.parnames=['alpha', 'theta', 'gamma']
        self.alpha=float(alpha)
        self.theta=float(theta)
        self.gamma=float(gamma)
        self.mu=self.alpha/self.gamma  # this is a reparametrisation 
        self.lambd=self.alpha**2
        self.name=name

    def pdf(self, t):
        t=np.maximum(t-self.theta, 1e-5) # absorbed into pdf
        r=self.alpha/(np.sqrt(2*np.pi*(t**3)))*np.exp(- ((self.alpha-self.gamma*t)**2)/(2*t))
        return np.maximum(0.0, r)
        
    def cdf(self,t):
        t=np.maximum(t-self.theta, 1e-5) # absorbed into cdf
        r=pnormP( (self.gamma*t-self.alpha)/np.sqrt(t))+np.exp(2*self.alpha*self.gamma)*pnormP(-(self.gamma*t+self.alpha)/(np.sqrt(t)))
        return np.minimum( np.maximum( 0., r ), 1.)

    def sample(self, n, upper_limit=np.infty):
        """draw n random samples from this accumulator's distribution"""
        nnot=n
        rts=[]
        while nnot>0:
            crts=self.theta+np.random.wald(self.mu, self.lambd, size=nnot)
            nnot=np.sum(crts>=upper_limit)
            rts+=list(crts[crts<upper_limit])
        return np.array(rts, dtype=np.double)
    
if __name__=="__main__":
    import pylab as pl
    # test LBA pdf
    acc=ShiftedWaldAccumulator(.2, .2, 1)
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
