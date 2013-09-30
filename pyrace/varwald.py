import numpy as np
import scipy
import scipy.stats as stats

from .tools import *
from .racemodel import *

class VarWaldAccumulator(Accumulator):
    """
    Wald Distribution with starting point variability
    """
    parnames=['alpha', 'theta', 'gamma', 'A']
    def __init__(self, alpha=np.nan, theta=np.nan, gamma=np.nan, A=np.nan, name='unknown'):
        """
        Wald Distribution with starting point variability.
        See Psychological Review. Logan, Van Zandt, Verbruggen, and Wagenmakers (in press).
        On the ability to inhibit thought and action: General and special theories
        of an act of control.


        Parametrization:

        alpha - boundary   (=b)
        gamma - drift-rate (=v)
        theta - non-decision-time (ter)
        A     - starting-point variability (uniform)

        the Wikipedia-parameters correspond to those as follows:

        mu = alpha/gamma
        lambda=alpha^2
        theta is just a shift (replace x with x-theta)

        The reparametrization to fit Logan's 2013 code is as follows:

        a=A/2
        k=alpha-A/2
        l=gamma
        """
        self.alpha=float(alpha)
        self.theta=float(theta)
        self.gamma=float(gamma)
        self.A=float(A)
        if self.gamma<1e-10:
            self.mu=self.alpha/1e-10
        else:
            self.mu=self.alpha/self.gamma  # this is a reparametrisation

        if self.alpha<self.A:
            raise ValueError("starting-point var exceeds threshold, b=%f, A=%f"%(self.alpha,self.A))
        self.lambd=self.alpha**2
        self.name=name

    def pdf(self, t):
        """
        Implementing the code distributed with Logan et al. 2013 using
        the reparametrization given above.

        Also, the constant theta is added as usual.
        """
        t=np.maximum(t-self.theta, 1e-5) # absorbed into pdf
        sqrt_t=np.sqrt(t)

        # reparametrization
        a=self.A/2.0
        k=self.alpha-self.A/2.0
        l=self.gamma

        if self.A<1e-10: # this is the solution without starting-point variability
            r=self.alpha/(np.sqrt(2*np.pi*(t**3)))*np.exp(- ((self.alpha-self.gamma*t)**2)/(2*t))
        elif self.gamma<1e-10:
            r=np.exp( -.5*( np.log(2)+np.log(np.pi)+np.log(t))
                      + np.log( np.exp(-( (k-a)**2/(2*t)))-np.exp(-( (k+a)**2/(2*t) )) )
                      - np.log(2) - np.log(a) )
        else:
            r=np.exp( np.log( (np.exp(- (a-k+t*l)**2/(2*t) )-np.exp(- (a+k-t*l)**2/(2*t) ))/np.sqrt(2*np.pi*t)
                              + np.exp(np.log(.5)+np.log(l))*( 2*pnormP( (-k+a)/sqrt_t + sqrt_t*l)-1
                                                               + 2*pnormP( (k+a)/sqrt_t - sqrt_t*l)-1) )
                      - np.log(2) - np.log(a))

        return np.maximum(0.0, np.where( np.isnan(r), 0, r))
        
    def cdf(self,t):
        t=np.maximum(t-self.theta, 1e-5) # absorbed into cdf

        sqrt_t=np.sqrt(t)

        # reparametrization
        a=self.A/2.0
        k=self.alpha-self.A/2.0
        l=self.gamma

        if self.A<1e-10: # this is the solution without starting-point variability
            r=pnormP( (self.gamma*t-self.alpha)/sqrt_t)+np.exp(2*self.alpha*self.gamma)*pnormP(-(self.gamma*t+self.alpha)/(sqrt_t))
        elif self.gamma<1e-10:
            r=(( -(k+a)*(2*pnormP( (k+a)/sqrt_t)-1)
                 -(k-a)*(2*pnormP(-(k-a)/sqrt_t)-1))/(2*a)) \
              + (1 + np.exp(-.5*(k-a)**2/t - .5*np.log(2) - .5*np.log(np.pi) + .5*np.log(t) - np.log(a))
                 -   np.exp(-.5*(k+a)**2/t - .5*np.log(2) - .5*np.log(np.pi) + .5*np.log(t) - np.log(a)))
        else:
            t1=np.exp( .5*np.log(t)-.5*np.log(2*np.pi) ) * (  np.exp( -((k-a-t*l)**2/t)/2.)
                                                            - np.exp( -((k+a-t*l)**2/t)/2.) ) # ok
            t2=a+(   np.exp(2*l*(k+a)+np.log(pnormP(-(k+a+t*l)/sqrt_t)))
                   - np.exp(2*l*(k-a)+np.log(pnormP(-(k-a+t*l)/sqrt_t))) )/(2*l) # ok
            t4= (.5*(t*l-a-k+.5/l)) * ( 2*pnormP((k+a)/sqrt_t-sqrt_t*l)-1) \
               + .5*(k-a-t*l-.5/l)*( 2*pnormP((k-a)/sqrt_t-sqrt_t*l)-1)
            r=(t4+t2+t1)/(2*a)


        return np.minimum( np.maximum( 0., np.where( np.isnan(r), 0, r) ), 1.)

    def sample(self, n, upper_limit=np.infty):
        """draw n random samples from this accumulator's distribution"""
        nnot=n
        rts=[]
        while nnot>0:
            crts=self.theta+np.random.wald(self.mu, self.lambd, size=nnot)-np.random.uniform(0,self.A,size=nnot)
            nnot=np.sum(crts>=upper_limit)
            rts+=list(crts[crts<upper_limit])
        return np.array(rts, dtype=np.double)
    

import rpy2.robjects.numpy2ri
import rpy2.rinterface as rinterface

if __name__=="__main__":
    import pylab as pl
    # test LBA pdf

    acc=VarWaldAccumulator(2, 0.1, 0, 1)
    print acc
    nsamples=10000

    x=np.linspace(0,3, nsamples)
    pl.plot( x, acc.pdf(x))
    pl.xlim(0,3)


    import rpy2.robjects as robjects

    r = robjects.r
    r['source']('./test/misc/UWA.R')


    rx = rinterface.FloatSexpVector(x-.1 )

    yr=r['digt'](rx, k=2-.5, l=0, a=.5)
    y=np.array(yr)
    pl.plot( x, y, color='red')
    assert np.any(y-acc.pdf(x) < 1e-10)
    #pl.plot( x, acc.cdf(x))
    


    #y=acc.sample(10000, upper_limit=3)
    #a=pl.hist(y, 100, normed=True)
    #print np.sum(a[1])
    
    #print "PDF Sum: ", (np.sum(acc.pdf(x))/nsamples*10)
    #print "PDF Integral: ",scipy.integrate.quad(acc.pdf, 0, np.infty)[0]
    #print "CDF(0,1000) :", acc.cdf(np.array([0,1000]))
    pl.show()
