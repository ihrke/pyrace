import numpy as np
import scipy.stats as stats

list_is_eq=lambda l1,l2: len(l1)==len(l2) and (len(set(l1).intersection(l2)) > 0)

## NOTE: is it a good idea to set a hard numerical threshold (i.e., 7)?
##  should maybe depend on mean/sd? E.g. crit=mean+7*sd ?
def pnormP(x, mean=0, sd=1):
    """standard normal CDF with numerical stability
    R: pnormP  <- function(x,mean=0,sd=1,lower.tail=T){ifelse(abs(x)<7,pnorm(x,mean,sd,lower.tail),ifelse(x<0,0,1))}
    """
    return np.where(np.abs(x-mean)<7.*sd, stats.norm.cdf(x, loc=mean,scale=sd), np.where(x<mean,0,1))
def dnormP(x, mean=0, sd=1):
    """standard normal PDF with numerical stability
    R: dnormP <- function(x,mean=0,sd=1,lower.tail=T){ifelse(abs(x)<7,dnorm(x,mean,sd),0)}
    """
    return np.where(np.abs(x-mean)<7.*sd,stats.norm.pdf(x, loc=mean, scale=sd),0)

def trans_logistic(x, a=0, b=1, inverse=False):
    """ goes from [a,b] to [-inf,+inf] and back;
    inverse=False: [a,b] -> [-inf, +inf]
    inverse=True:  [-inf,+inf] -> [a,b]
    """

    if inverse:
        return (1./(1+np.exp(-x)))*(b-a)+a 
    else:
        # take care of nan
        x=np.where(x<=a, a+1e-12, x)
        x=np.where(x>=b, b-1e-12, x)
        return -np.log( float(b-a)/(x-a)-1)

    
if __name__=="__main__":
    import pylab as pl
    x=np.linspace(-10,10,100)
    stats.norm.pdf(x)
    pl.plot(x,pnormP(x))
    pl.plot(x,dnormP(x))
    pl.show()
