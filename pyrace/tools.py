#from __future__ import print_function
import sys, time
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

class ProgressBar:
    """stolen from somewhere, I guess PyMC?"""
    def __init__(self, iterations):
        self.iterations = iterations
        self.prog_bar = '[]'
        self.fill_char = '*'
        self.width = 50
        self.__update_amount(0)

    def __call__(self, iter):
        self.animate(iter)

    def animate(self, iter):
        print "\r", self,
#        print('\r', self, end='')
        sys.stdout.flush()
        self.update_iteration(iter + 1)

    def update_iteration(self, elapsed_iter):
        self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
        self.prog_bar += '  %d of %s complete' % (elapsed_iter, self.iterations)

    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])

    def __str__(self):
        return str(self.prog_bar)    
    
if __name__=="__main__":
    import pylab as pl
    x=np.linspace(-10,10,100)
    stats.norm.pdf(x)
    pl.plot(x,pnormP(x))
    pl.plot(x,dnormP(x))
    pl.show()
