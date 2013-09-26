from common import *
from pyrace import ShiftedLogNormalAccumulator

class testLogNorm(PlottingEnabledTestCase):
    def test_lognorm_acc(self):
        acc=ShiftedLogNormalAccumulator(.2, .2, .1)
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        self.assertAlmostEqual(scipy.integrate.quad(acc.pdf, 0, np.infty)[0], 1)
        self.assertAlmostEqual(0, acc.cdf(0.01))
        assert abs( acc.cdf(10000)-1)<1e-4

    def test_sample(self):
        acc=ShiftedLogNormalAccumulator(.2, .2, .20)
        nsamples=100000
        x=np.linspace(0,10, nsamples)
        
        import pylab as pl
        samp=acc.sample(nsamples)
        #dens=scipy.stats.gaussian_kde(samp[samp<10])
        #pl.hist(acc.sample(nsamples),200, normed=True)
        h,hx=np.histogram(samp, density=True, bins=1000)
        hx=hx[:-1]+(hx[1]-hx[0])/2.
        assert np.any(np.abs(h-acc.pdf(hx))<0.5)
        
        if True:
            #pl.subplot(2,1,1)
            pl.hist(samp[samp<10],300, normed=True, alpha=.3)
            pl.xlim(0,3)

            #pl.subplot(2,1,2)                
            pl.plot(x,acc.pdf(x), color='red', label='analytical')
            #pl.plot(x,dens(x),    color='green', label='kde')
            pl.legend()
            self.savefig()
        
if __name__ == '__main__':
    unittest.main()
