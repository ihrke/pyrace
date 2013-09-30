from common import *
from pyrace import ShiftedWaldAccumulator, VarWaldAccumulator

class testWald(PlottingEnabledTestCase):
    def test_wald_acc(self):
        acc=ShiftedWaldAccumulator(.2, .2, 1.0)
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        self.assertAlmostEqual(scipy.integrate.quad(acc.pdf, 0, np.infty)[0], 1)
        self.assertAlmostEqual(0, acc.cdf(0.01))
        assert abs( acc.cdf(10000)-1)<1e-4

    def test_wald_sample(self):
        acc=ShiftedWaldAccumulator(.2, .2, 2.0)
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


import rpy2.rinterface as rinterface
import rpy2.robjects as robjects

## load the R-code by Logan
r = robjects.r
r['source']('./misc/UWA.R')

class testVarWald(PlottingEnabledTestCase):
    def cmp_pdf_against_R(self, pars):
        """
        check the plots, too
        """
        pars.update({'k':pars['b']-pars['A']/2.0,
                     'l':pars['v'],
                     'a':pars['A']/2.0})
        acc=VarWaldAccumulator( theta=pars['ter'],
                                gamma=pars['v'],
                                alpha=pars['b'],
                                A=pars['A'])
        nsamples=1000
        x=np.linspace(0,3,nsamples)
        rx=rinterface.FloatSexpVector(x - pars['ter'] )

        yr=r['digt'](rx, k=pars['k'], l=pars['l'], a=pars['a'])
        y=np.array(yr)
        assert np.all(y-acc.pdf(x) < 1e-10), str(pars)

        pl.clf()
        pl.plot(x, y, color='red', label='R', linewidth=3)
        pl.plot(x, acc.pdf(x), color='blue', label='python')

        pl.title(str(pars))
        pl.legend()
        self.savefig()

    def cmp_cdf_against_R(self, pars):
        """
        check the plots, too
        """
        pars.update({'k':pars['b']-pars['A']/2.0,
                     'l':pars['v'],
                     'a':pars['A']/2.0})
        acc=VarWaldAccumulator( theta=pars['ter'],
                                gamma=pars['v'],
                                alpha=pars['b'],
                                A=pars['A'])
        nsamples=1000
        x=np.linspace(0,3,nsamples)
        rx=rinterface.FloatSexpVector(x - pars['ter'] )

        yr=r['pigt'](rx, k=pars['k'], l=pars['l'], a=pars['a'])
        y=np.array(yr)
        assert np.all(y-acc.cdf(x) < 1e-10), str(pars)

        pl.clf()
        pl.plot(x, y, color='red', label='R', linewidth=3)
        pl.plot(x, acc.cdf(x), color='blue', label='python')
        pl.title(str(pars))
        pl.legend()
        self.savefig()


    def test_varwald_pdf_against_R(self):
        pars={'A':1.0, 'b':2.0, 'v':.5, 'ter':.1}
        self.cmp_pdf_against_R(pars)

        pars={'A':1.2, 'b':2.0, 'v':.6, 'ter':.5}
        self.cmp_pdf_against_R(pars)

        pars={'A':1.2, 'b':2.0, 'v':0, 'ter':.5}
        self.cmp_pdf_against_R(pars)

    def test_varwald_cdf_against_R(self):
        pars={'A':1.0, 'b':2.0, 'v':.5, 'ter':.1}
        self.cmp_cdf_against_R(pars)

        pars={'A':1.2, 'b':2.0, 'v':.6, 'ter':.5}
        self.cmp_cdf_against_R(pars)

        pars={'A':1.2, 'b':2.0, 'v':0, 'ter':.5}
        self.cmp_cdf_against_R(pars)

    def test_varwald_acc(self):
        pars={'A':1.0, 'b':2.0, 'v':.5, 'ter':.1}
        acc=VarWaldAccumulator( theta=pars['ter'],
                                gamma=pars['v'],
                                alpha=pars['b'],
                                A=pars['A'])
        nsamples=10000
        x=np.linspace(0,10, nsamples)
        self.assertAlmostEqual(scipy.integrate.quad(acc.pdf, 0, np.infty)[0], 1)
        self.assertAlmostEqual(0, acc.cdf(0.01))
        assert abs( acc.cdf(10000)-1)<1e-4

    def test_varwald_sample(self):
        pars={'A':1.0, 'b':2.0, 'v':.5, 'ter':.1}
        acc=VarWaldAccumulator( theta=pars['ter'],
                                gamma=pars['v'],
                                alpha=pars['b'],
                                A=pars['A'])
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
            pl.title(str(pars))
            #pl.xlim(0,3)

            #pl.subplot(2,1,2)
            pl.plot(x,acc.pdf(x), color='red', label='analytical')
            #pl.plot(x,dens(x),    color='green', label='kde')
            pl.legend()
            self.savefig()

if __name__ == '__main__':
    unittest.main()
