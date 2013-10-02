from common import *


class parspec(pr.Parameters):
    parnames=['ter', 'A', 'B', 'V', 'v']
    lower=[0, 0, 0, -5, -5]
    upper=[1, 4, 5, 5, 5]
    
class Testmod(pr.pSSLBA):
    paramspec=parspec;
    def __init__(self,design):
        self.design=design
        self.sv=1.0
        self.set_mixing_probabilities(0,0)
        self.set_params(self.paramspec().random())
    def set_params(self,pars):
        self.params=pars
        go_acc=[]
        stop_acc=[]
        for cond in range(self.design.nconditions()):
            go_acc.append( [pr.LBAAccumulator(pars.ter, pars.A, pars.V, self.sv, pars.A+pars.B)
                            for resp in self.design.get_responses()])
            stop_acc.append(pr.LBAAccumulator(pars.ter, pars.A, pars.V, self.sv, pars.A+pars.B))
        self.set_accumulators(go_acc, stop_acc)



def issorted(x):
    """Check if x is sorted"""
    return (np.diff(x) >= 0).all()

class testParameters(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        class mypar(pr.Parameters):
            parnames=['ich', 'du', 'er']

        pars=mypar(1,2,3, er=5)
        print pars
        assert pars[0]==1
        assert pars[1]==2
        assert pars[2]==5
        assert len(pars)==3
        assert pars.ich==1
        assert pars.du==pars[1]
        assert pars.er==5

        assert tuple(list(pars))==(1,2,5), str(tuple(list(pars)))

    def test_random(self):
        class whatever(pr.Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]

        pars=whatever()            
        for i in range(10):
            pars.random()
            assert pars.a>100 and pars.a<200
            assert pars.b>1 and pars.b<2
            assert pars.c>2 and pars.c<50

    def test_sub(self):
        class whatever(pr.Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]

        pars=whatever()
        pars.random()
        assert np.all(np.abs( (pars-pars) ))<1e-10

    def test_add(self):
        class whatever(pr.Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2.2]
            upper   =[200,   2,  50]

        pars=whatever()
        pars.random()
        assert np.all(np.abs( (pars+pars) - (2*pars) ))<1e-10

    def test_bounds(self):
        class whatever(pr.Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2.2]
            upper   =[200,   2,  50]

        pars=whatever()
        assert pars.bound_lower('a')==100
        assert pars.bound_upper('a')==200
        assert pars.bounds(2)==(2.2, 50)
        

    def test_in_range(self):
        class whatever(pr.Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]
        pars=whatever().random()
        assert pars.in_range()
        out=whatever(a=50, b=1.5, c=5)
        assert not out.in_range()
        assert not pars.in_range(out)
        assert pars.in_range(pars)
    
    def test_abs(self):
        class whatever(pr.Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]

        pars=whatever().random()
        pars2=whatever().random()
        pars3=abs(pars-pars2)
        assert np.all( pars3 )>=0

class testTools(unittest.TestCase):
    def test_trans_logistic(self):
        x=np.linspace(0,1,100)
        y=pr.trans_logistic(x)
        assert issorted(y)
        assert y[0]<0
        assert y[-1]>0
        assert pr.trans_logistic(0.4)<0
        assert pr.trans_logistic(0.6)>0
        assert np.abs(pr.trans_logistic(0.5))<1e-5

    def test_generalized_trans_logistic(self):
        x=np.linspace(-1000,1000,10000)
        y=pr.trans_logistic(x, -1000, +np.infty)
        x2=pr.trans_logistic(y, -1000, +np.infty, inverse=True)
        assert np.all( np.abs(x-x2)<1e-8 )

        x=np.linspace(-1000,1000,10000)
        y= pr.trans_logistic(x, -np.infty, +1000)
        x2=pr.trans_logistic(y, -np.infty, +1000, inverse=True)
        assert np.all( np.abs(x-x2)<1e-8 )

        x=np.linspace(-1000,1000,10000)
        y= pr.trans_logistic(x, -np.infty, +np.infty)
        x2=pr.trans_logistic(y, -np.infty, +np.infty, inverse=True)
        assert np.all( np.abs(x-x2)<1e-8 )


class testPlotting(PlottingEnabledTestCase):
    def test_plot_bar_parameters1(self):
        class whatever(pr.Parameters):
            parnames=['a', 'b', 'c']
            lower   =[100,   1,   2]
            upper   =[200,   2,  50]
        pl.figure()
        a=whatever().random()
        b=whatever().random()
        pr.plot_bar_parameters(a,b,whatever().random())

    def test_plot_bar_parameters(self):
        class whatever2(pr.Parameters):
            parnames=['a', 'b', 'c', 'd', 'ich', 'du']
            lower   =[100,   1,   2,   1,      1,   1]
            upper   =[200,   2,  50,   2,       2,   2]
        npar=10
        pars=[whatever2().random() for i in range(npar)]

        pl.figure()
        pr.plot_bar_parameters(*pars, title='%i tests'%npar)
#        self.savefig('test_plot_bar_parameters_1.png')
        self.savefig()

    def test_model_plot(self):
        design=pr.Design([{'fac1':['level11','level12', 'level13']},
                              {'fac2':['level21','level22', 'level23']},
                              {'fac3':['level31','level32']}],
                              ['level21', 'level22'], 'fac2')
        mod=Testmod(design)
        mod.plot_model()
        self.savefig('test_model_plot_1.png')
        
        
if __name__ == '__main__':
    unittest.main()
