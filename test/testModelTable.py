from common import *
from pyrace import ParMap

class testModelTable(PlottingEnabledTestCase):
    def test_model1(self):
        factors=[{'deprivation':['normal', 'sleep']},
                 {'stimulus':['left', 'right']}]
        responses=['left', 'right']
        design=pr.Design(factors, responses, 'stimulus', name='singlego')

        mt=pr.ModelTable('testModel', design, pr.SSWald,
                      fixed={'pgf':0, 'ptf':0},
                      ter=ParMap('theta'),
                      b=ParMap('alpha'),
                      V  =ParMap('gamma', correct=True, gostop='go'),
                      v  =ParMap('gamma', correct=False, gostop='go'),
                      Vs =ParMap('gamma', gostop='stop'))

        print mt
        modstr=mt.generate_model_str()
        print modstr
        loc={}
        exec( modstr, globals(), loc)
        modcl=loc['testModel']
        mod=modcl()

        pl.clf()
        mod.plot_model(lims=(.1,3))
        self.savefig()

    def test_model2(self):
        factors=[{'deprivation':['normal', 'sleep']},
                 {'stimulus':['left', 'right']}]
        responses=['left', 'right']
        design=pr.Design(factors, responses, 'stimulus', name='singlego')

        mt=pr.ModelTable('testModelLBA', design, pr.pSSLBA,
                      fixed={'sv':1.0, 'pgf':0, 'ptf':0},
                      ter=ParMap('ter'),
                      b  =ParMap('b'),
                      A  =ParMap('A'),
                      V  =ParMap('v', correct=True, gostop='go'),
                      v  =ParMap('v', correct=False, gostop='go'),
                      Vs =ParMap('v', gostop='stop'))

        print mt
        modstr=mt.generate_model_str()
        print modstr
        loc={}
        exec( modstr, globals(), loc)
        modcl=loc['testModelLBA']
        mod=modcl()

        pl.clf()
        mod.plot_model(lims=(.1,3))
        self.savefig()

    def test_model3(self):
        factors=[{'deprivation':['control', 'sleep']},
                 {'stimulus':['left', 'right']}]
        responses=['left', 'right']
        design=pr.Design(factors, responses, 'stimulus', name='singlego')

        mt=pr.ModelTable('testModelLBAmapping', design, pr.pSSLBA,
                      fixed={'sv':1, 'pgf':0, 'ptf':0},
                      ter=ParMap('ter'),
                      Bc =ParMap('b', mapping="Bc +A", deprivation='control', gostop='go'  ),
                      Bcs=ParMap('b', mapping="Bcs+As", deprivation='control', gostop='stop'),
                      Bd =ParMap('b', mapping="Bd +A", deprivation='sleep', gostop='go'  ),
                      Bds=ParMap('b', mapping="Bds+As", deprivation='sleep', gostop='stop'),
                      A  =ParMap('A', gostop='go'),
                      As =ParMap('A', gostop='stop'),
                      V  =ParMap('v', correct=True, gostop='go'),
                      v  =ParMap('v', correct=False, gostop='go'),
                      Vs =ParMap('v', gostop='stop'))

        print mt
        modstr=mt.generate_model_str()
        print modstr
        loc={}
        exec( modstr, globals(), loc)
        modcl=loc['testModelLBAmapping']
        mod=modcl()

        pl.clf()
        mod.plot_model(lims=(.1,3))
        self.savefig()

    def test_model4(self):
        factors=[{'deprivation':['control', 'sleep']},
                 {'stimulus':['left', 'right']}]
        responses=['left', 'right']
        design=pr.Design(factors, responses, 'stimulus', name='singlego')

        mt=pr.ModelTable('testModelLBA_gf', design, pr.pSSLBA,
                      fixed={'sv':1.0, 'ptf':0},
                      ter=ParMap('ter'),
                      gf =ParMap('pgf', deprivation='control'),
                      gfs=ParMap('pgf', deprivation='sleep'),
                      Bc =ParMap('b', mapping="Bc +A", deprivation='control', gostop='go'  ),
                      Bcs=ParMap('b', mapping="Bcs+As", deprivation='control', gostop='stop'),
                      Bd =ParMap('b', mapping="Bd +A", deprivation='sleep', gostop='go'  ),
                      Bds=ParMap('b', mapping="Bds+As", deprivation='sleep', gostop='stop'),
                      A  =ParMap('A', gostop='go'),
                      As =ParMap('A', gostop='stop'),
                      V  =ParMap('v', correct=True, gostop='go'),
                      v  =ParMap('v', correct=False, gostop='go'),
                      Vs =ParMap('v', gostop='stop'))

        print mt
        modstr=mt.generate_model_str()
        print modstr
        loc={}
        exec( modstr, globals(), loc)
        modcl=loc['testModelLBA_gf']
        mod=modcl()

        pl.clf()
        mod.plot_model(lims=(.1,3))
        self.savefig()

if __name__ == '__main__':
    unittest.main()
