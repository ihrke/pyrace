from common import *
from pyrace import ParMap

class testModelTable(PlottingEnabledTestCase):
    def test_model1(self):
        factors=[{'deprivation':['normal', 'sleep']},
                 {'stimulus':['left', 'right']}]
        responses=['left', 'right']
        design=pr.Design(factors, responses, 'stimulus', name='singlego')

        mt=pr.ModelTable('testModel', design, pr.SSWald,
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
        
if __name__ == '__main__':
    unittest.main()
