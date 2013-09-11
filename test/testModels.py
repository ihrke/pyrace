import unittest
import numpy as np
import pandas as pd
import scipy
import sys
from pyrace import Design, StopTaskDataSet
from pyrace import LBAAccumulator, pSSLBA
from pyrace.models.psslba import pSSLBA_modelB


class testpSSLBA_modelB(unittest.TestCase):
    def setUp(self):
        self.factors=[{'TUT':['tut','on-task']},
                      {'stimulus':['left', 'right']}]
        self.responses=['left','right']
        self.design=Design(self.factors,self.responses, 'stimulus')
        self.truepar =pSSLBA_modelB.paramspec(ster=.1, ter=.2, A=.2, Bs=.5,  B=.8, Vs=2.0, Vtut=.5, Vont=1.0, v=0.0)        
        self.mod=pSSLBA_modelB(self.design, self.truepar)

    def test_trans(self):
        for i in range(50):
            rpar=pSSLBA_modelB.paramspec().random()
            assert np.all(np.abs(rpar-self.mod.untrans(self.mod.trans(rpar))))<1e-5

    def test_setparam(self):
        tutcond=[i for i in range(self.design.nconditions()) if self.design.factorval(i,'TUT')=='tut' ]
        ontcond=[i for i in range(self.design.nconditions()) if self.design.factorval(i,'TUT')!='tut' ]
        rpar=self.mod.paramspec().random()
        self.mod.set_params(rpar)

        for tc in tutcond:
            correct=self.design.correct_response(tc)
            for ir,resp in enumerate(self.design.responses):
                if resp==correct:
                    assert self.mod.go_accumulators[tc][ir].v==rpar.Vtut
                else:
                    assert self.mod.go_accumulators[tc][ir].v==rpar.v

        for oc in ontcond:
            correct=self.design.correct_response(oc)
            for ir,resp in enumerate(self.design.responses):
                if resp==correct:
                    assert self.mod.go_accumulators[oc][ir].v==rpar.Vont
                else:
                    assert self.mod.go_accumulators[oc][ir].v==rpar.v
        
if __name__ == '__main__':
    unittest.main()
