import numpy as np
import math
import pyrace

class Rastrigin:
    """Defines the Rastrigin benchmark problem.
    
    This class defines the Rastrigin global optimization problem. This 
    is a highly multimodal minimization problem where the local minima
    are regularly distributed. It is defined as follows:
    
    .. math::
    
        f(x) = \sum_{i=1}^n (x_i^2 - 10\cos(2\pi x_i) + 10)
    
    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-5.12, 5.12]` for :math:`i=1,...,n`.
    
    .. figure:: _static/image12271.jpg
        :alt: Rastrigin function
        :align: center
        
        Two-dimensional Rastrigin function 
        (`image source <http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2607.htm>`__)
    
    Public Attributes:
    
    - *global_optimum* -- the problem input that produces the optimum output.
      Here, this corresponds to [0, 0, ..., 0].
    
    """
    def __init__(self, dimensions=2):
        self.dimensions=dimensions
        self.global_optimum = [0 for _ in range(self.dimensions)]
        
    def generator(self, npop):
        return [np.random.uniform(-5.12, 5.12, self.dimensions) for _ in range(npop)]
        
    def evaluator(self, cand):
        fitness=(sum([x**2 - 10 * math.cos(2 * math.pi * x) + 10 for x in cand]))
        print "Rastr(",cand,")=",fitness        
        return fitness

    def __call__(self, cand):
        return self.evaluator(cand)


npop=10
prob=Rastrigin(5)
pop_ini=prob.generator(npop)

final=pyrace.de_rand_1_bin( pop_ini, prob, gen_max=10)
print final


