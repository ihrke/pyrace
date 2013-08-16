import itertools

class Design:
    def __init__(self, factors, responses, response_to_factor):
        """
        factors : list of dicts
        responses : list of possible response; 
        response_to_factor : string matching one of the factors (used to determine 
                     correct response)
        """
        self.set_factors(factors)
        self.responses=responses
        self.response_to_factor=response_to_factor
        self.response_to_factor_idx=self.factors.index(response_to_factor)
        for resp in responses:
            if resp not in self.factor_dict[response_to_factor]:
                raise ValueError('response %s not in %s'%(resp, factors[response_to_factor]))
    
    def get_responses(self):
        return self.responses
    
    def set_factors(self, factors):
        its=[]
        for fac in factors:
            its+=fac.items()
        self.factor_dict=dict(its)
        self.factors=[v.keys()[0] for v in factors]
        l=list(itertools.product( *[v.values()[0] for v in factors]))
        self.factors_to_int={ ":".join([str(e) for e in v]):i for i,v in enumerate(l)}
        self.factors_from_int={ i:[e for e in v] for i,v in enumerate(l)}
    
    def nresponses(self):
        return len(self.responses)
    
    def naccumulators(self):
        """return number of accumulators in design (only GO)"""
        return len(self.factors_to_int)*len(self.responses)
    
    def nconditions(self):
        return len(self.factors_to_int)
    
    def correct_response(self, conditionidx):
        """return correct response in condition idx"""
        return self.factors_from_int[conditionidx][self.response_to_factor_idx]
    
    def __repr__(self):
        r="<Design>\n"
        for k,v in self.factor_dict.items():
            r+=" %s: %s\n"%(str(k),str(v))
        r+=" "+"-"*20+"\n"
        for i in range(self.nconditions()):
            r+= " %i <-> %s\n"%(i,self.factors_from_int[i])
        return r
    
    def condidx(self, condition):
        """for an index, return condition string,
           for a list/string return condition index"""
        if isinstance(condition,list):
            return self.factors_to_int[':'.join(condition)]
        elif isinstance(condition,str):
            return self.factors_to_int[condition]
        elif isinstance(condition,int):
            return self.factors_from_int[condition]
        else:
            raise TypeError("don't know what to do with condition %s"%str(condition))


if __name__=="__main__":
    factors=[{'sleepdep':['normal','deprived']},
        #         {'TUT':['on-task', 'tut']},
        {'stimulus':['left', 'right']}]
    responses=['left','right']
    design=Design(factors,responses, 'stimulus')
    print design.factors_to_int
    k=['deprived','left']
    print design.condidx(  k )
    print design.factor_dict
    print design
