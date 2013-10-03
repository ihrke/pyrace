import itertools

class Design:
    def __init__(self, factors, responses, response_to_factor, name='unknown'):
        """
        TODO: Need a mapping from responses to factors using a function instead
              of identity mapping to one of the factors.
        
        factors : list of dicts
        responses : list of possible response; 
        response_to_factor : string matching one of the factors (used to determine 
                     correct response)
        """
        self.name=name
        self.set_factors(factors)
#        if len(responses)<=1:
#            raise ValueError("Can't handle less than 2 responses, currently: responses=%s"%(str(responses)))
        if len(responses)>2:
            print "WARNING: you have a design with more than 2 responses"
            print "  though the code can handle this situation in principle, I've never tested it"
            print "  and there seem to be issues. That said, have fun!"
        self.responses=responses
        self.response_to_factor=response_to_factor
        self.response_to_factor_idx=self.factors.index(response_to_factor)
        for resp in responses:
            if resp not in self.factor_dict[response_to_factor]:
                raise ValueError('response %s not in %s'%(resp, factors[response_to_factor]))
    
    def get_responses(self):
        return self.responses
    
    def set_factors(self, factors):
        self.org_factors=factors
        its=[]
        for fac in factors:
            its+=fac.items()
        self.factor_dict=dict(its)
        self.factors=[v.keys()[0] for v in factors]
        self.factors_idx={ v:i for i,v in enumerate(self.factors)}
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
    
    def correct_response(self, conditionidx, as_index=False):
        """return correct response in condition idx"""
        if as_index:
            return self.responses.index(self.factors_from_int[conditionidx][self.response_to_factor_idx])
        else:
            return self.factors_from_int[conditionidx][self.response_to_factor_idx]

    def __repr__(self):
        return "Design(%s, %s, %s, name=%s)"%(repr(self.org_factors),
                                              repr(self.responses),
                                              repr(self.response_to_factor),
                                              repr(self.name))

    def __str__(self):
        r="<Design '%s'>\n"%self.name
        for k,v in self.factor_dict.items():
            r+=" %s: %s\n"%(str(k),str(v))
        r+=" "+"-"*20+"\n"
        for i in range(self.nconditions()):
            r+= " condition %i <-> %s\n"%(i,self.factors_from_int[i])
        for i in range(self.nresponses()):
            r+= " response  %i <-> %s\n"%(i,self.responses[i])
        return r

    def factorval(self, cond, factor):
        """
        return value of factor in condition cond
        """
        return self.factors_from_int[cond][self.factors_idx[factor]]

    def condidcs_factor(self, factor, value):
        """
        return all condition indices in which factor has value
        """
        idcs=[]
        for k,v in self.factors_from_int.items():
            if value in v[self.factors_idx[factor]]:
               idcs.append(k)
        if len(idcs)==0:
            print "WARNING: no conditions found for %s=%s"%(str(factor), str(value))
        return idcs
    
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

    def __eq__(self, other):
        return self.__dict__==other.__dict__

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
