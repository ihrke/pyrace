import numpy as np
import pandas as pd
from .tools import flatten

model_template="""
import pyrace as pr
class {modelname}({parentclass}):
    class {modelname}_paramspec(pr.Parameters):
        parnames=[{parnames}]
        lower   =[{lower}]
        upper   =[{upper}]

    paramspec={modelname}_paramspec

    def __init__(self, pars=None):
        self.design={design}
        if pars!=None:
            self.set_params(pars)
        else:
            self.set_params(self.__class__.paramspec().random())
        self.set_mixing_probabilities(0,0)

    def set_params(self, pars):
        self.params=pars
        go_acc=[ [None for resp in range(self.design.nresponses())] for cond in range(self.design.nconditions())]
        stop_acc=[ None for cond in range(self.design.nconditions())]

{pardef}

{go_accumulator_definition}

{stop_accumulator_definition}

        self.set_accumulators(go_acc, stop_acc)
    """


class ParMap(object):
    def __init__(self, accpar, **kwargs):
        self.accpar=accpar
        self.table_index=kwargs


class ModelTable():
    """Tabular representation of a race-model"""

    def __init__(self, modelname, design, parentcl, **modelspec):
        """
        modelname : str
            just a name for the model

        design : pyrace.Design

        parentcl : class
            the parentclass for the model (e.g., pyrace.SSVarWald)

        modelspec : kw-args
            dict of ParMap objects specifying the relationship between key
            and entry in table
        """
        self.name=modelname
        self.design=design
        self.parentcl=parentcl
        self.modelspec=modelspec
        self.init_modelspec()

    def init_empty_table(self):
        conditions=flatten([[i]*(self.design.nresponses()+1) for i in range(self.design.nconditions())])
        responses =flatten([self.design.get_responses()+['stop'] for i in range(self.design.nconditions())])
        correct   =[self.design.correct_response(cond)==resp for cond,resp in zip(conditions,responses)]
        gostop    =['go' if resp!='stop' else 'stop' for resp in responses]

        # factor columns in correct order
        factors=[(k,[]) for k in self.design.factors]
        for cond in conditions:
            condl=self.design.factors_from_int[cond]
            for i,(fac,val) in enumerate(zip(self.design.factors, condl)):
                factors[i][1].append(val)

        # accumulator columns in good order
        accpars=[(k,'*') for k in self.parentcl.accumulator_type.parnames]

        df=pd.DataFrame.from_items([('condition', conditions),
            ('response', responses), ('gostop',gostop), ('correct', correct)]+factors+accpars)
        self.table=df
        self.nrows=df.shape[0]

    def init_modelspec(self):
        self.init_empty_table()

        for parname, spec in self.modelspec.items():
            if spec.accpar not in self.table.columns:
                raise Exception("Accumulator parameter '%s' for Accumulator '%s' unknown"%(str(spec.accpar),
                                        str(self.parentcl.accumulator_type.__name__)))
            ind=np.ones(self.nrows, dtype=np.bool)*True
            for col,ix in spec.table_index.items():
                try:
                    self.table[col]==ix
                except:
                    raise Exception('index error: self.table[%s]==%s'%(str(col), str(ix)))
                ind=np.logical_and(ind, self.table[col]==ix)
            if np.any( np.array(self.table[spec.accpar][ind])!="*"):
                print "WARNING: overwriting previous specifications; Model may be misspecified."
                print "   Check resulting table wether the model does what you want!"
            self.table[spec.accpar][ind]=parname

    def check_table(self):
        accpars=self.parentcl.accumulator_type.parnames
        assert np.all([par in self.table.columns for par in accpars]), 'TableError: not all parameters present'
        assert np.all([fac in self.table.columns for fac in self.design.factors]), 'TableError: not all factors present'
        for fac in self.design.factors:
            levels=self.table[fac].unique()
            if set(levels)!=set(self.design.factor_dict[fac]):
                raise TypeError("column '%s' in table does not contain all levels from %s: %s"%(fac,str(self.design.factor_dict[fac]),str(levels)))


    def __repr__(self):
        return repr(self.table)

    def generate_model_str(self):
        self.check_table()

        accpars=self.parentcl.accumulator_type.parnames
        modpars=[]
        for accpar in accpars:
            modpars += list(self.table[accpar].unique())
        modpars=np.unique(modpars)

        tpl="""        go_acc[{cond}][{iresp}]=self.accumulator_type({parlist}, name='go-'+':'.join(self.design.condidx({cond})))\n"""
        go_acc_def=""
        for cond in range(self.design.nconditions()):
            for iresp, resp in enumerate(self.design.get_responses()):
                pard={k:self.table[k][(self.table.condition==cond) & (self.table.response==resp)
                                      & (self.table.gostop=='go')].iloc[0] for k in accpars}
                parlist=",".join(["%s=%s"%(k,v) for k,v in pard.items()])
                go_acc_def+=tpl.format(cond=cond,resp=resp,iresp=iresp, parlist=parlist)

        tpl="""        stop_acc[{cond}]=self.accumulator_type({parlist}, name='stop-'+':'.join(self.design.condidx({cond})))\n"""
        stop_acc_def=""
        for cond in range(self.design.nconditions()):
            pard={k:self.table[k][(self.table.condition==cond)
                                  & (self.table.gostop=='stop')].iloc[0] for k in accpars}
            parlist=",".join(["%s=%s"%(k,v) for k,v in pard.items()])
            stop_acc_def+=tpl.format(cond=cond, parlist=parlist)


        modelstr=model_template.format(
            modelname=self.name,
            pardef="\n".join(["        %s=pars.%s"%(par,par) for par in modpars]),
            parnames=",".join(['"%s"'%modpar for modpar in modpars]),
            parentclass="pr."+(self.parentcl.__name__.split(".")[-1]),
            lower=",".join(["0" for i in range(len(modpars))]),
            upper=",".join(["1" for i in range(len(modpars))]),
            design="pr."+repr(self.design),
            go_accumulator_definition=go_acc_def,
            stop_accumulator_definition=stop_acc_def
        )
        #loc_dict={}
        #print modelstr
        #exec(modelstr, globals(), loc_dict)
        #return loc_dict[modname]

        return modelstr



if __name__=="__main__":
    import pyrace as pr

    factors=[{'deprivation':['normal', 'sleep']},
             {'stimulus':['left', 'right']}]
    responses=['left', 'right']
    design=pr.Design(factors, responses, 'stimulus', name='singlego')

    mt=ModelTable('testModel', design, pr.SSWald,
                  ter=ParMap('theta'), b=ParMap('alpha'),
                  V  =ParMap('gamma', correct=True, gostop='go'),
                  v  =ParMap('gamma', correct=False, gostop='go'),
                  Vs =ParMap('gamma', gostop='stop'))

    print mt
    modstr=mt.generate_model_str()
    print modstr
    exec(modstr)
    print testModel()

#    modcl=generate_model( "TestModel", design, pr.SSWald, df)
#    mod=modcl(design)
#    print mod