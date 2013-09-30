import numpy as np
import pandas as pd
import pyrace as pr
import itertools

model_template="""
class {modelname}({parentclass}):
    class {modelname}_paramspec(pr.Parameters):
        parnames=[{parnames}]
        lower   =[{lower}]
        upper   =[{upper}]

    paramspec={modelname}_paramspec

    def __init__(self, design, pars=None):
        self.design=design
        if pars!=None:
            self.set_params(pars)
        else:
            self.set_params(self.__class__.paramspec().random())
        self.set_mixing_probabilities(0,0)

    def set_params(self, pars):
        self.params=pars
        go_acc=[ [None for resp in range(self.design.nresponses())] for cond in range(self.design.nconditions())]
        stop_acc=[ None for cond in range(self.design.nconditions())]

{go_accumulator_definition}

{stop_accumulator_definition}

        self.set_accumulators(go_acc, stop_acc)
    """


def generate_model( modname, design, modelclass, model_table, as_string=False):
    accpars=modelclass.accumulator_type.parnames
    assert np.all([par in model_table.columns for par in accpars]), 'TableError: not all parameters present'
    assert np.all([fac in model_table.columns for fac in design.factors]), 'TableError: not all factors present'
    for fac in design.factors:
        levels=model_table[fac].unique()
        if set(levels)!=set(design.factor_dict[fac]):
            raise TypeError("column '%s' in table does not contain all levels from %s: %s"%(fac,str(design.factor_dict[fac]),str(levels)))

    modpars=[]
    for accpar in accpars:
        modpars += list(model_table[accpar].unique())
    modpars=np.unique(modpars)

    tpl="""        go_acc[{cond}][{resp}]=self.accumulator_type({parlist}, name='go-'+':'.join(self.design.condidx({cond})))\n"""
    go_acc_def=""
    for cond,resp in itertools.product( range(design.nconditions()), design.get_responses()):
        print cond, resp
        pard={k:model_table[(model_table.condidx==cond) & (model_table.response==resp)][k][0] for k in accpars}
        print pard
        parlist=",".join(["%s=pars.%s"%(k,v[0]) for k,v in pard.items()])
        go_acc_def+=tpl.format(cond=cond,resp=resp, parlist=parlist)


    modelstr=model_template.format(
        modelname=modname,
        parnames=",".join(['"%s"'%modpar for modpar in modpars]),
        parentclass="pr."+(modelclass.__name__.split(".")[-1]),
        lower=",".join(["0" for i in range(len(modpars))]),
        upper=",".join(["1" for i in range(len(modpars))]),
        go_accumulator_definition=go_acc_def,
        stop_accumulator_definition=""
    )
    loc_dict={}
    print modelstr
    #exec(modelstr, globals(), loc_dict)

    return loc_dict[modname]


if __name__=="__main__":
    model_tab="""
    condidx, deprivation, stimulus, response, theta, alpha, gamma
          0,     normal,     left,     left,   ter,     b,     V
          0,     normal,     left,    right,   ter,     b,     v
          0,     normal,     left,     stop,   ter,     b,    Vs
          1,     normal,    right,     left,   ter,     b,     v
          1,     normal,    right,    right,   ter,     b,     V
          1,     normal,    right,     stop,   ter,     b,    Vs
          2,       sleep,     left,     left,   ter,     b,     V
          2,       sleep,     left,    right,   ter,     b,     v
          2,       sleep,     left,     stop,   ter,     b,    Vs
          3,       sleep,    right,     left,   ter,     b,     v
          3,       sleep,    right,    right,   ter,     b,     V
          3,       sleep,    right,     stop,   ter,     b,    Vs
    """
    import StringIO
    sio=StringIO.StringIO(model_tab.strip())
    #sio.write(model_tab)
    df=pd.read_csv(sio, skipinitialspace=True, )
    #print df

    #mocl2=modelspec("TestModel2", design, pr.SSWald,
    #                V=Par('gamma', varies=['gostop', 'correct']))

    factors=[{'deprivation':['normal', 'sleep']},
             {'stimulus':['left', 'right']}]
    responses=['left', 'right']
    design=pr.Design(factors, responses, 'stimulus', name='singlego')


    modcl=generate_model( "TestModel", design, pr.SSWald, df)
    mod=modcl(design)
    print mod