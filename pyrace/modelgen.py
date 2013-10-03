import numpy as np
from numpy import inf
import pandas as pd
import tempfile
from .tools import *

__all__=['ParMap', "ModelTable"]

model_template="""
import pyrace as pr
from pyrace.modelgen import _load_classinst_from_string
from numpy import inf

class {modelname}_paramspec(pr.Parameters):
    parnames=[{parnames}]
    lower   =[{lower}]
    upper   =[{upper}]
    modelstr=""

    def __reduce__(self):
        return (_load_classinst_from_string, (self.modelstr, self.__class__.__name__, self.__dict__))

class {modelname}({parentclass}):

    paramspec={modelname}_paramspec
    modelstr=""

    def __init__(self, pars=None):
        self.design={design}

        if pars!=None:
            if not isinstance(pars, self.paramspec):
                pars=self.paramspec(pars)
            self.set_params(pars)
        else:
            self.set_params(self.__class__.paramspec().random())
        self.set_mixing_probabilities(0,0)

    def copy(self):
        m=self.__class__(self.params)
        m.modelstr=self.modelstr
        return m

    def __reduce__(self):
        return (_load_classinst_from_string, (self.modelstr, self.name(), self.__dict__))

    def set_params(self, pars):
        self.params=pars
        go_acc=[ [None for resp in range(self.design.nresponses())] for cond in range(self.design.nconditions())]
        stop_acc=[ None for cond in range(self.design.nconditions())]

{pardef}

{go_accumulator_definition}

{stop_accumulator_definition}

        self.set_accumulators(go_acc, stop_acc)

        ## mixing probabilities
        pgf=[0 for cond in range(self.design.nconditions())]
        ptf=[0 for cond in range(self.design.nconditions())]

{prob_go_fail_definition}

{prob_trigger_fail_definition}

        self.set_mixing_probabilities(pgf, ptf)
    """


def _load_classinst_from_string(classstr, classname, objdict):
    fname=tempfile.mktemp(suffix='.py')
    with open(fname,'w') as f:
        f.write(classstr)

    ## TODO: ugly hack to enable resetting of model.paramspec to updated paramspec
    ## TODO: should really be removed in future versions
    if not classname.endswith("_paramspec"):
        classobj,co_paramspec=load_class_from_file(fname, [classname, classname+"_paramspec"])
        co_paramspec.modelstr=classstr
        co_paramspec.__module__="__main__"
        classobj.paramspec=co_paramspec
    else:
        classobj=load_class_from_file(fname, classname)
        co_paramspec=None
    if classobj==None:
        raise ValueError("reconstruction of '%s' failed. This is the classstr: %s"%(classname, classstr))

    classobj.modelstr=classstr
    classobj.__module__="__main__"

    inst=classobj()
    inst.__dict__=objdict

    return inst

class ParMap(object):
    def __init__(self, accpar, mapping=None, bounds=None, **kwargs):
        """
        accpar : str
            must be one of the accumulator's parameters or pgf/ptf

        bounds : (float, float)
            specify optimization bounds; if not specified, the maximum
            feasible bounds specified in the definition of each Accumulator
            are used

        mapping : str
            if mapping!=None (a string), it is used to replace the parameter value.
            I.e., normally when using

            >> ModelTable(..., driftrate=ParMap('v')),

            the generated code will have

            >> Accumulator( v=driftrate )

            using

            >> ModelTable(..., driftrate=ParMap('v', mapping='driftrate+A')),

            creates

            >> Accumulator( v=driftrate+A )
        """
        self.accpar=accpar
        self.mapping=mapping
        self.mybounds=bounds
        self.table_index=kwargs

    def bounds(self):
        if self.mybounds==None and self.accpar in ['pgf', 'ptf']:
            return (0,1)
        else:
            return self.mybounds


class ModelTable():
    """Tabular representation of a race-model.

    This should become the way to create a model. Currently not all possible cases are implemented.

    Stuff that is missing:

    DONE: specify constant parameters
    TODO: boundary specification
            * accumulators should propose maximal boundaries
            * ParMap should implement a tighter setting of those boundaries
    DONE: arbitrary parameter-mapping (e.g., map b to A+B)
    TODO: indexing better (currently only indexing of the sort dataframe[column==value] allowed in
          ParMap. Should allow, e.g., dataframe[column %in% list] etc

    """

    def __init__(self, modelname, design, parentcl, fixed={}, **modelspec):
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
        self.fixed=fixed
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

        # go/trigger failures
        pgf=["*" for i in conditions]
        ptf=["*" for i in conditions]

        df=pd.DataFrame.from_items([('condition', conditions),
            ('response', responses), ('gostop',gostop), ('correct', correct)]+factors+accpars+
                                   [('pgf',pgf), ('ptf',ptf)])
        self.table=df
        self.nrows=df.shape[0]

    def init_modelspec(self):
        self.init_empty_table()

        for k,v in self.fixed.items():
            self.table[k]=v
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
        self.check_table()

    def check_table(self):
        accpars=self.parentcl.accumulator_type.parnames
        assert np.all([par in self.table.columns for par in accpars]), 'TableError: not all parameters present'
        assert np.all([fac in self.table.columns for fac in self.design.factors]), 'TableError: not all factors present'
        for fac in self.design.factors:
            levels=self.table[fac].unique()
            if set(levels)!=set(self.design.factor_dict[fac]):
                raise TypeError("column '%s' in table does not contain all levels from %s: %s"%(fac,str(self.design.factor_dict[fac]),str(levels)))
        for par in accpars+['pgf','ptf']:
            assert np.all(np.array(self.table[par])!="*"), "some parameters are not set!\nHere is the model:\n"+str(self)


    def __repr__(self):
        return repr(self.table)

    def _repr_html_(self):
        """used by ipython notebook"""
        r="<table>\n"
        r+="<caption>%s</caption>\n"%self.name
        r+='<tr>\n'
        for col in self.table.columns:
            r+="<TH>%s</TH>"%(col)
        r+="</tr>\n"
        for i in range(self.nrows):
            r+="<tr>"
            r+="".join(["<TD>%s</TD>"%val for val in list(self.table.irow(i))])
            r+="\n</tr>\n"
        r+="</table>\n"
        return r

    def generate_model_str(self):
        self.check_table()

        accpars=self.parentcl.accumulator_type.parnames
        mod_accpars=accpars+["pgf","ptf"]
        modpars=[]
        for accpar in mod_accpars:
            if accpar in self.fixed.keys():
                continue
            modpars += list(self.table[accpar].unique())
        modpars=np.unique(modpars)


        ## go-accumulators
        tpl="""        go_acc[{cond}][{iresp}]=self.accumulator_type({parlist}, name='go-'+':'.join(self.design.condidx({cond})))\n"""
        go_acc_def=""
        for cond in range(self.design.nconditions()):
            for iresp, resp in enumerate(self.design.get_responses()):
                pard={k:self.table[k][(self.table.condition==cond) & (self.table.response==resp)
                                      & (self.table.gostop=='go')].iloc[0] for k in accpars}
                # replace identity mapping with non-trivial mapping
                for k in pard.keys():
                    if k not in self.fixed.keys() and self.modelspec[pard[k]].mapping!=None:
                        pard[k]=self.modelspec[pard[k]].mapping
                parlist=",".join(["%s=%s"%(k,v) for k,v in pard.items()])
                go_acc_def+=tpl.format(cond=cond,resp=resp,iresp=iresp, parlist=parlist)

        ## STOP-accumulators
        tpl="""        stop_acc[{cond}]=self.accumulator_type({parlist}, name='stop-'+':'.join(self.design.condidx({cond})))\n"""
        stop_acc_def=""
        for cond in range(self.design.nconditions()):
            pard={k:self.table[k][(self.table.condition==cond)
                                  & (self.table.gostop=='stop')].iloc[0] for k in accpars}
            # replace identity mapping with non-trivial mapping
            for k in pard.keys():
                if k not in self.fixed.keys() and self.modelspec[pard[k]].mapping!=None:
                    pard[k]=self.modelspec[pard[k]].mapping

            parlist=",".join(["%s=%s"%(k,v) for k,v in pard.items()])
            stop_acc_def+=tpl.format(cond=cond, parlist=parlist)

        ## mixing probabilities
        tpl="""        pgf[{cond}]={pgf}\n"""
        pgf_def=""
        for cond in range(self.design.nconditions()):
            pgf_def+=tpl.format(cond=cond, pgf=self.table['pgf'][self.table.condition==cond].iloc[0])

        tpl="""        ptf[{cond}]={ptf}\n"""
        ptf_def=""
        for cond in range(self.design.nconditions()):
            ptf_def+=tpl.format(cond=cond, ptf=self.table['ptf'][self.table.condition==cond].iloc[0])

        ## bounds
        lower=[]
        upper=[]
        for modpar in modpars:
            bounds=self.modelspec[modpar].bounds()
            if bounds==None:
                bounds=self.parentcl.accumulator_type().get_bounds(self.modelspec[modpar].accpar)
            lower.append(bounds[0])
            upper.append(bounds[1])

        ## build full model
        modelstr=model_template.format(
            modelname=self.name,
            pardef="\n".join(["        %s=pars.%s"%(par,par) for par in modpars]),
            parnames=",".join(['"%s"'%modpar for modpar in modpars]),
            parentclass="pr."+(self.parentcl.__name__.split(".")[-1]),
            lower=",".join([repr(i) for i in lower]),
            upper=",".join([repr(i) for i in upper]),
            design="pr."+repr(self.design),
            go_accumulator_definition=go_acc_def,
            stop_accumulator_definition=stop_acc_def,
            prob_go_fail_definition=pgf_def,
            prob_trigger_fail_definition=ptf_def
        )

        return modelstr

    def generate_model_class(self):
        """
        Return a class-object that can be used to instantiate model objects.

        Note: this implementation uses dynamic importing rather than exec()
              (which, believe me, can make a lot of trouble...)
        """
        modelstr=self.generate_model_str()

        fname=tempfile.mktemp(suffix='.py')
        with open(fname,'w') as f:
            f.write(modelstr)
        (classobj, co_paramspec)=load_class_from_file(fname, [self.name, self.name+"_paramspec"])

        classobj.modelstr=modelstr
        co_paramspec.modelstr=modelstr

        co_paramspec.__module__="__main__"
        classobj.__module__="__main__"
        classobj.paramspec=co_paramspec  ## required!! else we old paramspec
        return classobj

    def generate_model_obj(self, pars=None):
        modcl=self.generate_model_class()
        return modcl(pars=pars)

if __name__=="__main__":
    import pyrace as pr
    import pickle

    factors=[{'deprivation':['control', 'sleep']},
             {'stimulus':['left', 'right']}]
    responses=['left', 'right']
    design=pr.Design(factors, responses, 'stimulus', name='singlego')

    mt=ModelTable('testModelLBA_gf_bounds', design, pr.pSSLBA,
                  fixed={'sv':1.0, 'ptf':0},
                  ter=ParMap('ter', bounds=(0,3)),
                  gf =ParMap('pgf', deprivation='control'),
                  gfs=ParMap('pgf', deprivation='sleep'),
                  Bc =ParMap('b', mapping="Bc +A", deprivation='control', gostop='go', bounds=(0,5) ),
                  Bcs=ParMap('b', mapping="Bcs+As", deprivation='control', gostop='stop', bounds=(0,5) ),
                  Bd =ParMap('b', mapping="Bd +A", deprivation='sleep', gostop='go', bounds=(0,5) ),
                  Bds=ParMap('b', mapping="Bds+As", deprivation='sleep', gostop='stop', bounds=(0,5)),
                  A  =ParMap('A', gostop='go', bounds=(0,5)),
                  As =ParMap('A', gostop='stop', bounds=(0,5)),
                  V  =ParMap('v', correct=True, gostop='go', bounds=(-5,5)),
                  v  =ParMap('v', correct=False, gostop='go', bounds=(-5,5)),
                  Vs =ParMap('v', gostop='stop', bounds=(-5,5)))

    print mt
    modstr=mt.generate_model_str()
    print modstr

    modcl=mt.generate_model_class()
    mod=modcl()#testModel()
    print mod

    with open('test.pickle', 'w') as f:
        pickle.dump(mod, f)
    print mod.__module__

    #mod2 = pickle.loads(pickle.dumps(mod))

    #print mod2
    #print mod2.modelstr
    #import pylab as pl
    #mod.plot_model(lims=(.1,3))
    #pl.show()

#    modcl=generate_model( "TestModel", design, pr.SSWald, df)
#    mod=modcl(design)
#    print mod