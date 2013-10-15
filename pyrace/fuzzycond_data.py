"""
Fuzzy-condition dataset.

Highly experimental and undocumented.
"""
import pandas as pd
import numpy as np
from .data import StopTaskDataSet
from .design import Design


class FuzzyConditionStopTaskDataSet(StopTaskDataSet):
    def __init__(self, design, data=None, name=None, mapping=None):
        """
        Implements a data-set for the stop-task with fuzzy conditions.

        Internally, it stores the data using the following structure:

        condition SSD RT response

        where condition is an a list of probabilities that this trial corresponds
        to this index' condition.

        * if SSD==np.nan -> GO-trial
        * if RT==np.nan -> missed/successful STOP
        * if response==-1 -> miss/successful STOP
        * if condition<0 -> something is wrong with the dataset

        # TODO: implement from_dataset()

        **parameters**:

        design : Design
           description of the experiment
        data : pandas.DataFrame or dict
           with columns corresponding to the factors in design (=wide format)
           or with a single 'condition' column indexing design.conditions,
           or dict with fields 'condition', 'RT', 'SSD', 'response' (or a mapping thereof)
        mapping : dict
           must contain keys "RT", "SSD", 'response' and use the values to index into the
           pandas dataframe.
        """
        self.mapping={'RT':'RT', 'SSD':'SSD', 'response':'response','condition':'condition'}
        if mapping!=None:
            self.mapping.update(mapping)
        self.design=design
        self.name=name
        if isinstance(data, dict) or format=='dict':
            self.from_dict(data)
        else:
            raise NotImplementedError('data must be dict')

        if not self.check():
            print "WARNING: dataset '%s' is not consistent! Run dataset.check(verbose=True) to check!"%self.name

    def from_dict(self, data):
        self.org_data=data
        self.ntrials=len(data[self.mapping['RT']])

        # RT
        self.RT=np.array( data[self.mapping['RT']], dtype=np.float)
        self.RT[np.logical_not(np.isfinite(self.RT))]=np.nan

        # SSD
        self.SSD=np.array( data[self.mapping['SSD']], dtype=np.float)
        self.SSD[np.logical_not(np.isfinite(self.SSD))]=np.nan

        # responses
        self.response=np.array(data[self.mapping['response']], dtype=np.int)
        if np.any(self.response>=self.design.nresponses()):
            raise ValueError

        # conditions
        self.condition=np.array(data[self.mapping['condition']])
        if self.condition.shape!=(self.ntrials, self.design.nconditions()):
            raise ValueError

    def define_checks(self):
        checks=['r=len(self.RT)==self.ntrials',
            'r=len(self.SSD)==self.ntrials',
            'r=len(self.response)==self.ntrials',
            'r=len(self.condition)==self.ntrials',
            ('GO(NA)==responses<0', 'ix=np.isnan(self.SSD); r=np.all(np.where(np.isnan(self.RT[ix]))[0]==np.where(self.response[ix]<0)[0])'),
            'r=self.RT.dtype==np.float',
            'r=self.SSD.dtype==np.float',
            'r=self.response.dtype==np.int',
            'r=self.condition.shape==(self.ntrials,self.design.nconditions())',
            ('fuzzy sum to one?', 'r=np.all(np.sum(self.condition, axis=1)==1)'),
            ('RT  : any stop trials?', 'r=np.any(np.isnan(self.RT))'),
            ('SSD : any stop trials?', 'r=np.any(np.isnan(self.SSD))'),
            ('resp: any stop trials?', 'r=np.any(self.response<0)'),
            ('RT(nan) wherevever resp<0?', 'r=np.all( np.where(np.isnan(self.RT))[0]==np.where(self.response<0)[0])'),
            ]
        return checks

    def from_pandas_dataframe( self, data, format="wide"):
        raise NotImplementedError

    def as_dataframe(self, form='long', conditions_expanded=True):
        if form=='long':
            raise NotImplementedError
        elif form=='wide':
            df=pd.DataFrame({'SSD':self.SSD,
                             'RT':self.RT,
                             'response':self.response})
            for cond in range(self.design.nconditions()):
                col=":".join(self.design.condidx(cond))
                df[col]=self.condition[:,cond]
        else:
            raise ValueError("don't know how to handle format %s"%form)
        return df


if __name__=="__main__":
    import numpy as np
    factors=[{'tut':['on','off']},
             {'stimulus':['go']}]
    responses=['go']
    design=Design(factors,responses, 'stimulus')
    ntrials=100
    dat=dict(RT=np.random.rand(ntrials),
             SSD=np.random.rand(ntrials),
             response=[0 for _ in range(ntrials)],
             condition=[[0.3, 0.7] for _ in range(ntrials)])

    ds=FuzzyConditionStopTaskDataSet(design,dat)
    ds.check(verbose=True)

    print ds.head()