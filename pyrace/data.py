from tools import *
from design import Design
import pandas as pd
import numpy as np

class StopTaskDataSet(object):
    """
    Implements a data-set for the stop-task.

    
    
    Internally, it stores the data using the following structure:
    
    condition SSD RT response
    
    where condition is an integer number that is expanded by 
    the Design() as is response. 
    
    * if SSD==np.nan -> GO-trial
    * if RT==np.nan -> missed/successful STOP
    * if response==-1 -> miss/successful STOP
    * if condition<0 -> something is wrong with the dataset
    """
    def __init__(self, design, data=None, format='wide', mapping={'RT':'RT', 'SSD':'SSD', 'response':'response','condition':'condition'}):
        """
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
        self.mapping=mapping
        self.design=design
        if isinstance(data, dict) or format=='dict':
            self.from_dict(data)
        else:
            self.from_pandas_dataframe(data, format=format)

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
        self.condition=np.array(data[self.mapping['condition']], dtype=np.int)
        if np.any(self.condition<0) or np.any(self.condition>=self.design.nconditions()):
            raise ValueError
        
    def from_pandas_dataframe( self, data, format="wide"):
        self.org_data=data
        self.ntrials=data.shape[0]

        # RT
        self.RT=np.array( data[self.mapping['RT']], dtype=np.float)
        self.RT[np.logical_not(np.isfinite(self.RT))]=np.nan
        
        # SSD
        self.SSD=np.array( data[self.mapping['SSD']], dtype=np.float)
        self.SSD[np.logical_not(np.isfinite(self.SSD))]=np.nan
        
        # responses
        resp=pd.Categorical.from_array(data[self.mapping['response']])
        if not list_is_eq(resp.levels, self.design.responses):
            raise ValueError('response-array in data.frame does not match design: %s'%str(resp.levels))

        self.response=np.zeros( self.ntrials, dtype=np.int )-1
        for rix,r in enumerate(self.design.responses):
            self.response[resp==r]=rix

        if format=="wide":
            self.from_pandas_dataframe_wide(data)
        elif format=='long':
            self.from_pandas_dataframe_long(data)
        else:
            raise ValueError("unknown format '%s'"%format)
        
    def from_pandas_dataframe_long(self, data):
        # conditions
        self.condition=np.array(data['condition'], dtype=np.int)
        if np.any(self.condition<0) or np.any(self.condition>=self.design.nconditions()):
            raise ValueError
        
    def from_pandas_dataframe_wide(self, data):
        # conditions
        self.condition=np.zeros( self.ntrials, dtype=np.int)-1
        
        for i in range(self.ntrials):
            row=data.irow(i)
            cidx=[row[fac] for fac in self.design.factors]
            self.condition[i]=self.design.condidx(cidx)
            
    def as_dataframe(self, form='long', conditions_expanded=True):
        if form=='long':
            if conditions_expanded:
                df=pd.DataFrame({'condition':[":".join(self.design.condidx(c)) for c in self.condition],
                                 'SSD':self.SSD,
                                 'RT':self.RT,
                                 'response':self.response})
            else:
                df=pd.DataFrame({'condition':self.condition,
                                 'SSD':self.SSD,
                                 'RT':self.RT,
                                 'response':self.response})
        elif form=='wide':
            df=pd.DataFrame({'SSD':self.SSD,
                             'RT':self.RT,
                             'response':self.response})
            for cidx,col in enumerate(self.design.factors):
                df[col]=""
                for i,cond in enumerate(self.condition):
                    df[col].iloc[i]=self.design.condidx(cond)[cidx]
        else:
            raise ValueError("don't know how to handle format %s"%form)
        return df

if __name__=="__main__":
    factors=[{'sleepdep':['normal','deprived']},
             {'stimulus':['left', 'right']}]
    responses=['left','right']
    design=Design(factors,responses, 'stimulus')
    dat=pd.read_csv('./data/sleep_stop_onesubj_test.csv')
    print dat.shape[0]
    dat.columns=['sleepdep','stimulus','SSD','response','correct', 'RT']
    ds=StopTaskDataSet(design,dat)
    print ds.as_dataframe(conditions_expanded=False).head(50)
    print ds.as_dataframe(conditions_expanded=True).head(50)
    print ds.as_dataframe(form='wide', conditions_expanded=True).head(50)
