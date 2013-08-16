from tools import *
from design import Design
import pandas as pd
import numpy as np

class StopTaskDataSet(object):
    """
    Implements a data-set for the stop-task
    Internally, it stores the data using the following structure:
    
    condition SSD RT response
    
    where condition is an integer number that is expanded by 
    the Design() as is response. 
    
    * if SSD==np.nan -> GO-trial
    * if RT==np.nan -> missed/successful STOP
    * if response==-1 -> miss/successful STOP
    * if condition<0 -> something is wrong with the dataset
    """
    def __init__(self, design, data, mapping={'RT':'RT', 'SSD':'SSD', 'response':'response'}):
        """
        design : Design 
           description of the experiment
        data : pandas.DataFrame
           with columns corresponding to the factors in design
        mapping : dict
           must contain keys "RT", "SSD", 'response' and use the values to index into the
           pandas dataframe.
        """
        self.ntrials=data.shape[0]
        self.org_data=data
        self.mapping=mapping
        self.design=design
        
        # RT
        self.RT=np.array( data[mapping['RT']], dtype=np.float)
        self.RT[np.logical_not(np.isfinite(self.RT))]=np.nan
        
        # SSD
        self.SSD=np.array( data[mapping['SSD']], dtype=np.float)
        self.SSD[np.logical_not(np.isfinite(self.SSD))]=np.nan
        
        # responses
        resp=pd.Categorical.from_array(data[mapping['response']])
        if not list_is_eq(resp.levels, design.responses):
            raise ValueError('response-array in data.frame does not match design: %s'%str(resp.levels))

        self.response=np.zeros( self.ntrials, dtype=np.int )-1
        for rix,r in enumerate(self.design.responses):
            self.response[resp==r]=rix
        
        # conditions
        self.condition=np.zeros( self.ntrials, dtype=np.int)-1
        
        for i in range(self.ntrials):
            row=data.irow(i)
            cidx=[row[fac] for fac in self.design.factors]
            self.condition[i]=self.design.condidx(cidx)
    def as_dataframe(self):
        return pd.DataFrame({'condition':self.condition,
                             'SSD':self.SSD,
                             'RT':self.RT,
                             'response':self.response})

if __name__=="__main__":
    factors=[{'sleepdep':['normal','deprived']},
             {'stimulus':['left', 'right']}]
    responses=['left','right']
    design=Design(factors,responses, 'stimulus')
    dat=pd.read_csv('../data/sleep_stop_onesubj_test.csv')
    print dat.shape[0]
    dat.columns=['sleepdep','stimulus','SSD','response','correct', 'RT']
    ds=StopTaskDataSet(design,dat)
    print ds.as_dataframe().head(50)
