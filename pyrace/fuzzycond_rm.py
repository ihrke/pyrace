"""
highly experimental
"""
from .racemodel import StopTaskRaceModel
from .fuzzycond_data import FuzzyConditionStopTaskDataSet

class FuzzyConditionStopTaskRaceModel(StopTaskRaceModel):
    """
    overwrites

    sample()
    simulate()
    likelihood_trials()
    """

    def likelihood_trials(self, dat):
        if not isinstance(dat, FuzzyConditionStopTaskDataSet):
            raise ValueError('data needs to be a FuzzyConditionStopTaskDataSet')

        L=np.zeros(dat.ntrials,dtype=np.float)*np.nan

        for i in range(ntrials):
            print 'trial ',i
            


