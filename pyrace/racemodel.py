import scipy
import pylab as pl
from collections import Iterable
from itertools import cycle

from design import *
from tools import *
from data import *
import crace


class Accumulator:
    def __init__(self):
        self.parnames=[]
        self.name=""
    def __repr__(self):
        return self.__class__.__name__+"(%s; %s)"%( self.name, ",".join([ "=".join([k,str(self.__dict__[k])]) for k in self.parnames]))
    def pdf(self,t):
        raise NotImplementedError
    def cdf(self,t):
        raise NotImplementedError
  
class RaceModel(object):
    def __init__(self):
        pass


class Parameters(object):
    parnames=[]
    lower   =[]
    upper   =[]

    def __init__(self, *args, **kwargs):
        self.pars={}        
        if len(args)>len(self):
            raise ValueError('too many input args')
        for i,arg in enumerate(args):
            self.pars[self.__class__.parnames[i]]=arg

        for k,v in kwargs.items():
            self.__setattr__(k,v)

    def __repr__(self):
        r="{cname}({parlist})".format(cname=self.__class__.__name__,
                                      parlist=",".join(["%s=%s"%(k,self.pars[k]) for k in self.__class__.parnames]))
        return r
            
    def __getattr__(self, k):
        if k=='pars':
            return self.pars
        elif k in self.__class__.parnames:
            return self.pars[k]
        else:
            raise ValueError("don't have parameter '%s'"%k)

    def __setstate__(self,state):
        self.__dict__=state
        
    def __setattr__(self, k, v):
        if k=='__dict__':
            super(Parameters,self).__setattr__(k,v)
        elif k=='pars':
            self.__dict__[k]=v
        elif k in self.__class__.parnames:
            self.pars[k]=v
        else:
#            self.__dict__[k]=v
            raise ValueError("don't have parameter '%s'"%k)

    def __iter__(self):
        for k in self.__class__.parnames:
            yield self.pars[k]
    
    def __getitem__(self, i):
        if i<len(self.__class__.parnames) and i>=0:
            return self.pars[self.__class__.parnames[i]]
        else:
            raise ValueError("index out of bound: %i"%i)

    def __setitem__(self, i, v):
        if i<len(self.__class__.parnames) and i>=0:
            self.pars[self.__class__.parnames[i]]=v
        else:
            raise ValueError("index out of bound: %i"%i)

    def __len__(self):
        return len(self.__class__.parnames)


    def random(self):
        """
        Set the values to random elements within the borders
        set by lower and upper.
        """
        for i,(l,u) in enumerate(zip( self.__class__.lower, self.__class__.upper)):
            self.__setitem__(i, np.random.uniform( l, u ))
        
    
class StopTaskRaceModel(RaceModel):
    """
    This implements a generic Stop-Signal Task Race Model with the possibility
    to add a mixture of trials in which the GO or the STOP-signal has been missed.
    """
    
    def __init__(self, design, go_accumulators, stop_accumulators, 
                 prob_go_fail=1e-6, prob_trigger_fail=1e-6):
        """
        go_accumulators : list
          design.nconditions()-long list; each element is a
          design.nresponses()-long list containing the actual accumulators
          for each response
          
        stop_accumulators : list
          design.nconditions()-long list, each is an accumulator
          
        prob_go_fail : float [0,1]
          probability not to observe the GO-stimulus 
          
        prob_trigger_fail : float [0,1]
          probability not to observe the STOP-stimulus
        """
        self.design=design
        self.set_accumulators(go_accumulators, stop_accumulators)
        self.set_mixing_probabilities(prob_go_fail, prob_trigger_fail)
        
    def set_mixing_probabilities(self, prob_go_fail, prob_trigger_fail):
        if isinstance(prob_go_fail,Iterable):
            if len(prob_go_fail)!=self.design.nconditions():
                raise ValueError('prob_go_fail needs to have %i entries, have %i'%(self.design.nconditions(),
                                                                                   len(prob_go_fail)))
            else:
                self.prob_go_fail=prob_go_fail
        else:
            self.prob_go_fail=[prob_go_fail for i in range(self.design.nconditions())]

        if isinstance(prob_trigger_fail,Iterable):
            if len(prob_trigger_fail)!=self.design.nconditions():
                raise ValueError('prob_trigger_fail needs to have %i entries, have %i'%(self.design.nconditions(),
                                                                                       len(prob_trigger_fail)))
            else:
                self.prob_trigger_fail=prob_trigger_fail
        else:
            self.prob_trigger_fail=[prob_trigger_fail for i in range(self.design.nconditions())]
    
    def set_accumulators(self, go_accumulators, stop_accumulators):
        if len(go_accumulators)!=self.design.nconditions():
            raise ValueError('need GO-accumulators for each of the %i conditions'%(self.design.nconditions()))
        else:
            for goacc in go_accumulators:
                if len(goacc)!=self.design.nresponses():
                    raise ValueError('need %i GO-accumulators'%(self.design.nresponses()))
        if len(stop_accumulators)!=self.design.nconditions():
            raise ValueError('need %i STOP-accumulators'%(self.design.nconditions()))
        self.go_accumulators=go_accumulators
        self.stop_accumulators=stop_accumulators
    
    def dens_acc_go(self, t, condition, nacc):
        """
        Distribution of the nacc'th accumulator given
        that it is the first to reach threshold using
        only the GO-accumulators.
        
        In condition number 'condition'.
        
        equivalent to dfun.1 and dfun.2
        """
        accumulators=self.go_accumulators[condition]
        if nacc<0 or nacc>=len(accumulators):
            raise ValueError("nacc must be between 0 and %i"%(len(accumulators)-1))
        if np.any(t<1e-5):
            raise ValueError('need positive t here')
        out = np.ones_like(t, dtype=np.float)
        for i in range(len(accumulators)):
            if i==nacc:
                out *= accumulators[i].pdf(t)
            else:
                out *= (1-accumulators[i].cdf(t))
        if out.ndim>0:
            out[np.logical_not(np.isfinite(out)) | (out<0)]=0.0
        return out
    
    def dens_acc_stop(self, t, condition, SSD, nacc):
        """
        # Respond 2nd accumulator (1st response) in stop architecture

        equivalent to dfun.1s and dfun.2s

        except, when nacc<0, it returns p_s(1-F_1)(1-F_2)(...)
        """
        goacc=self.go_accumulators[condition]
        if nacc>=len(goacc):
            raise ValueError("nacc must be between 0 and %i"%(len(goacc)-1))
        
        stacc=self.stop_accumulators[condition]
        tstop=t-SSD
        out=(1-stacc.cdf(tstop)) if nacc>=0 else stacc.pdf(tstop)
        
        for i,acc in enumerate(goacc):
            if i==nacc:
                out*=acc.pdf(t)
            else:
                out*=(1-acc.cdf(t))
        if out.ndim>0:
            out[np.logical_not(np.isfinite(out)) | (out<0)]=0.0
        return out

    def dfun_stop(self, condition, SSD):
        """
        Successful stop probability when at SSD.
        
        equivalent to dfun.stop
        """
        go_accs=self.go_accumulators[condition]
        st_acc=self.stop_accumulators[condition]
        def tmpf(t,gos,st,SSD):
            r=st.pdf( t-SSD )
            for acc in gos:
                r*=(1-acc.cdf( t ))
            return r
        pstop=scipy.integrate.quad(tmpf, st_acc.ter+SSD, np.infty, args=(go_accs,st_acc,SSD))[0]
        return np.maximum( np.minimum( pstop,1), 0) # protect from numerical errors
        
    def likelihood_trials(self, dat):
        """generic and slow implementation (overwrite in child for performance)"""
        if not isinstance(dat,StopTaskDataSet):
            raise ValueError('data needs to be a StopTaskDataSet')
            
        L=np.zeros(dat.ntrials,dtype=np.float)*np.nan
        for cond in range(self.design.nconditions()):
            print cond, self.design.condidx(cond)
            pgf=self.prob_go_fail[cond]
            ptf=self.prob_trigger_fail[cond]
            
            # Failed go: GO(NA) = pgf
            idx=(dat.condition==cond) & (np.isnan(dat.SSD)) & (dat.response<0)
            L[idx]=pgf
            
            # Sucessful stops
            #if np.any(np.isfinite(dat.SSD))
            cssds=dat.SSD[(np.isfinite(dat.SSD)) & (dat.condition==cond) & (dat.response<0)]
            if len(cssds)>0:
                nstop=stats.itemfreq(cssds) #dat.SSD[np.isfinite(dat.SSD)])
                for j in nstop[:,0]:
                    idx=(dat.condition==cond) & (dat.response<0) & (dat.SSD==j)
                    # STOP(NA): pgf + (1-pgf)*(1-ptf)*pstop
                    L[idx]= (pgf + (1-pgf)*(1-ptf)*self.dfun_stop(cond, j))        
                       
            # trials with responses
            for resp in range(self.design.nresponses()):
                
                # Go-Responses
                idx=(dat.condition==cond) & (dat.response==resp) & (np.isnan(dat.SSD))
                if len(idx)>0:
                    L[idx]=(1-pgf)*self.dens_acc_go( dat.RT[idx], cond, resp )
                
                # STOP trials with a response
                idx=(dat.condition==cond) & (dat.response==resp) & (np.isfinite(dat.SSD))
                rts=dat.RT[idx]
                ssds=dat.SSD[idx]
                if len(rts)>0:
                    L[idx]=(1-pgf)*(ptf+(1-ptf)*self.dens_acc_stop(rts, cond, ssds, resp))
        return L

    def sample(self, n, condition, SSD=None, upper_limit=np.infty):
        """
        sample from condition accumulators.
        If SSD!=None, use the stop-distribution, else only go-accs.

        returns (response, RT) both are n-arrays
        """
        # sample RTs from all accumulators
        goaccs=self.go_accumulators[condition]
        go_rts=np.array([acc.sample(n, upper_limit=upper_limit) for acc in goaccs])

        rts=np.min(go_rts,axis=0)
        winner=np.argmin(go_rts,axis=0)

        # go-failures
        gof=stats.binom.rvs(1,self.prob_go_fail[condition],size=n).astype(np.bool) # 0/1 for go-failures
        winner[gof]=-1
        rts[gof]=np.nan

        # STOP trials
        if SSD!=None:
            stop_rt=SSD+self.stop_accumulators[condition].sample(n)
            # trigger-failures
            tf=stats.binom.rvs(1,self.prob_trigger_fail[condition],size=n).astype(np.bool) # 0/1 for trigger-failures
            winner[(stop_rt<rts) & np.logical_not(tf)]=-1            
            rts[(stop_rt<rts) & np.logical_not(tf)]=np.nan #stop_rt[stop_rt<rts]

        return winner,rts

    def simulate_ssd_dist(self, ngo, SSD, upper_limit=np.infty, name=None):
        """
        SSD is a (nssd x 2) array where
            SSD[:,0] are the SSDs and
            SSD[:,1] are the number of trials desired
        OR:

        SSD is a nconditions-long list of such (nssd x 2) arrays
            
        ngo is the number of go-trials, so that ntrials=ngo+np.sum(SSD[:,1])
        """
        if isinstance(SSD, list):
            if len(SSD)!=self.design.nconditions():
                raise ValueError('need distribution for each condition')
            for ssdc in SSD:
                if ssdc.shape[1]!=2:
                    raise ValueError(str(ssdc.shape))
        elif SSD.shape[1]!=2:
            raise ValueError


        conditions=[]
        RT=[]
        response=[]
        SSDs=[]

        for cond in range(self.design.nconditions()):
            if isinstance(SSD,list):
                cSSD=SSD[cond]
            else:
                cSSD=SSD
            nsim=int(ngo+np.sum(cSSD[:,1]))                
            conditions+=[cond]*nsim
            
            # go-trials
            resp,rt=self.sample(ngo, cond, upper_limit=upper_limit)
            RT+=list(rt)
            response+=list(resp)
            SSDs+=[np.nan]*ngo

            # stop-trials
            for ssd,ssdn in zip(cSSD[:,0], cSSD[:,1].astype(np.int)):
                SSDs+=[ssd]*ssdn
                resp,rt=self.sample(ssdn, cond, SSD=ssd)
                RT+=list(rt)
                response+=list(resp)

        ## wrap up in StopTaskDataSet
        dat={'condition':np.array(conditions, dtype=np.int),
             'RT':np.array(RT, dtype=np.double),
             'response':np.array(response, dtype=np.int),
             'SSD':np.array(SSDs, dtype=np.double)}
        ds=StopTaskDataSet(self.design, dat, format='dict', name=name)
        return ds
    
    def simulate(self, nsim, pstop=.25, SSD=np.array([.1, .2, .3, .4, .5]), upper_limit=np.infty, name=None):
        """Simulate a dataset of nsim*nconditions trials corresponding to self.design.

        There will be nsim*nconditions trials, a fraction pstop of which are stop-trials.
        There are an equal number of stop-trials for each of the provided SSDs
        if a single number pstop is used, else there can be a different
        pstop for each SSD.
        """
        nssd=len(SSD)
        if isinstance(pstop, (int, long, float, complex)) or len(pstop)==1:
            pstop=np.array([pstop/float(nssd) for i in range(nssd)], dtype=np.double)
        if len(pstop)!=nssd:
            raise ValueError("pstop and SSD array must be same length, have %i!=%i"%(len(pstop),nssd))

        ssdn = np.round(np.sum(pstop)*nsim/float(nssd))
        nstop = nssd*ssdn
        ngo = nsim-(nstop)
#        print "Simulating %i GO and %i STOP trials (%i of each of %i SSDs) (%i overall)"%(ngo, nstop, ssdn, nssd, ngo+nstop)

        conditions=[]
        RT=[]
        response=[]
        SSDs=[]
        for cond in range(self.design.nconditions()):
#            print "Condition: %s"%self.design.condidx(cond)
            conditions+=[cond]*nsim
            
            # go-trials
            resp,rt=self.sample(ngo, cond, upper_limit=upper_limit)
            RT+=list(rt)
            response+=list(resp)
            SSDs+=[np.nan]*ngo

            # stop-trials
            for ssd in SSD:
#                print " SSD=",ssd
                SSDs+=[ssd]*ssdn
                resp,rt=self.sample(ssdn, cond, SSD=ssd)
                RT+=list(rt)
                response+=list(resp)

        ## wrap up in StopTaskDataSet
        dat={'condition':np.array(conditions, dtype=np.int),
             'RT':np.array(RT, dtype=np.double),
             'response':np.array(response, dtype=np.int),
             'SSD':np.array(SSDs, dtype=np.double)}
        ds=StopTaskDataSet(self.design, dat, format='dict', name=name)
        return ds


    def plot_model(self,lims=(0.1,10),SSD=(None,0.2,0.5), pstop=True, pstop_SSD_lim=(0,1.0)):
        """
        Plot PDF for all accumulators (GO and stop) in each condition.
        """
        go_colors=['blue', 'green', 'yellow', 'magenta', 'cyan']
        stop_color='red'
        linestyles=['-', '--', '-.', ':']
        linecycler = cycle(linestyles)
        colorcycler=cycle(go_colors)

        lw=3
        a=int(np.sqrt(self.design.nconditions()))
        b=np.ceil(self.design.nconditions()/a)

        t=np.linspace(lims[0], lims[1], 1000)
        for cond in range(self.design.nconditions()):
            pl.subplot(a,b,cond)

            for ssdi,ssd in enumerate(SSD):
                lsty=linestyles[ssdi]
                for ri,resp in enumerate(self.design.responses):
                    gocol=go_colors[ri]
                    if ssd==None:
                        pl.plot(t, self.dens_acc_go(t, cond, ri), color=gocol, linewidth=lw, label='%s'%resp, ls=lsty)
                    else:
                        pl.plot(t, self.dens_acc_stop(t,cond,ssd,ri), color=gocol, linewidth=lw, label='%s,SSD=%.2f'%(resp,ssd), ls=lsty)
                if ssd!=None:
                    pl.plot(t, self.dens_acc_stop(t, cond, ssd, -1), color=stop_color, linewidth=lw, label='stop (SSD=%.2f)'%ssd, ls=lsty)
                pl.xlabel('t (s)')
                pl.ylabel('PDF')
            pl.title(":".join(self.design.condidx(cond)))
            if cond==0:
                pl.legend()


    def get_pstop_by_ssd(self, cond, SSD):
        return np.array([self.dfun_stop(cond, cssd) for cssd in SSD], dtype=np.double)
        
    def plot_pstop(self, data=None, SSD_lims=(0,1.0), npoints=10, subplots=False):
        """
        plot probability of successful stopping as a function of SSD.
        """
        lw=3
        bw=0.01 # bar-width
        a=int(np.sqrt(self.design.nconditions()))
        b=np.ceil(self.design.nconditions()/a)

        ssd=np.linspace(SSD_lims[0], SSD_lims[1], npoints)        
        for cond in range(self.design.nconditions()):
            if subplots:
                pl.subplot(a,b,cond)

            if data:
                tmp=stats.itemfreq(data.SSD[np.isfinite(data.SSD)])
                tmp2=[]
                for val in tmp[:,0]:
                    tmp2.append( np.sum(data.response[data.SSD==val]<0))
                pl.bar( tmp[:,0]-bw/2., tmp2/tmp[:,1], width=bw, alpha=.5)
            y=self.get_pstop_by_ssd(cond,ssd)
            pl.plot(ssd, y, linewidth=lw, label=":".join(self.design.condidx(cond)))
            pl.xlim(SSD_lims[0], SSD_lims[1])
            pl.xlabel('SSD')
            pl.ylabel('p(STOP)')
        if not subplots:
            pl.legend()
        
    def plot_fit_go(self, dat, lims=(0.1,10), nbins=20):
        """
        Plot histogram of GO-trials and the model-solution per condition.
        """
        colors=['red', 'blue', 'green', 'yellow', 'magenta', 'cyan']
        a=int(np.sqrt(self.design.nconditions()))
        b=np.ceil(self.design.nconditions()/a)

        t=np.linspace(lims[0], lims[1], 1000)
        for cond in range(self.design.nconditions()):
            pl.subplot(a,b,cond)
            goix=((dat.condition==cond) & np.isfinite(dat.RT) & np.isnan(dat.SSD))
            for ri,resp in enumerate(self.design.responses):
                d=dat.RT[goix & (dat.response==ri)]
                
                if len(d)>0:
                    h,hx=np.histogram(d, density=True, bins=nbins)
                    hx=hx[:-1]#-(hx[1]-hx[0])/2.
                    h=h*(len(d)/float((len(dat.RT[goix]))))
                    pl.bar(hx, h, width=(hx[1]-hx[0]), alpha=.3, color=colors[ri], label='resp=%s'%resp)
#                    pl.hist(d, nbins, normed=True, alpha=.3, color=colors[ri], label='resp=%s'%resp)
                pl.plot(t, self.dens_acc_go(t, cond, ri), color=colors[ri], linewidth=3)
            pl.title(":".join(self.design.condidx(cond)))
            if cond==0:
                pl.legend()            

    
    def loglikelihood(self,dat):
        return np.sum(np.log(np.maximum(self.likelihood_trials(dat),1e-10)))
    
    def deviance(self,dat):
        return -2.*self.loglikelihood(dat)
    
    def __repr__(self):
        return self.__class__.__name__+"(%s)"%self.parstring()
    
    def parstring(self, full=False):
        if not full:
            return str(self.params)
#            return ",".join(["=".join([str(k),str(v)]) for k,v in self.params.items()])
        else:
            r=""
            for cond in range(self.design.nconditions()):
                r+="{cidx} <-> {cstr}: {goaccstr}/{stopaccstr} / {pgf},{ptf}\n".format(
                        cidx=cond,
                        cstr=self.design.condidx(cond),
                        goaccstr=repr(self.go_accumulators[cond]),
                        stopaccstr=repr(self.stop_accumulators[cond]),
                        pgf=self.prob_go_fail[cond],
                        ptf=self.prob_trigger_fail[cond])
            return r
        
    def init_cmodule(self):
        crace.init()
