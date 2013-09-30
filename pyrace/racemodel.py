import scipy
import pylab as pl
from collections import Iterable
from textwrap import wrap
from itertools import cycle

from .design import *
from .tools import *
from .data import *
from .param import *

class Accumulator:
    def __init__(self):
        self.parnames=[]
        self.name=""
    def __repr__(self):
        return self.__class__.__name__+"(%s; %s)"%( self.name, ",".join([ "=".join([k,"%.2f"%(self.__dict__[k])])
                                                                          for k in self.parnames]))
    def pdf(self,t):
        raise NotImplementedError
    def cdf(self,t):
        raise NotImplementedError
  
class RaceModel(object):
    paramspec=Parameters;
    def __init__(self):
        pass

    def name(self):
        return (self.__class__.__name__).split(".")[-1]

    
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
    def copy(self):
        """for multiprocessing"""
        m=self.__class__(self.design, self.go_accumulators, self.stop_accumulators, self.prob_go_fail,self.prob_trigger_fail)
        return m

    def trans(self, pars):
        """generic logistic transformation"""
        x=np.array([trans_logistic(pars[i], a=pars.bound_lower(i), b=pars.bound_upper(i)) for i in range(len(pars))])
        return x
    
    def untrans(self, x):
        """generic logistic transformation"""        
        pars=self.paramspec([trans_logistic(x[i], a=self.paramspec.lower[i], b=self.paramspec.upper[i], inverse=True)
                             for i in range(len(x))])
        return pars

        
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
        if np.any(t<1e-15):
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
#        pstop=scipy.integrate.quad(tmpf, st_acc.ter+SSD, np.infty, args=(go_accs,st_acc,SSD))[0]
        pstop=scipy.integrate.quad(tmpf, SSD, np.infty, args=(go_accs,st_acc,SSD))[0]
        return np.maximum( np.minimum( pstop,1), 0) # protect from numerical errors
        
    def likelihood_trials(self, dat):
        """generic and slow implementation (overwrite in child for performance)

        # pgf  : probability of a go failure
        # ptf  : probability of a trigger failure
        # pstop: probability of a sucessful stop given no go or trigger fail 
        # Lg(t): likelihood of response at time t on go trial given no go fail
        # Ls(t): likelihood of response at time t on stop trial given no go or trigger fail
        # 
        # GO(NA)  : pgf
        # STOP(NA): pgf + (1-pgf)*(1-ptf)*pstop
        # GO(t)   : (1-pgf)*Lg(t)
        # STOP(t) : (1-pgf)*[ptf*Lg(t) + (1-ptf)*Ls(t)]
        """
        if not isinstance(dat,StopTaskDataSet):
            raise ValueError('data needs to be a StopTaskDataSet')
            
        L=np.zeros(dat.ntrials,dtype=np.float)*np.nan
        for cond in range(self.design.nconditions()):
            #print cond, self.design.condidx(cond)
            pgf=self.prob_go_fail[cond]
            ptf=self.prob_trigger_fail[cond]
            
            # Failed go: GO(NA) = pgfprint
            idx=(dat.condition==cond) & (np.isnan(dat.SSD)) & (dat.response<0)
            L[idx]=pgf
            
            # Sucessful stops: STOP(NA): pgf + (1-pgf)*(1-ptf)*pstop
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
                
                # Go-Responses: GO(t)   : (1-pgf)*Lg(t)
                idx=(dat.condition==cond) & (dat.response==resp) & (np.isnan(dat.SSD))
                if len(idx)>0:
                    L[idx]=(1-pgf)*self.dens_acc_go( dat.RT[idx], cond, resp )
                
                # STOP trials with a response: STOP(t) : (1-pgf)*[ptf*Lg(t) + (1-ptf)*Ls(t)]
                idx=(dat.condition==cond) & (dat.response==resp) & (np.isfinite(dat.SSD))
                rts=dat.RT[idx]
                ssds=dat.SSD[idx]
                if len(rts)>0:
                    L[idx]=(1-pgf)*(ptf*self.dens_acc_go(rts, cond, resp)+(1-ptf)*self.dens_acc_stop(rts, cond, ssds, resp))
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
        ngo is the number of go trials OR a list of ngo, one for each condition
        
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

        if isinstance(ngo, list):
            if len(ngo)!=self.design.nconditions():
                raise ValueError('need ngo-trials for each condition, ngo=%s'%(str(ngo)))

        conditions=[]
        RT=[]
        response=[]
        SSDs=[]

        for cond in range(self.design.nconditions()):
            if isinstance(SSD,list):
                cSSD=SSD[cond]
            else:
                cSSD=SSD
            if isinstance(ngo,list):
                cngo=ngo[cond]
            else:
                cngo=ngo
            nsim=int(cngo+np.sum(cSSD[:,1]))                
            conditions+=[cond]*nsim
            
            # go-trials
            resp,rt=self.sample(cngo, cond, upper_limit=upper_limit)
            RT+=list(rt)
            response+=list(resp)
            SSDs+=[np.nan]*cngo

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
        b=np.ceil(self.design.nconditions()/float(a))

        t=np.linspace(lims[0], lims[1], 1000)
        for cond in range(self.design.nconditions()):
            pl.subplot(a,b,cond+1)

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
                pl.xlim(lims)
            pl.title(":".join(self.design.condidx(cond)))
            if cond==0:
                pl.legend()


    def get_pstop_by_ssd(self, cond, SSD):
        return np.array([self.dfun_stop(cond, cssd) for cssd in SSD], dtype=np.double)
        
    def plot_pstop(self, data=None, SSD_lims=(0,1.0), npoints=8, subplots=True, conditions='all'):
        """
        plot probability of successful stopping as a function of SSD.
        """
        lw=3
        bw=0.01 # bar-width
        if conditions=='all':
            conditions=range(self.design.nconditions())
        elif conditions==None:
            conditions=[range(self.design.nconditions())]
        
        a=int(np.sqrt(len(conditions)))
        b=np.ceil(len(conditions)/float(a))

        ssd=np.linspace(SSD_lims[0], SSD_lims[1], npoints)

        for cix,cond in enumerate(conditions):
            if subplots:
                pl.subplot(a,b,cix+1)
            if isinstance(cond, int):
                condidx=(data.condition==cond)
            if isinstance(cond, Iterable):
                condidx=np.array( [dcond in cond for dcond in data.condition], dtype=np.bool)
    
            if data!=None:
                tmp=data.get_ssd_dist(condition=condidx)
                ssds=tmp[:,0]
                nssds=tmp[:,1].astype(np.int)
                pstop=[np.sum(np.isnan(data.RT[condidx & (data.SSD==cssd)]))/float(nssds[i]) for i,cssd in enumerate(ssds)]                    
                pl.bar( ssds-bw/2., pstop, width=bw, alpha=.5)
            if isinstance(cond, Iterable):
                for con in cond:
                    y=self.get_pstop_by_ssd(con,ssd)
                    pl.plot(ssd, y, linewidth=lw, label=":".join(self.design.condidx(con)))
            else:
                y=self.get_pstop_by_ssd(cond,ssd)
                pl.plot(ssd, y, linewidth=lw, label=":".join(self.design.condidx(cond)))

            pl.xlim(SSD_lims[0], SSD_lims[1])
            pl.xlabel('SSD')
            pl.ylabel('p(STOP)')
            if subplots:
                if isinstance(cond, Iterable):
                    titlab=(",".join( ":".join(self.design.condidx(con)) for con in cond))
                else:
                    titlab=(":".join(self.design.condidx(cond)))
                # long title, guaranteed to fit w matplotlib
                titlab=("\n".join(wrap(titlab,width=70)))
                pl.title(titlab, fontsize='x-small')
                
        if not subplots:
            pl.legend()
        
    def plot_fit_go(self, dat, lims=(0.001,5), res=100, nbins=20, conditions='all', split_by_response='correct', ylim_fixed=True):
        """
        Plot histogram of GO-trials and the model-solution per condition.

        nbins : int or list
           bins for the data-histogram
        
        res : int
           number of datapoints with which the density plots are resolved
        
        conditions : None, 'all' or list of items or lists
            this is used to collapse over design factors that are not modeled
            e.g., conditions=[ [0,1], [1,2] ] would plot the data of
            [0,1] together and [1,2] together. The density is plotted for
            both conditions on top of each other.

        split_by_response : None, 'correct' or 'response'
           plot separate histograms/densities depending on the response
             'correct'  : split into correct/incorrect responses
             'response' : split into the possible responses

        ylim_fixed : False, True
           should the y-axis be fixed across conditions?
        """
        colors=['red', 'blue', 'green', 'yellow', 'magenta', 'cyan']
        if isinstance(nbins, int):
            bins=np.linspace(lims[0], lims[1], nbins)
        else:
            bins=nbins
        
        if conditions=='all':
            conditions=range(self.design.nconditions())
        elif conditions==None:
            conditions=[range(self.design.nconditions())]

        a=int(np.sqrt(len(conditions)))
        b=np.ceil(len(conditions)/float(a))
        
        t=np.linspace(lims[0], lims[1], res)
        maxy=[]
        for cix,cond in enumerate(conditions):
            pl.subplot(a,b,cix+1)
            if isinstance(cond, int):
                condidx=(dat.condition==cond)
            if isinstance(cond, Iterable):
                condidx=np.array( [dcond in cond for dcond in dat.condition], dtype=np.bool)
                
            goix=((condidx) & np.isfinite(dat.RT) & np.isnan(dat.SSD))

            ## construct response indices
            if split_by_response not in ['response', 'correct']:
                split_by_response=='response'
                
            if split_by_response=='response':
                respidcs=[ (ri, resp, (dat.response==ri)) for ri,resp in enumerate(self.design.responses)]
            else:
                con=cond[0] if isinstance(cond,Iterable) else cond
                corr_ri=self.design.correct_response(con, as_index=True)
                inc_ri=((corr_ri + 1) % self.design.nresponses())
                respidcs=[(ri, resp, dat.correct==rbo) for rbo,resp,ri in zip([True, False], ['correct', 'incorrect'], [corr_ri,inc_ri])]

            for ri,resp,respix in respidcs:
                #resp, respix=resptup
                d=dat.RT[goix & respix]
                
                if len(d)>0:
                    h,hx=np.histogram(d, density=True, bins=bins)
                    hx=hx[:-1]
                    h=h*(len(d)/float((len(dat.RT[goix]))))
                    pl.bar(hx, h, width=(hx[1]-hx[0]), alpha=.3, color=colors[ri], label='resp=%s'%(resp))

                if isinstance(cond, Iterable):
                    for con in cond:
                        if split_by_response=='correct':
                            corr_ri=self.design.correct_response(con, as_index=True)
                            inc_ri =( (corr_ri+1) % self.design.nresponses())
                            cri=corr_ri if resp=='correct' else inc_ri
                        else:
                            cri=ri
                        pl.plot(t, self.dens_acc_go(t, con, cri), color=colors[ri], linewidth=3)
                else:
                    pl.plot(t, self.dens_acc_go(t, cond, ri), color=colors[ri], linewidth=3)
            if isinstance(cond, Iterable):
                titlab=(",".join( ":".join(self.design.condidx(con)) for con in cond))
            else:
                titlab=(":".join(self.design.condidx(cond)))
            # long title, guaranteed to fit w matplotlib
            titlab=("\n".join(wrap(titlab,width=70)))
            pl.title(titlab, fontsize='x-small')

            if cix==0:
                pl.legend()
            pl.xlim(lims)
            maxy.append(pl.ylim()[1])
        if ylim_fixed:
            maxy=np.max(maxy)
            for cix,cond in enumerate(conditions):
                pl.subplot(a,b,cix+1)
                pl.ylim(0,maxy)

    
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

    def npars(self):
        return len(self.params)
    
    def AIC(self, dat):
        return 2*self.npars()-2*self.loglikelihood(dat)
        
    def BIC(self, dat):
        return self.npars()*np.log(dat.ntrials)-2*self.loglikelihood(dat)

