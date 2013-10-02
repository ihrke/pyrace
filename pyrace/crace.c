#include <stdio.h>
#include <math.h>
#include "crace.h"

/** WARNING:

	 Gotta be really careful to include the correct GSL headers.
	 If you don't include them properly for any function you use, 
	 you will not get any warning but instead meaningless results... weird...

 */
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

/** \brief maximum of two elements. */
#define MAX(a,b) ((a) > (b) ? (a):(b))
/** \brief minimum of two elements. */
#define MIN(a,b) ((a) < (b) ? (a):(b))
/** \brief square of a number. */
#define SQR(a) ((a)*(a))
/** \brief absolute value of a number */
#define ABS(a) ( ((a)<0) ? (-1*(a)) : (a) )
#define ISODD(x)        ((((x)%2)==0)? (0) : (1))


#define dprintf(...) do{																\
		fprintf(stderr, "%s (%i), %s(): ", \
				  __FILE__, __LINE__, __FUNCTION__);								\
		fprintf(stderr, ## __VA_ARGS__);												\
	 } while(0)


double pnormP(double x, double mean, double sd){
  /* PDF and CDF 
  return np.where(np.abs(x-mean)<7.*sd, stats.norm.cdf(x, loc=mean,scale=sd), np.where(x<mean,0,1))
  TESTED
  */
  return (double)((ABS(x-mean)<(7*sd)) ? (gsl_cdf_gaussian_P(x-mean,sd)) : ( (x<mean) ? (0) : 1));
}

double dnormP(double x, double mean, double sd){
  /*     return np.where(np.abs(x-mean)<7.*sd,stats.norm.pdf(x, loc=mean, scale=sd),0) 
	TESTED
	*/
  return (double)(ABS(x-mean)<(7*sd)) ? (gsl_ran_gaussian_pdf(x-mean,sd)) : (0);
}

double slognorm_pdf(double t, double ter, double mu, double sigma){
  /* Python:
        t=np.maximum(t-self.ter, 1e-5) # absorbed into pdf
        r=1./(t*np.sqrt(2*np.pi)*self.sigma)*np.exp(-(np.log(t)-self.mu)**2/(2*self.sigma**2))
        return np.maximum(0.0, r)
  */
  t=MAX(t-ter, 1e-5);
  double r=1./(t*sqrt(2*M_PI)*sigma)*exp(- SQR(log(t)-mu)/(2*SQR(sigma)) );
  return (double)MAX(0,r);

}

double slognorm_cdf(double t, double ter, double mu, double sigma){
  /* Python:

        t=np.maximum(t-self.ter, 1e-5) # absorbed into cdf
        r=.5+.5*erf((np.log(t)-self.mu)/(np.sqrt(2)*self.sigma))
        return np.minimum( np.maximum( 0., r ), 1.)
   */
   t=MAX(t-ter, 1e-5);
   double r=.5+.5*gsl_sf_erf( (log(t)-mu)/(sqrt(2)*sigma));
   return (double)MIN( MAX( 0, r ), 1);
}


double lba_pdf(double t, double ter, double A, double v, double sv, double b){
  /* TESTED */
  t=MAX(t-ter, 1e-5);
  if( A<1e-10 ){ /* LATER */
	 return MAX( 0, (b/(SQR(t)))*dnormP(b/t, v, sv)/pnormP(v/sv,0,1));
  }
  double zs=t*sv;
  double zu=t*v;
  double bminuszu=b-zu;
  double bzu=bminuszu/zs;
  double bzumax=(bminuszu-A)/zs;
  return (double)MAX( 0, ( (v*(pnormP(bzu,0,1)-pnormP(bzumax,0,1))+(sv * (dnormP(bzumax,0,1)-dnormP(bzu,0,1))))/A)/pnormP(v/sv,0,1));
}
double lba_cdf(double t, double ter, double A, double v, double sv, double b){
  /* TESTED */
  t=MAX(t-ter, 1e-5);
  if( A<1e-10 ){
	 return MIN(1, MAX( 0, (pnormP(b/t,v,sv)/pnormP(v/sv,0,1)))); /* LATER */
  }
  double zs=t*sv;
  double zu=t*v;
  double bminuszu=b-zu;
  double xx=bminuszu-A;
  double bzu=bminuszu/zs;
  double bzumax=xx/zs;
  double tmp1=zs*(dnormP(bzumax,0,1)-dnormP(bzu,0,1));
  double tmp2=xx*pnormP(bzumax,0,1)-bminuszu*pnormP(bzu,0,1);
  
  return (double)MIN( 1, MAX( 0, (1+(tmp1+tmp2)/A)/pnormP(v/sv,0,1)));
}

double varwald_pdf(double t, double alpha, double gamma, double theta, double A){
  /* Python:

        t=np.maximum(t-self.theta, 1e-5) # absorbed into pdf
        sqrt_t=np.sqrt(t)

        # reparametrization
        a=self.A/2.0
        k=self.alpha-self.A/2.0
        l=self.gamma

        if self.A<1e-10: # this is the solution without starting-point variability
            r=self.alpha/(np.sqrt(2*np.pi*(t**3)))*np.exp(- ((self.alpha-self.gamma*t)**2)/(2*t))
        elif self.gamma<1e-10:
            r=np.exp( -.5*( np.log(2)+np.log(np.pi)+np.log(t))
                      + np.log( np.exp(-( (k-a)**2/(2*t)))-np.exp(-( (k+a)**2/(2*t) )) )
                      - np.log(2) - np.log(a) )
        else:
            r=np.exp( np.log( (np.exp(- (a-k+t*l)**2/(2*t) )-np.exp(- (a+k-t*l)**2/(2*t) ))/np.sqrt(2*np.pi*t)
                              + np.exp(np.log(.5)+np.log(l))*( 2*pnormP( (-k+a)/sqrt_t + sqrt_t*l)-1
                                                             + 2*pnormP(  (k+a)/sqrt_t - sqrt_t*l)-1) )
                      - np.log(2) - np.log(a))

        return np.maximum(0.0, np.where( np.isnan(r), 0, r))

  */
    double r;
    t=MAX(t-theta,1e-5);
    double sqrt_t=sqrt(t);
    double a=A/2.0;
    double k=alpha-a;
    double l=gamma;

    if( A<1e-10 ){ /* normal Wald */
        return (double)MAX(0, (alpha/(sqrt(2*M_PI*(t*t*t))))*exp( -( SQR(alpha-gamma*t))/(2*t)));
    }

    if( gamma<1e-10 ){ /* special solution for v=0 */
        r=exp( -.5*( log(2.0)+log(M_PI)+log(t))
                + log( exp( -( SQR(k-a)/(2*t))) - exp(-( SQR(k+a)/(2*t) )) )
                - log(2.0) - log(a));
    } else {
        r=exp( log( ( exp(-SQR(a-k+t*l)/(2*t) ) - exp(-SQR(a+k-t*l)/(2*t) ))/sqrt(2*M_PI*t)
            + exp( log(.5)+log(l))*( 2*pnormP( (-k+a)/sqrt_t + sqrt_t*l,0,1)-1
                                   + 2*pnormP(  (k+a)/sqrt_t - sqrt_t*l,0,1)-1))
            - log(2.0) -log(a));
    }

    if(!isfinite(r)) r=0.0;
    return (double)MAX(0, r);
}

double varwald_cdf(double t, double alpha, double gamma, double theta, double A){
  /* Python:

        t=np.maximum(t-self.theta, 1e-5) # absorbed into cdf

        sqrt_t=np.sqrt(t)

        # reparametrization
        a=self.A/2.0
        k=self.alpha-self.A/2.0
        l=self.gamma

        if self.A<1e-10: # this is the solution without starting-point variability
            r=pnormP( (self.gamma*t-self.alpha)/sqrt_t)+np.exp(2*self.alpha*self.gamma)*pnormP(-(self.gamma*t+self.alpha)/(sqrt_t))
        elif self.gamma<1e-10:
            r=(( -(k+a)*(2*pnormP( (k+a)/sqrt_t)-1)
                 -(k-a)*(2*pnormP(-(k-a)/sqrt_t)-1))/(2*a)) \
              + (1 + np.exp(-.5*(k-a)**2/t - .5*np.log(2) - .5*np.log(np.pi) + .5*np.log(t) - np.log(a))
                 -   np.exp(-.5*(k+a)**2/t - .5*np.log(2) - .5*np.log(np.pi) + .5*np.log(t) - np.log(a)))
        else:
            t1=np.exp( .5*np.log(t)-.5*np.log(2*np.pi) ) * (  np.exp( -((k-a-t*l)**2/t)/2.)
                                                            - np.exp( -((k+a-t*l)**2/t)/2.) ) # ok
            t2=a+(   np.exp(2*l*(k+a)+np.log(pnormP(-(k+a+t*l)/sqrt_t)))
                   - np.exp(2*l*(k-a)+np.log(pnormP(-(k-a+t*l)/sqrt_t))) )/(2*l) # ok
            t4= (.5*(t*l-a-k+.5/l)) * ( 2*pnormP((k+a)/sqrt_t-sqrt_t*l)-1) \
               +(.5*(k-a-t*l-.5/l)) * ( 2*pnormP((k-a)/sqrt_t-sqrt_t*l)-1)
            r=(t4+t2+t1)/(2*a)

        return np.minimum( np.maximum( 0., np.where( np.isnan(r), 0, r) ), 1.)
    */

    double r;
    t=MAX(t-theta,1e-5);
    double sqrt_t=sqrt(t);
    double a=A/2.0;
    double k=alpha-a;
    double l=gamma;

    if( A<1e-10 ){ /* normal Wald */
        return (double)MIN(1, MAX( 0, pnormP( (gamma*t-alpha)/sqrt(t), 0,1)+exp(2*alpha*gamma)*pnormP( -(gamma*t+alpha)/sqrt(t), 0,1) ));
    }

    if( gamma<1e-10){
        r=(( -(k+a)*(2*pnormP( (k+a)/sqrt_t,0,1)-1)
             -(k-a)*(2*pnormP(-(k-a)/sqrt_t,0,1)-1))/(2*a))
           + ( 1+ exp(-.5*SQR(k-a)/t - .5*log(2.0) - .5*log(M_PI) + .5*log(t) - log(a))
                - exp(-.5*SQR(k+a)/t - .5*log(2.0) - .5*log(M_PI) + .5*log(t) - log(a)));
    } else {
        double t1, t2, t4;
        t1=exp( .5*log(t) - .5*log(2*M_PI) ) * ( exp( -(SQR(k-a-t*l)/t)/2.)
                                                -exp( -(SQR(k+a-t*l)/t)/2.) );
        t2=a+( exp(2*l*(k+a)+log(pnormP(-(k+a+t*l)/sqrt_t,0,1)))
              -exp(2*l*(k-a)+log(pnormP(-(k-a+t*l)/sqrt_t,0,1))) )/(2*l);
        t4= (.5*(t*l-a-k+.5/l)) * (2*pnormP((k+a)/sqrt_t-sqrt_t*l,0,1)-1)
           +(.5*(k-a-t*l-.5/l)) * (2*pnormP((k-a)/sqrt_t-sqrt_t*l,0,1)-1);
        r=(t4+t2+t1)/(2*a);
    }
    if(!isfinite(r)) r=0.0;

    return (double)MIN( MAX( 0, r ), 1);
}

double wald_pdf(double t, double alpha, double gamma, double theta){
  /* Python-side:

    def pdf(self, t):
        t=np.maximum(t-self.theta, 1e-5) # absorbed into pdf
        r=self.alpha/(np.sqrt(2*np.pi*(t**3)))*np.exp(- ((self.alpha-self.gamma*t)**2)/(2*t))
        return np.maximum(0.0, r)
  */
  t=MAX(t-theta,1e-5);
  return (double)MAX(0, (alpha/(sqrt(2*M_PI*(t*t*t))))*exp( -( SQR(alpha-gamma*t))/(2*t)));
}
double wald_cdf(double t, double alpha, double gamma, double theta){
  /* Python-side:
        
    def cdf(self,t):
        t=np.maximum(t-self.theta, 1e-5) # absorbed into cdf
        r=pnormP( (self.gamma*t-self.alpha)/np.sqrt(t))+np.exp(2*self.alpha*self.gamma)*pnormP(-(self.gamma*t+self.alpha)/(np.sqrt(t)))
        return np.minimum( np.maximum( 0., r ), 1.)
  */
  t=MAX(t-theta,1e-5);
  return (double)MIN(1, MAX( 0, pnormP( (gamma*t-alpha)/sqrt(t), 0,1)+exp(2*alpha*gamma)*pnormP( -(gamma*t+alpha)/sqrt(t), 0,1) ));
}


double pstop_integrate( double t, void *params ){
  void *pptr=params;
  int nresponses=*((int*)pptr);   pptr+=sizeof(int);
  int stride=*((int*)pptr);       pptr+=sizeof(int);
  double *gv=*((double**)pptr);   pptr+=sizeof(double*);
  double *gter=*((double**)pptr); pptr+=sizeof(double*);
  double *gA=*((double**)pptr);   pptr+=sizeof(double*);
  double *gb=*((double**)pptr);   pptr+=sizeof(double*);
  double *gsv=*((double**)pptr);  pptr+=sizeof(double*);
  double sv=*((double*)pptr);     pptr+=sizeof(double);
  double ster=*((double*)pptr);   pptr+=sizeof(double);
  double sA=*((double*)pptr);     pptr+=sizeof(double);
  double sb=*((double*)pptr);     pptr+=sizeof(double);
  double ssv=*((double*)pptr);    pptr+=sizeof(double);
  double SSD=*((double*)pptr);

  int i, idx;
  double r=lba_pdf(t-SSD, ster, sA, sv, ssv, sb);

  for(i=0; i<nresponses; i++ ){
	 idx=i*stride;
	 r*=(1-lba_cdf(t, gter[idx], gA[idx], gv[idx], gsv[idx], gb[idx]));
  }
  return r;
}

struct pstop_params {
  int nresponses; int stride;
  double *gv; double *gter; double *gA; double *gb; double *gsv; 
  double sv; double ster; double sA; double sb; double ssv; double SSD;
};

double pstop_integrate2(double t, double *gv, double *gter, double *gA, double *gb, double *gsv, 
							  int nresponses, int stride,
							  double sv, double ster, double sA, double sb, double ssv, double SSD){
  int i, idx;
  double r=lba_pdf(t-SSD, ster, sA, sv, ssv, sb);

  for(i=0; i<nresponses; i++ ){
	 idx=i*stride;
	 r*=(1-lba_cdf(t, gter[idx], gA[idx], gv[idx], gsv[idx], gb[idx]));
  }
  return r;
}

/** Returns raw likelihoods for each trial in pointer L.

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
 */
void sslba_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
								  int *condition, int *response, double *RT, double *SSD, /* data */
								  double *go_v, double *go_ter, double *go_A, double *go_b, double *go_sv,
								  double *stop_v, double *stop_ter, double *stop_A, double *stop_b, double *stop_sv,
								  double *pgf, double *ptf, double *L ){
  
  /*
  dprintf("nconditions=%i, nres=%i, ntri=%i\n", nconditions, nresponses, ntrials);
  for( i=0; i<ntrials; i++ ){
	 dprintf("condition[%i]=%i, response=%i, RT=%f, SSD=%f\n", i, condition[i], response[i], RT[i], SSD[i]);
  }
  dprintf("L=%f\n",*L);
  */

  int i, j, idx;



  /* setup for numerical integration of pstop */
  gsl_integration_workspace * work =gsl_integration_workspace_alloc (1000);
  double result, error;
  int nerror=0;
  int errcode;
  int neval;
  gsl_function F;
  F.function = &pstop_integrate;
  struct pstop_params params;
  F.params=&params;
  params.nresponses=nresponses;
  params.stride=nconditions;
  double pstop;
  gsl_set_error_handler_off();

  
  double dens, densgo, densstop, tmp;

  for( i=0; i<ntrials; i++ ){ /* calc L for each datapoint */
	 /* failed GO */
	 if( isnan(SSD[i]) && response[i]<0 ){
		L[i]=pgf[ condition[i] ];
	 }

	 /* successful stop */
	 else if( response[i]<0 && isfinite(SSD[i]) ){
		//dprintf("succstop, trial=%i\n",i);
		/* set parameters for integration */
		params.gv  =&(go_v  [condition[i]]);
		params.gter=&(go_ter[condition[i]]);
		params.gA  =&(go_A  [condition[i]]);
		params.gb  =&(go_b  [condition[i]]);
		params.gsv =&(go_sv [condition[i]]);
		
		params.sv  =(stop_v  [condition[i]]);
		params.ster=(stop_ter[condition[i]]);
		params.sA  =(stop_A  [condition[i]]);
		params.sb  =(stop_b  [condition[i]]);
		params.ssv =(stop_sv [condition[i]]);
		params.SSD = SSD[i];

		/*		errcode=gsl_integration_qagiu (&F, stop_ter[condition[i]]+SSD[i], 1e-5, 1e-5, 1000, work, &result, &error); */
		errcode=gsl_integration_qagiu (&F, SSD[i], 1e-5, 1e-5, 1000, work, &result, &error); 
		if(errcode){
		  if(nerror==0){
			 dprintf("ERROR during integration, errcode=%i, abserr=%f, result=%f\n", errcode, error, result);
		  }
		  nerror++;
		}
		//dprintf("result=%f, error=%f, neval=%i\n",result, error, (int)(work->size));
		pstop=result;
		L[i]=pgf[condition[i]] + (1-pgf[condition[i]])*(1-ptf[condition[i]])*pstop;
	 }

	 /* go-trials with response */
	 else if( response[i]>=0 && isnan(SSD[i]) ){
		dens=1.0;
		//double lba_pdf(double t, double ter, double A, double v, double sv, double b){
		for( j=0; j<nresponses; j++ ){
		  idx=condition[i]+(j*nconditions);
		  if( j==response[i] ){
			 //			 dprintf("Trial %i, condition=%i, PDF target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 dens*=lba_pdf(RT[i], go_ter[idx], go_A[idx], go_v[idx], go_sv[idx], go_b[idx]);
		  } else {
			 //			 dprintf("Trial %i, condition=%i, CDF non-target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 dens*=(1-lba_cdf(RT[i], go_ter[idx], go_A[idx], go_v[idx], go_sv[idx], go_b[idx]));
		  }
		}
		if( (dens<0) || !isfinite(dens) )
		  dens=0.0;
		L[i]=(1-pgf[condition[i]])*dens;
	 }

	 /* STOP-trials with response 
		 STOP(t) : (1-pgf)*[ptf*Lg(t) + (1-ptf)*Ls(t)]
	  */
	 else if( response[i]>=0 && isfinite(SSD[i]) ){
		idx=condition[i];
		densgo=1.0;
		densstop=(1-lba_cdf(RT[i]-SSD[i], stop_ter[idx], stop_A[idx], stop_v[idx], stop_sv[idx], stop_b[idx]));

		for( j=0; j<nresponses; j++ ){
		  idx=condition[i]+(j*nconditions);
		  if( j==response[i] ){
			 //			 dprintf("Trial %i, condition=%i, PDF target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 tmp=lba_pdf(RT[i], go_ter[idx], go_A[idx], go_v[idx], go_sv[idx], go_b[idx]);
			 densstop*=tmp;
			 densgo*=tmp;
		  } else {
			 //			 dprintf("Trial %i, condition=%i, CDF non-target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 tmp=(1-lba_cdf(RT[i], go_ter[idx], go_A[idx], go_v[idx], go_sv[idx], go_b[idx]));
			 densstop*=tmp;
			 densgo*=tmp;
		  }
		}
		
		L[i]=(1-pgf[condition[i]])*(ptf[condition[i]]*densgo+(1-ptf[condition[i]])*densstop);
	 }

	 /* invalid -> error message */
	 else {
		dprintf("ERROR: trial %i has not been processed by loglikelihood!\n",i);
	 }
  }

  /* clean up */
  gsl_integration_workspace_free(work);
}


/*---------------------------------------------------------------------------*/
/* WALD Model (only Wald accs)
/*---------------------------------------------------------------------------*/

struct pstop_params_wald {
  int nresponses; int stride;
  double *ggamma; double *gtheta; double *galpha;
  double sgamma; double stheta; double salpha; double SSD;
};

double pstop_integrate_wald( double t, void *params ){
  void *pptr=params;
  int nresponses=*((int*)pptr);   pptr+=sizeof(int);
  int stride=*((int*)pptr);       pptr+=sizeof(int);
  double *ggamma=*((double**)pptr);   pptr+=sizeof(double*);
  double *gtheta=*((double**)pptr); pptr+=sizeof(double*);
  double *galpha=*((double**)pptr);   pptr+=sizeof(double*);
  double sgamma=*((double*)pptr);     pptr+=sizeof(double);
  double stheta=*((double*)pptr);   pptr+=sizeof(double);
  double salpha=*((double*)pptr);     pptr+=sizeof(double);
  double SSD=*((double*)pptr);

  int i, idx;
  double r=wald_pdf(t-SSD, salpha, sgamma, stheta);

  for(i=0; i<nresponses; i++ ){
	 idx=i*stride;
	 r*=(1-wald_cdf(t, galpha[idx], ggamma[idx], gtheta[idx]));
  }
  return r;
}

/** Returns raw likelihoods for each trial in pointer L.

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
 */
void sswald_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
			  					   int *condition, int *response, double *RT, double *SSD, /* data */
									double *go_gamma, double *go_theta, double *go_alpha,
									double *stop_gamma, double *stop_theta, double *stop_alpha, 
									double *pgf, double *ptf, double *L ){
  /*
  dprintf("nconditions=%i, nres=%i, ntri=%i\n", nconditions, nresponses, ntrials);
  for( i=0; i<ntrials; i++ ){
	 dprintf("condition[%i]=%i, response=%i, RT=%f, SSD=%f\n", i, condition[i], response[i], RT[i], SSD[i]);
  }
  dprintf("L=%f\n",*L);
  */

  int i, j, idx;

  /* setup for numerical integration of pstop */
  gsl_integration_workspace * work =gsl_integration_workspace_alloc (1000);
  double result, error;
  int nerror=0;
  int errcode;
  int neval;
  gsl_function F;
  F.function = &pstop_integrate_wald;
  struct pstop_params_wald params;
  F.params=&params;
  params.nresponses=nresponses;
  params.stride=nconditions;
  double pstop;
  gsl_set_error_handler_off();

  double dens, densgo, densstop, tmp;

  for( i=0; i<ntrials; i++ ){ /* calc L for each datapoint */
	 /* failed GO */
	 if( isnan(SSD[i]) && response[i]<0 ){
		L[i]=pgf[ condition[i] ];
	 }

	 /* successful stop */
	 else if( response[i]<0 && isfinite(SSD[i]) ){
		//dprintf("succstop, trial=%i\n",i);
		/* set parameters for integration */
		params.ggamma=&(go_gamma[condition[i]]);
		params.gtheta=&(go_theta[condition[i]]);
		params.galpha=&(go_alpha[condition[i]]);
		
		params.sgamma=(stop_gamma[condition[i]]);
		params.stheta=(stop_theta[condition[i]]);
		params.salpha=(stop_alpha[condition[i]]);
		params.SSD = SSD[i];

		/*		errcode=gsl_integration_qagiu (&F, stop_ter[condition[i]]+SSD[i], 1e-5, 1e-5, 1000, work, &result, &error); */
		errcode=gsl_integration_qagiu (&F, SSD[i], 1e-5, 1e-5, 1000, work, &result, &error); 
		if(errcode){
		  if(nerror==0){
			 dprintf("ERROR during integration, errcode=%i, abserr=%f, result=%f\n", errcode, error, result);
		  }
		  nerror++;
		}
		//dprintf("result=%f, error=%f, neval=%i\n",result, error, (int)(work->size));
		pstop=result;
		L[i]=pgf[condition[i]] + (1-pgf[condition[i]])*(1-ptf[condition[i]])*pstop;
	 }

	 /* go-trials with response */
	 else if( response[i]>=0 && isnan(SSD[i]) ){
		dens=1.0;
		//double wald_pdf(double t, double alpha, double gamma, double theta){
		for( j=0; j<nresponses; j++ ){
		  idx=condition[i]+(j*nconditions);
		  if( j==response[i] ){
			 //			 dprintf("Trial %i, condition=%i, PDF target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 dens*=wald_pdf(RT[i], go_alpha[idx], go_gamma[idx], go_theta[idx]);
		  } else {
			 //			 dprintf("Trial %i, condition=%i, CDF non-target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 dens*=(1-wald_cdf(RT[i], go_alpha[idx], go_gamma[idx], go_theta[idx]));
		  }
		}
		if( (dens<0) || !isfinite(dens) )
		  dens=0.0;
		L[i]=(1-pgf[condition[i]])*dens;
	 }

	 /* STOP-trials with response 
		 STOP(t) : (1-pgf)*[ptf*Lg(t) + (1-ptf)*Ls(t)]
	  */
	 else if( response[i]>=0 && isfinite(SSD[i]) ){
		idx=condition[i];
		densgo=1.0;
		densstop=(1-wald_cdf(RT[i]-SSD[i], stop_alpha[idx], stop_gamma[idx], stop_theta[idx]));

		for( j=0; j<nresponses; j++ ){
		  idx=condition[i]+(j*nconditions);
		  if( j==response[i] ){
			 //			 dprintf("Trial %i, condition=%i, PDF target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 tmp=wald_pdf(RT[i], go_alpha[idx], go_gamma[idx], go_theta[idx]);
			 densstop*=tmp;
			 densgo*=tmp;
		  } else {
			 //			 dprintf("Trial %i, condition=%i, CDF non-target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 tmp=(1-wald_cdf(RT[i], go_alpha[idx], go_gamma[idx], go_theta[idx]));
			 densstop*=tmp;
			 densgo*=tmp;
		  }
		}
		
		L[i]=(1-pgf[condition[i]])*(ptf[condition[i]]*densgo+(1-ptf[condition[i]])*densstop);
	 }

	 /* invalid -> error message */
	 else {
		dprintf("ERROR: trial %i has not been processed by loglikelihood!\n",i);
	 }
  }

  /* clean up */
  gsl_integration_workspace_free(work);
}

/*---------------------------------------------------------------------------*/
/* VarWALD Model (only VarWald accs)
/*---------------------------------------------------------------------------*/

struct pstop_params_varwald {
  int nresponses; int stride;
  double *ggamma; double *gtheta; double *galpha; double *gA;
  double sgamma; double stheta; double salpha; double sA; double SSD;
};

double pstop_integrate_varwald( double t, void *params ){
  void *pptr=params;
  int nresponses=*((int*)pptr);   pptr+=sizeof(int);
  int stride=*((int*)pptr);       pptr+=sizeof(int);
  double *ggamma=*((double**)pptr);   pptr+=sizeof(double*);
  double *gtheta=*((double**)pptr); pptr+=sizeof(double*);
  double *galpha=*((double**)pptr);   pptr+=sizeof(double*);
  double *gA=*((double**)pptr);   pptr+=sizeof(double*);
  double sgamma=*((double*)pptr);     pptr+=sizeof(double);
  double stheta=*((double*)pptr);   pptr+=sizeof(double);
  double salpha=*((double*)pptr);     pptr+=sizeof(double);
  double sA=*((double*)pptr);     pptr+=sizeof(double);
  double SSD=*((double*)pptr);

  int i, idx;
  double r=varwald_pdf(t-SSD, salpha, sgamma, stheta, sA);

  for(i=0; i<nresponses; i++ ){
	 idx=i*stride;
	 r*=(1-varwald_cdf(t, galpha[idx], ggamma[idx], gtheta[idx], gA[idx]));
  }
  return r;
}

/** Returns raw likelihoods for each trial in pointer L.

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
 */
void ssvarwald_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
		   				      int *condition, int *response, double *RT, double *SSD, /* data */
							  double *go_gamma, double *go_theta, double *go_alpha, double *go_A,
							  double *stop_gamma, double *stop_theta, double *stop_alpha, double *stop_A,
							  double *pgf, double *ptf, double *L ){
  /*
  dprintf("nconditions=%i, nres=%i, ntri=%i\n", nconditions, nresponses, ntrials);
  for( i=0; i<ntrials; i++ ){
	 dprintf("condition[%i]=%i, response=%i, RT=%f, SSD=%f\n", i, condition[i], response[i], RT[i], SSD[i]);
  }
  dprintf("L=%f\n",*L);
  */

  int i, j, idx;

  /* setup for numerical integration of pstop */
  gsl_integration_workspace * work =gsl_integration_workspace_alloc (1000);
  double result, error;
  int nerror=0;
  int errcode;
  int neval;
  gsl_function F;
  F.function = &pstop_integrate_varwald;
  struct pstop_params_varwald params;
  F.params=&params;
  params.nresponses=nresponses;
  params.stride=nconditions;
  double pstop;
  gsl_set_error_handler_off();

  double dens, densgo, densstop, tmp;

  for( i=0; i<ntrials; i++ ){ /* calc L for each datapoint */
	 /* failed GO */
	 if( isnan(SSD[i]) && response[i]<0 ){
		L[i]=pgf[ condition[i] ];
	 }

	 /* successful stop */
	 else if( response[i]<0 && isfinite(SSD[i]) ){
		//dprintf("succstop, trial=%i\n",i);
		/* set parameters for integration */
		params.ggamma=&(go_gamma[condition[i]]);
		params.gtheta=&(go_theta[condition[i]]);
		params.galpha=&(go_alpha[condition[i]]);
		params.gA    =&(go_A[condition[i]]);

		params.sgamma=(stop_gamma[condition[i]]);
		params.stheta=(stop_theta[condition[i]]);
		params.salpha=(stop_alpha[condition[i]]);
		params.sA    =(stop_A[condition[i]]);
		params.SSD = SSD[i];

		/*		errcode=gsl_integration_qagiu (&F, stop_ter[condition[i]]+SSD[i], 1e-5, 1e-5, 1000, work, &result, &error); */
		errcode=gsl_integration_qagiu (&F, SSD[i], 1e-5, 1e-5, 1000, work, &result, &error);
		if(errcode){
		  if(nerror==0){
			 dprintf("ERROR during integration, errcode=%i, abserr=%f, result=%f\n", errcode, error, result);
		  }
		  nerror++;
		}
		//dprintf("result=%f, error=%f, neval=%i\n",result, error, (int)(work->size));
		pstop=result;
		L[i]=pgf[condition[i]] + (1-pgf[condition[i]])*(1-ptf[condition[i]])*pstop;
	 }

	 /* go-trials with response */
	 else if( response[i]>=0 && isnan(SSD[i]) ){
		dens=1.0;

		for( j=0; j<nresponses; j++ ){
		  idx=condition[i]+(j*nconditions);
		  if( j==response[i] ){
			 //			 dprintf("Trial %i, condition=%i, PDF target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 dens*=varwald_pdf(RT[i], go_alpha[idx], go_gamma[idx], go_theta[idx], go_A[idx]);
		  } else {
			 //			 dprintf("Trial %i, condition=%i, CDF non-target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 dens*=(1-varwald_cdf(RT[i], go_alpha[idx], go_gamma[idx], go_theta[idx], go_A[idx]));
		  }
		}
		if( (dens<0) || !isfinite(dens) )
		  dens=0.0;
		L[i]=(1-pgf[condition[i]])*dens;
	 }

	 /* STOP-trials with response
		 STOP(t) : (1-pgf)*[ptf*Lg(t) + (1-ptf)*Ls(t)]
	  */
	 else if( response[i]>=0 && isfinite(SSD[i]) ){
		idx=condition[i];
		densgo=1.0;
		densstop=(1-varwald_cdf(RT[i]-SSD[i], stop_alpha[idx], stop_gamma[idx], stop_theta[idx], stop_A[idx]));

		for( j=0; j<nresponses; j++ ){
		  idx=condition[i]+(j*nconditions);
		  if( j==response[i] ){
			 //			 dprintf("Trial %i, condition=%i, PDF target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 tmp=varwald_pdf(RT[i], go_alpha[idx], go_gamma[idx], go_theta[idx], go_A[idx]);
			 densstop*=tmp;
			 densgo*=tmp;
		  } else {
			 //			 dprintf("Trial %i, condition=%i, CDF non-target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 tmp=(1-varwald_cdf(RT[i], go_alpha[idx], go_gamma[idx], go_theta[idx], go_A[idx]));
			 densstop*=tmp;
			 densgo*=tmp;
		  }
		}

		L[i]=(1-pgf[condition[i]])*(ptf[condition[i]]*densgo+(1-ptf[condition[i]])*densstop);
	 }

	 /* invalid -> error message */
	 else {
		dprintf("ERROR: trial %i has not been processed by loglikelihood!\n",i);
	 }
  }

  /* clean up */
  gsl_integration_workspace_free(work);
}



/*---------------------------------------------------------------------------*/
/* LOGNORMAL Model (only LN accs)
/*---------------------------------------------------------------------------*/


struct pstop_params_lognorm {
  int nresponses; int stride;
  double *gter; double *gmu; double *gsigma;
  double ster; double smu; double ssigma; double SSD;
};

double pstop_integrate_lognorm( double t, void *params ){
  void *pptr=params;
  int nresponses=*((int*)pptr);   pptr+=sizeof(int);
  int stride=*((int*)pptr);       pptr+=sizeof(int);
  double *gter=*((double**)pptr);   pptr+=sizeof(double*);
  double *gmu=*((double**)pptr); pptr+=sizeof(double*);
  double *gsigma=*((double**)pptr);   pptr+=sizeof(double*);
  double ster=*((double*)pptr);     pptr+=sizeof(double);
  double smu=*((double*)pptr);   pptr+=sizeof(double);
  double ssigma=*((double*)pptr);     pptr+=sizeof(double);
  double SSD=*((double*)pptr);

  int i, idx;
  double r=slognorm_pdf(t-SSD, ster, smu, ssigma);

  for(i=0; i<nresponses; i++ ){
	 idx=i*stride;
	 r*=(1-slognorm_cdf(t, gter[idx], gmu[idx], gsigma[idx]));
  }
  return r;
}


void sslognorm_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
			    			  int *condition, int *response, double *RT, double *SSD, /* data */
							  double *go_ter, double *go_mu, double *go_sigma,
							  double *stop_ter, double *stop_mu, double *stop_sigma,
							  double *pgf, double *ptf, double *L ){
  int i, j, idx;

  /* setup for numerical integration of pstop */
  gsl_integration_workspace * work =gsl_integration_workspace_alloc (1000);
  double result, error;
  int nerror=0;
  int errcode;
  int neval;
  gsl_function F;
  F.function = &pstop_integrate_lognorm;
  struct pstop_params_lognorm params;
  F.params=&params;
  params.nresponses=nresponses;
  params.stride=nconditions;
  double pstop;
  gsl_set_error_handler_off();

  double dens, densgo, densstop, tmp;

  for( i=0; i<ntrials; i++ ){ /* calc L for each datapoint */
	 /* failed GO */
	 if( isnan(SSD[i]) && response[i]<0 ){
		L[i]=pgf[ condition[i] ];
	 }

	 /* successful stop */
	 else if( response[i]<0 && isfinite(SSD[i]) ){
		//dprintf("succstop, trial=%i\n",i);
		/* set parameters for integration */
		params.gter  =&(go_ter  [condition[i]]);
		params.gmu   =&(go_mu   [condition[i]]);
		params.gsigma=&(go_sigma[condition[i]]);

		params.ster  =(stop_ter  [condition[i]]);
		params.smu   =(stop_mu   [condition[i]]);
		params.ssigma=(stop_sigma[condition[i]]);
		params.SSD = SSD[i];

		/*		errcode=gsl_integration_qagiu (&F, stop_ter[condition[i]]+SSD[i], 1e-5, 1e-5, 1000, work, &result, &error); */
		errcode=gsl_integration_qagiu (&F, SSD[i], 1e-5, 1e-5, 1000, work, &result, &error);
		if(errcode){
		  if(nerror==0){
			 dprintf("ERROR during integration, errcode=%i, abserr=%f, result=%f\n", errcode, error, result);
		  }
		  nerror++;
		}
		//dprintf("result=%f, error=%f, neval=%i\n",result, error, (int)(work->size));
		pstop=result;
		L[i]=pgf[condition[i]] + (1-pgf[condition[i]])*(1-ptf[condition[i]])*pstop;
	 }

	 /* go-trials with response */
	 else if( response[i]>=0 && isnan(SSD[i]) ){
		dens=1.0;
        //double slognorm_pdf(double t, double ter, double mu, double sigma);
		for( j=0; j<nresponses; j++ ){
		  idx=condition[i]+(j*nconditions);
		  if( j==response[i] ){
			 //			 dprintf("Trial %i, condition=%i, PDF target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 dens*=slognorm_pdf(RT[i], go_ter[idx], go_mu[idx], go_sigma[idx]);
		  } else {
			 //			 dprintf("Trial %i, condition=%i, CDF non-target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 dens*=(1-slognorm_cdf(RT[i], go_ter[idx], go_mu[idx], go_sigma[idx]));
		  }
		}
		if( (dens<0) || !isfinite(dens) )
		  dens=0.0;
		L[i]=(1-pgf[condition[i]])*dens;
	 }

	 /* STOP-trials with response
		 STOP(t) : (1-pgf)*[ptf*Lg(t) + (1-ptf)*Ls(t)]
	  */
	 else if( response[i]>=0 && isfinite(SSD[i]) ){
		idx=condition[i];
		densgo=1.0;
		densstop=(1-slognorm_cdf(RT[i]-SSD[i], stop_ter[idx], stop_mu[idx], stop_sigma[idx]));

		for( j=0; j<nresponses; j++ ){
		  idx=condition[i]+(j*nconditions);
		  if( j==response[i] ){
			 //			 dprintf("Trial %i, condition=%i, PDF target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 tmp=slognorm_pdf(RT[i], go_ter[idx], go_mu[idx], go_sigma[idx]);
			 densstop*=tmp;
			 densgo*=tmp;
		  } else {
			 //			 dprintf("Trial %i, condition=%i, CDF non-target-response=%i, idx=%i\n", i, condition[i], j, idx);
			 tmp=(1-slognorm_cdf(RT[i], go_ter[idx], go_mu[idx], go_sigma[idx]));
			 densstop*=tmp;
			 densgo*=tmp;
		  }
		}

		L[i]=(1-pgf[condition[i]])*(ptf[condition[i]]*densgo+(1-ptf[condition[i]])*densstop);
	 }

	 /* invalid -> error message */
	 else {
		dprintf("ERROR: trial %i has not been processed by loglikelihood!\n",i);
	 }
  }

  /* clean up */
  gsl_integration_workspace_free(work);
}
