#include <stdio.h>
#include <math.h>

/** WARNING:

	 Gotta be really careful to include the correct GSL headers.
	 If you don't include them properly for any function you use, 
	 you will not get any warning but instead meaningless results... weird...

 */
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

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



double lba_pdf(double t, double ter, double A, double v, double sv, double b){
  /* TESTED */
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
  
  return MIN( 1, MAX( 0, (1+(tmp1+tmp2)/A)/pnormP(v/sv,0,1)));
}



/** Returns raw likelihoods for each trial in pointer L.
 */
void sslba_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
								  int *condition, int *response, double *RT, double *SSD, /* data */
								  double *go_v, double *go_ter, double *go_A, double *go_b, double *go_sv,
								  double *stop_v, double *stop_ter, double *stop_A, double *stop_b, double *stop_sv,
								  double *pgf, double *ptf, double *L ){
  int i;
  dprintf("nconditions=%i, nres=%i, ntri=%i\n", nconditions, nresponses, ntrials);
  for( i=0; i<ntrials; i++ ){
	 dprintf("condition[%i]=%i, response=%i, RT=%f, SSD=%f\n", i, condition[i], response[i], RT[i], SSD[i]);
  }
  dprintf("L=%f\n",*L);

  for( i=0; i<ntrials; i++ ){ /* calc L for each datapoint */
	 /* failed GO */
	 if( isnan(SSD[i]) && response[i]<0 ){
		L[i]=pgf[ condition[i] ];
	 }

	 /* successful stop */
	 else if( response[i]<0 && isfinite(SSD[i]) ){
		L[i]=pgf[condition[i]] + (1-pgf[condition[i]])*(1-ptf[condition[i]])*pstop;
	 }
  }

}
