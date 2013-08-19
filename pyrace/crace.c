#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>

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


double lba_pdf(double t, double ter, double A, double v, double sv, double b){
  if( A<1e-10 ){ /* LATER */
	 return MAX( 0, (b/(SQR(t)))*dnormP(b/t, v, sv)/pnormP(v/sv,0,1));
  }
  return 0.0;
}
double lba_cdf(double t, double ter, double A, double v, double sv, double b){
  gsl_ran_gaussian_pdf(t, 1.0);
  return 0.0;
}


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
  }

  *L=1.0;
}
