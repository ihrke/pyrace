#ifndef CRACE_H
#define CRACE_H

void sslba_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
								  int *condition, int *response, double *RT, double *SSD, /* data */
								  double *go_v, double *go_ter, double *go_A, double *go_b, double *go_sv,
								  double *stop_v, double *stop_ter, double *stop_A, double *stop_b, double *stop_sv,
								  double *pgf, double *ptf, double *L );

void sswald_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
			  					   int *condition, int *response, double *RT, double *SSD, /* data */
									double *go_gamma, double *go_theta, double *go_alpha,
									double *stop_gamma, double *stop_theta, double *stop_alpha, 
									double *pgf, double *ptf, double *L );

double pnormP(double x, double mean, double sd);
double dnormP(double x, double mean, double sd);

double lba_pdf(double t, double ter, double A, double v, double sv, double b);
double lba_cdf(double t, double ter, double A, double v, double sv, double b);

double wald_pdf(double t, double alpha, double gamma, double theta);
double wald_cdf(double t, double alpha, double gamma, double theta);

double slognorm_pdf(double t, double ter, double mu, double sigma);
double slognorm_cdf(double t, double ter, double mu, double sigma);



#endif
