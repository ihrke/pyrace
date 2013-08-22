%define DOCSTRING
"
This is a C-module speeding up some of pyrace's functions.
"
%enddef

%module(docstring=DOCSTRING) crace
%{
#define SWIG_FILE_WITH_INIT
%}

%init %{
  import_array();
%}

%header %{
#include "numpy/arrayobject.h"
#include "crace.h"
%}

%feature("autodoc","1");

/** put in raw pointer
*/
%typemap(in) (double*){
  $1 = (double*)PyArray_DATA($input);
}
%typemap(in) (int*){
  $1 = (int*)PyArray_DATA($input);
}

void init(void);

void sslba_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
								  int *condition, int *response, double *RT, double *SSD, /* data */
								  double *go_v, double *go_ter, double *go_A, double *go_b, double *go_sv, /* GO-LBA params */
								  double *stop_v, double *stop_ter, double *stop_A, double *stop_b, double *stop_sv, /* STOP-LBA params */
								  double *pgf, double *ptf,  /* mixing probs */
								  double *L); /* output */

/* just exposing for debugging purposes */
double pnormP(double x, double mean, double sd);
double dnormP(double x, double mean, double sd);
double lba_pdf(double t, double ter, double A, double v, double sv, double b);
double lba_cdf(double t, double ter, double A, double v, double sv, double b);

%clear (double* );
%clear (int* );
