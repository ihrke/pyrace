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


void sslba_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
								  int *condition, int *response, double *RT, double *SSD, /* data */
								  double *go_v, double *go_ter, double *go_A, double *go_b, double *go_sv, /* GO-LBA params */
								  double *stop_v, double *stop_ter, double *stop_A, double *stop_b, double *stop_sv, /* STOP-LBA params */
								  double *pgf, double *ptf,  /* mixing probs */
								  double *L); /* output */


%clear (double* );
%clear (int* );
