%
% ------------- pyrace.i
%
 
%define DOCSTRING
"
This is a docstring
"
%enddef

%module(docstring=DOCSTRING) pyrace
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


double sslba_loglikelihood( int nconditions, int nresponses, int ntrials,           /* global pars */
									 int *condition, int *response, double *RT, double *SSD, /* data */
									 


%clear (double* );
