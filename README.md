pyrace
======

[![Build Status](https://travis-ci.org/ihrke/pyrace.png?branch=master)](https://travis-ci.org/ihrke/pyrace)

Race-Models for RT-data in python

This is work-in-progress and currently only for Stop-Tasks.

Have a look at these ipython-notebooks:

<https://github.com/ihrke/pyrace/tree/notebook_tests>

Mac OS X compilation
--------------------

Install libgsl (e.g. brew install gsl).

Use 

    ARCHFLAGS="-arch x86_64" python setup.py build 
	
or whatever architecture you used for compiling libgsl.



