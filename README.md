pyrace
======

Race-Models for RT-data in python

This is work-in-progress and currently only for Stop-Tasks.

Mac OS X compilation
--------------------

Install libgsl (e.g. brew install gsl).

Use 

    ARCHFLAGS="-arch x86_64" python setup.py build 
	
or whatever architecture you used for compiling libgsl.
