

Mac OS X compilation
--------------------

Install libgsl (e.g. brew install gsl).

Use 

    ARCHFLAGS="-arch x86_64" python setup.py build 
	
or whatever architecture you used for compiling libgsl.
