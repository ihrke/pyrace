Problem when running `python setup.py install` from fresh git-checkout
----------------------------------------------------------------------

Because `crace.py` is only generated in the call to swig, it is not
found by `build_py` which is executed before `build_ext` (in which it
is created).

Solution: always run

    python setup.py build
	python setup.py install

instead of `python setup.py install` alone.
