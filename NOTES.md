Problem when running `python setup.py install` from fresh git-checkout
----------------------------------------------------------------------

Because `crace.py` is only generated in the call to swig, it is not
found by `build_py` which is executed before `build_ext` (in which it
is created).

Solution: always run

    python setup.py build
	python setup.py install

instead of `python setup.py install` alone.


General Notes:
--------------

* It's awkward to have the ptf/pgf parameters as Model parameters rather than accumulator parameters
  - it would be much more elegant to assign the ptf/pgf to each accumulator (a probability of it starting)
    and let the model-definition layer deal with the mapping from pgf/ptf to those.
  - would resolve a lot of inconsistencies in the code