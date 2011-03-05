python setup.py build_ext -i
PYTHONPATH=. python bedtools/tests/tests.py
PYTHONPATH=. python ex1-intersect.py > /dev/null
nosetests --with-doctest --doctest-extension=.rst .
