python setup.py build_ext -i
PYTHONPATH=. python bedtools/tests/tests.py
PYTHONPATH=. python example-intersect.py > /dev/null
nosetests --with-doctest --doctest-extension=.rst .
