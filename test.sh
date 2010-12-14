python setup.py build_ext -i
PYTHONPATH=. python bedtools/tests/test.py
PYTHONPATH=. python test.py > /dev/null
nosetests --with-doctest --doctest-extension=.rst .
