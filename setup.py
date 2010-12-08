"""
    __Acknowledgments__: This module was made possible by 
       following Brent Pederson's very nice example of porting 
       FastaHack to Python using Cython.
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'bedtools-python',
  ext_modules=[
    Extension("btpython",
              sources=["bedtools-python/bedtools-python.pyx", "src/bedFile.cpp"],
              libraries=["stdc++"],
              include_dirs=["src/"],
              language="c++"),
    ],
    package_data = {'src': ['*.pyx', "*.c", "*.h", "README.rst"]},
    package_dir = {"bedtools-python": "bedtools-python"},
    cmdclass = {'build_ext': build_ext},
    packages = ['bedtools-python'],
    author = "Aaron Quinlan",
    author_email="aaronquinlan@gmail.com",
    #test_suite='nose.collector'
)