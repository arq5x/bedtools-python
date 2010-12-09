"""
    __Acknowledgments__: This module was made possible by 
       following Brent Pederson's very nice example of porting 
       FastaHack to Python using Cython.
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import glob

setup(
  name = 'bedtools',
  ext_modules=[
    Extension("bedtools.cbedtools",
              sources=["bedtools/cbedtools.pyx"] + glob.glob("src/*.cpp"),
              libraries=["stdc++", 'z'],
              include_dirs=["src/"],
              depends = glob.glob("src/*.h"),
              language="c++"),
    ],
    package_data = {'src': ['*.pyx', "*.c", "*.cpp", "*.h", "README.rst"]},
    package_dir = {"bedtools": "bedtools"},
    cmdclass = {'build_ext': build_ext},
    packages = ['bedtools'],
    author = "Aaron Quinlan",
    author_email="aaronquinlan@gmail.com",
    #test_suite='nose.collector'
)
