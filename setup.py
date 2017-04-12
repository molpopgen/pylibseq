#!/usr/bin/env python

from __future__ import print_function
from setuptools import setup,Extension

import platform, glob, sys, subprocess,os

#Dependency checks

##Check for libsequence version
try:
    out = subprocess.check_output('libsequenceConfig --version',shell=True)
    version = out.decode('utf-8').rstrip()
    print ("libsequence version",version," found.")
    if version < '1.9.0':
        print("libsequence >= ,'1.9.0' required, but ",version, "found.")
        sys.exit(2)
except:
    sys.exit("libsequenceConfig not found.  Please install libsequence (http://github.com/molpopgen/libsequence)")

##Can we compile a program based on libsequence?
print("Attempting to compile and link a test program using libsequence...")
try:
    proc = subprocess.check_output(['make','-f','check_deps/Makefile','clean'])
    proc = subprocess.check_output(['make','-f','check_deps/Makefile'])
    proc = subprocess.check_output(['make','-f','check_deps/Makefile','clean'])
except subprocess.CalledProcessError as e:
    print (e.returncode)
    sys.exit(2)
print("done")

#Are we gonna build using Cython or not?  Default is not to,
#which allows us to ship this in a standard way.
if '--use-cython' in sys.argv:
    USE_CYTHON = True
    sys.argv.remove('--use-cython')
else:
    USE_CYTHON = False


long_desc = open("README.rst").read()

EXTENSION = '.pyx' if USE_CYTHON else '.cpp'

extensions = []
pdata = {'libsequence':['*.pxd']}
provided = []
modules = ['polytable','summstats','windows','fst','extensions','parallel']
for i in modules:
    LIBS=["sequence"]
    if i == 'parallel':
        LIBS.append('tbb')
    extensions.append(Extension("libsequence."+i,
                                sources=["libsequence/"+i+EXTENSION],
                                language="c++",                  
                                extra_compile_args=["-std=c++11"],  
                                extra_link_args=["-std=c++11"],
                                libraries=LIBS,
                                ))
    provided.append('libsequence.'+i)
    pdata['libsequence.'+i]=['*.pxd']

#If using Cython, edit extensions here:
if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)


setup(name='pylibseq',
      version='0.1.9',
      author='Kevin R. Thornton',
      author_email='krthornt@uci.edu',
      maintainer='Kevin R. Thornton',
      maintainer_email='krthornt@uci.edu',
      url='http://github.com/molpopgen/pylibseq',
      description="Python interface to libsequence",
      long_description=long_desc,
      data_files=[('pylibseq', ['COPYING', 'README.rst'])],
      download_url='',
      classifiers=[],
      platforms=['Linux','OS X'],
      license='GPL >= 2',
      provides=provided,
      obsoletes=['none'],
      packages=['libsequence','libsequence.console'],
      py_modules=[],
      scripts=[],
      package_data=pdata,
      ext_modules=extensions,
      entry_points={
          'console_scripts':[
              'pymsstats=libsequence.console.msstats_cli:msstats_main'
              ]
          },
)
     
