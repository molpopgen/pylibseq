from distutils.core import setup, Extension
from Cython.Build import cythonize
import platform, glob, sys

extensions=[
     Extension("libsequence.polytable",
               sources=["libsequence/polytable.pyx"], # the Cython source and additional C++ source files
        language="c++",                        # generate and compile C++ code
        include_dirs=['.','include','..'], 
        extra_compile_args=["-std=c++11"],  
        extra_link_args=["-std=c++11"],
        libraries=["sequence"]),
     Extension("libsequence.summstats",
               sources=["libsequence/summstats.pyx"], # the Cython source and additional C++ source files
               language="c++",                        # generate and compile C++ code
               include_dirs=['.','include','..'], 
               extra_compile_args=["-std=c++11"],  
               extra_link_args=["-std=c++11"],
               libraries=["sequence"]),
    Extension("libsequence.windows",
               sources=["libsequence/windows.pyx"], # the Cython source and additional C++ source files
               language="c++",                        # generate and compile C++ code
               include_dirs=['.','include','..'], 
               extra_compile_args=["-std=c++11"],  
               extra_link_args=["-std=c++11"],
               libraries=["sequence"]),
]

setup(name='libsequence',
      version='0.0.1',      
      author='Kevin R. Thornton',
      author_email='krthornt@uci.edu',
      maintainer='Kevin R. Thornton',
      maintainer_email='krthornt@uci.edu',
      url='http://www.molpopgen.org',
      description="",
      download_url='',
      classifiers=['population genetics'],
      platforms=['Linux','OS X'],
      license='GPL >= 2',
      provides=['libsequence.polytable'],
      obsoletes=['none'],
      packages=['libsequence'],
      py_modules=[],
      scripts=[],
      package_data={'libsequence':['*.pxd'],'libsequence.polytable':['*.pxd'],'libsequence.summstats':['*.pxd'],'libsequence.windows':['*.pxd']},
      ext_modules=cythonize(extensions))

