from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import subprocess
import os
import glob

if sys.version_info < (3, 4):
    raise RuntimeError("Python >= 3.4 required")

try:
    libseq_version = subprocess.run(
        ['libsequenceConfig', '--version'], stdout=subprocess.PIPE)
except subprocess.CalledProcessError as error:
    print("Fatal error:", error)

if libseq_version.stdout.decode('utf8').rstrip() != "1.9.7":
    raise ValueError("libsequence == " + '1.9.7' + "required")

__version__ = '0.2.1.post0'

# clang/llvm is default for OS X builds.
# can over-ride darwin-specific options
# with setup.py --gcc install
if '--gcc' in sys.argv:
    USE_GCC = True
    sys.argv.remove('--gcc')
else:
    USE_GCC = False


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


PKGS = ['libsequence']

INCLUDES = [
    # Path to pybind11 headers
    get_pybind_include(),
    get_pybind_include(user=True),
    os.path.join(sys.prefix, 'include')
]

LIBRARY_DIRS = [
    os.path.join(sys.prefix, 'lib')
]

ext_modules = [
    Extension(
        'libsequence.polytable',
        ['libsequence/src/polytable.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['sequence'],
        language='c++'
    ),
    Extension(
        'libsequence.summstats',
        ['libsequence/src/summstats.cc','libsequence/src/omega_max.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['sequence'],
        language='c++'
    ),
    Extension(
        'libsequence.fst',
        ['libsequence/src/fst.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['sequence'],
        language='c++'
    ),
    Extension(
        'libsequence.windows_cpp',
        ['libsequence/src/windows_cpp.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['sequence'],
        language='c++'
    ),
    Extension(
        'libsequence.variant_matrix',
        ['libsequence/src/variant_matrix.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['sequence'],
        language='c++'
    ),

]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin' and USE_GCC is False:
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' %
                        self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden') and (sys.platform != 'darwin' or USE_GCC is True):
                opts.append('-fvisibility=hidden')
            if has_flag(self.compiler, '-g0'):
                opts.append('-g0')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' %
                        self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
            if sys.platform == 'darwin' and USE_GCC is False:
                ext.extra_link_args = [
                    '-stdlib=libc++', '-mmacosx-version-min=10.7']
        build_ext.build_extensions(self)


long_desc = open("README.rst").read()

setup(
    name='pylibseq',
    version=__version__,
    author='Kevin Thornton',
    author_email='krthornt@uci.edu',
    url='http://molpopgen.github.io/pylibseq',
    classifiers=['Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'],
    description='Python interface to libsequence.',
    license='GPL >= 3',
    requires=['pybind11', 'numpy'],
    provides=['pylibseq'],
    obsoletes=['none'],
    data_files=[('pylibseq', ['COPYING', 'README.rst'])],
    long_description=long_desc,
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.2.3', 'msprime>=0.5.0'],
    cmdclass={'build_ext': BuildExt},
    packages=PKGS,
    entry_points={
        'console_scripts': ['pymsstats = libsequence.msstats_cli:msstats_main']
    },
    zip_safe=False,
)
