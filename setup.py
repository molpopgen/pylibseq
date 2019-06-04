from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import subprocess
import platform
import os
import glob

if sys.version_info < (3, 4):
    raise RuntimeError("Python >= 3.4 required")

__version__ = '0.2.2'

# clang/llvm is default for OS X builds.
# can over-ride darwin-specific options
# with setup.py --gcc install
if '--gcc' in sys.argv:
    USE_GCC = True
    sys.argv.remove('--gcc')
else:
    USE_GCC = False

if '--no-weffcpp' in sys.argv:
    USE_WEFFCPP = False
    sys.argv.remove('--no-weffcpp')
else:
    USE_WEFFCPP = True

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(
            self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        # cfg = 'Debug' if DEBUG_MODE is True else 'Release'
        cfg = 'Release'

        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += [
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())

        if USE_WEFFCPP is False:
            cmake_args.append('-DUSE_WEFFCPP=OFF')
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        # if SKIP_BUILDING_TESTS is True:
        #     cmake_args.append('-DBUILD_UNIT_TESTS=OFF')
        subprocess.check_call(['cmake', ext.sourcedir] +
                              cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] +
                              build_args, cwd=self.build_temp)


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
    CMakeExtension('libsequence/_libsequence'),
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


# Figure out the headers we need to install:
generated_package_data = {}
for root, dirnames, filenames in os.walk('libsequence/src/libsequence/Sequence'):
    if 'testsuite' not in root and 'Coalescent' not in root:
        g = glob.glob(root + '/*.hpp')
        if len(g) > 0:
            replace = root.replace('/',  '.')
            # If there's a header file, we add the directory as a package
            if replace not in PKGS:
                PKGS.append(replace)
            try:
                if '*.hpp' not in generated_package_data[replace]:
                    generated_package_data[replace].append('*.hpp')
            except:
                generated_package_data[replace] = ['*.hpp']
        g = glob.glob(root + '/*.tcc')
        if len(g) > 0:
            replace = root.replace('/', '.')
            # If there's a template implementation file,
            # we add the directory as a package
            if replace not in PKGS:
                PKGS.append(replace)
            try:
                if '*.tcc' not in generated_package_data[replace]:
                    generated_package_data[replace].append('*.tcc')
            except:
                generated_package_data[replace] = ['*.tcc']

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
    install_requires=['pybind11>=2.2.3'],
    cmdclass={'build_ext': CMakeBuild},
    packages=PKGS,
    package_data=generated_package_data,
    # entry_points={
    #     'console_scripts': ['pymsstats = libsequence.msstats_cli:msstats_main']
    # },
    zip_safe=False,
)
