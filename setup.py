from setuptools import setup, find_packages, Extension
import sys
import traceback

try:
    from Cython.Build import cythonize
    extensions = cythonize([
        Extension(name="brainpy._speedup", sources=["brainpy/_speedup.pyx"]),
        ])
except ImportError:
    extensions = ([
        Extension(name="brainpy._speedup", sources=["brainpy/_speedup.c"]),
        ])


from distutils.command.build_ext import build_ext
from distutils.errors import (CCompilerError, DistutilsExecError,
                              DistutilsPlatformError)

ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError)
if sys.platform == 'win32':
    # 2.6's distutils.msvc9compiler can raise an IOError when failing to
    # find the compiler
    ext_errors += (IOError,)


class BuildFailed(Exception):

    def __init__(self):
        self.cause = sys.exc_info()[1]  # work around py 2/3 different syntax

    def __str__(self):
        return str(self.cause)


class ve_build_ext(build_ext):
    # This class allows C extension building to fail.

    def run(self):
        try:
            build_ext.run(self)
        except DistutilsPlatformError:
            traceback.print_exc()
            raise BuildFailed()

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except ext_errors:
            traceback.print_exc()
            raise BuildFailed()
        except ValueError:
            # this can happen on Windows 64 bit, see Python issue 7511
            traceback.print_exc()
            if "'path'" in str(sys.exc_info()[1]):  # works with both py 2/3
                raise BuildFailed()
            raise

cmdclass = {}

cmdclass['build_ext'] = ve_build_ext


def status_msgs(*msgs):
    print('*' * 75)
    for msg in msgs:
        print(msg)
    print('*' * 75)


def run_setup(include_cext=True):
    setup(
        name='brainpy',
        version='1.0.8',
        packages=find_packages(),
        description="Fast and efficient theoretical isotopic profile generation",
        long_description='''
        A Python Implementation of the Baffling Recursive Algorithm for Isotopic cluster distributioN
    ''',
        author=', '.join(["Joshua Klein", "Han Hu"]),
        author_email=["jaklein@bu.edu"],
        ext_modules=extensions if include_cext else None,
        cmdclass=cmdclass,
        classifiers=[
                'Development Status :: 4 - Beta',
                'Intended Audience :: Science/Research',
                'License :: OSI Approved :: BSD License',
                'Topic :: Scientific/Engineering :: Bio-Informatics']
    )

if __name__ == '__main__':
    try:
        run_setup(True)
    except Exception as exc:
        run_setup(False)

        status_msgs(
            "WARNING: The C extension could not be compiled, " +
            "speedups are not enabled.",
            "Plain-Python build succeeded."
        )
