# -*- coding: utf-8 -*-
from setuptools import setup, find_packages, Extension
import sys
import os
import traceback
import functools


def has_option(name):
    try:
        sys.argv.remove('--%s' % name)
        return True
    except ValueError:
        pass
    # allow passing all cmd line options also as environment variables
    env_val = os.getenv(name.upper().replace('-', '_'), 'false').lower()
    if env_val == "true":
        return True
    return False


include_diagnostics = has_option("include-diagnostics")
force_cythonize = has_option("force-cythonize")


try:
    from Cython.Build import cythonize
    cython_directives = {
        'embedsignature': True,
        "profile": include_diagnostics
    }
    if include_diagnostics:
        Extension = functools.partial(Extension, define_macros=[
            ("CYTHON_TRACE_NOGIL", "1"),
        ])
    extensions = cythonize([
        Extension(
            name="brainpy._speedup", sources=["brainpy/_speedup.pyx"],
            include_dirs=["brainpy/_c/"]),
        Extension(
            name="brainpy._c.composition", sources=["brainpy/_c/composition.pyx"],
            include_dirs=["brainpy/_c/"]),
        Extension(
            name="brainpy._c.double_vector", sources=["brainpy/_c/double_vector.pyx"],
            include_dirs=["brainpy/_c/"]),
        Extension(
            name="brainpy._c.isotopic_constants", sources=["brainpy/_c/isotopic_constants.pyx"],
            include_dirs=["brainpy/_c/"]),
        Extension(
            name="brainpy._c.isotopic_distribution", sources=["brainpy/_c/isotopic_distribution.pyx"],
            include_dirs=["brainpy/_c/"])
    ], compiler_directives=cython_directives, force=force_cythonize)
except ImportError:
    extensions = ([
        Extension(
            name="brainpy._speedup", sources=["brainpy/_speedup.c"],
            include_dirs=["brainpy/_c/"]),
        Extension(
            name="brainpy._c.composition", sources=["brainpy/_c/composition.c"],
            include_dirs=["brainpy/_c/"]),
        Extension(
            name="brainpy._c.double_vector", sources=["brainpy/_c/double_vector.c"],
            include_dirs=["brainpy/_c/"]),
        Extension(
            name="brainpy._c.isotopic_constants", sources=["brainpy/_c/isotopic_constants.c"],
            include_dirs=["brainpy/_c/"]),
        Extension(
            name="brainpy._c.isotopic_distribution", sources=["brainpy/_c/isotopic_distribution.c"],
            include_dirs=["brainpy/_c/"])
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
        name='brain-isotopic-distribution',
        version='1.5.7rc6',
        packages=find_packages(),
        description="Fast and efficient theoretical isotopic profile generation",
        long_description='''
:mod:`brainpy` is a small Python library implementing the *B* afflingly *R* ecursive
*A* lgorithm for *I* sotopic Patter *N* generation [Dittwald2014]_. It includes three implementations,
a pure-Python object oriented implementation, a :title-reference:`Cython` accelerated
version of the object oriented implementation, and a pure :title-reference:`C` implementation,
listed in order of ascending speed. The C implementation is used by default when available.

BRAIN, implemented in :func:`brainpy.isotopic_variants`, takes an elemental
composition represented by  any :class:`~.collections.abc.Mapping`-like Python object and uses
it to compute its aggregated isotopic distribution. All isotopic variants of the same number
of neutrons are collapsed into a single centroid peak, meaning it does not consider isotopic
fine structure.

Installing
----------
:mod:`brainpy` has three implementations, a pure Python implementation, a Cython translation
of that implementation, and a pure C implementation that releases the :title-reference:`GIL`.

To install from a package index, you will need to have a C compiler appropriate to your Python
version to build these extension modules. Additionally, there are prebuilt wheels for Windows
available on `PyPI <https://pypi.org/project/brain-isotopic-distribution/>`_:

.. code-block:: sh

    $ pip install brain-isotopic-distribution

To build from source, in addition to a C compiler you will also need to install a recent version
of `Cython <https://pypi.org/project/Cython/>`_ to transpile C code.

References
----------

This package is an implementation of the algorithm originally described in
    P. Dittwald, J. Claesen, T. Burzykowski, D. Valkenborg, and A. Gambin,
    "BRAIN: a universal tool for high-throughput calculations of the isotopic distribution for mass spectrometry.",
    Anal. Chem., vol. 85, no. 4, pp. 1991–4, Feb. 2013.

    H. Hu, P. Dittwald, J. Zaia, and D. Valkenborg,
    "Comment on 'Computation of isotopic peak center-mass distribution by fourier transform'.",
    Anal. Chem., vol. 85, no. 24, pp. 12189–92, Dec. 2013.
    ''',
        long_description_content_type='text/markdown',
        author=', '.join(["Joshua Klein", "Han Hu"]),
        author_email="jaklein@bu.edu",
        url="https://github.com/mobiusklein/brainpy",
        maintainer='Joshua Klein',
        keywords=["isotopic distribution", "isotopic pattern"],
        maintainer_email="jaklein@bu.edu",
        ext_modules=extensions if include_cext else None,
        include_package_data=True,
        cmdclass=cmdclass,
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
        project_urls={
            "Bug Tracker": "https://github.com/mobiusklein/brainpy/issues",
            "Source Code": "https://github.com/mobiusklein/brainpy",
            "Documentation": "http://mobiusklein.github.io/brainpy",
        }
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
