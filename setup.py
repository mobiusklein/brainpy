from setuptools import setup, find_packages

setup(
    name='brainpy',
    version='1.0.6',
    packages=find_packages(),
    description="Fast and efficient theoretical isotopic profile generation",
    long_description='''
    A Python Implementation of the Baffling Recursive Algorithm for Isotopic cluster distributioN
''',
    author=', '.join(["Joshua Klein", "Han Hu"]),
    author_email=["jaklein@bu.edu"],
    classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
)
