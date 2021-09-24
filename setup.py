#!/usr/bin/env python
import os

from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit tutorials/__version__.py
version = {}
with open(os.path.join(here, 'tutorials', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='tutorials',
    version=version['__version__'],
    description="A collection of tutorials for the analysis of nanocrystals.",
    long_description=readme + '\n\n',
    author=["Juliette Zito", "felipe zapata"],
    author_email='juliette.zito@hotmail.fr',
    url='https://github.com/nlesc-nano/Tutorials',
    packages=[
        'tutorials',
    ],
    include_package_data=True,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='tutorials',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Documentation :: Sphinx',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    test_suite='tests',
    install_requires=[],
    tests_require=[
        'pytest',
        'pytest-cov',
        'pycodestyle',
    ],
    extras_require={
        'docs':  ['sphinx', 'sphinx_rtd_theme>=1.0'],
        'test':  ['pytest', 'pytest-cov', 'pycodestyle'],
    }
)
