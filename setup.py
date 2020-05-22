#!/usr/bin/env python
from os.path import join
from setuptools import setup, find_packages
from versionpy import get_version, find_package_files, get_reqs

pkg = 'oncoviruses'
version = get_version(pkg)

setup(
    name=pkg,
    script_name=pkg,
    version=version,
    author='Vlad Saveliev',
    author_email='vladislav.sav@gmail.com',
    description='Detection of oncoviruses in tumor whole genome samples',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url='https://github.com/umccr/' + pkg,
    license='GPLv3',
    package_data={
        pkg: find_package_files('', pkg),
    },
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=get_reqs(),
    scripts=[
        join('scripts', 'oncoviruses'),
    ],
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ],
)
