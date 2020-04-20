#!/usr/bin/env python
from os.path import join
from setuptools import setup
import versionpy

import oncoviruses
pkg = oncoviruses.__name__

version = versionpy.get_version(pkg)

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
    packages=[pkg],
    include_package_data=True,
    zip_safe=False,
    install_requires=versionpy.get_reqs(),
    scripts=[
        join('scripts', 'oncoviruses'),
        # legacy:
        join('scripts', 'filter_bam_bridging_reads.py'),
        join('scripts', 'filter_bam_good_reads.py'),
        join('scripts', 'filter_bam_good_unmapped_or_mate.py'),
        join('scripts', 'filter_fastq.py'),
        join('scripts', 'filter_viral_bam.py'),
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
