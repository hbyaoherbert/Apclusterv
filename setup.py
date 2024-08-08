#!/usr/bin/env python
#-*- coding:utf-8 -*-
import os 
from setuptools import setup, find_packages

MAJOR =1
MINOR =2
PATCH =2
VERSION = f"{MAJOR}.{MINOR}.{PATCH}"

def get_install_requires():
    reqs = [
            'pandas>=0.18.0',
			'numpy>=1.9.2',
            'networkx>=2.8.4',
            'scipy>=1.8.1',
            'scikit-learn>=1.1.2'
            ]
    return reqs
setup(
	name = "apclusterv",
	version = VERSION,
    author ="haobinherbert",
    author_email = "haobinherbert@163.com",
    long_description_content_type="text/markdown",
	url = 'https://github.com/hbyaoherbert/Apclusterv.git',
	long_description = open('README.md',encoding="utf-8").read(),
    python_requires=">=3.7",
    install_requires=get_install_requires(),
	packages = find_packages(),
    entry_points={
        "console_scripts" : ['apclusterv = apclusterv.cluster_from_prot:main','prepare=apclusterv.prepare_from_contig:main']
    },
 	license = 'Apache',
   	classifiers = [
       'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
       'Natural Language :: English',
       'Operating System :: Unix',
       'Programming Language :: Python',       
       'Programming Language :: Python :: 3.7',
       'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    package_data={'': ['scripts/*.R', '*.txt','.toml']}, #这个很重要
    include_package_data=True #也选上

)
