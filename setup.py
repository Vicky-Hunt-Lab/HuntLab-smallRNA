#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='HuntLab-smallRNA',
    version='1.0.0',
    packages=[find_packages(include=['hlsmallrna', 'hlsmallrna.*'])],
    scripts=['bin/label_for_unitas'],
    entry_points = {
        'console_scripts': [
            'hlsmallrna = hlsmallrna:climain',                  
        ],              
    }
)