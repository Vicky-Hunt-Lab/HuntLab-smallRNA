#!/usr/bin/env python3

from distutils.core import setup

setup(
    name='HuntLab-smallRNA',
    version='1.0.0',
    description='Small RNA scripts developed for use in the Hunt Lab',
    author='Kieran Reynolds',
    packages=['hlsmallrna'],
    scripts=['bin/label_for_unitas', 'bin/build_coord_files'],
    entry_points = {
        'console_scripts': [
            'hlsmallrna = hlsmallrna:climain',                  
        ]          
    }
)