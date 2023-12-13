# Copyright 2022 Vicky Hunt Lab Members
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#!/usr/bin/env python3

from distutils.core import setup

setup(
    name='HuntLab-smallRNA',
    version='1.1.0',
    description='Small RNA scripts developed for use in the Hunt Lab',
    author='Vicky Hunt Lab',
    packages=['hlsmallrna'],
    scripts=['bin/label_for_unitas', 'bin/build_coord_files', 'bin/revcomp_rna', 'bin/overlap_ss.sh'],
    entry_points = {
        'console_scripts': [
            'hlsmallrna = hlsmallrna:climain',
            'overlap_ss = hlsmallrna:ssoverlap_main'               
        ]          
    }
)