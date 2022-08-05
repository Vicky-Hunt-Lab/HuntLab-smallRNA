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

# This file is a simple set of unit tests to work as a sanity check, 
# worth running to make sure everything is installed correctly, but
# not comprihensive

import os
import shutil
import unittest

from hlsmallrna import config
from hlsmallrna import __main__ as commands

PATH_TO_TEST_FILE = os.path.dirname(os.path.realpath(__file__))

class ConfigTest(unittest.TestCase):

    def test_validate_toml(self):
        '''
        Tests for the validate TOML method, without the loading step
        '''
        TOML_MATCH = {
            'key1': {
                'key1_2': 'string',
                'key1_3': 1,
                'key1_4': [1, 2, 3]
            },
            'key2': {
                'key2_1': ['1', '2', '3']
            }
        }

        result = config.validate_toml(
            {
                'key1': {
                    'key1_2': 'world',
                    'key1_3': 1,
                    'key1_4': [7, 2, 3]
                },
                'key2': {
                    'key2_1': ['1', '2', '7']
                }
            },
            TOML_MATCH
        )
        self.assertIsNone(result)

        result = config.validate_toml(
            {
                'key1': {
                    'key1_2': 'world',
                    'key1_3': '1',
                    'key1_4': [7, 2, 3]
                },
                'key2': {
                    'key2_1': ['1', '2', '7']
                }
            },
            TOML_MATCH
        )
        self.assertEqual(result, "key1_3 should be <class 'int'> but found <class 'str'>")

        result = config.validate_toml(
            {
                'key1': {
                    'key1_2': 7,
                    'key1_3': 1,
                    'key1_4': [7, 2, 3]
                },
                'key2': {
                    'key2_1': ['1', '2', '7']
                }
            },
            TOML_MATCH
        )
        self.assertEqual(result, "key1_2 should be <class 'str'> but found <class 'int'>")

        # Current behavior doesn't care what type in in an array
        result = config.validate_toml(
            {
                'key1': {
                    'key1_2': 'world',
                    'key1_3': 1,
                    'key1_4': [7, 2, 3]
                },
                'key2': {
                    'key2_1': ['1', 2, '7']
                }
            },
            TOML_MATCH
        )
        self.assertIsNone(result)

        result = config.validate_toml(
            {
                'key1': {
                    'key1_2': 'world',
                    'key1_3': 1,
                    'key1_4': [7, 2, 3]
                },
                'key2': {
                    'key2_1': 54321
                }
            },
            TOML_MATCH
        )
        self.assertEqual(result, "key2_1 should be <class 'list'> but found <class 'int'>")

    def test_load_toml(self):
        '''
        Test the code that loads the TOML from a file
        '''
        valid_config_path = os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'fullyValidConfig.toml')
        invalid_config_path = os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'invalidConfig.toml')

        self.assertIsNone(config.load_config(valid_config_path))
        assert config.CONFIG == {
            'general': {
                'output': './outputArray',
                'threads': 8
            },
            'cli-tools': {
                'trim': {
                    'trim_pass_threads': True,
                    'path_to_trim': 'cutadapt',
                    'trim_type': '',
                    'trim_params': []
                },
                'fastqc': {
                    'fastqc_pass_threads': True,
                    'path_to_fastqc': 'fastqc',
                    'fastqc_params': []
                },
                'bbmap': {
                    'bbmap_pass_threads': True,
                    'path_to_bbmap': 'bbmap.sh',
                    'bbmap_align_params': [],
                    'bbmap_index_params': ['--hello', '--world', '7']
                }
            }
        }

        with self.assertRaises(Exception):
            config.load_config(invalid_config_path)

    def test_get_config_key(self):
        '''
        Test retriving keys from a config set
        '''
        valid_config_path = os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'fullyValidConfig.toml')
        self.assertIsNone(config.load_config(valid_config_path))

        self.assertEqual(config.get_config_key('general', 'output'), './outputArray')
        self.assertTrue(config.get_config_key('cli-tools', 'bbmap', 'bbmap_pass_threads'))
        assert config.get_config_key('cli-tools', 'bbmap', 'bbmap_index_params') == ['--hello', '--world', '7']

        self.assertEqual(config.get_config_key('cli-tools', 'samtools', 'path_to_samtools'), 'samtools')

        with self.assertRaises(KeyError):
            config.get_config_key('some', 'rubbish', 'key')

class StepTests(unittest.TestCase):

    # def test_process_step(self):
    #     pass

    def test_unitas_step(self):
        '''
        Test the unitas command
        '''
        os.mkdir('./output')

        valid_config_path = os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'fullTestConfig.toml')
        self.assertIsNone(config.load_config(valid_config_path))

        commands.unitas_command(
            os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'unitasSmallRNA'),
            'x',
            [os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'randomSmallRNAGenome.fasta')],
            2
        )

        with open('./output/unitas_summery.csv') as f:
            self.assertNotEqual(f.read(), '')
    
    def test_sort_step(self):
        '''
        Test the sort command
        '''
        EXPECTED_RNA_LENGTH_FILE = '''RNA Length,A,C,G,U,All Frequency
20,1942,1968,1930,1933,7773
21,1826,1710,1803,2180,7519
22,1772,1776,1773,2383,7704
23,1681,1672,1766,2700,7819
24,1945,1942,1820,1939,7646
25,1779,1898,1867,2158,7702
26,1821,1802,1827,2360,7810
27,1672,1603,1731,2519,7525
28,1954,1909,1912,2021,7796
'''
        os.mkdir('./output')

        valid_config_path = os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'fullTestConfig.toml')
        self.assertIsNone(config.load_config(valid_config_path))

        commands.sort_command(
            os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'randomSmallRNAGenome.fasta'),
            os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'randomSmallRNA.fastq'),
            20,
            28,
            2
        )

        with open('./output/rna_length_report.csv') as f:
            self.assertEqual(f.read(), EXPECTED_RNA_LENGTH_FILE)

    def test_target_step(self):
        '''
        Test the targetid command function
        '''
        os.mkdir('./output')
        commands.targetid_command(
            os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'unitasAndTargetRNA.fastq'),
            [os.path.join(PATH_TO_TEST_FILE, 'testFiles', 'targetTest.fasta')],
            20,
            0,
            2
        )

        with open('./output/rna_target_list.tsv') as f:
            self.assertEqual(len(f.readlines()), 4)

    def tearDown(self):
        '''
        Delete any copies of the ouput directory
        '''
        shutil.rmtree('./output')

if __name__ == '__main__':
    unittest.main()