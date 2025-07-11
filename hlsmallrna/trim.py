# Copyright 2022 - 2025 Vicky Hunt Lab Members
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
import os
import sys
import shutil

from subprocess import run

def run_trim(
        file_to_trim, output, five_prime_adapter=None, three_prime_adapter=None, quality_cutoff=20, threads=4,
        cutadapt_path='cutadapt', keep_intermediate=False, verbose=False
    ):
    '''
    Call cutadapt to remove adapters
    '''
    current_input = file_to_trim
    files_to_remove = []

    # quality filter
    if quality_cutoff > 0:
        print('====> Removing low quality reads...', file=sys.stderr)
        cutadapt_quality = [
            cutadapt_path, '-q', str(quality_cutoff), '-o', os.path.join(output, 'quality_filtered.fq'), '-j', str(threads), current_input
        ]
        run(cutadapt_quality, capture_output=not verbose)

        current_input = os.path.join(output, 'quality_filtered.fq')
        files_to_remove.append(current_input)

    # trim five prime adapters
    if five_prime_adapter is not None:
        print('====> Removing complete 5\' adapters...', file=sys.stderr)
        cutadapt_five_prime = [
            cutadapt_path, '-g', five_prime_adapter, '-o', os.path.join(output, 'five_prime_removed.fq'), '-j', str(threads), current_input
        ]
        run(cutadapt_five_prime, capture_output=not verbose)

        current_input = os.path.join(output, 'five_prime_removed.fq')
        files_to_remove.append(current_input)

    # trim three prime adapters and remove reads without
    if three_prime_adapter is not None:
        print('====> Removing 3\' adapters and filtering reads...', file=sys.stderr)
        cutadapt_three_prime = [
            cutadapt_path, '-a', three_prime_adapter, '-m', '1', '--discard-untrimmed', '-o', os.path.join(output, 'three_prime_removed.fq'),
            '-j', str(threads), current_input
        ]
        run(cutadapt_three_prime, capture_output=not verbose)

        current_input = os.path.join(output, 'three_prime_removed.fq')
        files_to_remove.append(current_input)

    shutil.copy2(current_input, os.path.join(output, 'trimmed_reads.fq'))
    if not keep_intermediate:
        for filename in files_to_remove:
            os.remove(filename)