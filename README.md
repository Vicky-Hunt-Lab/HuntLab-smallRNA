# Small RNA Pipeline for the Hunt Lab

This repository contains the code for the pipeline used to carry out sRNA-seq ananlysis in Vicky Hunt's lab. For usage and the tasks the pipeline can run, see the documentation, avalible in PDF, DOCX and HTML format in the `docs/` directory in the repository. 

Note this software is very new, please report any bugs found to the GitHub bug tracker. Also note that it currently does not support GZipped input files, though support will be considered in the future if it is wanted enough.

## Configuration File

This pipeline can be configured with a "config.toml" file. This file is written in ["Tom's Obvious, Minimal Language"](https://toml.io/en/) (TOML) format, which is a common format for configeration. This file is split into multiple sections, with headers in square brackets. Valid keys are shown in the example file below. (Note lines starting with a # are comments like in python). An example is below, see the documentation for more comprihensive instructions.

```toml
# Configeration that doesn't fit anywhere else
[general]
# path to the directory to put the output in
output_directory = "/path/to/output"

# Number of threads to pass to any program that allows for multiprocessing
threads = 8

# This section contains subsitues for optional command line arguments
[command]
# for process command
adapter = "ACGTTTTTT"
front = "TTTTTACG"
anywhere = ""
cutoff = 25

# for sort command
min_length = 7
max_length = 35

# for unitas command
refseq = [
    "/path/to/refseq1",
    "/path/to/refseq2"
]
species = "x"

# for targetid command
target_files = [
    "/path/to/targetseq1",
    "/path/to/targetseq2"
]
min_seq_length = 9
num_mismatches = 0

# Section containing configeration for each of the exteranl programs called while in use
[cli-tools]
[cli-tools.trim]
trim_pass_threads = true
path_to_trim = "cutadapt"
trim_type = ""
trim_params = []

[cli-tools.fastqc]
fastqc_pass_threads = true
path_to_fastqc = "fastqc"
fastqc_params = []

[cli-tools.unitas]
unitas_pass_threads = true
path_to_unitas = "unitas"
unitas_params = []

[cli-tools.bowtie2]
bowtie2_pass_threads = true
path_to_bowtie2 = "bowtie2"
bowtie2_params = []
path_to_bowtie2_build = "bowtie2-build"
bowtie2_build_params = []

[cli-tools.samtools]
path_to_samtools = "samtools"
samtools_view_params = ["-h", "-F", "256", "-F", "4"]
samtools_fastq_params = []
```