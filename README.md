# Small RNA Pipeline for the Hunt Lab

This pipeline is designed to automate common tasks to do with small RNA analysis in the lab. It currently enables six different tasks:

1. "Process" - Takes the small RNA, removes specified adapter sequences (with cutadapt) and performs a quality check with FastQC. Automatically trims the 3' end to remove bases that fail to pass a minimum threshold.

2. "Sort" - Aligns the small RNA against a specified genome, removes ones that fail to align and splits them into files based on their length

3. "ExtractNC" - Extract the non-coding region to use with unitas. Needs a FASTA file containing a genome and a GFF containing feature labels

4. "Classify" - Uses unitas to try and classify the small RNA into groups and also runs a differential expression analysis with edgeR to allow visulising which RNA are up-regulated.

5. "TargetID" - Align small RNAs against a number of potential targets and creates a list of targets for each small RNA

6. "All" - run steps 1, 2 and 4 one after the other - piping the output of one into the next

## Configuration File

This pipeline can be configured with a "config.toml" file. This file is written in ["Tom's Obvious, Minimal Language"](https://toml.io/en/) (TOML) format, which is a common format for configeration. This file is split into multiple sections, with headers in square brackets. Valid keys are shown in the example file below. (Note lines starting with a # are comments like in python)

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

[cli-tools.bbmap]
bbmap_pass_threads = true
path_to_bbmap = "bbmap.sh"
bbmap_align_params = []
bbmap_index_params = []

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