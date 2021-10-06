# Small RNA Pipeline for the Hunt Lab

This pipeline is designed to automate commen tasks to do with small RNA analysis in the lab. It currently enables four different tasks:

1. "Process" - Takes the small RNA, removes specified adapter sequences (with cutadapt) and performs a quality check with FastQC. Automatically trims the 3' end to remove bases that fail to pass a minimun threashold.

2. "Sort" - Aligns the small RNA against a specified genome, removes ones that fail to align and splits them into files baised on their length

3. "Classify" - Uses unitas to try and classify the small RNA into groups and also runs a differntial expression analysis with edgeR to allow visulising which RNA are upregulated.

4. "TargetID" 

## Configeration File

This pipeline can be configered with a "config.toml" file. This file is written in ["Tom's Obvious, Minimal Language"](https://github.com/toml-lang/toml) (TOML) format, which is a common format for configeration. This file is split into multiple sections, with headers in square brackets. Valid keys are shown in the example file below. (Note lines starting with a # are comments like in python)

```toml
# Configeration that doesn't fit anywhere else
[general]
# path to the directory to put the output in
output_directory = "/path/to/output"

# Number of threads to pass to any program that allows for multiprocessing
threads = 8

# Section containing configeration for each of the exteranl programs called while in use
[cli-tools]

[cli-tools.trim]
trim_pass_threads = true
path_to_trim = "/path/to/cutadapt"
trim_params = []

[cli-tools.fastqc]
fastqc_pass_threads = true
path_to_fastqc = "/path/to/fastqc"
fastqc_params = []

[cli-tools.bowtie2]
path_to_bowtie2 = "/path/to/bowtie2"
path_to_bowtie2_build = "/path/to/bowtie2-build"
bowtie2_params = []
bowtie2_build_params = []

[cli-tools.samtools]
path_to_samtools = "/path/to/samtools"
samtools_view_params = []
samtools_sort_params = []

[cli-tools.bedtools]
path_to_bedtools = "/path/to/bedtools"
bedtools_bamToFastq_params = []

[cli-tools.unitas]
path_to_unitas = "/path/to/unitas"
unitas_params = []
```