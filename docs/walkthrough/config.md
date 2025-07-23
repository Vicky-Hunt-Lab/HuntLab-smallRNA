# The Config File

HuntLab-smallRNA 2 uses an overhauled config system that is different to HuntLab-smallRNA 1. This is based on a YAML file and the idea that anything that affects the analysis should be in the config file, while anything that doesn't is a command line argument. Alongside this it has been made easier to run multiple stages in one command. An annotated version of a minimal config file to run all 4 parts of the analysis is as follows:

```yaml
# The small RNA FASTQ file to use as input to the pipeline
smallRNA_fastq: smallrna.fastq

# The trim section, as it is present trim will be run
trim:
  # Use the bulit-in adapaters for the Qiagen kit
  kit: qiagen

# The sort section, as it is present sort will be run
sort:
  # The genome to align to filter small RNA 
  genome: genome.fasta

# The unitas section, as it is present unitas will be run
unitas:
  # files to pass to unitas refseq and the labels to give them
  refseq:
    - miRNA: test/miRNA.fasta
    - piRNA: test/piRNA.fasta
    - tRNA: test/tRNA.fasta
    - TE: test/transposable_elements.fasta

# The targetid section, as it is present targetid will be run
targetid:
  # Files containing the interesting targets to align to
  target_files:
    - test/file1.fasta
    - test/file2.fasta
```

!!! note "Absolute and Reletive Paths"

    When parsing the config file, the small RNA pipeline assumes that all paths in the config file are reletive to the directory the config file is in. So if you move the config file, things will break. If you need to move the config file, use absolute paths (that start with `/`) instead.

To run this, you could save the file as `hlsmallrna_config.yml` and run:

```
$ hlsmallrna hlsmallrna_config.yml
```

!!! tip
    Remember if you install via apptainer or docker you will need to prepend the command from the install page. For apptainer this command becomes:
    ```
    $ apptainer run -B $HOME:$HOME hlsmallrna.sif hlsmallrna hlsmallrna_config.yml
    ```

The results would then appear in the directory `hlsmallrna_output/`, presuming all the sequence data referenced is real. A more comprehensive file to do a similar analysis to the first one is:

```yaml
# The small RNA FASTQ file to use as input to the pipeline
smallRNA_fastq: smallrna.fastq
# CDS FASTA of the species of interest, can be used in sort, unitas and targetid
cds: cds.fasta
# unspliced transcriptome FASTA of the species of interest, can be used in unitas and targetid
unspliced_transcriptome: unspliced.fasta


# The trim section, as it is present trim will be run
trim:
  # The 5’ adapter sequence to trim
  5_prime: ACGTTTAG
  # The 3’ adapter sequence to trim
  3_prime: CGTAGGAT
  # The quality filter cutoff to use
  min_quality: 20

# The sort section, as it is present sort will be run
sort:
  # The genome to align to filter small RNA 
  genome: genome.fasta
  # If True, also align to the CDS file above
  align_to_cds: True
  # Minimum length of small RNA sequence to keep
  min_length: 15
  # Maximum length of small RNA sequence to keep
  max_length: 50
  # Max number to allow when aligning to the genome and CDS
  mismatches: 0

# The unitas section, as it is present unitas will be run
unitas:
  # files to pass to unitas refseq and the labels to give them
  refseq:
    # add genes to unitas, special keyword to add using cds and unspliced_transcriptome specified earlier
    - gene
    - miRNA: test/miRNA.fasta
    - piRNA: test/piRNA.fasta
    - tRNA: test/tRNA.fasta
    - TE: test/transposable_elements.fasta
  # species name to pass to unitas
  species: x

# The targetid section, as it is present targetid will be run
targetid:
  # Minimum length of small RNA, higher number will speed up runtime, used to calculate bowtie2 seed length
  min_seq_length: 5
  # Number of mismatches to allow when aligning to the targets
  mismatches: 0
  # Files containing the interesting targets to align to
  target_files:
    - test/file1.fasta
    - test/file2.fasta
  # If present also enrich GO terms, KEGG pathways and Pfams of targets
  enrich:
    # Path to the eggnog-mapper data dir
    eggnog_data_dir: /home/user/eggnog-mapper-data
    # Target files to ignore during enrichment
    exclude_files:
      - test/file2.fasta
```

If you want to skip a step, just leave it out of the file. e.g. if you only want to run sort and unitas, you could use the file:

```yaml
smallRNA_fastq: smallrna.fastq
sort:
  genome: genome.fasta
unitas:
  refseq:
    - miRNA: test/miRNA.fasta
    - piRNA: test/piRNA.fasta
    - tRNA: test/tRNA.fasta
    - TE: test/transposable_elements.fasta
```

!!! tip "Command Line Arguments"
    While most things are now specified in the config file, there are still a few CLI arguments that are as follows:

    `-o`, `--ouput` - directory to write pipeline output to (default: hlsmallrna_output)

    `-t`, `--threads` - number of threads to use when running analysis tools (default: 4)

    `-k`, `--keep-files` - If set, don't delete intermediate files, useful for debugging or if you need an intermediate file later

    `-v`, `--verbose` - makes the pipeline print out the output of the intermediate commands

The following sections will explain each parameter in the config file in detail.
