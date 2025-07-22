# Running Sort

## Basic Configuration

At minimum, sort requires a set of trimmed reads and a genome to align to from the species the small RNAs came from. If trim is in the same config file, the set to trimmed reads will be taken from the output of trim. Examples are as follows.

Without running trim.
```yaml
smallRNA_fastq: smallrna_trimmed.fastq
sort:
  genome: genome.fasta
```

With running trim.
```yaml
smallRNA_fastq: smallrna_untrimmed.fastq
trim:
  ...
sort:
  genome: genome.fasta
```

Once this config is written, it can be run as usual using:

```
$ hlsmallrna config.yml
```

!!! tip
    Remember if you install via apptainer or docker you will need to prepend the command from the install page. For apptainer this command becomes:
    ```
    $ apptainer run -B $HOME:$HOME hlsmallrna.sif hlsmallrna config.yml
    ```

## Optional Parameters

### align_to_cds

If true aligns reads that don't align to the genome to the specified CDS file, keeping any that align. This allows for small RNA that are fragments of genes and align accross introns to be kept as they may not align to the genome. 

```yaml
cds: cds_sequences.fasta
...
sort:
  ...
  align_to_cds: True
```

### min_length

Sets a minimum sequence length to look at (default is 18 nt). Anything shorter is filtered during sorting. e.g. to filter anything below 10 nt set the following:

```yaml
sort:
  min_length: 10
```

### max_length

Sets a maximum sequence length to look at (default is 30 nt). Anything longer is filtered during sorting. e.g. to filter anything above 50 nt set the following:

```yaml
sort:
  max_length: 50
```

### mismatches

Number of mismatches to allow in bowtie2 when aligning to the genome and CDS, defaults to 0. See the [Considerations](prereq.md) page for information about why you'd want to change this. The example below sets this to allow 1 mismatch:

```yaml
sort:
  mismatches: 1
```

!!! warning
    The small RNA pipeline supresses gaps by makeing them equal to 100 mismatches. If you set this over 100, some gaps may be included in the alignment.

## Output Files

### alignment_report.tsv

This is a table showing an overview of how many small RNA reads successfully aligned to the genome and optionally the CDS. It contains the following metrics:

- Total Reads
- Total Mapped
- Percentage Mapped
- Total Unmapped
- Percentage Unmapped
- Mapped to Genome
- Percentage Mapped to Genome
- Unmapped to Genome
- Percentage Unmapped to Genome
- Mapped to CDS
- Percentage Mapped to CDS (but not to genome)
- Unmapped to CDS or genome
- Percentage Unmapped to CDS or genome

### counts.tab

Read counts is TSV format for each small RNA sequence, designed as input to a differential expression analysis, if that wanted.

### rna_length_report.csv

Table of counts of small RNA for each length and first base, in a easily human readable format.

### baseplot.png

Plot of frequency of length and first base in the dataset.

### baseplot_data.csv

That length and first base data used in baseplot.png - normalised to percentages. Should be used if the user wants to replot that graph or feed into downstream analysis.

### binned_length_rna

Directory containing FASTQ files for each RNA length - named lengthX.fastq where X is the length of small RNA in the file.

### binned_group_rna

Directory containing FASTQ files for each RNA group - defined by length and first base - named {length}{first_base}.fastq (e.g. 22G.fastq) where {length} is the length of small RNA in the file and {first_base} is the letter corresponding to the first base.