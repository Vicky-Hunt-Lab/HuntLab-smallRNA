# Quick Start

To use this pipeline, you will need at least a FASTQ file containing raw small RNA-Seq reads. Then depending on what you want to achieve, you can run different steps of the pipeline. If you ever get stuck and need to see the options of a particular stage, all parts have their own help message, that can be seen by passing `--help` after the name of the stage. 

Before starting, you will need to create a configuration file that allows you to specify some parameters in advance, so you don’t have to keep retyping them. By default, it looks for these in a file in your current directory called `config.toml`, alternatives can be specified by passing the name to the `-C` argument before the name of the stage to run. More details on this can be found in the *Configuration file* section, but for now setup a simple file to specify an output directory of `output/` by putting the following in `config.toml`:

```toml
[general]
output_directory = "./output"
```

In addition, check that all of the commands can be run by just typing their name in the command line, as this is how the program tries to use them by default. If they can’t, you can set a custom path by setting the relevant `path_to_` variable in the configuration file. For example if unitas is run with `unitas.pl` instead of just `unitas` add the following to the bottom of the config file:

```toml
[cli-tools]
[cli-tools.unitas]
path_to_unitas = "unitas.pl"
```

## Process

The first stage is Process. This is for if you have raw data that hasn’t had it’s adapters trimmed and been quality filtered, so if your data has had this done to it already, this step can be skipped. To run you need a FASTQ file containing raw small RNA and the sequence and end of any adapters attached during sequencing.

Once you have this data, to include the trimming in the run, you need to set one of three flags with the adapter sequence, depending on the position of the particular adapter. 

If your adapter is on the 3’ end only, use `-a`; if it is from the  5’ end only, use `-g` and if it could be on either end use `-b`. If you don’t need adapters trimmed, just don’t specify any of them and the step will be skipped. For example, if your small RNAs are in `smallRNA.fastq` and you want to trim the adapter sequence `AGCATA` from the 3’ end only, you would run:

```
$ hlsmallrna process -a AGCATA smallRNA.fastq
```

In addition, this step runs FastQC to produce a report on the quality of the small RNA and automatically cuts anything with a lower quartile of quality below a cutoff. By default, this is set to 20, but can be changed by the user using the `-c` flag. If you do not want this, you can set `-c 0` and the FastQC step will be skipped for efficiency. For example if you want to change the cutoff for the last command to 5, you could run:

```
$ hlsmallrna process -c 5 -a AGCATA smallRNA.fastq
```

This produces a number of files in the output directory, importantly `cut_sequences.fastq` has the sequences that have adapters trimmed and low quality parts removed. (See *CLI Reference* section for information on other output for all of the commands).

## Sort

The second stage is Sort. This stage aligns the small RNAs to the genome of the organism and filters out any that don’t align, to remove contamination. Then it splits the remaining reads into FASTQ files by length, to allow for easier classification. 

To run this, you need a FASTQ file of processed small RNA and a FASTA file containing the genome of the organism of interest. For example, if your small RNA are in `smallRNA.fastq` and your genome is in `genome.fasta` you can run:

```
$ hlsmallrna sort smallRNA.fastq genome.fasta
```

If you are only interested in a subset of lengths of small RNA, this step can be set to restrict the lengths it outputs with `-l` for minimum and `-x` for maximum. For example, if you only want small RNA with lengths between 12 and 30, you could run:

```
$ hlsmallrna sort -l 12 -x 30 smallRNA.fastq genome.fasta
```

By default this aligns the small RNA to the genome with bowtie2’s default settings, which allow different amounts of mismatches depending on how long the RNA is. Details of this can be found in the score configuration options in the manual, that can be found at: [https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). To deal with cases where we may want to do a more stringent run, the `-m` or `--ref-mismatches` argument can be supplied with the number of mismatches for bowtie2 to allow. For example, if you want to perform an operation using the same files as above, but only include small RNA that map exactly to the genome, you could run:

```
$ hlsmallrna sort -m 0 smallRNA.fastq genome.fasta
```

This command produces a few useful files in the output directory: 

-	`binned_rna/` - directory containing one fastq file for each length of small RNA

-	`rna_length_report.csv` - table showing a summary of the RNAs by length and first base

-	`baseplot.png` - plot of the length and first base of the RNAs (this is sometimes wrong, replot from the csv values for now, will fix matplotlib code soon)

## Unitas

The third stage is Unitas. This runs unitas on a directory containing a range of lengths of small RNA to classify the origin of the RNA. To do this, you will need files containing reference sets of the genome targets you wish to find, with labels in their IDs and a directory containing small RNA of different lengths, similar to the output of sort.

You can set the reference file using the `-r` command line argument, but it is easier to set them in the configuration file, as it makes the command shorter to type. If your reference files are in the same directory you are running in and called `ref1.fasta`, `ref2.fasta` and `ref3.fasta`, you can add the following to the config file:

```toml
[command]
refseq = [
    "ref1.fasta",
    "ref2.fasta",
    "ref3.fasta"
]
```

Note that all of these reference sequences should be labelled with the name of the category they belong to and the | (bar or pipe character) before the start of their ID. If this hasn’t been done and all of the sequences in a FASTA file should be labelled the same, you can use the `label_for_unitas` script that was installed with the small RNA pipeline. For example if you want to label the sequences in `unlabelledTEs.fasta` with the label TE, with the result produced in `teRef.fasta`, you could run:

```
$ label_for_unitas "TE" unlabelledTEs.fasta -o teRef.fasta
```

Once the correct reference file are set, if your small RNA, split by length, are in `output/binned_rna`, you can run:

```
$  hlsmallrna unitas output/binned_rna
```

This produces a useful summary and graph in the output directory:

-	`unitasGraph.png` - graph showing the categories allocated by unitas against the length of the small RNA (this is sometimes wrong, replot from the csv values for now, will fix matplotlib code soon)

-	`unitas_graph_data.csv` - summary of the raw numbers produced from the unitas runs, used to produce the graph

If you want to get small RNA that are derived from genes, you will need to obtain both of the mRNA and CDS regions and pass them to speciled arguments in the pipeline so they can be combined. This can either be done using the config file as followers or with the corresponding command line arguments (`--cds` and `--unspliced-transcriptome`):

```toml
[command]
refseq = [
    "ref1.fasta",
    "ref2.fasta",
    "ref3.fasta"
]
cds = "cds_sequences.fasta"
unspliced_transcriptome = "unspliced_transcriptome.fasta"
```

## Extract Non-Coding

Sometimes you want to see if a small RNA comes from a non-coding region of the transcriptome, so the pipeline contains a utility command to automatically extract this using a genome and some annotations. At minimum, these annotations need to contain a transcribed region labelled with `mRNA` and a coding region within the transcribed region labelled `CDS`. Most of the ones on [Wormbase Parasite](https://parasite.wormbase.org/index.html) are fine for this purpose.

As an example, to extract the regions from `genome.fasta` using the annotations in `annotations.gff` you could run:

```
$ hlsmallrna extractnc genome.fasta annotations.gff
```

This will produce a FASTA file containing the non-coding regions, complete with unitas labels applied, called `noncoding.fasta` in the output directory.

## Target Identification

The fourth and final main step is Target Identification. This takes a set of small RNA and a set of potential targets and produces a list of genome features that are targeted by each small RNA. 

As with unitas, you can set target files using the command line option `-t` but it is easier to add them to the configuration file. For example if your target files are `target1.fasta`, `target2.fasta` and `target3.fasta`, they can be set by adding the following in the section below `[command]`:

```toml
target_files = [
    "target1.fasta",
    "target2.fasta",
    "target3.fasta"
]
```

Then, if your small RNA are in `smallRNA.fastq` you can run the following:

```
$ hlsmallrna targetid smallRNA.fastq
```

Make sure the sequences in both `smallRNA.fastq` and each of the reference files have unique sequence ids, otherwise the command will fail.

If you are working with very short sequences (<5 bases), you will need to reduce the minimum length with the -m option. Conversely, if you know your minimum length is longer than this, this step can be sped up by increasing it. For example, if you know the shortest sequence you are dealing with is 10 bases long, you could run:

```
$ hlsmallrna targetid -m 10 smallRNA.fastq
```

Though this is not necessary and the program will produce a correct result for sequences longer than 5 without it, so only add if you need the extra speed.

This produces `rna_target_list.csv` which contains a list of small RNA target pairs and some general information about each pair’s alignment.

If you are dealing with a small RNA that doesn’t need perfect complementarity, you can add the `--num-mismatches` option to set the number of mismatches the alignment should allow. For example, if you want to allow up to four mismatches, you could run:

```
$ hlsmallrna targetid --num-mismatches 4 smallRNA.fastq
```

## All
Since the first three steps are often run one after the other, a shorthand command `all` is provided for convenience. This takes all of the arguments from the Process, Sort and Unitas steps, then runs them, piping the result from one into the next. All files produced by those steps are then produced in the output directory. 

For example if you have small RNAs in `smallRNA.fastq`, a genome in `genome.fasta` and refseq set in the config file and you want to trim the adapter sequence `ACATA` from both ends, but not run quality filtering you could run:

```
$ hlsmallrna all -c 0 -b ACATA smallRNA.fastq genome.fasta
```