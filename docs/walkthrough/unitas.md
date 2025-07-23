# Running Unitas

## Basic Configuration

Unitas requires two parameters to work - a directory containing length sorted small RNA FASTQ and a set of refseq FASTA files to use a reference for the different types of small RNA or features that they may come from.

The former (directory containing length sorted small RNA FASTQ) can be specified in two different ways -
1. Directly though the `size_sorted_fastqs` parameter e.g.

```yaml
size_sorted_fastqs: sorted_small_rna
unitas:
  ...
```

2. Through having both `smallRNA_fastq` and `sort` specified, in which case it uses the output of `sort`. e.g.

```yaml
smallRNA_fastq: trimmed_smallrna.fastq
sort:
  ...
unitas:
  ...
```

The latter (set of refseq FASTA files) are specified as a list and can be written in three different forms:

1. If the FASTA file is not already labelled for unitas in the IDs (i.e. ID values do not start with `[class_name]|`) the class name for the pipeline to add can be specified as a key and the file as a value. For example to add a file of miRNA called `microRNA.fasta` that should be labelled `miRNA` specify the following:

```yaml
unitas:
  refseq:
    - miRNA: microRNA.fasta
```

2. If the file already has unitas labels in the IDs (i.e. ID values start with `[class_name]|`) you only need to specify the path to the FASTA file. e.g. if you have a labeled set of tRNA in a file called `transferRNA.fasta` specify the following:

```yaml
unitas:
  refseq:
    - transferRNA.fasta
```

3. If you want to add gene features from `cds` and/or `unspliced_transcriptome` specify the special string `gene`. e.g.

```yaml
cds: cds.fasta
unspliced_transcriptome: unspliced.fasta
unitas:
  refseq:
    - gene
```

These can be mixed and matched, so a full example could be:

```yaml
smallRNA_fastq: trimmed_smallrna.fastq
sort:
  ...
unitas:
  refseq:
    - gene
    - transferRNA.fasta
    - miRNA: microRNA.fasta
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
### species

If unitas allows, you can set a species to have unitas download existing data on it's small RNAs from the internet to help with classification. For example:

```yaml
unitas:
  species: caenorabditis_elegans
```

## Utility Programs

Two utility programs to help with creating input for unitas are contained in the small RNA pipeline.

### label_for_unitas

This is a simple program that adds a class label that unitas can read to a FASTA file of sequences. It is recommened to not use this directly most of the time and instead use the key-value labelling syntax shown above. However, if you need a unitas labelled file, you can run the following:

```
$ label_for_unitas [class] [FASTA_file]
```

e.g. to label `transposons.fasta` with TE and write the result to `transposons_labelled.fasta`, you could run:

```
$ label_for_unitas TE transposons.fasta -o transposons_labelled.fasta
```

!!! tip
    Remember if you install via apptainer or docker you will need to prepend the command from the install page. For apptainer this command becomes:
    ```
    $ apptainer run -B $HOME:$HOME hlsmallrna.sif label_for_unitas TE transposons.fasta -o transposons_labelled.fasta
    ```


### extract_nc

If you want to test what small RNA come form noncoding reigons - outside of genes - specifically you can use this tool to extract those from the genome. It also adds the unitas label noncoding to the ID of the sequences for you. 

!!! warning
    Most of the time you do not want to do this on the first run - as Unitas will likely double classify anything that aligns to something like a transposon. It is advised to take anything unitas has not classified on the first run and only run that on the noncoding if you want the information.

It can be run as follows:

```
$ extract_nc [genome_FASTA] [annotation_GFF] -o [output_FASTA]
```

!!! tip
    Remember if you install via apptainer or docker you will need to prepend the command from the install page. For apptainer this command becomes:
    ```
    $ apptainer run -B $HOME:$HOME hlsmallrna.sif extract_nc [genome_FASTA] [annotation_GFF] -o [output_FASTA]
    ```

## Output Files
### unitas

Directory containing the raw output of unitas for each length file run. Can be inspected to extract details for particular small RNAs.

### unitas_summary.csv

Human-readable unitas classification summary by length. Conatins the classifications for each RNA length on a row, allowing for inspection of all the details at once.

### unitasGraph.png

Graph of unitas classification by length - uses only top level unitas classes.

### unitas_graph_data.csv

Data for `unitasGraph.png`, in case it needs to be replotted, misses the more granular classes from `unitas_summary.csv`.