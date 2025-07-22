# Prerequisites and Considerations

Before we get on to how to run the pipeline this page covers what files are needed to run the pipeline and things to think about to inform the parameters to set. Read each section you want to run carefully.

## Trim
The trimming in the pipeline is what we think is best practice for the small RNA data we use - remove low quality ends and trim any 5' adapters first, then trim the 3' adapter and remove any sequences that lack it - as they are likely incomplete or noise. If you want to do something different, you can do this and provide the already trimmed small RNAs to the pipeline to run the other analysis. To run this, you will need the following data:

-	Small RNA-seq reads in FASTQ or FASTA format
-	3' adapter sequence from small RNA library prep
-	Optionally, a 5' adapter sequence and quality cutoff to use

!!! tip "Small RNA Library Prep Kits"
    For the Quiagen and NEB kits the best adapters to use are built into the pipeline. You just need to know the kit and if we have data in the forward or reverse direction. 

!!! note "Paired End Sequencing"
    Most small RNA-seq provides a single sequence, but if you have paired end data, you should only provide the read that contains the sequence of the small RNA to the program (Usually the forward read), not it's reverse complement.

## Sort

This takes trimmed reads and aligns them to the genome and optionally a set of CDS sequences to catch small RNA derived from spliced genes. Then sorts them into length and plots length and first base to give you an overview of the small RNAs you have. It also outputs a summary table on the number and percentage of small RNA reads aligned and FASTQs of small RNAs of each class (length and first base) for easier downstream analysis. For this you'll need:

-	A FASTQ or FASTA of trimmed small RNA sequences (can be automatically retrieved if trim is run)
-	A genome assembly to align the reads to
-	Optionally, a set of CDS sequences to also align reads to 
-	Optionally, A minimum and maximum length of small RNA to consider, otherwise defaults to 18 - 30 nt
-	Optionally, A number of mismatches to use when aligning, otherwise defaults to 0

!!! tip "Mismatches"
    By default, sort allows no mismatches or gaps, which is fine if you have a good reference genome that is of a close relative to the sample small RNA-seq was performed on and want to be conservative about what you retain. If you know your genome is likely to have lots of SNPs from the reference or are uncertain about how good you assembly is, you may want to up the amount of mismatches allowed to compensate for this.

!!! warning "Downstream use of BLAST"
    While it seems like a good to run small RNAs through BLAST downstream, if you do this you need to be careful of the parameters used, particularly word size and E value. This is because BLAST will not return any hits if your sequence is shorter than the word size set (often above 20nt) and the short size of small RNA may mean even close to exact marches have a high E value. Due to this a basic search may return no hits, even if huntlab-smallrna has found an exact match in your genome.

## Unitas

Unitas runs unitas on a number of files representing small RNA lengths and plots a graphs of the origin or class for each small RNA length. This requires the following files:

-	A directory of small RNAs of different lengths named in the form length.fastq - where X is the length of the small RNA. This can be automatically derived if sort has been run
-	One or both of:
    - A species name in the form supported by unitas - see examples in [https://www.smallrnagroup.uni-mainz.de/software/unitas_documentation_1.7.0.pdf](https://www.smallrnagroup.uni-mainz.de/software/unitas_documentation_1.7.0.pdf) 
    - One or more files containing known sequences of a particular class (e.g. miRNA, rRNA, tRNA, Transposable Elements etc.)
-	Optionally a CDS and/or unspliced transcriptome files to supplement the above

!!! tip "Labelling for Unitas"
    Unitas expects files that have the class of each sequence labelled in their ID e.g. for TEs with ids te1 and te2, the file would look like:
    ```
    >TE|te1
    ACGT...
    >TE|te2
    CGAT...
    ```
    Unless you have mixed files, the pipeline can label this for you, so don't worry about adding that yet. See the unitas section for information on this.

## TargetID
This stage aligns the reverse complement exactly to lists of potential targets, to find likely targets of a small RNA, outputting both and overall list and by class (e.g. 22Us). Then it can optionally enrich the targets GO terms, KEGG pathways and Pfams to help understand the tragets. For this you will need:

-	FASTQ file of small RNAs or a directory or sorted small RNAs - can be derived from trim or sort if run
-	One or more FASTA files containing sequences of potential targets (e.g. genes, TEs, host genes etc.)
-	Optionally, the number of mismatches to tolerate when aligning, if not 0
-	If you want to run the enrichment, an eggnog-mapper data directiory (see targetid section on how to create this)
-	Optionally, a list of target files to exclude from the enrichment

Now you know what you need you can continue to the next page of the tutorial to understand the config file you need to write.