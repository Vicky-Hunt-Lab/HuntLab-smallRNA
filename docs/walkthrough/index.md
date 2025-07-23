# Walkthrough

This walkthrough covers all of the major tasks of the pipeline and how they can be performed, either together or separately and considerations for your analysis. To start, here is an overview of the analysis the pipeline can perform:

## Trim

This filters low quality ends, trims 5' adapters and then trims 3' adapters. Any reads without a detectable 3' adapter are removed, as they may not cover a full small RNA.

## Sort

Firstly, aligns the small RNA to the genome and removes any that don’t align to remove contamination. Then sorts into FASTQ files of each length within a specified range. Finally, plots a graph of frequency of each length and first base of small RNA to give an overview of your sample.

## Unitas

This step takes the files of each small RNA length from sort and runs unitas on them to classify the origins of the small RNA. Then plots a graph of the origins of each length of small RNA.

## TargetID

This aligns the reverse complement of the small RNA exactly to FASTAs of interesting targets (e.g. genes in the organism). These confident targets can then optionally be enriched for GO terms, KEGG pathways and Pfams using eggnog-mapper to allocate.

!!! note "Pipeline Stage Ordering"
    While you do not have to run all of these parts of the pipeline, they will always run in the order:

    Trim -> Sort -> Unitas -> TargetID -> Enrich 

    If you want to do a different order you will need to write multiple config files and run the pipeline multiple times.
