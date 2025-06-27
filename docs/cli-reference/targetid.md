# TargetID

```
targetid [-h] [-m MIN_SEQ_LENGTH] [-t TARGET_FILES [TARGET_FILES ...]]
                       	[--num-mismatches NUM_MISMATCHES]
                       	small_rna

positional arguments:
  small_rna         	Path to the FASTQ containing the small RNA to find targets of

optional arguments:
  -h, --help        	show this help message and exit
  -m MIN_SEQ_LENGTH, --min-seq-length MIN_SEQ_LENGTH
                    	Minimum sequence length to properly align
  -t TARGET_FILES [TARGET_FILES ...], --target-files TARGET_FILES [TARGET_FILES ...]
                    	Files containing genome features that could be targeted
  --num-mismatches NUM_MISMATCHES
                    	Number of mismatches to allow in the alignment, defaults to 0
```

Input: 

- `small_rna` - fastq file containing small RNA to look for targets of

- `target_files` - one or more fastq files containing a list of potential target sequences

Output:

- `bowtie_indexes/` - indexes produced of the target files with bowtie2

- `target_alignments/` - SAM files containing the results of alignment attempts and FASTQ files containing all the small RNA successfully aligned against the targets. One set of files are produced for each target file provided

- `rna_target_list.csv` - list of RNA target pairs

*Note that the samtools view command used in this step is the only one that has defaults pre-set in the `samtools_view_params` config key, namly `-h -F 256 -F 4`. They will be removed if you set the value differently in the config file. To only add parameters, you will need to set these ones again before the extra ones you want to set. E.g. `samtools_view_params = ["-h", "-F", "256", "-F", "4", ...]`*

*Make sure all of your target files are labelled with what they are, by prepending the name and | (the bar or pipe character) to the ID in the fasta file. If you don’t do this, it will not be nicely sorted in the results correctly.*
