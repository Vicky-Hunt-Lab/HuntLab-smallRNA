# Process

```
process [-h] [-a ADAPTER] [-g FRONT] [-b ANYWHERE]
                       	[-c CUTOFF]
                       	small_rna

positional arguments:
  small_rna         	Path to FASTQ containing the small RNA

optional arguments:
  -h, --help        	show this help message and exit
  -a ADAPTER, --adapter ADAPTER
                    	Sequence of the adapter to remove from the 3' end
  -g FRONT, --front FRONT
                    	Sequence of the adapter to remove from the 5' end
  -b ANYWHERE, --anywhere ANYWHERE
                    	Sequence of the adapters to remove from both ends
  -c CUTOFF, --cutoff CUTOFF
                    	Quality cutoff to trim RNA sequences at
```
Input files:

- `small_rna` - fastq file containing raw RNA-seq data

Output files:

- `Fastqc/` - directory containing the raw output of FastQC

- `Trimmed_rna.fastq` - file containing the raw RNA with adapters removed (only produced if an adapter sequence is provided)

- `Cut_sequences.fastq` - sequences with low quality parts removed

*Note that the parts of this step can easily be skipped. If you want to skip adapter trimming, don’t specify any adapters with `-a`, `-g` or `-c`. If you don’t want to run fastQC set `-c 0`. The program knows that in these cases, the steps do not have to be run.*
