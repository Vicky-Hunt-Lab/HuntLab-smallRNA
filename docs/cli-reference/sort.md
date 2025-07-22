# Sort

```
sort [-h] [-d CDS] [-l MIN_LENGTH] [-x MAX_LENGTH]
                   	[-m REF_MISMATCHES] [--disable-alignment]
                   	small_rna [genome]

positional arguments:
  small_rna         	Path to FASTQ containing the small RNA
  genome            	Genome to align against

optional arguments:
  -h, --help        	show this help message and exit
  -d CDS, --cds CDS 	Optional CDS region, also align this to the CDS region
                    	as well as the genome
  -l MIN_LENGTH, --min-length MIN_LENGTH
                    	Minimum length to bin
  -x MAX_LENGTH, --max-length MAX_LENGTH
                    	Maximum length to bin
  -m REF_MISMATCHES, --ref-mismatches REF_MISMATCHES
                    	Number of mismatches to use in bowtie2, None for
                    	default behaviour
  --disable-alignment   Skip the alignment to the reference genome step
```

Input files:

- `small_rna` - fastq file containing RNA with adapters removed
	
- `genome` - Reference genome of the species you are using to align against in fasta format

Output files:

- `Bbmap_index/` - contains the index of the reference genome created by bowtie2

- `Mapped_sequences.fastq` - sequences successfully mapped to the reference by bowtie2

- `binned_rna/` - directory containing one fastq file for each length of sequence contained in the bowtie2 output

- `rna_length_report.csv` - table showing a summary of the RNAs by length and first base

- `Baseplot.png` - plot of the length and first base of the RNAs, using the data in rna_length_report.csv

- `baseplot_data.csv` - raw data used to make baseplot.png to allow for easy regraphing

*Note BBMap has been replaced with bowtie2 for this step, but file names haven't been changed.*