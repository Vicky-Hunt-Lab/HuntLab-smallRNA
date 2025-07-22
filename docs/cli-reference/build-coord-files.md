# Build Coord Files

```
usage: build_coord_files [-h] [-a] [-q] [-s] [-b] [-g] [-r] [-c]
                     	[--feature FEATURE] [--genome-coords GENOME_COORDS]
                     	[--band-width BAND_WIDTH] [-o OUTPUT]
                     	input_file

Convert common bioinformatics file formats into coordinate file to plot

positional arguments:
  input_file        	File to convert, attempts to autodetect type

optional arguments:
  -h, --help        	show this help message and exit
  -a, --fasta       	Treat input as a FASTA file
  -q, --fastq       	Treat input as a FASTQ file
  -s, --sam         	Treat input as a SAM file
  -b, --bam         	Treat input as a BAM file
  -g, --gff         	Treat input as a GFF file
  -r, --rm-fa-out   	Treat input as a RepeatMasker .fa.out file
  -c, --scaffold-aware  Merge scaffolds into one chromosome
  --feature FEATURE 	Select a gff feature to use
  --genome-coords GENOME_COORDS
                    	File containing the coordinates of the genome
  --band-width BAND_WIDTH
                    	Size of the bands showing change in scaffolds
  -o OUTPUT, --output OUTPUT
                    	File to output to
```