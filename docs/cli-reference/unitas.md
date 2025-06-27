# Unitas

```
unitas [-h] [-d CDS] [-u UNSPLICED_TRANSCRIPTOME]
                     	[-r [REFSEQ [REFSEQ ...]]] [-s SPECIES]
                     	path_to_rnas

positional arguments:
  path_to_rnas      	Path to the folder with varying length RNAs in

optional arguments:
  -h, --help        	show this help message and exit
  -d CDS, --cds CDS 	Optional CDS region, passed to unitas
  -u UNSPLICED_TRANSCRIPTOME, --unspliced-transcriptome UNSPLICED_TRANSCRIPTOME
                    	Optional, unspliced transcriptome, passed to unitas
  -r [REFSEQ [REFSEQ ...]], --refseq [REFSEQ [REFSEQ ...]]
                    	References for use with unitas
  -s SPECIES, --species SPECIES
                    	Species to set in unitas arguments
```

Input files:

- `path_to_rnas` - path to the directory containing fastq files of each set of  small RNAs you want to procesus in unitas

*Note: to get unitas to work, you need to set either the species to a compatible species or provide one or more files of reference small RNAs. These can either be done with the -s and -r optional arguments or by adding them to your configuration file (see below).*

- `CDS` - Coding sequence, pipeline automatically labels with Gene and combines with the unspliced transcriptome, optional

- `unspliced_transcriptome` - Unspliced transcriptome, pipeline automatically labels with `Gene` and combines with the `CDS`, optional

Output files:

- `unitas/` - directory containing the results for all the unitas runs performed by the program

- `unitas_summery.csv` - CSV file containing the combined summaries of all of the unitas runs (plus some extra numbers). It is recommended to open it in a spreadsheet package (e.g. microsoft excel, libreoffice calc) as there is a lot of data and it may look a mess in a text editor.

- `UnitasGraph.png` - graph showing the categories allocated by unitas against the length of the small RNA

- `unitas_graph_data.csv` - raw data used for creating unitasGraph.png to allow for easy regraphing
