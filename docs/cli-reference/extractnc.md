# ExtractNC

```
usage: extract_nc [-h] [-o OUTPUT] genome gff_file

positional arguments:
  genome                FASTA containing the genome to extract from
  gff_file              GFF file containing annotations of CDS and mRNA
                        regions

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        FASTA file to write output to
```

Input files:

- `genome` - file containing the whole genome of the species of interest

- `gff_file` - GFF3 file containing annotations at least for the mRNA and coding region (labelled CDS) for the genome

Output files:

- `noncoding.fasta` - FASTA file containing only the transcribed noncoding regions of the DNA

!!! warning
    The method used here will not work if there a coding region is labelled outside a single mRNA region. The program does a check before running to make sure this is true and throws an exception before trying to run if so, to prevent nonsense results from being produced.
