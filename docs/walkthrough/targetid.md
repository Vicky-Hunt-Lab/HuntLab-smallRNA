# Running TargetID

## Basic Configuration

To run TargetID, you need two things:

1. A small RNA file - if sort is specified the output of that is used - if not `smallRNA_fastq` is used.

2. One or more files containing potential targets of small RNAs (e.g. genes, transposons, lncRNA etc.) specified in the `targetid` section as `target_files`.

For example:

```yaml
smallRNA_fastq: small_rna.fastq
targetid:
  target_files:
    - targets/lncRNA.fasta
    - targets/genes.fasta
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

### min_seq_length

This is used to set the seed length in bowtie2 when aligning to targets (defaults to 5). Larger numbers makes it fasta but less accurate for sequences close to this length. For example, to set to 18 add the following:

```yaml
targetid:
  min_seq_length: 18
```

### mismatches

Number of mismatches to allow in bowtie2 when aligning to the target files, defaults to 0. See the [Considerations](prereq.md) page for information about why you'd want to change this. The example below sets this to allow 1 mismatch:

```yaml
targetid:
  mismatches: 1
```

### enrich - eggnog_data_dir

The pipeline can optionally enrich GO terms, KEGG pathways and PFams for the targets of the small RNA. This is a rough analysis based on `emapper.py --itype CDS` mode that uses `diamond blastx` to map, so is likely only a rough analysis. May be useful if you want to have a quick look at if there are any common annotations within the targets though should be redone probably later. To run this you need to add the eggnog_data_dir to the config file. This directory can be created with the command:

```
$ mkdir -p /home/user/path/to/eggnog-data
$ download_eggnog_data.py -y --data_dir /home/user/path/to/eggnog-data
```

!!! tip
    Remember if you install via apptainer or docker you will need to prepend the command from the install page. For apptainer this command becomes:
    ```
    $ mkdir -p /home/user/path/to/eggnog-data
    $ apptainer run -B $HOME:$HOME hlsmallrna.sif download_eggnog_data.py -y --data_dir /home/user/path/to/eggnog-data
    ```
Then add the following to the config file

```yaml
targetid:
  ...
  enrich:
    eggnog_data_dir: /home/user/path/to/eggnog-data
```

### enrich - exclude

If you have run targetid with targets it doesn't make sense to map to the eggnog database (e.g. lncRNA), you can exclude them with the exclude key, but must ensure the path exactly matches the one in the target_files list. e.g.

```yaml
targetid:
  target_files:
    - targets/genes.fasta
    - targets/lncRNA.fasta
  enrich:
    eggnog_data_dir: /home/user/path/to/eggnog-data
    exclude:
      - targets/lncRNA.fasta
```

## Output Files
### rna_target_list.tsv

List of small RNA and there targets identified from the files specified.

### targets_by_group

Same information as `rna_target_list.tsv`, but split by small RNA group (length and first base).

## Output Files - Enrich

These files are only produced if the `enrich` key is present.

### eggnog_mapper

Results of the eggnog-mapper run on the target files.

### targets_enriched_go_terms.tsv

Enriched GO terms in the targets of small RNA vs the rest of the target files.

### targets_enriched_kegg_pathways.tsv

Enriched KEGG pathways in the targets of small RNA vs the rest of the target files.

### targets_enriched_pfams.tsv

Enriched PFams in the targets of small RNA vs the rest of the target files.