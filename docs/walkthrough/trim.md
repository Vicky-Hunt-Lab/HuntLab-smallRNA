# Running Trim

## Basic Configuration

Trim requires a file of untrimmed small RNAs and the adpaters to trim - either explicitly or through the kit used - to be specified. Examples of config files for  possible usage is as follows.

Trimming adapters with the 5' sequence ACGTTTAGAT and the 3' sequence GCTAGGATA:

```yaml
smallRNA_fastq: smallrna_untrimmed.fastq
trim:
  5_prime: ACGTTTAGAT
  3_prime: GCTAGGATA
```

Trimming adapters from the qiagen kit in forward direction:

```yaml
smallRNA_fastq: smallrna_untrimmed.fastq
trim:
  kit: qiagen
```

Trimming adapters from the qiagen kit in reverse direction:

```yaml
smallRNA_fastq: smallrna_untrimmed.fastq
trim:
  kit: qiagen-rev
```

Currently the valid values for kit are:

- `qiagen` - qiagen kit, forward direction
- `qiagen-rev` - qiagen kit, reverse direction
- `neb` - NEB kit, forward direction 
- `neb-rev` - NEB kit, reverse direction

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
### min_quality
This sets the quality cutoff in cutadapt when filtering low quality ends. Defaults to 20 - example below sets it to 15:

```yaml
trim:
  ...
  min_quality: 15
```

## Output Files
### trimmed_reads.fq

This is a FASTQ file containing all the trimmed reads that passed the filters of quaility and had a 3' adapter present.