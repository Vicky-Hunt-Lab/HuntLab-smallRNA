# Configuration and Troubleshooting Notes

## Configuration File

For ease of use this program requires a configuration file in TOML format to run successfully. At minimum it can be a blank file, but it is strongly recommended setting at least the `output_directory` setting to be named something more useful than `output`, the default setting.

By default the program tries to use config.toml in your current directory, but it can be set to a different file using the `-C` option before specifying which task to run, e.g.:

```
$ hlsmallrna -C alternate.toml sort ...
```

The layout of a TOML file has a number of sections each with a number of configuration keys. Each key has an expected type (e.g. string, integer, floating point number etc.) and the program will check they are as expected before running. If one is wrong, the program will crash and alert the user of the problem. A list of these keys and example values can be found in the README.md included with the program.

Manually setting parameters of internal commands
For advanced uses, you may want to add parameters to the internal commands run. This can be done through variables in the cli-tools section of the configuration file. A full list of these parameters can be found in the README.md provided with the code. 

There is a main one that you can set, called `<name>_params`. It is a list of custom parameters to be added on to the command when called. For example if you wanted to set min-length to 10 for FastQC you could add the following to the configuration file:

```toml
[cli-tools]
[cli-tools.fastqc]
fastqc_params = ["--min_length", "10"]
```

## Suppressing excess output

If you do not want to view the output from each of the commands, you can pass the `-q`option to the main script. (Before the task selection, like with `-C`). It is recommended to test the command first however, as this option sometimes suppresses useful error messages.

*Note: this doesn’t stop all output, only that of internal commands. If you don’t want any progress messages whatsoever, only errors and warnings from python, set `-q` twice. E.g. `-q -q` or `-qq`.*

## Troubleshooting

If it is taking too long to run or producing odd data, it is strongly recommended to delete or move the output directory somewhere else, then do a clean run. This is because some steps use the data in the output directory as input, so may be combining old and new data.

When using the `targetid` command, `samtools` is used to filter out duplicates and create the FASTQ file outputs. `samtools` is known to be picky about how the input file is layed out, so will fail. The program should keep running, but if you want these files, check the following in 
your input files:

-	For FASTA files, there is no `@` at the start of any of the header strings
-	There are no duplicate headers in any of the input files
