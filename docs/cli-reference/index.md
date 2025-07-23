# CLI Reference
References of how to use the command line interface of the program and notes on some ways particular commands work that could trip you up.

General help message is as follows:

```
usage: hlsmallrna [-h] [-o OUTPUT] [-t THREADS] [-k] [-v] config_file

Pipeline to process small RNAs - version 2

positional arguments:
  config_file           Path to YAML config file

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Directory to write output to (will be created if it
                        doesn't exist, default: hlsmallrna_output)
  -t THREADS, --threads THREADS
                        Number of threads to run the pipeline on (default: 4)
  -k, --keep-files      If set, keep the intermediate files, to help debug the
                        application
  -v, --verbose         If set, print the output for the intermediate commands
```

!!! note "The Config File"
	If you are looking for instructions on writing the config file, use [the relevent page of the walkthrough](../walkthrough/config.md) to help.