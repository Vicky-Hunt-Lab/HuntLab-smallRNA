# CLI Reference
References of how to use the command line interface of the program and notes on some ways particular commands work that could trip you up.

General help message is as follows:

```
usage: hlsmallrna [-h] [-q] [-C PATH_TO_CONFIG]
              	{process,sort,extractnc,unitas,targetid,all} ...

Pipeline to process small RNAs

positional arguments:
  {process,sort,extractnc,unitas,targetid,all}
	process         	Preprocessing for the RNA
	sort            	Find RNAs that align to a genome and sort them by
                    	length
	extractnc       	Extract the noncoding region from a fasta with a GFF
                    	file
	unitas          	Run unitas on split files and merge results
	targetid        	Align small RNA to a number of genome features to find
                    	out what is targeted
	all             	Run process, sort and unitas one after the other

optional arguments:
  -h, --help        	show this help message and exit
  -q, --quiet       	Suppress output from intermediate commands
  -C PATH_TO_CONFIG, --path-to-config PATH_TO_CONFIG
                    	Path to the TOML format config file to use
```