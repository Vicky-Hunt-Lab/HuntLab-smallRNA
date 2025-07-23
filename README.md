# Small RNA Pipeline for the Hunt Lab

This repository contains the code for the pipeline used to carry out sRNA-seq ananlysis in Vicky Hunt's lab. This is version 2 for version 1 see the `v1.x` branch.

Documentation and information on how to install and run this pipeline can be found [here](https://vicky-hunt-lab.github.io/HuntLab-smallRNA/latest). Please report any bugs found or ask any other questions relating to the pipeline on the GitHub bug tracker.

## Installation

A containerised version of the pipeline is provided to allow for easy installation. Apptainer is recommended, as it allows for running with relative paths.

### Installing with Apptainer
Create a SIF file with the code by running the following:

```
$ apptainer pull hlsmallrna.sif docker://kieranr51/huntlab-smallrna
```

To test it works try displaying the help message with the following command:

```
$ apptainer run hlsmallrna.sif hlsmallrna --help
```

If this displays the help message, you are good to go. All subsequent commands can be run with:

```
$ apptainer run -B $HOME:$HOME hlsmallrna.sif [COMMAND]
```

If this is not possible, other methods of installation are covered in the documentation at:
[https://vicky-hunt-lab.github.io/HuntLab-smallRNA/latest/installation/](https://vicky-hunt-lab.github.io/HuntLab-smallRNA/latest/installation/)