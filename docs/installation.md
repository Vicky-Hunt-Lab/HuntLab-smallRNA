# Installation

You can obtain the code for the pipeline from [https://github.com/kieranr51/HuntLab-smallRNA](https://github.com/kieranr51/HuntLab-smallRNA). To download, click the code button and copy the URL. You can then download it to a computer with git installed on it by opening the command line and running:

```
$ git clone https://github.com/Vicky-Hunt-Lab/HuntLab-smallRNA.git
```

It is recommended to use [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) to install this pipeline. To perform the installation, ensure mamba is installed on the machine, navigate to the directory containing the code and environment.yml file, then run:

```
$ mamba env create -f environment.yml
```

This will create a conda environment called `huntlab-smallrna`. This can then be activated each time you want to use the software with the command:

```
$ mamba activate huntlab-smallrna
```

By default, this provides the `hlsmallrna` executable. This provides all of the tasks listed in this documentation, each one can be run by typing `hlsmallrna` followed by the name and options for the task.