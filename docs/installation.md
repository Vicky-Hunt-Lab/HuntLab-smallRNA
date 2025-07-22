# Installation

There are currently three different ways to install the pipeline, to allow for flexibility. The apptainer method is recommended as it is the quickest and most reliable, but the other two work if that is not possible.

## Installing using Apptainer or Docker
A containerised version of the pipeline is provided to allow for easy installation. Apptainer is recommended, as it allows for running with relative paths, but it will run in Docker as well. If you do not have these installed, see their websites for installation instructions:

-	Apptainer:  [https://apptainer.org/docs/admin/main/installation.html](https://apptainer.org/docs/admin/main/installation.html) 
-	Docker: [https://docs.docker.com/](https://docs.docker.com/)

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

!!! note "Binding Directories in Apptainer"
    To allow access to your files to the pipeline in apptainer you will need to bind the directories for you input and output. The above command binds your home directory, but you may need to change this if your files are outside of that directory. See [https://apptainer.org/docs/user/main/bind_paths_and_mounts.html](https://apptainer.org/docs/user/main/bind_paths_and_mounts.html) for more information.

### Installing with Docker

First pull the container with:

```
$ docker pull kieranr51/huntlab-smallrna
```

To test it works try displaying the help message with the following command:

```
$ docker run kieranr51/huntlab-smallrna hlsmallrna --help
```

If this displays the help message, you are good to go. All subsequent commands can be run with:

```
$ docker run -v $HOME:$HOME kieranr51/huntlab-smallrna [COMMAND]
```

!!! note "Paths and Mounting Directories in Docker"
    To allow access to your files to the pipeline in apptainer you will need to mount the directories for you input and output. The above command mounts your home directory, but you may need to change this if your files are outside of that directory. See [https://docs.docker.com/engine/storage/bind-mounts/](https://docs.docker.com/engine/storage/bind-mounts/) for more information.

    Docker hides the current path you are working in from the program as well, so you will have to specify the absolute path to your input and output files to get it to work. E.g. if your config file is in /home/user/smallrna/config.yml you would have to run:

    ```
    $ docker run -v $HOME:$HOME kieranr51/huntlab-smallrna hlsmallrna -o /home/user/smallrna/hlsmallrna_output /home/user/smallrna/config.yml
    ```

## Installing using Mamba or Conda

If you can't use the containers, using mamba or conda to set up dependencies is the next supported route. Mamba is preferred. If you do not have these installed see the following links:

-	Mamba: [https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
-	Conda: [https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

Then download the latest version of HuntLab-smallRNA from [https://github.com/Vicky-Hunt-Lab/HuntLab-smallRNA/releases](https://github.com/Vicky-Hunt-Lab/HuntLab-smallRNA/releases). Extract the zipped file, then navigate into the directory. To install run one of the following commands:

```
$ mamba env create -f environment.yml
```

```
$ conda env create -f environment.yml
```

This will create an environment called huntlab-smallrna with the pipeline and dependencies installed into it. To test activate it with one of these commands:

```
$ mamba activate huntlab-smallrna
```

```
$ conda activate huntlab-smallrna
```

Note you will have to run the above command each time you start your computer, before you use this software. Then you should be able to display the help message with the command:

```
$ hlsmallrna --help
```

## Manual Installation
If none of the above are possible, you can retrieve each dependency and manually install it. See the following list of dependencies to install:

-	Python 3 - [https://www.python.org/](https://www.python.org/)
-	Pyyaml - [https://pyyaml.org/](https://pyyaml.org/)
-	Pysam - [https://pysam.readthedocs.io/en/latest/api.html](https://pysam.readthedocs.io/en/latest/api.html)
-	Biopython - [https://biopython.org/](https://biopython.org/)
-	Gffutils - [https://daler.github.io/gffutils/](https://daler.github.io/gffutils/)
-	Scipy - [https://scipy.org/](https://scipy.org/) 
-	Matplotlib - [https://matplotlib.org/](https://matplotlib.org/)
-	Cutadapt - [https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)
-	Unitas - [https://www.smallrnagroup.uni-mainz.de/software.html](https://www.smallrnagroup.uni-mainz.de/software.html)
-	Bowtie2 - [https://bowtie-bio.sourceforge.net/bowtie2/index.shtml](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
-	Samtools - [http://www.htslib.org/](http://www.htslib.org/)
-	Bedtools - [https://bedtools.readthedocs.io/en/latest/](https://bedtools.readthedocs.io/en/latest/)
-	Eggnog-mapper - [https://github.com/eggnogdb/eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)
-	R -[https://www.r-project.org/](https://www.r-project.org/)
-	topGO - [https://bioconductor.org/packages/3.22/bioc/html/topGO.html](https://bioconductor.org/packages/3.22/bioc/html/topGO.html)
-	pip - [https://packaging.python.org/en/latest/tutorials/installing-packages/](https://packaging.python.org/en/latest/tutorials/installing-packages/)

Once all of these are installed you need to download and extract the latest version of HuntLab-smallRNA from [https://github.com/Vicky-Hunt-Lab/HuntLab-smallRNA/releases](https://github.com/Vicky-Hunt-Lab/HuntLab-smallRNA/releases). Navigate into the extracted directory and run the following command to install:

```
$ pip install .
```

Then you should be able to display the help message with the command:

```
$ hlsmallrna --help
```