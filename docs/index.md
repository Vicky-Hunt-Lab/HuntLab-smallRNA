---
hide:
  - navigation
  - toc
---
# Welcome to the Hunt Lab Small RNA Pipeline Documentation

This documentation is for version 2 of the pipeline, select 1.1.0 from the version dropdown for the version 1 documentation.

This is a pipeline developed in [Vicky Hunt's Lab](https://www.vickyhuntlab.org/) at the [University of Bath](https://www.bath.ac.uk), UK to analyse small RNA data in a quick and consistent way between projects. It is deliberately designed to be highly configurable so it can be easily applied to a wide range of scenarios, though most of our work is on nematode and nematomorph worms. 

# Background

Most of the original methodology was devised by Mona Suleiman, Dominika Lastik and Vicky Hunt for the papers below. This was then developed into a usable pipeline by Kieran Reynolds in 2022 and updated to version 2 in 2025. 

> Suleiman, M., Kounosu, A., Murcott, B. et al. piRNA-like small RNAs target transposable elements in a Clade IV parasitic nematode. Sci Rep 12, 10156 (2022). [https://doi.org/10.1038/s41598-022-14247-1](https://doi.org/10.1038/s41598-022-14247-1)

> Lastik, D., Kounosu, A., Dayi, M. et al. Small non-coding RNAs have predicted roles in reproductive biology and transposable element regulation in the parasitic worm Strongyloides venezuelensis. Sci Rep 15, 20608 (2025). [https://doi.org/10.1038/s41598-025-01968-2](https://doi.org/10.1038/s41598-025-01968-2)

# Citation

This pipeline is mainly a combination of existing tools so, doesn't have a paper attached to it but can be referenced with the repository URL and version number. For example:

> This analysis was performed using HuntLab-smallRNA v2.0.0 ([https://github.com/Vicky-Hunt-Lab/HuntLab-smallRNA](https://github.com/Vicky-Hunt-Lab/HuntLab-smallRNA)).

For the internal steps please cite the papers of the tools used:

**Trim**:

> Martin, Marcel. (2011). CUTADAPT removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 17. 10.14806/ej.17.1.200.

**Sort**:

> Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357-359 (2012). https://doi.org/10.1038/nmeth.1923 

**Unitas**:

> Gebert, D., Hewel, C. & Rosenkranz, D. unitas: the universal tool for annotation of small RNAs. BMC Genomics 18, 644 (2017). https://doi.org/10.1186/s12864-017-4031-9

**Targetid**:

> Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357-359 (2012). https://doi.org/10.1038/nmeth.1923

**Enrich**:

> Carlos P Cantalapiedra, Ana Hernández-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas, eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale, Molecular Biology and Evolution, Volume 38, Issue 12, December 2021, Pages 5825-5829, https://doi.org/10.1093/molbev/msab293
