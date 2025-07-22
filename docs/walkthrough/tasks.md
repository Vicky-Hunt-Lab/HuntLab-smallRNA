# Tasks

Currently the following tasks have been implemented in it:

-	Process - Trim adapters off raw RNA-Seq data, run a quality control report and trim of any bases that fall below a specified threshold from all the RNA

-	Sort - Align the RNA against a given reference, removing any sequences that fail to align, then split them into files based on sequence length to allow for easy further processing. Can also remove sequences below a minimum length and above a maximum length if desired

-	ExtractNC - Using a genome and GFF file containing annotations, extract the region that is transcribed, but not a coding region

-	Unitas - Run a directory of RNA files (e.g. output of sort step) though unitas to classify the types of RNA, then combine the summaries into one spreadsheet for easy viewing and manipulation

-	TargetID - Identify which of a set of potential targets a set of small RNA are targeting by reverse complimenting them and aligning them to the potential targets

More in depth documentation on each of these steps is provided, starting with examples and suggestions for further processing. There is also a CLI Reference, that contains general help for each command, including what files are required and what outputs you get from each step.
