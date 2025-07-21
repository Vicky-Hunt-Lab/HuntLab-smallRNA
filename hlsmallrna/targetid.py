# Copyright 2022 - 2025 Vicky Hunt Lab Members
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os
import os.path
import csv

from subprocess import run
from collections import defaultdict, Counter

import pysam

from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import fisher_exact

from .kegg_pathway_names import KEGG_PATHWAY_NAMES

GO_ENRICHMENT_SCRIPT = '''
interesting_genes <- c({{INTERSTING_GENE_LIST}})

library(topGO)

all_data <- readMappings("{{PATH_TO_GO_TXT}}")
universe <- names(all_data)
duplicated_list <- factor(as.integer(universe %in% interesting_genes))
names(duplicated_list) <- universe

go_ids <- c()
go_signif <- c()
go_names <- c()

for (ont in c('BP', 'CC', 'MF')) {
  myGOdata <- new("topGOdata", description="My project", ontology=ont, allGenes=duplicated_list, annot = annFUN.gene2GO, gene2GO = all_data)
  significant <- sigGenes(myGOdata)
  resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")
  mysummary <- summary(attributes(resultFisher)$score <= 0.01)
  
  if (length(mysummary) > 2) {
    numsignif <- as.integer(mysummary[[3]]) 
    
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = numsignif)
    go_ids <- c(go_ids, allRes$GO.ID)
    go_names <- c(go_names, allRes$Term)
    go_signif <- c(go_signif, allRes$classicFisher)
  }
}

lines <- c()
for (i in 1:length(go_ids)) {
  lines <- append(lines, c(paste(go_ids[i], go_names[i], go_signif[i], sep = "\t")))
}
fileConn <- file("{{PATH_TO_OUTPUT_TSV}}")
writeLines(lines, fileConn)
close(fileConn)
'''

def check_ids_unique(path_to_file, format):
    '''
    Make sure that the sequences in a FASTA or FASTQ file are unique
    '''
    seen_ids = set()
    for seq in SeqIO.parse(path_to_file, format):
        if seq.id in seen_ids:
            return False
        else:
            seen_ids.add(seq.id)

    return True

def seqio_merge_multiple(files, filetype):
    '''
    Merge SeqIO.parse over multiple files
    '''
    for filepath in files:
        for seq in SeqIO.parse(filepath, filetype):
            yield seq

def seqio_revcomp_multiple(fastq_files):
    '''
    Merge SeqIO.parse over multiple files and reverse complement sequences
    '''
    seen_ids = set()

    for seq in seqio_merge_multiple(fastq_files, 'fastq'):
        if seq.id in seen_ids:
            print(f'Error: duplicate ID {seq.id} found in your small RNA files, make IDs unique and try again')
            exit(1)
        else:
            seen_ids.add(seq.id)

        revcomp = seq.reverse_complement()
        revcomp.id = seq.id

        yield revcomp

def merge_and_revcomp_input_files(smallRNA_files, output):
    '''
    Create a file containing the reverse complement of a file of small RNA
    '''

    print('====> Reverse complimenting RNA...')
    PATH_TO_REVCOMP = os.path.join(output, 'revcomp_rna.fastq')

    SeqIO.write(seqio_revcomp_multiple(smallRNA_files), PATH_TO_REVCOMP, 'fastq')

    return PATH_TO_REVCOMP

def find_targets(smallRNA, possible_target_list, output, program_paths, threads=4, min_seq_length=2, mismatches_allowed=0, verbose=False):
    '''
    Align the small RNA against the lists of possible targets with bowtie2 and
    analyse the output
    '''

    CWD = os.getcwd()
    INDEX_DIR = os.path.join(output, 'bowtie_indexes')
    SAM_DIR = os.path.join(CWD, output, 'target_alignments')

    sam_files = []

    if not os.path.exists(SAM_DIR):
        os.makedirs(SAM_DIR)

    if not os.path.exists(INDEX_DIR):
        os.makedirs(INDEX_DIR)

    os.chdir(INDEX_DIR)

    for target in possible_target_list:
        if not check_ids_unique(os.path.join(CWD, target), 'fasta'):
            print(f'Error: duplicate ID found in your target file: {target}. Make IDs unique and try again')
            exit(1)

        print(f'====> Building index for {target}...')
        INDEX_NAME = os.path.basename(target) + '_index'

        bowtie2_build_command = [
            program_paths['bowtie2-build'],
            '--threads', str(threads),
            os.path.join(CWD, target),
            INDEX_NAME
        ]

        run(bowtie2_build_command, capture_output=not verbose)

        print(f'====> Aligning small RNA against {target}...')

        sam_files.append(os.path.join(SAM_DIR, INDEX_NAME + '.sam'))

        if mismatches_allowed > 0:
            bowtie2_align_command = [
                program_paths['bowtie2'],
                '--threads', str(threads),
                '-L', str(min_seq_length),
                '--score-min', 'L,-' + str(mismatches_allowed) + ',0',
                '--end-to-end',
                '--norc',
                '--mp', '1,1',
                '--ignore-quals',
                '--rdg', '100,100',
                '--rfg', '100,100',
                '-x', INDEX_NAME,
                '-U', os.path.join(CWD, smallRNA),
                '-S', sam_files[-1],
                '-a'
            ]
        else:
            bowtie2_align_command = [
                program_paths['bowtie2'],
                '--threads', str(threads),
                '-L', str(min_seq_length),
                '--no-1mm-upfront',
                '--score-min', 'L,0,0',
                '--end-to-end',
                '--norc',
                '-M', '0',
                '-x', INDEX_NAME,
                '-U', os.path.join(CWD, smallRNA),
                '-S', sam_files[-1],
                '-a'
            ]

        run(bowtie2_align_command, capture_output=not verbose)

    os.chdir(CWD)

    return sam_files

def build_summary_files(sam_files, output, program_paths, verbose=False):
    '''
    Takes the alignment SAM and extracts the sequences into a summary file and
    a set of FASTA files
    '''
    print('====> Removing unaligend sequences and building summary files...')

    target_rows = defaultdict(lambda: [])
    # Create a set of targets to analyse for enrichment later
    target_set = set()
    for filename in sam_files:
            samfile = pysam.AlignmentFile(filename, 'r')
            fileheader = samfile.header

            for read in samfile:
                if not read.is_unmapped:
                    query_name = read.query_name
                    query_sequence = Seq(read.query_sequence).reverse_complement()
                    r_name = fileheader.get_reference_name(read.reference_id)
                    start = read.reference_start
                    end = read.reference_end

                    if read.is_reverse:
                        strand = 'antisense'
                    else:
                        strand = 'sense'

                    length = len(query_sequence)
                    first_base = query_sequence[0].replace('T', 'U')
                    target_rows[f'{length}{first_base}'].append([query_name, query_sequence, os.path.abspath(filename), r_name, start, end, strand])
                    target_set.add(r_name)

    TRAGETS_DIR = os.path.join(output, 'targets_by_group')

    if not os.path.exists(TRAGETS_DIR):
        os.makedirs(TRAGETS_DIR)

    for target in target_rows:
        with open(os.path.join(TRAGETS_DIR, f'{target}_targets.tsv'), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['Query Name', 'Query Sequence', 'Target File', 'Target Sequence', 'Start Coordinate', 'End Coordinate', 'Strand'])
            writer.writerows(target_rows[target])
    
    with open(os.path.join(output, 'rna_target_list.tsv'), 'w') as tablefile:
        writer = csv.writer(tablefile, delimiter='\t')
        writer.writerow(['Query Name', 'Query Sequence', 'Target File', 'Target Sequence', 'Start Coordinate', 'End Coordinate', 'Strand'])
        for target in target_rows:
            writer.writerows(target_rows[target])
    try:
        remove_duplicate_command = [
            program_paths['samtools'],
            'view',
            '-h', 
            '-F', '256',
            '-F', '4',
            '-o', filename + '.nodup.sam',
            filename
        ]
    
        result = run(remove_duplicate_command, capture_output=not verbose)
        result.check_returncode()
        filename = filename + '.nodup.sam'

        bam_to_fastq_command = [
            program_paths['samtools'],
            'fastq',
            filename
        ]

        result = run(bam_to_fastq_command, capture_output=True)
        result.check_returncode()
        with open(filename + '.fastq', 'wb') as f:
            f.write(result.stdout)

    except:
        print('====> SAMTOOLS convertion failed, skipping creating filtered and FASTQ files...')

    return target_set

def do_term_enrichment(files_to_annotate, target_list, path_to_emapper_data, output, program_paths, threads=4, verbose=False):
    '''
    Use eggnog-mapper to annotate GO, KEGG and PFam with diamond blastx to provide a quick overview enrichment of terms
    '''
    EGGNOG_MAPPER_OUTPUT = os.path.join(output, 'eggnog_mapper')
    COMBINED_TARGETS_FILE = os.path.join(output, 'eggnog_mapper', 'targets.fasta')

    if not os.path.exists(EGGNOG_MAPPER_OUTPUT):
        os.makedirs(EGGNOG_MAPPER_OUTPUT)

    print('====> Merging target files...')
    SeqIO.write(seqio_merge_multiple(files_to_annotate, 'fasta'), COMBINED_TARGETS_FILE, 'fasta')

    print('====> Annotating with eggnog-mapper...')
    eggnog_mapper_command = [
        program_paths['eggnog-mapper'], '-i', COMBINED_TARGETS_FILE, '-o', 'eggnog_mapper', '--output_dir', EGGNOG_MAPPER_OUTPUT, '--cpu', str(threads), 
        '--data_dir', path_to_emapper_data, '--itype', 'CDS'
    ]
    run(eggnog_mapper_command, capture_output=not verbose)

    print('====> Reading eggnog-mapper output...')
    go_terms = {}

    kegg_pathways = {}
    kegg_counts = defaultdict(lambda: 0)
    total_kegg_annot = 0

    pfams = {}
    pfam_counts = defaultdict(lambda: 0)
    total_pfam_annot = 0
    with open(os.path.join(EGGNOG_MAPPER_OUTPUT, 'eggnog_mapper.emapper.annotations')) as f:
        reader = csv.reader(f, delimiter='\t')

        for line in reader:
            if not line[0].startswith('#'):
                if line[9] == '-':
                    go_terms[line[0]] = set()
                else:
                    go_terms[line[0]] = set(line[9].split(','))

                kegg_pathways[line[0]] = set()
                for kegg in line[12].split(','):
                    if 'map' in kegg:
                        kegg_pathways[line[0]].add(kegg)
                        kegg_counts[kegg] += 1
                        total_kegg_annot += 1

                if line[-1] == '-':
                    pfams[line[0]] = set()
                else:
                    pfam_set = set(line[-1].split(','))
                    pfams[line[0]] = pfam_set

                    total_pfam_annot += len(pfam_set)

                    for pfam in pfam_set:
                        pfam_counts[pfam] += 1

    with open(os.path.join(EGGNOG_MAPPER_OUTPUT, 'go_map.txt'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for name in go_terms:
            writer.writerow([name, ','.join(go_terms[name])])

    with open(os.path.join(EGGNOG_MAPPER_OUTPUT, 'kegg_map.txt'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for name in kegg_pathways:
            writer.writerow([name, ','.join(kegg_pathways[name])])

    with open(os.path.join(EGGNOG_MAPPER_OUTPUT, 'pfam_map.txt'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for name in pfams:
            writer.writerow([name, ','.join(pfams[name])])

    genes_to_analyse_string = ''
    for gene in target_list:
        genes_to_analyse_string += f',"{gene}"'

    print('====> Enriching GO terms using TopGO...')
    customised_go_enrichment = GO_ENRICHMENT_SCRIPT.replace(
        '{{INTERSTING_GENE_LIST}}', genes_to_analyse_string[1:],  1
    ).replace(
        '{{PATH_TO_GO_TXT}}', os.path.join(EGGNOG_MAPPER_OUTPUT, 'go_map.txt'), 1
    ).replace(
        '{{PATH_TO_OUTPUT_TSV}}', os.path.join(output, 'targets_enriched_go_terms.tsv'), 1
    )

    GO_ENRICHMENT_PATH = os.path.join(EGGNOG_MAPPER_OUTPUT, 'enrichGOTerms.R')
    with open(GO_ENRICHMENT_PATH, 'w') as f:
        f.write(customised_go_enrichment)

    go_enrich_command = [
        program_paths['rscript'], GO_ENRICHMENT_PATH
    ]
    run(go_enrich_command, capture_output=not verbose)

    print('====> Enriching KEGG pathways with scipy fisher_exact...')
    kegg_output_table = []

    kegg_in_set = []
    for interesting_gene in target_list:
        if interesting_gene in kegg_pathways:
            kegg_in_set += kegg_pathways[interesting_gene]

    kegg_in_total = len(kegg_in_set)
    kegg_counter = Counter(kegg_in_set)

    for kegg in kegg_counter:
        in_intersting = kegg_counter[kegg]
        out_intesting = kegg_counts[kegg] - in_intersting

        in_not_intersting = kegg_in_total - in_intersting
        out_not_intersting = total_kegg_annot - out_intesting - in_intersting - in_not_intersting

        fisher_table = [[in_intersting, out_intesting], [in_not_intersting, out_not_intersting]]
        result = fisher_exact(fisher_table, alternative='greater')

        if result.pvalue < 0.01:
            kegg_output_table.append([kegg, KEGG_PATHWAY_NAMES[kegg], result.pvalue])
        
    with open(os.path.join(output, 'targets_enriched_kegg_pathways.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(kegg_output_table)

    print('====> Enriching PFams with scipy fisher_exact...')
    pfam_output_table = []

    pfam_in_set = []
    for interesting_gene in target_list:
        if interesting_gene in pfams:
            pfam_in_set += pfams[interesting_gene]

    pfam_in_total = len(pfam_in_set)
    pfam_counter = Counter(pfam_in_set)

    for pfam in pfam_counter:
        in_intersting = pfam_counter[pfam]
        out_interesting = pfam_counts[pfam] - in_intersting

        in_not_intersting = pfam_in_total - in_intersting
        out_not_intersting = total_pfam_annot - in_intersting - out_interesting - in_not_intersting

        fisher_table = [[in_intersting, out_interesting], [in_not_intersting, out_not_intersting]]
        result = fisher_exact(fisher_table, alternative='greater')

        if result.pvalue < 0.01:
            pfam_output_table.append([pfam, result.pvalue])

    with open(os.path.join(output, 'targets_enriched_pfams.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(pfam_output_table)