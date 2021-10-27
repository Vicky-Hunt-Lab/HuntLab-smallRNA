from Bio import SeqIO
from gffutils.iterators import DataIterator

from sys import argv
from collections import defaultdict

from .config import do_log

def get_pairs(array):
    '''
    Return pairs of values from an array
    '''
    for i in range(0, len(array), 2):
        yield array[i], array[i + 1]

def extract_fragments(sequences, coordinates, start_dict, end_dict):
    '''
    Pull out fragments from a set of sequences based on a 
    dictionary of coordinates
    '''
    for seq in sequences:
        if seq.id + '+' in coordinates.keys():
            current_labels = []

            starts = sorted(start_dict[seq.id + '+'].keys())
            ends = sorted(end_dict[seq.id + '+'].keys())

            for i, j in get_pairs(coordinates[seq.id + '+']):
                while len(starts) > 0 and starts[0] <= j:
                    current_labels += start_dict[seq.id + '+'][starts.pop(0)]

                while len(ends) > 0 and ends[0] <= i:
                    for item in end_dict[seq.id + '+'][ends.pop(0)]:
                        current_labels.remove(item)

                fragement = seq[i:j - 1]
                fragement.id = 'noncoding|' + '|'.join(current_labels)

                if len(fragement) > 0:
                    yield fragement

        if seq.id + '-' in coordinates.keys():
            current_labels = []

            starts = sorted(start_dict[seq.id + '-'].keys())
            ends = sorted(end_dict[seq.id + '-'].keys())

            for i, j in get_pairs(coordinates[seq.id + '-']):
                while len(starts) > 0 and starts[0] <= j:
                    current_labels += start_dict[seq.id + '-'][starts.pop(0)]

                while len(ends) > 0 and ends[0] <= i:
                    for item in end_dict[seq.id + '-'][ends.pop(0)]:
                        current_labels.remove(item)

                fragement = seq[i:j - 1].reverse_complement()
                fragement.id = 'noncoding|' + '|'.join(current_labels)
                fragement.description = seq.description

                if len(fragement) > 0:
                    yield fragement

def validate_gff(mRNAs, coding_reigon):
    '''
    Check the assumption that the coding reigon is within the mRNA
    Crash if it isn't true as the algorithm will not work
    '''
    for key in coding_reigon.keys():
        for coding in coding_reigon[key]:
            valid = False

            for coord in mRNAs[key]:
                if coding[0] >= coord[0] and coding[1] <= coord[1]:
                    valid = True
                    break

            if not valid:
                raise Exception(f'GFF File Error : coding reigon outside mRNA coords. coding: {coding[0]} {coding[1]}')

def merge_cds(coding_reigon):
    '''
    Sometime CDS overlap, this just merges them into the longest reigon
    '''
    n_merged = 0
    for key in coding_reigon.keys():
        new_cds = []

        for reigon in coding_reigon[key]:
            merged = False
            
            for i in range(len(new_cds)):
                if (reigon[0] >= new_cds[i][0] and reigon[0] <= new_cds[i][1]) or (reigon[1] >= new_cds[i][0] and reigon[1] <= new_cds[i][1]):
                    new_cds[i] = [min(new_cds[i][0], reigon[0]), max(new_cds[i][1], reigon[1])]
                    merged = True

                    n_merged += 1

            if not merged:
                new_cds.append(reigon)

        coding_reigon[key] = new_cds

    return n_merged

def extract_noncoding(genome, gff_path, quiet=0, output='result.fasta'):
    '''
    Extract the noncoding reigon from the genome, basied on a GFF file
    '''
    gff_iter = DataIterator(gff_path)
    genome_data = SeqIO.parse(genome, 'fasta')

    do_log(quiet, '====> Calculating coordinates')

    coordinates = defaultdict(lambda: [])
    mRNAs = defaultdict(lambda: [])
    coding_reigon = defaultdict(lambda: [])

    mRNA_start = defaultdict(lambda: defaultdict(lambda: []))
    mRNA_end = defaultdict(lambda: defaultdict(lambda: []))

    for item in gff_iter:
        if item[2] == 'mRNA':
            mRNAs[item[0] + item[6]].append([int(item[3]), int(item[4])])

            mRNA_start[item[0] + item[6]][int(item[3])].append(item[8]['ID'][0])
            mRNA_end[item[0] + item[6]][int(item[4])].append(item[8]['ID'][0])
        
        if item[2] == 'CDS':
            coding_reigon[item[0] + item[6]].append([int(item[3]), int(item[4])])

    do_log(quiet, '====> Merging and validating coordinates')

    cds_merged = merge_cds(coding_reigon)
    mRNA_merged = merge_cds(mRNAs)
    validate_gff(mRNAs, coding_reigon)

    do_log(quiet, f'Merged {cds_merged} coding reigons and {mRNA_merged} mRNAs')

    for key in mRNAs.keys():
        for item in mRNAs[key]:
            coordinates[key].append(item[0])
            coordinates[key].append(item[1] + 1)

    for key in coding_reigon.keys():
        for item in coding_reigon[key]:
            coordinates[key].append(item[0])
            coordinates[key].append(item[1])

    do_log(quiet, '====> Extarcting fragments')
    
    for key in coordinates.keys():
        coordinates[key].sort()

    SeqIO.write(extract_fragments(genome_data, coordinates, mRNA_start, mRNA_end), output, 'fasta')

if __name__ == '__main__':
    extract_noncoding(argv[1], argv[2])