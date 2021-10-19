#!/usr/bin/env python3

from argparse import ArgumentParser

from Bio import SeqIO

def append_label_to_id(label, seqs):
    for seq in seqs:
        seq.id = label + '|' + seq.id

        yield seq

path_to_file = '/home/kmhr20/Reference files/TE_sequences.fasta'
data = SeqIO.parse(path_to_file, 'fasta')
labelled_data = append_label_to_id('TE', data)

SeqIO.write(labelled_data, 'TE_labelled.fasta', 'fasta')