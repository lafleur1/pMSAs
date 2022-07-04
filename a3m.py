#reformat stockholms to a3m files
#reformat to a3m format - which is what rosettafold expects formatwise

import os
from Bio import SeqIO, AlignIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import chain


def convert_sto_to_a3m(sto_file_name, a3m_file_name, path_to_reformat = "~/hh-suite/scripts/reformat.pl" ):
    #use hhsuite reformat.pl to transform it to a3m format
    os.system(path_to_reformat + " sto a3m " + sto_file_name + " " + a3m_file_name)

def reformat_all_ortholog_pipeline_msas(project_name):
    pipeone_name = project_name + 'filteredFastaFinalAlign.sto'
    pipeone_reformat = pipeone_name.replace('.sto', '.a3m')
    pipetwo_name = project_name + 'filteredOrigMSAFinalAlign.sto'
    pipetwo_reformat = pipetwo_name.replace('.sto', '.a3m')
    pipethree_name = project_name + 'startMSAAlign.sto'
    pipethree_reformat = pipethree_name.replace('.sto', '.a3m')

    #freformat names
    convert_sto_to_a3m(pipeone_name, pipeone_reformat)
    convert_sto_to_a3m(pipetwo_name, pipetwo_reformat)
    convert_sto_to_a3m(pipethree_name, pipethree_reformat)

def check_if_gaps_in_query(a3m_format):
    s1_msa = SeqIO.parse(a3m_format, "fasta")
    first = next(s1_msa)
    for gap_symbol in ['-', '*', '.']:
        if str(first.seq).count('-') != 0:
            return False
    return True

def check_if_gaps_in_pipeline_msas(project_name):
    pipeone_name = project_name + 'filteredFastaFinalAlign.sto'
    pipeone_reformat = pipeone_name.replace('.sto', '.a3m')
    pipetwo_name = project_name + 'filteredOrigMSAFinalAlign.sto'
    pipetwo_reformat = pipetwo_name.replace('.sto', '.a3m')
    pipethree_name = project_name + 'startMSAAlign.sto'
    pipethree_reformat = pipethree_name.replace('.sto', '.a3m')

    #freformat names
    print (pipeone_reformat, check_if_gaps_in_query(pipeone_reformat))
    filter_high_gap_lines(pipeone_reformat)
    print(pipetwo_reformat, check_if_gaps_in_query(pipetwo_reformat))
    filter_high_gap_lines(pipetwo_reformat)
    print(pipethree_reformat, check_if_gaps_in_query(pipethree_reformat))
    filter_high_gap_lines(pipethree_reformat)

#remove any lines from msa with high gaps (> 50)
def filter_high_gap_lines(a3m_file_name):
    keep_fasta_lines = []
    lens_all = []
    s1_msa = SeqIO.parse(a3m_file_name, "fasta")
    first = next(s1_msa)
    keep_fasta_lines.append(first)
    lens_all.append(len(first.seq))
    total_seqs = 0
    distrib_gaps = []
    for record in s1_msa:
        total_seqs += 1
        len_seq = len(str(record.seq))
        lens_all.append(len_seq)
        count_gaps = str(record.seq).count('-')
        distrib_gaps.append(count_gaps / len_seq)
        if count_gaps / len_seq <= 0.5:
            keep_fasta_lines.append(record)

    print ('total seqs: ', total_seqs, ' kept: ', len(keep_fasta_lines) - 1)
    #save filtered a3m
    print ('max percent gaps: ', max(distrib_gaps))
    print ('min percent gaps: ', min(distrib_gaps))
    print ('mean percent gaps: ', sum(distrib_gaps)/len(distrib_gaps))

    print('max percent gaps: ', max(lens_all))
    print('min percent gaps: ', min(lens_all))
    print('mean percent gaps: ', sum(lens_all) / len(lens_all))
    new_name = a3m_file_name.replace('.a3m', '_filtered50.a3m')
    SeqIO.write(keep_fasta_lines, new_name, 'fasta')
