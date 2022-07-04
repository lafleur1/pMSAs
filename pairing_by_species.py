#actual pairing by species step of pMSA generation

import os
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product

#procedure:
#1) split into species seq /species seq for each msa
#2) ID shared species pairs
#3) Estimate output size MSA
#4) Generate dummy gap sequences per protein
#5) Assemble final pmsa

def save_a3m_format(seq_record_list, savename):
    #save as correct foramt
    print ('saving a3m with ', len(seq_record_list), ' lines')
    with open(savename, 'w') as f:
        for record in seq_record_list:
            f.write('>' + record.id + '\n')
            f.write(str(record.seq) + '\n')

def make_dictionary_of_alignments_by_species(a3m_file_name):
    #first 5 letters are the UniprotKB species code + unique 5 digit number for the protein
    species_dict = {}
    s1_msa = SeqIO.parse(a3m_file_name, "fasta")
    wt_record = next(s1_msa)
    for record in s1_msa:
        species_id = record.id
        if species_id[0:5] in species_dict:
            species_dict[species_id[0:5]].append(record)
        else:
            species_dict[species_id[0:5]] = [record]
    return wt_record, species_dict

def id_shared_species(dict_spec_aligns1, dict_spec_aligns2):
    #returns estimate of msa to be generated
    size_pairs = {}
    spec_one = set(dict_spec_aligns1.keys())
    spec_two = set(dict_spec_aligns2.keys())
    print ('species in one: ', len(spec_one))
    print ('species in two: ', len(spec_two))
    ids_all = spec_one.union(spec_two)
    print ('total species (and number pMSA lines): ', len(ids_all))
    for id in ids_all:
        l1 = 0
        l2 = 0
        if id in dict_spec_aligns1:
            l1 = len(dict_spec_aligns1[id])
        if id in dict_spec_aligns2:
            l2 = len(dict_spec_aligns2[id])
        size_pairs[id] = [l1, l2, l1 * l2]
    #report
    print ("rep 1 > 0, rep 2 > 0") #full complex per species
    print (sum([1 for x in size_pairs if size_pairs[x][0] >0 and size_pairs[x][1]>0]))
    #half complex cases 1 & 2
    print ("rep 1 > 0, rep 2 == 0")
    print(sum([1 for x in size_pairs if size_pairs[x][0] > 0 and size_pairs[x][1] == 0]))
    print ("rep 1 == 0, rep 2 > 0")
    print(sum([1 for x in size_pairs if size_pairs[x][0] == 0 and size_pairs[x][1] > 0]))
    return size_pairs

def construct_pmsa(a3m_file_name1, a3m_file_name2, pmsa_name):
    #constructs in order of 1-2
    #using dummy gaps of '-' of len wt for inserts
    wt_one, dict_one = make_dictionary_of_alignments_by_species(a3m_file_name1)
    wt_two, dict_two = make_dictionary_of_alignments_by_species(a3m_file_name2)
    spec_pairs_dict = id_shared_species(dict_one, dict_two)

    #save wt string for robetta input
    #print (str(wt_one.seq) + "/" + str(wt_two.seq))
    with open(pmsa_name.replace('.a3m', '_wt.txt'), 'w') as f:
        f.write(str(wt_one.seq) + "/" + str(wt_two.seq))
    dummy_seq_one = '-' * len(wt_one.seq)
    dummy_seq_two = '-' * len(wt_two.seq)

    #constrct list of seqrecords to save the a3m with biopython
    final_fasta_seqs = []
    #add wt line first
    wt_pair = SeqIO.SeqRecord(Seq(str(wt_one.seq) + str(wt_two.seq)), id = 'wt1_wt2')
    final_fasta_seqs.append(wt_pair)
    for id in spec_pairs_dict:
        #add all possible combinations to the msa
        if id in dict_one:
            records_1 = dict_one[id]
        else:
            records_1 = [SeqIO.SeqRecord(Seq(dummy_seq_one), id = 'dummy')]

        if id in dict_two:
            records_2 = dict_two[id]
        else:
            records_2 = [SeqIO.SeqRecord(Seq(dummy_seq_two), id = 'dummy')]

        #use itertools to get all the combinations
        lines_to_add = list(product(records_1, records_2))
        for pair in lines_to_add:
            new_seq = Seq(str(pair[0].seq) + str(pair[1].seq))
            new_id = pair[0].id + '_' + pair[1].id
            final_fasta_seqs.append(SeqIO.SeqRecord(new_seq, id = new_id))
    save_a3m_format(final_fasta_seqs, pmsa_name)

def generate_all_pmsa_pairs(project_name, project_name_2, complex_output_folder):
    #generate pmsa for all 3 pipeline combos between both projects

    if not os.path.isdir(complex_output_folder):
        os.mkdir(complex_output_folder)

    names = ["startMSAAlign_filtered50.a3m", "filteredOrigMSAFinalAlign_filtered50.a3m", "filteredFastaFinalAlign_filtered50.a3m"]
    for n in names:
        a3m_p1 = project_name + n
        a3m_p2 = project_name_2 + n
        pmsa_name = n.replace('.a3m', 'paired.a3m')
        construct_pmsa(a3m_p1, a3m_p2,  complex_output_folder + pmsa_name)
        #quality filter the msa
        quality_filter_with_hhfilter(complex_output_folder + pmsa_name)
        print ("-----------------------")



def quality_filter_with_hhfilter(a3m_path):
    #filtering paired msa with suggested pmsa filtering criteria from roesttafold complex example page
    #id = 90, 95 and cov 75, 50
    #adding in orig 'high quality' alignment filter
    #gaps_percent_seq_ids(a3m_path)
    for id in [90,95,99,100]:
        for cov in [75,50,35]:
            new_name = a3m_path.replace('.a3m', '_id' + str(id) + '_cov' + str(cov) + '_max_hhfilered.a3m')
            os.system("hhfilter -i " + a3m_path + ' -o ' + new_name + ' -id ' + str(id) + ' -cov ' + str(cov)  + ' -v 0')
            #look at number of lines in filtered file
            wt_one, dict_one = make_dictionary_of_alignments_by_species(new_name)
            print('max: ', id, cov, len(dict_one))
            new_name = a3m_path.replace('.a3m', '_id' + str(id) + '_cov' + str(cov) + '_min40_hhfilered.a3m')
            os.system("hhfilter -i " + a3m_path + ' -o ' + new_name + ' -id ' + str(id) + ' -cov ' + str(cov) + ' -v 0 + -qid 40 ')
            # look at number of lines in filtered file
            wt_one, dict_one = make_dictionary_of_alignments_by_species(new_name)
            print ('min 40 qid: ', id, cov, len(dict_one))
            print ("**************")
    #jsut min 40 with query & min 35 seq id
    new_name = a3m_path.replace('.a3m', '_min40idmin35gap_hhfilered.a3m')
    os.system("hhfilter -i " + a3m_path + ' -o ' + new_name + ' -cov 35 -v 0 + -qid 40 -id 100 -qsc 0')
    wt_one, dict_one = make_dictionary_of_alignments_by_species(new_name)
    print('min 40 qid & min 35 gap: ', id, cov, len(dict_one))
    #making one minimum quality filtered alignment (40 percent seq id & 35% gaps


def gaps_percent_seq_ids(a3m_file):
    #calculate percent id & gaps for a a3m alignment line
    s1_msa = SeqIO.parse(a3m_file, "fasta")
    wt_record = next(s1_msa)
    gaps = []
    seqids = []
    seq_len = len(str(wt_record.seq))
    for record in s1_msa:
        gaps.append(str(record.seq).count('-')/seq_len)
        seqids.append(sum([1 for char in str(record.seq) if char.isupper()]) / seq_len)
    plt.scatter(gaps, seqids)
    plt.xlabel('gaps')
    plt.ylabel('seq id')
    plt.show()
    return gaps, seqids
