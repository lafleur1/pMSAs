#Stockholm format MSA generation per heterodimer complex protein

#Use single-copy fasta db -> generate seed alignment -> align with seed alignment -> stockhom format msa

from Bio import SeqIO, AlignIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from itertools import chain

#hmmsuite functions for all stockholm generation steps
#install hhsuite for these to work
def run_phmmer(search_fasta, db_fasta, output_sto_name):
    #make seed alignment from search fasta and orthologs db file from OMA
    os.system("phmmer -A " + output_sto_name + ' ' + search_fasta + ' ' + db_fasta)

def run_hmmbuild(hmm_out_name, sto_name_in):
    #run hmmbuild to convert sto name to the seed alignment
    os.system('hmmbuild '  + hmm_out_name + ' ' + sto_name_in)

def run_hmmalign(seed_hmm, db_fasta, output_sto):
    os.system("hmmalign -o " + output_sto + ' ' + seed_hmm + ' ' + db_fasta)

def check_sto_is_complete(sto_lines):
    #basic formatting check for if this stockholm is complete at the end of the file
    # should have version number at front and // at end
    return 'STOCKHOLM' in sto_lines[0] and '//' == sto_lines[-1]

#makes rf string
#RF string is the reference alignment string
def get_gc_rf(stock_file_name):
    file_lines = open(stock_file_name, 'r')
    all_lines = [x.strip() for x in file_lines.readlines()]
    ongoing_line = ""
    for i in range(0, len(all_lines)):# in all_lines:
        if '#=GC RF' == all_lines[i][0:len('#=GC RF')]:
            #print (i)
            ongoing_line += all_lines[i].replace(" ","").split("RF")[1]
    return ongoing_line

#makes wt dict
def make_lookup_alignment_dict(wt_seq, rf_line):
    #For manually computing percent identity using the RF line in the stockholm file
    #makes a lookup dictionary of the wildtype and alignment line positions to remove wild-type gaps
    dict_lookup = {}
    on_wt_index = 0
    for i in range(0, len(rf_line)):
        if rf_line[i] == 'x':
            #is a wt position
            dict_lookup[i] = wt_seq[on_wt_index]
            on_wt_index += 1 #look for next wt position
        else:
            dict_lookup[i] = rf_line[i]
    return dict_lookup

def get_all_indices(list_search, search_val):
    all_inds = []
    for i in range(0, len(list_search)):
        if list_search[i] == search_val:
            all_inds.append(i)
    return all_inds

def parse_gs_tags(tag_block, block_one_lines):
    # make default gs tags block and look to see what the other rf tags are
    tags = []
    descriptions = []
    for line in tag_block:
        # #GS NAME DE Description
        # print (line)
        if line[0:5] != '#=GS ' and len(line) > 0:
            print("WTF")
            print(line)
        elif len(line) > 0:
            # print (line)
            line = line[5:]  # skip GS start
            split_gs = line.split(' DE ')
            tags.append(split_gs[0].replace(" ", ""))
            descriptions.append(split_gs[1][1:])

    weird_tags = []
    # now find the weird lines (RF, PP_cons)
    for line in block_one_lines:
        # #GS NAME DE Description
        if '#=GC ' in line and len(line) > 0:
            # check to see if id is in the tags
            if line[0].split(" ")[0] not in tags:
                weird_tags.append(line.split(" ")[1])

    return tags, descriptions, weird_tags

#class of per-line info in the stockholm file
#represents one orthologs alignment in the file
class alignment_line:
    def __init__(self, name, description, alignment_str=""):
        self.name = name
        self.description = description
        self.alignment = alignment_str
        self.gap_score = -1
        self.id_score = -1
        if '/' in self.name:
            self.species = self.name.split("/")[0]
        else:
            self.species = self.name

    def gap_score_alignment(self, wt_dict):
        # return gap score (count '-' - count '*')
        return 100 * (self.alignment.count('-') / (len(self.alignment) - self.alignment.count('.')))

    def percent_seq_id(self, lookup_wt_dict):
        total_aligned = len(self.alignment)
        sum_aligned = 0
        for i in range(0, len(self.alignment)):
            if self.alignment[i] == lookup_wt_dict[i]:
                sum_aligned += 1
        return 100 * (sum_aligned - self.alignment.count('.')) / (total_aligned - self.alignment.count('.'))

    def score_seq(self, wt_dict):
        self.gap_score = self.gap_score_alignment(wt_dict)
        self.id_score = self.percent_seq_id(wt_dict)


class stockholm_file:
    # the biopython stockholm opener does weird stuff to the gaps - I want to keep the exact gap info
    # plus I need the RF string for alignment calculations later
    def __init__(self, file_name, wt_fasta_name):
        # wt setup
        wt_seq = next(SeqIO.parse(wt_fasta_name, "fasta"))
        self.wt_seq = str(wt_seq.seq)
        self.lookup_dict = make_lookup_alignment_dict(wt_seq, get_gc_rf(file_name))

        sto_file = open(file_name, 'r')
        sto_file_lines = [x.strip() for x in sto_file.readlines()]

        #sto_version_number = sto_file_lines[0].replace("# ", "")
        self.complete = check_sto_is_complete(sto_file_lines)
        self.size = 0
        self.rf = None
        self.alignments = {}

        if self.complete:
            breaks = get_all_indices(sto_file_lines, '')
            tags, descriptions, weird_tags = parse_gs_tags(sto_file_lines[breaks[0]:breaks[1] + 1],
                                                       sto_file_lines[breaks[1] + 1:breaks[2] + 1])
            for i in range(0, len(tags)):
                self.alignments[tags[i]] = alignment_line(tags[i], descriptions[i])
            for w in weird_tags:
                self.alignments[w] = alignment_line(w, 'special')

            # break up all blocks in order
            for line in sto_file_lines[breaks[1]:]:
                # add to alignment_line in tracker
                if len(line) > 2:  # has something in the line
                    split_lines = line.split()
                    if split_lines[0][0] != '#':
                        # is not the posterior prob
                        self.alignments[split_lines[0]].alignment += split_lines[1]
                    elif split_lines[0] == '#=GC':
                        self.alignments[split_lines[1]].alignment += split_lines[2]


    def score_all_alignments(self):
        for tag in self.alignments:
            if self.alignments[tag].description != 'special':
                self.alignments[tag].score_seq(self.lookup_dict)


    def produce_heatmap_percent_covered(self):
        gap_scores = [self.alignments[tag].gap_score for tag in self.alignments if
                      self.alignments[tag].description != 'special']
        id_score = [self.alignments[tag].id_score for tag in self.alignments if
                    self.alignments[tag].description != 'special']
        species = [self.alignments[tag].species for tag in self.alignments if
                    self.alignments[tag].description != 'special']
        dummy_df = pd.DataFrame({'gap': gap_scores, 'seq_id': id_score, 'species': species})
        #dummy_df['spec'] = dummy_df['id'].apply(lambda x: x.split('/')[0])

        vals = np.arange(0, 105, 5)

        z = np.zeros((vals.shape[0], vals.shape[0])) #percent covered
        import math
        for i in range(0, vals.shape[0]):# in vals:
            for j in range(0, vals.shape[0]):#seqs_val in vals:
                gap_val = vals[i]
                seqs_val = vals[j]
                z[i,j] = dummy_df[(dummy_df.gap <= gap_val) & (dummy_df.seq_id >= seqs_val)].shape[0] / dummy_df.shape[0]

        plt.matshow(z)
        plt.colorbar()
        plt.xlabel('ID score')
        plt.ylabel('Gap score')
        plt.show()

        #looking at the three criteria from the rosettafold paper

        #criteria 1
        print (dummy_df[(dummy_df.gap <= 20) & (dummy_df.seq_id >= 55)].shape[0]/dummy_df.shape[0])
        print(dummy_df[(dummy_df.gap <= 20) & (dummy_df.seq_id >= 55)].shape[0],
              dummy_df[(dummy_df.gap <= 20) & (dummy_df.seq_id >= 55)].species.nunique())

                #criteria 2
        print(dummy_df[(dummy_df.gap <= 35) & (dummy_df.seq_id >= 40)].shape[0] / dummy_df.shape[0])
        print(dummy_df[(dummy_df.gap <= 35) & (dummy_df.seq_id >= 40)].shape[0],
              dummy_df[(dummy_df.gap <= 35) & (dummy_df.seq_id >= 40)].species.nunique())

        #criteria 3
        print(dummy_df[(dummy_df.gap <= 50) & (dummy_df.seq_id >= 25)].shape[0] / dummy_df.shape[0])
        print(dummy_df[(dummy_df.gap <= 50) & (dummy_df.seq_id >= 25)].shape[0],
              dummy_df[(dummy_df.gap <= 50) & (dummy_df.seq_id >= 25)].species.nunique())

    def scatter_plot_scores(self):
        gap_scores = [self.alignments[tag].gap_score for tag in self.alignments if self.alignments[tag].description != 'special']
        id_score = [self.alignments[tag].id_score for tag in self.alignments if self.alignments[tag].description != 'special']
        plt.scatter(gap_scores, id_score)
        plt.xlabel('Gap score')
        plt.ylabel("id score")
        plt.show()

    def fitler_tags(self, gap_filter, id_filter):
        tags_len = [self.alignments[tag].name for tag in self.alignments if self.alignments[tag].gap_score <= gap_filter and self.alignments[tag].id_score >= id_filter]
        spec_len = [self.alignments[tag].species for tag in self.alignments if self.alignments[tag].gap_score <= gap_filter and self.alignments[tag].id_score >= id_filter]
        if len(tags_len) == len(spec_len):
            return tags_len, spec_len
        else:
            #not a case for septin proteins but with others will need to do the fancy mutli-alignment phmmer correction from the paper
            print ("ERROR! Need to do combo lines.... :(")
            print (len(tags_len), len(spec_len))
            return None, None

#functions for running differnt options of msa generation (differnet ways of creating seed alignment + final alignmnet)
def reduce_fasta(fasta_name, species_filtered, filtered_fasta_name):
    #reduce the fasta to just the proteins in the filtered tags list
    # checking how many orthologs from OMA are not paralogs - for S1, S5, S1ds (same as S1- will have to manually do the reciprocal check for the S1 isoform vs the proteome of the target species to verify that the selected is the best one)
    seqs = []
    ids = []

    seqs_ortho = SeqIO.parse(fasta_name, "fasta")
    for record in seqs_ortho:
        seqs.append(record)
        ids.append(str(record.id))

    seq_record_list = []
    # make new fasta of non-paralog seqs for phmmer step from rosettafold euk. paper
    for i in range(0, len(ids)):
        if ids[i] in species_filtered:
            seq_record_list.append(seqs[i])
    # now write to a single filtered fasta
    SeqIO.write(seq_record_list, filtered_fasta_name, "fasta")

def reduce_sto(tags_filtered, intermed_msa, filtered_orig_msa_name):
    #remove sto lines that are trashy
    #use SeqIO alignment reader for this - should be okay for the gaps here.....?
    old_sto = open(intermed_msa, 'r') #delete entries from readlines then rejoin with '/n' and save the string
    old_sto = [x.strip() for x in old_sto.readlines()]
    new_msa = []
    for line in old_sto:
        #check if the line should be stripped (kind or reverse of init?)
        splits = line.split()
        if len(splits) > 0 and '#=' in splits[0]: #alignment line
            if splits[1] not in tags_filtered:
                new_msa.append(line)
        elif len(splits) > 0 and splits[0] not in tags_filtered:
            new_msa.append(line)
        else:
            new_msa.append(line)
    #write out new msa
    with open(filtered_orig_msa_name, 'w') as f:
        f.write('\n'.join(new_msa))

def concatenate_fastas(fasta_1, fasta_2, new_name):
    #adds first fasta lines before the second and saves the concatenated fasta
    fasta_1 = SeqIO.parse(fasta_1, 'fasta')
    fasta_2 = SeqIO.parse(fasta_2, 'fasta')
    long_iterator = chain(fasta_1, fasta_2)
    SeqIO.write(long_iterator, new_name, 'fasta')


def all_options_stockholm_msa_generation(search_fasta_name, db_fasta_name, project_name):
    if not os.path.isdir(project_name):
        os.mkdir(project_name)
    #)0
    #make full db fasta (needs to have query as first sequence to do successful reformat.pl for .sto to .a3m conversion)
    concatenate_fastas(search_fasta_name, db_fasta_name, project_name + 'fullDB.fasta')

    #1) make initial msa
    intermed_msa = project_name + 'first_msa.sto'
    run_phmmer(search_fasta_name, db_fasta_name, intermed_msa)

    #2) Analyze the intermediate fasta
    s1 = stockholm_file(intermed_msa, search_fasta_name)
    s1.score_all_alignments()
    tags_filtered, spec_filtered = s1.fitler_tags(35, 40) #gap and percent id filters recommended for initial MSA seed creation from rosettafold paper

    #3) process initial fasta
    #Approach 1: make filtered fasta, remake a sto from that, use that sto to generate a seed msa, then align all the sequences
    filtered_fastsa_name = project_name + 'filteredFasta.fasta'
    reduce_fasta(db_fasta_name, spec_filtered, filtered_fastsa_name)
    #run phmmer on the filtered fasta
    filtered_fastsa_sto = project_name + 'filteredFastaAlign.sto'
    run_phmmer(search_fasta_name, filtered_fastsa_name, filtered_fastsa_sto)
    #generate seed hmm
    filtered_fastsa_seed = project_name + 'filteredFastaSeed.hmm'
    run_hmmbuild(filtered_fastsa_seed, filtered_fastsa_sto)
    #align all sequences
    pipeone_name = project_name + 'filteredFastaFinalAlign.sto'
    run_hmmalign(filtered_fastsa_seed, project_name + 'fullDB.fasta', pipeone_name)
    #analyze the scores for this alignment
    s1 = stockholm_file(pipeone_name, search_fasta_name)
    #s1.score_all_alignments()
    #s1.scatter_plot_scores()

    #Approach 2: filter the original msa down to just those lines
    filtered_orig_msa = project_name + 'filteredOrigMSA.sto'
    reduce_sto(tags_filtered, intermed_msa, filtered_orig_msa)
    # generate seed hmm with reduced sto
    filtered_orig_seed = project_name + 'filteredOrigMSA.hmm'
    run_hmmbuild(filtered_orig_seed, filtered_orig_msa)
    # align all sequences
    pipetwo_name = project_name + 'filteredOrigMSAFinalAlign.sto'
    run_hmmalign(filtered_orig_seed, project_name + 'fullDB.fasta', pipetwo_name)
    # analyze the scores for this alignment
    s1 = stockholm_file(pipetwo_name, search_fasta_name)
    #s1.score_all_alignments()
    #s1.scatter_plot_scores()

    #Approach 3: use the start msa as the seed for the alignment
    start_msa_seed = project_name + 'startMSASeed.hmm'
    run_hmmbuild(start_msa_seed, intermed_msa)
    # align all sequences
    pipethree_name = project_name + 'startMSAAlign.sto'
    run_hmmalign(start_msa_seed, project_name + 'fullDB.fasta', pipethree_name)
    # analyze the scores for this alignment
    s1 = stockholm_file(pipethree_name, search_fasta_name)
    #s1.score_all_alignments()
    #s1.scatter_plot_scores()








