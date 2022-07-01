#make sure blast is installed to run these

import os
import numpy as np
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def run_pblast_self_v_self(fasta_file, save_name, output_dir = './blastp_outputs/'):
    #runs pBLAST for a fasta vs. its self (all by all alignment) and saves alignment stats as a csv in the output_dir
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    _ = NcbiblastpCommandline(query=fasta_file,
                          subject=fasta_file,
                          out=output_dir + save_name + '.tab', outfmt=6, num_threads = 8)()[0]
    #adding columns to the output
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
            'bitscore']
    output = pd.read_csv(output_dir + save_name  + '.tab', names=cols, sep="\t")
    #removing self-self blast
    output = output[output.qseqid != output.sseqid]
    # fitler out same len duplicate pairs
    output.to_csv(output_dir + save_name + '.tab', index=False, sep="\t")
    return output

def run_blastp_all_v_all(species_id, ortholog_db, output_dir = './blastp_outputs/', fasta_dir = './single_family_fastas/'):
    #make a temp fasta and align it with itself with blastp
    dummy_list = []
    ids = []
    lens = []
    #set up outputs
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(fasta_dir):
        os.mkdir(fasta_dir)

    for ortho in ortholog_db[species_id]:
        dummy_list.append(SeqRecord(Seq(ortho.sequence),
                   id=ortho.oma_id))
        ids.append(ortho.oma_id)
        lens.append(len(ortho.sequence))
    SeqIO.write(dummy_list, fasta_dir+ species_id + "ava.fasta", "fasta")
    _ = NcbiblastpCommandline(query= fasta_dir + species_id + "ava.fasta",
                              subject= fasta_dir + species_id + "ava.fasta",
                              out=output_dir + species_id +'ava.tab', outfmt=6, num_threads = 8)()[0]
    #add column names
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
            'bitscore']
    output = pd.read_csv(output_dir + species_id +'ava.tab', names = cols, sep = "\t")
    #lengths of query & target sequences
    q_len_df = pd.DataFrame({'qseqid':ids, 'q_seq_len':lens})
    s_len_df = pd.DataFrame({'sseqid':ids, 's_seq_len':lens})
    #remove self alignments
    output = output[output.qseqid != output.sseqid]
    #add lengths of sequences
    output = output.merge(q_len_df, on = 'qseqid', how = 'left')
    output = output.merge(s_len_df, on ='sseqid', how = 'left')
    #keep only alignments where the query is larger than the subject - we are trying to decide
    #if we should keep the query or not
    output = output[output.q_seq_len >= output.s_seq_len]
    #fitler out same len duplicate pairs
    #add a unqique pair id to the pairs and drop duplicate pair ids
    output['pair_id'] = output.apply(lambda row: '_'.join(sorted([row['qseqid'], row['sseqid']])), axis = 1)
    output = output.drop_duplicates(subset = 'pair_id')
    output.to_csv(output_dir + species_id +'ava.tab', index = False, sep="\t")
    return output


def run_blastp_target_v_all(species_id, orig_fasta_name, ortholog_db, output_dir = './blastp_outputs/', fasta_dir = './single_family_fastas/'):
    #make a temp fasta and align it with itself with blastp
    dummy_list = []
    ids = []
    lens = []

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(fasta_dir):
        os.mkdir(fasta_dir)

    for ortho in ortholog_db[species_id]:
        dummy_list.append(SeqRecord(Seq(ortho.sequence),
                   id=ortho.oma_id))
        ids.append(ortho.oma_id)
        lens.append(len(ortho.sequence))
    SeqIO.write(dummy_list, fasta_dir+ species_id + "tva.fasta", "fasta")

    #running blastp on orig vs the family to see if there are diff overlap regions
    _ = NcbiblastpCommandline(query= orig_fasta_name,
                              subject=fasta_dir + species_id + "tva.fasta",
                              out=output_dir + species_id + '_overlaps.tab', outfmt=6, num_threads = 8)()[0]
    #add columns to results
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
            'bitscore']
    output = pd.read_csv(output_dir + species_id + '_overlaps.tab', names = cols, sep = "\t")
    #drop smaller subalignmnets produced by blastp for the longest alignment possible for each ortholog seq
    output = output.sort_values('length', ascending=False).drop_duplicates('sseqid').sort_index()
    s_len_df = pd.DataFrame({'sseqid':ids, 's_seq_len':lens})
    output = output.merge(s_len_df, on ='sseqid', how = 'left')
    output.to_csv(output_dir + species_id + '_overlaps.tab', index = False, sep="\t")
    return output

def calc_id90(x, blast_results):
    vals = blast_results[blast_results.qseqid == x].pident.to_numpy()
    return np.percentile(vals, 90)

def run_blastp_species_v_singles(species_id, single_family_fasta, ortholog_db, output_dir = './blastp_outputs/', fasta_dir = './single_family_fastas/'):
    #make a temp fasta and align it with itself with blastp
    dummy_list = []
    ids = []
    lens = []
    seqs = []

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(fasta_dir):
        os.mkdir(fasta_dir)

    for ortho in ortholog_db[species_id]:
        dummy_list.append(SeqRecord(Seq(ortho.sequence),
                   id=ortho.oma_id))
        ids.append(ortho.oma_id)
        lens.append(len(ortho.sequence))
        seqs.append(ortho.sequence)
    SeqIO.write(dummy_list, fasta_dir+ species_id + "fvs.fasta", "fasta")

    #running blastp on orig vs the family to see if there are diff overlap regions
    output = \
    NcbiblastpCommandline(query=fasta_dir+ species_id + "fvs.fasta", subject=single_family_fasta, out= output_dir + species_id +'_id90.tab', outfmt=6, num_threads = 8)()[0]

    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
            'bitscore']
    output = pd.read_csv(output_dir + species_id +'_id90.tab', names = cols, sep = "\t")
    #calculate id90 value
    frame_id90 = pd.DataFrame({'oma_id': ids})
    frame_id90['id_90'] = frame_id90.oma_id.apply(lambda x: calc_id90(x, output))
    #sort by best
    frame_id90 = frame_id90.sort_values(by = 'id_90', ascending= False)
    return frame_id90


