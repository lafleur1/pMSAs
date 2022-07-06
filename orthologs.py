#need to trim down the number of orthologs from one species - in case it's ambiguous.....
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from Bio import SeqIO
from omadb import Client

from pblast import *

class ortholog_with_isoforms:
    def __init__(self, description, sequence):
        #break up description to useful info
        #OMA_ID | OTHER_ID [species name]
        self.oma_id = re.search('[^|]*', description).group().strip(" ")
        self.has_isoform = False
        self.ref_iso = None
        if len(self.oma_id.split(" ")) > 1:
            #is an isoform
            self.has_isoform = True
            self.ref_iso = self.oma_id.split(" ")[1]
            self.oma_id = self.oma_id.split(" ")[0]
        self.other_id = re.search('(?<=\|).*?(?=\[)', description).group().replace(" ","")
        self.sequence = sequence
        self.species_code = self.oma_id[0:5]
        self.species_desc = re.search('(?<=\[).*?(?=\])', description).group()

    def return_seqrecord(self):
        #returns it as a seq record
        description = self.oma_id + " | " + self.other_id + '[' + self.species_desc + ']'
        return SeqRecord(id = self.oma_id, description = description, seq= Seq(self.sequence))

    def retrieve_isofroms(self):
        c = Client()
        query_info = c.proteins[self.oma_id]
        isoform_query = query_info.isoforms
        isoforms = [] #list of other omaids for the gene
        for iso in isoform_query:
            isoforms.append(iso['omaid'])
        return isoforms

class ortholog_database:
    def __init__(self, fasta_name):
        self.species_other_dict = {}
        for record in SeqIO.parse(fasta_name, 'fasta'):
            ortho_temp = ortholog_with_isoforms(record.description, str(record.seq))
            #assign each with same species to dict list
            if ortho_temp.species_code in self.species_other_dict:
                self.species_other_dict[ortho_temp.species_code].append(ortho_temp)
            else:
                self.species_other_dict[ortho_temp.species_code] = [ortho_temp]
        #if it's a single ortholog, add it to this dict
        self.single_orthologs_per_species_dict = {}
        for key in self.species_other_dict:
            if len(self.species_other_dict[key]) == 1:
                self.single_orthologs_per_species_dict[key] = self.species_other_dict[key]
        self.size = sum([len(self.species_other_dict[x]) for x in self.species_other_dict])
        self.ref_isoforms_list = []#[x.oma_id for x in self.species_other_dict[y] for y in self.species_other_dict if x.has_isoform]
        self.update_reference_isoforms()


    def get_unique_species(self):
        print ('number spec: ', len(self.species_other_dict))
        print ("with 1 ortho: ", len(self.single_orthologs_per_species_dict))


    def get_oma_species(self):
        return list(self.species_other_dict.keys())

    def update_reference_isoforms(self):
        #updates list of all refernece isoforms in the db
        self.ref_isoforms_list = []
        for spec in self.species_other_dict:
            self.ref_isoforms_list = self.ref_isoforms_list + [x.ref_iso for x in self.species_other_dict[spec] if x.has_isoform]
        self.ref_isoforms_list = list(set(self.ref_isoforms_list))

def genome_fragment_filter(filtered_sept1_orthodb):

    #Align multi-cipy orthologs from proteome with BLAST
    #If two are > 95% identical in over 80% of residues in the shorter, the shorter is redundant
    #for more than one do all-by-all pairwise alignment of longest with the shortest ones and see if just the longest should be kept

    species_mult_entries = []

    for key in filtered_sept1_orthodb.species_other_dict:
        if len(filtered_sept1_orthodb.species_other_dict[key]) > 1:
            species_mult_entries.append(key)

    #print ('spec with mult counts: ', len(species_mult_entries))
    to_discard = {}
    filtered_dict = filtered_sept1_orthodb.species_other_dict.copy()
    # for these, make a little fasta and run blastp all v all

    print('start filtering....')
    for mult_key in species_mult_entries:
        #print (mult_key)
        #print ( [x.oma_id for x in filtered_sept1_orthodb.species_other_dict[mult_key]])
        outputs = run_blastp_all_v_all(mult_key, filtered_sept1_orthodb.species_other_dict)
        if outputs[(outputs.pident >= 95) & (outputs.length >= 0.8 * outputs.s_seq_len)].shape[0] > 0:
            # print (outputs[(outputs.pident >= 95) & (outputs.length >= 0.8 * outputs.s_seq_len)]) #shorter is redundant - add strip list
            for ind, row in outputs[(outputs.pident >= 95) & (outputs.length >= 0.8 * outputs.s_seq_len)].iterrows():
                if mult_key in to_discard:
                    to_discard[mult_key].append(row['sseqid'])
                else:
                    to_discard[mult_key] = [row['sseqid']]

    print('species with duplicates: ', len(to_discard))

    # remove them and return the dictionary
    num_items = sum(
        [len(filtered_sept1_orthodb.species_other_dict[key]) for key in filtered_sept1_orthodb.species_other_dict])
    print('num items orig: ', num_items)

    for mult_key in to_discard:
        origs = filtered_sept1_orthodb.species_other_dict[mult_key]
        new_list = []
        for thing in origs:
            if thing.oma_id not in to_discard[mult_key]:
                new_list.append(thing)
        filtered_dict[mult_key] = new_list

    num_items = sum([len(filtered_dict[key]) for key in filtered_dict])
    print('filtered: ', num_items)
    return filtered_dict

def calc_overlap_two_orthologs(start1, end1, start2, end2):
    #figure out value of overlap between longer (start1, end1) and shorter (start2, end2) proteins in the family
    if start2 >= start1 and start2 <= end1:
        return end1 - start2 + 1
    elif end2 >= start1 and end2 <= end1:
        return end2 - start1 + 1

def region_merge_filter(spec_dict, target_fasta_name):
    #seeing if proteins can be merged into one
    #looking to see if there are diff regions of overlap to orig for the family
    #do mini blastp with orig and the family only
    species_mult_entries = []
    for key in spec_dict:
        if len(spec_dict[key]) > 1:
            species_mult_entries.append(key)
    print('spec with mult counts: ', len(species_mult_entries))

    #make indiv blastp
    print('start filtering....')
    merged_pairs = {}
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
            'bitscore']
    for i in range(len(species_mult_entries)):
        mult_key = species_mult_entries[i]
        if i % 25 == 0:
            print(i, 'out of ', len(species_mult_entries))
        outputs = run_blastp_target_v_all(mult_key, target_fasta_name, spec_dict)
        #figure out if length overlap < 25% of longer protein of the two alignments to the target
        #output tables are ordered by length - should only have to look at following row pairs to see if they should be merged
        for i in range(0, outputs.shape[0]):
            for j in range(i, outputs.shape[0]):
                #compare row i to row j
                rowj = outputs.iloc[j]
                rowi = outputs.iloc[i]
                overlap_current = calc_overlap_two_orthologs(int(rowj['qstart']), int(rowj['qend']), int(rowi['qstart']), int(rowj['qend']))
                if overlap_current <= 0.25 * int(rowj['length']):
                    if mult_key in merged_pairs:
                        merged_pairs[mult_key].append((rowj['sseqid'], rowi['sseqid']))
                    else:
                        merged_pairs[mult_key] = [(rowj['sseqid'], rowi['sseqid'])]
    #merge proteins into one
    print ("Orthologs to merge: ", len(merged_pairs))

def save_db_fasta(db_dict, savename, only_single = True):
    #generate a fasta of orthologs per species
    #use seqrecord option to maintain all species info
    dummy_list = []
    for family in db_dict:
        if only_single: # == 'single':
            if len(db_dict[family]) == 1:
                for ortho in db_dict[family]:
                    dummy_list.append(ortho.return_seqrecord())
        else:
            for ortho in db_dict[family]:
                dummy_list.append(ortho.return_seqrecord())
    SeqIO.write(dummy_list, savename, "fasta")
    print ('number single copy orthologs saved: ', len(dummy_list))

def make_single_copy_fasta(db_dict, savename = './single_family_fastas/single_copy_db.fasta'):
    #generate a fasta of every single copy HOG sequence for ortholog selection per species
    dummy_list = []
    ids = []
    lens = []

    for family in db_dict:
        if len(db_dict[family]) == 1:
            for ortho in db_dict[family]:
                dummy_list.append(SeqRecord(Seq(ortho.sequence),
                                            id=ortho.oma_id))
                ids.append(ortho.oma_id)
                lens.append(len(ortho.sequence))
    SeqIO.write(dummy_list, savename, "fasta")

def ancient_gene_duplication_check(filtered_sept1_orthodb, cutoff_val = 5):
    #figure out highest id90 value to the other single-copyy orthologs
    #make single copy fasta db for blasting against with a given species family
    make_single_copy_fasta(filtered_sept1_orthodb)

    #find each multi-ex species
    species_mult_entries = []
    for key in filtered_sept1_orthodb:
        if len(filtered_sept1_orthodb[key]) > 1:
            species_mult_entries.append(key)
    print('spec with mult counts: ', len(species_mult_entries))

    # make indiv blastp
    print('start filtering....')
    select_one = {}
    bad_groupings = {}
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
            'bitscore']
    for i in range(len(species_mult_entries)):
        mult_key = species_mult_entries[i]
        if i % 25 == 0:
            print (i, 'out of ', len(species_mult_entries))
        #run_blastp_species_v_singles(species_id, single_family_fasta, ortholog_db)
        outputs = run_blastp_species_v_singles(mult_key, './single_family_fastas/single_copy_db.fasta', filtered_sept1_orthodb)
        first_val = outputs.iloc[0].id_90
        outputs['adjusted_90'] = outputs['id_90'] - first_val
        others = outputs.iloc[1:]
        #print (others[others.adjusted_90  > -5], others[others.adjusted_90  > -5].shape)
        if others[others.adjusted_90  > -1 * cutoff_val].shape[0] != 0:
            #not selection one....
            bad_groupings[mult_key] = outputs
            #print(others.iloc[0].adjusted_90)
        else:
            #print (others)
            select_one[mult_key] = outputs.iloc[0].oma_id
        #print ("------------------------")
    print (len(select_one), len(bad_groupings))

    #filter down multispecies
    filtered_dict = filtered_sept1_orthodb.copy()
    for mult_key in select_one:
        origs = filtered_sept1_orthodb[mult_key]
        new_list = []
        for thing in origs:
            if thing.oma_id in select_one[mult_key]:
                new_list.append(thing)
        filtered_dict[mult_key] = new_list

    num_items = sum([len(filtered_dict[key]) for key in filtered_dict])
    print('filtered: ', num_items)
    return filtered_dict


