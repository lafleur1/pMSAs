from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from Bio import SeqIO
from omadb import Client
from pblast import *

def filter_for_paralogs(ortholog_fasta, paralog_fasta, new_fasta_name, output_dir = "./fastas/"):
    seqs = []
    ids = []
    seqs_ortho = SeqIO.parse(ortholog_fasta, "fasta")
    for record in seqs_ortho:
        seqs.append(record)
        ids.append(str(record.id))
    seqs_para = []
    ids_para = []
    seqs_p = SeqIO.parse(paralog_fasta, "fasta")
    for record in seqs_p:
        seqs_para.append(record)
        ids_para.append(str(record.id))
    # look at overlap
    print('orthologs: ', len(ids), 'unique orthologs: ', len(set(ids)))
    print('paralogs ', len(ids_para), len(set(ids_para)))
    # no overlap if equal to orig
    print('Unique orthologs overlap: ', len(set(ids).difference(set(ids_para))))

    seq_record_list = []
    for i in range(0, len(ids)):
        if ids[i] not in ids_para:
            seq_record_list.append(seqs[i])
    SeqIO.write(seq_record_list, output_dir + new_fasta_name, "fasta")

#for isoform specific ortholog filtering
def retrieve_all_oma_isoforms(ortholog_db, name_orthology_fasta):
    c = Client()
    species = ortholog_db.get_oma_species()
    seqrecord_list = []
    for i in range(0, len(species)):
        # df_one = next(gen_df)
        print (species[i], len(ortholog_db.species_other_dict[species[i]]))
        for ortho in ortholog_db.species_other_dict[species[i]]:
            print ("    ", ortho.oma_id)
            retrieve_seq_ids = ortho.retrieve_isofroms()
            #if isoform list is just 1
            if len(retrieve_seq_ids) == 1:
                seqrecord_list.append(ortho.return_seqrecord())
            else:
                #if it has isoforms, retrieve them all from oma to see if they match the current isoform better than the others
                print("    num isoforms: ", len(retrieve_seq_ids))
                retrieve_seqs = c.proteins[retrieve_seq_ids]
                for entry in retrieve_seqs:
                    id_entry = entry['omaid']
                    fasta_desc = retrieve_seqs[0]['omaid'] + ' | ' + retrieve_seqs[0]['canonicalid'] + ' [' + \
                                 retrieve_seqs[0]['species']['species'] + ']'
                    sequence = entry['sequence']
                    seqrecord_list.append(SeqRecord(id=id_entry, description=fasta_desc, seq=Seq(sequence)))
    SeqIO.write(seqrecord_list, name_orthology_fasta, "fasta")


def reblast_and_filter_isoforms_fasta(full_isoform_db, orig_fasta_name, new_fasta_db_name):
    #reblast and filter each family against the isoform - if one of the isoforms is higher than the orig, keep that instead
    fasta_lines = []
    species = full_isoform_db.get_oma_species()
    full_isoform_db.update_reference_isoforms()
    select_new_isoforms = 0
    num_isoforms2 = 0
    replaced_isoforms = []
    for spec in species:
        #is isoform if they have same cannonical id for the othorlogs
        iso_organizer2 = {}
        #print (spec)
        for ortho in full_isoform_db.species_other_dict[spec]:
            if ortho.has_isoform:
                num_isoforms2 += 1
                if ortho.ref_iso in iso_organizer2:
                    iso_organizer2[ortho.ref_iso].append(ortho.oma_id)
                else:
                    iso_organizer2[ortho.ref_iso] = [ortho.oma_id]
            elif ortho.oma_id not in full_isoform_db.ref_isoforms_list:
                fasta_lines.append(ortho.return_seqrecord())
            #elif ortho.oma_id in full_isoform_db.ref_isoforms_list:
                #print ("ref iso? ")
                #print (ortho.oma_id)

        for ref_isoform in iso_organizer2:
            ids_list = [ref_isoform] + iso_organizer2[ref_isoform] #all isoforms ot get from the species list
            orthos_to_blast = [ortho for ortho in full_isoform_db.species_other_dict[spec] if ortho.oma_id in ids_list]
            stupid_lookup = {}
            for orth in orthos_to_blast:
                stupid_lookup[orth.oma_id] = orth
            #make a fasta to blast against
            blast_results = run_blastp_target_v_isoforms(orthos_to_blast, orig_fasta_name, full_isoform_db)
            #making sure it's not a fake isoform / annotation error isoform in oma
            if (ref_isoform != blast_results.iloc[0].sseqid) and (len(blast_results.s_seq.unique()) == blast_results.shape[0]):
                select_new_isoforms +=1
                print (ref_isoform)
                print (blast_results[['sseqid', 'pident', 'length', 's_seq_len', 'evalue', 'bitscore']])
                #print (blast_results.evalue.to_list())
                print ('-----------------------')
                #replace reference protein
                fasta_lines.append(stupid_lookup[blast_results.iloc[0].sseqid].return_seqrecord())
                replaced_isoforms.append(stupid_lookup[blast_results.iloc[0].sseqid].return_seqrecord())
            else:
                fasta_lines.append(stupid_lookup[ref_isoform].return_seqrecord())

    #make new fastadb with the replaced isoforms
    SeqIO.write(fasta_lines, new_fasta_db_name, "fasta")
    SeqIO.write(replaced_isoforms, new_fasta_db_name.replace('.fasta', '_replacement_isoforms.fasta'), "fasta")

    print ('Total proteins: ', full_isoform_db.size, " pros with isoforms: ", num_isoforms2, " difference isoform selected for db ", select_new_isoforms)


def refilter_for_paralogs(isoform_swapped_fasta, paralog_fasta = "HUMAN27832_paralogs_s1.txt"):
    #make sure that the selected new isoform isn't in the paralog file (shouldn't be, since origianls weren't paralogous to query but I'm paranoid)
    seqs = []
    ids = []
    seqs_ortho = SeqIO.parse(isoform_swapped_fasta, "fasta")
    for record in seqs_ortho:
        seqs.append(record)
        ids.append(str(record.id))
    seqs_para = []
    ids_para = []
    seqs_p = SeqIO.parse(paralog_fasta, "fasta")
    for record in seqs_p:
        seqs_para.append(record)
        ids_para.append(str(record.id))

    if len(set(ids)) != len(set(ids).difference(set(ids_para))):
        print ("warning! reintroduced paralogs")

#for option 2 of isoform specific paired MSA generation, the entire OMA proteome needs to be downloaded and parsed into constituent proteomes
def get_all_oma_proteomes():
    #7/4/22 - this link works, but probably need to update it with oma
    #zipped = 4GB, unzipped 8GB
    os.system("wget https://omabrowser.org/All/oma-seqs.fa.gz")
    os.system("gunzip oma-seqs.fa.gz")
    os.remove("oma-seqs.fa.gz")


def break_up_into_proteomes():
    #assumes you have oma-seqs.fa unzipped in this directory
    #the master proteome list is organized by species id

    proteo_dir = "./oma_proteomes/"
    if not os.path.isdir(proteo_dir):
        os.mkdir(proteo_dir)

    current_records = []
    fasta_sequences = SeqIO.parse(open("oma-seqs.fa"), 'fasta')
    #parse first sequence
    curr_seq = next(fasta_sequences)
    current_spec_id = curr_seq.id[0:5]
    current_file_name = proteo_dir + current_spec_id +'_oma_proteome.fasta'
    current_records.append(curr_seq)
    for fasta in fasta_sequences:
        #for all the other lines
        spec_id = fasta.id[0:5]
        if spec_id == current_spec_id:
            #if the species is still the same, append the new line
            current_records.append(fasta)
        else:
            #new species, reset to make a new fasta file
            #save old records
            print ('writing out genome for ', currnet_spec_id, ' with ', len(current_records), ' proteins ')
            SeqIO.write(current_records, current_file_name, 'fasta')
            #set up new genome file
            current_records = []
            currnet_spec_id = spec_id
            current_file_name = proteo_dir + current_spec_id + '_oma_proteome.fasta'
            current_records.append(fasta)
            #
    #write out last proteome
    print('writing out genome for ', currnet_spec_id, ' with ', len(current_records), ' proteins ')
    SeqIO.write(current_records, current_file_name, 'fasta')

