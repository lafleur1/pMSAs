from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from Bio import SeqIO
#from omadb import Client

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

