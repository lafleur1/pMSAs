#processing the s1ds isoform complex

#option 1 - swapping out isoforms for a different isoform if its better than the existing one
#option 2 - reblast the isoform against the target proteome

#option 1
from oma_access import *
from orthologs import *
from stockholm import *
from a3m import *
from pairing_by_species import *

fasta_dir = "./fastas/"
msas = "./alignments/"
if not os.path.isdir(fasta_dir):
    os.mkdir(fasta_dir)
if not os.path.isdir(msas):
    os.mkdir(msas)

s12 = "septin12.fasta" #Uniprot Q8IYM1
s12_filt_ortho = 's12_filtered_orthologs.fasta'
s1 = "septin1.fasta" #Uniprot Q8WYJ6
s1_filt_ortho = 's1_filtered_orthologs.fasta'
s1ds = "septin1ds.fasta" #Uniprot J3kNL2
s1_filt_ortho_w_iso = 'all_orthologs_and_isoforms_septin1.fasta'
s1ds_filt_ortho =  's1ds_isoform_filt_orthologs.fasta'

filtered_sept1_orthodb = ortholog_database(fasta_dir + s1_filt_ortho) #created in original septin complex demo
species = filtered_sept1_orthodb.get_oma_species()
#retrieve_all_oma_isoforms(filtered_sept1_orthodb, fasta_dir + s1_filt_ortho_w_iso) #get all annotated isoforms of ortholog seqeuences where relevant

#count diffs in size
orig_filtered = ortholog_database(fasta_dir + s1_filt_ortho)
print ('No isoform orthologs: ', orig_filtered.size)
#load isoform db
all_isoforms_orthologs_septin1 = ortholog_database(fasta_dir + s1_filt_ortho_w_iso)
print ('Orthologs + isoforms: ', all_isoforms_orthologs_septin1.size)
#for every ortholog with a isoform, run pBLAST with the isoform sequence vs original sequence and its isoforms and see if a better match is possible (by evalyue, bit score)
reblast_and_filter_isoforms_fasta(all_isoforms_orthologs_septin1, fasta_dir + s1ds, fasta_dir +s1ds_filt_ortho)

#now do the same steps as previous
