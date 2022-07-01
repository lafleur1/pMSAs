#I'm interesting in predicting complexes for S12-S1, S12-S1DS, S12-S5
#S1DS is a splice variant of S1 with an additional exon

from pblast import *
from oma_access import *
from orthologs import *

# Basic info

#Looking at the sequence identity they have with each other

fasta_dir = "./fastas/"

s12 = "septin12.fasta" #Uniprot Q8IYM1
s12_filt_ortho = 's12_filtered_orthologs.fasta'
s1 = "septin1.fasta" #Uniprot Q8WYJ6
s1_filt_ortho = 's1_filtered_orthologs.fasta'
s1ds = "septin1ds.fasta" #Uniprot J3kNL2
s1ds_filt_ortho = 's1ds_filtered_orthologs.fasta'
s5 = "septin5.fasta" #Uniprot Q99719
s5_filt_ortho = 's5_filtered_orthologs.fasta'
all_septins = 'all_septin_fastas.fasta' #all wt above in one fasta for all v all blasting

run_pblast_self_v_self(fasta_dir + all_septins, 'septin_v_septin_blast', output_dir = './blastp_outputs/')
"""

% Seq ID longest pBLAST alignment 

        S1 | S1DS | S5 | S12
| S1   | -   | 100 | 63 | 45 |
| S1DS | 100 | -   | 63 | 45 |
| S5   | 64  | 64  | -  | 46 |
| S12  | 45  | 45  | 46 | -  |
"""

#filter OMA orthologs and paralogs
print ("SEPTIN 1")
filter_for_paralogs(fasta_dir + 's1_orthologs.txt', fasta_dir + 's1_paralogs.txt', s1_filt_ortho)
print ("SEPTIN 5")
filter_for_paralogs(fasta_dir + 's5_orthologs.txt', fasta_dir + 's5_paralogs.txt', s5_filt_ortho)
print ("SEPTIN 12")
filter_for_paralogs(fasta_dir + 's12_orthologs.txt', fasta_dir + 's12_paralogs.txt', s12_filt_ortho)

"""
       |# orthologs | # paralogs| # filtered orthologs |
| S1   | 1003   | 2123 | 1003 | 
| S5   | 1143  | 2060  | 1066 | 
| S12  | 1034  | 2126  | 995  :w|
"""

#For each filtered db, select orthologs to use for alignment

#s1
#SEPTIN 1
print ("SEPTIN 1")
print ("Genome fragment filter....")
filtered_sept1_orthodb = ortholog_database(fasta_dir + s1_filt_ortho )
filtered_sept1_orthodb.get_unique_species()
print ('Size db: ', filtered_sept1_orthodb.size)
#2) filter by genome duplication issues - removing shorter orthologs from a species which align highly to a longer one
fragmented_filtered_dict = genome_fragment_filter(filtered_sept1_orthodb)
print ("multi-alignment filter - check if it's an issue (if it is, implement this....")
#region_merge_filter(fragmented_filtered_dict, fasta_dir + s1_filt_ortho)
#filter on 5%, 1% and 0% (just keep top seq90) values for selecting orthologs from each species
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict)
save_db_fasta(filtered_dict, fasta_dir + 'sept1_filt_id90_5.fasta')
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 1)
save_db_fasta(filtered_dict, fasta_dir + 'sept1_filt_id90_1.fasta')
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 0)
save_db_fasta(filtered_dict, fasta_dir + 'sept1_filt_id90_all.fasta')


print ("SEPTIN 5")
print ("Genome fragment filter....")
filtered_sept5_orthodb = ortholog_database(fasta_dir + s5_filt_ortho )
filtered_sept5_orthodb.get_unique_species()
print ('Size db: ', filtered_sept5_orthodb.size)
#2) filter by genome duplication issues - removing shorter orthologs from a species which align highly to a longer one
fragmented_filtered_dict = genome_fragment_filter(filtered_sept5_orthodb)
print ("multi-alignment filter - check if it's an issue (if it is, implement this....")
region_merge_filter(fragmented_filtered_dict, fasta_dir + s5_filt_ortho)
#filter on 5%, 1% and 0% (just keep top seq90) values for selecting orthologs from each species
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict)
save_db_fasta(filtered_dict, 'sept5_filt_id90_5.fasta')
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 1)
save_db_fasta(filtered_dict, 'sept5_filt_id90_1.fasta')
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 0)
save_db_fasta(filtered_dict, 'sept5_filt_id90_all.fasta')


print ("SEPTIN 12")
print ("Genome fragment filter....")
filtered_sept12_orthodb = ortholog_database(fasta_dir + s12_filt_ortho )
filtered_sept12_orthodb.get_unique_species()
print ('Size db: ', filtered_sept12_orthodb.size)
#2) filter by genome duplication issues - removing shorter orthologs from a species which align highly to a longer one
fragmented_filtered_dict = genome_fragment_filter(filtered_sept12_orthodb)
print ("multi-alignment filter - check if it's an issue (if it is, implement this....")
region_merge_filter(fragmented_filtered_dict, fasta_dir + s12_filt_ortho)
#filter on 5%, 1% and 0% (just keep top seq90) values for selecting orthologs from each species
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict)
save_db_fasta(filtered_dict, 'sept12_filt_id90_5.fasta')
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 1)
save_db_fasta(filtered_dict, 'sept12_filt_id90_1.fasta')
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 0)
save_db_fasta(filtered_dict, 'sept12_filt_id90_all.fasta')

#sizes of filtered databases
