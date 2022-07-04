#I'm interesting in predicting complexes for S12-S1, S12-S1DS, S12-S5
#S1DS is a splice variant of S1 with an additional exon

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

#Looking at the sequence identity they have with each other


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
print ("-----------------------------")
print ("SEPTIN 5")
filter_for_paralogs(fasta_dir + 's5_orthologs.txt', fasta_dir + 's5_paralogs.txt', s5_filt_ortho)
print ("-----------------------------")
print ("SEPTIN 12")
filter_for_paralogs(fasta_dir + 's12_orthologs.txt', fasta_dir + 's12_paralogs.txt', s12_filt_ortho)
print ("-----------------------------")

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
filtered_sept1_orthodb = ortholog_database(fasta_dir + s1_filt_ortho )
filtered_sept1_orthodb.get_unique_species()
print ('Size db: ', filtered_sept1_orthodb.size)
#2) filter by genome duplication issues - removing shorter orthologs from a species which align highly to a longer one
fragmented_filtered_dict = genome_fragment_filter(filtered_sept1_orthodb)
#save filtered db
save_db_fasta(fragmented_filtered_dict, fasta_dir + 'sept1_retain_all.fasta', only_single = False)
#print ("multi-alignment filter - check if it's an issue (if it is, implement this....")
#region_merge_filter(fragmented_filtered_dict, fasta_dir + s1_filt_ortho)
print ("-----------------------------")
#Option 1- MSAs with all species orthologs
all_options_stockholm_msa_generation(fasta_dir + s1, fasta_dir + 'sept1_retain_all.fasta', msas +'sept1_retain_all/')
#filter on 5%, 1% and 0% (just keep top seq90) values for selecting orthologs from each species
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict)
save_db_fasta(filtered_dict, fasta_dir + 'sept1_filt_id90_5.fasta')
all_options_stockholm_msa_generation(fasta_dir + s1, fasta_dir + 'sept1_filt_id90_5.fasta', msas +'sept1_filt_id90_5/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 1)
save_db_fasta(filtered_dict, fasta_dir + 'sept1_filt_id90_1.fasta')
all_options_stockholm_msa_generation(fasta_dir + s1, fasta_dir + 'sept1_filt_id90_1.fasta', msas +'sept1_filt_id90_1/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 0)
save_db_fasta(filtered_dict, fasta_dir + 'sept1_filt_id90_all.fasta')
all_options_stockholm_msa_generation(fasta_dir + s1, fasta_dir + 'sept1_filt_id90_all.fasta', msas + 'sept1_filt_id90_all/')
print ("-----------------------------")


print ("SEPTIN 5")
filtered_sept5_orthodb = ortholog_database(fasta_dir + s5_filt_ortho )
filtered_sept5_orthodb.get_unique_species()
print ('Size db: ', filtered_sept5_orthodb.size)
#2) filter by genome duplication issues - removing shorter orthologs from a species which align highly to a longer one
fragmented_filtered_dict = genome_fragment_filter(filtered_sept5_orthodb)
#save filtered db
save_db_fasta(fragmented_filtered_dict, fasta_dir + 'sept5_retain_all.fasta', only_single = False)
#print ("multi-alignment filter - check if it's an issue (if it is, implement this....")
#region_merge_filter(fragmented_filtered_dict, fasta_dir + s5_filt_ortho)
print ("-----------------------------")
#Option 1- MSAs with all species orthologs
all_options_stockholm_msa_generation(fasta_dir + s5, fasta_dir + 'sept5_retain_all.fasta', msas +'sept5_retain_all/')
#filter on 5%, 1% and 0% (just keep top seq90) values for selecting orthologs from each species
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict)
save_db_fasta(filtered_dict, fasta_dir + 'sept5_filt_id90_5.fasta')
all_options_stockholm_msa_generation(fasta_dir + s5, fasta_dir + 'sept5_filt_id90_5.fasta', msas +'sept5_filt_id90_5/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 1)
save_db_fasta(filtered_dict, fasta_dir + 'sept5_filt_id90_1.fasta')
all_options_stockholm_msa_generation(fasta_dir + s5, fasta_dir + 'sept5_filt_id90_1.fasta', msas +'sept5_filt_id90_1/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 0)
save_db_fasta(filtered_dict, fasta_dir + 'sept5_filt_id90_all.fasta')
all_options_stockholm_msa_generation(fasta_dir + s5, fasta_dir + 'sept5_filt_id90_all.fasta', msas +'sept5_filt_id90_all/')
print ("-----------------------------")



print ("SEPTIN 12")
print ("Genome fragment filter....")
filtered_sept12_orthodb = ortholog_database(fasta_dir + s12_filt_ortho )
filtered_sept12_orthodb.get_unique_species()
print ('Size db: ', filtered_sept12_orthodb.size)
#2) filter by genome duplication issues - removing shorter orthologs from a species which align highly to a longer one
fragmented_filtered_dict = genome_fragment_filter(filtered_sept12_orthodb)
save_db_fasta(fragmented_filtered_dict, fasta_dir + 'sept12_retain_all.fasta', only_single = False)
all_options_stockholm_msa_generation(fasta_dir + s12, fasta_dir + 'sept12_retain_all.fasta', msas +'sept12_retain_all/')

#print ("multi-alignment filter - check if it's an issue (if it is, implement this....")
#region_merge_filter(fragmented_filtered_dict, fasta_dir + s12_filt_ortho)
print ("-----------------------------")
#filter on 5%, 1% and 0% (just keep top seq90) values for selecting orthologs from each species
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict)
save_db_fasta(filtered_dict, fasta_dir + 'sept12_filt_id90_5.fasta')
all_options_stockholm_msa_generation(fasta_dir + s12, fasta_dir + 'sept12_filt_id90_5.fasta', msas +'sept12_filt_id90_5/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 1)
save_db_fasta(filtered_dict, fasta_dir + 'sept12_filt_id90_1.fasta')
all_options_stockholm_msa_generation(fasta_dir + s12, fasta_dir + 'sept12_filt_id90_1.fasta', msas +'sept12_filt_id90_1/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 0)
save_db_fasta(filtered_dict, fasta_dir + 'sept12_filt_id90_all.fasta')
all_options_stockholm_msa_generation(fasta_dir + s12, fasta_dir + 'sept12_filt_id90_all.fasta', msas +'sept12_filt_id90_all/')
print ("-----------------------------")



#reformat all sto files to a3m files, filter out high gap lines
proteins = ['sept1', 'sept5', 'sept12']
for protein in proteins:
    for alignment_folder in ['_retain_all', '_filt_id90_5', '_filt_id90_1', '_filt_id90_all']:
        full_path = msas + protein + alignment_folder + '/'
        reformat_all_ortholog_pipeline_msas(full_path)
        check_if_gaps_in_pipeline_msas(full_path)

#making actual pMSAs now
#septin12 & septin1 complex pmsas
for alignment_folder in ['_retain_all', '_filt_id90_5', '_filt_id90_1', '_filt_id90_all']:
    septin12_folder = msas + 'sept12' + alignment_folder +'/'
    septin1_folder =  msas + 'sept1' + alignment_folder +'/'
    septin5_folder = msas + 'sept5' + alignment_folder +'/'
    output_folder = msas + 'sept12_sept1_complex' + alignment_folder +'/'
    generate_all_pmsa_pairs(septin12_folder, septin1_folder, output_folder)
    print ("------------------")
    output_folder = msas + 'sept12_sept5_complex' + alignment_folder +'/'
    generate_all_pmsa_pairs(septin12_folder, septin5_folder, output_folder)
    print("------------------")
