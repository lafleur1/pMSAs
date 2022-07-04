#tryign to make pmsas with omadb for the same example complex they have on the trRosetta github

#p1 (P0AAV4 · PXPB_ECOLI)
#MQRARCYLIGETAVVLELEPPVTLASQKRIWRLAQRLVDMPNVVEAIPGMNNITVILRNPESLALDAIERLQRWWEESEALEPESRFIEIPVVYGGAGGPDLAVVAAHCGLSEKQVVELHSSVEYVVWFLGFQPGFPYLGSLPEQLHTPRRAEPRLLVPAGSVGIGGPQTGVYPLATPGGWQLIGHTSLSLFDPARDEPILLRPGDSVRFVPQKEGVC
#p2 (P75745 · PXPC_ECOLI)
#MLKIIRAGMYTTVQDGGRHGFRQSGISHCGALDMPALRIANLLVGNDANAPALEITLGQLTVEFETDGWFALTGAGCEARLDDNAVWTGWRLPMKAGQRLTLKRPQHGMRSYLAVAGGIDVPPVMGSCSTDLKVGIGGLEGRLLKDGDRLPIGKSKRDSMEAQGVKQLLWGNRIRALPGPEYHEFDRASQDAFWRSPWQLSSQSNRMGYRLQGQILKRTTDRELLSHGLLPGVVQVPHNGQPIVLMNDAQTTGGYPRIACIIEADMYHLAQIPLGQPIHFVQCSLEEALKARQDQQRYFEQLAWRLHNEN

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

P0AAV4 = "P0AAV4.fasta" #Uniprot Q8IYM1
P0AAV4_filt_ortho = 'P0AAV4_filtered_orthologs.fasta'
P75745 = "P75745.fasta" #Uniprot Q8WYJ6
P75745_filt_ortho = 'P75745_filtered_orthologs.fasta'

print ("P0AAV4 (1)")
filter_for_paralogs(fasta_dir + 'P0AAV4_orthologs.txt', fasta_dir + 'P0AAV4_paralogs.txt', P0AAV4_filt_ortho)
print ("P75745 (2)")
filter_for_paralogs(fasta_dir + 'P75745_orthologs.txt', fasta_dir + 'P75745_paralogs.txt', P75745_filt_ortho)

print ("P0AAV4 (1)")
filtered_P0AAV4_orthodb = ortholog_database(fasta_dir + P0AAV4_filt_ortho )
filtered_P0AAV4_orthodb.get_unique_species()
print ('Size db: ', filtered_P0AAV4_orthodb.size)
#2) filter by genome duplication issues - removing shorter orthologs from a species which align highly to a longer one
fragmented_filtered_dict = genome_fragment_filter(filtered_P0AAV4_orthodb)
#save filtered db
save_db_fasta(fragmented_filtered_dict, fasta_dir + 'P0AAV4_retain_all.fasta')
#print ("multi-alignment filter - check if it's an issue (if it is, implement this....")
#region_merge_filter(fragmented_filtered_dict, fasta_dir + P0AAV4_filt_ortho)
print ("-----------------------------")
#Option 1- MSAs with all species orthologs
all_options_stockholm_msa_generation(fasta_dir + P0AAV4, fasta_dir + 'P0AAV4_retain_all.fasta', msas +'P0AAV4_retain_all/')
#filter on 5%, 1% and 0% (just keep top seq90) values for selecting orthologs from each species
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict)
save_db_fasta(filtered_dict, fasta_dir + 'P0AAV4_filt_id90_5.fasta')
all_options_stockholm_msa_generation(fasta_dir + P0AAV4, fasta_dir + 'P0AAV4_filt_id90_5.fasta', msas +'P0AAV4_filt_id90_5/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 1)
save_db_fasta(filtered_dict, fasta_dir + 'P0AAV4_filt_id90_1.fasta')
all_options_stockholm_msa_generation(fasta_dir + P0AAV4, fasta_dir + 'P0AAV4_filt_id90_1.fasta', msas +'P0AAV4_filt_id90_1/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 0)
save_db_fasta(filtered_dict, fasta_dir + 'P0AAV4_filt_id90_all.fasta')
all_options_stockholm_msa_generation(fasta_dir + P0AAV4, fasta_dir + 'P0AAV4_filt_id90_all.fasta', msas + 'P0AAV4_filt_id90_all/')
print ("-----------------------------")

print ("P75745 (2)")
filtered_P75745_orthodb = ortholog_database(fasta_dir + P75745_filt_ortho )
filtered_P75745_orthodb.get_unique_species()
print ('Size db: ', filtered_P75745_orthodb.size)
#2) filter by genome duplication issues - removing shorter orthologs from a species which align highly to a longer one
fragmented_filtered_dict = genome_fragment_filter(filtered_P75745_orthodb)
#save filtered db
save_db_fasta(fragmented_filtered_dict, fasta_dir + 'P75745_retain_all.fasta')
#print ("multi-alignment filter - check if it's an issue (if it is, implement this....")
#region_merge_filter(fragmented_filtered_dict, fasta_dir + P75745_filt_ortho)
print ("-----------------------------")
#Option 1- MSAs with all species orthologs
all_options_stockholm_msa_generation(fasta_dir + P75745, fasta_dir + 'P75745_retain_all.fasta', msas +'P75745_retain_all/')
#filter on 5%, 1% and 0% (just keep top seq90) values for selecting orthologs from each species
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict)
save_db_fasta(filtered_dict, fasta_dir + 'P75745_filt_id90_5.fasta')
all_options_stockholm_msa_generation(fasta_dir + P75745, fasta_dir + 'P75745_filt_id90_5.fasta', msas +'P75745_filt_id90_5/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 1)
save_db_fasta(filtered_dict, fasta_dir + 'P75745_filt_id90_1.fasta')
all_options_stockholm_msa_generation(fasta_dir + P75745, fasta_dir + 'P75745_filt_id90_1.fasta', msas +'P75745_filt_id90_1/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 0)
save_db_fasta(filtered_dict, fasta_dir + 'P75745_filt_id90_all.fasta')
all_options_stockholm_msa_generation(fasta_dir + P75745, fasta_dir + 'P75745_filt_id90_all.fasta', msas + 'P75745_filt_id90_all/')
print ("-----------------------------")

#reformat all sto files to a3m files, filter out high gap lines
proteins = ['P0AAV4', 'P75745']
for protein in proteins:
    for alignment_folder in ['_retain_all', '_filt_id90_5', '_filt_id90_1', '_filt_id90_all']:
        full_path = msas + protein + alignment_folder + '/'
        reformat_all_ortholog_pipeline_msas(full_path)
        check_if_gaps_in_pipeline_msas(full_path)


for alignment_folder in ['_retain_all', '_filt_id90_5', '_filt_id90_1', '_filt_id90_all']:
    septin12_folder = msas + 'P0AAV4' + alignment_folder +'/'
    septin1_folder =  msas + 'P75745' + alignment_folder +'/'
    output_folder = msas + 'P0AAV4_P75745_complex' + alignment_folder +'/'
    generate_all_pmsa_pairs(septin12_folder, septin1_folder, output_folder)
    print ("------------------")
