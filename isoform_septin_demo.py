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

#now do the same steps as previous - but with the 9 swapped isoforms
print ("SEPTIN 1SDS")
filtered_sept1ds_orthodb = ortholog_database(fasta_dir + s1ds_filt_ortho )
filtered_sept1ds_orthodb.get_unique_species()
print ('Size db: ', filtered_sept1ds_orthodb.size)
#2) filter by genome duplication issues - removing shorter orthologs from a species which align highly to a longer one
fragmented_filtered_dict = genome_fragment_filter(filtered_sept1ds_orthodb)
#save filtered db
save_db_fasta(fragmented_filtered_dict, fasta_dir + 'sept1ds_retain_all.fasta', only_single = False)
#print ("multi-alignment filter - check if it's an issue (if it is, implement this....")
#region_merge_filter(fragmented_filtered_dict, fasta_dir + s1_filt_ortho)
print ("-----------------------------")
#Option 1- MSAs with all species orthologs
all_options_stockholm_msa_generation(fasta_dir + s1ds, fasta_dir + 'sept1ds_retain_all.fasta', msas +'sept1ds_retain_all/')
#filter on 5%, 1% and 0% (just keep top seq90) values for selecting orthologs from each species
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict)
save_db_fasta(filtered_dict, fasta_dir + 'sept1ds_filt_id90_5.fasta')
all_options_stockholm_msa_generation(fasta_dir + s1ds, fasta_dir + 'sept1ds_filt_id90_5.fasta', msas +'sept1ds_filt_id90_5/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 1)
save_db_fasta(filtered_dict, fasta_dir + 'sept1ds_filt_id90_1.fasta')
all_options_stockholm_msa_generation(fasta_dir + s1ds, fasta_dir + 'sept1ds_filt_id90_1.fasta', msas +'sept1ds_filt_id90_1/')
print ("-----------------------------")
filtered_dict = ancient_gene_duplication_check(fragmented_filtered_dict, 0)
save_db_fasta(filtered_dict, fasta_dir + 'sept1ds_filt_id90_all.fasta')
all_options_stockholm_msa_generation(fasta_dir + s1ds, fasta_dir + 'sept1ds_filt_id90_all.fasta', msas + 'sept1ds_filt_id90_all/')
print ("-----------------------------")


#reformat all sto files to a3m files, filter out high gap lines
proteins = ['sept1ds']
for protein in proteins:
    for alignment_folder in ['_retain_all', '_filt_id90_5', '_filt_id90_1', '_filt_id90_all']:
        full_path = msas + protein + alignment_folder + '/'
        reformat_all_ortholog_pipeline_msas(full_path)
        check_if_gaps_in_pipeline_msas(full_path)

#making actual pMSAs now
#septin12 & septin1 complex pmsas
for alignment_folder in ['_retain_all', '_filt_id90_5', '_filt_id90_1', '_filt_id90_all']:
    septin12_folder = msas + 'sept12' + alignment_folder +'/'
    septin1ds_folder =  msas + 'sept1ds' + alignment_folder +'/'
    output_folder = msas + 'sept12_sept1ds_complex' + alignment_folder +'/'
    generate_all_pmsa_pairs(septin12_folder, septin1ds_folder, output_folder)
    print ("------------------")
