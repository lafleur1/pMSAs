#showing how the different fitlering options for pMSA generation result in different sizes/qualities of pMSAs
import numpy as np
import matplotlib.pyplot as plt
from orthologs import *
import seaborn as sns
import pandas as pd
import os

def convert_a3m_to_sto_score_against_query(a3m_file_name, path_to_reformat = "~/hh-suite/scripts/reformat.pl" ):
    #use hhsuite reformat.pl to a3m to sto in a temp file
    #score with esp-apiquid
    #save the scores as a csv
    #return just query vs other alignment scores
    #0) get wt name
    fasta_lines = SeqIO.parse(a3m_file_name, 'fasta')
    first_line = next(fasta_lines)

    id_first = first_line.id
    os.system(path_to_reformat + " a3m sto " + a3m_file_name + " " + "temp.sto")
    os.system("esl-alipid temp.sto > temp_out.txt")
    #header is messed up by # at beginning, fixing it while loading
    pandas_temp = pd.read_csv("temp_out.txt", skiprows= 1, names = ["seqname1", "seqname2", "%id", "nid", "denomid", "%match", "nmatch", "denommatch"], delim_whitespace= True ) #uses variable sized whitespace as the iterator
    pandas_temp = pandas_temp[pandas_temp.seqname1 == id_first]
    #clean up temps
    os.remove("temp.sto")
    os.remove("temp_out.txt")
    #get number gaps

    #gaps_count = []
    #for l in fasta_lines:
    #    gaps_count.append(str(l.seq).count("-")/len(str(l.seq)))
    return pandas_temp

def size_pMSA():
    return None

#0) ID some good test heterocomplexes to try out (Look at those with known structures to compare predicted complexes against)


#1) evaluate what DB is the best (unfiltered, all, 5, or 1)(IE - how many sequences are lost)
proteins = ['sept1', 'sept5', 'sept12', 'sept1ds', 'P0AAV4', 'P75745']
sizes_df = pd.DataFrame({'protein':[], 'type':[], 'size':[], 'spec':[]})
for pro in proteins:
    all_db = current_db = './fastas/' + pro + '_retain_all' + '.fasta'
    ortho_db = ortholog_database(current_db)
    total_size = ortho_db.size
    total_spec = len(ortho_db.get_oma_species())
    for db_type in  ['_filt_id90_5', '_filt_id90_1', '_filt_id90_all']:
        current_db = './fastas/' + pro + db_type + '.fasta'
        ortho_db = ortholog_database(current_db)
        sizes_df.loc[len(sizes_df.index)] = [pro, db_type[1:], ortho_db.size/total_size, len(ortho_db.get_oma_species())/total_spec]
print (sizes_df)
sizes_df.to_csv('sizes_df.csv', index = False)
#sns.scatterplot(data = sizes_df, x = 'size', y = 'spec', hue = 'protein', style = 'type')
#plt.show()
#lose a good portion of the orthologs filtering by species....this is unsuprising but needed to be confirmed
#TODO: verify that dropping them at the DB stage results in different pMSAs then filtering the full DB alignments

#2) Evaluate which msa generation method is the best (full seed, partial seed, or subset seed)
#reformat temp. with reformat.pl then use esl-apid from hhsuite to score & drop irrelevant data - look at score distributions with query
quality_df = pd.DataFrame({'protein':[], 'type':[], 'delta_id_filt_v_unfilt':[], 'delta_match_filt_v_unfilt':[], 'size':[]})
proteins = ['sept1', 'sept5', 'sept12', 'sept1ds', 'P0AAV4', 'P75745']
#izes_df = pd.DataFrame({'protein':[], 'type':[], 'size':[], 'spec':[]})
for pro in proteins:
    for db_type in  ['_retain_all', '_filt_id90_5', '_filt_id90_1', '_filt_id90_all']:
        project_name = './alignments/' + pro + db_type + '/'
        #looking at all 3 start msa type (before filtering)
        pipeone_name = project_name + 'filteredFastaFinalAlign.sto' #starting from a filtered seed
        pipeone_reformat = pipeone_name.replace('.sto', '.a3m')
        pipetwo_name = project_name + 'filteredOrigMSAFinalAlign.sto' #starting from subset seed
        pipetwo_reformat = pipetwo_name.replace('.sto', '.a3m')
        pipethree_name = project_name + 'startMSAAlign.sto' #seed of full orig MSA
        pipethree_reformat = pipethree_name.replace('.sto', '.a3m')
        #full and subset seeds = same final alignments.....
        #open all 3 options and compare
        df_one = convert_a3m_to_sto_score_against_query(pipeone_reformat)[["seqname1", "seqname2", "%id","%match"]]
        df_one['type'] = 'filt_seed'
        #print (df_one)
        #print (df_one.columns)
        df_two = convert_a3m_to_sto_score_against_query(pipetwo_reformat)[["seqname1", "seqname2", "%id","%match"]]
        #print (df_two)
        df_two['type'] = 'subset_seed'
        df_three = convert_a3m_to_sto_score_against_query(pipethree_reformat)[["seqname1", "seqname2", "%id","%match"]]
        df_three['type'] = 'full_seed'
        concat_df = pd.concat([df_one, df_two, df_three])
        #print (concat_df)
        #print (df_three)
        merged_onetwo = df_one.merge(df_two, on = ['seqname1', 'seqname2'], suffixes = ['_1', '_2'])
        merged_onetwothree = merged_onetwo.merge(df_three, on = ['seqname1', 'seqname2'])
        merged_onetwothree['%id_delta_filt_full'] = merged_onetwothree['%id_1'] - merged_onetwothree['%id_2']
        merged_onetwothree['%match_delta_filt_full'] = merged_onetwothree['%match_1'] - merged_onetwothree['%match_2']
        #print (merged_onetwothree[["%id_1","%match_1","%id_2","%match_2","%id","%match"]])
        #print ('starting from a filtered seed: ', df_one['%id'].mean(), df_one['%match'].mean() )
        #print('starting from subset seed: ', df_two['%id'].mean(), df_two['%match'].mean())
        #print('seed of full orig MSA: ', df_three['%id'].mean(), df_three['%match'].mean())
        print ('avg_diff_filt_seed & subset and/or full: ', merged_onetwothree['%id_delta_filt_full'].mean())
        print('avg_diff_filt_seed & subset and/or full: ', merged_onetwothree['%match_delta_filt_full'].mean())
        quality_df.loc[len(quality_df.index)] = [pro, db_type[1:], merged_onetwothree['%id_delta_filt_full'].mean(),
                                             merged_onetwothree['%match_delta_filt_full'].mean(), merged_onetwothree.shape[0]]
        #sns.scatterplot(x = '%id', y = '%match', hue = 'type', data = concat_df)
        #mark which are lost during pre-filtering v post-filtering

        #plt.show()
        print ("-------------------")
print (quality_df)
quality_df.to_csv('quality_df.csv', index = False)

#it appears that using the filtered HMM seed has slight benefits / very small detriments for the examined proteins so far
#Also - there is unsuprisingly no difference between the subset original MSA seed and the full original MSA seed for making the alingnment HMM

#3)
#TODO Evaluate how many lines are lost with different hhfilter citeria
#4)
#TODO Evaluate double vs single info retained with prefitler vs. postfiltering with hhfilter (prefilter = filter individual MSAs THEN pair, postfilter = pair MSAs then filter the pMSA)
