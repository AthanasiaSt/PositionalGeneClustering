#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 18:17:02 2022

@author: Athanasia
"""

#libraries used in the code below 

import pandas as pd
import numpy as np
import re 
import pickle 
from itertools import combinations
import seaborn as sns
import ast
import matplotlib.pyplot as plt


#---------------------------------------------------------------------------------------------------------------------------------------------------

# Converting a gtf file to bed-like file, extracting all transcript coordinates 

# These will serve as the background gene 'assembly' on which permutation tests will be applied// 
# this should work for every organism as long as the gtf file is of the right format
# ensGene.gtf files with specific columns // Otherwise you should modify it or provide a bed-like file with chr start end gene_name columns
# you should choose what kind of coordinates you are using for each gene -- here we use transcript coordinates

#THE USER PROVIDES THIS
path_to_gtf = "/home/labuser/Athanasia/master_thesis_files/TFs_sacc/various_data/sacCer2.ensGene.gtf"
 
def gtf_to_bed(path_to_gtf):
    
    # All S cerevisiae genes in a dataframe. The start and end of transcript were kept. The assembly is saccer2.
    # columns : 0-> chromosome , 1-> ensgene, 2-> part of gene, 3-> start of part, 4-> end of part, 8-> ids 
    matrix_genes_sac = pd.read_csv(path_to_gtf,sep='\t', header=None )
    matrix_genes_sac_transcript = matrix_genes_sac[matrix_genes_sac[2] == 'transcript']
    matrix_genes_sac_transcript.dropna(subset=[8], inplace=True) #make sure that all transcripts have a name in the corresponding column 
    
    #finding the gene names that will be put in the same order as in the dataframe --> be careful to not change the order 
    gene_names = []
    for index, row in matrix_genes_sac_transcript.iterrows():
        if type(row[8]) == str:
            lst = re.findall(r' gene_name "(.+)";', row[8])
        gene_names += lst
    
    #matrix_genes_sac_final contains 5 columns: index, chr, start, end of transcript and gene_name 
    matrix_genes_sac_final = pd.concat([matrix_genes_sac_transcript[[0,3,4]].reset_index(), pd.DataFrame(gene_names, columns=['gene_name'])], axis=1).rename(columns={0: 'chr', 3: "start", 4:'end'})
    matrix_genes_sac_final = matrix_genes_sac_final[matrix_genes_sac_final.chr != 'chrM'] # discard the mitochondrial genes

    # finding the trna genes and keeping the protein coding genes 
    trnas = matrix_genes_sac_final[matrix_genes_sac_final['gene_name'].str.contains(r'\(.+\)') == True]['gene_name'].to_list() #discard the tRNA genes
    matrix_genes_sac_final = matrix_genes_sac_final.loc[~matrix_genes_sac_final['gene_name'].isin(trnas)]

    
    matrix_genes_sac_final=matrix_genes_sac_final.reset_index()
    # drop the column index     
    matrix_genes_sac_final = matrix_genes_sac_final.drop(columns=['index'] ,axis=0)
    matrix_genes_sac_final = matrix_genes_sac_final.drop(columns=['level_0'])

    #export a bed file 
    matrix_genes_sac_final.to_csv("./Scerevisiae_All_genes.bed", sep = '\t', index = False)
    
    return matrix_genes_sac_final

gtf_all_genes=gtf_to_bed(path_to_gtf)

#----------------------------------------------------------------------------------------------------------------------------------------------------


# Computing and storing all the intergenic distances // 
# Overllaping genes have zero distance 
# The algorithm does not consider strands but the coordinates should all be based on a single strand 

# coordinate_df is the dataframe that contains the genes and their coordinates produced at the previous step or provided by the user// 
# BE CAREFUL the columns should be named chr start end gene_name //
# Better be sorted based on chr and start columns

def intergenic_distances(coordinate_df):  
    
    # sort the dataframes based on chr name and start coordinate  // 
    coordinate_df = coordinate_df.sort_values(by=['chr','start'])
    
    #1) for each of the chromosomes and all the combinations of genes compute the distance between them // a pair of genes is not taken twice // no directionality
     #here the distances will be stored for each pair of genes and then for each chromosome
    intergenic_dist_all_genes = {}

    lst_chr=coordinate_df['chr'].unique()
    
    for chrom in lst_chr:
        print('Computing distances on chr: ',chrom)

        df = coordinate_df[coordinate_df['chr'] == chrom]
        df=df.reset_index() #index has to be reset because it will be used in order to access the genes
    
        intergenic_dist_all_genes_chrom = {}

        
        assert df['start'].is_monotonic # assert that the coordinates are sorted 


        for i in range(len(df)): # for all the combinations of genes but only once 
            for j in range (i+1,len(df)):
            
               gene_downstream_start = df.at[j,'start']
               gene_downstream_end = df.at[j, 'end']
    
               gene_upstream_start = df.at[i,'start']
               gene_upstream_end = df.at[i,'end']
    
               #Different distance choices
               
               #initial
               #distance = gene_downstream_start - gene_upstream_end
               
               #take the median base as coordinate
               #distance= (gene_downstream_end-(gene_downstream_end-gene_downstream_start)/2)-(gene_upstream_end-(gene_upstream_end-gene_upstream_start)/2)

               #normalize the initial distance with the genes length
               distance = (gene_downstream_start - gene_upstream_end)/((gene_downstream_end-gene_downstream_start)+(gene_upstream_end-gene_upstream_start))
           
               if distance < 0: #if genes are overlapping then distance = 0 // depending on the distance metric some overlapping genes may not have zero distance 
                   distance = 0
           
               #storing the results in a dictionary --> the order of genes in the name are the same as the order on the chromosomes
               intergenic_dist_all_genes_chrom[df.at[i,'gene_name']+'_'+df.at[j,'gene_name']] = {'chr':chrom, 'gene1':df.at[i,'gene_name'],'gene2':df.at[j,'gene_name'],'distance': float(distance)}
               intergenic_dist_all_genes_chrom[df.at[j,'gene_name']+'_'+df.at[i,'gene_name']] = {'chr':chrom, 'gene1':df.at[j,'gene_name'],'gene2':df.at[i,'gene_name'],'distance': float(distance)}

        #store the results for each chromosome also in a dictionary 
        intergenic_dist_all_genes[chrom] = intergenic_dist_all_genes_chrom

    return intergenic_dist_all_genes

import time
start = time.time()
intergenic_all_distances=intergenic_distances(gtf_all_genes)
end = time.time()
print(end - start)

# -----------------------------------------------------------------------------------------------------------------------------
    # Compute z-score–normalized inter-gene distances for each chromosome.
    # 
    # For each chromosome present in df_significant_zscores, this function:
    # 1. Uses the distribution of consecutive inter-gene distances.
    # 2. Computes the mean and standard deviation of all distances on that chromosome.
    # 3. Converts each raw distance into a z-score relative to the chromosome-specific background.
    # 
    # This allows comparison of observed clustering patterns to the expected 
    # chromosome-specific distance distribution.

def compute_z_distances(df_significant_zscores, consecutive_distances):
      #find the intra cluster z-score for positionally clustered cases
      z_scores_in_cluster_comparison_d = {}
      for index, row  in df_significant_zscores.iterrows():
          chrom = row['chr']
  
          mean_all_genes = np.mean(consecutive_distances[chrom])
          std_all_genes = np.std(consecutive_distances[chrom])
         
          z_scores = []
          for dist in consecutive_distances[chrom]:
              z_scores += [(dist-mean_all_genes)/std_all_genes]
                 
          z_scores_in_cluster_comparison_d[chrom]=z_scores
      return z_scores_in_cluster_comparison_d

# -----------------------------------------------------------------------------------------------------------------------------

# Identify sub-clusters of positionally clustered genes based on intra-cluster
#     z-score distances and summarize their genomic and statistical properties.
# 
#     For each chromosome:
#     1. Break large linear clusters into sub-clusters using a z-score threshold
#        on inter-gene distances.
#     2. Assign cluster labels to consecutive genes.
#     3. Remove singleton clusters (clusters with only one gene).
#     4. Compute genomic coordinates, density, and preference/avoidance statistics
#        for each resulting sub-cluster.

def subclustering(z_scores_in_cluster_comparison_d,neighbouring_pairs_genes,significant_z_scores,coordinates_df,dict_preference_chromosome,threshold):
    d_clusters ={}
    d_clusters_final={}
    for chrom in z_scores_in_cluster_comparison_d.keys(): 
        
        d_clusters[chrom] = pd.DataFrame(columns=['chr','start', 'end','gene_name','z_score_intraCluster','cluster'])

        array_zscores = np.array(z_scores_in_cluster_comparison_d[chrom]) 
        indices_breaks = np.where(array_zscores>=threshold)[0] #setting the threshold of significance to break into two clusters// 
        lst_of_genes = np.array(neighbouring_pairs_genes[chrom]) # genes 

        #set the number corresponding to each cluster 
        clusters = []
        count=0
        previous_elem=-1
        for elem in indices_breaks:
            clusters+=[int(count)]*(elem-previous_elem) # make the cluster labels 
            count+=1
            previous_elem = elem 
            
        clusters += [int(count)]*(len(lst_of_genes) - previous_elem-1)
        z_scores = z_scores_in_cluster_comparison_d[chrom] + ['-']# the last element is '-' beacause there is no other genes below
        
        # make the final dataframe 
        df_temp = coordinates_df.loc[coordinates_df['gene_name'].isin(lst_of_genes)].reset_index()
        df_temp = df_temp.drop('index', axis=1)
        d_temp = {'z_score_intraCluster':z_scores, 'cluster':clusters}
        df_temp = pd.concat([df_temp, pd.DataFrame(d_temp)], axis=1)
        d_clusters[chrom] = pd.concat([d_clusters[chrom],df_temp], axis=0) # almost final dataframes in a dictionary 
        
        # discard the genes that are alone in a cluster 
        d_clusters_grouped_all =  d_clusters[chrom].groupby(['cluster']).size().reset_index().rename(columns={0:'count'})
        clusters_exclude= d_clusters_grouped_all[d_clusters_grouped_all['count'] == 1]['cluster'].to_list()
        d_clusters_final[chrom] = d_clusters[chrom].loc[~d_clusters[chrom]['cluster'].isin(clusters_exclude)]  
        
    # find the density of sub-clusters and their coordinates from first to last gene // 
    clusters_final_df = pd.DataFrame(columns=['chr','start','end','cluster_label','density','perc_of_clust','genes','preference_avoidance', 'preference_pvalue', 'linear_clustering_z_score'])
    for chrom in d_clusters_final.keys():
        df_temp = d_clusters_final[chrom]
        z_linear = significant_z_scores[significant_z_scores['chr'] == chrom]['z_score'].item() 
        count=0
        for elem in df_temp['cluster'].unique(): # find the start gene and end gene of the cluster in the whole genome dataset 
            genes_coordinates_start = coordinates_df[coordinates_df['gene_name'] == df_temp[df_temp['cluster'] == elem]['gene_name'].iloc[0]].index[0]
            genes_coordinates_end = coordinates_df[coordinates_df['gene_name'] == df_temp[df_temp['cluster'] == elem]['gene_name'].iloc[-1]].index[0] 
            density = len(df_temp[df_temp['cluster'] == elem]['gene_name'].to_list())/len(coordinates_df.loc[genes_coordinates_start:genes_coordinates_end])
    
            if dict_preference_chromosome[dict_preference_chromosome['chr'] == chrom]['observed_number'].item() > dict_preference_chromosome[dict_preference_chromosome['chr'] == chrom]['expected_number'].item():
                preference = 'preferred'
            else:
                preference='avoided'
                
            d_temp = {'chr':chrom, 'start':df_temp[df_temp['cluster'] == elem]['start'].iloc[0],'end':df_temp[df_temp['cluster'] == elem]['end'].iloc[-1],'genes':'kati','cluster_label':count,'density':density,'preference_avoidance':preference,'preference_pvalue':dict_preference_chromosome[dict_preference_chromosome['chr'] == chrom]['pvalue'].item(),'linear_clustering_z_score':z_linear, 'perc_of_clust':significant_z_scores['genes_per_chr'].sum()/significant_z_scores['total_genes'].to_list()[0]}
            df_d_temp = pd.DataFrame(d_temp, index=[0])
            df_d_temp.at[0, 'genes'] = df_temp[df_temp['cluster'] == elem]['gene_name'].to_list()

            clusters_final_df = pd.concat([clusters_final_df, df_d_temp], axis=0)
            count+=1
      
    clusters_final_df.to_csv('subclusters_positionallyClust_Results.bed',sep='\t',index=False)
    return clusters_final_df
     
# -----------------------------------------------------------------------------------------------------------------------------

# instructions
# linearly_clust is function that takes as input a dataframe (chr, start, end, gene_name) with the query genes and their coordinates, 
# The dataframe storing the coordinates of all genes of the query organism 
# The output of the intergenic_distances function 

# Cases of only one target on chromosomes are discarded from the linear - positional clustering test 

# This function evaluates if the input genes are clustered positionally (per chromosome) comparing the real intergenic distances
# to random ones. The output --> dataframe with all z-scores and pvalues per chromosome, dataframe with the preference or avoidance per chromosome
# dataframe with the significant z-scores per chromosome indicating positional clustering and finally a dataframe with the subclusters of positionally clustered genes per chromosome 
 
#be careful. the input dataframe and coordinates_df should have the columns chr, start, end, gene_name
coordinates_df = gtf_all_genes
dataframe = coordinates_df.sample(frac=1)
dataframe = coordinates_df.loc[coordinates_df.index[:1000]] #Take the first *len("total_genes") rows

def linearly_clust(dataframe, coordinates_df, intergenic_all_distances):

    #sort the coordinates
    coordinates_df = coordinates_df.sort_values(by=['chr','start']) # a dataframe with all gene coordinates of the query organism
    
    # store gene number per chromosome
    genes_per_chromosome_counts = coordinates_df.groupby(['chr']).size().reset_index().rename(columns={0:'count'})

    # Firstly compute the intergenic distances between consecutive genes of the input dataframe, for each chromosome --> resulting dictionary with mean intergenic distances per chromosome 
    
    # Dictionaries to store the results 
    neighbouring_distances_all_genes = {} # All interegenic distances per chromosome
    neighbouring_distances_genes_wo_outliers={} # Intergenic distances but without the 5% highest percentile of distances
    neighbouring_distances_genes_mean = {} # Mean intergenic distances per chromosomes, based on neighbouring_distances_genes_wo_outliers
    neighbouring_pairs_genes={} # all pairs of genes per chromosome 

    # Counting the genes per chromosome, of the input
    genes_per_chromosome_df = dataframe.groupby(['chr']).size().reset_index().rename(columns={0:'count'})  
    chrom_with_targets = genes_per_chromosome_df[genes_per_chromosome_df['count'] != 1]['chr'].unique()  # keeping chromosomes with more than one gene

    for chrom in chrom_with_targets:
        
        df = dataframe[dataframe['chr'] ==  chrom] # keep the genes for each chromosome separately
    
        #extracting the coordinates of the genes from coordinates_df, sorted as positined on the linear chromosome 
        df_genes = coordinates_df.loc[coordinates_df['gene_name'].isin(df['gene_name'].to_list())]
        assert df_genes['start'].is_monotonic #checking that they are sorted
        
        
        dic_chr = intergenic_all_distances[chrom] # loading the dictionary of all gene-pair distances for the specific chromosome 

        genes_lst = df_genes['gene_name'].to_list() #gene names based on their order on the chromosome 
        neighbouring_pairs_genes[chrom] = genes_lst #store them in a dataframe

        # store the intergenic distances of the input genes 
        average_distance = [dic_chr[genes_lst[i]+'_'+genes_lst[i+1]]['distance'] for i in range (len(genes_lst)-1)] #finding all consecutive pairs of genes

        neighbouring_distances_all_genes[chrom] = average_distance #store all intergenic distances in a dataframe
        
        average_distance = np.array(average_distance) # make the list into an numpy array to handle better 
        if len(average_distance) <= 4: #when <= 5 genes on the chromosome don't discard the outliers
            neighbouring_distances_genes_wo_outliers[chrom] = average_distance # store again the distances
            average_distance = average_distance.mean() # copmute the mean intergenic distance for this chromosome
        else: # if there are more than 5 genes then discard the 5% percentile 
            average_distance=average_distance[np.where(average_distance < np.percentile(average_distance,95))] #keeping only distances less than the 5% percentile of distances 
            neighbouring_distances_genes_wo_outliers[chrom] = average_distance # store again the distances
            average_distance = average_distance.mean() # computing the mean intergenic distance
    
       
        # store the results in a final dictionary "neighbouring_distances_genes_mean" which contains the mean intergenic distances per chromosome
        neighbouring_distances_genes_mean[chrom] = average_distance
        
    #---------------------------------------------------------------------------------------------------------
    # Here starts the random permutations (1000). 1000 random intergenic mean distances will be computed for each chromosome.
    # The same number of genes as in the dataframe will be randomly chosen and random integenic distances of the consecutive random genes will be used to create a random distribution of mean distances
    # Finally z-scores and p values indicating the difference between the observed mean intergenic distances and the expected, will be computed
    
    shuffled_mean_distances_of_genes = {} # here will be stored the average distances
    dictionary_of_genes_counts_random_chr = pd.DataFrame(columns=['chr','count']) # storing the number of genes per chromosomes for each loop    

    total_genes = genes_per_chromosome_df['count'].sum() # number of genes in the input dataframe 
    df_shuffled = coordinates_df.copy() # the dataframe that will be used to shuffle in each permutation
    
    count_loop = 0 # to keep the count of the permutations 
    for j in range(1000):
        count_loop +=1
        print(count_loop)
        # shuffling for the procedure of chromosome preference
        df_shuffled_temp = df_shuffled.sample(frac=1) # shuffling the rows of the dataset 
        
        shuffled_genes = df_shuffled_temp.loc[df_shuffled_temp.index[:total_genes]] #Take the first *len("total_genes") rows
        dictionary_of_genes_counts_random_chr = pd.concat([dictionary_of_genes_counts_random_chr, shuffled_genes.groupby(['chr']).size().reset_index().rename(columns={0:'count'})], axis=0) # count genes per chromosome
        
        for chrom in chrom_with_targets:
            
            chr_genes = genes_per_chromosome_df[genes_per_chromosome_df['chr']==chrom]['count'].item() #keep the number of real genes per chromosome
            shuffled_genes_part = df_shuffled_temp[df_shuffled_temp['chr'] == chrom] # keep genes per chromosome separately
            shuffled_genes_part = shuffled_genes_part.loc[shuffled_genes_part.index[:chr_genes]]
            shuffled_genes_part =  shuffled_genes_part.sort_values(by=['start']) #sort the genes based on their order on the chromosome
            #assert shuffled_genes_part['start'].is_monotonic # assert they are sorted 
            genes_lst = coordinates_df.loc[coordinates_df['gene_name'].isin(shuffled_genes_part['gene_name'].to_list())]['gene_name'].to_list()

            dic_chr = intergenic_all_distances[chrom] # loading the dictionary of distances for the specific chromosome 

            # store the intergenic distances of the input genes 
            average_distance = [dic_chr[genes_lst[i]+'_'+genes_lst[i+1]]['distance'] for i in range (len(genes_lst)-1)] #finding all consecutive pairs of genes
        
            average_distance = np.array(average_distance) # make the list into an numpy array to handle better 
            if len(average_distance) <= 4: #when <= 5 genes on the chromosome don't discard the outliers
                average_distance = average_distance.mean() # copmute the mean intergenic distance for this chromosome
            else: # if there are more than 5 genes then discard the 5% percentile 
                average_distance=average_distance[np.where(average_distance < np.percentile(average_distance,95))] #keeping only distances less than the 5% percentile of distances 
                average_distance = average_distance.mean() # computing the mean intergenic distance
        
            if chrom in shuffled_mean_distances_of_genes.keys():
                shuffled_mean_distances_of_genes[chrom]=np.append(shuffled_mean_distances_of_genes[chrom], average_distance)
            else:
                shuffled_mean_distances_of_genes[chrom] = average_distance

        #print(count_loop)
    #save the shuffled resutls in a pickle file 
    #a_file = open("shuffled_mean_distances_1000perm_{name}.pkl".format(name=name_1), "wb")
    #pickle.dump(shuffled_mean_distances_of_genes, a_file)
    #a_file.close()
    
    # Here starts the evaluation 
    #------------------------------------------------------------------------------------------------------------
    dict_preference_chromosome = pd.DataFrame(columns=['chr','observed_number','expected_number','pvalue'])
    
    #chromosome preference evaluation 
    chr_shuffled = dictionary_of_genes_counts_random_chr['chr'].unique()
    for chrom in chr_shuffled:
        
        # find the preference or avoidance per chromosome // pvalue computation 
        df_shuffled = dictionary_of_genes_counts_random_chr[dictionary_of_genes_counts_random_chr['chr'] == chrom] # random genes per chromosome (1000)

        if chrom in genes_per_chromosome_df['chr'].unique():
            df_real = genes_per_chromosome_df[genes_per_chromosome_df['chr'] == chrom] # real genes per chromosome
        
            if df_real['count'].item() > df_shuffled['count'].mean(): #computing pvalues
                pvalue = [len(df_shuffled[df_shuffled['count'] >= df_real['count'].item()])/1000]
            else:
                pvalue = [len(df_shuffled[df_shuffled['count'] <= df_real['count'].item()])/1000]
            preference_chrom = {'chr':chrom, 'observed_number':df_real['count'].item(),'expected_number':df_shuffled['count'].mean(),'pvalue':pvalue}

        else:#if there are no genes on that chromosome 
            pvalue = [len(df_shuffled[df_shuffled['count'] <= 0])/1000]
            preference_chrom = {'chr':chrom, 'observed_number':0,'expected_number':df_shuffled['count'].mean(),'pvalue':pvalue}
            
        df_0 = pd.DataFrame(preference_chrom, index=[0])
        dict_preference_chromosome = pd.concat([dict_preference_chromosome, df_0], axis=0)
       
    dict_preference_chromosome.to_csv('chromosome_preference_avoidance_test.bed', index=False, sep='\t')
    
    ###############################################################################################
    # computing the z-scores and the p-values between the real and random mean intergenic distances 
    
    total_z_scores_df = pd.DataFrame(columns=['chr','z_score','pvalue','average_observed_distance','observed_percentage','average_random_distance','std','genes_per_chr','total_genes'])
    for key in neighbouring_distances_genes_wo_outliers.keys():

        # for computing the z-scores
        x = neighbouring_distances_genes_mean[key]
        mean_sample = shuffled_mean_distances_of_genes[key].mean()
        std_sample =  shuffled_mean_distances_of_genes[key].std()
    
        # for computing the pvalues 
        df_shuffled = shuffled_mean_distances_of_genes[key]

        if x > mean_sample:
            pvalue = [((df_shuffled >= x).sum())/1000]
        else:
            pvalue = [((df_shuffled <= x).sum())/1000]
        
        z_score = {'chr':key, 'z_score':((x-mean_sample)/std_sample),'pvalue':pvalue,'average_observed_distance': x,'observed_percentage':(len(neighbouring_distances_all_genes[key])+1)/genes_per_chromosome_counts[genes_per_chromosome_counts['chr'] == key]['count'].item(),'average_random_distance': mean_sample,'std': std_sample, 'genes_per_chr': (len(neighbouring_distances_all_genes[key])+1),'total_genes':len(dataframe)}
        df_0 = pd.DataFrame(z_score, index=[0])
        total_z_scores_df = pd.concat([total_z_scores_df, df_0], axis=0)
    
    #total_z_scores_df.to_csv("all_z_scores_df_1000perm_{name}.bed".format(name=name), sep = '\t', index = False) #store all the results 

    total_z_scores_df_more5 = total_z_scores_df[(total_z_scores_df['genes_per_chr'] > 5 )& (abs(total_z_scores_df['z_score']) >= 1.96)] # keep only cases with more than 5 genes per chromosome
    
    #total_z_scores_df_more5.to_csv("significant_z_scores_morethan5genes_df_{name}.bed".format(name=name), sep = '\t', index = False)
    
    ##############################################################################################
    ##############################################################################################
    # Identify positionally clustered chromosomes and derive sub-clusters
    
    # Keep only chromosomes with sufficient genes and significant clustering
    total_z_scores_df_more5 = total_z_scores_df[
        (total_z_scores_df['genes_per_chr'] > 5) &
        (abs(total_z_scores_df['z_score']) >= 1.96)
    ]
    
    # Select only significantly positionally clustered cases
    # Negative z-score = smaller-than-expected distances (clustering)
    significant_z_scores = total_z_scores_df_more5[
        total_z_scores_df_more5['z_score'] < 0
    ]
    
    # Initialize outputs to avoid UnboundLocalError
    positionally_clustered_subclusters = pd.DataFrame()
    positionally_distanced_subclusters = pd.DataFrame()
    
    # ------------------------------------------
    # Sub-clustering for positionally clustered chromosomes
    # ------------------------------------------
    if len(significant_z_scores) >= 1:
    
        # Compute intra-cluster distance z-scores
        z_scores_in_cluster_comparison_d = compute_z_distances(
            significant_z_scores,
            neighbouring_distances_all_genes
        )
    
        # Divide linear clusters into sub-clusters
        positionally_clustered_subclusters = subclustering(
            z_scores_in_cluster_comparison_d,
            neighbouring_pairs_genes,
            significant_z_scores,
            coordinates_df,
            dict_preference_chromosome,
            threshold=2
        )
    
    else:
        print("No significantly positionally clustered chromosomes detected.")
    
    # ------------------------------------------
    # OPTIONAL: Sub-clustering for positionally dispersed cases
    # (currently disabled but kept for future use)
    # ------------------------------------------
    """
    significant_z_scores_pos = total_z_scores_df_more5[
        total_z_scores_df_more5['z_score'] > 0
    ]
    
    if len(significant_z_scores_pos) >= 1:
    
        z_scores_in_cluster_comparison_d = compute_z_distances(
            significant_z_scores_pos,
            neighbouring_distances_all_genes
        )
    
        positionally_distanced_subclusters = subclustering(
            z_scores_in_cluster_comparison_d,
            neighbouring_pairs_genes,
            significant_z_scores_pos,
            coordinates_df,
            dict_preference_chromosome,
            threshold=2
        )
    """
    
    # ------------------------------------------
    # Return all main results
    # ------------------------------------------
    return  total_z_scores_df,total_z_scores_df_more5,dict_preference_chromosome,positionally_clustered_subclusters, positionally_distanced_subclusters


def sub_clusters_elaborate(df, df_coordinates, category, name_output):
    """
    Expand subcluster-level results into gene-level annotations.

    For each subcluster entry (containing a list of genes), this function:
    - Extracts gene coordinates
    - Assigns the cluster label
    - Assigns the provided gene category
    - Returns a gene-level DataFrame suitable for genome browser tracks or R plotting

    Parameters
    ----------
    df : pd.DataFrame
        Subcluster-level results containing a 'genes' column (list of gene names)
        and 'cluster_label'.

    df_coordinates : pd.DataFrame
        Gene coordinate table with columns ['chr', 'start', 'end', 'gene_name'].

    category : str
        Column name in df indicating the gene category (e.g. TF class, feature type).

    name_output : str
        String used for output filename.

    Returns
    -------
    pd.DataFrame
        Gene-level dataframe with coordinates, category, and cluster labels.
    """

    df_final = pd.DataFrame(
        columns=['chr', 'start', 'end', 'gene_name', category, 'cluster_label']
    )

    for _, row in df.iterrows():

        # Extract coordinates for genes in this subcluster
        genes = df_coordinates.loc[
            df_coordinates['gene_name'].isin(row['genes'])
        ].copy()

        # Assign category and cluster label to all genes in this subcluster
        genes[category] = row[category]
        genes['cluster_label'] = row['cluster_label']

        df_final = pd.concat([df_final, genes], axis=0)

    # Save gene-level subcluster annotations
    output_file = f'gene_categories_subclusters_{name_output}.bed'
    df_final.to_csv(output_file, sep='\t', index=False)

    return df_final

def compute_perc_of_clustering(dataframe, z_score_df, column_name, output_name):
    """
    Compute percentage of clustered genes per category.

    For each category, computes:
        sum(genes_per_chr) / total_genes

    Parameters
    ----------
    dataframe : pd.DataFrame
        Subcluster dataframe containing category column.

    z_score_df : pd.DataFrame
        Z-score results dataframe containing genes_per_chr and total_genes.

    column_name : str
        Column defining gene category or feature.

    output_name : str
        Output filename.

    Returns
    -------
    pd.DataFrame
        Updated dataframe with 'perc_of_clust' column.
    """

    perc_map = {}

    for cat in dataframe[column_name].unique():
        df_spec = z_score_df[z_score_df[column_name] == cat]

        if len(df_spec) == 0:
            perc_map[cat] = 0
        else:
            perc_map[cat] = (
                df_spec['genes_per_chr'].sum() /
                df_spec['total_genes'].iloc[0]
            )

    dataframe['perc_of_clust'] = dataframe[column_name].map(perc_map)

    dataframe.to_csv(output_name, sep='\t', index=False)

    return dataframe


def subclusters_plot_across_categories(
    significant_z_scores,
    clusters_final_df,
    coordinate_df,
    name,
    assembly,
    category
):
    """
    Plot subclusters across binned chromosomes as heatmaps.

    For each chromosome:
    - Create genomic bins (1 kb)
    - Map subcluster genes to bins
    - Build a matrix encoding subcluster membership
    - Plot heatmap of subclusters across chromosome

    This visualization highlights spatial organization of subclusters.

    Parameters
    ----------
    significant_z_scores : pd.DataFrame
        Chromosome-level significant clustering results.

    clusters_final_df : pd.DataFrame
        Subcluster-level results including 'genes' and 'cluster_label'.

    coordinate_df : pd.DataFrame
        Gene coordinate table.

    name : str
        Name used for output files.

    assembly : str
        Genome assembly ('saccer2' or 'saccer3').

    category : str
        Feature/category column name.
    """

    # Load chromosome sizes
    if assembly == 'saccer2':
        chromSizes_df = pd.read_csv(
            "/home/labuser/Athanasia/ChromHMM/CHROMSIZES/saccer2.txt",
            sep='\t',
            names=["Chr_name", "size"]
        )
    else:
        chromSizes_df = pd.read_csv(
            "/home/labuser/Athanasia/TFs_sacc/saccer3_chromSizes.bed",
            sep='\t',
            names=["Chr_name", "size"]
        )

    # Keep only nuclear chromosomes
    chromSizes_df = chromSizes_df.iloc[:16]

    for chrom in significant_z_scores['chr'].unique():

        # Define 1 kb bins along chromosome
        chr_size = chromSizes_df.loc[
            chromSizes_df['Chr_name'] == chrom, 'size'
        ].item()

        chr_coordinates = np.arange(0, chr_size + 1, 1000)
        chr_coordinates = np.append(chr_coordinates, chr_size)

        df_temp = significant_z_scores[
            significant_z_scores['chr'] == chrom
        ]

        cluster_ids = np.arange(1, len(df_temp) + 1)

        # Initialize matrix: rows = clusters, cols = bins
        N = cluster_ids.shape[0]
        D = chr_coordinates.shape[0]

        matrix = np.zeros((N, D), dtype=int)
        lst_features = []

        max_cluster_label = 0

        for i, (_, row) in enumerate(df_temp.iterrows()):

            df_cluster_temp = clusters_final_df[
                (clusters_final_df[category] == row[category]) &
                (clusters_final_df['chr'] == chrom)
            ]

            lst_features.append(row[category])

            for cluster in df_cluster_temp['cluster_label'].unique():
                max_cluster_label = max(max_cluster_label, cluster)

                genes = df_cluster_temp[
                    df_cluster_temp['cluster_label'] == cluster
                ]['genes'].item()

                for gene in genes:
                    start = coordinate_df.loc[
                        coordinate_df['gene_name'] == gene, 'start'
                    ].item()

                    end = coordinate_df.loc[
                        coordinate_df['gene_name'] == gene, 'end'
                    ].item()

                    # Find bin indices
                    idx_start = np.searchsorted(chr_coordinates, start) - 1
                    idx_end = np.searchsorted(chr_coordinates, end)

                    matrix[i, idx_start:idx_end] = cluster + 1

        # ---------------- Heatmap ----------------
        sns.set_theme(style="darkgrid", font_scale=1)
        plt.figure(figsize=(20, 10))

        ax = sns.heatmap(
            matrix,
            vmin=1,
            vmax=max_cluster_label + 1,
            cmap=sns.color_palette("Paired", max_cluster_label + 1),
            mask=matrix == 0,
            yticklabels=lst_features
        )

        ax.set(
            title=f'{name} genes — subclusters on chromosome {chrom}',
            xlabel="Genomic bins (1 kb)"
        )

        plt.tight_layout()
        plt.savefig(f'{name}_subclusters_on_{chrom}.png', dpi=300)
        plt.close()
