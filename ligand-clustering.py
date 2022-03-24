#!/usr/bin/python3

import numpy as np
import json
import requests
import csv
import pickle
import os
import sys
import re
import statistics
import scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy import cluster
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_samples, silhouette_score
from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import rdFMCS, AllChem, Draw, Lipinski
from itertools import chain, combinations


# import ligand smiles strings
ligand_file = 'lig_1118.json'

with open(ligand_file) as ff:
    smiles_dict = json.load(ff)
    
    
# check for ligands with bad or missing smiles strings
bad_smiles = []

for key,value in smiles_dict.items():
    if 'smiles_cactus' not in value or Chem.MolFromSmiles(smiles_dict[key]['smiles_cactus']) is None:
        bad_smiles.append(key)
    

# for each viral protein, make list of residues in consensus pockets
def pocket_residues(consensus_pockets,all_consensus_residues,directory,protnow):
    file = open(directory+'clusters_'+protnow+'.txt','r')
    line_list = file.readlines()
    consensus_pockets[protnow] = {}
    all_consensus_residues[protnow] = []
    
    for line in line_list:
        pocket = line.split()[0].split(':')[0]
        residues = line.split()[1].split(',')[0:-1]
        consensus_pockets[protnow][pocket] = residues
        for res in residues:
            if res not in all_consensus_residues[protnow]:
                all_consensus_residues[protnow].append(res)

    file.close()
    return consensus_pockets, all_consensus_residues


# for each viral protein pocket, make list of pocket ligands each residue contacts
def pocket_residue_ligand_pairs(directory,filenames,consensus_pockets,consensus_pocket_ligands,protnow):
    consensus_pocket_reslig_pairs[protnow] = {}
    
    for pocket,residues in consensus_pockets[protnow].items():
        consensus_pocket_reslig_pairs[protnow][pocket] = {}
        for res in residues:
            consensus_pocket_reslig_pairs[protnow][pocket][res] = []
     
    for fl in filenames:
        file = open(directory+fl,'r')
        line_list = file.readlines()
    
        for line in line_list:
            interaction = line.split()[0].split(':')[0]
            binding_residues = line.split()[-1].split(',')[0:-1]
            ligand = line.split()[0].split('.')[6]

            # viral protein
            if line.split()[0].split('.')[0].split('_')[0]=='nCoV':
                protein = line.split()[0].split('.')[0].split('_')[1]
                if protein=='Spike':
                    protein = 'S'
                           
                for pocket,residues in consensus_pockets[protnow].items():
                    for res in residues:
                        if protein==protnow and res in binding_residues and ligand in consensus_pocket_ligands[protnow][pocket] and ligand not in consensus_pocket_reslig_pairs[protnow][pocket][res]:
                            consensus_pocket_reslig_pairs[protnow][pocket][res].append(ligand)
                            consensus_pocket_reslig_pairs[protnow][pocket][res].sort()

        file.close()
    return consensus_pocket_reslig_pairs


# for each viral protein pocket, make list of filtered ligands that bind (require ligand size>=8)
def filtered_pocket_ligands(directory,filenames,consensus_pockets,protnow,ligs_leaveout):
    consensus_pocket_ligands[protnow] = {}
    
    for pocket,residues in consensus_pockets[protnow].items():
        consensus_pocket_ligands[protnow][pocket] = []
     
    for fl in filenames:
        file = open(directory+fl,'r')
        line_list = file.readlines()
    
        for line in line_list:
            interaction = line.split()[0].split(':')[0]
            binding_residues = line.split()[-1].split(',')[0:-1]
            ligand = line.split()[0].split('.')[6]
            lig_size = line.split()[0].split('.')[7]

            # viral protein
            if line.split()[0].split('.')[0].split('_')[0]=='nCoV':
                protein = line.split()[0].split('.')[0].split('_')[1]
                if protein=='Spike':
                    protein = 'S'
                           
                for pocket,residues in consensus_pockets[protnow].items():
                    if protein==protnow and (set(binding_residues) & set(residues)) and ligand not in ligs_leaveout[protnow]:
                        if len(ligand)<4:
                            if float(lig_size)>=8:
                                if ligand not in consensus_pocket_ligands[protnow][pocket] and ligand in smiles_dict and ligand not in bad_smiles:
                                    consensus_pocket_ligands[protnow][pocket].append(ligand)
                                    consensus_pocket_ligands[protnow][pocket].sort()
        file.close()
    return consensus_pocket_ligands


# calculate normalized Tanimoto distance between all pairs of ligands
def calc_Tanimoto_dist_norm(fp_radius,nBits,chemtax_dict):
    Tdist_dict = {}
    Tdistnorm_dict = {}
    Tdistlist = []
    for i1 in range(0,len(pocket_ligs)):
        lig1 = pocket_ligs[i1]
        if lig1 not in bad_smiles:
            m1 = Chem.MolFromSmiles(smiles_dict[lig1]['smiles_cactus'])
            fp1 = AllChem.GetMorganFingerprintAsBitVect(m1,fp_radius,nBits)
            for i2 in range(i1+1,len(pocket_ligs)):
                lig2 = pocket_ligs[i2]
                if lig2 not in bad_smiles:
                    m2 = Chem.MolFromSmiles(smiles_dict[lig2]['smiles_cactus'])
                    fp2 = AllChem.GetMorganFingerprintAsBitVect(m2,fp_radius,nBits)
                    Tsim = DataStructs.FingerprintSimilarity(fp1,fp2)
                    Tdist = 1-Tsim
                    Tdistlist.append(Tdist)
                    Tdist_dict[(lig1,lig2)] = Tdist
                    
                    if lig1 in chemtax_dict and lig2 in chemtax_dict:
                        if Tdist==0.0 and chemtax_dict[lig1]['kingdom']!='' and chemtax_dict[lig2]['kingdom']=='':
                            chemtax_dict[lig2] = chemtax_dict[lig1]
                        elif Tdist==0.0 and chemtax_dict[lig1]['kingdom']=='' and chemtax_dict[lig2]['kingdom']!='':
                            chemtax_dict[lig1] = chemtax_dict[lig2]
                    
    Tdistavg = sum(Tdistlist)/float(len(Tdistlist))
    Tdiststd = statistics.stdev(Tdistlist)
    
    print(Tdistavg,Tdiststd,min(Tdistlist))
    
    for ligpair in Tdist_dict.keys():
        Tdistnorm_dict[ligpair] = ((Tdist_dict[ligpair]-Tdistavg)/Tdiststd)+15
     
    return Tdistnorm_dict,chemtax_dict,Tdistavg,Tdiststd
                        
    
# calculate Tanimoto distance for pairs of ligands in each pocket
def get_Tanimoto_dist(Tdistlist_dict,consensus_pocket_ligands,protnow,Tdistnorm_dict):
    Tdistlist_dict[protnow] = {}
    for pocket,ligands in consensus_pocket_ligands[protnow].items():
        Tdistlist = []
        for i1 in range(0,len(ligands)):
            lig1 = ligands[i1]
            for i2 in range(i1+1,len(ligands)):
                lig2 = ligands[i2]
                if (lig1,lig2) in Tdistnorm_dict:
                    Tdist = Tdistnorm_dict[(lig1,lig2)]
                    Tdistlist.append(Tdist)
                elif (lig2,lig1) in Tdistnorm_dict:
                    Tdist = Tdistnorm_dict[(lig2,lig1)]
                    Tdistlist.append(Tdist)
        Tdistlist_dict[protnow][pocket] = Tdistlist 
    return Tdistlist_dict


# calculate normalized chemical taxonomy distance between all pairs of ligands
def calc_chemtax_dist_norm():
    Cdist_dict = {}
    Cdistnorm_dict = {}
    Cdistlist = []
    for i1 in range(0,len(pocket_ligs)):
        lig1 = pocket_ligs[i1]
        if lig1 not in bad_smiles:
            for i2 in range(i1+1,len(pocket_ligs)):
                lig2 = pocket_ligs[i2]
                Csim = 0
                if lig2 not in bad_smiles:
                    if chemtax_dict[lig1]['kingdom']!='' and chemtax_dict[lig2]['kingdom']!='':
                        for level in ['kingdom','superclass','class','subclass']:
                            if chemtax_dict[lig1][level]==chemtax_dict[lig2][level]:
                                Csim = Csim + 1
                        Csim = float(Csim)/float(4)
                        Cdist = 1-Csim
                        Cdistlist.append(Cdist)
                        Cdist_dict[(lig1,lig2)] = Cdist 
                    
    Cdistavg = sum(Cdistlist)/float(len(Cdistlist))
    Cdiststd = statistics.stdev(Cdistlist)
    
    print(Cdistavg,Cdiststd,min(Cdistlist))
    
    for ligpair in Cdist_dict.keys():
        Cdistnorm_dict[ligpair] = ((Cdist_dict[ligpair]-Cdistavg)/Cdiststd)+15
     
    return Cdistnorm_dict,Cdistavg,Cdiststd


# calculate distance based on chemical taxonomy for pairs of ligands in each pocket
def get_chemtax_dist(Cdistlist_dict,consensus_pocket_ligands,protnow,Cdistnorm_dict):
    Cdistlist_dict[protnow] = {}
    for pocket,ligands in consensus_pocket_ligands[protnow].items():
        Cdistlist = []
        for i1 in range(0,len(ligands)):
            lig1 = ligands[i1]
            for i2 in range(i1+1,len(ligands)):
                lig2 = ligands[i2]
                if (lig1,lig2) in Cdistnorm_dict:
                    Cdist = Cdistnorm_dict[(lig1,lig2)]
                    Cdistlist.append(Cdist)
                elif (lig2,lig1) in Cdistnorm_dict:
                    Cdist = Cdistnorm_dict[(lig2,lig1)]
                    Cdistlist.append(Cdist)
                else:
                    Cdistlist.append('NA')
        Cdistlist_dict[protnow][pocket] = Cdistlist  
    return Cdistlist_dict


# calculate normalized word context distance between all pairs of ligands
def calc_ligname_dist_norm(ligname_dist_dict):
    LNdistnorm_dict = {}
    LNdist_dict = {}
    LNdistlist = []
    for i1 in range(0,len(pocket_ligs)):
        lig1 = pocket_ligs[i1]
        if lig1 not in bad_smiles:
            for i2 in range(i1+1,len(pocket_ligs)):
                lig2 = pocket_ligs[i2]
                if lig2 not in bad_smiles:
                    if (lig1,lig2) in ligname_dist_dict or (lig2,lig1) in ligname_dist_dict:
                        if (lig1,lig2) in ligname_dist_dict:
                            LNdist = ligname_dist_dict[(lig1,lig2)]
                        elif (lig2,lig1) in ligname_dist_dict:
                            LNdist = ligname_dist_dict[(lig2,lig1)]
                        LNdistlist.append(LNdist)
                        LNdist_dict[(lig1,lig2)] = LNdist 
                    
    LNdistavg = sum(LNdistlist)/float(len(LNdistlist))
    LNdiststd = statistics.stdev(LNdistlist)
    
    print(LNdistavg,LNdiststd,min(LNdistlist))
    
    for ligpair in LNdist_dict.keys():
        LNdistnorm_dict[ligpair] = ((LNdist_dict[ligpair]-LNdistavg)/LNdiststd)+15
     
    return LNdistnorm_dict,LNdistavg,LNdiststd


# make distance matrix from ligand name distance dictionary
def get_ligname_dist(LNdistlist_dict,consensus_pocket_ligands,protnow,LNdistnorm_dict):
    LNdistlist_dict[protnow] = {}
    for pocket,ligands in consensus_pocket_ligands[protnow].items():
        LNdistlist = []
        for i1 in range(0,len(ligands)):
            lig1 = ligands[i1]
            for i2 in range(i1+1,len(ligands)):
                lig2 = ligands[i2]
                if (lig1,lig2) in LNdistnorm_dict:
                    LNdist = LNdistnorm_dict[(lig1,lig2)]
                    LNdistlist.append(LNdist)
                elif (lig2,lig1) in LNdistnorm_dict:
                    LNdist = LNdistnorm_dict[(lig2,lig1)]
                    LNdistlist.append(LNdist)
                else:
                    LNdistlist.append('NA')
        LNdistlist_dict[protnow][pocket] = LNdistlist  
    return LNdistlist_dict


# take weighted average of distance matrices 
def weighted_dist(protnow,consensus_pocket_ligands,Tdistlist_dict,Cdistlist_dict,LNdistlist_dict,Wdistlist_dict):
    Wdistlist_dict[protnow] = {}
    for pocket,ligands in consensus_pocket_ligands[protnow].items():
        Wdistlist_dict[protnow][pocket] = []
        if len(Tdistlist_dict[protnow][pocket])!=len(Cdistlist_dict[protnow][pocket]):
            print('distance problem')
        if len(Cdistlist_dict[protnow][pocket])!=len(LNdistlist_dict[protnow][pocket]):
            print('distance problem')
        for ind in range(0,len(Tdistlist_dict[protnow][pocket])):
            Tdist = Tdistlist_dict[protnow][pocket][ind]
            Cdist = Cdistlist_dict[protnow][pocket][ind]
            LNdist = LNdistlist_dict[protnow][pocket][ind]
            if Cdist!='NA' and LNdist!='NA':
                Wdist = (float(1)/float(3))*Tdist + (float(1)/float(3))*Cdist + (float(1)/float(3))*LNdist
            elif Cdist!='NA' and LNdist=='NA':
                Wdist = 0.5*Tdist + 0.5*Cdist
            elif Cdist=='NA' and LNdist!='NA':
                Wdist = 0.5*Tdist + 0.5*LNdist
            elif Cdist=='NA' and LNdist=='NA':
                Wdist = Tdist
            else:
                print(ind,'no Wdist assignment')                
            Wdistlist_dict[protnow][pocket].append(Wdist) 
    return Wdistlist_dict


# get distance matrix from distance list
def get_dist_matrix(Wdistlist_dict,Wdistmat_dict,consensus_pocket_ligands):
    Wdistmat_dict[protnow] = {}
    for pocket,ligands in consensus_pocket_ligands[protnow].items():
        Wdistmat = -1*np.ones((len(ligands),len(ligands)))
        ind = 0
        for i1 in range(0,len(ligands)):
            for i2 in range(i1,len(ligands)):
                if i1==i2:
                    Wdistmat[i1,i2] = 0
                else:
                    Wdistmat[i1,i2] = Wdistlist_dict[protnow][pocket][ind]
                    Wdistmat[i2,i1] = Wdistlist_dict[protnow][pocket][ind]
                    ind = ind + 1 
        Wdistmat_dict[protnow][pocket] = Wdistmat
    return Wdistmat_dict


# cluster ligands in each pocket using DBSCAN
def dbscan_cluster_pocket_ligands(cluster_dict,protnow,consensus_pocket_ligands,Wdistmat_dict):    
    cluster_dict[protnow] = {}
    for pocket,ligands in consensus_pocket_ligands[protnow].items():
        best_params = {}
        if len(ligands)>2:
            best_params['sscore'] = -1
            eps_vec = np.linspace(1,15,14*20)
            minsamp_vec = np.arange(3,11)
            for eps_val in eps_vec:
                for minsamp_val in minsamp_vec:
                    db = DBSCAN(eps=eps_val, min_samples=minsamp_val, metric='precomputed').fit(Wdistmat_dict[protnow][pocket])
                    labels = db.labels_
                    n_clusters = len(set(labels)) - (1 if -1 in labels else 0) # Number of clusters in labels, ignoring noise if present
                    n_noise = list(labels).count(-1)
            
                    if len(ligands)>n_clusters and n_clusters>=2:
                        ss = silhouette_score(Wdistmat_dict[protnow][pocket],labels,metric='precomputed')
                        if ss > best_params['sscore']:
                            best_params = {'eps': eps_val, 'min_samples': minsamp_val, 'sscore': ss, 'n_clusters': n_clusters, 'n_noise': n_noise}               
            
            if 'eps' in best_params.keys():
                db = DBSCAN(eps=best_params['eps'], min_samples=best_params['min_samples'], metric='precomputed').fit(Wdistmat_dict[protnow][pocket])
                labels = db.labels_
                n_clusters = len(set(labels)) - (1 if -1 in labels else 0) # Number of clusters in labels, ignoring noise if present
                n_noise = list(labels).count(-1)
            
                labels_array = np.empty((len(consensus_pocket_ligands[protnow][pocket]),1),dtype=np.int64)
                for k,ligand in enumerate(consensus_pocket_ligands[protnow][pocket]):
                    labels_array[k] = np.empty((1,),dtype=np.int64)
                    labels_array[k][0] = np.int64(labels[k])
                            
                cluster_dict[protnow][pocket] = labels_array
                
            else:
                labels_array = np.empty((len(consensus_pocket_ligands[protnow][pocket]),1),dtype=np.int64)
                for k,ligand in enumerate(consensus_pocket_ligands[protnow][pocket]):
                    labels_array[k] = np.empty((1,),dtype=np.int64)
                    labels_array[k][0] = np.int64(k)
                            
                cluster_dict[protnow][pocket] = labels_array
              
        elif len(ligands)==2:
            cluster_dict[protnow][pocket] = [[0], [0]]
            
    return cluster_dict


# silhouette plot
def silhouette(Wdistmat_dict,cluster_dict,consensus_pocket_ligands):
    for pocket, clusters in cluster_dict[protnow].items():
        clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
        clustall=[]
        for k in range(max([x[1] for x in clusout])+1):
            clustall.append([x[0] for x in clusout if x[1]==k])
        n_clusters=len(clustall)
        dist_matrix = Wdistmat_dict[protnow][pocket]

        try:
            cluster_labels = np.empty((len(consensus_pocket_ligands[protnow][pocket]),))
            max_clust_ind = np.amax(cluster_dict[protnow][pocket])
            for k in range(0,len(cluster_dict[protnow][pocket])):
                if cluster_dict[protnow][pocket][k][0]==-1:
                    cluster_labels[k] = max_clust_ind+1
                    max_clust_ind = max_clust_ind+1
                else:
                    cluster_labels[k] = cluster_dict[protnow][pocket][k][0]

            # Create a subplot with 1 row and 1 column
            fig = plt.figure()
            fig.set_size_inches(9, 6)
            ax=fig.add_subplot(111)
            
            # The (n_clusters+1)*10 is for inserting blank space between silhouette plots of individual clusters, to demarcate them clearly.
            ax.set_ylim([0, len(cluster_labels) + (n_clusters + 1) * 10])

            # The silhouette_score gives the average value for all the samples.
            silhouette_avg = silhouette_score(dist_matrix, cluster_labels, metric="precomputed", sample_size=None)
            print("There are ",n_clusters," clusters and the average silhouette_score is : ",silhouette_avg)

            # Compute the silhouette scores for each sample
            sample_silhouette_values = silhouette_samples(dist_matrix, cluster_labels, metric="precomputed")

            y_lower = 10
            for i in range(n_clusters):
                # Aggregate the silhouette scores for samples belonging to cluster i, and sort them
                ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

                ith_cluster_silhouette_values.sort()

                size_cluster_i = ith_cluster_silhouette_values.shape[0]
                y_upper = y_lower + size_cluster_i

                color = cm.nipy_spectral(float(i)/n_clusters)
                ax.fill_betweenx(np.arange(y_lower, y_upper),0, ith_cluster_silhouette_values,facecolor=color, edgecolor=color, alpha=0.7)

                # Compute the new y_lower for next plot
                y_lower = y_upper + 10  # 10 for the 0 samples

            plt.title(protnow+', Pocket '+pocket,fontsize=16)
            plt.xlabel("Silhouette coefficient",fontsize=16)
            plt.ylabel("Cluster label",fontsize=16)

            # vertical line for average silhouette score of all the values
            ax.axvline(x=silhouette_avg, color="red", linestyle="--")

            ax.set_yticks([])  # Clear the yaxis labels / ticks
            plt.xlim((-0.1,0.6))
            ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
            plt.xticks(fontsize=14)

            plt.show()
            #plt.savefig('figures/silhouette_plot_'+protnow+'_'+pocket+'.png')
        
        except:
            if len(cluster_dict[protnow][pocket])==n_clusters:
                print('All ligands clustered separately')
            elif n_clusters==1:
                print('All ligands clustered together')
    
    return 

    
# find max common substructure for ligand cluster in each pocket
def pocket_mcs(cluster_dict,consensus_pocket_ligands):
    os.system('rm images/CCC-15-10-'+gdccut+'-4-0-ligs-8-current-resall/'+protnow+'_pocket*')
    with open('ligand-cluster-key-CCC-15-10-'+gdccut+'-4-0-ligs-8-current-resall-'+protnow+'.csv','w') as f:
        writeCSV = csv.DictWriter(f,fieldnames=['Pocket','Cluster Index','Ligands in Cluster'])
        writeCSV.writeheader()
        for pocket, clusters in cluster_dict[protnow].items():
            clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
            clustall=[]
            for k in range(max([x[1] for x in clusout])+1):
                clustall.append([x[0] for x in clusout if x[1]==k])
            n_clusters=len(clustall)
            for cind,clust in enumerate(clustall,1):
                newrow = {'Pocket': pocket, 'Cluster Index': cind, 'Ligands in Cluster': clust}
                writeCSV.writerow(newrow)
                
                subclasses = []
                for lig in clust:
                    if chemtax_dict[lig]['subclass'] not in subclasses:
                        subclasses.append(chemtax_dict[lig]['subclass'])
                if len(clust)>1:
                    print('Pocket',pocket,'\t','Cluster',cind,'\t','# Ligands',len(clust),'\t',\
                          'Class',chemtax_dict[clust[0]]['class'],'\t','# Subclasses',len(subclasses))
                
                if len(clust)>=3:
                    molecules = []
                    for lig in clust:
                        if lig not in bad_smiles:
                            molecules.append(Chem.MolFromSmiles(smiles_dict[lig]['smiles_cactus']))
                    if len(molecules)>1:
                        mcs = Chem.rdFMCS.FindMCS(molecules)
                        if mcs.numAtoms>=4:
                            mcs_smiles = Chem.MolToSmiles(Chem.MolFromSmarts(mcs.smartsString))
                            mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                            mcs_coords = AllChem.Compute2DCoords(mcs_mol)
                            image_file = 'images/CCC-15-10-'+gdccut+'-4-0-ligs-8-current-resall/'\
                            +protnow+'_pocket'+pocket+'_cluster'+str(cind)+'_mcs.png'
                            Draw.MolToFile(mcs_mol,image_file)
    return


def fraction_cluster_contacts_heatmap(cluster_dict,consensus_pocket_ligands,consensus_pocket_reslig_pairs,protnow,fraction_ligand_contacts_matrix_dict):
    fraction_ligand_contacts_matrix_dict[protnow] = {}
    for pocket, clusters in cluster_dict[protnow].items():
        sorted_residues = sorted(consensus_pocket_reslig_pairs[protnow][pocket].keys(), key = lambda r: int(r[1:-2]))
        clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
        clustall=[]
        for k in range(max([x[1] for x in clusout])+1):
            clustall.append([x[0] for x in clusout if x[1]==k])
        n_clusters=len(clustall)
        flc_matrix = -1*np.ones((n_clusters,len(consensus_pocket_reslig_pairs[protnow][pocket].keys())))
        for cind,clust in enumerate(clustall,1):
            for res in sorted_residues:
                overlap = set(consensus_pocket_reslig_pairs[protnow][pocket][res]).intersection(set(clust))
                resind = sorted_residues.index(res)
                flc_matrix[cind-1,resind] = len(overlap)/float(len(clust))   
                
        fraction_ligand_contacts_matrix_dict[protnow][pocket] = (flc_matrix, sorted_residues)
        
        ## heatmap using matplotlib (colorbar has same range for all plots)
        cb_viridis = cm.get_cmap('viridis', 100)
        plt.figure()
        plt.pcolor(np.arange(len(sorted_residues)), np.arange(n_clusters), flc_matrix, cmap=cb_viridis, vmin=0, vmax=1, shading='auto')
        plt.title(protnow+', Pocket '+str(pocket))
        plt.xlabel('Residues')
        plt.ylabel('Clusters') 
        plt.xticks(ticks=np.arange(len(sorted_residues)), labels=sorted_residues, rotation=90)
        plt.yticks(ticks=list(np.arange(n_clusters)), labels=list(np.arange(1,n_clusters+1)))
        plt.colorbar(label='Fraction of Cluster Ligands in Contact')
        plt.show()

    return fraction_ligand_contacts_matrix_dict


# save files with PDB IDs for ligands in each cluster
def save_ligand_clusters(cluster_dict,consensus_pocket_ligands):
    with open('ligand-cluster-key-CCC-15-10-'+gdccut+'-4-0-ligs-8-current-resall-'+protnow+'.csv','w') as f:
        writeCSV = csv.DictWriter(f,fieldnames=['Pocket','Cluster Index','Ligands in Cluster'])
        writeCSV.writeheader()
        for pocket, clusters in cluster_dict[protnow].items():
            clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
            clustall=[]
            for k in range(max([x[1] for x in clusout])+1):
                clustall.append([x[0] for x in clusout if x[1]==k])
            n_clusters=len(clustall)
            for cind,clust in enumerate(clustall,1):
                newrow = {'Pocket': pocket, 'Cluster Index': cind, 'Ligands in Cluster': clust}
                writeCSV.writerow(newrow)


# compare experimentally screened positive compounds with ligand clusters
def compare_exp_pos_compounds(protnow,smiles_dict,cluster_dict,consensus_pocket_ligands):
    cluster_distance_pos = {}
    if protnow=='nsp5':
        with open('./fret_crys_test1.csv','r') as pos_smiles_file:
            readCSV = csv.DictReader(pos_smiles_file)
            for row in readCSV:
                exp_smiles = str(row['smiles'])
                exp_id = row['compound_id']
                m1 = Chem.MolFromSmiles(exp_smiles)
                fp1 = AllChem.GetMorganFingerprintAsBitVect(m1,fp_radius,nBits)
                image_file = 'images/experimental_compounds_positive/'+exp_id+'.png'
                Draw.MolToFile(m1,image_file)
                cluster_distance_pos[exp_id] = {}
                for pocket, clusters in cluster_dict[protnow].items():
                    cluster_distance_pos[exp_id][pocket] = {}
                    clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
                    clustall=[]
                    for k in range(max([x[1] for x in clusout])+1):
                        clustall.append([x[0] for x in clusout if x[1]==k])
                    n_clusters=len(clustall)
                    for cind,clust in enumerate(clustall,1):
                        Tdist_avg = 0
                        for lig in clust:
                            # calculate Tanimoto distance
                            m2 = Chem.MolFromSmiles(smiles_dict[lig]['smiles_cactus'])
                            fp2 = AllChem.GetMorganFingerprintAsBitVect(m2,fp_radius,nBits)
                            Tsim = DataStructs.FingerprintSimilarity(fp1,fp2)
                            Tdist = 1-Tsim
                            Tdist_avg = Tdist_avg + Tdist
                        Tdist_avg = Tdist_avg/float(len(clust))
                        cluster_distance_pos[exp_id][pocket][cind] = Tdist_avg
                        
    return cluster_distance_pos


# compare experimentally screened negative compounds with ligand clusters
def compare_exp_neg_compounds(protnow,smiles_dict,cluster_dict,consensus_pocket_ligands):
    cluster_distance_neg = {}
    if protnow=='nsp5':
        with open('./all_postera_match_neg_10132021.csv','r') as neg_smiles_file:
            readCSV = csv.DictReader(neg_smiles_file)
            for row in readCSV:
                exp_smiles = str(row['SMILES'])
                exp_id = row['compound_id']
                m1 = Chem.MolFromSmiles(exp_smiles)
                try: 
                    fp1 = AllChem.GetMorganFingerprintAsBitVect(m1,fp_radius,nBits)
                    image_file = 'images/experimental_compounds_negative/'+exp_id+'.png'
                    Draw.MolToFile(m1,image_file)
                    cluster_distance_neg[exp_id] = {}
                    for pocket, clusters in cluster_dict[protnow].items():
                        cluster_distance_neg[exp_id][pocket] = {}
                        clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
                        clustall=[]
                        for k in range(max([x[1] for x in clusout])+1):
                            clustall.append([x[0] for x in clusout if x[1]==k])
                        n_clusters=len(clustall)
                        for cind,clust in enumerate(clustall,1):
                            Tdist_avg = 0
                            for lig in clust:
                                # calculate Tanimoto distance
                                m2 = Chem.MolFromSmiles(smiles_dict[lig]['smiles_cactus'])
                                fp2 = AllChem.GetMorganFingerprintAsBitVect(m2,fp_radius,nBits)
                                Tsim = DataStructs.FingerprintSimilarity(fp1,fp2)
                                Tdist = 1-Tsim
                                Tdist_avg = Tdist_avg + Tdist
                            Tdist_avg = Tdist_avg/float(len(clust))
                            cluster_distance_neg[exp_id][pocket][cind] = Tdist_avg
                except:
                    print(exp_id)
                    pass
                        
    return cluster_distance_neg


# find most similar ligand cluster for experimentally screened compounds 
def find_closest_cluster(cluster_distance,protnow,cluster_dict,consensus_pocket_ligands):
    closest_cluster = {}
    closest_dist_list = {}
    if protnow=='nsp5':
        for pocket, clusters in cluster_dict[protnow].items():
            closest_dist_list[pocket] = {}
            clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
            clustall=[]
            for k in range(max([x[1] for x in clusout])+1):
                clustall.append([x[0] for x in clusout if x[1]==k])
            n_clusters=len(clustall)
            for cind,clust in enumerate(clustall,1):
                closest_dist_list[pocket][cind] = []
        
        for exp_id,itm in cluster_distance.items():
            closest_cluster[exp_id] = {}
            for pocket,cind in cluster_distance[exp_id].items():
                min_dist = 1
                for cind,Tdist in cluster_distance[exp_id][pocket].items():
                    if Tdist<min_dist:
                        min_dist = Tdist
                        closest_cluster[exp_id][pocket] = (cind,Tdist)
                    
        for pocket, clusters in cluster_dict[protnow].items():
            clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
            clustall=[]
            for k in range(max([x[1] for x in clusout])+1):
                clustall.append([x[0] for x in clusout if x[1]==k])
            n_clusters=len(clustall)
            for cind,clust in enumerate(clustall,1):
                for exp_id,itm in closest_cluster.items():
                    for pckt,tup in closest_cluster[exp_id].items():
                        if str(pckt)==str(pocket) and str(cind)==str(closest_cluster[exp_id][pckt][0]):
                            closest_dist_list[pocket][cind].append(closest_cluster[exp_id][pckt][1])
                    
    return closest_cluster, closest_dist_list


# for each nsp5 ligand cluster, make overlapping histograms of Tanimoto distances for positive and negative compounds closest to that cluster
def closest_cluster_hist(protnow,closest_dist_list_pos,closest_dist_list_neg,cluster_dict,consensus_pocket_ligands):
    if protnow=='nsp5':
        for pocket, clusters in cluster_dict[protnow].items():
            clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
            clustall=[]
            for k in range(max([x[1] for x in clusout])+1):
                clustall.append([x[0] for x in clusout if x[1]==k])
            n_clusters=len(clustall)
            for cind,clust in enumerate(clustall,1):
                print('ligand cluster index :',cind)
                plt.figure()
                
                if cind==1:
                    plt.figtext(0.14,0.8,'p-value = 0.702\nno difference')
                elif cind==2:
                    plt.figtext(0.14,0.8,'p-value = 4.87e-17\ndifferent')
                elif cind==4:
                    plt.figtext(0.14,0.8,'p-value = 1.83e-29\ndifferent')
                elif cind==5:
                    plt.figtext(0.14,0.8,'p-value = 0.0825\nno difference')
                elif cind==7:
                    plt.figtext(0.14,0.8,'p-value = 1.36e-32\ndifferent')
                elif cind==10:
                    plt.figtext(0.14,0.8,'p-value = 5.90e-30\ndifferent')

                if len(closest_dist_list_pos[pocket][cind])>0 and len(closest_dist_list_neg[pocket][cind])>0:
                    print('average distance positive',statistics.mean(closest_dist_list_pos[pocket][cind]))
                    print('average distance negative',statistics.mean(closest_dist_list_neg[pocket][cind]))
                    plt.hist(closest_dist_list_pos[pocket][cind],range=(0.5,1.0),bins=20,density=True,alpha=0.5,label='positive',color='red')
                    plt.hist(closest_dist_list_neg[pocket][cind],range=(0.5,1.0),bins=20,density=True,alpha=0.5,label='negative',color='blue')
                    
                elif len(closest_dist_list_pos[pocket][cind])==0 and len(closest_dist_list_neg[pocket][cind])>0:
                    print('average distance negative',statistics.mean(closest_dist_list_neg[pocket][cind]))
                    plt.hist(closest_dist_list_neg[pocket][cind],range=(0.5,1.0),bins=20,density=True,alpha=0.5,label='negative',color='blue')
                    
                elif len(closest_dist_list_pos[pocket][cind])>0 and len(closest_dist_list_neg[pocket][cind])==0:
                    print('average distance negative',statistics.mean(closest_dist_list_pos[pocket][cind]))
                    plt.hist(closest_dist_list_pos[pocket][cind],range=(0.5,1.0),bins=20,density=True,alpha=0.5,label='positive',color='red')
                
                plt.legend(loc='upper right')
                plt.title(protnow+', Pocket '+pocket+', Cluster '+str(cind),fontsize=16)
                plt.xlabel('Tanimoto distance',fontsize=16)
                plt.xlim(0.5,1.0)
                plt.ylim(0,30)
                plt.ylabel('Normalized count',fontsize=16)
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                plt.show()
                #plt.savefig('figures/exp_compounds_distance_distr_pocket'+pocket+'_cluster'+str(cind)+'.png')  
                
    return


# for each nsp5 ligand cluster, perform t-test for positive and negative compounds closest to that cluster
def closest_cluster_ttest(protnow,closest_dist_list_pos,closest_dist_list_neg,cluster_dict,consensus_pocket_ligands):
    if protnow=='nsp5':
        for pocket, clusters in cluster_dict[protnow].items():
            clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
            clustall=[]
            for k in range(max([x[1] for x in clusout])+1):
                clustall.append([x[0] for x in clusout if x[1]==k])
            n_clusters=len(clustall)
            for cind,clust in enumerate(clustall,1):
                ttest_output = scipy.stats.ttest_ind(closest_dist_list_pos[pocket][cind],closest_dist_list_neg[pocket][cind])
                pval = ttest_output[1]
                print(pocket,cind,ttest_output)            
    return 


# for each cluster, compile list of experimentally screened compounds that were found to be most similar to cluster
def list_closest_ligands(closest_cluster,cluster_dict):
    closest_ligands = {}
    for pocket, clusters in cluster_dict[protnow].items():
        closest_ligands[pocket] = {}
        clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
        clustall=[]
        for k in range(max([x[1] for x in clusout])+1):
            clustall.append([x[0] for x in clusout if x[1]==k])
        n_clusters=len(clustall)
        for cind,clust in enumerate(clustall,1):
            closest_ligands[pocket][cind] = []
            for exp_id in closest_cluster.keys():
                if closest_cluster[exp_id][pocket][0]==cind:
                    closest_ligands[pocket][cind].append(exp_id)
                    
    return closest_ligands



prot_list_focus = ['S','nsp5','nsp12']

consensus_pockets = {}
all_consensus_residues = {}
consensus_pocket_ligands = {}
consensus_pocket_reslig_pairs = {}
fraction_lig_contacts = {}
fraction_ligand_contacts_matrix_dict = {}
Tdistlist_dict = {}
Cdistlist_dict = {}
LNdistlist_dict = {}
Wdistlist_dict = {}
Wdistmat_dict = {}
Tweight = 0.5
cluster_dict = {}
comout_dict = {}
fp_radius = 2
nBits = 1024
gdccut = '60'


ligs_leaveout = pickle.load(open('ligs_leaveout.p','rb'))
chemtax_dict = pickle.load(open('chemtax_dict.p', 'rb')) 
ligname_dist_dict_notscaled = pickle.load(open('ligname_dist_dict_notscaled.p','rb'))

directory = 'cluster-output-ncov-residues-shortestpath-CCC-15-10-'+gdccut+'-4-0.ligs_8/date_current_resall/'
for protnow in prot_list_focus:
    consensus_pockets = pocket_residues(consensus_pockets,all_consensus_residues,directory,protnow)[0]
    all_consensus_residues = pocket_residues(consensus_pockets,all_consensus_residues,directory,protnow)[1]

directory = './'
filenames = ['CCC.confidence_centroid_contacts.15_10_'+gdccut+'_4_0.ligs_8.nCoV.current.resall']

pocket_ligs = []
for protnow in prot_list_focus:
    consensus_pocket_ligands = filtered_pocket_ligands(directory,filenames,consensus_pockets,protnow,ligs_leaveout)
    consensus_pocket_reslig_pairs = pocket_residue_ligand_pairs(directory,filenames,consensus_pockets,consensus_pocket_ligands,protnow)
    for pocket,ligands in consensus_pocket_ligands[protnow].items():
        for lig in ligands:
            if lig not in pocket_ligs:
                pocket_ligs.append(lig)


Tdist_output = calc_Tanimoto_dist_norm(fp_radius,nBits,chemtax_dict)
Tdistnorm_dict = Tdist_output[0]
pickle.dump(Tdistnorm_dict,open('normalized-Tanimoto-dist.p','wb'))
#Tdistnorm_dict = pickle.load(open('normalized-Tanimoto-dist.p','rb'))

# chemical taxonomy dictionary
chemtax_dict = Tdist_output[1]
pickle.dump(chemtax_dict,open('chemtax_dict_updated.p','wb'))
#chemtax_dict = pickle.load(open('chemtax_dict_updated.p','rb'))

Cdist_output = calc_chemtax_dist_norm()
Cdistnorm_dict = Cdist_output[0]
pickle.dump(Cdistnorm_dict,open('normalized-chemtax-dist.p','wb'))
#Cdistnorm_dict = pickle.load(open('normalized-chemtax-dist.p','rb'))

LNdist_output = calc_ligname_dist_norm(ligname_dist_dict_notscaled)
LNdistnorm_dict = LNdist_output[0]
pickle.dump(LNdistnorm_dict,open('normalized-wordvec-dist.p','wb'))
#LNdistnorm_dict = pickle.load(open('normalized-wordvec-dist.p','rb'))

for protnow in prot_list_focus:
    Tdistlist_dict = get_Tanimoto_dist(Tdistlist_dict,consensus_pocket_ligands,protnow,Tdistnorm_dict)
    Cdistlist_dict = get_chemtax_dist(Cdistlist_dict,consensus_pocket_ligands,protnow,Cdistnorm_dict)
    LNdistlist_dict = get_ligname_dist(LNdistlist_dict,consensus_pocket_ligands,protnow,LNdistnorm_dict)
    
    Wdistlist_dict = weighted_dist(protnow,consensus_pocket_ligands,Tdistlist_dict,Cdistlist_dict,LNdistlist_dict,Wdistlist_dict)
    Wdistmat_dict = get_dist_matrix(Wdistlist_dict,Wdistmat_dict,consensus_pocket_ligands)
    
    cluster_dict = dbscan_cluster_pocket_ligands(cluster_dict,protnow,consensus_pocket_ligands,Wdistmat_dict)
    
    save_ligand_clusters(cluster_dict,consensus_pocket_ligands)
    
    fraction_ligand_contacts_matrix_dict = fraction_cluster_contacts_heatmap(cluster_dict,consensus_pocket_ligands,consensus_pocket_reslig_pairs,protnow,fraction_ligand_contacts_matrix_dict)
    
    silhouette(Wdistmat_dict,cluster_dict,consensus_pocket_ligands)
    
    pocket_mcs(cluster_dict,consensus_pocket_ligands)
    
    cluster_distance_pos = compare_exp_pos_compounds(protnow,smiles_dict,cluster_dict,consensus_pocket_ligands)
    cluster_distance_neg = compare_exp_neg_compounds(protnow,smiles_dict,cluster_dict,consensus_pocket_ligands)
    closest_cluster_pos = find_closest_cluster(cluster_distance_pos,protnow,cluster_dict,consensus_pocket_ligands)[0]
    closest_dist_list_pos = find_closest_cluster(cluster_distance_pos,protnow,cluster_dict,consensus_pocket_ligands)[1]
    closest_cluster_neg = find_closest_cluster(cluster_distance_neg,protnow,cluster_dict,consensus_pocket_ligands)[0]
    closest_dist_list_neg = find_closest_cluster(cluster_distance_neg,protnow,cluster_dict,consensus_pocket_ligands)[1]
    
    #print(closest_dist_list_pos)
    #print(closest_dist_list_neg)
    
    closest_cluster_hist(protnow,closest_dist_list_pos,closest_dist_list_neg,cluster_dict,consensus_pocket_ligands)
    closest_cluster_ttest(protnow,closest_dist_list_pos,closest_dist_list_neg,cluster_dict,consensus_pocket_ligands)
       
    #closest_ligands_pos = list_closest_ligands(closest_cluster_pos,cluster_dict)
    #closest_ligands_neg = list_closest_ligands(closest_cluster_neg,cluster_dict)
    #print(closest_ligands_pos)


            

    
    
    
        
        
   
        

                                                     




