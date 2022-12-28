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
from rdkit.Chem import rdFMCS, AllChem, Draw, Lipinski, Fragments, Descriptors
from itertools import chain, combinations


# import ligand smiles strings
ligand_file = 'lig_031022.json'

with open(ligand_file) as ff:
    smiles_dict = json.load(ff)
    
    
# check for ligands with bad or missing smiles strings
bad_smiles = []

for key,value in smiles_dict.items():
    if 'smiles_cactus' not in value or Chem.MolFromSmiles(smiles_dict[key]['smiles_cactus']) is None:
        bad_smiles.append(key)
        
good_smiles = []
        
for key,value in smiles_dict.items():
    if key not in bad_smiles:
        good_smiles.append(key)
    

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
                        
    
# calculate Tanimoto distance for pairs of ligands in each pocket
def get_Tanimoto_dist_withmissingvalues(Tdistlist_dict,consensus_pocket_ligands,protnow,Tdistnorm_dict):
    Tdistlist_dict[protnow] = {}
    ligands = consensus_pocket_ligands[protnow][pocket]
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
            else:
                Tdistlist.append('NA')
    Tdistlist_dict[protnow][pocket] = Tdistlist 
        
    return Tdistlist_dict


# calculate normalized chemical taxonomy distance between all pairs of ligands
def calc_chemtax_dist_norm_revised():
    Cdist_dict = {}
    Cdistnorm_dict = {}
    Cdistlist = []
    for i1 in range(0,len(consensus_pocket_ligands[protnow][pocket])):
        lig1 = consensus_pocket_ligands[protnow][pocket][i1]
        if lig1 not in bad_smiles:
            for i2 in range(i1+1,len(consensus_pocket_ligands[protnow][pocket])):
                lig2 = consensus_pocket_ligands[protnow][pocket][i2]
                Csim = 0
                denom = 0
                if lig2 not in bad_smiles:
                    if lig1 in chemtax_dict and lig2 in chemtax_dict:
                        for level in ['kingdom','superclass','class','subclass']:
                            if chemtax_dict[lig1][level]!='None' and chemtax_dict[lig2][level]!='None':
                                denom = denom + 1
                                if chemtax_dict[lig1][level]==chemtax_dict[lig2][level]:
                                    Csim = Csim + 1
                        Csim = float(Csim)/float(denom)
                        Cdist = 1-Csim
                        Cdistlist.append(Cdist)
                        Cdist_dict[(lig1,lig2)] = Cdist 
                    
    Cdistavg = sum(Cdistlist)/float(len(Cdistlist))
    Cdiststd = statistics.stdev(Cdistlist)
    
    print(Cdistavg,Cdiststd,min(Cdistlist),max(Cdistlist))
    
    for ligpair in Cdist_dict.keys():
        Cdistnorm_dict[ligpair] = ((Cdist_dict[ligpair]-Cdistavg)/Cdiststd)+15
     
    return Cdistnorm_dict,Cdistavg,Cdiststd




# calculate distance based on chemical taxonomy for pairs of ligands in each pocket
def get_chemtax_dist(Cdistlist_dict,consensus_pocket_ligands,protnow,Cdistnorm_dict):
    Cdistlist_dict[protnow] = {} 
    ligands = consensus_pocket_ligands[protnow][pocket]
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
    for i1 in range(0,len(consensus_pocket_ligands[protnow][pocket])):
        lig1 = consensus_pocket_ligands[protnow][pocket][i1]
        if lig1 not in bad_smiles:
            for i2 in range(i1+1,len(consensus_pocket_ligands[protnow][pocket])):
                lig2 = consensus_pocket_ligands[protnow][pocket][i2]
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
    
    print(LNdistavg,LNdiststd,min(LNdistlist),max(LNdistlist))
    
    for ligpair in LNdist_dict.keys():
        LNdistnorm_dict[ligpair] = ((LNdist_dict[ligpair]-LNdistavg)/LNdiststd)+15
     
    return LNdistnorm_dict,LNdistavg,LNdiststd


# make distance matrix from ligand name distance dictionary
def get_ligname_dist(LNdistlist_dict,consensus_pocket_ligands,protnow,LNdistnorm_dict):
    LNdistlist_dict[protnow] = {}
    ligands = consensus_pocket_ligands[protnow][pocket]
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
def weighted_dist_revised(protnow,consensus_pocket_ligands,Tdistlist_dict,Cdistlist_dict,LNdistlist_dict,Wdistlist_dict):
    Wdistlist_dict[protnow] = {}
    Wdistlist = []
    ligands = consensus_pocket_ligands[protnow][pocket]
    Wdistlist_dict[protnow][pocket] = []
    if len(Tdistlist_dict[protnow][pocket])!=len(Cdistlist_dict[protnow][pocket]):
        print('distance problem')
    if len(Cdistlist_dict[protnow][pocket])!=len(LNdistlist_dict[protnow][pocket]):
        print('distance problem')
    for ind in range(0,len(Tdistlist_dict[protnow][pocket])):
        Tdist = Tdistlist_dict[protnow][pocket][ind]
        Cdist = Cdistlist_dict[protnow][pocket][ind]
        LNdist = LNdistlist_dict[protnow][pocket][ind]
        if Tdist!='NA' and Cdist!='NA' and LNdist!='NA':
            Wdist = (float(1)/float(3))*Tdist + (float(1)/float(3))*Cdist + (float(1)/float(3))*LNdist
        elif Tdist!='NA' and Cdist!='NA' and LNdist=='NA':
            Wdist = 0.5*Tdist + 0.5*Cdist
        elif Tdist!='NA' and Cdist=='NA' and LNdist!='NA':
            Wdist = 0.5*Tdist + 0.5*LNdist
        elif Tdist!='NA' and Cdist=='NA' and LNdist=='NA':
            Wdist = Tdist
        elif Tdist=='NA' and Cdist!='NA' and LNdist!='NA':
            Wdist = 0.5*Cdist + 0.5*LNdist
        elif Tdist=='NA' and Cdist!='NA' and LNdist=='NA':
            Wdist = Cdist
        elif Tdist=='NA' and Cdist=='NA' and LNdist!='NA':
            Wdist = LNdist
        else:
            print(ind,'no Wdist assignment')                
        Wdistlist_dict[protnow][pocket].append(Wdist) 
        Wdistlist.append(Wdist)
            
    Wdistavg = sum(Wdistlist)/float(len(Wdistlist))
    Wdiststd = statistics.stdev(Wdistlist)
    
    print(Wdistavg,Wdiststd,min(Wdistlist),max(Wdistlist))
            
    return Wdistlist_dict


# get distance matrix from distance list
def get_dist_matrix(Wdistlist_dict,Wdistmat_dict,consensus_pocket_ligands):
    Wdistmat_dict[protnow] = {}
    ligands = consensus_pocket_ligands[protnow][pocket]
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
    ligands = consensus_pocket_ligands[protnow][pocket]
    best_params = {}
    if len(ligands)>2:
        best_params['sscore'] = -1
            
        eps_vec = np.linspace(0.05,15,15*20) 
            
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
                        
                #elif len(ligands)==n_clusters:
                    #print('All ligands clustered separately')
                #elif n_clusters==1:
                    #print('All ligands clustered together')
                #elif n_clusters==0:
                    #print('No ligands clustered')
                
        print(protnow,pocket)
        print(len(ligands))
        print(best_params) 
            
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

            fig = plt.figure()
            fig.set_size_inches(9, 6)
            ax=fig.add_subplot(111)
            
            ax.set_ylim([0, len(cluster_labels) + (n_clusters + 1) * 10])

            silhouette_avg = silhouette_score(dist_matrix, cluster_labels, metric="precomputed", sample_size=None)
            print("There are ",n_clusters," clusters and the average silhouette_score is : ",silhouette_avg)

            sample_silhouette_values = silhouette_samples(dist_matrix, cluster_labels, metric="precomputed")

            y_lower = 10
            for i in range(n_clusters):
                ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

                ith_cluster_silhouette_values.sort()

                size_cluster_i = ith_cluster_silhouette_values.shape[0]
                y_upper = y_lower + size_cluster_i

                color = cm.nipy_spectral(float(i)/n_clusters)
                ax.fill_betweenx(np.arange(y_lower, y_upper),0, ith_cluster_silhouette_values,facecolor=color, edgecolor=color, alpha=0.7)

                y_lower = y_upper + 10  

            plt.title(protnow+', Pocket '+pocket,fontsize=16)
            plt.xlabel("Silhouette coefficient",fontsize=16)
            plt.ylabel("Clusters",fontsize=16)

            ax.axvline(x=silhouette_avg, color="red", linestyle="--")

            ax.set_yticks([])  
            plt.xlim((-0.15,0.25))
            ax.set_xticks([-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25])
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
    contact_substr_chemotypes = {}
    contact_substr_chemotypes[protnow] = {}
    with open('ligand-cluster-key-CCC-15-10-'+gdccut+'-4-0-ligs-8-current-resall-'+protnow+'.csv','w') as f:
        writeCSV = csv.DictWriter(f,fieldnames=['Pocket','Cluster Index','Ligands in Cluster'])
        writeCSV.writeheader()
        for pocket, clusters in cluster_dict[protnow].items():
            contact_substr_chemotypes[protnow][pocket] = {}
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
                    if lig in chemtax_dict:
                        if chemtax_dict[lig]['subclass'] not in subclasses:
                            subclasses.append(chemtax_dict[lig]['subclass'])

                if len(clust)>=2:
                    molecules = []
                    for lig in clust:
                        if lig not in bad_smiles:
                            molecules.append(Chem.MolFromSmiles(smiles_dict[lig]['smiles_cactus']))
                    if len(molecules)>1:
                        mcs = Chem.rdFMCS.FindMCS(molecules)
                        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                        contact_substr_chemotypes[protnow][pocket][cind] = mcs_mol
                        if mcs.numAtoms>=4:
                            mcs_smiles = Chem.MolToSmiles(Chem.MolFromSmarts(mcs.smartsString))
                            mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                            mcs_coords = AllChem.Compute2DCoords(mcs_mol)
                            image_file = 'images/CCC-15-10-'+gdccut+'-4-0-ligs-8-current-resall/'\
                            +protnow+'_pocket'+pocket+'_cluster'+str(cind)+'_mcs.png'
                            Draw.MolToFile(mcs_mol,image_file)
    
    return contact_substr_chemotypes


# save info about each ligand cluster in each pocket
def pocket_cluster_info(cluster_dict,consensus_pocket_ligands):
    with open('ligand-cluster-info-CCC-15-10-'+gdccut+'-4-0-ligs-8-current-resall-'+protnow+'.csv','w') as f:
        writeCSV = csv.DictWriter(f,fieldnames=['Pocket','Cluster Index','Number of Ligands','Class','Number of Subclasses'])
        writeCSV.writeheader()
        
        for pocket, clusters in cluster_dict[protnow].items():
            clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
            clustall=[]
            for k in range(max([x[1] for x in clusout])+1):
                clustall.append([x[0] for x in clusout if x[1]==k])
            n_clusters=len(clustall)
            for cind,clust in enumerate(clustall,1):
                subclasses = []
                for lig in clust:
                    if lig in chemtax_dict:
                        if chemtax_dict[lig]['subclass'] not in subclasses:
                            subclasses.append(chemtax_dict[lig]['subclass'])
                    else:
                        if 'NA' not in subclasses:
                            subclasses.append('NA')
                if len(clust)>1:
                    print('Cluster',cind)
                    print('# Ligands',len(clust))
                    print('# Subclasses',len(subclasses))
                    
                    newrow = {'Pocket':pocket, 'Cluster Index':cind, 'Number of Ligands':len(clust),'Class':chemtax_dict[clust[0]]['class'], 'Number of Subclasses':len(subclasses)}
                    writeCSV.writerow(newrow)

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
    with open('ligand-cluster-key-CCC-15-10-'+gdccut+'-4-0-ligs-8-current-resall-'+protnow+'-contactsubstr-withsingletons.csv','w') as f:
        writeCSV = csv.DictWriter(f,fieldnames=['Pocket','Cluster Index','Ligands in Cluster'])
        writeCSV.writeheader()
        for pocket, clusters in cluster_dict[protnow].items():
            clusout=[(x,clusters[k][0]) for k,x in enumerate(consensus_pocket_ligands[protnow][pocket])]
            clustall=[]
            for k in range(-1,max([x[1] for x in clusout])+1):
                clustall.append([x[0] for x in clusout if x[1]==k])
            n_clusters=len(clustall)
            for cind,clust in enumerate(clustall,0):
                newrow = {'Pocket': pocket, 'Cluster Index': cind, 'Ligands in Cluster': clust}
                writeCSV.writerow(newrow)


def make_lig_pdbfile(pdbdirectory,lig):
    for root, dirs, files in os.walk(pdbdirectory+lig):
        for name in files:
            usefile = pdbdirectory+lig+'/'+name
            template = name.split('.')[3]
            
            break
    
    with open(usefile,'r') as origfile:
        atomlist = []
        lines = origfile.readlines()
        with open('./ligand_pdb_files/ligand.'+lig+'.'+template+'.pdb','w') as newfile:
            for line in lines:
                if 'HETATM' not in line.split()[0] and line.split()[0]!='ATOM' and line.split()[0]!='FORMUL':
                    newfile.write(line)
                elif 'HETATM' in line.split()[0] and line.split()[0]!='HETATM':
                    if line.split()[2].strip()==lig and line.split()[1].strip() not in atomlist:
                        newfile.write(line) 
                        atomlist.append(line.split()[1].strip())
                elif 'HETATM' in line.split()[0] and line.split()[0]=='HETATM':
                    if line.split()[3].strip()==lig and line.split()[2].strip() not in atomlist:
                        newfile.write(line) 
                        atomlist.append(line.split()[2].strip())
    return template


# make line plot showing percent of experimental compounds assigned to each chemotype that are active, going from strictest to least strict 
def line_plot_exp_chemotypes_bylevel(protnow,pos_exp_cmpd_assigned_chemotype,neg_exp_cmpd_assigned_chemotype,ordered_chemotypes,number_cmpds_per_chemotype,total_number_cmpds_good_smiles,tot_num_pos_exp_cmpds,tot_num_neg_exp_cmpds):
    num_pos_cmpds = tot_num_pos_exp_cmpds
    num_neg_cmpds = tot_num_neg_exp_cmpds
    
    chemotype_counts = {}
    for chemotype in ordered_chemotypes:
        chemotype_counts[chemotype] = {'pos':0, 'neg':0}
        
        for pcmpd in pos_exp_cmpd_assigned_chemotype:
            if pos_exp_cmpd_assigned_chemotype[pcmpd]==chemotype:
                chemotype_counts[chemotype]['pos'] = chemotype_counts[chemotype]['pos'] + 1
                
        for ncmpd in neg_exp_cmpd_assigned_chemotype:
            if neg_exp_cmpd_assigned_chemotype[ncmpd]==chemotype:
                chemotype_counts[chemotype]['neg'] = chemotype_counts[chemotype]['neg'] + 1
                   
    bar_counts = []
    for chemotype in chemotype_counts:
        bar_counts.append((chemotype,chemotype_counts[chemotype]['pos'],chemotype_counts[chemotype]['neg'],number_cmpds_per_chemotype[chemotype]))
        
    chemos = list(np.arange(1,len(bar_counts)+1))
    x_axis = np.arange(len(chemos)) 
    pos_counts = [bc[1] for bc in bar_counts]
    neg_counts = [bc[2] for bc in bar_counts]
    percent_positive = [float(pc)/float(pc+nc) for (pc,nc) in zip(pos_counts,neg_counts) if (pc+nc)!=0]
    sum_counts = [pc+nc for (pc,nc) in zip(pos_counts,neg_counts)]
    x_axis_defined = [k+1 for k,v in enumerate(sum_counts) if v!=0]
    
    plt.scatter(x_axis_defined,percent_positive)
    plt.axhline(y = tot_num_pos_exp_cmpds/(tot_num_pos_exp_cmpds+tot_num_neg_exp_cmpds), color = 'r', linestyle = '--')
    
    plt.xticks(x_axis_defined,fontsize=14)
    plt.xlabel("Chemotype",fontsize=16)
    plt.yticks([0,0.2,0.4,0.6,0.8,1.0],fontsize=14)
    plt.ylabel("Fraction of experimentally tested compounds \nassigned to chemotype that are active",fontsize=16)  
    plt.show()
                      
    return


# for each ligand, create substructure including only atoms in contact with protein
def extract_substructure_contact_atoms(protnow,pocket,pdbspheres_lig_atom_contacts,pdbdirectory,consensus_pocket_ligands):
    contact_atom_substructures = {}
    ligands_rdkit_pdb_issues = []
    
    for ligand in pdbspheres_lig_atom_contacts[protnow][pocket]: 
        if len(ligand)<=3 and ligand in good_smiles:
            # make ligand PDB file from PDB file of protein-ligand complex
            template = make_lig_pdbfile(pdbdirectory,ligand)
            ligmol = Chem.MolFromPDBFile('./ligand_pdb_files/ligand.'+ligand+'.'+template+'.pdb')
            
            try:
                atom_indlab = []
   
                for i, atom in enumerate(ligmol.GetAtoms()):
                    aind = atom.GetIdx()
                    at = ligmol.GetAtomWithIdx(i)
                    atlab = at.GetPDBResidueInfo().GetName()
                    atom_indlab.append((aind,atlab.strip()))
                
                edit_ligmol = ligmol
                
                for ail in reversed(atom_indlab):
                    aind = ail[0]
                    alab = ail[1]
                    edit_ligmol_temp = edit_ligmol
                    if alab not in pdbspheres_lig_atom_contacts[protnow][pocket][ligand]:
                        edit_ligmol = Chem.RWMol(edit_ligmol_temp)
                        edit_ligmol.RemoveAtom(aind)
                        try:
                            AllChem.GetMorganFingerprintAsBitVect(edit_ligmol,fp_radius,nBits)
                        except:
                            edit_ligmol = edit_ligmol_temp
                    else:
                        edit_ligmol = edit_ligmol_temp
                    
                contact_atom_substructures[ligand] = edit_ligmol
                
            except AttributeError:
                print(ligand,template,'problem with PDB file')
                ligands_rdkit_pdb_issues.append(ligand)
    
    for ligand in ligands_rdkit_pdb_issues:
        consensus_pocket_ligands[protnow][pocket].remove(ligand)
                                   
    return contact_atom_substructures, consensus_pocket_ligands


# calculate normalized Tanimoto distances between ligands based on their contact substructures
def calc_Tanimoto_dist_norm_substr(fp_radius,nBits,contact_atom_substructures):
    Tdist_dict = {}
    Tdistnorm_dict = {}
    Tdistlist = []
    ligs_missing_contacts = []
    for i1 in range(0,len(consensus_pocket_ligands[protnow][pocket])):
        lig1 = consensus_pocket_ligands[protnow][pocket][i1]
        if lig1 not in bad_smiles and lig1 in contact_atom_substructures:    
            m1 = contact_atom_substructures[lig1]
            fp1 = AllChem.GetMorganFingerprintAsBitVect(m1,fp_radius,nBits)
            for i2 in range(i1+1,len(consensus_pocket_ligands[protnow][pocket])):
                lig2 = consensus_pocket_ligands[protnow][pocket][i2]
                if lig2 not in bad_smiles and lig2 in contact_atom_substructures:
                    m2 = contact_atom_substructures[lig2]
                    fp2 = AllChem.GetMorganFingerprintAsBitVect(m2,fp_radius,nBits)
                    Tsim = DataStructs.FingerprintSimilarity(fp1,fp2)
                    Tdist = 1-Tsim
                    Tdistlist.append(Tdist)
                    Tdist_dict[(lig1,lig2)] = Tdist
        
        elif lig1 not in bad_smiles and lig1 not in contact_atom_substructures:
            ligs_missing_contacts.append(lig1)
                       
    Tdistavg = sum(Tdistlist)/float(len(Tdistlist))
    Tdiststd = statistics.stdev(Tdistlist)
    
    print(Tdistavg,Tdiststd,min(Tdistlist),max(Tdistlist))
    
    for ligpair in Tdist_dict.keys():
        Tdistnorm_dict[ligpair] = ((Tdist_dict[ligpair]-Tdistavg)/Tdiststd)+15
     
    return Tdistnorm_dict,Tdistavg,Tdiststd,ligs_missing_contacts


# compare experimentally screened compounds with ligand chemotypes
def assign_exp_compounds_to_contactsubstr_chemotype(protnow,datafile,smiles_header,cmpd_id_header,contact_substr_chemotypes):
    exp_cmpd_assigned_chemotype = {}
    chemotype_chemicalweights = {}
    all_unique_chemotypes = []
    
    for pocket in contact_substr_chemotypes[protnow]:
        for cind in contact_substr_chemotypes[protnow][pocket]:
            chemotype = contact_substr_chemotypes[protnow][pocket][cind]
            mw = Descriptors.HeavyAtomMolWt(chemotype)
            chemotype_chemicalweights[chemotype] = mw
            if chemotype not in all_unique_chemotypes:
                all_unique_chemotypes.append(chemotype)
    
    ordered_chemotypes = sorted(all_unique_chemotypes, key=lambda chemo: chemotype_chemicalweights[chemo], reverse=True)  
     
    with open(datafile,'r') as smiles_file:
        readCSV = csv.DictReader(smiles_file)
        tot_num_exp_cmpds = 0
        for row in readCSV:
            tot_num_exp_cmpds = tot_num_exp_cmpds + 1
            exp_smiles = str(row[smiles_header])
            exp_id = row[cmpd_id_header]
            exp_mol = Chem.MolFromSmiles(exp_smiles)
            
            for chemotype in ordered_chemotypes:
                substr_check = exp_mol.HasSubstructMatch(chemotype)
                if substr_check==True:
                    exp_cmpd_assigned_chemotype[exp_id] = chemotype
                    break
                        
    return exp_cmpd_assigned_chemotype, ordered_chemotypes, tot_num_exp_cmpds


# assign all ligands in file to chemotype; determine order of strictest to least strict chemotype based on number of ligands assigned to each
def assign_all_ligands_to_contactsubstr_chemotype(protnow,contact_substr_chemotypes,smiles_dict,good_smiles):
    exp_cmpd_assigned_chemotype = {} 
    chemotype_chemicalweights = {}
    all_unique_chemotypes = []
    
    all_unique_chemotypes = []
    for pocket in contact_substr_chemotypes[protnow]:
        for cind in contact_substr_chemotypes[protnow][pocket]:
            chemotype = contact_substr_chemotypes[protnow][pocket][cind]
            mw = Descriptors.HeavyAtomMolWt(chemotype)
            chemotype_chemicalweights[chemotype] = mw
            if chemotype not in all_unique_chemotypes:
                all_unique_chemotypes.append(chemotype)
                
    ordered_chemotypes = sorted(all_unique_chemotypes, key=lambda chemo: chemotype_chemicalweights[chemo], reverse=True) 
    
    for cmpd in good_smiles:
        smiles = smiles_dict[cmpd]['smiles_cactus']
        exp_mol = Chem.MolFromSmiles(smiles)
        exp_id = cmpd
            
        for chemotype in ordered_chemotypes:
            substr_check = exp_mol.HasSubstructMatch(chemotype)
            if substr_check==True:
                exp_cmpd_assigned_chemotype[exp_id] = chemotype
                break

    total_number_cmpds_good_smiles = len(good_smiles)
    
    number_cmpds_per_chemotype = {}
    for chemotype in all_unique_chemotypes:
        number_cmpds_per_chemotype[chemotype] = 0
        for exp_id in exp_cmpd_assigned_chemotype:
            if exp_cmpd_assigned_chemotype[exp_id]==chemotype:
                number_cmpds_per_chemotype[chemotype] = number_cmpds_per_chemotype[chemotype] + 1
    
    ordered_chemotypes_by_prevalence = sorted(all_unique_chemotypes, key=lambda chemo: number_cmpds_per_chemotype[chemo], reverse=False)  
    
    chemotype_compounds = {}
    for chemotype in ordered_chemotypes:
        chemotype_compounds[chemotype] = []
    
    for exp_id in exp_cmpd_assigned_chemotype:
        for chemotype in ordered_chemotypes:
            if exp_cmpd_assigned_chemotype[exp_id]==chemotype:
                chemotype_compounds[chemotype].append(exp_id)
                        
    return exp_cmpd_assigned_chemotype, ordered_chemotypes, number_cmpds_per_chemotype, ordered_chemotypes_by_prevalence, total_number_cmpds_good_smiles, chemotype_compounds





prot_list_focus = ['nsp5']


consensus_pockets = {}
all_consensus_residues = {}
consensus_pocket_ligands = {}
consensus_pocket_reslig_pairs = {}
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
chemtax_dict = pickle.load(open('chemtax_dict_221018_augmented.p', 'rb')) 
ligname_dist_dict_notscaled = pickle.load(open('ligname_dist_dict_notscaled.p','rb'))
pdbspheres_lig_atom_contacts = pickle.load(open('pdbspheres_lig_atom_contacts_221111_nsp5.p','rb'))

directory = 'cluster-output-ncov-residues-shortestpath-CCC-15-10-'+gdccut+'-4-0.ligs_8/date_current_resall/'
for protnow in prot_list_focus:
    consensus_pockets = pocket_residues(consensus_pockets,all_consensus_residues,directory,protnow)[0]
    all_consensus_residues = pocket_residues(consensus_pockets,all_consensus_residues,directory,protnow)[1]

directory = './'
filenames = ['CCC.confidence_centroid_contacts.15_10_'+gdccut+'_4_0.ligs_8.nCoV.current.resall']

for protnow in prot_list_focus:
    consensus_pocket_ligands = filtered_pocket_ligands(directory,filenames,consensus_pockets,protnow,ligs_leaveout)
    consensus_pocket_reslig_pairs = pocket_residue_ligand_pairs(directory,filenames,consensus_pockets,consensus_pocket_ligands,protnow)               
            
pdbdirectory = '/usr/workspace/kpath/BindingSites/PREDICT_binding_sites/PDB_SPHERES/Pdb_HET.SPHERES_12/'

for protnow in prot_list_focus:
    for pocket in consensus_pocket_ligands[protnow]:
        extract_substructures_output = extract_substructure_contact_atoms(protnow,pocket,pdbspheres_lig_atom_contacts,pdbdirectory,consensus_pocket_ligands)
        contact_atom_substructures = extract_substructures_output[0]
        consensus_pocket_ligands = extract_substructures_output[1]
        Tdist_substr_output = calc_Tanimoto_dist_norm_substr(fp_radius,nBits,contact_atom_substructures)
        Tdistnorm_dict = Tdist_substr_output[0]
        ligs_missing_contacts = Tdist_substr_output[3]

        Cdist_output = calc_chemtax_dist_norm_revised()
        Cdistnorm_dict = Cdist_output[0]

        LNdist_output = calc_ligname_dist_norm(ligname_dist_dict_notscaled)
        LNdistnorm_dict = LNdist_output[0]

        Tdistlist_dict = get_Tanimoto_dist_withmissingvalues(Tdistlist_dict,consensus_pocket_ligands,protnow,Tdistnorm_dict)
        Cdistlist_dict = get_chemtax_dist(Cdistlist_dict,consensus_pocket_ligands,protnow,Cdistnorm_dict)
        LNdistlist_dict = get_ligname_dist(LNdistlist_dict,consensus_pocket_ligands,protnow,LNdistnorm_dict)
    
        Wdistlist_dict = weighted_dist_revised(protnow,consensus_pocket_ligands,Tdistlist_dict,Cdistlist_dict,LNdistlist_dict,Wdistlist_dict)
        Wdistmat_dict = get_dist_matrix(Wdistlist_dict,Wdistmat_dict,consensus_pocket_ligands)
    
        cluster_dict = dbscan_cluster_pocket_ligands(cluster_dict,protnow,consensus_pocket_ligands,Wdistmat_dict)
    
        pocket_cluster_info(cluster_dict,consensus_pocket_ligands)
    
        save_ligand_clusters(cluster_dict,consensus_pocket_ligands)
    
        fraction_ligand_contacts_matrix_dict = fraction_cluster_contacts_heatmap(cluster_dict,consensus_pocket_ligands,consensus_pocket_reslig_pairs,protnow,fraction_ligand_contacts_matrix_dict)
    
        silhouette(Wdistmat_dict,cluster_dict,consensus_pocket_ligands)
    
        contact_substr_chemotypes = pocket_mcs(cluster_dict,consensus_pocket_ligands)
    
        if protnow=='nsp5':
            assign_all_ligands_output = assign_all_ligands_to_contactsubstr_chemotype(protnow,contact_substr_chemotypes,smiles_dict,good_smiles)
            ordered_chemotypes = assign_all_ligands_output[1]
            number_cmpds_per_chemotype = assign_all_ligands_output[2]
            ordered_chemotypes_by_prevalence = assign_all_ligands_output[3]
            total_number_cmpds_good_smiles = assign_all_ligands_output[4]
            chemotype_compounds = assign_all_ligands_output[5]
        
            pos_compounds_output = assign_exp_compounds_to_contactsubstr_chemotype(protnow,'./protease_A0A6H1PK90_scores_with_labels_positivesonly.csv','SMILES','Compound ID',contact_substr_chemotypes)
            pos_exp_cmpd_assigned_chemotype = pos_compounds_output[0]
            tot_num_pos_exp_cmpds = pos_compounds_output[2]
            neg_compounds_output = assign_exp_compounds_to_contactsubstr_chemotype(protnow,'./protease_A0A6H1PK90_scores_with_labels_negativesonly.csv','SMILES','Compound ID',contact_substr_chemotypes)
            neg_exp_cmpd_assigned_chemotype = neg_compounds_output[0]
            tot_num_neg_exp_cmpds = neg_compounds_output[2]
            
            line_plot_exp_chemotypes_bylevel(protnow,pos_exp_cmpd_assigned_chemotype,neg_exp_cmpd_assigned_chemotype,ordered_chemotypes,number_cmpds_per_chemotype,total_number_cmpds_good_smiles,tot_num_pos_exp_cmpds,tot_num_neg_exp_cmpds)
    
            sars2_drug_compounds_output = assign_exp_compounds_to_contactsubstr_chemotype(protnow,'sars-cov2-pdb-drugs.csv','SMILES','Drug PDB ID',contact_substr_chemotypes)
            sars2_drug_cmpd_assigned_chemotype = sars2_drug_compounds_output[0]
            #print(sars2_drug_cmpd_assigned_chemotype)
    
            hepC_drug_compounds_output = assign_exp_compounds_to_contactsubstr_chemotype(protnow,'hepC-pdb-drugs.csv','SMILES','Drug PDB ID',contact_substr_chemotypes)
            hepC_drug_cmpd_assigned_chemotype = hepC_drug_compounds_output[0]
            #print(hepC_drug_cmpd_assigned_chemotype)
            
            nsp5_inhib_compounds_output = assign_exp_compounds_to_contactsubstr_chemotype(protnow,'nsp5_inhibitors.csv','SMILES','Drug PDB ID',contact_substr_chemotypes)
            nsp5_inhib_cmpd_assigned_chemotype = nsp5_inhib_compounds_output[0]
            #print(nsp5_inhib_cmpd_assigned_chemotype)
        

                                                 

