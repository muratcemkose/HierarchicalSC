
"""
Created on Mon Dec  3 14:45:33 2018

@author: Murat Cem Köse
"""

import scanpy.api as sc
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import collections

def getDEgenes(refDataset,annot,level,n=None):
    """ 
    Creates a dictionary for differentially highly expressed genes for all pairwise cell types in the a reference data set.
    
    Parameters
    ----------
    refDataset : DataFrame
        Reference dataset gene expression matrix.
        
    annot : DataFrame
        Annotations for samples in reference dataset.
        
    level : String
        Level of annotation. Should be one of the columns indicated in reference annotation dataset.
        
    n : Int
        Number of top differential genes to select.
        
    Returns
    -------
    deGenes : dict with multiple index
        Dictionary containing differentially highly expressed genes for each combination of cell types. 
    """
    
    if n is None:
        n = int(500*np.power(2/3,np.log2(len([i for i in annot.groupby(level).groups]))))

    types=[i for i in annot.groupby(level).groups]
    median=refDataset.groupby(annot[level].values,axis=1).apply(np.median,axis=1)
    deGenes={}
    [deGenes.update({(i,j):median[i]-median[j]}) for i in median.index  for j in median.index if i!=j]
    for i in deGenes.keys():
        deGenes[i]=refDataset.iloc[deGenes.get(i).argsort()[-n:]].index.values.tolist()
    return deGenes
    
def readSCData(path,min_genes):
    """ 
    Reads, precesses and returns single cell data. 
    
    Parameters
    ----------
    path : Str
        Diractory path, the location of single cell data.
        
    min_genes : Int
        The minimum number of genes for a cell to have in order to participate the analysis.
        
    Returns
    -------
    scData : AnnData
        Single cell data. 
        
    """
    result = sc.read(path + 'matrix.mtx').transpose() #, cache=True
    result.var_names = np.genfromtxt(path + 'genes.tsv', dtype=str)[:, 1]
    result.obs_names = np.genfromtxt(path + 'barcodes.tsv', dtype=str)
    result.var_names_make_unique()
    result.obs['n_counts'] = np.sum(result.X, axis=1).A1
    sc.pp.filter_cells(result, min_genes=min_genes)
    
    return result

def convertAnnDataToDf(scData):
    """ 
    Converts AnnData object obtained from scanpy into a pandas dataframe.
    
    Parameters
    ----------
    scData : AnnData
        Single cell data. 
        
    Returns
    -------
    sc_data : DataFrame
        Single cell data as pandas dataframe. 
        
    """
    try:
        result = pd.DataFrame(scData.X.toarray()) # If data is 10x data
    except:
        result = pd.DataFrame(scData.X[:]) # If data is Digital Gene expression matrix
        
    result.index = scData.obs_names.values
    result.columns = scData.var_names.values
    return result.T

def findMarkers(sc_annot, de_dict):
    """ 
    Finds markers genes for each cell type at each level.
    
    Parameters
    ----------
    sc_annot : DataFrame
        Single cell annotations. 
    de_dict : Dict
        Dictionary containing DE genes generated while annotating cell types.
    Returns
    -------
    marker_genes : Dict
        A dictionary including marker genes for cell types based on level. 
        
    """
    marker_genes = {}
    k = 1
    while "level"+str(k) in de_dict.keys():
        level = "level"+str(k) 
        de = de_dict.get(level)
        unique_types = [i for i in sc_annot.groupby(level).groups]
        if k >1:
            combinations = [(i,j) for i in unique_types  for j in unique_types if i != j and i[:i.rfind(":")] == j[:j.rfind(":")]]
        else:
            combinations = [(i,j) for i in unique_types  for j in unique_types if i != j]
        markers_lvl = pd.DataFrame()
        for type1 in unique_types:
            keys = [i for i in combinations if type1 in i]
            de_all = []
            for key in keys:
                try:
                    de_all.extend(de.get(key))
                except:
                    de_all=de_all
            type1_markers = pd.DataFrame(collections.Counter(de_all),index=[type1]).T.sort_values(by=type1, ascending = False)/(len(keys)/2) 
            markers_lvl = pd.concat([markers_lvl, type1_markers[type1_markers[type1] == 1]],axis =1)
        markers_lvl = markers_lvl.replace(float("NaN"),False).replace(1,True)
        marker_genes.update({level : markers_lvl})
        k += 1
    return marker_genes


#this function will be revised
def plotTree(final_annotations):
    """ 
    Plots tree representation of the cell type annotations.
    
    Parameters
    ----------
    final_annotations : DataFrame
        Single cell annotations. 
        
    """
    G = nx.DiGraph()
    pos = {}
    G.add_node("Cells")
    location = 1
    step = 0
    pos_first = 0
    pos_last = 0
    for level in final_annotations.columns:
        groups = [i for i in final_annotations.groupby(level).groups]
        groups.sort()
        if level[-1] == "1":
            for i in groups:
                G.add_node(i)
                pos.update({i:(location,1)})
                location = location + 12
                G.add_edge("Cells",i)
            pos.update({"Cells":((pos.get(groups[0])[0] + pos.get(groups[-1])[0])/2,2)})
            pos_first = pos.get(groups[0])[0]
            pos_last = pos.get(groups[-1])[0]
        else:
            step = (pos_first + pos_last+12) / (len(groups) + 1)
            location = pos_first-6
            for i in groups:
                l1 = i.split(":")[0] #these designed for level2!!!
                l2 = i.split(":")[1]
                G.add_node(l2)
                G.add_edge(l1, l2)
                pos.update({l2:(location,2-int(level[-1]))})
                location = location + step
    tree=nx.dfs_tree(G)
    plt.figure(1,figsize=(40,10),facecolor="white")
    plt.axis('off')
    nx.draw_networkx(tree,pos,node_color="w",font_size=12)
    plt.show()