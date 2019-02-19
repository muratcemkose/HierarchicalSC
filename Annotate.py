"""
Created on Mon Feb  19 14:22:31 2019

@author: Murat Cem Köse
"""
import Utils
import numpy as np
import pandas as pd
import scipy

def annotateTree(sc_data, refDataset, annot):
    """ 
    The wrapper function to obtain cell types at all levels.
    
    Parameters
    ----------
    sc_data : DataFrame
        Single cell dataset gene expression matrix.
    refDataset : DataFrame
        Reference dataset gene expression matrix.
    annot : DataFrame
        Annotations for samples in reference dataset.
    Returns
    -------
    final_annotations : DataFrame
        A matrix containing cell type annotations at all levels for all single cells.
    de_dict : Dict
        Dictionary containing DE genes generated while annotating cell types.
    
    """
    intersect = np.intersect1d(refDataset.index.values,sc_data.index)
    sc_data = sc_data.loc[intersect]
    refDataset = refDataset.loc[intersect]
    final_annotations=pd.DataFrame(index=sc_data.columns)
    de_dict={}
    #Column names are checked
    for i in ["level" in i for i in annot.columns]:
        if i == False:
            print("Please indicate annotation columns as levels. Such as 'level1', 'level2'.")
            break

    k = 1
    while "level"+str(k) in annot.columns:
        final_annotations, de = annotateLevel(sc_data, refDataset, annot, k, final_annotations)
        de_dict.update({"level"+str(k) : de})
        k += 1 
    return final_annotations, de_dict

def annotateLevel(sc_data, refDataset, annot, k, final_annotations):
    """ 
    A function to obtain cell types for a certain level.
    
    Parameters
    ----------
    sc_data : DataFrame
        Single cell dataset gene expression matrix.
    refDataset : DataFrame
        Reference dataset gene expression matrix.
    annot : DataFrame
        Annotations for samples in reference dataset.
    k : int
        The number of the level.
    final_annotations : DataFrame
        A matrix containing cell type annotations at all levels for all single cells (an incomplete version).
    Returns
    -------
    final_annotations : DataFrame
        An updated version of the current annotation data matrix.
    de_lvl : Dict
        Dictionary containing DE genes generated at certain level.
    """
    level = "level" + str(k)
    de_lvl={}
    if k == 1:
        de = Utils.getDEgenes(refDataset,annot,level)
        unique_types = [i for i in annot.groupby(level).groups]
        round_res = annotateGroup(sc_data, refDataset, annot, unique_types, de, level)
        final_annotations[level] = round_res
        de_lvl = de
        return final_annotations, de_lvl
    else:
        prev_level = "level"+str(k-1)
        prev_groups = [i for i in final_annotations.groupby(prev_level).groups]
        level_annotations = pd.DataFrame()

        for j in prev_groups:
            sc_cells = final_annotations.groupby(prev_level).groups.get(j)
            ref_cells = annot[annot[prev_level] == j].index
            unique_types = [i for i in annot[annot[prev_level] == j].groupby(level).groups]
            # Should I assign in these two cases???
            if len(unique_types) == 0:
                round_res = pd.DataFrame([float("Nan")]*len(sc_cells),index=sc_cells,columns=["type"])
                level_annotations = pd.concat([level_annotations,round_res],axis=0)
            elif len(unique_types) == 1:
                round_res = pd.DataFrame(unique_types*len(sc_cells),index=sc_cells,columns=["type"])
                level_annotations = pd.concat([level_annotations,round_res],axis=0)
            else:
                de = Utils.getDEgenes(refDataset.loc[:,ref_cells],annot.loc[ref_cells],level)
                de_lvl.update(de)
                round_res = annotateGroup(sc_data.loc[:,sc_cells], refDataset.loc[:,ref_cells], annot.loc[ref_cells],unique_types, de, level)
                level_annotations = pd.concat([level_annotations,round_res],axis=0)
        final_annotations[level] = level_annotations
        return final_annotations, de_lvl

def annotateGroup(sc_data, refDataset, annot, unique_types, de, level):
    """ 
    A function to annotate subtypes of a specific cell type.
    
    Parameters
    ----------
    sc_data : DataFrame
        Single cell dataset gene expression matrix.
    refDataset : DataFrame
        Reference dataset gene expression matrix.
    annot : DataFrame
        Annotations for samples in reference dataset.
    unique_types : List
        A list containing unique cell types (subtypes) belonging to group of interest (maintype).
    de : Dict
        A dictionarycontaining pairwise differentially expressed genes between subgroups.
    level : Str
        A string indicating which level the analysis is taking place.
    Returns
    -------
    round_res : DataFrame
        A matrix, including annotations of the subtypes.
        
    """
    n = 1
    round_res = pd.DataFrame()
    combinations = [(i,j) for i in unique_types  for j in unique_types if i != j]

    while(n < len(unique_types)):
        if n == 1:
            type1 = unique_types[n-1]
            type2 = unique_types[n]
            round_res = paiwiseComparison(sc_data, refDataset, annot, de, level, combinations, type1, type2)
            n+=1
        else:
            round_res = comparisonRound(sc_data, refDataset, annot, unique_types, de, level, combinations, round_res, n)
            n+=1
    return round_res

# This function first defines groups that are obtained until this round. It takes cells with same annotation and compares 
# defined annotation with the next one. Does this for all current annotation types and returns this rounds result.
def comparisonRound(sc_data, refDataset, annot, unique_types, de, level, combinations, round_res, n):
    """ 
    A function to make comparison of correlations between the current annotation and next possible annotation.
    
    Parameters
    ----------
    sc_data : DataFrame
        Single cell dataset gene expression matrix.
    refDataset : DataFrame
        Reference dataset gene expression matrix.
    annot : DataFrame
        Annotations for samples in reference dataset.
    unique_types : List
        A list containing unique cell types (subtypes) belonging to group of interest (maintype).
    de : Dict
        A dictionarycontaining pairwise differentially expressed genes between subgroups.
    level : Str
        A string indicating which level the analysis is taking place.
    combinations : List
        A list containing possible cell type combinations.
    round res : DataFrame
        The results form the previous round
    n : Int
        The number of the current round.
    Returns
    -------
    temp : DataFrame
        A temporary matrix, including annotations of the subtypes for the current round.
        
    """
    groups = [i for i in round_res.groupby("type").groups]
    temp = pd.DataFrame()
    for i in groups:
        type1 = i
        type2 = unique_types[n]
        comparison_res = paiwiseComparison(sc_data.loc[:,round_res.groupby("type").groups.get(i)], refDataset, annot, de, level, combinations, type1, type2)
        temp = pd.concat([temp,comparison_res],axis = 0)
    return temp

def paiwiseComparison(sc_data, refDataset, annot, de, level, combinations, type1, type2):
    """ 
    A function to make pairwise comparison of correlations between cell types.
    
    Parameters
    ----------
    sc_data : DataFrame
        Single cell dataset gene expression matrix.
    refDataset : DataFrame
        Reference dataset gene expression matrix.
    annot : DataFrame
        Annotations for samples in reference dataset.
    de : Dict
        A dictionarycontaining pairwise differentially expressed genes between subgroups.
    level : Str
        A string indicating which level the analysis is taking place.
    combinations : List
        A list containing possible cell type combinations.
    type1 : Str
        First type to compare.
    type2 : Str
        Second type to compare.
    Returns
    -------
    comparison_res : DataFrame
        A temporary matrix, including annotations of the cell types for this comparison.
        
    """
    keys = [i for i in combinations if type1 in i and type2 in i]
    de_merged = []
    de_merged.extend(de.get(keys[0]))
    de_merged.extend(de.get(keys[1]))
    
    cells = []
    cells.extend(annot[annot[level] == type1].index)
    cells.extend(annot[annot[level] == type2].index)
    
    cor = scipy.stats.spearmanr(sc_data.loc[de_merged], refDataset.loc[de_merged, cells])
    cor = pd.DataFrame(cor[0]).iloc[:,0:len(sc_data.columns)][-len(cells):]
    cor.columns = sc_data.columns
    cor[level] = annot[level].loc[cells].values
    comparison_res = pd.DataFrame(cor.groupby(level).quantile(q = 0.8).idxmax(), columns = ["type"])
    return comparison_res

