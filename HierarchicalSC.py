"""
Created on Mon Feb  19 14:22:31 2019

@author: Murat Cem Köse
"""

import numpy as np
import pandas as pd
import Utils
import Annotate

class HierarchicalSC:
    def __init__(self, scData, refDataset, refAnnot):
        """
        Contructor function for HierarchicalSC class.
    
        Parameters
        ----------
        scData : DataFrame
            Single cell data matrix. Columns representing sample names, rows representing gene symbols.
            
        refDataset : DataFrame
            The reference dataset gene expression matrix. Columns representing sample names, rows representing gene symbols.
            
        refAnnot : DataFrame
            Annotations for samples in reference dataset. Columns of annotations should be cell type levels. Additionally, 
            higher levels should contain lower levels bound with ':'. Example structure; level1 (including B-cells),
            level2 (including B-cells:Naive)
            
        """
        
        self.sc_data = scData
        self.refDataset = refDataset.astype(float)
        self.refAnnot = refAnnot
        
    def getCellTypes(self):
        """
        Annotates single cell types at each level and adds the result to the object.
        
        """
        sc_data = Utils.convertAnnDataToDf(self.sc_data)
        try:
            self.sc_annot, self.de_dict = Annotate.annotateTree(sc_data, self.refDataset, self.refAnnot)
        except:
            print("Columns of annotations should be cell type levels. Additionally, higher levels should contain lower levels bound with ':'. Example structure; level1 (including B-cells), level2 (including B-cells:Naive)")
            
    def getMarkerGenes(self):
        """
        Finds marker genes and adds the result to the object.
        
        """
        try:
            self.marker_genes = Utils.findMarkers(self.sc_annot, self.de_dict)
        except:
            print("Please run getCellTypes first to get cell annotations. This step is needed for marker gene finding.")
            
    def writeMarkerGenes(self, location):
        """
        Writes marker genes to an excel file. Each sheet indicates different levels.
        
        Parameters
        ----------
        location : Str
            Location of the file to be written. Example location, './new_analysis/marker_genes/'.
        
        """
        try:
            writer = pd.ExcelWriter(location+"marker_genes.xlsx", engine="xlsxwriter")
            for key in self.marker_genes:
                self.marker_genes.get(key).to_excel(writer, sheet_name=key)
            writer.save()
        except:
            print("Please run getMarkerGenes first to get marker genes. This step is needed to write them to excel.")
            
    def getTreePlot(self):
        """
        Plots tree representation of single cell annotations.
        
        """
        try:
            Utils.plotTree(self.sc_annot)
        except:
            print("Please run getCellTypes first to get cell annotations. This step is needed for plotting.")

