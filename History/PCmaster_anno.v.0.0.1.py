# -*- coding: utf-8 -*-
# Some references:
# D2l. https://zh-v2.d2l.ai/
# Pandas. https://pandas.pydata.org/docs/index.html
# PyTorch. https://pytorch.org/
# Scikit-learn. https://scikit-learn.org/stable/

import sys
import itertools
import statistics
import shutil
import random
import math
from functools import reduce
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import gc
import os
import psutil
import inspect
import re
import pickle
import torch
import torchvision.models
import torch.nn.functional
import datetime
from torch import nn
from d2l import torch as d2l
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from torch.utils.data import DataLoader, Dataset, SequentialSampler, RandomSampler
from PIL import Image
from sklearn.preprocessing import MinMaxScaler
import time
from IPython.core.interactiveshell import InteractiveShell 
from torchvision import transforms  
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from collections import Counter

InteractiveShell.ast_node_interactivity = "all"

now = datetime.datetime.now()
current_time = now.strftime("%Y-%m-%d %H:%M:%S")
print(current_time)

gc.collect()

available_memory = psutil.virtual_memory().available

available_memory_gb = available_memory / (1024 * 1024 * 1024)

print(f"available_memory: {available_memory_gb} GB")

used_memory = psutil.virtual_memory().used

used_memory_gb = used_memory / (1024 * 1024 * 1024)

print(f"used_memory: {used_memory_gb} GB")


class PCmaster_anno_0_Error(Exception):
    pass

class trainDataset(Dataset):  
    def __init__(self,train_data,train_label):
        self.Data = torch.from_numpy(np.float32(np.array(train_data)))
        self.Data = torch.reshape(self.Data, (train_data.shape[0], 224, 224))
        self.Data = torch.unsqueeze(self.Data, 1)
        self.Data = torch.cat((self.Data, self.Data, self.Data), dim=1)
        self.Label = torch.from_numpy(np.int64(np.array(train_label)).reshape(1,-1)[0])
    def __getitem__(self, index):
        x = self.Data[index]
        y = self.Label[index]
        return x, y  
    def __len__(self):
        return len(self.Data)
            
class validDataset(Dataset):  
    def __init__(self,valid_data,valid_label):
        self.Data = torch.from_numpy(np.float32(np.array(valid_data)))
        self.Data = torch.reshape(self.Data, (valid_data.shape[0], 224, 224))
        self.Data = torch.unsqueeze(self.Data, 1)
        self.Data = torch.cat((self.Data, self.Data, self.Data), dim=1)
        self.Label = torch.from_numpy(np.int64(np.array(valid_label)).reshape(1,-1)[0])
    def __getitem__(self, index):
        x = self.Data[index]
        y = self.Label[index]
        return x, y  
    def __len__(self):
        return len(self.Data)

class CustomDataset(Dataset):  
    def __init__(self, root_dir, transform=None):  
        self.root_dir = root_dir  
        self.classes = sorted(os.listdir(root_dir))  
        self.class_to_idx = {cls_name: i for i, cls_name in enumerate(self.classes)}  
        self.samples = self._make_dataset()  
        self.transform = transforms.Compose([  
            transforms.Resize((224, 224)),  # 224x224  
            transforms.ToTensor()  # to tensor
        ])   
  
    def _make_dataset(self):  
        samples = []  
        for class_name in self.classes:  
            class_dir = os.path.join(self.root_dir, class_name)  
            if not os.path.isdir(class_dir):  
                continue  
            for file_name in os.listdir(class_dir):  
                file_path = os.path.join(class_dir, file_name)  
                samples.append((file_path, self.class_to_idx[class_name]))  
        return samples  
  
    def __len__(self):  
        return len(self.samples)  
  
    def __getitem__(self, idx):  
        file_path, label = self.samples[idx]  
        image = Image.open(file_path)  
        if self.transform:  
            image = self.transform(image)  
        return image, label  

class PCmaster_anno_0():
    result_file = None
    shape = None
    
    obj = None 
    
    species = None
    
    organ = None
    dd = None  # dense dataframe, pandas object
    genes = None
    cells = None
    cellnamelabel = None
    
    objs = []
    
    batch_species = []
    
    batch_organs = []
    markergenepd_head5 = None
    markergenepd_head10 = None
    markergenepd_head15 = None
    markergenepd_head20 = None
    markergenepd_head = None
    
    marker_grade_1 = None
    anno_re_using_rank_1_markers = None
    anno_re_using_rank_1and2_markers = None
    dir_gene_species = None
    dir_gene_organ_celltype = None
    dir_gene_score = None
    anno_re_using_marker_scores = None
    
    # deep learning
    ref_obj = None
    test_obj = None
    num_types = None
    mapping_1 = None
    mapping_2 = None
    model_weight = None
    accuracys = None
    precisions = None
    recalls = None
    predicted_results = None
    dl_cellnames = None
    original_celltypes = None
    train_accuracy = None
    train_macro_f1_score = None

    def memory_collect_0(self):
        gc.collect()

        available_memory = psutil.virtual_memory().available

        available_memory_gb = available_memory / (1024 * 1024 * 1024)
        
        print(f"max_memory: {available_memory_gb} GB")
        
        used_memory = psutil.virtual_memory().used
        
        used_memory_gb = used_memory / (1024 * 1024 * 1024)
        
        print(f"used_memory: {used_memory_gb} GB")

    def set_result_file_0(self, result_file=None):
        if result_file is not None:
            self.result_file = result_file
        else:
            raise PCmaster_anno_0_Error("The result_file is None.")

    def set_default_0(self, verbosity=None, dpi=None, facecolor=None):
        # default settings of scanpy
        if verbosity is None:
            verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
        if dpi is None:
            dpi = 80
        if facecolor is None:
            facecolor = 'white'
        sc.settings.verbosity = verbosity  # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=dpi, facecolor=facecolor)

    def data_input_0(self, filepath=None, species=None, organ=None, organs=None):
        if filepath is not None:
            if filepath[-6:] == "tsv.gz":    
                self.obj = sc.read_umi_tools(filepath)
            elif filepath[-6:] == "mtx.gz":
                self.obj = sc.read_10x_mtx(filepath[:-13], var_names='gene_symbols', cache=True)
            elif 'filtered_feature_bc_matrix' in filepath:    
                self.obj = sc.read_10x_mtx(filepath, var_names='gene_symbols', cache=True)
            else:
                self.obj = sc.read_10x_mtx(filepath, var_names='gene_symbols', cache=True)
                
        if organs is not None:
            self.batch_organs = organs
        if species is not None or organ is not None:
            self.species = species
            self.organ = organ
        if self.obj is not None:
            self.genes = self.obj.var_names
            self.cells = self.obj.obs_names
        else:
            raise PCmaster_anno_0_Error("The expression matrix is not successfully loaded.")

    def batch_data_input_0(self, filepaths=None, species=None, organs=None):
        if filepaths is None:
            raise PCmaster_anno_0_Error("The filepaths are None.")
        for i in filepaths:
            if i[-6:] == "tsv.gz":
                self.objs.append(sc.read_umi_tools(i))
            if i[-6:] == "mtx.gz":
                self.objs.append(sc.read_10x_mtx(i[:-13], var_names='gene_symbols', cache=True))
            tmpstr1 = 'filtered_feature_bc_matrix'
            if tmpstr1 in i:
                self.objs.append(sc.read_10x_mtx(i, var_names='gene_symbols', cache=True))
        if species is not None:
            for i in species:
                self.batch_species.append(i)
        if organs is not None:
            for i in organs:
                self.batch_organs.append(i)
        if self.objs is None:
            raise PCmaster_anno_0_Error("The expression matrices are not successfully loaded.")

    def data_input_0_from_anndata(self, anndata=None, species=None, organ=None, organs=None):
        if anndata is not None:
            self.obj = anndata
        else:
            raise PCmaster_anno_0_Error("The anndata is None.")
        if organs is not None:
            self.batch_organs = organs
        if species is not None or organ is not None:
            self.species = species
            self.organ = organ
        if self.obj is not None:
            self.genes = self.obj.var_names
            self.cells = self.obj.obs_names
        else:
            raise PCmaster_anno_0_Error("The expression matrix is not successfully loaded.")
            
    def batch_data_input_0_from_anndata(self, anndatas=None, species=None, organs=None):
        if anndatas is not None:
            for i in anndatas:
                self.objs.append(i)
        else:
            raise PCmaster_anno_0_Error("The anndatas are None.")
        if species is not None:
            for i in species:
                self.batch_species.append(i)
        if organs is not None:
            for i in organs:
                self.batch_organs.append(i)
        if self.objs is None:
            raise PCmaster_anno_0_Error("The expression matrices are not successfully loaded.")

    def concat_0(self):
        batch_categories = []
        if self.objs is None:
            raise PCmaster_anno_0_Error("The expression matrices are not successfully loaded.")
        for i in range(len(self.objs)):
            batch_categories.append("data_number_" + str(i))
        if self.obj is not None:
            print("Warning: self.obj is not None and will be changed.")
        self.obj = self.objs[0].concatenate(self.objs[1:], batch_categories=batch_categories)

    def batch_data_input_0_and_concat_0(self, filepaths=None, species=None, organs=None):
        pass

    def get_dense_0(self):
        if self.obj is None:
            raise PCmaster_anno_0_Error("The expression matrix is not successfully loaded.")
        self.dd = pd.DataFrame(self.obj.X.todense(), index=self.obj.obs_names, columns=self.obj.var_names)

    # Pre-processing, filtering correlation, using basic parameters of scanpy.
    # Added the function of filtering by absolute median difference MAD.
    # MAD = the median of absolute values of all values after the median of Xi-X.
    # The median of # 1, 2, and 3 is 2. If you subtract 2 from each item, you will get -1,0,1. After taking the absolute value, sort it. If it is 0,1,1, then MAD = 1.
    # This function is to keep the data whose values are within the range of+-mad _ coefficient * mad.
    # In general, it doesn't matter if you don't need to filter through the absolute median difference MAD.
    def filter_0(self, n_top=None, num_min_genes=None, num_max_genes=None, num_min_cells=None,
                 max_percent_of_specific_genes=None,
                 genelist=None, mad_coeffcient=None):
        # tmp = self.obj.copy()
        if n_top is None:
            n_top = 20
        if num_min_genes is None:
            num_min_genes = 200
        if num_min_cells is None:
            num_min_cells = 3
        # The following parameter is used to filter specific gene sets, such as mitochondrial genes and so on.
        if max_percent_of_specific_genes is None:
            max_percent_of_specific_genes = 0.05
        self.genes = self.obj.var.index
        
        mad_min_genes = 0
        mad_max_genes = 0
        mad_min_counts = 0
        mad_max_counts = 0
        mad_max_pct = 0
        
        print("\n")
        print("Before filter")
        print("\n")
        sc.pl.highest_expr_genes(self.obj, n_top, )
        sc.pp.filter_cells(self.obj, min_genes=-1)
        sc.pp.filter_genes(self.obj, min_cells=-1)
        if genelist is None:
            self.obj.var['sp'] = self.obj.var_names.str.startswith('meaninglessmeaningless-')
        else:
            genelist = genelist.copy()
            tmp = []
            for i in self.genes:
                if i in genelist:
                    tmp.append(True)
                else:
                    tmp.append(False)
            # Here, add a column of my own identification of genes in the var of scanpy object. For example, if it is a mitochondrial gene, it will be true, and if it is not, it will be false.
            self.obj.var['sp'] = np.asarray(tmp)
        sc.pp.calculate_qc_metrics(self.obj, qc_vars=['sp'], percent_top=None, log1p=False, inplace=True)
        sc.pl.violin(self.obj, ['n_genes_by_counts'])
        sc.pl.violin(self.obj, ['total_counts'])
        sc.pl.violin(self.obj, ['pct_counts_sp'])
        sc.pl.scatter(self.obj, x='total_counts', y='pct_counts_sp')
        sc.pl.scatter(self.obj, x='total_counts', y='n_genes_by_counts')
        if mad_coeffcient is not None:
            mad_min_genes = self.obj.obs.n_genes_by_counts.median() \
                            - mad_coeffcient * (
                                    self.obj.obs.n_genes_by_counts - self.obj.obs.n_genes_by_counts.median()).abs().median()
            mad_max_genes = self.obj.obs.n_genes_by_counts.median() \
                            + mad_coeffcient * (
                                    self.obj.obs.n_genes_by_counts - self.obj.obs.n_genes_by_counts.median()).abs().median()
            mad_min_counts = self.obj.obs.total_counts.median() \
                             - mad_coeffcient * (
                                     self.obj.obs.total_counts - self.obj.obs.total_counts.median()).abs().median()
            mad_max_counts = self.obj.obs.total_counts.median() \
                             + mad_coeffcient * (
                                     self.obj.obs.total_counts - self.obj.obs.total_counts.median()).abs().median()
            # mad_max_pct = self.obj.obs.pct_counts_sp.median() \
            #     + mad_coeffcient*(self.obj.obs.pct_counts_sp - self.obj.obs.pct_counts_sp.median()).abs().median()
        total_counts_describe_before = self.obj.obs.total_counts.describe()
        n_genes_by_counts_describe_before = self.obj.obs.n_genes_by_counts.describe()
        pct_counts_sp_describe_before = self.obj.obs.pct_counts_sp.describe()
        # self.obj = tmp
        # del tmp
        print("\n")
        print("After filter")
        print("\n")
        
        # filter
        if num_min_genes is not None:
            sc.pp.filter_cells(self.obj, min_genes=num_min_genes)
        if num_min_cells is not None:
            sc.pp.filter_genes(self.obj, min_cells=num_min_cells)
        if num_max_genes is not None:
            self.obj = self.obj[self.obj.obs.n_genes_by_counts < num_max_genes, :]
        if max_percent_of_specific_genes is not None:
            self.obj = self.obj[self.obj.obs.pct_counts_sp < max_percent_of_specific_genes, :]
        if mad_coeffcient is not None:
            self.obj = self.obj[self.obj.obs.n_genes_by_counts < mad_max_genes, :]
            self.obj = self.obj[self.obj.obs.n_genes_by_counts > mad_min_genes, :]
            self.obj = self.obj[self.obj.obs.total_counts < mad_max_counts, :]
            self.obj = self.obj[self.obj.obs.total_counts > mad_min_counts, :]
            # if genelist !=None:
            #     self.obj = self.obj[self.obj.obs.pct_counts_sp < mad_max_pct, :]
        
        # show
        sc.pl.violin(self.obj, ['n_genes_by_counts'])
        sc.pl.violin(self.obj, ['total_counts'])
        sc.pl.violin(self.obj, ['pct_counts_sp'])
        sc.pl.scatter(self.obj, x='total_counts', y='pct_counts_sp')
        sc.pl.scatter(self.obj, x='total_counts', y='n_genes_by_counts')
        total_counts_describe_after = self.obj.obs.total_counts.describe()
        n_genes_by_counts_describe_after = self.obj.obs.n_genes_by_counts.describe()
        pct_counts_sp_describe_after = self.obj.obs.pct_counts_sp.describe()
        print("\n")
        print("Before filter")
        print("\n")
        print(total_counts_describe_before)
        print(n_genes_by_counts_describe_before)
        print(pct_counts_sp_describe_before)
        print("\n")
        print("After filter")
        print("\n")
        print(total_counts_describe_after)
        print(n_genes_by_counts_describe_after)
        print(pct_counts_sp_describe_after)

    # Normalization, screening hypervariable genes, etc., are roughly completed with scanpy's own functions.
    def norm_log_hvg_filter_regress_scale_0(self, norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,
                               hvg_min_dispersions=None, n_top_hvg_genes=None,
                               regresslist=None, scale_max_value=None):
        if norm_target_sum is None:
            norm_target_sum = 1e4
        if hvg_min_mean is None:
            hvg_min_mean = 0.0125
        if hvg_max_mean is None:
            hvg_max_mean = 3
        if hvg_min_dispersions is None:
            hvg_min_dispersions = 0.5
        if regresslist is None:
            regresslist = ['total_counts']
        if scale_max_value is None:
            scale_max_value = 10

        if norm_target_sum is not None:
            sc.pp.normalize_total(self.obj, target_sum=norm_target_sum)
            sc.pp.log1p(self.obj)
        if hvg_min_mean is not None and hvg_max_mean is not None and hvg_min_dispersions is not None:
            if n_top_hvg_genes is not None:
                sc.pp.highly_variable_genes(self.obj, min_mean=hvg_min_mean, max_mean=hvg_max_mean,
                                   min_disp=hvg_min_dispersions,n_top_genes=n_top_hvg_genes,)
            else:
                sc.pp.highly_variable_genes(self.obj, min_mean=hvg_min_mean, max_mean=hvg_max_mean,
                                   min_disp=hvg_min_dispersions,)
            sc.pl.highly_variable_genes(self.obj)
            self.obj.raw = self.obj
            self.obj = self.obj[:, self.obj.var.highly_variable]
        if regresslist is not None:
            regresslist = regresslist.copy()
            sc.pp.regress_out(self.obj, regresslist)
        if scale_max_value is not None:
            sc.pp.scale(self.obj, max_value=scale_max_value)

    # Carry out PCA to screen the differentially expressed gene DEGs.
    def pca_cluster_markerSelect_0(self, pca_svd_algorithm=None, # 'randomized','auto','full' # is ok too
                            pca_genelist=None,
                            default_cluster_method=None,leiden_resolution=None,louvain_resolution=None,
                            cluster_n_neighbors=None, cluster_n_pcs=None, use_rep=None, perform_umap=None,
                            umap_genelist=None, perform_paga=None,
                            markerSelect_method=None,  # 'wilcoxon','logreg' # is ok too
                            markerSelect_n_genes=None, show_n=None):
        if pca_svd_algorithm is None:
            pca_svd_algorithm = 'arpack'
        if default_cluster_method is None:
            default_cluster_method = 'leiden'
        if leiden_resolution is None:
            leiden_resolution = 1.0
        if louvain_resolution is None:
            louvain_resolution = 1.0
        if cluster_n_neighbors is None:
            cluster_n_neighbors = 10
        if cluster_n_pcs is None:
            cluster_n_pcs = 40
        if perform_umap is None:
            perform_umap = True
        if perform_paga is None:
            perform_paga = False
        if markerSelect_method is None:
            markerSelect_method = 't-test'  # 'wilcoxon','logreg' # is ok too
        if markerSelect_n_genes is None:
            markerSelect_n_genes = 25
        if show_n is None:
            show_n = 25

        if pca_svd_algorithm is not None:
            sc.tl.pca(self.obj, svd_solver=pca_svd_algorithm)
            if pca_genelist is not None:
                pca_genelist = pca_genelist.copy()
                sc.pl.pca(self.obj, color=pca_genelist)
            sc.pl.pca_variance_ratio(self.obj, log=True)
        if cluster_n_neighbors is not None and cluster_n_pcs is not None and use_rep is None:
            sc.pp.neighbors(self.obj, n_neighbors=cluster_n_neighbors, n_pcs=cluster_n_pcs)
        
        # harmony for batch effect correction
        if cluster_n_neighbors is not None and cluster_n_pcs is not None and use_rep == "X_pca_harmony":
            sce.pp.harmony_integrate(self.obj, 'batch')
            sc.pp.neighbors(self.obj, n_neighbors=cluster_n_neighbors, n_pcs=cluster_n_pcs, use_rep="X_pca_harmony")
        if perform_umap:
            sc.tl.umap(self.obj)
            sc.tl.leiden(self.obj,resolution = leiden_resolution)
            if default_cluster_method == 'louvain':
                sc.tl.louvain(self.obj,resolution = louvain_resolution)
            if default_cluster_method == 'louvain':
                tmplist = ['louvain']
            else:
                tmplist = ['leiden']
            if umap_genelist is not None:
                umap_genelist = umap_genelist.copy()
                for i in umap_genelist:
                    tmplist.append(i)
            sc.pl.umap(self.obj, color=tmplist)
        
        # paga
        if perform_paga:
            sc.tl.paga(self.obj)
            sc.pl.paga(self.obj)  # remove `plot=False` if you want to see the coarse-grained graph
            sc.tl.umap(self.obj, init_pos='paga')
            if default_cluster_method == 'louvain':
                tmplist = ['louvain']
            else:
                tmplist = ['leiden']
            if use_rep == "X_pca_harmony" and default_cluster_method == 'louvain':
                tmplist = ['louvain', 'batch']
            if use_rep == "X_pca_harmony" and default_cluster_method == 'leiden':
                tmplist = ['leiden', 'batch']
            if umap_genelist is not None:
                umap_genelist = umap_genelist.copy()
                for i in umap_genelist:
                    tmplist.append(i)
            for i in tmplist:
                sc.pl.umap(self.obj, color=i)
        if markerSelect_method is not None:
            if default_cluster_method == 'leiden':
                sc.tl.rank_genes_groups(self.obj, 'leiden', method=markerSelect_method)
            else:
                sc.tl.rank_genes_groups(self.obj, 'louvain', method=markerSelect_method)
            sc.pl.rank_genes_groups(self.obj, n_genes=markerSelect_n_genes, sharey=False)
        print("\n")
        print("The below is the top marker gene list for each cluster.\n\
            You can change the annotation result according to this.")
        print("\n")
        print("The list of top marker genes \n col : clusters")
        print("\n")
        print(pd.DataFrame(self.obj.uns['rank_genes_groups']['names']).head(show_n))
        # print( pd.DataFrame(
        #     {group + '_' + key: self.obj.uns['rank_genes_groups'][key][group]
        #     for group in self.obj.uns['rank_genes_groups']['names'].dtype.names \
        # for key in ['names', 'pvals']}).head(show_n) )
        # result = self.obj.uns['rank_genes_groups']
        # groups = result['names'].dtype.names
        # print(result)
        # print("\n")
        # print("row : marker_genes ; col : n , names , p , pvals")
        # print("\n")
        # pd.DataFrame(
        #     {group + '_' + key[:1]: result[key][group]
        #     for group in groups for key in ['names', 'pvals']}).head(show_n)
        result = self.obj.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        
        # save DEGs and p-values
        self.markergenepd_head5 = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(5)
        self.markergenepd_head10 = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(10)
        self.markergenepd_head15 = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(15)
        self.markergenepd_head20 = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(20)
        print("\n")
        print("Top marker genes and their p-values of each cluster: ")
        print("\n")
        if show_n is not None:
            self.markergenepd_head = pd.DataFrame(
                {group + '_' + key[:1]: result[key][group]
                 for group in groups for key in ['names', 'pvals']}).head(show_n)
            print(self.markergenepd_head)
        else:
            print(self.markergenepd_head20)
        # print(self.markergenepd_head5)
        # print(self.markergenepd_head10)
        # print(self.markergenepd_head15)
        # print(self.markergenepd_head20)

    def markergenepd_make_0(self):
        if self.obj is None:
            raise PCmaster_anno_0_Error("The expression matrix is not successfully loaded.")
        result = self.obj.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        self.markergenepd_head5 = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(5)
        self.markergenepd_head10 = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(10)
        self.markergenepd_head15 = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(15)
        self.markergenepd_head20 = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(20)
        # print(self.markergenepd_head5)
        # print(self.markergenepd_head10)
        # print(self.markergenepd_head15)
        # print(self.markergenepd_head20)
        
    def print_organ_types(self):
        ref_txt = open('./plant_marker_gene_list.txt', 'r', encoding='utf-8')
        tmpset_for_print = set()
        for line in ref_txt:
                # if len(line.split('\t')) < 5:
                #     print(line)
                if len(line.split('\t')) >= 5:
                    species = line.split('\t')[0]
                    organ = line.split('\t')[1]
                    celltype = line.split('\t')[2]
                    gene = line.split('\t')[3]
                    if organ not in tmpset_for_print:
                        tmpset_for_print.add(organ)
        tmpset_for_print = sorted(tmpset_for_print)
        for i in tmpset_for_print:
            print(i)
        
    # for generate dict{gene-others} in ./plant_marker_gene_list.txt
    def make_marker_gene_list_0(self,marker_select_flag = None,marker_gene_list = None):
        # organ = None
        # batch_organs = []
        if (marker_select_flag is None) or (marker_select_flag != '#1'):
            marker_select_flag = '#1 and #2'
        ref_txt = None
        if marker_gene_list is None:
            ref_txt = open('./plant_marker_gene_list.txt', 'r', encoding='utf-8')
        else:
            ref_txt = open(marker_gene_list, 'r', encoding='utf-8')
        
        # Here, a filter is carried out first. If it is a single organ type, only the marker gene of that organ type will be left, and others will be deleted.
        ref1 = []
        for line in ref_txt:
            if len(line.split('\t')) >= 5:
                if line.split('\t')[4] == 'Marker #1\n' or line.split('\t')[4] == 'Cell marker #1\n' \
                or line.split('\t')[4] == 'Marker #2\n' or line.split('\t')[4] == 'Cell marker #2\n':
                    the_organ = line.split('\t')[1]
                    if self.organ is not None:
                        if self.organ.upper() in the_organ.upper() or self.organ.upper() == 'ALL':
                            ref1.append(line)
                    elif len(self.batch_organs) != 0:
                        for i in self.batch_organs:
                            if i.upper() in the_organ.upper() or i.upper() == 'ALL':
                                ref1.append(line)
                    else:
                        ref1.append(line)
        
        # The following are the steps to make two dictionary for comparing DEG with known marker.
        dir_gene_species = {}
        dir_gene_organ_celltype = {}
        count = 0
        if marker_select_flag == '#1 and #2':
            forbidden_set = set()
            set_rank1 = set()
            for line in ref1:
                # print(line.split('\t'))
                if len(line.split('\t')) >= 5:
                    if line.split('\t')[4] == 'Marker #1\n' or line.split('\t')[4] == 'Cell marker #1\n':
                        species = line.split('\t')[0]
                        organ = line.split('\t')[1]
                        celltype = line.split('\t')[2]
                        gene = line.split('\t')[3]
                        # print(species)
                        if gene not in forbidden_set:
                            if gene in dir_gene_organ_celltype:
                                # The recognition of'/'indicates that the gene points to at least three different cell types.
                                if '/' in dir_gene_organ_celltype[gene]:
                                    # Delete all the information previously left by this gene.
                                    del dir_gene_organ_celltype[gene]
                                    del dir_gene_species[gene]
                                    # Blacken the gene
                                    forbidden_set.add(gene)
                                else:
                                    tmpinfo_1 = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                    # If the same item is added, it will be skipped and no changes will be made.
                                    if dir_gene_organ_celltype[gene] == tmpinfo_1:
                                        pass
                                    # Only, the same gene is allowed to point to at most two cell types.
                                    else:
                                        dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + \
                                        ' / ' + 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                        dir_gene_species[gene] = species
                            else:
                                # first add
                                dir_gene_organ_celltype[gene] = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                dir_gene_species[gene] = species
                                # record marker#1 gene
                                set_rank1.add(gene)
                        # print(organ)
                        # print(celltype)
                        # print(gene)
                count = count + 1
            for line in ref1:
                # print(line.split('\t'))
                if len(line.split('\t')) < 5:
                    print(line)
                if len(line.split('\t')) >= 5 and ( line.split('\t')[4] == 'Marker #2\n' or \
                                                   line.split('\t')[4] == 'Cell marker #2\n' ):
                    species = line.split('\t')[0]
                    organ = line.split('\t')[1]
                    celltype = line.split('\t')[2]
                    gene = line.split('\t')[3]
                    # organ_yesno = 1
                    # if self.organ is not None:
                    #     if self.organ.upper() not in organ.upper() and self.organ.upper() != 'ALL':
                    #         organ_yesno = 0
                    # elif self.batch_organs is not None:
                    #     tmp_organs = []
                    #     for tmp_organ in self.batch_organs:
                    #         tmp_organs.append(tmp_organ.upper())
                    #     if 'ALL' in tmp_organs:
                    #         organ_yesno = 1
                    #     elif organ.upper() in tmp_organs:
                    #         organ_yesno = 1
                    #     else:
                    #         organ_yesno = 0
                    # if organ_yesno == 0:
                    #     continue
                    # print(species)
                    
                    # marker#1 and marker#2 will not interact each other
                    if (gene not in forbidden_set) and (gene not in set_rank1):
                        if gene in dir_gene_organ_celltype:
                            if '/' in dir_gene_organ_celltype[gene]:
                                del dir_gene_organ_celltype[gene]
                                del dir_gene_species[gene]
                                forbidden_set.add(gene)
                            else:
                                tmpinfo_1 = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                if dir_gene_organ_celltype[gene] == tmpinfo_1:
                                    pass
                                else:
                                    dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + \
                                    ' / ' + 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                    dir_gene_species[gene] = species
                        else:
                            dir_gene_organ_celltype[gene] = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                            dir_gene_species[gene] = species
                    # print(organ)
                    # print(celltype)
                    # print(gene)
                count = count + 1
            self.marker_grade_1 = set_rank1
        else:
            # the below only process for marker#1
            forbidden_set = set()
            for line in ref1:
                # print(line.split('\t'))
                if len(line.split('\t')) < 5:
                    print(line)
                if len(line.split('\t')) >= 5:
                    if line.split('\t')[4] == 'Marker #1\n' or line.split('\t')[4] == 'Cell marker #1\n':
                        species = line.split('\t')[0]
                        organ = line.split('\t')[1]
                        celltype = line.split('\t')[2]
                        gene = line.split('\t')[3]
                        # print(species)
                        if gene not in forbidden_set:
                            if gene in dir_gene_organ_celltype:
                                if '/' in dir_gene_organ_celltype[gene]:
                                    del dir_gene_organ_celltype[gene]
                                    del dir_gene_species[gene]
                                    forbidden_set.add(gene)
                                else:
                                    tmpinfo_1 = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                    if dir_gene_organ_celltype[gene] == tmpinfo_1:
                                        pass
                                    else:
                                        dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + \
                                        ' / ' + 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                        dir_gene_species[gene] = species
                            else:
                                dir_gene_organ_celltype[gene] = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                dir_gene_species[gene] = species
                        # print(organ)
                        # print(celltype)
                        # print(gene)
                count = count + 1
        return dir_gene_species, dir_gene_organ_celltype

    def annotation_with_marker_gene_list_0(self,marker_gene_list = None,marker_select_flag = None,
                                           marker_grade_1_coefficient = None):
        # Update the unchanged default situation, which is the weight of marker#1. Only when marker#1 and marker#2 are mixed,
        # and one gene is marker#1 gene and only corresponds to one cell type.
        if marker_grade_1_coefficient is None:
            marker_grade_1_coefficient = 1
        
        # get dictionary
        self.dir_gene_species, self.dir_gene_organ_celltype = self.make_marker_gene_list_0(marker_select_flag = '#1 and #2',
                                                         marker_gene_list = marker_gene_list)
        dir_gene_species, dir_gene_organ_celltype = self.make_marker_gene_list_0(marker_select_flag = '#1 and #2',
                                                         marker_gene_list = marker_gene_list)
        count_1 = 0
        count_2 = 0
        anno_re = []

        # use self.markergenepd_head20 as a temporary carrier, it will be changed back later.
        tmpmarkergenepd = self.markergenepd_head20.copy()
        if self.markergenepd_head is not None:
            self.markergenepd_head20 = self.markergenepd_head
        for i in range(0, int(self.markergenepd_head20.shape[1] / 2)):
            tmp = {}
            # According to the original organization form of DEG table, DEG is extracted independently for each cluster, 
            # and then compared with the known marker table one by one.
            for j in self.markergenepd_head20[str(i) + '_n'].values:
                j = str(j).strip()
                # print(j)
                # break
                if j in dir_gene_species:
                    organ_celltype = dir_gene_organ_celltype[j]
                    organ_celltype_1 = None
                    organ_celltype_2 = None
                    if '/' in organ_celltype:
                        organ_celltype_1 = organ_celltype.split('/')[0].strip()
                        organ_celltype_2 = organ_celltype.split('/')[1].strip()
                        
                    if organ_celltype in tmp:
                        if j in self.marker_grade_1 and ('/' in organ_celltype) is False:
                            tmp[organ_celltype] = tmp[organ_celltype] + 1*marker_grade_1_coefficient
                        else:
                            tmp[organ_celltype] = tmp[organ_celltype] + 1
                        if organ_celltype_1 is not None and organ_celltype_1 in tmp:
                            tmp[organ_celltype_1] = tmp[organ_celltype_1] + 1
                        if organ_celltype_2 is not None and organ_celltype_2 in tmp:
                            tmp[organ_celltype_2] = tmp[organ_celltype_2] + 1
                    else:
                        if j in self.marker_grade_1 and ('/' in organ_celltype) is False:
                            tmp[organ_celltype] = 1*marker_grade_1_coefficient
                        else:
                            tmp[organ_celltype] = 1
                        if organ_celltype_1 is not None and (organ_celltype_1 in tmp) is False:
                            tmp[organ_celltype_1] = 1
                        if organ_celltype_2 is not None and (organ_celltype_2 in tmp) is False:
                            tmp[organ_celltype_2] = 1
            if tmp != {}:
                anno_re.append(max(tmp, key=tmp.get))
            else:
                anno_re.append('Unknown')
        self.anno_re_using_rank_1and2_markers = anno_re

        # for i in anno_re:
        #     print(i)
        # print(len(anno_re))
        anno_re_1and2 = anno_re

        dir_gene_species, dir_gene_organ_celltype = self.make_marker_gene_list_0(marker_select_flag = '#1',
                                                                                 marker_gene_list = marker_gene_list)

        count_1 = 0
        count_2 = 0
        anno_re = []
        for i in range(0, int(self.markergenepd_head20.shape[1] / 2)):
            tmp = {}
            for j in self.markergenepd_head20[str(i) + '_n'].values:
                j = str(j).strip()
                # print(j)
                # break
                if j in dir_gene_species:
                    organ_celltype = dir_gene_organ_celltype[j]
                    organ_celltype_1 = None
                    organ_celltype_2 = None
                    if '/' in organ_celltype:
                        organ_celltype_1 = organ_celltype.split('/')[0].strip()
                        organ_celltype_2 = organ_celltype.split('/')[1].strip()
                        
                    if organ_celltype in tmp:
                        tmp[organ_celltype] = tmp[organ_celltype] + 1
                        if organ_celltype_1 is not None and organ_celltype_1 in tmp:
                            tmp[organ_celltype_1] = tmp[organ_celltype_1] + 1
                        if organ_celltype_2 is not None and organ_celltype_2 in tmp:
                            tmp[organ_celltype_2] = tmp[organ_celltype_2] + 1
                    else:
                        tmp[organ_celltype] = 1
                        if organ_celltype_1 is not None and (organ_celltype_1 in tmp) is False:
                            tmp[organ_celltype_1] = 1
                        if organ_celltype_2 is not None and (organ_celltype_2 in tmp) is False:
                            tmp[organ_celltype_2] = 1
            if tmp != {}:
                anno_re.append(max(tmp, key=tmp.get))
            else:
                anno_re.append('Unknown')
        self.anno_re_using_rank_1_markers = anno_re

        # for i in anno_re:
        #     print(i)
        # print(len(anno_re))
        anno_re_1 = anno_re

        # print("The below is the top marker gene list for each cluster.\n\
        #     You can change the annotation result according to this.")
        # print(pd.DataFrame(self.obj.uns['rank_genes_groups']['names']).head(20))

        tmplist = []
        for i in self.obj.obs['leiden']:
            tmplist.append(self.anno_re_using_rank_1and2_markers[int(i)])
        # print(tmplist)
        # print(len(self.anno_re_using_rank_1and2_markers))
        # add annotation results
        self.obj.obs['anno_re_using_rank_1and2_markers'] = tmplist.copy()

        for i, j in enumerate(anno_re_1and2):
            print('cluster ' + str(i) + ': ' + str(j))
        sc.pl.umap(self.obj, color='anno_re_using_rank_1and2_markers')

        tmplist = []
        for i in self.obj.obs['leiden']:
            tmplist.append(self.anno_re_using_rank_1_markers[int(i)])
        # print(tmplist)
        # print(len(self.anno_re_using_rank_1_markers))
        # add annotation results
        self.obj.obs['anno_re_using_rank_1_markers'] = tmplist.copy()
        for i, j in enumerate(anno_re_1):
            print('cluster ' + str(i) + ': ' + str(j))
        sc.pl.umap(self.obj, color='anno_re_using_rank_1_markers')

        self.markergenepd_head20 = tmpmarkergenepd.copy()
        del tmpmarkergenepd
        self.memory_collect_0()
        
    def get_dir_gene_score(self,marker_gene_list = None):
        self.dir_gene_score = {}
        if marker_gene_list is None:
            ref_txt = open('./plant_marker_gene_list.txt', 'r', encoding='utf-8')
        else:
            ref_txt = open(marker_gene_list, 'r', encoding='utf-8')
            
        ref1 = []
        for line in ref_txt:
            if len(line.split('\t')) == 6:
                ref1.append(line.strip())
        for line in ref1:
            gene = line.split('\t')[3].strip()
            score = float(line.split('\t')[5].strip())
            self.dir_gene_score[gene] = score
        
    def annotation_with_marker_gene_scores(self,the_least_score = None):
        if the_least_score is None:
            the_least_score = 0.0
        count_1 = 0
        count_2 = 0
        anno_re = []

        tmpmarkergenepd = self.markergenepd_head20.copy()
        if self.markergenepd_head is not None:
            self.markergenepd_head20 = self.markergenepd_head
        for i in range(0, int(self.markergenepd_head20.shape[1] / 2)):
            tmp = {}
            for j in self.markergenepd_head20[str(i) + '_n'].values:
                j = str(j).strip()
                # print(j)
                # break
                if j in self.dir_gene_species and j in self.dir_gene_score and float(self.dir_gene_score[j])>=the_least_score:
                    organ_celltype = self.dir_gene_organ_celltype[j]
                    organ_celltype_1 = None
                    organ_celltype_2 = None
                    if '/' in organ_celltype:
                        organ_celltype_1 = organ_celltype.split('/')[0].strip()
                        organ_celltype_2 = organ_celltype.split('/')[1].strip()
                        
                    if organ_celltype in tmp:
                        tmp[organ_celltype] = tmp[organ_celltype] + 1.0*float(self.dir_gene_score[j])
                        if organ_celltype_1 is not None and organ_celltype_1 in tmp:
                            tmp[organ_celltype_1] = tmp[organ_celltype_1] + 1.0*float(self.dir_gene_score[j])
                        if organ_celltype_2 is not None and organ_celltype_2 in tmp:
                            tmp[organ_celltype_2] = tmp[organ_celltype_2] + 1.0*float(self.dir_gene_score[j])
                    else:
                        tmp[organ_celltype] = 1.0*float(self.dir_gene_score[j])
                        if organ_celltype_1 is not None and (organ_celltype_1 in tmp) is False:
                            tmp[organ_celltype_1] = 1.0*float(self.dir_gene_score[j])
                        if organ_celltype_2 is not None and (organ_celltype_2 in tmp) is False:
                            tmp[organ_celltype_2] = 1.0*float(self.dir_gene_score[j])
            if tmp != {}:
                anno_re.append(max(tmp, key=tmp.get))
            else:
                anno_re.append('Unknown')
            print(i)
            for tmp_key in tmp:
                print(tmp_key)
                print(tmp[tmp_key])
        self.anno_re_using_marker_scores = anno_re
        
        tmplist = []
        for i in self.obj.obs['leiden']:
            tmplist.append(self.anno_re_using_marker_scores[int(i)])
        # print(tmplist)
        # print(len(self.anno_re_using_rank_1and2_markers))
        self.obj.obs['anno_re_using_marker_scores'] = tmplist.copy()

        for i, j in enumerate(anno_re):
            print('cluster ' + str(i) + ': ' + str(j))
        sc.pl.umap(self.obj, color='anno_re_using_marker_scores')

        self.markergenepd_head20 = tmpmarkergenepd.copy()
        del tmpmarkergenepd
        self.memory_collect_0()

    # integrated function
    def to_clusters_0(self, verbosity=None, dpi=None, facecolor=None, # set_default_0
            anndata=None,filepath=None, species=None, organ=None, organs=None, # data_input_0
            n_top=None, num_min_genes=None, num_max_genes=None, num_min_cells=None, # filter_0
                max_percent_of_specific_genes=None,
                genelist=None, mad_coeffcient=None,
            norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,  # norm_log_hvg_filter_regress_scale_0
                hvg_min_dispersions=None, regresslist=None, scale_max_value=None,n_top_hvg_genes=None,
            pca_svd_algorithm=None, pca_genelist=None, # pca_cluster_markerSelect_0
                default_cluster_method=None,leiden_resolution=None,louvain_resolution=None,
                cluster_n_neighbors=None, cluster_n_pcs=None, use_rep=None, perform_umap=None,
                umap_genelist=None, perform_paga=None,
                markerSelect_method=None,  # 'wilcoxon','logreg' # is ok too
                markerSelect_n_genes=None, show_n=None):
        self.set_default_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor)
        if filepath is not None:
            self.data_input_0(filepath=filepath, species=species, organ=organ, organs=organs)
        if anndata is not None:
            self.data_input_0_from_anndata(anndata=anndata, species=species, organ=organ, organs=organs)
        self.filter_0(n_top=n_top, num_min_genes=num_min_genes, 
            num_max_genes=num_max_genes, num_min_cells=num_min_cells,
            max_percent_of_specific_genes=max_percent_of_specific_genes,
            genelist=genelist, mad_coeffcient=mad_coeffcient)
        self.norm_log_hvg_filter_regress_scale_0(norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, 
            hvg_max_mean=hvg_max_mean, hvg_min_dispersions=hvg_min_dispersions, n_top_hvg_genes=n_top_hvg_genes,
            regresslist=regresslist, scale_max_value=scale_max_value)
        self.pca_cluster_markerSelect_0(pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist,
            default_cluster_method=default_cluster_method,
            leiden_resolution=leiden_resolution,louvain_resolution=louvain_resolution,
            cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, use_rep=use_rep, 
            perform_umap=perform_umap,umap_genelist=umap_genelist, perform_paga=perform_paga,
            markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
            markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)

    # integrated function
    def batch_to_clusters_0(self, verbosity=None, dpi=None, facecolor=None, # set_default_0
            anndatas=None,filepaths=None, species=None, organs=None, # batch_data_input_0
            n_top=None, num_min_genes=None, num_max_genes=None, num_min_cells=None, # filter_0
                max_percent_of_specific_genes=None,
                genelist=None, mad_coeffcient=None,
            norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,  # norm_log_hvg_filter_regress_scale_0
                hvg_min_dispersions=None, regresslist=None, scale_max_value=None,n_top_hvg_genes=None,
            pca_svd_algorithm=None, pca_genelist=None, # pca_cluster_markerSelect_0
                default_cluster_method=None,leiden_resolution=None,louvain_resolution=None,
                cluster_n_neighbors=None, cluster_n_pcs=None, use_rep=None, perform_umap=None,
                umap_genelist=None, perform_paga=None,
                markerSelect_method=None,  # 'wilcoxon','logreg' # is ok too
                markerSelect_n_genes=None, show_n=None):
        self.set_default_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor)
        if filepaths is not None:
            self.batch_data_input_0(filepaths=filepaths, species=species, organs=organs)
        if anndatas is not None:
            self.batch_data_input_0_from_anndata(anndatas=anndatas, species=species, organs=organs)
        self.concat_0()
        self.filter_0(n_top=n_top, num_min_genes=num_min_genes, 
            num_max_genes=num_max_genes, num_min_cells=num_min_cells,
            max_percent_of_specific_genes=max_percent_of_specific_genes,
            genelist=genelist, mad_coeffcient=mad_coeffcient)
        self.norm_log_hvg_filter_regress_scale_0(norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, 
            hvg_max_mean=hvg_max_mean, hvg_min_dispersions=hvg_min_dispersions, n_top_hvg_genes=n_top_hvg_genes,
            regresslist=regresslist, scale_max_value=scale_max_value)
        self.pca_cluster_markerSelect_0(pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist,
            default_cluster_method=default_cluster_method,
            leiden_resolution=leiden_resolution,louvain_resolution=louvain_resolution,
            cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, use_rep=use_rep, 
            perform_umap=perform_umap,umap_genelist=umap_genelist, perform_paga=perform_paga,
            markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
            markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)

    # integrated function
    def to_annotation_0(self, verbosity=None, dpi=None, facecolor=None, # set_default_0
            species=None,anndata=None,filepath=None,organ=None,anndatas=None,filepaths=None,organs=None, 
                        # batch_data_input_0
            n_top=None, num_min_genes=None, num_max_genes=None, num_min_cells=None, # filter_0
                max_percent_of_specific_genes=None,
                genelist=None, mad_coeffcient=None,
            norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,  # norm_log_hvg_filter_regress_scale_0
                hvg_min_dispersions=None, regresslist=None, scale_max_value=None,n_top_hvg_genes=None,
            pca_svd_algorithm=None, pca_genelist=None, # pca_cluster_markerSelect_0
                default_cluster_method=None,leiden_resolution=None,louvain_resolution=None,
                cluster_n_neighbors=None, cluster_n_pcs=None, use_rep=None, perform_umap=None,
                umap_genelist=None, perform_paga=None,
                markerSelect_method=None,  # 'wilcoxon','logreg' # is ok too
                markerSelect_n_genes=None, show_n=None,
            marker_select_flag = None,marker_grade_1_coefficient = None,marker_gene_list = None,
            the_least_score = None):
        if filepath is not None:
            self.to_clusters_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor, # set_default_0
                anndata=anndata,filepath=filepath, species=species, organ=organ, organs=organs, # data_input_0
                n_top=n_top, num_min_genes=num_min_genes, num_max_genes=num_max_genes, # filter_0
                    num_min_cells=num_min_cells, 
                    max_percent_of_specific_genes=max_percent_of_specific_genes,
                    genelist=genelist, mad_coeffcient=mad_coeffcient,
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  
                               # norm_log_hvg_filter_regress_scale_0
                    hvg_min_dispersions=hvg_min_dispersions, regresslist=regresslist, scale_max_value=scale_max_value,
                               n_top_hvg_genes=n_top_hvg_genes,
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
                    default_cluster_method=default_cluster_method,
                    leiden_resolution=leiden_resolution,louvain_resolution=louvain_resolution,
                    cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, 
                    use_rep=use_rep, perform_umap=perform_umap,
                    umap_genelist=umap_genelist, perform_paga=perform_paga,
                    markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
                    markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)
        elif anndata is not None:
            self.to_clusters_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor, # set_default_0
                anndata=anndata,filepath=filepath, species=species, organ=organ, organs=organs, # batch_data_input_0
                n_top=n_top, num_min_genes=num_min_genes, num_max_genes=num_max_genes, # filter_0
                    num_min_cells=num_min_cells, 
                    max_percent_of_specific_genes=max_percent_of_specific_genes,
                    genelist=genelist, mad_coeffcient=mad_coeffcient,
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  
                               # norm_log_hvg_filter_regress_scale_0
                    hvg_min_dispersions=hvg_min_dispersions, regresslist=regresslist, scale_max_value=scale_max_value,
                               n_top_hvg_genes=n_top_hvg_genes,
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
                    default_cluster_method=default_cluster_method,
                    leiden_resolution=leiden_resolution,louvain_resolution=louvain_resolution,
                    cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, 
                    use_rep=use_rep, perform_umap=perform_umap,
                    umap_genelist=umap_genelist, perform_paga=perform_paga,
                    markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
                    markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)
        elif filepaths is not None:
            if use_rep is None:
                use_rep="X_pca_harmony"
            self.batch_to_clusters_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor, # set_default_0
                anndatas=anndatas,filepaths=filepaths, species=species, organs=organs, # batch_data_input_0
                n_top=n_top, num_min_genes=num_min_genes, num_max_genes=num_max_genes, # filter_0
                    num_min_cells=num_min_cells, 
                    max_percent_of_specific_genes=max_percent_of_specific_genes,
                    genelist=genelist, mad_coeffcient=mad_coeffcient,
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  
                                     # norm_log_hvg_filter_regress_scale_0
                    hvg_min_dispersions=hvg_min_dispersions, regresslist=regresslist, scale_max_value=scale_max_value,
                                     n_top_hvg_genes=n_top_hvg_genes,
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
                    default_cluster_method=default_cluster_method,
                    leiden_resolution=leiden_resolution,louvain_resolution=louvain_resolution,
                    cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, 
                    use_rep=use_rep, perform_umap=perform_umap,
                    umap_genelist=umap_genelist, perform_paga=perform_paga,
                    markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
                    markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)
        elif anndatas is not None:
            if use_rep is None:
                use_rep="X_pca_harmony"
            self.batch_to_clusters_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor, # set_default_0
                anndatas=anndatas,filepaths=filepaths, species=species, organs=organs, # batch_data_input_0
                n_top=n_top, num_min_genes=num_min_genes, num_max_genes=num_max_genes, # filter_0
                    num_min_cells=num_min_cells, 
                    max_percent_of_specific_genes=max_percent_of_specific_genes,
                    genelist=genelist, mad_coeffcient=mad_coeffcient,
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  
                                     # norm_log_hvg_filter_regress_scale_0
                    hvg_min_dispersions=hvg_min_dispersions, regresslist=regresslist, scale_max_value=scale_max_value,
                                     n_top_hvg_genes=n_top_hvg_genes,
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
                    default_cluster_method=default_cluster_method,
                    leiden_resolution=leiden_resolution,louvain_resolution=louvain_resolution,
                    cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, 
                    use_rep=use_rep, perform_umap=perform_umap,
                    umap_genelist=umap_genelist, perform_paga=perform_paga,
                    markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
                    markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)
        elif self.obj is not None:
            self.to_clusters_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor, # set_default_0
                filepath=filepath, species=species, organ=organ, organs=organs, # data_input_0
                n_top=n_top, num_min_genes=num_min_genes, num_max_genes=num_max_genes, # filter_0
                    num_min_cells=num_min_cells, 
                    max_percent_of_specific_genes=max_percent_of_specific_genes,
                    genelist=genelist, mad_coeffcient=mad_coeffcient,
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  
                               # norm_log_hvg_filter_regress_scale_0
                    hvg_min_dispersions=hvg_min_dispersions, regresslist=regresslist, scale_max_value=scale_max_value,
                               n_top_hvg_genes=n_top_hvg_genes,
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
                    default_cluster_method=default_cluster_method,
                    leiden_resolution=leiden_resolution,louvain_resolution=louvain_resolution,
                    cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, 
                    use_rep=use_rep, perform_umap=perform_umap,
                    umap_genelist=umap_genelist, perform_paga=perform_paga,
                    markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
                    markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)
        else:
            raise PCmaster_anno_0_Error("No input.")
        if self.obj is not None:
            self.markergenepd_make_0()
            self.annotation_with_marker_gene_list_0(marker_select_flag = marker_select_flag,
                                                    marker_grade_1_coefficient = marker_grade_1_coefficient,
                                                    marker_gene_list = marker_gene_list)
            self.get_dir_gene_score(marker_gene_list = marker_gene_list)
            self.annotation_with_marker_gene_scores(the_least_score = the_least_score)

    # not used
    def batch_data_input_for_CRA004082_PRJCA004855_0(self,
                                                     filepaths=None, species=None,
                                                     organs=None, labelpath=None):
        if labelpath is not None:
            self.cellnamelabel = pd.read_csv(labelpath)
        for i in filepaths:
            if i[-6:] == "tsv.gz":
                tmp = sc.read_umi_tools(i)
                print(tmp.obs_names)
                cellnamelisttmp = []
                for j in tmp.obs_names:
                    j = j[0:-2] + '-' + i.split('/')[7].replace('-', '').replace('_', '')
                    # print(i)
                    # break
                    cellnamelisttmp.append(j)

                tmp.obs_names = cellnamelisttmp
                print(tmp.obs_names)
                tmp_any_tmp_list = []
                label_set = set()
                tflist = []
                for i in self.cellnamelabel['cellname']:
                    label_set.add(str(i))
                for i in tmp.obs_names:
                    tmp_any_tmp_list.append(str(i))
                for i in tmp_any_tmp_list:
                    if i in label_set:
                        tflist.append(1)
                    else:
                        tflist.append(0)

                ls = {'cellname': tmp_any_tmp_list,
                      'tf': tflist}
                tfpd = pd.DataFrame(ls)
                tfpd.describe()

                tmp.obs['tf'] = np.array(tfpd['tf']).tolist()
                tmp = tmp[tmp.obs.tf == 1, :]
                print(tmp.obs_names)

                # self.objs.append( tmp )
                del tmp
                del cellnamelisttmp
                del tmp_any_tmp_list
                del label_set
                del tflist
                del ls
                del tfpd
                gc.collect()
                info = psutil.virtual_memory()
                print(psutil.Process(os.getpid()).memory_info().rss)
                print(info.total)
                print(info.percent)
                print(psutil.cpu_count())
            if i[-6:] == "mtx.gz":
                tmp = sc.read_10x_mtx(i[:-13], var_names='gene_symbols', cache=False)
                print(tmp.obs_names)
                cellnamelisttmp = []
                for j in tmp.obs_names:
                    j = j[0:-2] + '-' + i.split('/')[7].replace('-', '').replace('_', '')
                    # print(i)
                    # break
                    cellnamelisttmp.append(j)

                tmp.obs_names = cellnamelisttmp
                print(tmp.obs_names)
                tmp_any_tmp_list = []
                label_set = set()
                tflist = []
                for i in self.cellnamelabel['cellname']:
                    label_set.add(str(i))
                for i in tmp.obs_names:
                    tmp_any_tmp_list.append(str(i))
                for i in tmp_any_tmp_list:
                    if i in label_set:
                        tflist.append(1)
                    else:
                        tflist.append(0)

                ls = {'cellname': tmp_any_tmp_list,
                      'tf': tflist}
                tfpd = pd.DataFrame(ls)
                tfpd.describe()

                tmp.obs['tf'] = np.array(tfpd['tf']).tolist()
                tmp = tmp[tmp.obs.tf == 1, :]
                print(tmp.obs_names)

                # self.objs.append( tmp )
                del tmp
                del cellnamelisttmp
                del tmp_any_tmp_list
                del label_set
                del tflist
                del ls
                del tfpd
                gc.collect()
                info = psutil.virtual_memory()
                print(psutil.Process(os.getpid()).memory_info().rss)
                print(info.total)
                print(info.percent)
                print(psutil.cpu_count())
        if species is not None:
            for i in species:
                self.batch_species.append(i)
        if organs is not None:
            for i in organs:
                self.batch_organs.append(i)

    # not used
    def process_for_CRA004082_PRJCA004855_0(self, filepaths=None, use_rep="X_pca_harmony",
                                            labelpath=None):
        # if labelpath is not None:
        #     self.cellnamelabel = pd.read_csv(labelpath)
        self.set_default_0()
        self.batch_data_input_for_CRA004082_PRJCA004855_0(filepaths, labelpath)
        self.concat_0()

        self.objs = []
        gc.collect()
        info = psutil.virtual_memory()
        print(psutil.Process(os.getpid()).memory_info().rss)
        print(info.total)
        print(info.percent)
        print(psutil.cpu_count())

        self.filter_0()
        self.norm_log_hvg_filter_regress_scale_0()
        self.pca_cluster_markerSelect_0(use_rep=use_rep)
        
    # input ref and test, annotation by deep learning model
    def auto_anno_mldl_train_0(self,ref_obj=None,test_obj=None,
                               random_sample_num=1000,
                               gpu_code = None,
                               learning_rate = None,epochs = None,
                               batch_size = None,dropout = None,
                               last_linear_n = None,num_workers = None,
                               model_weight_save_path = None,
                               mapping_1_save_path = None,
                               mapping_2_save_path = None,
                               num_types_save_path = None):
        # which GPU
        if gpu_code is None:
            gpu_code = 0
        self.accuracys = []
        self.precisions = []
        self.recalls = []
        lr = self.set_the_value(input_value = learning_rate,default = 0.01)
        num_epochs = self.set_the_value(input_value = epochs,default = 10)
        the_batch_size = self.set_the_value(input_value = batch_size,default = 16)
        the_num_workers = self.set_the_value(input_value = num_workers,default = 4)
        the_dropout = self.set_the_value(input_value = dropout,default = 0.3)
        the_last_linear_n = self.set_the_value(input_value = last_linear_n,default = 1000)
        print('lr='+str(lr))
        print('num_epochs='+str(num_epochs))
        print('the_batch_size='+str(the_batch_size))
        print('the_num_workers='+str(the_num_workers))
        print('the_dropout='+str(the_dropout))
        print('the_last_linear_n='+str(the_last_linear_n))

        try:   
            df,mapping_1,mapping_2 = self.auto_annotation_with_deep_learning_transfer_adata_to_df_0(adata_input = ref_obj)  
        except Exception:   
            df,mapping_1,mapping_2 = self.auto_annotation_with_deep_learning_transfer_adata_to_df_without_todense_0(adata_input = ref_obj) 
        
        # displays the number of categories and the number of each category.
        category_counts = df['celltype'].value_counts()
        print("the number of categories: ", len(category_counts))
        print("the number of each category: ")
        print(category_counts)
        print(mapping_1)
        print(mapping_2)
        
        category_counts = df['celltype'].value_counts()
        train_data = pd.DataFrame()
        valid_data = pd.DataFrame()
        test_data = pd.DataFrame()
        # for each celltype, catch random_sample_num cells
        for category in category_counts.index:
            category_data = df[df['celltype'] == category]
            print(category)
            print(category_data.shape[0])
            # print(type(category_data.shape[0]))
            category_data = category_data.sample(n=min(random_sample_num,int(category_data.shape[0])), random_state=1, replace=True)
            # from sklearn.model_selection import train_test_split
            train, temp = train_test_split(category_data, train_size=0.8, random_state=1)
            valid, test = train_test_split(temp, train_size=0.5, random_state=1)
            train_data = pd.concat([train_data, train])
            valid_data = pd.concat([valid_data, valid])
            test_data = pd.concat([test_data, test])
        category_counts = train_data['celltype'].value_counts()
        train_data_balanced = pd.DataFrame()
        for category in category_counts.index:
            category_data = train_data[train_data['celltype'] == category]
            print(category)
            print(category_data.shape[0])
            # print(type(category_data.shape[0]))
            category_data = category_data.sample(n=random_sample_num, random_state=1, replace=True)
            train_data_balanced = pd.concat([train_data_balanced, category_data])
        train_data = train_data_balanced
            
        train = train_data
        valid = valid_data
        test = test_data
        
        # divide data and label
        # +0.001
        train_data = train.iloc[:, :-2]
        train_label = train.iloc[:, -1:]
        train_data = train_data + 0.001
        valid_data = valid.iloc[:, :-2]
        valid_label = valid.iloc[:, -1:]
        valid_data = valid_data + 0.001
        test_data = test.iloc[:, :-2]
        test_label = test.iloc[:, -1:]
        test_data = test_data + 0.001
            
        import torch
            
        traindata = trainDataset(train_data,train_label)
        validdata = validDataset(valid_data,valid_label)
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        
        num_types = train_label.nunique().values[0]

        import torchvision.models
        
        class PCMA_MYGO_RC(nn.Module):
            def __init__(self):
                super(PCMA_MYGO_RC, self).__init__()
                # change "pretrained=False" to "weights=None" may be better
                self.resnet18 = torchvision.models.resnet18(weights=None)   
                self.PCMA_MYGO_FC = nn.Sequential(
                                nn.ReLU(),
                                nn.Dropout(the_dropout),
                                nn.Linear(the_last_linear_n, num_types))
                
            def forward(self, x):
                x = self.resnet18(x)
                x = self.PCMA_MYGO_FC(x)
                return x
        
        pcma_mygo_rc = PCMA_MYGO_RC()
        # load pth
        pretrained_dict = torch.load('resnet18.pth')
        model_dict = pcma_mygo_rc.resnet18.state_dict()
        
        pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
        
        model_dict.update(pretrained_dict)
        pcma_mygo_rc.resnet18.load_state_dict(model_dict)
        
        import torch
        torch.cuda.init()
                    
        # GPU
        def evaluate_accuracy_gpu(net, data_iter, device=None): 
            if isinstance(net, nn.Module):
                net.eval()
                if not device:
                    device = next(iter(net.parameters())).device
            metric = d2l.Accumulator(2)
            with torch.no_grad():
                for X, y in data_iter:
                    if isinstance(X, list):
                        X = [x.to(device) for x in X]
                    else:
                        X = X.to(device)
                    y = y.to(device)
                    metric.add(d2l.accuracy(net(X), y), y.numel())
            return metric[0] / metric[1]
        
        # GPU
        def train_with_GPU(net, train_iter, valid_iter, num_epochs, lr, device):
            def init_weights(m):
                if type(m) == nn.Linear or type(m) == nn.Conv2d:
                    nn.init.xavier_uniform_(m.weight)
            net.apply(init_weights)
            print('training on', device)
            net.to(device)
            optimizer = torch.optim.SGD(net.parameters(), lr=lr, momentum=0.9)
            scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=150, eta_min=0.00001)
            loss = nn.CrossEntropyLoss()
            animator = d2l.Animator(xlabel='epoch', xlim=[1, num_epochs],ylim=[0,1],
                                    legend=['train loss', 'train acc', 'valid acc'])
            timer, num_batches = d2l.Timer(), len(train_iter)
            for epoch in range(num_epochs):
                metric = d2l.Accumulator(3)
                net.train()
                for i, (X, y) in enumerate(train_iter):
                    timer.start()
                    optimizer.zero_grad()
                    X, y = X.to(device), y.to(device)
                    y_hat = net(X)
                    l = loss(y_hat, y)
                    l.backward()
                    optimizer.step()
                    scheduler.step()
                    with torch.no_grad():
                        metric.add(l * X.shape[0], d2l.accuracy(y_hat, y), X.shape[0])
                    timer.stop()
                    train_l = metric[0] / metric[2]
                    train_acc = metric[1] / metric[2]
                    if (i + 1) % (num_batches // 5) == 0 or i == num_batches - 1:
                        animator.add(epoch + (i + 1) / num_batches,
                                     (train_l, train_acc, None))
                valid_acc = evaluate_accuracy_gpu(net, valid_iter)
                animator.add(epoch + 1, (None, None, valid_acc))
            print(f'loss {train_l:.3f}, train acc {train_acc:.3f}, '
                  f'valid acc {valid_acc:.3f}')
            print(f'{metric[2] * num_epochs / timer.sum():.1f} examples/sec '
                  f'on {str(device)}')
        train_with_GPU(pcma_mygo_rc, train_iter, valid_iter, num_epochs, lr, d2l.try_all_gpus()[gpu_code])
        self.model_weight = pcma_mygo_rc.state_dict()
        self.mapping_1 = mapping_1
        self.mapping_2 = mapping_2
        self.num_types = num_types
        if model_weight_save_path is not None and \
                mapping_1_save_path is not None and \
                mapping_2_save_path is not None and \
                num_types_save_path is not None:
            torch.save(pcma_mygo_rc.state_dict(), model_weight_save_path)
            tmp_txt = open(mapping_1_save_path, "w", encoding='utf-8')
            for key, value in mapping_1.items():
                tmp_txt.write(str(key) + ": " + str(value) + "\n")
            tmp_txt.close()
            tmp_txt = open(mapping_2_save_path, "w", encoding='utf-8')
            for key, value in mapping_2.items():
                tmp_txt.write(str(key) + ": " + str(value) + "\n")
            tmp_txt.close()
            tmp_txt = open(num_types_save_path, "w", encoding='utf-8')
            tmp_txt.write(str(num_types) + "\n")
            tmp_txt.close()
        else:
            now = datetime.datetime.now()
            current_time = now.strftime("%Y-%m-%d %H:%M:%S")
            print(current_time)
            torch.save(pcma_mygo_rc.state_dict(), 'model_weight_'+current_time.replace(':','-').replace(' ','-')+'.pth')
            tmp_txt = open('model_'+current_time.replace(':','-').replace(' ','-')+"_mapping_1.txt", "w", encoding='utf-8')
            for key, value in mapping_1.items():
                tmp_txt.write(str(key) + ": " + str(value) + "\n")
            tmp_txt.close()
            tmp_txt = open('model_'+current_time.replace(':','-').replace(' ','-')+"_mapping_2.txt", "w", encoding='utf-8')
            for key, value in mapping_2.items():
                tmp_txt.write(str(key) + ": " + str(value) + "\n")
            tmp_txt.close()
            tmp_txt = open('model_'+current_time.replace(':','-').replace(' ','-')+"_num_types.txt", "w", encoding='utf-8')
            tmp_txt.write(str(num_types) + "\n")
            tmp_txt.close()
        
        # test_data
        tmp_shape = test_data.shape[0]
        test_data = torch.from_numpy(np.float32(np.array(test_data)))
        test_data = torch.reshape(test_data, (tmp_shape, 224, 224))
        test_data = torch.unsqueeze(test_data, 1)
        test_data = torch.cat((test_data, test_data, test_data), dim=1)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        
        pcma_mygo_rc = pcma_mygo_rc.to('cpu')
        pcma_mygo_rc.eval()
        
        # close
        with torch.no_grad():
            outputs = pcma_mygo_rc(test_data)
        
        _, predicted = torch.max(outputs, 1)
        
        # from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
        
        # to NumPy
        predicted = predicted.numpy()
        test_label = test_label.numpy()
        
        # accuracy = accuracy_score(test_label, predicted)
        # 
        # precision = precision_score(test_label, predicted, average='macro', zero_division=1)
        # 
        # recall = recall_score(test_label, predicted, average='macro', zero_division=1)
        # 
        # print('accuracy')
        # print(accuracy)
        # print('precision')
        # print(precision)
        # print('recall')
        # print(recall)
        # print('macro f1 score')
        # print( 2*precision*recall / (precision+recall) )
        
        # accuracy
        accuracy = accuracy_score(test_label, predicted)  
        
        # macro precision 
        macro_precision = precision_score(test_label, predicted, average='macro', zero_division=1)  
        
        # macro recall  
        macro_recall = recall_score(test_label, predicted, average='macro', zero_division=1)  
        
        # macro F1 score
        macro_f1_score = f1_score(test_label, predicted, average='macro', zero_division=1)  
        
        print("Accuracy:", accuracy)  
        print("Macro Precision:", macro_precision)  
        print("Macro Recall:", macro_recall)  
        print("Macro F1 Score:", macro_f1_score)  
        self.train_accuracy = accuracy
        self.train_macro_f1_score = macro_f1_score
        
    def auto_anno_mldl_predict_0(self,ref_obj=None,test_obj=None,gpu_code = None,
                                 model_weight_load_path = None,
                                 dropout = None,
                                 mapping_1_load_path = None,
                                 mapping_2_load_path = None,
                                 num_types_load_path = None,
                                 predicted_results_save_path = None):
        the_dropout = self.set_the_value(input_value = dropout,default = 0.3)
        if self.model_weight is None and model_weight_load_path is None:
            raise PCmaster_anno_0_Error('No model weight! \n self.model_weight is None and model_weight_load_path is None!')
        if self.model_weight is not None and model_weight_load_path is not None:
            raise PCmaster_anno_0_Error('Only one is needed! \n self.model_weight and model_weight_load_path both exist!')
        if self.model_weight is None:
            self.model_weight = torch.load(model_weight_load_path)
            print('load model weight from '+ model_weight_load_path)
        if self.mapping_1 is None:
            self.mapping_1 = self.auto_annotation_with_deep_learning_load_mapping_0(mapping_load_path = mapping_1_load_path)
            print('load mapping_1 from '+mapping_1_load_path)
        if self.mapping_2 is None:
            self.mapping_2 = self.auto_annotation_with_deep_learning_load_mapping_0(mapping_load_path = mapping_2_load_path)
            print('load mapping_2 from '+mapping_2_load_path)
        if self.num_types is None:
            self.num_types = self.auto_annotation_with_deep_learning_load_num_types_0(num_types_load_path = num_types_load_path)
            print('load num_types from '+num_types_load_path)
        # print(self.model_weight)
        print(self.mapping_1)
        print(self.mapping_2)
        print(self.num_types)
        aaa_num_types = self.num_types
        
        import torchvision.models

        class PCMA_MYGO_RC(nn.Module):
            def __init__(self):
                super(PCMA_MYGO_RC, self).__init__()
                self.resnet18 = torchvision.models.resnet18(pretrained=False)   
                self.PCMA_MYGO_FC = nn.Sequential(
                                nn.ReLU(),
                                nn.Dropout(the_dropout),
                                nn.Linear(1000, aaa_num_types))
                
            def forward(self, x):
                x = self.resnet18(x)
                x = self.PCMA_MYGO_FC(x)
                return x
        
        pcma_mygo_rc = PCMA_MYGO_RC()
        pretrained_dict = self.model_weight
        model_dict = pcma_mygo_rc.state_dict()
        
        pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
        
        model_dict.update(pretrained_dict)
        pcma_mygo_rc.load_state_dict(model_dict)
        
        try:   
            df,_,_ = self.auto_annotation_with_deep_learning_transfer_adata_to_df_0(adata_input = test_obj)  
        except Exception:   
            df,_,_ = self.auto_annotation_with_deep_learning_transfer_adata_to_df_without_todense_0(adata_input = test_obj) 
        
        df_label = df.iloc[:, -1:]
        df_label_full = df.iloc[:, -2:]
        df = df.iloc[:, :-2]
        
        import math
        
        batch_size = 32 
        total_size = df.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = pcma_mygo_rc
        net_try = net_try.to('cpu')
        net_try.eval()
        
        with torch.no_grad():
            predicted_results = []
            for i in range(total_batches):
                start = i * batch_size
                end = min(start + batch_size, total_size)
                batch_data = df[start:end]
                batch_label = df_label[start:end]
                batch_data = torch.from_numpy(np.float32(np.array(batch_data)))
                batch_data = torch.reshape(batch_data, (batch_data.shape[0], 224, 224))
                batch_data = torch.unsqueeze(batch_data, 1)
                batch_data = torch.cat((batch_data, batch_data, batch_data), dim=1)
                batch_label = torch.from_numpy(np.int64(np.array(batch_label)).reshape(1,-1)[0])
                outputs = net_try(batch_data)
                _, predicted = torch.max(outputs, 1)
                predicted_results.append(predicted.numpy())
        
        predicted_results = np.concatenate(predicted_results)
        dl_cellnames = df_label_full['cellname'].tolist()
        original_celltypes = df_label_full['celltype'].tolist()
        # print(predicted_results)
        print(len(predicted_results))
        # print(dl_cellnames)
        print(len(dl_cellnames))
        # print(original_celltypes)
        print(len(original_celltypes))
        
        list2 = predicted_results.tolist()
        
        reverse_mapping_2 = {v: k for k, v in self.mapping_2.items()}
        
        list2_0_mapped = [reverse_mapping_2[value] for value in list2]
        
        list2_0 = [next(key for key, value in self.mapping_1.items() if value == mapped_value) for mapped_value in list2_0_mapped]
        
        # print(list2_0)
        
        predicted_results_true = list2_0
        
        self.predicted_results = predicted_results_true
        self.dl_cellnames = dl_cellnames
        self.original_celltypes = original_celltypes
        now = datetime.datetime.now()
        current_time = now.strftime("%Y-%m-%d %H:%M:%S")
        print(current_time)
        tmp_txt = open('model_'+current_time.replace(':','-').replace(' ','-')+"_predicted_results.txt", "w", encoding='utf-8')
        for list_item in self.predicted_results:
            tmp_txt.write(str(list_item) + "\n")
        tmp_txt.close()
        tmp_txt = open('model_'+current_time.replace(':','-').replace(' ','-')+"_dl_cellnames.txt", "w", encoding='utf-8')
        for list_item in self.dl_cellnames:
            tmp_txt.write(str(list_item) + "\n")
        tmp_txt.close()
        tmp_txt = open('model_'+current_time.replace(':','-').replace(' ','-')+"_original_celltypes.txt", "w", encoding='utf-8')
        for list_item in self.original_celltypes:
            tmp_txt.write(str(list_item) + "\n")
        tmp_txt.close()
        
    def auto_annotation_with_deep_learning_0(self,original_obj=None,ref_test_label=None,
                                           ref_test_label_for_ref=None,ref_test_label_for_test=None,
                                           random_sample_num=1000,gpu_code = None,
                                           learning_rate = None,epochs = None,
                                           batch_size = None,dropout = None,
                                           last_linear_n = None,num_workers = None,
                                           model_weight_save_path = None,
                                           mapping_1_save_path = None,
                                           mapping_2_save_path = None,
                                           num_types_save_path = None,
                                           model_weight_load_path = None,
                                           mapping_1_load_path = None,
                                           mapping_2_load_path = None,
                                           num_types_load_path = None,
                                           predicted_results_save_path = None,
                                           run_model = None):
        if run_model is None:
            run_model = 'train_and_predict'
        if model_weight_load_path is not None:
            run_model = 'predict'
        if ref_test_label is None:
            ref_test_label = 'ref_test_label'
        if ref_test_label_for_ref is None:
            ref_test_label_for_ref = 'ref'
        if ref_test_label_for_test is None:
            ref_test_label_for_test = 'test'
        ref_obj,test_obj = self.auto_annotation_with_deep_learning_split_dataset_0(original_obj=original_obj,ref_test_label=ref_test_label,
                                                                            ref_test_label_for_ref=ref_test_label_for_ref,
                                                                            ref_test_label_for_test=ref_test_label_for_test)
        if run_model == 'train':
            self.auto_anno_mldl_train_0(ref_obj=ref_obj,test_obj=test_obj,random_sample_num=random_sample_num,gpu_code = gpu_code,
                                        learning_rate = learning_rate,epochs = epochs,
                                        batch_size = batch_size,dropout = dropout,
                                        last_linear_n = last_linear_n,num_workers = num_workers,
                                        model_weight_save_path = model_weight_save_path,
                                        mapping_1_save_path = mapping_1_save_path,
                                        mapping_2_save_path = mapping_2_save_path,
                                        num_types_save_path = num_types_save_path)
        if run_model == 'train_and_predict':
            self.auto_anno_mldl_train_0(ref_obj=ref_obj,test_obj=test_obj,random_sample_num=random_sample_num,gpu_code = gpu_code,
                                        learning_rate = learning_rate,epochs = epochs,
                                        batch_size = batch_size,dropout = dropout,
                                        last_linear_n = last_linear_n,num_workers = num_workers,
                                        model_weight_save_path = model_weight_save_path,
                                        mapping_1_save_path = mapping_1_save_path,
                                        mapping_2_save_path = mapping_2_save_path,
                                        num_types_save_path = num_types_save_path)
            self.auto_anno_mldl_predict_0(ref_obj=ref_obj,test_obj=test_obj,gpu_code = gpu_code,
                                          model_weight_load_path = model_weight_load_path,
                                          dropout = dropout,
                                          mapping_1_load_path = mapping_1_load_path,
                                          mapping_2_load_path = mapping_2_load_path,
                                          num_types_load_path = num_types_load_path,
                                          predicted_results_save_path = predicted_results_save_path)
        if run_model == 'predict':
            self.auto_anno_mldl_predict_0(ref_obj=ref_obj,test_obj=test_obj,gpu_code = gpu_code,
                                          model_weight_load_path = model_weight_load_path,
                                          dropout = dropout,
                                          mapping_1_load_path = mapping_1_load_path,
                                          mapping_2_load_path = mapping_2_load_path,
                                          num_types_load_path = num_types_load_path,
                                          predicted_results_save_path = predicted_results_save_path)
        
    def auto_annotation_with_deep_learning_split_dataset_0(self,original_obj=None,ref_test_label=None,
                                                         ref_test_label_for_ref=None,ref_test_label_for_test=None):
        ref_obj = original_obj[original_obj.obs[ref_test_label] == ref_test_label_for_ref]
        test_obj = original_obj[original_obj.obs[ref_test_label] == ref_test_label_for_test]
        return ref_obj,test_obj
    
    def auto_annotation_with_deep_learning_load_mapping_0(self,mapping_load_path = None):
        tmp_txt = open(mapping_load_path, "r", encoding='utf-8')
        tmp_dict = {}
        for line in tmp_txt:
            tmp_dict[int(line.strip().split(": ")[0])] = line.strip().split(": ")[1]
        return tmp_dict
    
    def auto_annotation_with_deep_learning_load_num_types_0(self,num_types_load_path = None):
        tmp_txt = open(num_types_load_path, "r", encoding='utf-8')
        tmp_num_types = None
        for line in tmp_txt:
            tmp_num_types = line.strip()
        tmp_num_types = int(tmp_num_types)
        return tmp_num_types
    
    def set_the_value(self,input_value = None,default = None):
        if input_value is None:
            return default
        else:
            return input_value

    def auto_annotation_with_deep_learning_transfer_adata_to_df_for_all_0(self,adata_input = None,the_length = None):
        if the_length is None:
            the_length = 224
        try:   
            df,mapping_1,mapping_2 = self.auto_annotation_with_deep_learning_transfer_adata_to_df_0(
                adata_input = adata_input,the_length = the_length,
            )  
        except Exception:   
            df,mapping_1,mapping_2 = self.auto_annotation_with_deep_learning_transfer_adata_to_df_without_todense_0(
                adata_input = adata_input,the_length = the_length,
            )
        return df,mapping_1,mapping_2
        
    def auto_annotation_with_deep_learning_transfer_adata_to_df_0(self,adata_input = None,the_length = None):
        # normalize_total
        # sc.pp.normalize_total(adata_input, target_sum=1e4)
        if the_length is None:
            the_length = 224
        adata_input.obs['cellname'] = adata_input.obs_names
        
        tmp_pcm = adata_input
        data=tmp_pcm.X.todense()
        pcm_pd=pd.DataFrame(data,index=tmp_pcm.obs_names,columns=tmp_pcm.var_names)
        pcm_pd.loc[:,'cellname']=pcm_pd.index
        df = pcm_pd
        df.drop_duplicates(subset='cellname', keep='first', inplace=True)
        
        duplicate_rows = df.duplicated(subset='cellname', keep=False)
        df = df[~duplicate_rows]
        pcm_pd = df
        # print('pcm')
        # print(pcm_pd)
        # print('\n')
        
        label=pd.concat([adata_input.obs['cellname'], adata_input.obs['celltype']], axis=1)
        df = label
        num_categories = df['celltype'].nunique()
        mapping_1 = {category: i for i, category in enumerate(df['celltype'].unique())}
        # print('mapping_1')
        # print(mapping)
        # print('\n')
        df['celltype'] = df['celltype'].map(mapping_1)
        duplicate_rows = df.duplicated(subset='cellname', keep=False)
        df = df[~duplicate_rows]
        # print('df')
        # print(df)
        # print('\n')
        label = df
        
        merged = pd.merge(pcm_pd, label, on='cellname', how='inner')
        # print("merged")
        # print(merged)
        df = merged
        num_categories = df['celltype'].nunique()
        mapping_2 = {category: i for i, category in enumerate(df['celltype'].unique())}
        # print('mapping_2')
        # print(mapping)
        # print('\n')
        df['celltype'] = df['celltype'].map(mapping_2)
        
        # Create an all-zero DataFrame with the same number of rows as the original DataFrame, and the number of columns is the number of target columns to be expanded.
        # PCMA_MYGO_RC is used by default here. To change the data shape to 224*224=50176.
        # Here +2 is because the columns of celltype and cellname should be removed.
        # Because others' PCMA_MYGO_RC has trained tens of thousands of epoch, we use others' weights, and the convergence speed will be much faster than training from scratch.
        # The newly added line will be filled with 0 because it is meaningless, and then it will become 0.0001.
        # Here, making every value in the whole world not 0 takes into account the statement that data of 0 is not conducive to updating parameters.
        target_columns = the_length*the_length-df.shape[1]+2
        zeros_df = pd.DataFrame(np.zeros((df.shape[0], target_columns)), columns=[f'new_col_{i}' for i in range(target_columns)])
        
        expanded_df = pd.concat([df, zeros_df], axis=1)
        
        # print(expanded_df.shape) 
        
        df = expanded_df
        # print(df)
        columns = df.columns.tolist()
        
        columns.remove('cellname')
        columns.append('cellname')
        
        columns.remove('celltype')
        columns.append('celltype')
        
        # reorder
        df = df[columns]
        
        return df,mapping_1,mapping_2

    def auto_annotation_with_deep_learning_transfer_adata_to_df_without_todense_0(self,adata_input = None,the_length = None):
        # sc.pp.normalize_total(adata_input, target_sum=1e4)
        if the_length is None:
            the_length = 224
        
        adata_input.obs['cellname'] = adata_input.obs_names
        
        tmp_pcm = adata_input
        data=tmp_pcm.X
        pcm_pd=pd.DataFrame(data,index=tmp_pcm.obs_names,columns=tmp_pcm.var_names)
        pcm_pd.loc[:,'cellname']=pcm_pd.index
        df = pcm_pd
        df.drop_duplicates(subset='cellname', keep='first', inplace=True)
        
        duplicate_rows = df.duplicated(subset='cellname', keep=False)
        df = df[~duplicate_rows]
        pcm_pd = df
        # print('pcm')
        # print(pcm_pd)
        # print('\n')
        
        label=pd.concat([adata_input.obs['cellname'], adata_input.obs['celltype']], axis=1)
        df = label
        num_categories = df['celltype'].nunique()
        mapping_1 = {category: i for i, category in enumerate(df['celltype'].unique())}
        # print('mapping_1')
        # print(mapping)
        # print('\n')
        df['celltype'] = df['celltype'].map(mapping_1)
        duplicate_rows = df.duplicated(subset='cellname', keep=False)
        df = df[~duplicate_rows]
        # print('df')
        # print(df)
        # print('\n')
        label = df
        
        merged = pd.merge(pcm_pd, label, on='cellname', how='inner')
        # print("merged")
        # print(merged)
        df = merged
        num_categories = df['celltype'].nunique()
        mapping_2 = {category: i for i, category in enumerate(df['celltype'].unique())}
        # print('mapping_2')
        # print(mapping)
        # print('\n')
        df['celltype'] = df['celltype'].map(mapping_2)

        target_columns = the_length*the_length-df.shape[1]+2
        zeros_df = pd.DataFrame(np.zeros((df.shape[0], target_columns)), columns=[f'new_col_{i}' for i in range(target_columns)])
        
        expanded_df = pd.concat([df, zeros_df], axis=1)
        
        # print(expanded_df.shape)  
        
        df = expanded_df
        # print(df)
        columns = df.columns.tolist()
        
        columns.remove('cellname')
        columns.append('cellname')
        
        columns.remove('celltype')
        columns.append('celltype')
        
        df = df[columns]
        
        return df,mapping_1,mapping_2

    def save_df_to_images(self,root_dir = None,df = None,mapping_1 = None,mapping_2 = None,the_length = None):
        if the_length is None:
            the_length = 224
        
        start_time = time.time()
        if root_dir is None or df is None or mapping_1 is None or mapping_2 is None:
            raise PCmaster_anno_0_Error("root_dir == None or df == None or mapping_1 == None or mapping_2 == None.")
        # root_dir = './data_to_image_test/'
        reverse_mapping_1 = {v: k for k, v in mapping_1.items()}  
        reverse_mapping_2 = {v: k for k, v in mapping_2.items()} 
        # Read the DataFrame containing gene expression data  
        # df = pd.read_csv('your_dataframe.csv')  
        the_max = df.iloc[:, :the_length*the_length].values.max()
        print(the_max)
        # Create a MinMaxScaler object  
        # scaler = MinMaxScaler(feature_range=(0, 255))  
        print('ok')
        # Iterate over each cell  
        for index in range(df.shape[0]):
            row = df.iloc[index]
            # Extract the gene expression data  
            gene_expression = row[:the_length*the_length].values  
            # print(gene_expression)
            # the_max = gene_expression.max()
            # gene_expression = gene_expression*255.0/the_max
            
            gene_expression = gene_expression.reshape(1,-1)
            # print(gene_expression)
            # break

            # Scale the gene_expression to have values between 0 and 255  
            # scaled_gene_expression = scaler.fit_transform(gene_expression)  
            scaled_gene_expression = gene_expression  
            
            # Reshape the gene expression data into a 224x224 matrix  
            matrix = scaled_gene_expression.reshape(the_length, the_length)  
            scaled_matrix = matrix
            
            # Convert the scaled matrix to uint8 data type  
            scaled_matrix = scaled_matrix.astype(np.float32)  
            
            # Copy and expand the matrix into a 224x224x3 matrix  
            tmp_matrix = np.repeat(scaled_matrix[:, :, np.newaxis], 3, axis=2)  
            
            # Create an image object  
            image = Image.fromarray(tmp_matrix)  
            
            # Get the celltype and cellname of the cell  
            cell_type = row[the_length*the_length+1]
            tmp_re_1 = reverse_mapping_2.get(cell_type)
            cell_type_origin = reverse_mapping_1.get(tmp_re_1)
            cell_type = cell_type_origin.replace(' ','___').replace('/','---')
            # print('celltype')
            # print(cell_type)
            
            cell_name = str(row[the_length*the_length])  
            # print('cellname')
            # print(cell_name)
            
            # Create the directory to save the image (if it doesn't exist)  
            save_dir = root_dir + str(cell_type)  
            os.makedirs(save_dir, exist_ok=True)  
            
            # Save the image as a JPEG file  
            image.save(f'{save_dir}/{cell_name}.jpg')
            # break
        
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"time cost：{execution_time} s")

    def save_df_to_npys(self,root_dir = None,df = None,mapping_1 = None,mapping_2 = None,the_length = None):
        if the_length is None:
            the_length = 224
        
        start_time = time.time()
        if root_dir is None or df is None or mapping_1 is None or mapping_2 is None:
            raise PCmaster_anno_0_Error("root_dir == None or df == None or mapping_1 == None or mapping_2 == None.")
        # root_dir = './data_to_image_test/'
        reverse_mapping_1 = {v: k for k, v in mapping_1.items()}  
        reverse_mapping_2 = {v: k for k, v in mapping_2.items()} 
        # Read the DataFrame containing gene expression data  
        # df = pd.read_csv('your_dataframe.csv')  
        the_max = df.iloc[:, :the_length*the_length].values.max()
        print(the_max)
        # Create a MinMaxScaler object  
        # scaler = MinMaxScaler(feature_range=(0, 255))  
        print('ok')
        # Iterate over each cell  
        for index in range(df.shape[0]):
            row = df.iloc[index]
            # Extract the gene expression data  
            gene_expression = row[:the_length*the_length].values  
            # print(gene_expression)
            # the_max = gene_expression.max()
            # gene_expression = gene_expression*255.0/the_max
            
            # gene_expression = gene_expression.reshape(1,-1)
            # print(gene_expression)
            # break

            # Scale the gene_expression to have values between 0 and 255  
            # scaled_gene_expression = scaler.fit_transform(gene_expression)  
            # scaled_gene_expression = gene_expression  
            
            # Reshape the gene expression data into a 224x224 matrix  
            # matrix = scaled_gene_expression.reshape(the_length, the_length)  
            # scaled_matrix = matrix
            
            # Convert the scaled matrix to uint8 data type  
            # scaled_matrix = scaled_matrix.astype(np.uint8)  
            
            # Copy and expand the matrix into a 224x224x3 matrix  
            # tmp_matrix = np.repeat(scaled_matrix[:, :, np.newaxis], 3, axis=2)  
            
            # Create an image object  
            # image = Image.fromarray(tmp_matrix)  
            
            # Get the celltype and cellname of the cell  
            cell_type = row[the_length*the_length+1]
            tmp_re_1 = reverse_mapping_2.get(cell_type)
            cell_type_origin = reverse_mapping_1.get(tmp_re_1)
            cell_type = cell_type_origin.replace(' ','___').replace('/','---')
            # print('celltype')
            # print(cell_type)
            
            cell_name = str(row[the_length*the_length])  
            # print('cellname')
            # print(cell_name)
            
            # Create the directory to save the image (if it doesn't exist)  
            save_dir = root_dir + str(cell_type)  
            os.makedirs(save_dir, exist_ok=True)  
            
            # Save the image as a JPEG file  
            # image.save(f'{save_dir}/{cell_name}.jpg')
            np.save(f'{save_dir}/{cell_name}.npy', gene_expression)
            # break
        
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"time cost：{execution_time} s")







    # other functions
    
    def calculate_median_f1_score(self,int_list1=None, int_list2=None):
        # from sklearn.metrics import f1_score  
        # import numpy as np  
        f1_scores = f1_score(int_list1, int_list2, average=None)  
    
        # Calculate median F1 score  
        median_f1_score = np.median(f1_scores)  
          
        print('f1-scores')
        print(f1_scores)
        print('median_f1_score')
        print(median_f1_score)
        
        return f1_scores,median_f1_score
    
    def get_adata_gene_list_irgsp(self,the_adata = None):
        adata_gene_list_zh11 = []
        adata_gene_list_irgsp = []
        for i in the_adata.var_names:
            adata_gene_list_zh11.append(i.split('.')[0].strip())
        for i in adata_gene_list_zh11:
            if i in zh11_irgsp_dict.keys():
                adata_gene_list_irgsp.append(zh11_irgsp_dict[i])
            else:
                adata_gene_list_irgsp.append(i)
        # for i in adata_gene_list_irgsp:
        #     print(i)
        return adata_gene_list_irgsp
    
    def compare_two_adata_gene_names(self,adata1 = None,adata2 = None):
        set1 = set()
        set2 = set()
        for i in adata1.var_names:
            set1.add(str(i))
        for i in adata2.var_names:
            set2.add(str(i))
        count = 0
        for i in set1:
            if i in set2:
                count = count+1
        print(count)
    
    # GPU
    def evaluate_accuracy_gpu(self,net, data_iter, device=None): 
        if isinstance(net, nn.Module):
            net.eval()
            if not device:
                device = next(iter(net.parameters())).device
        metric = d2l.Accumulator(2)
        with torch.no_grad():
            for X, y in data_iter:
                if isinstance(X, list):
                    X = [x.to(device) for x in X]
                else:
                    X = X.to(device)
                y = y.to(device)
                metric.add(d2l.accuracy(net(X), y), y.numel())
        return metric[0] / metric[1]
    
    # GPU
    def train_with_GPU(self,net, train_iter, valid_iter, num_epochs, lr, device):
        def init_weights(m):
            if type(m) == nn.Linear or type(m) == nn.Conv2d:
                nn.init.xavier_uniform_(m.weight)
        net.apply(init_weights)
        print('training on', device)
        net.to(device)
        optimizer = torch.optim.SGD(net.parameters(), lr=lr, momentum=0.9)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=150, eta_min=0.00001)
        loss = nn.CrossEntropyLoss()
        animator = d2l.Animator(xlabel='epoch', xlim=[1, num_epochs],ylim=[0,1],
                                legend=['train loss', 'train acc', 'valid acc'])
        timer, num_batches = d2l.Timer(), len(train_iter)
        for epoch in range(num_epochs):
            metric = d2l.Accumulator(3)
            net.train()
            for i, (X, y) in enumerate(train_iter):
                timer.start()
                optimizer.zero_grad()
                X, y = X.to(device), y.to(device)
                y_hat = net(X)
                l = loss(y_hat, y)
                l.backward()
                optimizer.step()
                scheduler.step()
                with torch.no_grad():
                    metric.add(l * X.shape[0], d2l.accuracy(y_hat, y), X.shape[0])
                timer.stop()
                train_l = metric[0] / metric[2]
                train_acc = metric[1] / metric[2]
                if (i + 1) % (num_batches // 5) == 0 or i == num_batches - 1:
                    animator.add(epoch + (i + 1) / num_batches,
                                 (train_l, train_acc, None))
            valid_acc = evaluate_accuracy_gpu(net, valid_iter)
            animator.add(epoch + 1, (None, None, valid_acc))
        print(f'loss {train_l:.3f}, train acc {train_acc:.3f}, '
              f'valid acc {valid_acc:.3f}')
        print(f'{metric[2] * num_epochs / timer.sum():.1f} examples/sec '
              f'on {str(device)}')
    
    def get_sequence_index_dict(self,n = None, input_dict = None):
        # Given sequence length  
        # Given input dictionary   
          
        # Getting the original index  
        origin_index = len(input_dict.keys())  
          
        # Generating all possible sequences  
        bases = ['A', 'T', 'C', 'G']  # Assuming the sequence can consist of A, T, C, G  
        all_sequences = [''.join(seq) for seq in itertools.product(bases, repeat=n)]  
          
        # Creating the new dictionary with sequences and indices  
        new_dict = {sequence: origin_index + index for index, sequence in enumerate(all_sequences)}  
          
        # Merging the new dictionary with the input dictionary  
        input_dict.update(new_dict)  
    
    def get_gene_sequence_folder(self,input_gene_fa_path = None,output_folder = None):
        if output_folder is None:
            output_folder = input_gene_fa_path.replace('.fa','')
        if not os.path.exists(output_folder):  
            os.makedirs(output_folder) 
        with open(input_gene_fa_path, 'r',encoding='utf-8') as the_fa_file:   
            the_dict = {}
            tmp_gene = ''
            tmp_sequence = ''
            for line in the_fa_file: 
                if line.startswith('>'):
                    if tmp_gene != '': # save old,reset
                        if os.path.exists(output_folder+'/'+tmp_gene+'.txt'): # add a new line
                            with open(output_folder+'/'+tmp_gene+'.txt', 'a',encoding='utf-8') as tmp_file:  
                                # Write the string to the file  
                                tmp_file.write('\n'+tmp_sequence)
                        else:
                            with open(output_folder+'/'+tmp_gene+'.txt', 'w',encoding='utf-8') as tmp_file:  
                                # Write the string to the file  
                                tmp_file.write(tmp_sequence) 
                        # the_dict[tmp_gene]=tmp_sequence
                        tmp_gene = line.split(' ')[1].strip().replace('\n','')
                        tmp_sequence = ''
                    else: # the_start
                        tmp_gene = line.split(' ')[1].strip().replace('\n','')
                        tmp_sequence = ''
                else: # append the gene sequence
                    tmp_sequence = tmp_sequence+str(line.strip().replace('\n',''))
                    # print(line)  
            # return the_dict
        
    def get_gene_sequence_vector(self,input_folder=None, output_folder=None, vector_len=None, seq_dict=None):  
        tmp_dict = {i: 0 for i in range(vector_len)}  # Step 1  
        seq_count_dict = {j:0.0 for j in seq_dict.keys()}
        gene_seq_num_list = []
        count_count = 0
        if not os.path.exists(output_folder):  
            os.makedirs(output_folder) 
        # import os  
        for file_name in os.listdir(input_folder):  
            if file_name.endswith(".txt"):  
                file_path = os.path.join(input_folder, file_name)  
                with open(file_path, 'r') as file:  
                    lines = file.readlines()  
                    line_count = len(lines)
                    for line in lines:
                        # print(line)
                        for i in range(len(line) - 5):  # Step 2.1  
                            sequence = line[i:i+5]  
                            if sequence in seq_dict:  
                                pos = seq_dict[sequence]  
                                tmp_dict[pos] = tmp_dict[pos] +1  
                                seq_count_dict[sequence] = seq_count_dict[sequence]+1.0/float(line_count)
                                # seq_count_dict[sequence] = tmp_dict[pos]
                        for i in range(len(line) - 6):  # Step 2.1  
                            sequence = line[i:i+6]  
                            if sequence in seq_dict:  
                                pos = seq_dict[sequence]  
                                tmp_dict[pos] = tmp_dict[pos] +1  
                                seq_count_dict[sequence] = seq_count_dict[sequence]+1.0/float(line_count)
                                # seq_count_dict[sequence] = tmp_dict[pos]
                        for i in range(len(line) - 7):  # Step 2.1  
                            sequence = line[i:i+7]  
                            if sequence in seq_dict:  
                                pos = seq_dict[sequence]  
                                tmp_dict[pos] = tmp_dict[pos] +1  
                                seq_count_dict[sequence] = seq_count_dict[sequence]+1.0/float(line_count)
                                # seq_count_dict[sequence] = tmp_dict[pos]
                line_count = len(lines)  # Step 2.2  
                gene_seq_num_list.append(line_count)
                # print(len(lines))
                # for seqseq in seq_count_dict.keys():
                #     print(seqseq)
                #     print(seq_count_dict[seqseq])
                tmp_dict_values = list(tmp_dict.values())  
                tmp_dict_values = [float(val)/float(line_count) for val in tmp_dict_values]  
                output_file_path = os.path.join(output_folder, file_name)  
                with open(output_file_path, 'w') as output_file:  
                    output_file.write('\n'.join(map(str, tmp_dict_values)))  # Step 2.3  
                count_count = count_count+1
            # break
            # if count_count >=50:
            #     break
        # print('len of gene_seq_num_list')
        # gene_seq_num_list_len = len(gene_seq_num_list)
        # print(gene_seq_num_list_len)
        # print('max min mean median q1 q3 of gene_seq_num_list')
        # max_value = max(gene_seq_num_list)  
        # min_value = min(gene_seq_num_list)  
        # mean_value = statistics.mean(gene_seq_num_list)  
        # median_value = statistics.median(gene_seq_num_list)  
        # q1 = statistics.quantiles(gene_seq_num_list, n=4)[0]  
        # q3 = statistics.quantiles(gene_seq_num_list, n=4)[-1]  
        # print(max_value)  
        # print(min_value)  
        # print(mean_value)  
        # print(median_value)  
        # print(q1)  
        # print(q3)  
        # for seqseq in seq_count_dict.keys():
        #     print(seqseq)
        #     print(seq_count_dict[seqseq]/gene_seq_num_list_len)
        return gene_seq_num_list,seq_count_dict
    
    def get_all_files_under_a_folder(self,input_folder = None,end_str = None):
        file_list = [os.path.join(input_folder, file) for file in os.listdir(input_folder) if file.endswith(end_str)]  
        return file_list
    
    def get_all_biotypes_from_gtf(self,input_gtf = None):
        tmp_set = set()
        with open(input_gtf, 'r',encoding='utf-8') as the_gtf_file:
            for line in the_gtf_file:
                tmp_set.add(line.split('biotype ')[1].split(';')[0].replace('\"','').strip())
        print(input_gtf)
        for i in tmp_set:
            print(i)
        print('\n')
        
    def change_gtf_to_specific_biotype(self,input_gtf = None,output_gtf = None,the_biotype = None):
        tmp_set = set()
        with open(input_gtf, 'r',encoding='utf-8') as the_gtf_file:
            for line in the_gtf_file:
                if str(line.split('biotype ')[1].split(';')[0].replace('\"','').strip()) == the_biotype:
                    tmp_set.add(line.strip())
        print(input_gtf)
        print(len(tmp_set))
        
    def save_big_adata_to_small_adata(self,input_adata_path=None, input_adata_name=None, small_adata_size=None, output_folder=None):
        if not os.path.exists(output_folder):  
            os.makedirs(output_folder) 
        max_count = 0
        big_adata = sc.read(input_adata_path + input_adata_name)
        try:
            del big_adata.raw
        except Exception:
            max_count = 0
        # Define the start and end cell indices  
        tmp_list = []
        max_count = 0
        for i in range(0,len(big_adata.obs_names)):
            tmp_count = i//small_adata_size+1
            tmp_list.append(tmp_count)
            if tmp_count>max_count:
                max_count = tmp_count
        big_adata.obs['split_group_number'] = tmp_list
        # start_cell_index = 0  
        # end_cell_index = start_cell_index + small_adata_size
        # Loop through the big Anndata and save smaller parts 
        end_count = 1
        while end_count <= max_count - 1:  
            # Extract a smaller part of the big Anndata  
            tmp_adata = big_adata[big_adata.obs['split_group_number'] == end_count, :]  
            tmp_adata.obs.reset_index(inplace=True)  
            tmp_adata.var.reset_index(inplace=True) 
            # Save the smaller Anndata  
            tmp_adata.write(output_folder + input_adata_name.replace('.h5ad','') + \
                            f"_from_{small_adata_size*(end_count-1)}_to_{small_adata_size*(end_count)}.h5ad")  
            end_count = end_count+1
            # Update indices for the next iteration  
            # start_cell_index = end_cell_index  
            # end_cell_index = min(start_cell_index + small_adata_size, big_adata.shape[0])  
            
    def change_adata_to_cell_vector(self,input_adata=None,gene_seq_vector_folder=None, 
                                    output_folder=None, if_change_var_names=None,
                                    the_length=None,):
        max_count = 0
        gene_seq_vector_dict = get_gene_seq_vector_dict(gene_seq_vector_folder=gene_seq_vector_folder)
        the_adata = sc.read(input_adata)
        # try:
        #     del the_adata.raw
        # except Exception:
        #     max_count = 0
        if if_change_var_names is True:
            the_adata.var_names = the_adata.var['features']
        if the_length is None:
            the_length = 50176
        print(the_adata.var_names)
        # print(the_adata.var['features'])
        the_adata.obs['cellname_tmp'] = the_adata.obs_names
        the_adata.var['genename_tmp'] = the_adata.var_names
        for i in the_adata.obs_names:
            tmp_vector = np.zeros(the_length)
            tmp_adata = the_adata[the_adata.obs['cellname_tmp']==i,:]
            print(len(tmp_adata.var_names))
            sc.pp.filter_genes(tmp_adata, min_cells=1)
            print(len(tmp_adata.var_names))
            for j in tmp_adata.var_names:
                # print(j)
                if j in gene_seq_vector_dict.keys():
                    # print(j)
                    # with open(gene_seq_vector_dict[j], 'r',encoding='utf-8') as file:
                    #     value_list = []
                    #     for line in file:
                    #         value_list.append(float(line.strip()))
                    #     data_array = np.array(value_list)
                        # print(data_array)
                        # print(type(tmp_adata[:,tmp_adata.var['genename_tmp']==j].X[0,0])) 
                    data_array = np.load(gene_seq_vector_dict[j])
                    mask_log10 = (np.arange(the_length) >= 1) & (np.arange(the_length) <= 1024)  
                    mask_loge = (np.arange(the_length) >= 1025) & (np.arange(the_length) <= 5120)  
                      
                    # Apply the logarithm operations selectively  
                    data_array[mask_log10] = data_array[mask_log10]+1
                    data_array[mask_loge] = data_array[mask_loge]+1
                    data_array[mask_log10] = np.log10(data_array[mask_log10])  
                    data_array[mask_loge] = np.log(data_array[mask_loge])  
                      
                    # Print the resulting array to verify  
                    # print(array)  
                    tmp_value = tmp_adata[:,tmp_adata.var['genename_tmp']==j].X[0,0]
                    if tmp_value>0.0:
                        tmp_vector = tmp_vector+data_array*tmp_value
                # break
            array_sum = np.sum(tmp_vector)
            tmp_vector = tmp_vector * 100000 / array_sum
            tmp_vector = tmp_vector+1
            tmp_vector = np.log(tmp_vector)
            # print(tmp_vector)
            break
            
    def change_adata_to_cell_vector_2(self,input_adata=None,gene_seq_vector_folder=None, 
                                    output_folder=None, if_change_var_names=None,
                                    the_length=None,):
        max_count = 0
        gene_seq_vector_dict = get_gene_seq_vector_dict(gene_seq_vector_folder=gene_seq_vector_folder)
        the_adata = sc.read(input_adata)
        # try:
        #     del the_adata.raw
        # except Exception:
        #     max_count = 0
        if if_change_var_names is True:
            the_adata.var_names = the_adata.var['features']
        if the_length is None:
            the_length = 50176
        print(the_adata.var_names)
        # print(the_adata.var['features'])
        the_adata.obs['cellname_tmp'] = the_adata.obs_names
        the_adata.var['genename_tmp'] = the_adata.var_names
        for i in the_adata.obs_names:
            tmp_vector = np.zeros(the_length)
            tmp_adata = the_adata[the_adata.obs['cellname_tmp']==i,:]
            print(len(tmp_adata.var_names))
            # sc.pp.filter_genes(tmp_adata, min_cells=1)
            print(len(tmp_adata.var_names))
            for j in tmp_adata.var_names:
                # print(j)
                tmp_value = tmp_adata[:,tmp_adata.var['genename_tmp']==j].X[0,0]
                if j in gene_seq_vector_dict.keys() and tmp_value>0.0:
                    # print(j)
                    # with open(gene_seq_vector_dict[j], 'r',encoding='utf-8') as file:
                    #     value_list = []
                    #     for line in file:
                    #         value_list.append(float(line.strip()))
                    #     data_array = np.array(value_list)
                        # print(data_array)
                        # print(type(tmp_adata[:,tmp_adata.var['genename_tmp']==j].X[0,0])) 
                    data_array = np.load(gene_seq_vector_dict[j])
                    mask_log10 = (np.arange(the_length) >= 1) & (np.arange(the_length) <= 1024)  
                    mask_loge = (np.arange(the_length) >= 1025) & (np.arange(the_length) <= 5120)  
                      
                    # Apply the logarithm operations selectively  
                    data_array[mask_log10] = data_array[mask_log10]+1
                    data_array[mask_loge] = data_array[mask_loge]+1
                    data_array[mask_log10] = np.log10(data_array[mask_log10])  
                    data_array[mask_loge] = np.log(data_array[mask_loge])  
                      
                    # Print the resulting array to verify  
                    # print(array)  
                    # tmp_value = tmp_adata[:,tmp_adata.var['genename_tmp']==j].X[0,0]
                    # if tmp_value>0.0:
                    tmp_vector = tmp_vector+data_array*tmp_value
                # break
            array_sum = np.sum(tmp_vector)
            tmp_vector = tmp_vector * 100000 / array_sum
            tmp_vector = tmp_vector+1
            tmp_vector = np.log(tmp_vector)
            # print(tmp_vector)
            break
    
    def get_gene_seq_vector_dict(self,gene_seq_vector_folder=None):
        tmp_dict = {}
        files = os.listdir(gene_seq_vector_folder)
        for file in files:
            if file.endswith('.npy'):
                tmp_dict[file.replace('.npy','')] = gene_seq_vector_folder+file
        return tmp_dict
    
    def folder_merge(self,folders=None, output_folder=None):
        for folder in folders:  
            files = os.listdir(folder)  
            for file in files:  
                shutil.copy(os.path.join(folder, file), os.path.join(output_folder, file))  
                
    def change_txt_vector_to_npy_vector(self,input_folder=None,output_folder=None):
        files = os.listdir(input_folder)
        for the_file in files:
            with open(os.path.join(input_folder, the_file), 'r',encoding='utf-8') as file:
                value_list = []
                for line in file:
                    value_list.append(float(line.strip()))
                data_array = np.array(value_list)
                np.save(os.path.join(output_folder, the_file.replace('.txt','.npy')),data_array)
                
    def compare_load_content_and_speed(self,input_txt=None, input_npy=None):
        with open(input_txt, 'r',encoding='utf-8') as file:
            value_list = []
            for line in file:
                value_list.append(float(line.strip()))
            data_array1 = np.array(value_list)
            data_array2 = np.load(input_npy)
            # Check if the arrays have the same shape and elements  
            result = np.array_equal(data_array1, data_array2)  
            print(result)  # This will print True if the arrays have the same values, and False otherwise 
            
        # Process 1  
        start_time_p1 = time.time()  
        with open(input_txt, 'r',encoding='utf-8') as file:
            value_list = []
            for line in file:
                value_list.append(float(line.strip()))
            data_array1 = np.array(value_list)  
        end_time_p1 = time.time()  
        time_taken_p1 = end_time_p1 - start_time_p1  
        print(f"Time taken by process 1: {time_taken_p1} seconds")  
          
        # Process 2  
        start_time_p2 = time.time()  
        data_array2 = np.load(input_npy) 
        end_time_p2 = time.time()  
        time_taken_p2 = end_time_p2 - start_time_p2  
        print(f"Time taken by process 2: {time_taken_p2} seconds")  
        
    def change_npy_vector_length(self,input_folder=None,the_length=None, start_pos=None, end_pos=None, output_folder=None):
        if not os.path.exists(output_folder):  
            os.makedirs(output_folder)
        if the_length is None:
            the_length = 50176
        files = os.listdir(input_folder)
        for the_file in files:
            if the_file.endswith('.npy'):
                data_array = np.load(os.path.join(input_folder, the_file))
                the_mask = (np.arange(the_length) >= start_pos) & (np.arange(the_length) <= end_pos)
                np.save(os.path.join(output_folder, the_file),data_array[the_mask])
                
    def change_adata_to_cell_vector_3(self,input_adata=None,gene_seq_vector_folder=None, 
                                      output_folder=None, if_change_var_names=None,
                                      the_length=None,# all_gene_seq_possibilities
                                      gene_batch_size=None,gene_seq_dict=None,obs_celltype_label=None):
        max_count = 0
        gene_seq_vector_dict = get_gene_seq_vector_dict(gene_seq_vector_folder=gene_seq_vector_folder)
        the_adata = sc.read(input_adata)
        # try:
        #     del the_adata.raw
        # except Exception:
        #     max_count = 0
        if if_change_var_names is True:
            the_adata.var_names = the_adata.var['features']
        if the_length is None:
            the_length = 50176
        if gene_batch_size is None:
            gene_batch_size = 1000
        print(the_adata.var_names)
        # print(the_adata.var['features'])
        # the_adata.obs['cellname_tmp'] = the_adata.obs_names
        # the_adata.var['genename_tmp'] = the_adata.var_names
        sc.pp.filter_genes(the_adata, min_cells=1)
        the_number_of_cells = len(the_adata.obs_names)
        var_len = len(the_adata.var_names)
        start_len = 0
        result_matrix = np.zeros((the_number_of_cells, the_length))
        matched_count = 0
        miss_count = 0
        miss_list = []
        while start_len<var_len:
            matrix2_blocks = []
            tmp_list = []
            # tmp_vector = np.zeros(the_length)
            for i in range(start_len,min(start_len+gene_batch_size,var_len)):
                tmp_gene = the_adata.var_names[i]
                tmp_list.append(tmp_gene)
                if tmp_gene in gene_seq_vector_dict.keys(): 
                    matrix2_blocks.append(np.load(gene_seq_vector_dict[tmp_gene]))
                    matched_count = matched_count+1
                else:
                    matrix2_blocks.append(np.zeros(the_length))
                    miss_count = miss_count+1
                    miss_list.append(tmp_gene)
            matrix1 = the_adata[:,tmp_list].X.todense()
            matrix2 = np.vstack(matrix2_blocks) 
            print(matrix1.shape)
            print(matrix2.shape)
            result_matrix = result_matrix + the_adata[:,tmp_list].X.todense() @ matrix2
            print(result_matrix.shape)
            start_len = start_len+gene_batch_size
        print(matched_count)
        print(miss_count)
        # for i in miss_list:
        #     print(i)
        result_matrix = result_matrix[:, :21504]
        the_df = pd.DataFrame(result_matrix, index=the_adata.obs_names, columns=list(the_dict.keys()))
        return the_df,miss_list,the_adata.obs[obs_celltype_label]
    
    def change_adata_to_cell_vector_4(self,input_adata_dir=None,input_adata_name=None,gene_seq_vector_folder=None, 
                                      output_folder=None, if_change_var_names=None,if_change_obs_names=None,
                                      the_length=None,result_matrix_length=None,# all_gene_seq_possibilities,
                                      gene_batch_size=None,gene_seq_dict=None,obs_celltype_label=None):
        max_count = 0
        gene_seq_vector_dict = self.get_gene_seq_vector_dict(gene_seq_vector_folder=gene_seq_vector_folder)
        the_adata = sc.read(input_adata_dir+input_adata_name)
        # try:
        #     del the_adata.raw
        # except Exception:
        #     max_count = 0
        if result_matrix_length is None:
            result_matrix_length=21504
        if output_folder is None:
            output_folder=input_adata_dir
        if not os.path.exists(output_folder):  
            os.makedirs(output_folder)
        if if_change_var_names is True:
            the_adata.var_names = the_adata.var['features']
        if if_change_obs_names is True:
            the_adata.obs_names = list(the_adata.obs['index'])
        if the_length is None:
            the_length = 50176
        if gene_batch_size is None:
            gene_batch_size = 1000
        print(the_adata.var_names)
        # print(the_adata.var['features'])
        # the_adata.obs['cellname_tmp'] = the_adata.obs_names
        # the_adata.var['genename_tmp'] = the_adata.var_names
        sc.pp.filter_genes(the_adata, min_cells=1)
        the_number_of_cells = len(the_adata.obs_names)
        var_len = len(the_adata.var_names)
        start_len = 0
        result_matrix = np.zeros((the_number_of_cells, the_length))
        matched_count = 0
        miss_count = 0
        miss_list = []
        while start_len<var_len:
            matrix2_blocks = []
            tmp_list = []
            # tmp_vector = np.zeros(the_length)
            for i in range(start_len,min(start_len+gene_batch_size,var_len)):
                tmp_gene = the_adata.var_names[i]
                tmp_list.append(tmp_gene)
                if tmp_gene in gene_seq_vector_dict.keys(): 
                    matrix2_blocks.append(np.load(gene_seq_vector_dict[tmp_gene]))
                    matched_count = matched_count+1
                else:
                    matrix2_blocks.append(np.zeros(the_length))
                    miss_count = miss_count+1
                    miss_list.append(tmp_gene)
            matrix1 = the_adata[:,tmp_list].X.todense()
            matrix2 = np.vstack(matrix2_blocks) 
            print(matrix1.shape)
            print(matrix2.shape)
            result_matrix = result_matrix + the_adata[:,tmp_list].X.todense() @ matrix2
            print(result_matrix.shape)
            start_len = start_len+gene_batch_size
        print(matched_count)
        print(miss_count)
        # for i in miss_list:
        #     print(i)
        result_matrix = result_matrix[:, :result_matrix_length]
        the_df = pd.DataFrame(result_matrix, index=the_adata.obs_names, columns=list(gene_seq_dict.keys()))
        the_df.to_csv(output_folder+input_adata_name.replace('.h5ad','_changed.csv'))  
        del the_df
        gc.collect()
        new_anndata = sc.read_csv(output_folder+input_adata_name.replace('.h5ad','_changed.csv'))  
        true_gene_list = list(new_anndata.var_names)[0:]
        print(true_gene_list[0])
        new_anndata.var['seq_names']=new_anndata.var_names
        new_anndata=new_anndata[:,new_anndata.var['seq_names'].isin(true_gene_list)]
        new_anndata.obs['celltype']=the_adata.obs[obs_celltype_label]
        
        print(new_anndata.obs_names)
        print(new_anndata.var_names)
        print(new_anndata.X.max())
        print(new_anndata.obs['celltype'])
        new_anndata.write(output_folder+input_adata_name.replace('.h5ad','_changed.h5ad'))
        
        # return the_df,miss_list,the_adata.obs[obs_celltype_label]

    def downsample_celltypes(self, new_anndata=None, n=None, celltype_label='celltype',random_seed=42):
        # downsample each celltype, max is n
        np.random.seed(random_seed)
        celltype_counts = new_anndata.obs[celltype_label].value_counts()
        re_anndata = None
        for celltype, count in celltype_counts.items():
            if count > n:
                indices_to_keep = np.random.choice(
                    new_anndata.obs.index[new_anndata.obs[celltype_label] == celltype],
                    size=n,  
                    replace=False,)
                if re_anndata is None:
                    re_anndata = new_anndata[new_anndata.obs.index.isin(indices_to_keep)]
                else:
                    re_anndata = re_anndata.concatenate(new_anndata[new_anndata.obs.index.isin(indices_to_keep)], index_unique=None)
            else:
                if re_anndata is None:
                    re_anndata = new_anndata[new_anndata.obs.index.isin(indices_to_keep)]
                else:
                    re_anndata = re_anndata.concatenate(new_anndata[new_anndata.obs[celltype_label] == celltype], index_unique=None)
                    
        return re_anndata

    def split_anndata_by_celltype_to_train_valid_test(self, new_anndata=None, celltype_label='celltype', split_ratio=[0.8,0.1,0.1], random_seed=42):  
        # split adata to three part by given ratio
        # import scanpy as sc  
        # import numpy as np 
        np.random.seed(random_seed)
        train_anndata = None
        valid_anndata = None
        test_anndata = None
        
        for celltype in new_anndata.obs[celltype_label].unique(): 
            celltype_indices = new_anndata.obs[celltype_label] == celltype  
            celltype_data = new_anndata[celltype_indices]  
            
            total_cells = len(celltype_data)  
            x, y, z = split_ratio  
            x_count = int(total_cells * x / (x + y + z))  
            y_count = int(total_cells * y / (x + y + z))  
            z_count = total_cells - x_count - y_count  
            
            indices = np.random.permutation(total_cells)  
            train_indices = indices[:x_count]  
            valid_indices = indices[x_count:x_count + y_count]  
            test_indices = indices[x_count + y_count:]  
            
            train_data = celltype_data[train_indices]  
            valid_data = celltype_data[valid_indices]  
            test_data = celltype_data[test_indices]  
            
            if train_anndata is None:
                train_anndata = train_data
            else:
                train_anndata = train_anndata.concatenate(train_data, index_unique=None)  
            if valid_anndata is None:
                valid_anndata = valid_data
            else:
                valid_anndata = valid_anndata.concatenate(valid_data, index_unique=None)   
            if test_anndata is None:
                test_anndata = test_data
            else:
                test_anndata = test_anndata.concatenate(test_data, index_unique=None)   
            
        return train_anndata, valid_anndata, test_anndata  

    def save_to_pkl(obj=None, name=None):
        with open(name,'wb') as f:
            pickle.dump(obj,f)

    def load_pkl(name=None):
        with open(name,'rb') as f:
            obj = pickle.load(f)
            return obj
    

import sys
import itertools
import statistics
import shutil
import random
import math
from functools import reduce
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import gc
import os
import psutil
import inspect
import re
import pickle
import torch
import torch.nn.functional
import datetime
from torch import nn
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from torch.utils.data import DataLoader, Dataset, SequentialSampler, RandomSampler
from PIL import Image
from sklearn.preprocessing import MinMaxScaler
import time
from IPython.core.interactiveshell import InteractiveShell 
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from collections import Counter

InteractiveShell.ast_node_interactivity = "all"

now = datetime.datetime.now()
current_time = now.strftime("%Y-%m-%d %H:%M:%S")
print(current_time)

gc.collect()

available_memory = psutil.virtual_memory().available

available_memory_gb = available_memory / (1024 * 1024 * 1024)

print(f"available_memory: {available_memory_gb} GB")

used_memory = psutil.virtual_memory().used

used_memory_gb = used_memory / (1024 * 1024 * 1024)

print(f"used_memory: {used_memory_gb} GB")




if the_method == 'cell-blast':
    import os
    import time
    import pandas as pd
    import warnings
    warnings.filterwarnings("ignore")
    
    import Cell_BLAST as cb
    import numpy as np
    from numpy import genfromtxt as gft
    
    import scanpy as sc
    import psutil
    import anndata
    
    warnings.simplefilter("ignore")
    cb.config.RANDOM_SEED = 0
elif the_method == 'scvi-tools':
    import os
    import tempfile
    import warnings
    
    import anndata
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    # import scrublet as scr
    import scvi
    import torch
    import psutil
    
    scvi.settings.seed = 0
    print("Last run with scvi-tools version:", scvi.__version__)
    
    sc.set_figure_params(figsize=(4, 4), frameon=False)
    torch.set_float32_matmul_precision("high")
    save_dir = './'
    warnings.simplefilter(action="ignore", category=FutureWarning)
else:
    import torchvision.models
    from d2l import torch as d2l
    from torchvision import transforms

#!/usr/bin/env python
# coding: utf-8

# Some references:
# D2l. https://zh-v2.d2l.ai/
# Pandas. https://pandas.pydata.org/docs/index.html
# PyTorch. https://pytorch.org/
# Scikit-learn. https://scikit-learn.org/stable/

# import sys
# import itertools
# import statistics
# import shutil
# import random
# import math
# from functools import reduce
# import numpy as np
# import pandas as pd
# import scanpy as sc
# import scanpy.external as sce
# import seaborn as sns
# import gc
# import os
# import psutil
# import inspect
# import re
# import pickle
# import torch
# import torchvision.models
# import torch.nn.functional
# import datetime
# from torch import nn
# from d2l import torch as d2l
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
# from torch.utils.data import DataLoader, Dataset, SequentialSampler, RandomSampler
# from PIL import Image
# from sklearn.preprocessing import MinMaxScaler
# import time
# from IPython.core.interactiveshell import InteractiveShell 
# from torchvision import transforms  
# from sklearn.svm import LinearSVC
# from sklearn.calibration import CalibratedClassifierCV
# from collections import Counter

InteractiveShell.ast_node_interactivity = "all"

now = datetime.datetime.now()
current_time = now.strftime("%Y-%m-%d %H:%M:%S")
print(current_time)

gc.collect()

available_memory = psutil.virtual_memory().available

available_memory_gb = available_memory / (1024 * 1024 * 1024)

print(f"available_memory: {available_memory_gb} GB")

used_memory = psutil.virtual_memory().used

used_memory_gb = used_memory / (1024 * 1024 * 1024)

print(f"used_memory: {used_memory_gb} GB")



if the_method == 'pcma_mygo_rc' or the_method == 'pcma_mygo_rc2':
    class Dataset_pcma_mygo_rc_load_from_hardware(Dataset):  
        def __init__(self, root_dir, transform=None):  
            self.root_dir = root_dir  
            self.classes = sorted(os.listdir(root_dir))  
            self.class_to_idx = {cls_name: i for i, cls_name in enumerate(self.classes)}  
            self.samples = self._make_dataset()  
            self.transform = transforms.Compose([  
                transforms.Resize((224, 224)),  # 224x224  
                transforms.ToTensor()  # to tensor
            ])   
      
        def _make_dataset(self):  
            samples = []  
            for class_name in self.classes:  
                class_dir = os.path.join(self.root_dir, class_name)  
                if not os.path.isdir(class_dir):  
                    continue  
                for file_name in os.listdir(class_dir):  
                    file_path = os.path.join(class_dir, file_name)  
                    samples.append((file_path, self.class_to_idx[class_name]))  
            return samples  
      
        def __len__(self):  
            return len(self.samples)  
      
        def __getitem__(self, idx):  
            file_path, label = self.samples[idx]  
            image = torch.load(file_path)  
            return image, label  
            
    class CustomDataset2(Dataset):  
        def __init__(self, root_dir, transform=None):  
            self.root_dir = root_dir  
            self.classes = sorted(os.listdir(root_dir))  
            self.class_to_idx = {cls_name: i for i, cls_name in enumerate(self.classes)}  
            self.samples = self._make_dataset()  
            self.transform = transforms.Compose([  
                transforms.Resize((224, 224)),  # 224x224  
                transforms.ToTensor()  # to tensor
            ])   
      
        def _make_dataset(self):  
            samples = []  
            for class_name in self.classes:  
                class_dir = os.path.join(self.root_dir, class_name)  
                if not os.path.isdir(class_dir):  
                    continue  
                for file_name in os.listdir(class_dir):  
                    file_path = os.path.join(class_dir, file_name)  
                    samples.append((file_path, self.class_to_idx[class_name]))  
            return samples  
      
        def __len__(self):  
            return len(self.samples)  
      
        def __getitem__(self, idx):  
            file_path, label = self.samples[idx]  
            image = torch.load(file_path)  
            return image, label  

class Dataset_PCMA_MYGO_RC_2(Dataset):  
    def __init__(self,data,label):
        self.Data = torch.from_numpy(np.float32(np.array(data)))
        new_tensor = torch.zeros(self.Data.size(0), 3*224*224)
        new_tensor[:, :self.Data.size(1)] = self.Data
        self.Data = new_tensor
        self.Data = self.Data.view(self.Data.size(0), 3, 224, 224)
        self.Label = torch.from_numpy(np.int64(np.array(label)).reshape(1,-1)[0])
    def __getitem__(self, index):
        x = self.Data[index]
        y = self.Label[index]
        return x, y  
    def __len__(self):
        return len(self.Data)

if the_method == 'pcma_mygo_rc2':
    class PCMA_MYGO_RC2(nn.Module):
        def __init__(self, dropout_ratio1=0.1, dropout_ratio2=0.1, dropout_ratio3=0.1, num_types=None):
            super(PCMA_MYGO_RC2, self).__init__()
            self.dropout_layer1 = nn.Dropout(dropout_ratio1)
            self.resnet18 = torchvision.models.resnet18(weights=None)   
            self.PCMA_MYGO_FC = nn.Sequential(
                            nn.ReLU(),
                            nn.Dropout(dropout_ratio2),
                            nn.Linear(1000, 512),
                            nn.ReLU(),
                            nn.Dropout(dropout_ratio3),
                            nn.Linear(512, num_types),
            )
            
        def forward(self, x):
            x = self.dropout_layer1(x)
            x = self.resnet18(x)
            x = self.PCMA_MYGO_FC(x)
            return x

class Dataset_PCMA_MYGO_RC(Dataset):  
    def __init__(self,data,label):
        self.Data = torch.from_numpy(np.float32(np.array(data)))
        self.Data = torch.reshape(self.Data, (data.shape[0], 224, 224))
        self.Data = torch.unsqueeze(self.Data, 1)
        self.Data = torch.cat((self.Data, self.Data, self.Data), dim=1)
        self.Label = torch.from_numpy(np.int64(np.array(label)).reshape(1,-1)[0])
    def __getitem__(self, index):
        x = self.Data[index]
        y = self.Label[index]
        return x, y  
    def __len__(self):
        return len(self.Data)

if the_method == 'pcma_mygo_rc':
    class PCMA_MYGO_RC(nn.Module):
        def __init__(self, dropout_ratio1=0.3, dropout_ratio2=0.1, dropout_ratio3=0.1, num_types=None):
            super(PCMA_MYGO_RC, self).__init__()
            self.resnet18 = torchvision.models.resnet18(pretrained=False)   
            self.PCMA_MYGO_FC = nn.Sequential(
                            nn.ReLU(),
                            nn.Dropout(dropout_ratio1),
                            nn.Linear(1000, num_types))
            
        def forward(self, x):
            x = self.resnet18(x)
            x = self.PCMA_MYGO_FC(x)
            return x

class Dataset_PCMA_MYGO_FC(Dataset):  
    def __init__(self,data,label):
        self.Data = torch.from_numpy(np.float32(np.array(data)))
        self.Label = torch.from_numpy(np.int64(np.array(label)).reshape(1,-1)[0])
    def __getitem__(self, index):
        x = self.Data[index]
        y = self.Label[index]
        return x, y  
    def __len__(self):
        return len(self.Data)

class PCMA_MYGO_FC(nn.Module):
    def __init__(self, feature_len=None, dropout_ratio1=0.1, dropout_ratio2=0.1, dropout_ratio3=0.1, num_classes=None):
        super().__init__() 
        self.PCMA_MYGO_FC = nn.Sequential(nn.Flatten(),
                               nn.Dropout(dropout_ratio1),
                               nn.Linear( feature_len, feature_len//4 ),
                               nn.ReLU(),
                               nn.Dropout(dropout_ratio2),
                               nn.Linear( feature_len//4, feature_len//16 ),
                               nn.ReLU(),
                               nn.Dropout(dropout_ratio3),
                               nn.Linear(feature_len//16, num_classes))
        
    def forward(self, x):
        x = self.PCMA_MYGO_FC(x)
        return x

class Dataset_transformer(Dataset):  
    def __init__(self,data,label):
        self.Data = torch.from_numpy(np.float32(np.array(data)))
        self.Data = torch.reshape(self.Data, (data.shape[0], 2, 1024))
        self.Label = torch.from_numpy(np.int64(np.array(label)).reshape(1,-1)[0])
    def __getitem__(self, index):
        x = self.Data[index]
        y = self.Label[index]
        return x, y  
    def __len__(self):
        return len(self.Data)

class PositionalEncoding(nn.Module):
    def __init__(self, dim1, dropout=0.1, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)
        pe = torch.zeros(max_len, dim1)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, dim1, 2).float() * (-math.log(10000.0) / dim1))
        
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[:x.size(0), :]
        return self.dropout(x)

class PCMA_MYGO_SA(nn.Module):
    def __init__(self, input_dim=1024, num_head=2, dim1=1024, num_classes=2, dropout=0.2):
        super().__init__()
        self.positionalEncoding = PositionalEncoding(dim1=dim1, dropout=dropout)
        self.encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim1, dim_feedforward=dim1, nhead=num_head, dropout=dropout*0.75,
        )
        
        self.pred_layer = nn.Sequential(
            nn.Dropout(dropout*0.5),
            nn.Linear(dim1, dim1),
            nn.ReLU(),
            nn.Dropout(dropout*0.25),
            nn.Linear(dim1, num_classes)
        )

    def forward(self, x):
        out = x.permute(1, 0, 2)
        out = self.positionalEncoding(out)
        out = self.encoder_layer(out)
        out = out.transpose(0, 1)
        out = out.mean(dim=1)
        out = self.pred_layer(out)
        return out

def same_seeds(seed=None):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

def print_memory_used():
    process = psutil.Process()  
    memory_usage = process.memory_info().rss / 1024 / 1024 /1024 
    print(f"memory used: {memory_usage} GB")  

def print_anndata(input_adata=None):
    new_anndata = input_adata
    print(new_anndata)
    print(new_anndata.obs_names)
    print(new_anndata.var_names)
    print(new_anndata.X.max())
    print(new_anndata.obs)
    print(new_anndata.var)

def preprocessing(input_adata=None,if_filter_genes=True,if_norm=True,if_hvg=False,n_top_hvgs=None):
    if isinstance(input_adata, list) and len(input_adata)==2:
        new_anndata1 = input_adata[0]
        new_anndata1.obs['ref_test'] = 'ref'
        new_anndata2 = input_adata[1]
        new_anndata2.obs['ref_test'] = 'test'
        new_anndata = new_anndata1.concatenate(new_anndata2)
    else:
        new_anndata = input_adata
    if if_filter_genes==True:
        sc.pp.filter_genes(new_anndata, min_cells=1)
        print('Genes filtered.')
    if if_norm==True:
        if new_anndata.X.max() > 50:
            sc.pp.normalize_total(new_anndata, target_sum=1e4)
            sc.pp.log1p(new_anndata)
            print('Data normalized.')
    new_anndata.var_names_make_unique()
    new_anndata.obs_names_make_unique()
    print('Names are changed to be unique.')
    if if_hvg==True:
        if n_top_hvgs is not None:
            sc.pp.highly_variable_genes(new_anndata,n_top_genes=n_top_hvgs)
        else:
            sc.pp.highly_variable_genes(new_anndata)
        sc.pl.highly_variable_genes(new_anndata)
        new_anndata = new_anndata[:, new_anndata.var.highly_variable]
        print('Hvgs selected.')
    return new_anndata

def split_anndata_by_celltype_to_train_valid_test(new_anndata=None, celltype_label='celltype', split_ratio=[0.8,0.1,0.1], random_seed=42):  
    # split adata to three part by given ratio
    np.random.seed(random_seed)
    train_anndata = None
    valid_anndata = None
    test_anndata = None

    if len(split_ratio)==2:
        for celltype in new_anndata.obs[celltype_label].unique(): 
            celltype_indices = new_anndata.obs[celltype_label] == celltype  
            celltype_data = new_anndata[celltype_indices]  
            
            total_cells = len(celltype_data)  
            x, y = split_ratio  
            x_count = int(total_cells * x / (x + y))   
            y_count = total_cells - x_count  
            
            indices = np.random.permutation(total_cells)  
            train_indices = indices[:x_count]  
            valid_indices = indices[x_count:]  
            
            train_data = celltype_data[train_indices]  
            valid_data = celltype_data[valid_indices]  
            
            if train_anndata is None:
                train_anndata = train_data
            else:
                train_anndata = train_anndata.concatenate(train_data, index_unique=None)  
            if valid_anndata is None:
                valid_anndata = valid_data
            else:
                valid_anndata = valid_anndata.concatenate(valid_data, index_unique=None)   
        return train_anndata, valid_anndata
    
    for celltype in new_anndata.obs[celltype_label].unique(): 
        celltype_indices = new_anndata.obs[celltype_label] == celltype  
        celltype_data = new_anndata[celltype_indices]  
        
        total_cells = len(celltype_data)  
        x, y, z = split_ratio  
        x_count = int(total_cells * x / (x + y + z))  
        y_count = int(total_cells * y / (x + y + z))  
        z_count = total_cells - x_count - y_count  
        
        indices = np.random.permutation(total_cells)  
        train_indices = indices[:x_count]  
        valid_indices = indices[x_count:x_count + y_count]  
        test_indices = indices[x_count + y_count:]  
        
        train_data = celltype_data[train_indices]  
        valid_data = celltype_data[valid_indices]  
        test_data = celltype_data[test_indices]  
        
        if train_anndata is None:
            train_anndata = train_data
        else:
            train_anndata = train_anndata.concatenate(train_data, index_unique=None)  
        if valid_anndata is None:
            valid_anndata = valid_data
        else:
            valid_anndata = valid_anndata.concatenate(valid_data, index_unique=None)   
        if test_anndata is None:
            test_anndata = test_data
        else:
            test_anndata = test_anndata.concatenate(test_data, index_unique=None)   
        
    return train_anndata, valid_anndata, test_anndata

def generate_datas_and_labels(input_adatas=None, the_length=None, if_balance_train_data=None, random_sample_num=None,if_need_test_mapping=None):
    new_anndata_train = input_adatas[0]
    new_anndata_valid = input_adatas[1]
    new_anndata_test = input_adatas[2]
    pcma = PCmaster_anno_0()
    train_data,mapping_1,mapping_2 = pcma.auto_annotation_with_deep_learning_transfer_adata_to_df_for_all_0(
        adata_input = new_anndata_train, the_length=the_length)
    valid_data,_,_ = pcma.auto_annotation_with_deep_learning_transfer_adata_to_df_for_all_0(
        adata_input = new_anndata_valid, the_length=the_length)
    test_data,mapping_1_test,mapping_2_test = pcma.auto_annotation_with_deep_learning_transfer_adata_to_df_for_all_0(
        adata_input = new_anndata_test, the_length=the_length)
    
    # displays the number of categories and the number of each category.
    category_counts = train_data['celltype'].value_counts()
    print("the number of categories: ", len(category_counts))
    print("the number of each category: ")
    print(category_counts)
    print(mapping_1)
    print(mapping_2)
    print('mapping_test_dataset')
    print(mapping_1_test)
    print(mapping_2_test)
    
    if if_balance_train_data is True:
        train_data_balanced = pd.DataFrame()
        for category in category_counts.index:
            category_data = train_data[train_data['celltype'] == category]
            print(category)
            print(category_data.shape[0])
            # print(type(category_data.shape[0]))
            category_data = category_data.sample(n=random_sample_num, random_state=1, replace=True)
            train_data_balanced = pd.concat([train_data_balanced, category_data])
        train_data = train_data_balanced
    
    print_memory_used() 
    
    train_label = train_data.iloc[:, -1:]
    train_data = train_data.iloc[:, :-2]
    valid_label = valid_data.iloc[:, -1:]
    valid_data = valid_data.iloc[:, :-2]
    test_label = test_data.iloc[:, -1:]
    test_data = test_data.iloc[:, :-2]
    
    if if_need_test_mapping is True:
        return train_data,train_label,valid_data,valid_label,test_data,test_label,mapping_1,mapping_2,category_counts,mapping_1_test,mapping_2_test
    else:
        return train_data,train_label,valid_data,valid_label,test_data,test_label,mapping_1,mapping_2,category_counts

# GPU
def evaluate_accuracy_gpu(net, data_iter, device=None): 
    if isinstance(net, nn.Module):
        net.eval()
        if not device:
            device = next(iter(net.parameters())).device
    metric = d2l.Accumulator(2)
    with torch.no_grad():
        for X, y in data_iter:
            if isinstance(X, list):
                X = [x.to(device) for x in X]
            else:
                X = X.to(device)
            y = y.to(device)
            metric.add(d2l.accuracy(net(X), y), y.numel())
    return metric[0] / metric[1]

# GPU
def train_with_GPU(net=None, train_iter=None, valid_iter=None, num_epochs=None, lr=None, device=None, optimizer='Adam', momentum=None):
    def init_weights(m):
        if type(m) == nn.Linear or type(m) == nn.Conv2d:
            nn.init.xavier_uniform_(m.weight)
    net.apply(init_weights)
    print('training on', device)
    net.to(device)
    if optimizer == 'SGD':
        if momentum is not None:
            optimizer = torch.optim.SGD(net.parameters(), lr=lr, momentum=momentum)
        else:
            optimizer = torch.optim.SGD(net.parameters(), lr=lr)
    if optimizer == 'Adam':
        optimizer = torch.optim.Adam(net.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=150, eta_min=0.00001)
    loss = nn.CrossEntropyLoss()
    animator = d2l.Animator(xlabel='epoch', xlim=[1, num_epochs],ylim=[0,1],
                            legend=['train loss', 'train acc', 'valid acc'])
    timer, num_batches = d2l.Timer(), len(train_iter)
    best_model_weights = None
    best_valid_acc = 0.0
    best_epoch = 0
    for epoch in range(num_epochs):
        metric = d2l.Accumulator(3)
        net.train()
        for i, (X, y) in enumerate(train_iter):
            timer.start()
            optimizer.zero_grad()
            X, y = X.to(device), y.to(device)
            y_hat = net(X)
            l = loss(y_hat, y)
            l.backward()
            optimizer.step()
            scheduler.step()
            with torch.no_grad():
                metric.add(l * X.shape[0], d2l.accuracy(y_hat, y), X.shape[0])
            timer.stop()
            train_l = metric[0] / metric[2]
            train_acc = metric[1] / metric[2]
            if (i + 1) % (num_batches // 5) == 0 or i == num_batches - 1:
                animator.add(epoch + (i + 1) / num_batches,
                             (train_l, train_acc, None))
        valid_acc = evaluate_accuracy_gpu(net, valid_iter)
        if valid_acc > best_valid_acc:
            best_valid_acc = valid_acc
            best_epoch = epoch
            best_model_weights = net.state_dict()
        animator.add(epoch + 1, (None, None, valid_acc))
    net.load_state_dict(best_model_weights)
    print(f'At the epoch {best_epoch}, the model got the best valid accuracy {best_valid_acc}.')
    print(f'loss {train_l:.3f}, train acc {train_acc:.3f}, '
          f'valid acc {valid_acc:.3f}')
    print(f'{metric[2] * num_epochs / timer.sum():.1f} examples/sec '
          f'on {str(device)}')

def merge_dicts(dict1, dict2):    
    common_keys = set(dict1.keys()) & set(dict2.keys())    
    dict1_only_keys = [key for key in dict1.keys() if key not in common_keys]  
    dict2_only_keys = [key for key in dict2.keys() if key not in common_keys]   
    merged_dict = {key: i for i, key in enumerate(common_keys)}   
    for i, key in enumerate(dict1_only_keys, start=len(merged_dict)):  
        merged_dict[key] = i   
    for i, key in enumerate(dict2_only_keys, start=len(merged_dict)):  
        merged_dict[key] = i  
    return merged_dict, set(common_keys), set(dict1_only_keys), set(dict2_only_keys) 

def back_to_unmapped(obj=None,mapping_1=None,mapping_2=None):
    if not isinstance(obj, list):
        list2 = obj.tolist()
    else:
        list2 = obj
    reverse_mapping_2 = {v: k for k, v in mapping_2.items()}
    list2_0_mapped = [reverse_mapping_2[value] for value in list2]
    list2_0 = [next(key for key, value in mapping_1.items() if value == mapped_value) for mapped_value in list2_0_mapped]
    return list2_0

def generate_mapping_for_list(input_list=None):
    df = pd.DataFrame(input_list, columns=['celltype'])
    mapping_1 = {category: i for i, category in enumerate(df['celltype'].unique())}
    df['celltype'] = df['celltype'].map(mapping_1)
    mapping_2 = {category: i for i, category in enumerate(df['celltype'].unique())}
    df['celltype'] = df['celltype'].map(mapping_2)
    celltype_list = df['celltype'].tolist()
    
    return celltype_list,mapping_1,mapping_2

def print_list(input_list=None):
    for i in range(0,min(10,len(input_list))):
        print(input_list[i])

def change_list_list_to_list(input_list=None):
    if isinstance(input_list[0],list):
        print('change list_list to list')
        tmp_list = []
        for i in input_list:
            tmp_list.append(i[0])
        return tmp_list
    else:
        return input_list

if the_method == 'cell-blast':
    def processing_for_cell_blast(test_label=None,predicted=None,
                                    mapping_1=None,mapping_2=None,
                                    mapping_1_test=None,mapping_2_test=None,):
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        if mapping_1 is not None and mapping_2 is not None and \
        mapping_1_test is not None and mapping_2_test is not None and mapping_1 != mapping_1_test:
            test_label = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
            predicted = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
            test_label_unmapped = test_label
            predicted_unmapped = predicted
            merged_mapping,set_common_keys,set_dict1_only_keys,set_dict2_only_keys = merge_dicts(dict1=mapping_1,dict2=mapping_1_test)
            # print(merged_mapping)
            tmp_pre = []
            tmp_test = []
            test_total_count = 0
            true_common = 0
            predicted_common_and_true_common = 0
            for i,j in enumerate(predicted):
                test_total_count = test_total_count+1
                if test_label[i] in set_common_keys:
                    true_common = true_common+1
                    if j in set_common_keys:
                        predicted_common_and_true_common = predicted_common_and_true_common+1
                        tmp_pre.append(j)
                        tmp_test.append(test_label[i])
            test_label = tmp_test
            predicted = tmp_pre
            test_label = [mapping_1_test[key] for key in test_label]
            test_label = [mapping_2_test[key] for key in test_label]
            predicted = [mapping_1[key] for key in predicted]
            predicted = [mapping_2[key] for key in predicted]
        
        return test_label,predicted

def calculate_f1_for_each_class(true_list, pred_list):  
    f1_scores = {}  
    classes = set(true_list + pred_list)  # Get all unique classes  
    for class_ in classes:  
        true_class = [1 if x == class_ else 0 for x in true_list]  # Convert true labels to binary array  
        pred_class = [1 if x == class_ else 0 for x in pred_list]  # Convert predicted labels to binary array  
        f1_scores[class_] = f1_score(true_class, pred_class)  # Calculate F1 score and store  
    return f1_scores  

def get_accuracy_precision_recall_f1_score(test_label=None,predicted=None,
                                           mapping_1=None,mapping_2=None,
                                           mapping_1_test=None,mapping_2_test=None,
                                           unknown_model=True,
                                          ):
    model_result_list = []
    unknown_label = 'Unknown'
    unknown_accuracy = -1.0
    unknown_precision = -1.0
    unknown_recall = -1.0
    unknown_f1_score = -1.0
    true_common_in_test_total_rate = -1.0
    predicted_common_in_true_common_rate = -1.0
    if unknown_model is True:
        mapping_1[unknown_label]=-100
        mapping_2[-100]=-100
    test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
    predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
    if mapping_1 is not None and mapping_2 is not None and \
    mapping_1_test is not None and mapping_2_test is not None and mapping_1 != mapping_1_test:
        if unknown_model is True:
            mapping_1[unknown_label]=-100
            mapping_2[-100]=-100
            mapping_1_test[unknown_label]=-100
            mapping_2_test[-100]=-100
        test_label = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        test_label_unmapped = test_label
        predicted_unmapped = predicted
        merged_mapping,set_common_keys,set_dict1_only_keys,set_dict2_only_keys = merge_dicts(dict1=mapping_1,dict2=mapping_1_test)
        print('merged_mapping')
        print(merged_mapping)
        print('set_common_keys')
        print(set_common_keys)
        print('set_dict1_only_keys')
        print(set_dict1_only_keys)
        print('set_dict2_only_keys')
        print(set_dict2_only_keys)
        tmp_pre = []
        tmp_test = []
        pre_unknown = []
        test_unknown = []
        for i,j in enumerate(predicted):
            if j==unknown_label:
                pre_unknown.append(j)
                if test_label[i] in set_common_keys:
                    test_unknown.append('Not_'+unknown_label)
                elif test_label[i] in set_dict2_only_keys:
                    test_unknown.append(unknown_label)
            elif j!=unknown_label and test_label[i] in set_dict2_only_keys:
                pre_unknown.append('Not_'+unknown_label)
                test_unknown.append(unknown_label)
            elif j!=unknown_label and test_label[i] in set_common_keys:
                pre_unknown.append('Not_'+unknown_label)
                test_unknown.append('Not_'+unknown_label)
        test_total_count = 0
        true_common = 0
        predicted_common_and_true_common = 0
        for i,j in enumerate(predicted):
            test_total_count = test_total_count+1
            if test_label[i] in set_common_keys: # when in set_common_keys, normally calculate acc f1, is that appropriate?
                true_common = true_common+1
                tmp_pre.append(j)
                tmp_test.append(test_label[i])
                if j != unknown_label and j in set_common_keys:
                    predicted_common_and_true_common = predicted_common_and_true_common+1
                    # tmp_pre.append(j)
                    # tmp_test.append(test_label[i])
        print('test_total_count')
        print(test_total_count)
        print('true_common')
        print(true_common)
        print('predicted_common_and_true_common')
        print(predicted_common_and_true_common)
        # true_common_in_test_total_rate = true_common / test_total_count
        # predicted_common_in_true_common_rate = predicted_common_and_true_common / true_common
        # print("true_common_in_test_total_rate:", true_common_in_test_total_rate)  
        # print("predicted_common_in_true_common_rate:", predicted_common_in_true_common_rate)
        test_label = tmp_test
        predicted = tmp_pre
        test_label = [merged_mapping[key] for key in test_label]
        predicted = [merged_mapping[key] for key in predicted]
        unknown_mapping = {unknown_label:0,('Not_'+unknown_label):1}
        pre_unknown = [unknown_mapping[key] for key in pre_unknown]
        test_unknown = [unknown_mapping[key] for key in test_unknown]
        unknown_accuracy = accuracy_score(test_unknown, pre_unknown)
        unknown_precision = precision_score(test_unknown, pre_unknown, average='macro', zero_division=1)   
        unknown_recall = recall_score(test_unknown, pre_unknown, average='macro', zero_division=1)  
        unknown_f1_score = f1_score(test_unknown, pre_unknown, average='macro', zero_division=1) 
        if len(merged_mapping.keys()) == len(set_common_keys):
            print('mapping_1 is just different from mapping_1_test by order, the keys are the same.')
            unknown_accuracy = -1.0
            unknown_precision = -1.0   
            unknown_recall = -1.0  
            unknown_f1_score = -1.0
        print("Unknown Accuracy:", unknown_accuracy)  
        print("Unknown Macro Precision:", unknown_precision)  
        print("Unknown Macro Recall:", unknown_recall)  
        print("Unknown Macro F1 Score:", unknown_f1_score)
    model_result_list.append('test_total_count: '+str(test_total_count)+' []')
    model_result_list.append('true_common: '+str(true_common)+' []')
    model_result_list.append('predicted_common_and_true_common: '+str(predicted_common_and_true_common)+' []')
    model_result_list.append('unknown_accuracy: '+str(unknown_accuracy)+' []')
    model_result_list.append('unknown_precision: '+str(unknown_precision)+' []')
    model_result_list.append('unknown_recall: '+str(unknown_recall)+' []')
    model_result_list.append('unknown_f1_score: '+str(unknown_f1_score)+' []')
    accuracy = accuracy_score(test_label, predicted)  
    macro_precision = precision_score(test_label, predicted, average='macro', zero_division=1)   
    macro_recall = recall_score(test_label, predicted, average='macro', zero_division=1)  
    macro_f1_score = f1_score(test_label, predicted, average='macro', zero_division=1)  
    model_result_list.append('accuracy: '+str(accuracy)+' []')
    model_result_list.append('macro_precision: '+str(macro_precision)+' []')
    model_result_list.append('macro_recall: '+str(macro_recall)+' []')
    model_result_list.append('macro_f1_score: '+str(macro_f1_score)+' []')
    dict_label_count = {}
    for label in test_label:
        if label not in dict_label_count:
            dict_label_count[label]=1
        else:
            dict_label_count[label] = dict_label_count[label]+1
    print('test_label_count(celltype-count):')
    print(dict_label_count)
    sorted_labels = sorted(dict_label_count, key=dict_label_count.get)  
    top3_labels = sorted_labels[:3]  
    top1_labels = sorted_labels[:1]  
    print('three cell types with the least number:')
    print(top3_labels)
    all_f1_scores = f1_score(test_label, predicted, average=None)  
    all_f1_scores_dict = calculate_f1_for_each_class(test_label, predicted)
    median_f1_score = np.median(all_f1_scores)  
    model_result_list.append('median_f1_score: '+str(median_f1_score)+' []')
    if len(top3_labels) == 3:
        f1_score_least_1 = all_f1_scores_dict[top3_labels[0]]
        f1_score_least_2 = all_f1_scores_dict[top3_labels[1]]
        f1_score_least_3 = all_f1_scores_dict[top3_labels[2]]
        print('f1 score of the cell type with the least number:')
        print(f1_score_least_1)
        print('f1 score of the cell type with the second least number:')
        print(f1_score_least_2)
        print('f1 score of the cell type with the third least number:')
        print(f1_score_least_3)
        average_f1_score_three_least = (f1_score_least_1+f1_score_least_2+f1_score_least_3)/3.0
    else:
        print('The num of common cell types between the ref dataset and the unknown dataset is less than three!')
        average_f1_score_three_least = -1
    model_result_list.append('average_f1_score_three_least: '+str(average_f1_score_three_least)+' []')
    if len(top1_labels) == 1:
        f1_score_least_1 = all_f1_scores_dict[top3_labels[0]]
        print('f1 score of the cell type with the least number:')
        print(f1_score_least_1)
    else:
        print('The num of common cell types between the ref dataset and the unknown dataset is less than one!')
        f1_score_least_1 = -1
    model_result_list.append('f1_score_least_1: '+str(f1_score_least_1)+' []')
    
    print("Accuracy:", accuracy)  
    print("Macro Precision:", macro_precision)  
    print("Macro Recall:", macro_recall)  
    print("Macro F1 Score:", macro_f1_score)
    print("All F1 Scores:", all_f1_scores_dict)
    print("Median F1 Score:", median_f1_score)
    print("Average f1 score of three least cell types:",average_f1_score_three_least)
    
    return model_result_list,test_label_unmapped,predicted_unmapped

if the_method == 'cell-blast':
    def train_Cell_BLAST_with_unknown(
        new_anndata=None,new_anndata2=None, # new_anndata = new_anndata_train; new_anndata2 = new_anndata_test
        epochs=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        start_time=time.time()
        models = []
        if epochs is None:
            for i in range(4):
                models.append(cb.directi.fit_DIRECTi(
                    new_anndata, genes=new_anndata.var_names,
                    latent_dim=10, cat_dim=20, random_seed=i,
                ))
        else:
            for i in range(4):
                models.append(cb.directi.fit_DIRECTi(
                    new_anndata, genes=new_anndata.var_names,
                    latent_dim=10, cat_dim=20, random_seed=i, epoch=epochs,
                ))
        print("Time elapsed: %.1fs" % (time.time() - start_time))
        
        blast = cb.blast.BLAST(models, new_anndata)
        new_anndata2_hits = blast.query(new_anndata2)
        print("Time per query: %.1fms" % (
            (time.time() - start_time) * 1000 / new_anndata2.shape[0]
        ))
        new_anndata2_hits = new_anndata2_hits.reconcile_models().filter(by="pval", cutoff=0.05)
        new_anndata2_predictions = new_anndata2_hits.annotate("celltype")
        truelist = new_anndata2.obs['celltype'].tolist()
        predlist = new_anndata2_predictions['celltype'].tolist()
        
        # print('truelist')
        # print_list(input_list=truelist)
        # print('predlist')
        # print_list(input_list=predlist)
        model = blast
        
        gc.collect()
        print_memory_used()
        
        test_label,mapping_1_test,mapping_2_test = generate_mapping_for_list(input_list=truelist)
        predicted,mapping_1,mapping_2 = generate_mapping_for_list(input_list=predlist)
        # print('test_label')
        # print_list(input_list=test_label)
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        test_label_unmapped = truelist
        predicted_unmapped = predlist
        
        return model,test_label_unmapped,predicted_unmapped
        
elif the_method == 'scvi-tools':
    def train_scvi_tools_no_unknown(
        new_anndata_train=None,new_anndata_test=None, 
        epochs=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):  
        start_time=time.time()
        scvi.model.SCVI.setup_anndata(new_anndata_train)
            
        scvi_ref = scvi.model.SCVI(
            new_anndata_train,
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
            dropout_rate=0.2,
            n_layers=2,
        )
        scvi_ref.train()
        # scvi_ref.train(max_epochs=5)
        SCVI_LATENT_KEY = "X_scVI"
        new_anndata_train.obsm[SCVI_LATENT_KEY] = scvi_ref.get_latent_representation()
        
        SCANVI_LABELS_KEY = "labels_scanvi"
        new_anndata_train.obs[SCANVI_LABELS_KEY] = new_anndata_train.obs[celltype_label].values
        scanvi_ref = scvi.model.SCANVI.from_scvi_model(
            scvi_ref,
            unlabeled_category="Unknown",
            labels_key=SCANVI_LABELS_KEY,
        )
        scanvi_ref.train(max_epochs=epochs[0], n_samples_per_label=100)
        SCANVI_LATENT_KEY = "X_scANVI"
        new_anndata_train.obsm[SCANVI_LATENT_KEY] = scanvi_ref.get_latent_representation()
        sc.pp.neighbors(new_anndata_train, use_rep=SCANVI_LATENT_KEY)
        sc.tl.leiden(new_anndata_train)
        sc.tl.umap(new_anndata_train)
        print("Train finished. Time elapsed: %.1fs" % (time.time() - start_time))
        
        model = scanvi_ref
        
        gc.collect()
        print_memory_used()
        new_anndata_test_copy = new_anndata_test.copy()
        scvi.model.SCANVI.prepare_query_anndata(new_anndata_test, scanvi_ref)
        scanvi_query = scvi.model.SCANVI.load_query_data(new_anndata_test_copy, scanvi_ref)
        
        scanvi_query.train(
            max_epochs=epochs[1],
            plan_kwargs={"weight_decay": 0.0},
            check_val_every_n_epoch=10,
        )
        
        SCANVI_PREDICTIONS_KEY = "predictions_scanvi"
        
        new_anndata_test.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
        new_anndata_test.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()
        truelist = new_anndata_test.obs['celltype'].tolist()
        predlist = new_anndata_test.obs[SCANVI_PREDICTIONS_KEY].tolist()
        
        test_label,mapping_1_test,mapping_2_test = generate_mapping_for_list(input_list=truelist)
        predicted,mapping_1,mapping_2 = generate_mapping_for_list(input_list=predlist)
        # print('test_label')
        # print_list(input_list=test_label)
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        test_label_unmapped = truelist
        predicted_unmapped = predlist
        
        return model,test_label_unmapped,predicted_unmapped
else:
    def train_PCMA_MYGO_FC_with_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_batch_size=None,num_classes=None,lr=None,dropout_ratio1=None,dropout_ratio2=None,dropout_ratio3=None,feature_len=None,
        epochs=None,gpu_code=None,the_num_workers=None,optimizer=None, momentum=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        traindata = Dataset_PCMA_MYGO_FC(train_data,train_label) 
        validdata = Dataset_PCMA_MYGO_FC(valid_data,valid_label)
        
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        
        del train_data,train_label,valid_data,valid_label,traindata,validdata
        gc.collect()
        print_memory_used()
        print(f'num_classes = {num_classes}')
        
        pcma_mygo_fc = PCMA_MYGO_FC(
            feature_len=feature_len,
            dropout_ratio1=dropout_ratio1, dropout_ratio2=dropout_ratio2, dropout_ratio3=dropout_ratio3, 
            num_classes=num_classes,
        )
        
        torch.cuda.init()
        # net=None, train_iter=None, valid_iter=None, num_epochs=None, lr=None, device=None, optimizer='Adam', momentum=None
        train_with_GPU(net=pcma_mygo_fc, train_iter=train_iter, valid_iter=valid_iter, 
                       num_epochs=epochs, lr=lr, device=d2l.try_all_gpus()[gpu_code], optimizer=optimizer, momentum=momentum)
        model = pcma_mygo_fc
        del train_sampler,train_iter,valid_sampler,valid_iter
        gc.collect()
        print_memory_used()
        
        import math
        
        batch_size = 32 
        total_size = test_data.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = model
        net_try = net_try.to('cpu')
        net_try.eval()
        
        with torch.no_grad():
            predicted_results = []
            for i in range(total_batches):
                start = i * batch_size
                end = min(start + batch_size, total_size)
                batch_data = test_data[start:end]
                batch_label = test_label[start:end]
                batch_data = torch.from_numpy(np.float32(np.array(batch_data)))
                batch_label = torch.from_numpy(np.int64(np.array(batch_label)).reshape(1,-1)[0])
                outputs = net_try(batch_data)
                outputs = torch.nn.functional.softmax(outputs, dim=1)
                max_values, predicted = torch.max(outputs, dim=1)  
                max_values_list = max_values.numpy().tolist()  
                max_values_list = [round(num, 4) for num in max_values_list]  
                the_median = np.median(max_values_list)
                _, indices = torch.topk(outputs, 2, dim=1)  
                second_largest_values = torch.gather(outputs, 1, indices[:, 1].unsqueeze(1))   
                second_largest_values_list = second_largest_values.numpy().flatten().tolist()  
                second_largest_values_list = [round(num, 4) for num in second_largest_values_list]  
                second_median = np.median(second_largest_values_list)   
                predicted_list = predicted.numpy().flatten().tolist()
                predicted_list_with_unknown = []
                for i,j in enumerate(predicted_list):
                    if max_values_list[i] > 3*second_largest_values_list[i]:
                        predicted_list_with_unknown.append(j)
                    else:
                        predicted_list_with_unknown.append(-100)
                predicted_results.append(np.array(predicted_list_with_unknown))
                # break
                
        predicted = np.concatenate(predicted_results)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        test_label = test_label.numpy()
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped
        
    def train_PCMA_MYGO_FC_no_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_batch_size=None,num_classes=None,lr=None,dropout_ratio1=None,dropout_ratio2=None,dropout_ratio3=None,feature_len=None,
        epochs=None,gpu_code=None,the_num_workers=None,optimizer=None, momentum=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        traindata = Dataset_PCMA_MYGO_FC(train_data,train_label) 
        validdata = Dataset_PCMA_MYGO_FC(valid_data,valid_label)
        
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        
        del train_data,train_label,valid_data,valid_label,traindata,validdata
        gc.collect()
        print_memory_used()
        print(f'num_classes = {num_classes}')
        
        pcma_mygo_fc = PCMA_MYGO_FC(
            feature_len=feature_len,
            dropout_ratio1=dropout_ratio1, dropout_ratio2=dropout_ratio2, dropout_ratio3=dropout_ratio3, 
            num_classes=num_classes,
        )
        
        torch.cuda.init()
        # net=None, train_iter=None, valid_iter=None, num_epochs=None, lr=None, device=None, optimizer='Adam', momentum=None
        train_with_GPU(net=pcma_mygo_fc, train_iter=train_iter, valid_iter=valid_iter, 
                       num_epochs=epochs, lr=lr, device=d2l.try_all_gpus()[gpu_code], optimizer=optimizer, momentum=momentum)
        model = pcma_mygo_fc
        del train_sampler,train_iter,valid_sampler,valid_iter
        gc.collect()
        print_memory_used()
        
        import math
        
        batch_size = 32 
        total_size = test_data.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = model
        net_try = net_try.to('cpu')
        net_try.eval()
        
        with torch.no_grad():
            predicted_results = []
            for i in range(total_batches):
                start = i * batch_size
                end = min(start + batch_size, total_size)
                batch_data = test_data[start:end]
                batch_label = test_label[start:end]
                batch_data = torch.from_numpy(np.float32(np.array(batch_data)))
                batch_label = torch.from_numpy(np.int64(np.array(batch_label)).reshape(1,-1)[0])
                outputs = net_try(batch_data)
                _, predicted = torch.max(outputs, 1)
                predicted_results.append(predicted.numpy())
                # break
                
        predicted = np.concatenate(predicted_results)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        test_label = test_label.numpy()
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped

    def train_PCMA_MYGO_HP_with_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_threshold=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        if the_threshold is None:
            the_threshold = 0.7
        Classifier = LinearSVC(C=1.2)
        clf = CalibratedClassifierCV(Classifier)
        start_time = time.time()
        clf.fit(train_data, train_label)
        
        del train_data,train_label,valid_data,valid_label
        gc.collect()
        print_memory_used()
        # print(f'num_classes = {num_classes}')
        
        model = clf
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"time cost：{execution_time} s")
        start_time = time.time()
        
        predicted = clf.predict(test_data)
        prob = np.max(clf.predict_proba(test_data), axis = 1)
        unlabeled = np.where(prob < the_threshold)
        predicted[unlabeled] = -100
        
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"time cost：{execution_time} s")
        
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        test_label = test_label.values.tolist()
        predicted = predicted.tolist()
        test_label = change_list_list_to_list(input_list=test_label)
        predicted = change_list_list_to_list(input_list=predicted)
        print(type(test_label))
        print(type(predicted))
        print_list(input_list=test_label)
        print_list(input_list=predicted)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped

    def train_PCMA_MYGO_HP_no_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_threshold=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        # if the_threshold is None:
        #     the_threshold = 0.7
        Classifier = LinearSVC(C=1.2)
        clf = CalibratedClassifierCV(Classifier)
        start_time = time.time()
        clf.fit(train_data, train_label)
        
        del train_data,train_label,valid_data,valid_label
        gc.collect()
        print_memory_used()
        # print(f'num_classes = {num_classes}')
        
        model = clf
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"time cost：{execution_time} s")
        start_time = time.time()
        
        predicted = clf.predict(test_data)
        
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"time cost：{execution_time} s")
        
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        test_label = test_label.values.tolist()
        predicted = predicted.tolist()
        test_label = change_list_list_to_list(input_list=test_label)
        predicted = change_list_list_to_list(input_list=predicted)
        print(type(test_label))
        print(type(predicted))
        print_list(input_list=test_label)
        print_list(input_list=predicted)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped

    def train_PCMA_MYGO_RC_with_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_batch_size=None,num_classes=None,lr=None,dropout_ratio1=None,dropout_ratio2=None,dropout_ratio3=None,
        epochs=None,gpu_code=None,the_num_workers=None,optimizer=None, momentum=None, the_threshold1=None, the_threshold2=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        # if the_threshold1 is None:
        #     the_threshold1 = 0.5
        # if the_threshold2 is None:
        #     the_threshold2 = 0.7
        traindata = Dataset_PCMA_MYGO_RC(train_data,train_label) 
        validdata = Dataset_PCMA_MYGO_RC(valid_data,valid_label)
        
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        
        del train_data,train_label,valid_data,valid_label,traindata,validdata
        gc.collect()
        print_memory_used()
        print(f'num_classes = {num_classes}')
        
        pcma_mygo_rc = PCMA_MYGO_RC(
            dropout_ratio1=dropout_ratio1, dropout_ratio2=dropout_ratio2, dropout_ratio3=dropout_ratio3, num_types=num_classes,
        )
        pretrained_dict = torch.load('resnet18.pth')
        model_dict = pcma_mygo_rc.resnet18.state_dict()
        pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
        model_dict.update(pretrained_dict)
        pcma_mygo_rc.resnet18.load_state_dict(model_dict)
        torch.cuda.init()
        
        # net=None, train_iter=None, valid_iter=None, num_epochs=None, lr=None, device=None, optimizer='Adam', momentum=None
        train_with_GPU(net=pcma_mygo_rc, train_iter=train_iter, valid_iter=valid_iter, 
                       num_epochs=epochs, lr=lr, device=d2l.try_all_gpus()[gpu_code], optimizer=optimizer, momentum=momentum)
        model = pcma_mygo_rc
        del train_sampler,train_iter,valid_sampler,valid_iter
        gc.collect()
        print_memory_used()
        
        import math
        
        batch_size = 32 
        total_size = test_data.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = model
        net_try = net_try.to('cpu')
        net_try.eval()
        
        with torch.no_grad():
            predicted_results = []
            for i in range(total_batches):
                start = i * batch_size
                end = min(start + batch_size, total_size)
                batch_data = test_data[start:end]
                batch_label = test_label[start:end]
                batch_data = torch.from_numpy(np.float32(np.array(batch_data)))
                batch_data = torch.reshape(batch_data, (batch_data.shape[0], 224, 224))
                batch_data = torch.unsqueeze(batch_data, 1)
                batch_data = torch.cat((batch_data, batch_data, batch_data), dim=1)
                batch_label = torch.from_numpy(np.int64(np.array(batch_label)).reshape(1,-1)[0])
                outputs = net_try(batch_data)
                outputs = torch.nn.functional.softmax(outputs, dim=1)
                # print(outputs)
                # print('type(outputs)')
                # print(type(outputs))
                max_values, predicted = torch.max(outputs, dim=1)  
                max_values_list = max_values.numpy().tolist()  
                max_values_list = [round(num, 4) for num in max_values_list]  
                the_median = np.median(max_values_list)
                _, indices = torch.topk(outputs, 2, dim=1)  
                second_largest_values = torch.gather(outputs, 1, indices[:, 1].unsqueeze(1))   
                second_largest_values_list = second_largest_values.numpy().flatten().tolist()  
                second_largest_values_list = [round(num, 4) for num in second_largest_values_list]  
                second_median = np.median(second_largest_values_list)   
                # print('predicted')
                # print(predicted)
                # print('max_values')
                # print(max_values_list)
                # print('the_median')
                # print(the_median)
                # print('second_largest_values')
                # print(second_largest_values_list)
                # print('second_median')
                # print(second_median)
                predicted_list = predicted.numpy().flatten().tolist()
                predicted_list_with_unknown = []
                for i,j in enumerate(predicted_list):
                    if max_values_list[i] > 3*second_largest_values_list[i]:
                        predicted_list_with_unknown.append(j)
                    else:
                        predicted_list_with_unknown.append(-100)
                # print('predicted_list_with_unknown')
                # print(predicted_list_with_unknown)
                # print(predicted)
                # the_threshold1 = the_threshold1*the_median
                # the_threshold2 = the_threshold2*the_median
                # predicted[max_values < the_threshold1] = -100 
                # num_over_09 = (outputs > the_threshold2).sum(dim=1) 
                # predicted[num_over_09 >= 2] = -100
                # print(predicted)
                predicted_results.append(np.array(predicted_list_with_unknown))
                # break
                
        predicted = np.concatenate(predicted_results)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        test_label = test_label.numpy()
        
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        test_label = test_label.tolist()
        predicted = predicted.tolist()
        test_label = change_list_list_to_list(input_list=test_label)
        predicted = change_list_list_to_list(input_list=predicted)
        print(type(test_label))
        print(type(predicted))
        print_list(input_list=test_label)
        print_list(input_list=predicted)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped

    def train_PCMA_MYGO_RC_no_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_batch_size=None,num_classes=None,lr=None,dropout_ratio1=None,dropout_ratio2=None,dropout_ratio3=None,
        epochs=None,gpu_code=None,the_num_workers=None,optimizer=None, momentum=None, the_threshold1=None, the_threshold2=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        # if the_threshold1 is None:
        #     the_threshold1 = 0.5
        # if the_threshold2 is None:
        #     the_threshold2 = 0.7
        traindata = Dataset_PCMA_MYGO_RC(train_data,train_label) 
        validdata = Dataset_PCMA_MYGO_RC(valid_data,valid_label)
        
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        
        del train_data,train_label,valid_data,valid_label,traindata,validdata
        gc.collect()
        print_memory_used()
        print(f'num_classes = {num_classes}')
        
        pcma_mygo_rc = PCMA_MYGO_RC(
            dropout_ratio1=dropout_ratio1, dropout_ratio2=dropout_ratio2, dropout_ratio3=dropout_ratio3, num_types=num_classes,
        )
        pretrained_dict = torch.load('resnet18.pth')
        model_dict = pcma_mygo_rc.resnet18.state_dict()
        pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
        model_dict.update(pretrained_dict)
        pcma_mygo_rc.resnet18.load_state_dict(model_dict)
        torch.cuda.init()
        
        # net=None, train_iter=None, valid_iter=None, num_epochs=None, lr=None, device=None, optimizer='Adam', momentum=None
        train_with_GPU(net=pcma_mygo_rc, train_iter=train_iter, valid_iter=valid_iter, 
                       num_epochs=epochs, lr=lr, device=d2l.try_all_gpus()[gpu_code], optimizer=optimizer, momentum=momentum)
        model = pcma_mygo_rc
        del train_sampler,train_iter,valid_sampler,valid_iter
        gc.collect()
        print_memory_used()
        
        import math
        
        batch_size = 32 
        total_size = test_data.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = model
        net_try = net_try.to('cpu')
        net_try.eval()
        
        with torch.no_grad():
            predicted_results = []
            for i in range(total_batches):
                start = i * batch_size
                end = min(start + batch_size, total_size)
                batch_data = test_data[start:end]
                batch_label = test_label[start:end]
                batch_data = torch.from_numpy(np.float32(np.array(batch_data)))
                batch_data = torch.reshape(batch_data, (batch_data.shape[0], 224, 224))
                batch_data = torch.unsqueeze(batch_data, 1)
                batch_data = torch.cat((batch_data, batch_data, batch_data), dim=1)
                batch_label = torch.from_numpy(np.int64(np.array(batch_label)).reshape(1,-1)[0])
                outputs = net_try(batch_data)
                _, predicted = torch.max(outputs, 1)
                predicted_results.append(predicted.numpy())
                
        predicted = np.concatenate(predicted_results)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        test_label = test_label.numpy()
        
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        test_label = test_label.tolist()
        predicted = predicted.tolist()
        test_label = change_list_list_to_list(input_list=test_label)
        predicted = change_list_list_to_list(input_list=predicted)
        print(type(test_label))
        print(type(predicted))
        print_list(input_list=test_label)
        print_list(input_list=predicted)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped

    def train_PCMA_MYGO_RC2_with_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_batch_size=None,num_classes=None,lr=None,dropout_ratio1=None,dropout_ratio2=None,dropout_ratio3=None,
        epochs=None,gpu_code=None,the_num_workers=None,optimizer=None, momentum=None, the_threshold1=None, the_threshold2=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        # if the_threshold1 is None:
        #     the_threshold1 = 0.8
        # if the_threshold2 is None:
        #     the_threshold2 = 0.8
        traindata = Dataset_PCMA_MYGO_RC_2(train_data,train_label) 
        validdata = Dataset_PCMA_MYGO_RC_2(valid_data,valid_label)
        
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        
        del train_data,train_label,valid_data,valid_label,traindata,validdata
        gc.collect()
        print_memory_used()
        print(f'num_classes = {num_classes}')
        
        pcma_mygo_rc = PCMA_MYGO_RC2(
            dropout_ratio1=dropout_ratio1, dropout_ratio2=dropout_ratio2, dropout_ratio3=dropout_ratio3, num_types=num_classes,
        )
        
        torch.cuda.init()
        # net=None, train_iter=None, valid_iter=None, num_epochs=None, lr=None, device=None, optimizer='Adam', momentum=None
        train_with_GPU(net=pcma_mygo_rc, train_iter=train_iter, valid_iter=valid_iter, 
                       num_epochs=epochs, lr=lr, device=d2l.try_all_gpus()[gpu_code], optimizer=optimizer, momentum=momentum)
        model = pcma_mygo_rc
        del train_sampler,train_iter,valid_sampler,valid_iter
        gc.collect()
        print_memory_used()
        
        import math
        
        batch_size = 32 
        total_size = test_data.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = model
        net_try = net_try.to('cpu')
        net_try.eval()
        
        with torch.no_grad():
            predicted_results = []
            for i in range(total_batches):
                start = i * batch_size
                end = min(start + batch_size, total_size)
                batch_data = test_data[start:end]
                batch_label = test_label[start:end]
                batch_data = torch.from_numpy(np.float32(np.array(batch_data)))
                new_tensor = torch.zeros(batch_data.size(0), 3*224*224)
                new_tensor[:, :batch_data.size(1)] = batch_data
                batch_data = new_tensor
                batch_data = batch_data.view(batch_data.size(0), 3, 224, 224)
                batch_label = torch.from_numpy(np.int64(np.array(batch_label)).reshape(1,-1)[0])
                outputs = net_try(batch_data)
                outputs = torch.nn.functional.softmax(outputs, dim=1)
                # print(outputs)
                # print('type(outputs)')
                # print(type(outputs))
                max_values, predicted = torch.max(outputs, dim=1)  
                max_values_list = max_values.numpy().tolist()  
                max_values_list = [round(num, 4) for num in max_values_list]  
                the_median = np.median(max_values_list)
                _, indices = torch.topk(outputs, 2, dim=1)  
                second_largest_values = torch.gather(outputs, 1, indices[:, 1].unsqueeze(1))   
                second_largest_values_list = second_largest_values.numpy().flatten().tolist()  
                second_largest_values_list = [round(num, 4) for num in second_largest_values_list]  
                second_median = np.median(second_largest_values_list)   
                # print('predicted')
                # print(predicted)
                # print('max_values')
                # print(max_values_list)
                # print('the_median')
                # print(the_median)
                # print('second_largest_values')
                # print(second_largest_values_list)
                # print('second_median')
                # print(second_median)
                predicted_list = predicted.numpy().flatten().tolist()
                predicted_list_with_unknown = []
                for i,j in enumerate(predicted_list):
                    if max_values_list[i] > 3*second_largest_values_list[i]:
                        predicted_list_with_unknown.append(j)
                    else:
                        predicted_list_with_unknown.append(-100)
                # print('predicted_list_with_unknown')
                # print(predicted_list_with_unknown)
                # print(predicted)
                # the_threshold1 = the_threshold1*the_median
                # the_threshold2 = the_threshold2*the_median
                # predicted[max_values < the_threshold1] = -100 
                # num_over_09 = (outputs > the_threshold2).sum(dim=1) 
                # predicted[num_over_09 >= 2] = -100
                # print(predicted)
                predicted_results.append(np.array(predicted_list_with_unknown))
                # break
                
        predicted = np.concatenate(predicted_results)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        test_label = test_label.numpy()
        
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        test_label = test_label.tolist()
        predicted = predicted.tolist()
        test_label = change_list_list_to_list(input_list=test_label)
        predicted = change_list_list_to_list(input_list=predicted)
        print(type(test_label))
        print(type(predicted))
        print_list(input_list=test_label)
        print_list(input_list=predicted)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped

    def train_PCMA_MYGO_RC2_no_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_batch_size=None,num_classes=None,lr=None,dropout_ratio1=None,dropout_ratio2=None,dropout_ratio3=None,
        epochs=None,gpu_code=None,the_num_workers=None,optimizer=None, momentum=None, the_threshold1=None, the_threshold2=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        # if the_threshold1 is None:
        #     the_threshold1 = 0.8
        # if the_threshold2 is None:
        #     the_threshold2 = 0.8
        traindata = Dataset_PCMA_MYGO_RC_2(train_data,train_label) 
        validdata = Dataset_PCMA_MYGO_RC_2(valid_data,valid_label)
        
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        
        del train_data,train_label,valid_data,valid_label,traindata,validdata
        gc.collect()
        print_memory_used()
        print(f'num_classes = {num_classes}')
        
        pcma_mygo_rc = PCMA_MYGO_RC2(
            dropout_ratio1=dropout_ratio1, dropout_ratio2=dropout_ratio2, dropout_ratio3=dropout_ratio3, num_types=num_classes,
        )
        
        torch.cuda.init()
        # net=None, train_iter=None, valid_iter=None, num_epochs=None, lr=None, device=None, optimizer='Adam', momentum=None
        train_with_GPU(net=pcma_mygo_rc, train_iter=train_iter, valid_iter=valid_iter, 
                       num_epochs=epochs, lr=lr, device=d2l.try_all_gpus()[gpu_code], optimizer=optimizer, momentum=momentum)
        model = pcma_mygo_rc
        del train_sampler,train_iter,valid_sampler,valid_iter
        gc.collect()
        print_memory_used()
        
        import math
        
        batch_size = 32 
        total_size = test_data.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = model
        net_try = net_try.to('cpu')
        net_try.eval()
        
        with torch.no_grad():
            predicted_results = []
            for i in range(total_batches):
                start = i * batch_size
                end = min(start + batch_size, total_size)
                batch_data = test_data[start:end]
                batch_label = test_label[start:end]
                batch_data = torch.from_numpy(np.float32(np.array(batch_data)))
                new_tensor = torch.zeros(batch_data.size(0), 3*224*224)
                new_tensor[:, :batch_data.size(1)] = batch_data
                batch_data = new_tensor
                batch_data = batch_data.view(batch_data.size(0), 3, 224, 224)
                batch_label = torch.from_numpy(np.int64(np.array(batch_label)).reshape(1,-1)[0])
                outputs = net_try(batch_data)
                _, predicted = torch.max(outputs, 1)
                predicted_results.append(predicted.numpy())
                
        predicted = np.concatenate(predicted_results)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        test_label = test_label.numpy()
        
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        test_label = test_label.tolist()
        predicted = predicted.tolist()
        test_label = change_list_list_to_list(input_list=test_label)
        predicted = change_list_list_to_list(input_list=predicted)
        print(type(test_label))
        print(type(predicted))
        print_list(input_list=test_label)
        print_list(input_list=predicted)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped

    def train_PCMA_MYGO_SA_with_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_batch_size=None,num_classes=None,lr=None,dropout_ratio1=None,dropout_ratio2=None,dropout_ratio3=None,feature_len=None,
        epochs=None,gpu_code=None,the_num_workers=None,the_num_heads=None,optimizer=None, momentum=None, the_threshold1=None, the_threshold2=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        if the_threshold1 is None:
            the_threshold1 = 0.8
        if the_threshold2 is None:
            the_threshold2 = 0.8
        traindata = Dataset_transformer(train_data,train_label) 
        validdata = Dataset_transformer(valid_data,valid_label)
        
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        del traindata,train_data,train_label,validdata,valid_data,valid_label
        gc.collect() 
        print_memory_used()
        
        s = feature_len//2
        heads = the_num_heads
        dp = dropout_ratio1
        batch_sizes = the_batch_size
        num_classes = num_classes
        print(f'num_classes = {num_classes}')
        
        pcma_mygo_sa = PCMA_MYGO_SA(input_dim=s, num_head=heads, dim1=s,
                           num_classes=num_classes,dropout=dp)
        
        torch.cuda.init()
        # net=None, train_iter=None, valid_iter=None, num_epochs=None, lr=None, device=None, optimizer='Adam', momentum=None
        train_with_GPU(net=pcma_mygo_sa, train_iter=train_iter, valid_iter=valid_iter, 
                       num_epochs=epochs, lr=lr, device=d2l.try_all_gpus()[gpu_code], optimizer=optimizer, momentum=momentum)
        model = pcma_mygo_sa
        del train_sampler,train_iter,valid_sampler,valid_iter
        gc.collect()
        print_memory_used()
        
        import math
        
        batch_size = 32 
        total_size = test_data.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = model
        net_try = net_try.to('cpu')
        net_try.eval()
        
        with torch.no_grad():
            predicted_results = []
            for i in range(total_batches):
                start = i * batch_size
                end = min(start + batch_size, total_size)
                batch_data = test_data[start:end]
                batch_label = test_label[start:end]
                batch_data = torch.from_numpy(np.float32(np.array(batch_data)))
                batch_data = torch.reshape(batch_data, (batch_data.shape[0], 2, 1024))
                batch_label = torch.from_numpy(np.int64(np.array(batch_label)).reshape(1,-1)[0])
                outputs = net_try(batch_data)
                # print(outputs)
                # outputs = torch.sigmoid(outputs)
                outputs = torch.nn.functional.softmax(outputs, dim=1)
                # print(outputs)
                _, predicted = torch.max(outputs, 1)
                # print(predicted)
                max_vals, _ = torch.max(outputs, dim=1)
                predicted[max_vals < the_threshold1] = -100 
                num_over_09 = (outputs > the_threshold2).sum(dim=1) 
                predicted[num_over_09 >= 2] = -100
                # print(predicted)
                predicted_results.append(predicted.numpy())
                # break
        predicted = np.concatenate(predicted_results)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        test_label = test_label.numpy()
        
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        test_label = test_label.tolist()
        predicted = predicted.tolist()
        test_label = change_list_list_to_list(input_list=test_label)
        predicted = change_list_list_to_list(input_list=predicted)
        print(type(test_label))
        print(type(predicted))
        print_list(input_list=test_label)
        print_list(input_list=predicted)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped

    def train_PCMA_MYGO_SA_no_unknown(
        train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
        the_batch_size=None,num_classes=None,lr=None,dropout_ratio1=None,dropout_ratio2=None,dropout_ratio3=None,feature_len=None,
        epochs=None,gpu_code=None,the_num_workers=None,the_num_heads=None,optimizer=None, momentum=None, the_threshold1=None, the_threshold2=None,
        mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
    ):
        # if the_threshold1 is None:
        #     the_threshold1 = 0.8
        # if the_threshold2 is None:
        #     the_threshold2 = 0.8
        traindata = Dataset_transformer(train_data,train_label) 
        validdata = Dataset_transformer(valid_data,valid_label)
        
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        del traindata,train_data,train_label,validdata,valid_data,valid_label
        gc.collect() 
        print_memory_used()
        
        s = feature_len//2
        heads = the_num_heads
        dp = dropout_ratio1
        batch_sizes = the_batch_size
        num_classes = num_classes
        print(f'num_classes = {num_classes}')
        
        pcma_mygo_sa = PCMA_MYGO_SA(input_dim=s, num_head=heads, dim1=s,
                           num_classes=num_classes,dropout=dp)
        
        torch.cuda.init()
        # net=None, train_iter=None, valid_iter=None, num_epochs=None, lr=None, device=None, optimizer='Adam', momentum=None
        train_with_GPU(net=pcma_mygo_sa, train_iter=train_iter, valid_iter=valid_iter, 
                       num_epochs=epochs, lr=lr, device=d2l.try_all_gpus()[gpu_code], optimizer=optimizer, momentum=momentum)
        model = pcma_mygo_sa
        del train_sampler,train_iter,valid_sampler,valid_iter
        gc.collect()
        print_memory_used()
        
        import math
        
        batch_size = 32 
        total_size = test_data.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = model
        net_try = net_try.to('cpu')
        net_try.eval()
        
        with torch.no_grad():
            predicted_results = []
            for i in range(total_batches):
                start = i * batch_size
                end = min(start + batch_size, total_size)
                batch_data = test_data[start:end]
                batch_label = test_label[start:end]
                batch_data = torch.from_numpy(np.float32(np.array(batch_data)))
                batch_data = torch.reshape(batch_data, (batch_data.shape[0], 2, 1024))
                batch_label = torch.from_numpy(np.int64(np.array(batch_label)).reshape(1,-1)[0])
                outputs = net_try(batch_data)
                _, predicted = torch.max(outputs, 1)
                predicted_results.append(predicted.numpy())
                
        predicted = np.concatenate(predicted_results)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        test_label = test_label.numpy()
        
        print('mapping_1_test,mapping_2_test')
        print(mapping_1_test)
        print(mapping_2_test)
        # print('predicted')
        # print_list(input_list=predicted)
        print('mapping_1,mapping_2')
        print(mapping_1)
        print(mapping_2)
        
        test_label = test_label.tolist()
        predicted = predicted.tolist()
        test_label = change_list_list_to_list(input_list=test_label)
        predicted = change_list_list_to_list(input_list=predicted)
        print(type(test_label))
        print(type(predicted))
        print_list(input_list=test_label)
        print_list(input_list=predicted)
        
        mapping_1['Unknown']=-100
        mapping_2[-100]=-100
        test_label_unmapped = back_to_unmapped(obj=test_label,mapping_1=mapping_1_test,mapping_2=mapping_2_test)
        predicted_unmapped = back_to_unmapped(obj=predicted,mapping_1=mapping_1,mapping_2=mapping_2)
        
        return model,test_label_unmapped,predicted_unmapped

def calculate_f1_for_each_class(true_list, pred_list):  
    f1_scores = {}  
    classes = set(true_list + pred_list)  # Get all unique classes  
    for class_ in classes:  
        true_class = [1 if x == class_ else 0 for x in true_list]  # Convert true labels to binary array  
        pred_class = [1 if x == class_ else 0 for x in pred_list]  # Convert predicted labels to binary array  
        f1_scores[class_] = f1_score(true_class, pred_class)  # Calculate F1 score and store  
    return f1_scores  

def get_accuracy_precision_recall_f1_score_adata1_adata2(
    test_label_unmapped=None,predicted_unmapped=None,
    ref_mapping=None,# maybe aaaa is in ref mapping, but not in pred!
    unknown_model=True,unknown_label=None,unknown_labels_set=None):
    
    model_result_list = []
    if unknown_label is None:
        unknown_label = 'Unknown'
    if unknown_labels_set is None:
        unknown_labels_set = {'unknown','Unknown','rejected','ambiguous'}
    unknown_accuracy = -1.0
    unknown_precision = -1.0
    unknown_recall = -1.0
    unknown_f1_score = -1.0
    true_common_in_test_total_rate = -1.0
    predicted_common_in_true_common_rate = -1.0
    
    _,test_mapping,_ = generate_mapping_for_list(input_list=test_label_unmapped)
    merged_mapping,set_common_keys,set_dict1_only_keys,set_dict2_only_keys = merge_dicts(dict1=ref_mapping,dict2=test_mapping)
    print('merged_mapping')
    print(merged_mapping)
    print('set_common_keys')
    print(set_common_keys)
    print('set_ref_only_keys')
    print(set_dict1_only_keys)
    print('set_test_only_keys')
    print(set_dict2_only_keys)
    
    pre_unknown = []
    test_unknown = []
    for i,j in enumerate(predicted_unmapped):
        if j in unknown_labels_set and test_label_unmapped[i] in set_dict2_only_keys:
            pre_unknown.append(unknown_label)
            test_unknown.append(unknown_label)
        elif j in unknown_labels_set and test_label_unmapped[i] in set_common_keys:
            pre_unknown.append(unknown_label)
            test_unknown.append('Not_'+unknown_label)
        elif j not in unknown_labels_set and test_label_unmapped[i] in set_dict2_only_keys:
            pre_unknown.append('Not_'+unknown_label)
            test_unknown.append(unknown_label)
        elif j not in unknown_labels_set and test_label_unmapped[i] in set_common_keys:
            pre_unknown.append('Not_'+unknown_label)
            test_unknown.append('Not_'+unknown_label)
    unknown_mapping = {unknown_label:0,('Not_'+unknown_label):1}
    pre_unknown = [unknown_mapping[key] for key in pre_unknown]
    test_unknown = [unknown_mapping[key] for key in test_unknown]
    unknown_accuracy = accuracy_score(test_unknown, pre_unknown)
    unknown_precision = precision_score(test_unknown, pre_unknown, average='macro')
    unknown_recall = recall_score(test_unknown, pre_unknown, average='macro')
    unknown_f1_score = f1_score(test_unknown, pre_unknown, average='macro')
    if_ref_common = 0
    if_test_common = 0
    for i in ref_mapping.keys():
        if i not in set_common_keys:
            if_ref_common = 1
    for i in test_mapping.keys():
        if i not in set_common_keys:
            if_test_common = 1
    if if_ref_common == 0 and if_test_common == 0:
        print('ref labels and test labels are the same.')
        unknown_accuracy = -1.0
        unknown_precision = -1.0   
        unknown_recall = -1.0  
        unknown_f1_score = -1.0
    print("Unknown Accuracy:", unknown_accuracy)  
    print("Unknown Macro Precision:", unknown_precision)  
    print("Unknown Macro Recall:", unknown_recall)  
    print("Unknown Macro F1 Score:", unknown_f1_score)
    model_result_list.append('unknown_accuracy: '+str(unknown_accuracy)+' []')
    model_result_list.append('unknown_precision: '+str(unknown_precision)+' []')
    model_result_list.append('unknown_recall: '+str(unknown_recall)+' []')
    model_result_list.append('unknown_f1_score: '+str(unknown_f1_score)+' []')
    
    test_total_count = 0
    true_common = 0
    predicted_common_and_true_common = 0
    for i,j in enumerate(predicted_unmapped):
        test_total_count = test_total_count+1
        if test_label_unmapped[i] in set_common_keys: 
            true_common = true_common+1
            if j not in unknown_labels_set and j in set_common_keys:
                predicted_common_and_true_common = predicted_common_and_true_common+1
    print('test_total_count')
    print(test_total_count)
    print('true_common')
    print(true_common)
    print('predicted_common_and_true_common')
    print(predicted_common_and_true_common)
    model_result_list.append('test_total_count: '+str(test_total_count)+' []')
    model_result_list.append('true_common: '+str(true_common)+' []')
    model_result_list.append('predicted_common_and_true_common: '+str(predicted_common_and_true_common)+' []')
    
    # merge test and pred, when test_label is not in set_common_keys and pred is in unknown_labels_set, consider pred as test_label
    # # when in set_common_keys, normally calculate acc f1, is that appropriate?
    tmp_test = []
    tmp_pred = []
    for i,j in enumerate(predicted_unmapped):
        testtest = test_label_unmapped[i]
        predpred = j
        if testtest in set_common_keys and predpred in unknown_labels_set: 
            tmp_test.append(testtest)
            tmp_pred.append(unknown_label) # change different unknown labels to the same unknown label
        elif testtest in set_common_keys and predpred not in unknown_labels_set:
            tmp_test.append(testtest)
            tmp_pred.append(predpred)
        elif testtest not in set_common_keys and predpred in unknown_labels_set:
            tmp_test.append(testtest)
            tmp_pred.append(testtest)
        elif testtest not in set_common_keys and predpred not in unknown_labels_set:
            tmp_test.append(testtest)
            tmp_pred.append(predpred)
    _,pred_mapping,_ = generate_mapping_for_list(input_list=tmp_pred)
    merged_mapping_test_and_pred,_,_,_ = merge_dicts(dict1=pred_mapping,dict2=test_mapping)
    test_label = [merged_mapping_test_and_pred[key] for key in tmp_test]
    predicted = [merged_mapping_test_and_pred[key] for key in tmp_pred]
    accuracy = accuracy_score(test_label, predicted)  
    macro_precision = precision_score(test_label, predicted, average='macro', zero_division=0)   
    macro_recall = recall_score(test_label, predicted, average='macro', zero_division=0)  
    macro_f1_score = f1_score(test_label, predicted, average='macro', zero_division=0)  
    model_result_list.append('accuracy: '+str(accuracy)+' []')
    model_result_list.append('macro_precision: '+str(macro_precision)+' []')
    model_result_list.append('macro_recall: '+str(macro_recall)+' []')
    model_result_list.append('macro_f1_score: '+str(macro_f1_score)+' []')
    print("Accuracy:", accuracy)  
    print("Macro Precision:", macro_precision)  
    print("Macro Recall:", macro_recall)  
    print("Macro F1 Score:", macro_f1_score)
    
    all_f1_scores = f1_score(test_label, predicted, average=None)  
    median_f1_score = np.median(all_f1_scores)  
    model_result_list.append('median_f1_score: '+str(median_f1_score)+' []')
    print("Median F1 Score:", median_f1_score)
    
    dict_label_count = {}
    for label in test_label_unmapped:
        if label not in dict_label_count:
            dict_label_count[label]=1
        else:
            dict_label_count[label] = dict_label_count[label]+1
    print('test_label_count(celltype-count):')
    print(dict_label_count)
    
    sorted_labels = sorted(dict_label_count, key=dict_label_count.get)  
    top3_labels = sorted_labels[:3]  
    top1_labels = sorted_labels[:1]  
    print('three cell types with the least number:')
    print(top3_labels) 
    
    # all_f1_scores_dict = calculate_f1_for_each_class(test_label_unmapped, predicted_unmapped)
    all_f1_scores_dict = calculate_f1_for_each_class(true_list=tmp_test, pred_list=tmp_pred)
    print("All F1 Scores:", all_f1_scores_dict)

    all_f1_scores_dict_filtered = {}
    for tmpkey in all_f1_scores_dict:
        if tmpkey in test_mapping:
            all_f1_scores_dict_filtered[tmpkey] = all_f1_scores_dict[tmpkey]
    print("All F1 Scores filtered:", all_f1_scores_dict_filtered)
    tmpf1scorelist = []
    for tmpv in all_f1_scores_dict_filtered.values():
        tmpf1scorelist.append(tmpv)
    macro_f1_adjusted = np.average(tmpf1scorelist)
    median_f1_adjusted = np.median(tmpf1scorelist)
    model_result_list.append('macro_f1_score_adjusted: '+str(macro_f1_adjusted)+' []')
    print("Macro F1 Score adjusted:", macro_f1_adjusted)
    model_result_list.append('median_f1_score_adjusted: '+str(median_f1_adjusted)+' []')
    print("Median F1 Score adjusted:", median_f1_adjusted)
    
    if len(top3_labels) == 3:
        f1_score_least_1 = all_f1_scores_dict[top3_labels[0]]
        f1_score_least_2 = all_f1_scores_dict[top3_labels[1]]
        f1_score_least_3 = all_f1_scores_dict[top3_labels[2]]
        print('f1 score of the cell type with the least number:')
        print(top3_labels[0])
        print(f1_score_least_1)
        print('f1 score of the cell type with the second least number:')
        print(top3_labels[1])
        print(f1_score_least_2)
        print('f1 score of the cell type with the third least number:')
        print(top3_labels[2])
        print(f1_score_least_3)
        average_f1_score_three_least = (f1_score_least_1+f1_score_least_2+f1_score_least_3)/3.0
        print("Average f1 score of three least cell types:",average_f1_score_three_least)
    else:
        print('The num of cell types in the unknown dataset is less than three!')
        # print('The num of common cell types between the ref dataset and the unknown dataset is less than three!')
        average_f1_score_three_least = -1
    model_result_list.append('average_f1_score_three_least: '+str(average_f1_score_three_least)+' []')
    
    if len(top1_labels) == 1:
        f1_score_least_1 = all_f1_scores_dict[top3_labels[0]]
        print('f1 score of the cell type with the least number:')
        print(top3_labels[0])
        print(f1_score_least_1)
    else:
        # print('The num of common cell types between the ref dataset and the unknown dataset is less than one!')
        print('The num of cell types in the unknown dataset is less than one!')
        f1_score_least_1 = -1
    model_result_list.append('f1_score_least_1: '+str(f1_score_least_1)+' []')
    
    return model_result_list

def cal_acc_f1_simply(test_label_unmapped,predicted_unmapped):
    test_label,mapping_1_test,mapping_2_test = generate_mapping_for_list(input_list=test_label_unmapped)
    predicted = [mapping_1_test[key] for key in predicted_unmapped]
    from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score  
    from sklearn.metrics import precision_recall_fscore_support  
      
    # Assuming y_true and y_pred are your actual and predicted labels  
    accuracy = accuracy_score(test_label, predicted)  
    precision, recall, f1, _ = precision_recall_fscore_support(test_label, predicted, average='macro')  
      
    print("Accuracy:", accuracy)  
    print("Macro Precision:", precision)  
    print("Macro Recall:", recall)  
    print("Macro F1 Score:", f1)  

    