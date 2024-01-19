# -*- coding: utf-8 -*-
# Some references：
# D2l. https://zh-v2.d2l.ai/
# Pandas. https://pandas.pydata.org/docs/index.html
# PyTorch. https://pytorch.org/
# Scikit-learn. https://scikit-learn.org/stable/

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
import datetime
from torch import nn
from d2l import torch as d2l
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Dataset, SequentialSampler, RandomSampler

from IPython.core.interactiveshell import InteractiveShell

InteractiveShell.ast_node_interactivity = "all"

now = datetime.datetime.now()
current_time = now.strftime("%Y-%m-%d %H:%M:%S")
print(current_time)

gc.collect()

# 获取当前节点的可用物理内存大小（单位为字节）
available_memory = psutil.virtual_memory().available

# 将字节转换为 GB
available_memory_gb = available_memory / (1024 * 1024 * 1024)

print(f"当前可用的最大内存为: {available_memory_gb} GB")

# 获取当前节点已使用的物理内存大小（单位为字节）
used_memory = psutil.virtual_memory().used

# 将字节转换为 GB
used_memory_gb = used_memory / (1024 * 1024 * 1024)

print(f"当前已使用的内存为: {used_memory_gb} GB")


class PCmaster_anno_0_Error(Exception):
    pass


class PCmaster_anno_0():
    # 下面是PCmaster_anno_0的现有属性,emmm,其实有一些是暂时用不到的
    result_file = None
    shape = None
    # 放单数据输入时的scanpy对象,或者放多数据输入时整合后的scanpy对象 
    obj = None 
    # 放单数据输入的物种
    species = None
    # 放单数据输入的器官
    organ = None
    dd = None  # dense dataframe, pandas object
    genes = None
    cells = None
    cellnamelabel = None
    # 放多数据输入时整合前的scanpy对象
    objs = []
    # 放多数据输入时的物种
    batch_species = []
    # 放多数据输入时的器官
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
    
    # 机器学习模型相关
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

    # 看看内存情况
    def memory_collect_0(self):
        gc.collect()

        # 获取当前节点的可用物理内存大小（单位为字节）
        available_memory = psutil.virtual_memory().available
        
        # 将字节转换为 GB
        available_memory_gb = available_memory / (1024 * 1024 * 1024)
        
        print(f"当前可用的最大内存为: {available_memory_gb} GB")
        
        # 获取当前节点已使用的物理内存大小（单位为字节）
        used_memory = psutil.virtual_memory().used
        
        # 将字节转换为 GB
        used_memory_gb = used_memory / (1024 * 1024 * 1024)
        
        print(f"当前已使用的内存为: {used_memory_gb} GB")

    # 目前用不太到
    def set_result_file_0(self, result_file=None):
        if result_file is not None:
            self.result_file = result_file
        else:
            raise PCmaster_anno_0_Error("The result_file is None.")

    # 设置dpi等参数
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

    # 数据导入
    # 这里包含的就是字符串从后往前截取,以及字符串a是否在字符串b之中
    # 这里的sc就是scanpy
    def data_input_0(self, filepath=None, species=None, organ=None, organs=None):
        if filepath is not None:
            if filepath[-6:] == "tsv.gz":    # 字符串从后往前截取
                self.obj = sc.read_umi_tools(filepath)
            elif filepath[-6:] == "mtx.gz":
                self.obj = sc.read_10x_mtx(filepath[:-13], var_names='gene_symbols', cache=True)
            elif 'filtered_feature_bc_matrix' in filepath:    # 字符串a是否在字符串b之中
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

    # 批量数据导入
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

    # 这里是后来补的,通过scanpy的数据格式anndata导入
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

    # 拼接批量数据,用scanpy的方法,加了batch信息
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

    # 目前用不到
    def get_dense_0(self):
        if self.obj is None:
            raise PCmaster_anno_0_Error("The expression matrix is not successfully loaded.")
        self.dd = pd.DataFrame(self.obj.X.todense(), index=self.obj.obs_names, columns=self.obj.var_names)

    # 预处理,过滤相关,用的基本是scanpy的参数
    # 增加了通过绝对中位差MAD过滤的功能
    # MAD = 所有xi-x中位数之后的值的绝对值的中位数       
    # 1,2,3 中位数为2,各项减去2,得-1,0,1,取绝对值后排序,0,1,1,则MAD = 1
    # 该功能就是保留各项数值在+-mad_coeffcient*MAD范围内的数据
    # 一般情况下不用通过绝对中位差MAD过滤的功能也没关系
    def filter_0(self, n_top=None, num_min_genes=None, num_max_genes=None, num_min_cells=None,
                 max_percent_of_specific_genes=None,
                 genelist=None, mad_coeffcient=None):
        # tmp = self.obj.copy()
        # 参数初始化
        if n_top is None:
            n_top = 20
        if num_min_genes is None:
            num_min_genes = 200
        if num_min_cells is None:
            num_min_cells = 3
        # 下面这个参数是用来过滤特定基因集合的,比如线粒体基因等等
        if max_percent_of_specific_genes is None:
            max_percent_of_specific_genes = 0.05
        # 获取scanpy对象的基因名字,用var_names也可以
        self.genes = self.obj.var.index
        
        mad_min_genes = 0
        mad_max_genes = 0
        mad_min_counts = 0
        mad_max_counts = 0
        mad_max_pct = 0
        
        # 这里是分成了过滤前数据情况和过滤后数据情况,分别展示了一下
        print("\n")
        print("Before filter")
        print("\n")
        sc.pl.highest_expr_genes(self.obj, n_top, )
        sc.pp.filter_cells(self.obj, min_genes=-1)
        sc.pp.filter_genes(self.obj, min_cells=-1)
        if genelist is None:
            # 这里其实不应该这样写,万一别人也这样标识了一下基因,我就会执行过滤
            # 下面这一行命令会识别所有以'meaningless-'开头的字符串
            # self.obj.var['sp'] = self.obj.var_names.str.startswith('meaningless-')
            # 改一下,让发生这事的概率更小一点
            self.obj.var['sp'] = self.obj.var_names.str.startswith('meaninglessmeaningless-')
        else:
            # 拷贝一个对象用.copy()
            genelist = genelist.copy()
            tmp = []
            for i in self.genes:
                if i in genelist:
                    tmp.append(True)
                else:
                    tmp.append(False)
            # 这里在scanpy对象的var里面添加了一列自己对基因的标识,比如,如果是线粒体基因就会是true,不是则为false
            self.obj.var['sp'] = np.asarray(tmp)
        sc.pp.calculate_qc_metrics(self.obj, qc_vars=['sp'], percent_top=None, log1p=False, inplace=True)
        sc.pl.violin(self.obj, ['n_genes_by_counts', 'total_counts', 'pct_counts_sp'], jitter=0.4, multi_panel=True)
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
        
        # 下面执行了过滤操作
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
        
        # 展示结果环节
        sc.pl.violin(self.obj, ['n_genes_by_counts', 'total_counts', 'pct_counts_sp'], jitter=0.4, multi_panel=True)
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

    # 归一化,筛选高变基因等,大致用scanpy自带的功能完成
    def norm_log_hvg_filter_regress_scale_0(self, norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,
                                            hvg_min_dispersions=None, regresslist=None, scale_max_value=None):
        if norm_target_sum is None:
            norm_target_sum = 1e4
        # 高变基因表达量最小值,跟scanpy相同,一般不用改
        if hvg_min_mean is None:
            hvg_min_mean = 0.0125
        # 高变基因表达量最大值,跟scanpy相同,一般不用改
        if hvg_max_mean is None:
            hvg_max_mean = 3
        # 高变基因存在表达的细胞数最小值,跟scanpy相同,一般不用改
        if hvg_min_dispersions is None:
            hvg_min_dispersions = 0.5
        # 回归掉一些要素的影响,这里是总表达量
        if regresslist is None:
            regresslist = ['total_counts']
        # 标准化后的最大值为10
        if scale_max_value is None:
            scale_max_value = 10

        if norm_target_sum is not None:
            # 首先把所有细胞的表达量缩放为同一个值
            sc.pp.normalize_total(self.obj, target_sum=norm_target_sum)
            # 取log
            sc.pp.log1p(self.obj)
        if hvg_min_mean is not None and hvg_max_mean is not None and hvg_min_dispersions is not None:
            # 筛hvg
            sc.pp.highly_variable_genes(self.obj, min_mean=hvg_min_mean, max_mean=hvg_max_mean,
                                        min_disp=hvg_min_dispersions)
            sc.pl.highly_variable_genes(self.obj)
            self.obj.raw = self.obj
            self.obj = self.obj[:, self.obj.var.highly_variable]
        if regresslist is not None:
            regresslist = regresslist.copy()
            sc.pp.regress_out(self.obj, regresslist)
        if scale_max_value is not None:
            sc.pp.scale(self.obj, max_value=scale_max_value)

    # 进行PCA,筛选差异表达基因DEGs,这里markerSelect应该改成DEGSelect的,之前没注意
    def pca_cluster_markerSelect_0(self, pca_svd_algorithm=None, pca_genelist=None,
                                   cluster_n_neighbors=None, cluster_n_pcs=None, use_rep=None, perform_umap=None,
                                   umap_genelist=None, perform_paga=None,
                                   markerSelect_method=None,  # 'wilcoxon','logreg' # is ok too
                                   markerSelect_n_genes=None, show_n=None):
        # 这里的值尽量都设置得比较大了
        if pca_svd_algorithm is None:
            pca_svd_algorithm = 'arpack'
        if cluster_n_neighbors is None:
            cluster_n_neighbors = 10
        if cluster_n_pcs is None:
            cluster_n_pcs = 40
        if perform_umap is None:
            perform_umap = True
        if perform_paga is None:
            perform_paga = True
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
        
        # 用harmony来进行去批次
        if cluster_n_neighbors is not None and cluster_n_pcs is not None and use_rep == "X_pca_harmony":
            sce.pp.harmony_integrate(self.obj, 'batch')
            sc.pp.neighbors(self.obj, n_neighbors=cluster_n_neighbors, n_pcs=cluster_n_pcs, use_rep="X_pca_harmony")
        if perform_umap:
            sc.tl.umap(self.obj)
            sc.tl.leiden(self.obj)
            tmplist = ['leiden']
            if umap_genelist is not None:
                umap_genelist = umap_genelist.copy()
                for i in umap_genelist:
                    tmplist.append(i)
            sc.pl.umap(self.obj, color=tmplist)
        
        # paga的原文说,这是一种能让聚类结果变得更合理的方法
        if perform_paga:
            sc.tl.paga(self.obj)
            sc.pl.paga(self.obj)  # remove `plot=False` if you want to see the coarse-grained graph
            sc.tl.umap(self.obj, init_pos='paga')
            tmplist = ['leiden']
            if use_rep == "X_pca_harmony":
                tmplist = ['leiden', 'batch']
            if umap_genelist is not None:
                umap_genelist = umap_genelist.copy()
                for i in umap_genelist:
                    tmplist.append(i)
            for i in tmplist:
                sc.pl.umap(self.obj, color=i)
        if markerSelect_method is not None:
            sc.tl.rank_genes_groups(self.obj, 'leiden', method=markerSelect_method)
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
        
        # 这里其实就是把选出来的DEG和对应的p值存起来,方便后面回过头来看
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

    # 结果这个函数最后好像也没用到的样子
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
        
    # 好像没用到
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
        
    # 这个函数,写得就没那么美观了,而且看起来乱糟糟的,可能还有隐藏的漏洞
    # 用来生成./plant_marker_gene_list.txt中所有的基因-其他的dict
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
        
        # 这里先进行了一个过滤,如果是单个器官类型,则会只留下该器官类型的marker gene,删去其他
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
        
        # 下面是制作两个dictionary的步骤,用于DEG和已知marker比对
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
                                # 识别到'/',说明该基因,指向了至少三种不同的细胞类型
                                if '/' in dir_gene_organ_celltype[gene]:
                                    # 删除该基因先前留下的所有信息
                                    del dir_gene_organ_celltype[gene]
                                    del dir_gene_species[gene]
                                    # 将该基因拉黑
                                    forbidden_set.add(gene)
                                else:
                                    tmpinfo_1 = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                    # 如果相同的一条被添加,则跳过,不进行改动
                                    if dir_gene_organ_celltype[gene] == tmpinfo_1:
                                        pass
                                    # 只允许,相同基因,指向最多两个细胞类型
                                    else:
                                        dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + \
                                        ' / ' + 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                        dir_gene_species[gene] = species
                            else:
                                # 首次添加
                                dir_gene_organ_celltype[gene] = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                dir_gene_species[gene] = species
                                # 记录一下marker#1的基因
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
                    
                    # 这里是直接写成marker#1和marker#2不会相互影响了
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
            # 这里就是只弄marker#1的
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
        # 不变默认情况的更新,这是marker#1的权重,只在marker#1和marker#2混用时,
        # 且一个基因是marker#1基因且只对应一个细胞类型时起作用
        if marker_grade_1_coefficient is None:
            marker_grade_1_coefficient = 1
        
        # 获取关于marker基因的两个dictionary
        self.dir_gene_species, self.dir_gene_organ_celltype = self.make_marker_gene_list_0(marker_select_flag = '#1 and #2',
                                                         marker_gene_list = marker_gene_list)
        dir_gene_species, dir_gene_organ_celltype = self.make_marker_gene_list_0(marker_select_flag = '#1 and #2',
                                                         marker_gene_list = marker_gene_list)
        count_1 = 0
        count_2 = 0
        anno_re = []

        # 这里不知道为什么当初选用了self.markergenepd_head20作为一个暂时的载体,虽然后面改回去了
        tmpmarkergenepd = self.markergenepd_head20.copy()
        if self.markergenepd_head is not None:
            self.markergenepd_head20 = self.markergenepd_head
        for i in range(0, int(self.markergenepd_head20.shape[1] / 2)):
            tmp = {}
            # 根据当初DEG的表格的组织形式,针对每个cluster,独立地分别地提取DEG,然后一一跟已知marker表比对
            for j in self.markergenepd_head20[str(i) + '_n'].values:
                j = str(j).strip()
                # print(j)
                # break
                if j in dir_gene_species:
                    # 这里存在可能不合理的地方,如果是有'/'的话,最后细胞类型的呈现可能也是有'/'的
                    # 我得把这里改一下
                    # 改之后,如果一个基因是A/B,则A+1分,B+1分,(A/B)+1分
                    # 虽然还是有不合理的地方,不过先这样好了
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
                    # 这里存在可能不合理的地方,如果是有'/'的话,最后细胞类型的呈现可能也是有'/'的
                    # 我得把这里改一下
                    # 改之后,如果一个基因是A/B,则A+1分,B+1分,(A/B)+1分
                    # 虽然还是有不合理的地方,不过先这样好了
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
        # 这里就是把注释结果添加进去
        self.obj.obs['anno_re_using_rank_1and2_markers'] = tmplist.copy()

        for i, j in enumerate(anno_re_1and2):
            print('cluster ' + str(i) + ': ' + str(j))
        sc.pl.umap(self.obj, color='anno_re_using_rank_1and2_markers')

        tmplist = []
        for i in self.obj.obs['leiden']:
            tmplist.append(self.anno_re_using_rank_1_markers[int(i)])
        # print(tmplist)
        # print(len(self.anno_re_using_rank_1_markers))
        # 这里就是把注释结果添加进去
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

        # 这里不知道为什么当初选用了self.markergenepd_head20作为一个暂时的载体,虽然后面改回去了
        tmpmarkergenepd = self.markergenepd_head20.copy()
        if self.markergenepd_head is not None:
            self.markergenepd_head20 = self.markergenepd_head
        for i in range(0, int(self.markergenepd_head20.shape[1] / 2)):
            tmp = {}
            # 根据当初DEG的表格的组织形式,针对每个cluster,独立地分别地提取DEG,然后一一跟已知marker表比对
            for j in self.markergenepd_head20[str(i) + '_n'].values:
                j = str(j).strip()
                # print(j)
                # break
                if j in self.dir_gene_species and j in self.dir_gene_score and float(self.dir_gene_score[j])>=the_least_score:
                    # 这里存在可能不合理的地方,如果是有'/'的话,最后细胞类型的呈现可能也是有'/'的
                    # 我得把这里改一下
                    # 改之后,如果一个基因是A/B,则A+1分,B+1分,(A/B)+1分
                    # 虽然还是有不合理的地方,不过先这样好了
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
        # 这里就是把注释结果添加进去
        self.obj.obs['anno_re_using_marker_scores'] = tmplist.copy()

        for i, j in enumerate(anno_re):
            print('cluster ' + str(i) + ': ' + str(j))
        sc.pl.umap(self.obj, color='anno_re_using_marker_scores')

        self.markergenepd_head20 = tmpmarkergenepd.copy()
        del tmpmarkergenepd
        self.memory_collect_0()

    # 整合后的函数
    def to_clusters_0(self, verbosity=None, dpi=None, facecolor=None, # set_default_0
            anndata=None,filepath=None, species=None, organ=None, organs=None, # data_input_0
            n_top=None, num_min_genes=None, num_max_genes=None, num_min_cells=None, # filter_0
                max_percent_of_specific_genes=None,
                genelist=None, mad_coeffcient=None,
            norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,  # norm_log_hvg_filter_regress_scale_0
                hvg_min_dispersions=None, regresslist=None, scale_max_value=None,
            pca_svd_algorithm=None, pca_genelist=None, # pca_cluster_markerSelect_0
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
            hvg_max_mean=hvg_max_mean, hvg_min_dispersions=hvg_min_dispersions, 
            regresslist=regresslist, scale_max_value=scale_max_value)
        self.pca_cluster_markerSelect_0(pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist,
            cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, use_rep=use_rep, 
            perform_umap=perform_umap,umap_genelist=umap_genelist, perform_paga=perform_paga,
            markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
            markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)

    # 整合后的函数
    def batch_to_clusters_0(self, verbosity=None, dpi=None, facecolor=None, # set_default_0
            anndatas=None,filepaths=None, species=None, organs=None, # batch_data_input_0
            n_top=None, num_min_genes=None, num_max_genes=None, num_min_cells=None, # filter_0
                max_percent_of_specific_genes=None,
                genelist=None, mad_coeffcient=None,
            norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,  # norm_log_hvg_filter_regress_scale_0
                hvg_min_dispersions=None, regresslist=None, scale_max_value=None,
            pca_svd_algorithm=None, pca_genelist=None, # pca_cluster_markerSelect_0
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
            hvg_max_mean=hvg_max_mean, hvg_min_dispersions=hvg_min_dispersions, 
            regresslist=regresslist, scale_max_value=scale_max_value)
        self.pca_cluster_markerSelect_0(pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist,
            cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, use_rep=use_rep, 
            perform_umap=perform_umap,umap_genelist=umap_genelist, perform_paga=perform_paga,
            markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
            markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)

    # 整合后的函数
    def to_annotation_0(self, verbosity=None, dpi=None, facecolor=None, # set_default_0
            species=None,anndata=None,filepath=None,organ=None,anndatas=None,filepaths=None,organs=None, 
                        # batch_data_input_0
            n_top=None, num_min_genes=None, num_max_genes=None, num_min_cells=None, # filter_0
                max_percent_of_specific_genes=None,
                genelist=None, mad_coeffcient=None,
            norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,  # norm_log_hvg_filter_regress_scale_0
                hvg_min_dispersions=None, regresslist=None, scale_max_value=None,
            pca_svd_algorithm=None, pca_genelist=None, # pca_cluster_markerSelect_0
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
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
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
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
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
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
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
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
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
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
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

    # 暂时用不到
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

    # 暂时用不到
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
        
    # 输入参考集合待注释集,自动训练模型
    def auto_anno_mldl_train_0(self,ref_obj=None,test_obj=None,gpu_code = None,
                               learning_rate = None,epochs = None,
                               batch_size = None,dropout = None,
                               last_linear_n = None,num_workers = None,
                               model_weight_save_path = None,
                               mapping_1_save_path = None,
                               mapping_2_save_path = None,
                               num_types_save_path = None):
        # 要在哪块GPU上跑
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
        
        df,mapping_1,mapping_2 = self.auto_annotation_with_deep_learning_transfer_adata_to_df_0(adata_input = ref_obj)
        
        # 显示类别数和每个类别的数量
        category_counts = df['celltype'].value_counts()
        print("类别数：", len(category_counts))
        print("每个类别的数量：")
        print(category_counts)
        print(mapping_1)
        print(mapping_2)
        
        # 统计每个细胞类型的数目
        category_counts = df['celltype'].value_counts()
        # 划分数据集
        # 创建空数据框
        train_data = pd.DataFrame()
        valid_data = pd.DataFrame()
        test_data = pd.DataFrame()
        # 对每个细胞类型随机采样1000个细胞
        for category in category_counts.index:
            category_data = df[df['celltype'] == category]
            category_data = category_data.sample(n=1000, random_state=1, replace=True)
            # from sklearn.model_selection import train_test_split
            train, temp = train_test_split(category_data, train_size=0.8, random_state=1)
            valid, test = train_test_split(temp, train_size=0.5, random_state=1)
            train_data = pd.concat([train_data, train])
            valid_data = pd.concat([valid_data, valid])
            test_data = pd.concat([test_data, test])
            
        train = train_data
        valid = valid_data
        test = test_data
        
        # 分离注释和数据
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
        # 这里继承了Dataset,封装、继承、多态,面向对象编程三大特性,可以去了解一下
        # 这里应该算是固定招式了,也拿别人定好的模板来用,实现初始化、根据索引拿数据、求长度,三个功能
        class trainDataset(Dataset):  # 这是一个Dataset子类
            def __init__(self):
                # 这里的数据转换过程涉及到数据从cpu转移到gpu上的准备工作
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
            
        class validDataset(Dataset):  # 这是一个Dataset子类
            def __init__(self):
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
            
        traindata = trainDataset()
        validdata = validDataset()
        # 随机打乱,设置batch,设置取的时候的最大cpu线/进程数(?有点不知道这里怎么形容)
        train_sampler = RandomSampler(traindata)
        train_iter = DataLoader(traindata,batch_size=the_batch_size,shuffle=False,sampler=train_sampler,
                                num_workers=the_num_workers)
        valid_sampler = RandomSampler(validdata)
        valid_iter = DataLoader(validdata,batch_size=the_batch_size,shuffle=False,sampler=valid_sampler,
                                num_workers=the_num_workers)
        
        # 获取细胞类型总数
        num_types = train_label.nunique().values[0]

        import torchvision.models
        
        class PscResNet(nn.Module):
            def __init__(self):
                super(PscResNet, self).__init__()
                # change "pretrained=False" to "weights=None" may be better,"pretrained=False"比较旧了
                self.resnet18 = torchvision.models.resnet18(pretrained=False)   
                self.pscMLP = nn.Sequential(
                                nn.ReLU(),#隐藏层激活函数
                                nn.Dropout(the_dropout),
                                nn.Linear(the_last_linear_n, num_types))#输出层
                
            def forward(self, x):
                x = self.resnet18(x)
                x = self.pscMLP(x)
                return x
        
        pscResNet1 = PscResNet()
        # 加载预训练模型的权重
        pretrained_dict = torch.load('resnet18.pth')
        model_dict = pscResNet1.resnet18.state_dict()
        
        # 从预训练模型的权重中提取对应的部分
        pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
        
        # 更新 self.resnet18 的权重
        model_dict.update(pretrained_dict)
        pscResNet1.resnet18.load_state_dict(model_dict)
        
        import torch
        torch.cuda.init()
                    
        # 使用GPU计算模型在数据集上的精度
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
        
        # 用GPU训练模型
        def train_with_GPU(net, train_iter, valid_iter, num_epochs, lr, device):
            # 初始化权重
            def init_weights(m):
                if type(m) == nn.Linear or type(m) == nn.Conv2d:
                    nn.init.xavier_uniform_(m.weight)
            net.apply(init_weights)
            print('training on', device)
            # 移到gpu上
            net.to(device)
            # 优化器,momentum越大越倾向于保持旧的移动方向
            optimizer = torch.optim.SGD(net.parameters(), lr=lr, momentum=0.9)
            # 学习率变动器
            scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=150, eta_min=0.00001)
            loss = nn.CrossEntropyLoss()
            # 动画效果
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
        # 训练
        train_with_GPU(pscResNet1, train_iter, valid_iter, num_epochs, lr, d2l.try_all_gpus()[gpu_code])
        # 保存一下权重、mapping_1和mapping_2
        self.model_weight = pscResNet1.state_dict()
        self.mapping_1 = mapping_1
        self.mapping_2 = mapping_2
        self.num_types = num_types
        if model_weight_save_path is not None and \
                mapping_1_save_path is not None and \
                mapping_2_save_path is not None and \
                num_types_save_path is not None:
            torch.save(pscResNet1.state_dict(), model_weight_save_path)
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
            torch.save(pscResNet1.state_dict(), 'model_weight_'+current_time.replace(':','-').replace(' ','-')+'.pth')
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
        
        # 将测试数据和真实标签转换为 Tensor 对象
        tmp_shape = test_data.shape[0]
        test_data = torch.from_numpy(np.float32(np.array(test_data)))
        test_data = torch.reshape(test_data, (tmp_shape, 224, 224))
        test_data = torch.unsqueeze(test_data, 1)
        test_data = torch.cat((test_data, test_data, test_data), dim=1)
        test_label = torch.from_numpy(np.int64(np.array(test_label)).reshape(1,-1)[0])
        
        pscResNet1 = pscResNet1.to('cpu')
        pscResNet1.eval()
        
        # 关闭梯度跟踪
        with torch.no_grad():
            outputs = pscResNet1(test_data)
        
        # 获取预测结果
        _, predicted = torch.max(outputs, 1)
        
        from sklearn.metrics import accuracy_score, precision_score, recall_score
        
        # 将预测结果和真实标签转换为 NumPy 数组
        predicted = predicted.numpy()
        test_label = test_label.numpy()
        
        # 计算准确率
        accuracy = accuracy_score(test_label, predicted)
        
        # 计算精确度
        precision = precision_score(test_label, predicted, average='macro', zero_division=1)
        
        # 计算召回率
        recall = recall_score(test_label, predicted, average='macro', zero_division=1)
        
        print('accuracy')
        print(accuracy)
        print('precision')
        print(precision)
        print('recall')
        print(recall)
        
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

        class PscResNet(nn.Module):
            def __init__(self):
                super(PscResNet, self).__init__()
                self.resnet18 = torchvision.models.resnet18(pretrained=False)   # change "pretrained=False" to "weights=None" may be better
                self.pscMLP = nn.Sequential(
                                nn.ReLU(),#隐藏层激活函数
                                nn.Dropout(the_dropout),
                                nn.Linear(1000, aaa_num_types))#输出层
                
            def forward(self, x):
                x = self.resnet18(x)
                x = self.pscMLP(x)
                return x
        
        pscResNet1 = PscResNet()
        # 加载预训练模型的权重
        pretrained_dict = self.model_weight
        model_dict = pscResNet1.state_dict()
        
        # 从预训练模型的权重中提取对应的部分
        pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
        
        # 更新 self.resnet18 的权重
        model_dict.update(pretrained_dict)
        pscResNet1.load_state_dict(model_dict)
        
        df,_,_ = self.auto_annotation_with_deep_learning_transfer_adata_to_df_0(adata_input = test_obj)
        
        df_label = df.iloc[:, -1:]
        df_label_full = df.iloc[:, -2:]
        df = df.iloc[:, :-2]
        
        import math
        
        batch_size = 32 # 每个批次的大小
        total_size = df.shape[0]
        total_batches = math.ceil(total_size / batch_size)
        
        net_try = pscResNet1
        net_try = net_try.to('cpu')
        net_try.eval()
        
        # 关闭梯度跟踪
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
        
        # 创建逆映射的字典
        reverse_mapping_2 = {v: k for k, v in self.mapping_2.items()}
        
        # 还原 list1 中的值为 list2_0 中的映射后的值
        list2_0_mapped = [reverse_mapping_2[value] for value in list2]
        
        # 还原 list2_0 中的映射后的值为原始值
        list2_0 = [next(key for key, value in self.mapping_1.items() if value == mapped_value) for mapped_value in list2_0_mapped]
        
        # 打印还原后的 list2_0
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
                                           gpu_code = None,
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
            self.auto_anno_mldl_train_0(ref_obj=ref_obj,test_obj=test_obj,gpu_code = gpu_code,
                                        learning_rate = learning_rate,epochs = epochs,
                                        batch_size = batch_size,dropout = dropout,
                                        last_linear_n = last_linear_n,num_workers = num_workers,
                                        model_weight_save_path = model_weight_save_path,
                                        mapping_1_save_path = mapping_1_save_path,
                                        mapping_2_save_path = mapping_2_save_path,
                                        num_types_save_path = num_types_save_path)
        if run_model == 'train_and_predict':
            self.auto_anno_mldl_train_0(ref_obj=ref_obj,test_obj=test_obj,gpu_code = gpu_code,
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
        
    def auto_annotation_with_deep_learning_transfer_adata_to_df_0(self,adata_input = None):
        # 这里就进行一个normalize_total的处理
        # sc.pp.normalize_total(adata_input, target_sum=1e4)
        adata_input.obs['cellname'] = adata_input.obs_names
        
        # 将scanpy对象的表达矩阵,转成数据框,去除重复细胞名的细胞
        tmp_pcm = adata_input
        data=tmp_pcm.X.todense()
        pcm_pd=pd.DataFrame(data,index=tmp_pcm.obs_names,columns=tmp_pcm.var_names)
        pcm_pd.loc[:,'cellname']=pcm_pd.index
        df = pcm_pd
        df.drop_duplicates(subset='cellname', keep='first', inplace=True)
        
        duplicate_rows = df.duplicated(subset='cellname', keep=False)
        # 过滤掉重复行
        df = df[~duplicate_rows]
        pcm_pd = df
        # print('pcm')
        # print(pcm_pd)
        # print('\n')
        
        label=pd.concat([adata_input.obs['cellname'], adata_input.obs['celltype']], axis=1)
        df = label
        # 获取"celltype"中值的类别数
        num_categories = df['celltype'].nunique()
        # 构建映射字典
        mapping_1 = {category: i for i, category in enumerate(df['celltype'].unique())}
        # print('mapping_1')
        # print(mapping)
        # print('\n')
        # 使用映射字典将"celltype"中的每个值转换为对应的数字
        df['celltype'] = df['celltype'].map(mapping_1)
        # 根据"cellname"列判断重复行
        duplicate_rows = df.duplicated(subset='cellname', keep=False)
        # 过滤掉重复行
        df = df[~duplicate_rows]
        # print('df')
        # print(df)
        # print('\n')
        label = df
        
        # 保留交集部分
        merged = pd.merge(pcm_pd, label, on='cellname', how='inner')
        # 打印结果
        # print("merged")
        # print(merged)
        df = merged
        # 获取"celltype"中值的类别数
        num_categories = df['celltype'].nunique()
        # 构建映射字典
        mapping_2 = {category: i for i, category in enumerate(df['celltype'].unique())}
        # print('mapping_2')
        # print(mapping)
        # print('\n')
        # 使用映射字典将"celltype"中的每个值转换为对应的数字
        df['celltype'] = df['celltype'].map(mapping_2)
        
        # 创建一个全零的 DataFrame,行数与原始 DataFrame 相同,列数为要扩充的目标列数
        # 这里是默认使用ResNet了,要更改数据形状为224*224=50176
        # 这里+2是因为celltype和cellname列之后要去掉
        # 因为别人的ResNet训练了可能几千几万个epoch,所以我们拿别人的权重来用,比自己从头开始训练,收敛速度会快很多
        # 新添加的行因为没有意义就全填0,后面会变成0.0001
        # 这里使全局每个值都不是0是考虑到一个说法,数据为0不利于更新参数
        target_columns = 50176-df.shape[1]+2
        zeros_df = pd.DataFrame(np.zeros((df.shape[0], target_columns)), columns=[f'new_col_{i}' for i in range(target_columns)])
        
        # 将原始 DataFrame 与全零 DataFrame 进行列拼接
        expanded_df = pd.concat([df, zeros_df], axis=1)
        
        # 扩充后的 DataFrame 大小
        # print(expanded_df.shape)  # 输出: (5020, 50176)
        
        df = expanded_df
        # print(df)
        # 获取列名列表
        columns = df.columns.tolist()
        
        # 移动 'cellname' 列到倒数第一列
        columns.remove('cellname')
        columns.append('cellname')
        
        # 再移动 'celltype' 列到倒数第一列
        columns.remove('celltype')
        columns.append('celltype')
        
        # 重新排列列顺序
        df = df[columns]
        
        return df,mapping_1,mapping_2