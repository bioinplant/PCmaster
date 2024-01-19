import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import gc
import os
import psutil

from IPython.core.interactiveshell import InteractiveShell

InteractiveShell.ast_node_interactivity = "all"

print("Memory collecting...")
gc.collect()
print("Memory information: ")
info = psutil.virtual_memory()
print("Used: ")
print(psutil.Process(os.getpid()).memory_info().rss)
print("Total: ")
print(info.total)
print("Used (%): ")
print(info.percent)


class PCmaster_anno_0_Error(Exception):
    pass


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
    anno_re_using_rank_1_markers = None
    anno_re_using_rank_1and2_markers = None

    def memory_collect_0(self):
        print("Memory collecting...")
        gc.collect()
        print("Memory information: ")
        info = psutil.virtual_memory()
        print("Used: ")
        print(psutil.Process(os.getpid()).memory_info().rss)
        print("Total: ")
        print(info.total)
        print("Used (%): ")
        print(info.percent)

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

    def data_input_0(self, filepath=None, species=None, organ=None):
        if filepath is not None:
            if filepath[-6:] == "tsv.gz":
                self.obj = sc.read_umi_tools(filepath)
            if filepath[-6:] == "mtx.gz":
                self.obj = sc.read_10x_mtx(filepath[:-13], var_names='gene_symbols', cache=True)
            tmpstr1 = 'filtered_feature_bc_matrix'
            if tmpstr1 in filepath:
                self.obj = sc.read_10x_mtx(filepath, var_names='gene_symbols', cache=True)
                
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

    def data_input_0_from_anndata(self, anndata=None, species=None, organ=None):
        if anndata is not None:
            self.obj = anndata
        else:
            raise PCmaster_anno_0_Error("The anndata is None.")
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
            batch_categories.append("batch_" + str(i))
        if self.obj is not None:
            print("Warning: self.obj is not None and will be changed.")
        self.obj = self.objs[0].concatenate(self.objs[1:], batch_categories=batch_categories)

    def batch_data_input_0_and_concat_0(self, filepaths=None, species=None, organs=None):
        pass

    def get_dense_0(self):
        if self.obj is None:
            raise PCmaster_anno_0_Error("The expression matrix is not successfully loaded.")
        self.dd = pd.DataFrame(self.obj.X.todense(), index=self.obj.obs_names, columns=self.obj.var_names)

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
            self.obj.var['sp'] = self.obj.var_names.str.startswith('meaningless-')
        else:
            genelist = genelist.copy()
            tmp = []
            for i in self.genes:
                if i in genelist:
                    tmp.append(True)
                else:
                    tmp.append(False)
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

    def norm_log_hvg_filter_regress_scale_0(self, norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,
                                            hvg_min_dispersions=None, regresslist=None, scale_max_value=None):
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

    def pca_cluster_markerSelect_0(self, pca_svd_algorithm=None, pca_genelist=None,
                                   cluster_n_neighbors=None, cluster_n_pcs=None, use_rep=None, perform_umap=None,
                                   umap_genelist=None, perform_paga=None,
                                   markerSelect_method=None,  # 'wilcoxon','logreg' # is ok too
                                   markerSelect_n_genes=None, show_n=None):
        if pca_svd_algorithm is None:
            pca_svd_algorithm = 'arpack'
        if cluster_n_neighbors is None:
            cluster_n_neighbors = 20
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
            show_n = 20

        if pca_svd_algorithm is not None:
            sc.tl.pca(self.obj, svd_solver=pca_svd_algorithm)
            if pca_genelist is not None:
                pca_genelist = pca_genelist.copy()
                sc.pl.pca(self.obj, color=pca_genelist)
            sc.pl.pca_variance_ratio(self.obj, log=True)
        if cluster_n_neighbors is not None and cluster_n_pcs is not None and use_rep is None:
            sc.pp.neighbors(self.obj, n_neighbors=cluster_n_neighbors, n_pcs=cluster_n_pcs)
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
        #     for group in self.obj.uns['rank_genes_groups']['names'].dtype.names for key in ['names', 'pvals']}).head(show_n) )
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
        
    def make_marker_gene_list_0(self,marker_select_flag = None):
        # organ = None
        # batch_organs = []
        if (marker_select_flag is None) or (marker_select_flag != '#1'):
            marker_select_flag = '#1 and #2'
        ref_txt = open('./plant_marker_gene_list.txt', 'r', encoding='utf-8')
        ref1 = []
        for line in ref_txt:
            if self.organ is not None:
                if self.organ.upper() in line.upper() or self.organ.upper() == 'ALL':
                    ref1.append(line)
            elif len(self.batch_organs) != 0:
                # tmp = self.batch_organs[0]
                # allsame = 'yes'
                # for i in self.batch_organs:
                #     if i != tmp:
                #         allsame = 'no'
                # if (tmp.upper()!='ALL' and allsame == 'yes' and self.batch_organs[0].upper() in line.upper()) or (allsame == 'no') or (tmp.upper()=='ALL'):
                #     ref1.append(line)
                for i in self.batch_organs:
                    if i.upper() in line.upper() or i.upper() == 'ALL':
                        ref1.append(line)
            else:
                ref1.append(line)
        
        dir_gene_species = {}
        dir_gene_organ_celltype = {}
        count = 0
        if marker_select_flag == '#1 and #2':
            forbidden_set = set()
            set_rank1 = set()
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
                                    tmpinfo_2 = ''
                                    if self.organ is not None:
                                        tmpinfo_2 = 'organ: ' + self.organ + '  ' + 'celltype: ' + celltype
                                    if self.batch_organs is not None and len(self.batch_organs)==1:
                                        tmpinfo_2 = 'organ: ' + self.batch_organs[0] + '  ' + 'celltype: ' + celltype
                                    if dir_gene_organ_celltype[gene] == tmpinfo_1 or dir_gene_organ_celltype[gene] == tmpinfo_2:
                                        pass
                                    else:
                                        if self.organ is not None and self.organ.upper() in organ.upper():
                                            dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + ' / ' + 'organ: ' + self.organ + '  ' + 'celltype: ' + celltype
                                            dir_gene_species[gene] = species
                                        elif self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper() in organ.upper():
                                            dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + ' / ' + 'organ: ' + self.batch_organs[0] + '  ' + 'celltype: ' + celltype
                                            dir_gene_species[gene] = species
                                        elif (self.organ is not None and self.organ.upper()=='ALL') or (self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper()=='ALL') or ( self.batch_organs is not None and len(self.batch_organs)>1 ):
                                            dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + ' / ' + 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                            dir_gene_species[gene] = species
                                        else:
                                            pass
                            else:
                                if self.organ is not None and self.organ.upper() in organ.upper():
                                    dir_gene_organ_celltype[gene] = 'organ: ' + self.organ + '  ' + 'celltype: ' + celltype
                                    dir_gene_species[gene] = species
                                elif self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper() in organ.upper():
                                    dir_gene_organ_celltype[gene] = 'organ: ' + self.batch_organs[0] + '  ' + 'celltype: ' + celltype
                                    dir_gene_species[gene] = species
                                elif (self.organ is not None and self.organ.upper()=='ALL') or (self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper()=='ALL') or ( self.batch_organs is not None and len(self.batch_organs)>1 ):
                                    dir_gene_organ_celltype[gene] = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                    dir_gene_species[gene] = species
                                else:
                                    pass
                        # print(organ)
                        # print(celltype)
                        # print(gene)
                count = count + 1
            for line in ref1:
                # print(line.split('\t'))
                if len(line.split('\t')) < 5:
                    print(line)
                if len(line.split('\t')) >= 5 and ( line.split('\t')[4] == 'Marker #2\n' or line.split('\t')[4] == 'Cell marker #2\n' ):
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
                    
                    if (gene not in forbidden_set) and (gene not in set_rank1):
                        if gene in dir_gene_organ_celltype:
                            if '/' in dir_gene_organ_celltype[gene]:
                                del dir_gene_organ_celltype[gene]
                                del dir_gene_species[gene]
                                forbidden_set.add(gene)
                            else:
                                tmpinfo_1 = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                tmpinfo_2 = ''
                                if self.organ is not None:
                                    tmpinfo_2 = 'organ: ' + self.organ + '  ' + 'celltype: ' + celltype
                                if self.batch_organs is not None and len(self.batch_organs)==1:
                                    tmpinfo_2 = 'organ: ' + self.batch_organs[0] + '  ' + 'celltype: ' + celltype
                                if dir_gene_organ_celltype[gene] == tmpinfo_1 or dir_gene_organ_celltype[gene] == tmpinfo_2:
                                    pass
                                else:
                                    if self.organ is not None and self.organ.upper() in organ.upper():
                                        dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + ' / ' + 'organ: ' + self.organ + '  ' + 'celltype: ' + celltype
                                        dir_gene_species[gene] = species
                                    elif self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper() in organ.upper():
                                        dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + ' / ' + 'organ: ' + self.batch_organs[0] + '  ' + 'celltype: ' + celltype
                                        dir_gene_species[gene] = species
                                    elif (self.organ is not None and self.organ.upper()=='ALL') or (self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper()=='ALL') or ( self.batch_organs is not None and len(self.batch_organs)>1 ):
                                        dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + ' / ' + 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                        dir_gene_species[gene] = species
                                    else:
                                        pass
                        else:
                            if self.organ is not None and self.organ.upper() in organ.upper():
                                dir_gene_organ_celltype[gene] = 'organ: ' + self.organ + '  ' + 'celltype: ' + celltype
                                dir_gene_species[gene] = species
                            elif self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper() in organ.upper():
                                dir_gene_organ_celltype[gene] = 'organ: ' + self.batch_organs[0] + '  ' + 'celltype: ' + celltype
                                dir_gene_species[gene] = species
                            elif (self.organ is not None and self.organ.upper()=='ALL') or (self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper()=='ALL') or ( self.batch_organs is not None and len(self.batch_organs)>1 ):
                                dir_gene_organ_celltype[gene] = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                dir_gene_species[gene] = species
                            else:
                                pass
                    # print(organ)
                    # print(celltype)
                    # print(gene)
                count = count + 1
        else:
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
                                    tmpinfo_2 = ''
                                    if self.organ is not None:
                                        tmpinfo_2 = 'organ: ' + self.organ + '  ' + 'celltype: ' + celltype
                                    if self.batch_organs is not None and len(self.batch_organs)==1:
                                        tmpinfo_2 = 'organ: ' + self.batch_organs[0] + '  ' + 'celltype: ' + celltype
                                    if dir_gene_organ_celltype[gene] == tmpinfo_1 or dir_gene_organ_celltype[gene] == tmpinfo_2:
                                        pass
                                    else:
                                        if self.organ is not None and self.organ.upper() in organ.upper():
                                            dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + ' / ' + 'organ: ' + self.organ + '  ' + 'celltype: ' + celltype
                                            dir_gene_species[gene] = species
                                        elif self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper() in organ.upper():
                                            dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + ' / ' + 'organ: ' + self.batch_organs[0] + '  ' + 'celltype: ' + celltype
                                            dir_gene_species[gene] = species
                                        elif (self.organ is not None and self.organ.upper()=='ALL') or (self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper()=='ALL') or ( self.batch_organs is not None and len(self.batch_organs)>1 ):
                                            dir_gene_organ_celltype[gene] = dir_gene_organ_celltype[gene] + ' / ' + 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                            dir_gene_species[gene] = species
                                        else:
                                            pass
                            else:
                                if self.organ is not None and self.organ.upper() in organ.upper():
                                    dir_gene_organ_celltype[gene] = 'organ: ' + self.organ + '  ' + 'celltype: ' + celltype
                                    dir_gene_species[gene] = species
                                elif self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper() in organ.upper():
                                    dir_gene_organ_celltype[gene] = 'organ: ' + self.batch_organs[0] + '  ' + 'celltype: ' + celltype
                                    dir_gene_species[gene] = species
                                elif (self.organ is not None and self.organ.upper()=='ALL') or (self.batch_organs is not None and len(self.batch_organs)==1 and self.batch_organs[0].upper()=='ALL') or ( self.batch_organs is not None and len(self.batch_organs)>1 ):
                                    dir_gene_organ_celltype[gene] = 'organ: ' + organ + '  ' + 'celltype: ' + celltype
                                    dir_gene_species[gene] = species
                                else:
                                    pass
                        # print(organ)
                        # print(celltype)
                        # print(gene)
                count = count + 1
        return dir_gene_species, dir_gene_organ_celltype

    def annotation_with_marker_gene_list_0(self,marker_select_flag = None):
        dir_gene_species, dir_gene_organ_celltype = self.make_marker_gene_list_0(marker_select_flag = '#1 and #2')
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
                if j in dir_gene_species:
                    organ_celltype = dir_gene_organ_celltype[j]
                    if organ_celltype in tmp:
                        tmp[organ_celltype] = tmp[organ_celltype] + 1
                    else:
                        tmp[organ_celltype] = 1
            if tmp != {}:
                anno_re.append(max(tmp, key=tmp.get))
            else:
                anno_re.append('Unknown')
        self.anno_re_using_rank_1and2_markers = anno_re

        # for i in anno_re:
        #     print(i)
        # print(len(anno_re))
        anno_re_1and2 = anno_re

        dir_gene_species, dir_gene_organ_celltype = self.make_marker_gene_list_0(marker_select_flag = '#1')

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
                    if organ_celltype in tmp:
                        tmp[organ_celltype] = tmp[organ_celltype] + 1
                    else:
                        tmp[organ_celltype] = 1
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
        self.obj.obs['anno_re_using_rank_1and2_markers'] = tmplist.copy()

        for i, j in enumerate(anno_re_1and2):
            print('cluster ' + str(i) + ': ' + str(j))
        sc.pl.umap(self.obj, color='anno_re_using_rank_1and2_markers')

        tmplist = []
        for i in self.obj.obs['leiden']:
            tmplist.append(self.anno_re_using_rank_1_markers[int(i)])
        # print(tmplist)
        # print(len(self.anno_re_using_rank_1_markers))
        self.obj.obs['anno_re_using_rank_1_markers'] = tmplist.copy()
        for i, j in enumerate(anno_re_1):
            print('cluster ' + str(i) + ': ' + str(j))
        sc.pl.umap(self.obj, color='anno_re_using_rank_1_markers')

        self.markergenepd_head20 = tmpmarkergenepd.copy()
        del tmpmarkergenepd
        self.memory_collect_0()

    def to_clusters_0(self, verbosity=None, dpi=None, facecolor=None, # set_default_0
            anndata=None,filepath=None, species=None, organ=None, # batch_data_input_0
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
            self.data_input_0(filepath=filepath, species=species, organ=organ)
        if anndata is not None:
            self.data_input_0_from_anndata(anndata=anndata, species=species, organ=organ)
        self.filter_0()
        self.norm_log_hvg_filter_regress_scale_0()
        self.pca_cluster_markerSelect_0()

    def batch_to_clusters_0(self, verbosity=None, dpi=None, facecolor=None, # set_default_0
            anndatas=None,filepaths=None, species=None, organs=None, # batch_data_input_0
            n_top=None, num_min_genes=None, num_max_genes=None, num_min_cells=None, # filter_0
                max_percent_of_specific_genes=None,
                genelist=None, mad_coeffcient=None,
            norm_target_sum=None, hvg_min_mean=None, hvg_max_mean=None,  # norm_log_hvg_filter_regress_scale_0
                hvg_min_dispersions=None, regresslist=None, scale_max_value=None,
            pca_svd_algorithm=None, pca_genelist=None, # pca_cluster_markerSelect_0
                cluster_n_neighbors=None, cluster_n_pcs=None, use_rep="X_pca_harmony", perform_umap=None,
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
            cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, use_rep="X_pca_harmony", 
            perform_umap=perform_umap,umap_genelist=umap_genelist, perform_paga=perform_paga,
            markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
            markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)

    def to_annotation_0(self, verbosity=None, dpi=None, facecolor=None, # set_default_0
            species=None,anndata=None,filepath=None,organ=None,anndatas=None,filepaths=None,organs=None, # batch_data_input_0
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
            marker_select_flag = None):
        if filepath is not None:
            self.to_clusters_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor, # set_default_0
                anndata=anndata,filepath=filepath, species=species, organ=organ, # batch_data_input_0
                n_top=n_top, num_min_genes=num_min_genes, num_max_genes=num_max_genes, # filter_0
                    num_min_cells=num_min_cells, 
                    max_percent_of_specific_genes=max_percent_of_specific_genes,
                    genelist=genelist, mad_coeffcient=mad_coeffcient,
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  # norm_log_hvg_filter_regress_scale_0
                    hvg_min_dispersions=hvg_min_dispersions, regresslist=regresslist, scale_max_value=scale_max_value,
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
                    cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, 
                    use_rep=use_rep, perform_umap=perform_umap,
                    umap_genelist=umap_genelist, perform_paga=perform_paga,
                    markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
                    markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)
        elif anndata is not None:
            self.to_clusters_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor, # set_default_0
                anndata=anndata,filepath=filepath, species=species, organ=organ, # batch_data_input_0
                n_top=n_top, num_min_genes=num_min_genes, num_max_genes=num_max_genes, # filter_0
                    num_min_cells=num_min_cells, 
                    max_percent_of_specific_genes=max_percent_of_specific_genes,
                    genelist=genelist, mad_coeffcient=mad_coeffcient,
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  # norm_log_hvg_filter_regress_scale_0
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
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  # norm_log_hvg_filter_regress_scale_0
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
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  # norm_log_hvg_filter_regress_scale_0
                    hvg_min_dispersions=hvg_min_dispersions, regresslist=regresslist, scale_max_value=scale_max_value,
                pca_svd_algorithm=pca_svd_algorithm, pca_genelist=pca_genelist, # pca_cluster_markerSelect_0
                    cluster_n_neighbors=cluster_n_neighbors, cluster_n_pcs=cluster_n_pcs, 
                    use_rep=use_rep, perform_umap=perform_umap,
                    umap_genelist=umap_genelist, perform_paga=perform_paga,
                    markerSelect_method=markerSelect_method,  # 'wilcoxon','logreg' # is ok too
                    markerSelect_n_genes=markerSelect_n_genes, show_n=show_n)
        elif self.obj is not None:
            self.to_clusters_0(verbosity=verbosity, dpi=dpi, facecolor=facecolor, # set_default_0
                filepath=filepath, species=species, organ=organ, # batch_data_input_0
                n_top=n_top, num_min_genes=num_min_genes, num_max_genes=num_max_genes, # filter_0
                    num_min_cells=num_min_cells, 
                    max_percent_of_specific_genes=max_percent_of_specific_genes,
                    genelist=genelist, mad_coeffcient=mad_coeffcient,
                norm_target_sum=norm_target_sum, hvg_min_mean=hvg_min_mean, hvg_max_mean=hvg_max_mean,  # norm_log_hvg_filter_regress_scale_0
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
            self.annotation_with_marker_gene_list_0(marker_select_flag = marker_select_flag)

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


