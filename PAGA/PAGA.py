#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns 

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80)
sc.logging.print_versions()
results_file = '/home/PAGA/PAGA.h5ad'

adata = sc.read_loom("/home/yoursample.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')

adata

adata.var_names_make_unique()
adata.X = adata.X.astype('float64')
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=20)
sc.pl.paga(adata, color=['louvain'], edge_width_scale=0.2, threshold=0.2)
sc.tl.umap(adata)
adata.obsm['X_umap'][:, 1] = np.clip(adata.obsm['X_umap'][:, 1], a_min=-20, a_max=None)
sc.pl.umap(adata, color=['louvain'], legend_loc='on data')
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=31, use_rep='X_diffmap')
sc.tl.louvain(adata, resolution=0.2)
sc.pl.umap(adata, color='louvain', legend_loc='on data')
sc.tl.paga(adata, groups='clusters')
sc.pl.paga(adata, threshold=0.02, edge_width_scale=1, layout='fr', random_state=0)

pos = adata.uns['paga']['pos']
sc.pl.paga(adata, threshold=0.02, edge_width_scale=1, layout='fr', pos=pos)
sc.tl.umap(adata, init_pos='paga')
sc.pl.umap(adata, color='clusters', legend_loc='on data')

adata
#sc.pp.pca(adata, svd_solver='arpack') 
#sc.pp.neighbors(adata, n_neighbors=20, use_rep='harmony_cell_embeddings',n_pcs=20)
#sc.tl.diffmap(adata)
#sc.pp.neighbors(adata, n_neighbors=20, use_rep='X_diffmap',n_pcs=20)
sc.pp.neighbors(adata, n_neighbors=20, use_rep='umap_cell_embeddings',n_pcs=20)
#adata.uns['iroot'] = np.flatnonzero(adata.obs['seurat_clusters']  == '2')[0]

#change cluster names if need
#adata.obs['seurat_clusters_anno'] = adata.obs['sub_types']
#adata.obs['seurat_clusters_anno'].cat.categories = []
sc.tl.paga(adata, groups='new_cluster_id')
sc.pl.paga(adata, color=['new_cluster_id'],threshold=0.2,edge_width_scale=0.5, layout='fr', random_state=0)

pos = adata.uns['paga']['pos']
pos

#positions abstracted from seurat-umap 
pos[0]=[-2.60362,4.271616]
pos[1]=[-1.456266,-4.241566]
pos[2]=[-2.730585,0.9757567]
pos[3]=[2.014696,-1.116894]
pos[4]=[2.087086,2.428047]
pos[5]=[0.02982375,2.882216]
pos[6]=[3.366552,5.750513]
pos[7]=[-0.6422363,6.289839]
pos[8]=[4.570878,1.133407]
pos[9]=[5.348038,5.977872]
pos[10]=[-1.101751,-1.788876]
pos[11]=[3.809338,-2.370346]
pos[12]=[0.1240367,-8.425795]
pos[13]=[1.617643,6.250533]

fig = plt.figure(frameon = False)
for i in np.arange(0.0, 0.51, 0.01):
    plt.clf()
    sc.pl.paga(adata,color=['new_cluster_id'], threshold=i,edge_width_scale=0.2,layout='fr',pos=pos,fontsize=0,node_size_scale=0,cmap='black',frameon = None)
    fig.set_size_inches(100, 75)
    plt.savefig("/home/PAGA/PAGA_MKPMK_new_{}.pdf".format(i),dpi=800,bbox_inches = 'tight')
    plt.clf()
    plt.close("all")
