#!/usr/bin/env python
# coding: utf-8

# import dependencies
import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import loompy as lp
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss

from IPython.display import HTML, display
sc.settings.njobs = 32
from MulticoreTSNE import MulticoreTSNE as TSNE

RESOURCES_FOLDERNAME = "/home/SCENIC/resources/"
AUXILLIARIES_FOLDERNAME = "/home/SCENIC/auxilliaries/"
RESULTS_FOLDERNAME = "/home/SCENIC/results/"
FIGURES_FOLDERNAME = "/home/SCENIC/figures/"
sc.settings.figdir = FIGURES_FOLDERNAME
BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"
def savesvg(fname: str, fig, folder: str=FIGURES_FOLDERNAME) -> None:
    """
    Save figure as vector-based SVG image format.
    """
    fig.tight_layout()
    fig.savefig(os.path.join(folder, fname), format='svg')
def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', -1)
    display(HTML(df.head().to_html(escape=False)))
    pd.set_option('display.max_colwidth', MAX_COL_WIDTH)
# Downloaded fromm pySCENIC github repo: https://github.com/aertslab/pySCENIC/tree/master/resources
HUMAN_TFS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'lambert2018.txt')
# Ranking databases. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
RANKING_DBS_FNAMES = list(map(lambda fn: os.path.join(AUXILLIARIES_FOLDERNAME, fn),
                       ['hg19-500bp-upstream-10species.mc9nr.feather',
                       'hg19-tss-centered-5kb-10species.mc9nr.feather',
                        'hg19-tss-centered-10kb-10species.mc9nr.feather']))
# Motif annotations. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
MOTIF_ANNOTATIONS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'motifs-v9-nr.hgnc-m0.001-o0.0.tbl')
DATASET_ID = "yoursample"
METADATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.metadata.csv'.format(DATASET_ID))
EXP_MTX_QC_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.qc.tpm.csv'.format(DATASET_ID))
ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.adjacencies.tsv'.format(DATASET_ID))
MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.motifs.csv'.format(DATASET_ID))
REGULONS_DAT_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.regulons.dat'.format(DATASET_ID))
AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.auc.csv'.format(DATASET_ID))
BIN_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.bin.csv'.format(DATASET_ID))
THR_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.thresholds.csv'.format(DATASET_ID))
ANNDATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.h5ad'.format(DATASET_ID))
LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.loom'.format( DATASET_ID))

adata = sc.read_loom("/home/SCENIC/yoursample.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')

#change the cluster names 
#adata.obs.seurat_clusters= adata.obs.seurat_clusters.replace({'1': ''...})
adata.to_df().to_csv(EXP_MTX_QC_FNAME)

os.system("pyscenic grn /home/SCENIC/results/yoursample.qc.tpm.csv /home/SCENIC/auxilliaries/lambert2018-edit.txt -o /home/SCENIC/results/yoursample.adjacencies.tsv --num_workers 15 --seed 777")

DBS_PARAM = ' '.join(RANKING_DBS_FNAMES)
os.system("pyscenic ctx /home/SCENIC/results/yoursample.adjacencies.tsv \
                 /home/SCENIC/auxilliaries/hg19-500bp-upstream-10species.mc9nr.feather  \
                 /home/SCENIC/auxilliaries/hg19-tss-centered-5kb-10species.mc9nr.feather \
                 /home/SCENIC/auxilliaries/hg19-tss-centered-10kb-10species.mc9nr.feather \
                 --annotations_fname /home/SCENIC/auxilliaries/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
                 --expression_mtx_fname /home/SCENIC/results/yoursample.qc.tpm.csv \
                 --output /home/SCENIC/results/yoursample.motifs.csv \
                 --nes_threshold 2 \
                 --num_workers 26")

regulons = load_signatures("/home/SCENIC/results/yoursample.motifs.csv")
auc_mtx = aucell(adata.to_df(), regulons, num_workers=26)
auc_mtx.to_csv("/home/SCENIC/results/yoursample.auc.csv")
#auc_mtx = pd.read_csv("/home/SCENIC/results/yoursample.auc.csv", index_col=0)
bin_mtx, thresholds = binarize(auc_mtx,num_workers=26)
bin_mtx.to_csv("/home/SCENIC/results/yoursample.bin.csv")
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv("/home/SCENIC/results/yoursample.thresholds.csv")
bin_mtx = pd.read_csv(BIN_MTX_FNAME, index_col=0)
thresholds = pd.read_csv(THR_FNAME, index_col=0).threshold

def savesvg(fname: str, fig, folder: str='/data/tmp_data/zxy/SCENIC/results/') -> None:
    """
    Save figure as vector-based SVG image format.
    """
    fig.tight_layout()
    fig.savefig(os.path.join(folder, fname), format='svg')
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(8, 4), dpi=100)
plot_binarization(auc_mtx, 'FOXA3(+)', thresholds['FOXA3(+)'], ax=ax1)
plot_binarization(auc_mtx, 'IKZF3(+)', thresholds['IKZF3(+)'], ax=ax2)
plot_binarization(auc_mtx, 'MSC(+)', thresholds['MSC(+)'], ax=ax3)
plot_binarization(auc_mtx, 'SPIC(+)', thresholds['SPIC(+)'], ax=ax4)
plot_binarization(auc_mtx, 'ZNF675(+)', thresholds['ZNF675(+)'], ax=ax5)
plot_binarization(auc_mtx, 'ZNF879(+)', thresholds['ZNF879(+)'], ax=ax6)
plt.tight_layout()
savesvg('yoursample.test.binarization.svg', fig)

#RSS
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import matplotlib as mpl
import matplotlib.pyplot as plt
from scanpy.plotting._tools.scatterplots import plot_scatter
import seaborn as sns
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import seaborn as sns
from pyscenic.binarization import binarize
import matplotlib.pyplot as plt

rss_cellType = regulon_specificity_scores( auc_mtx, adata.obs['cell_types'] )
rss_cellType
from adjustText import adjust_text

cats = sorted(list(set(adata.obs['cell_types'])))
topreg = []
for i,c in enumerate(cats):
    topreg.extend(
        list(rss_cellType.T[c].sort_values(ascending=False)[:10].index)
    )
topreg = list(set(topreg))
