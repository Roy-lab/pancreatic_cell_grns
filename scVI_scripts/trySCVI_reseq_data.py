#This is just an example redoing what is in t github link https://docs.scvi-tools.org/en/stable/tutorials/notebooks/harmonization.html
#We do the following to play with the data

##Before running this install scvitools, scanpy, leiden and  pymde

import scanpy as sc
import pandas as pd
import os

import scvi
input_file='/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/scVI_application/input_data/resequenced_monocle3_normalized_matrix_nodup_with_genotypes_reordered_noheader_t.txt'
anndata = sc.read_csv(input_file,delimiter='\t',dtype='float32')

samplenames='/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/scVI_application/input_data/samplenames.txt'
batchnames='/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/scVI_application/input_data/batchnames.txt'
genenames='/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/scVI_application/input_data/genenames.txt'
sdata = pd.read_table(samplenames, index_col=0, header=None, names=['samples'], dtype=str)
bdata = pd.read_table(batchnames, index_col=0, header=None, names=['batch'], dtype=str)
gdata=pd.read_table(genenames,header=None,names=['genenames'], dtype=str);
anndata.obs.index=sdata.index
anndata.obs['samples']=sdata['samples']
anndata.obs['batch']=bdata['batch']
anndata.var_names=gdata['genenames']

anndata.layers["counts"] = anndata.X.copy()

##make a copy
anndata.raw=anndata
#sc.pp.highly_variable_genes(
#    anndata,
#    flavor="seurat_v3",
#    n_top_genes=2000,
#    layer="counts",
#    batch_key="samples",
#    subset=True
#)
##After this step, anndata becomes a 81k by 2000 matrix. With the backup staying in the anndata.raw object

#try scvi tools:
#scvi.model.SCVI.setup_anndata(anndata, layer="counts", batch_key="samples")
scvi.model.SCVI.setup_anndata(anndata, layer="counts", batch_key="batch")
vae = scvi.model.SCVI(anndata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()

anndata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(anndata, use_rep="X_scVI")
#sc.tl.leiden(anndata)


resolutions_to_try = [0.2,0.4,0.6,0.8,1.0]

for resolution in resolutions_to_try:
    sc.tl.leiden(anndata, resolution=resolution, key_added="leiden_" + str(resolution))
    num_clusters = len(anndata.obs["leiden_" + str(resolution)].unique())
    cluster_info = anndata.obs[["leiden_" + str(resolution)]]
    cluster_info.to_csv(f"cluster_info_resolution_{resolution}.txt")


import matplotlib.pyplot as plt
os.mkdir("/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/scVI_application/results2")
os.chdir("/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/scVI_application/results2")
from scvi.model.utils import mde
import pymde
anndata.obsm["X_mde"] = mde(anndata.obsm["X_scVI"])
sc.pl.embedding(anndata, basis="X_mde", color=["samples"], frameon=False, ncols=1)
plt.savefig('vae_reseq.pdf', transparent=True, format = "pdf", bbox_inches='tight')

##Save vae mode
vae.save('vae_reseq',save_anndata=True);
##Save ann data
anndata.write('anndata_scvi_reseq')

##Save leiden clusters
#anndata.obs.leiden.to_csv('scvi_clusters.txt',sep="\t")
import numpy as np
#anndata.obsm["normcounts"] = vae.get_normalized_expression()
#anndata.obsm["normcounts_scaled"] = vae.get_normalized_expression(library_size=100000)
##New, based on KE's query
anndata.obsm["normcounts_corrected_scaled"] = vae.get_normalized_expression(transform_batch=['ire1','xbp1'],library_size=100000,n_samples=7)
anndata.obsm["normcounts_corrected_scaled_10k"] = vae.get_normalized_expression(transform_batch=['ire1','xbp1'],library_size=10000,n_samples=7)
anndata.obsm["normcounts_corrected_scaled_1k"] = vae.get_normalized_expression(transform_batch=['ire1','xbp1'],library_size=1000,n_samples=7)
##Unfortunately anndata could not be saved with such a big matrix.
np.savetxt("resequenced_scVI_mde_coord.txt",anndata.obsm["X_mde"],delimiter="\t");
#np.savetxt("scRNAseqalone_and_scRNAseq_MO_D3_scvi_normalized_batch.txt",anndata.obsm["normcounts"],delimiter="\t")
#np.savetxt("scRNAseqalone_and_scRNAseq_MO_D3_scvi_normalized_scaled_batch.txt",anndata.obsm["normcounts_scaled"],delimiter="\t")
np.savetxt("resequenced_scVI_normalized_corrected_scaled_batch_t.txt",anndata.obsm["normcounts_corrected_scaled"],delimiter="\t")
np.savetxt("resequenced_scVI_normalized_corrected_scaled_10K_batch_t.txt",anndata.obsm["normcounts_corrected_scaled_10K"],delimiter="\t")
np.savetxt("resequenced_scVI_normalized_corrected_scaled_1K_batch_t.txt",anndata.obsm["normcounts_corrected_scaled_1K"],delimiter="\t")

