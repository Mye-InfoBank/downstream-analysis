import scanpy as sc
import pandas as pd
import anndata as ad
import os
import infercnvpy as cnv
import matplotlib.pyplot as plt
import seaborn as sns
# to get python and libraries installed version
# import sys
# import platform
# import pkg_resources


## READ SINGLE-CELL DATA   

# %%
adata_path = "PATH/toy_atlas.h5ad" # replace by your PATH and read .h5ad file

# %%
# Caching the mtx as anndata object saves time
if os.path.exists(adata_path):
    adata = ad.read_h5ad(adata_path)
else:
    adata = sc.read_mtx("count_matrix_sparse.mtx").transpose()
    
    barcodes = pd.read_csv("count_matrix_barcodes.tsv", header=None, sep='\t', names=['barcodes'])
    features = pd.read_csv("count_matrix_genes.tsv", header=None, sep='\t', names=['gene_names'])

    adata.obs_names = barcodes['barcodes']
    adata.var_names = features['gene_names']
    adata.obs = pd.read_csv('metadata.csv', header=0, index_col=0).loc[adata.obs_names]

    adata.write_h5ad(adata_path)

adata

## DATA PREPROCESSING: Skip this part if data pre-processing is already done !!

# %%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.obs['percent_mt'] = adata[:, adata.var['mt']].X.sum(1) / adata.X.sum(1)
adata = adata[adata.obs['percent_mt'] < 0.05, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata

## DATA FILTERING: Keep only epithelial cells
adata_Epc = adata[adata.obs['label'] == 'Epithelial_cell', :] 
# in this step you filter data keeping only cells annotated as 'Epitelial_cell' in 'label' column


## READ CHROMOSOME INFORMATION

# %%
df_loc = pd.read_csv('PATH/Chrom_info_toyAtlas_v2.csv', index_col=0) # replace PATH by your directory and read chromosome information (structure: Chromosome (Col1) | Start_Position (Col2)| End_Position (Col3)
df_loc['chromosome'] = df_loc['chromosome'].apply(lambda x: 'chr' + str(x))
var_merged = pd.merge(adata.var, df_loc, how='inner', left_index=True, right_index=True)
adata_Epc = adata_Epc[:, var_merged.index]
adata_Epc.var = var_merged


## RUN infercnvpy

# %%
cnv.tl.infercnv(adata_Epc, # AnnData of Epithelial cells
                # reference_key='label',  # Column representing the column of the proper cell type annotation contained in adata.obs (optional)
                # reference_cat='Normal Epithelial',  # Cells used as reference (optional)
                window_size=100 # Window size for analysis. Adjust it
                )
adata_Epc


# %%
cnv.pl.chromosome_heatmap(adata_Epc, groupby="label", dendrogram = True) # groupby is the appropiate cell type annotation column

# Clustering by CNV profiles and identifying tumor cells
cnv.tl.pca(adata_Epc)
cnv.pp.neighbors(adata_Epc)
cnv.tl.leiden(adata_Epc)
cnv.pl.chromosome_heatmap(adata_Epc, groupby="cnv_leiden", dendrogram=True)

# UMAP plot of CNV profiles
cnv.tl.umap(adata_Epc)
cnv.tl.cnv_score(adata_Epc)

# %%
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap(
    adata_Epc,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False,
)
cnv.pl.umap(adata_Epc, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata_Epc, color="label", ax=ax3)
plt.gcf().savefig("OUTPUT_PATH/UMAP_CLUSTER_CNV.pdf", format="pdf", dpi=300) # change OUTPUT_PATH by your output directory

## CLASSIFYING TUMOR CELLS
adata_Epc.obs["cnv_status"] = "normal" # labelled as normal all a priori normal cells
adata_Epc.obs.loc[adata_Epc.obs["cnv_leiden"].isin(["2", "9", "6"]), "cnv_status"] = (
    "tumor"
)
# change ([‘2’, ‘9’, ‘6’]) by clusters with high values of cnv scores and labelled as 'tumor'
# Yo should change tumor clusters IDs based on your results

fig, (ax1) = plt.subplots(1, figsize=(10, 5), gridspec_kw={"wspace": 0.5})
cnv.pl.umap(adata_Epc, color="cnv_status", ax=ax1, show=True)
plt.savefig("OUTPUT_PATH/UMAP_tumor_vs_normal", format="pdf") # change OUTPUT_PATH by your output directory
plt.close(fig) 

cnv.pl.chromosome_heatmap(adata_Epc[adata_Epc.obs["cnv_status"] == "tumor", :])
cnv.pl.chromosome_heatmap(adata_Epc[adata_Epc.obs["cnv_status"] == "normal", :])

# Write new metadata as .csv
metadata_cnv = adata_Epc.obs
metadata_cnv.to_csv("/OUTPUT_PATH/metadata_finalcnvToy.csv", sep="\t", index=False, )
metadata_cnv.reset_index().rename(columns={"index": "RowNames"}).to_csv("OUTPUT_PATH/metadata_finalcnvToy.csv", index=False) # change OUTPUT_PATH by your output directory

# Show CNV values
# Create dataframe containing variables of our interest
cnv_scores = metadata_cnv[["cnv_score", "cnv_leiden"]]  # cnv_score represents infercnvpy results and cnv_leiden is the cathegory variable of cnv cluster. Both of them from pipeline results
cnv_scores = cnv_scores.rename(columns={"cnv_leiden": "Group", "cnv_score": "Value"})

# Show CNV by Dot Plot
plt.figure(figsize=(10, 6))
sns.stripplot(x="Group", y="Value", data=cnv_scores, color="#9370DB", jitter=True, size=5)
plt.title("Distribution of CNV Scores per CNV clusters")
plt.xlabel("CNV cluster")
plt.ylabel("CNV Score")
plt.xticks(rotation=45)
plt.show()

# Show CNV by Swarm Plot: Dot Plot alternative but with non-overlapping dots
plt.figure(figsize=(10, 6))
sns.swarmplot(x="Group", y="Value", data=cnv_scores, color="#9370DB", size=5)
plt.title("Distribution of CNV Scores per CNV clusters")
plt.xlabel("CNV cluster")
plt.ylabel("CNV Score")
plt.xticks(rotation=45)
plt.show()
# %%

# *** Python version info ***
# print(f"Python version: {sys.version}")
# Python version: 3.10.12 (main, Mar 22 2024, 16:50:05) [GCC 11.4.0]

# print(f"Python executable: {sys.executable}")
# Python executable: /bin/python3

# print(f"Platform: {platform.system()} {platform.release()}")
# Platform: Linux 5.19.0-38-generic

# Info about python libraries installed
# installed_packages = pkg_resources.working_set
# for dist in installed_packages:
#    print(f"{dist.project_name}=={dist.version}")

# Libraries of 
# scanpy==1.10.3
# pandas==2.2.3
# anndata==0.10.9
# infercnvpy==0.5.0
# matplotlib==3.9.2
# matplotlib-inline==0.1.7
# seaborn==0.13.2