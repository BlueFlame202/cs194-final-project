import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
#DE - anirudh
#load
adata_gene = sc.read_h5ad("post_scvi_scanvi_integrated_gene.h5ad")
adata = sc.read_h5ad("post_scvi_scanvi_integrated_isoform.h5ad")
#leiden clustering from aathreya
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, key_added="X_scVI_umap")
sc.tl.leiden(adata, resolution=0.6, key_added="scvi_leiden_class")

h1975_nuclei = (adata.obs['scvi_leiden_class'].isin(['3', '5'])) & (adata.obs['dataset'].str.contains('nuclei', case=False, na=False))

adata_sub = adata[h1975_nuclei].copy()

sc.tl.rank_genes_groups(
    adata_sub,
    groupby='scvi_leiden_class',
    groups=['3' and '5'],     
    method='wilcoxon',
    use_raw=False,    
    layer='counts',   
    n_genes=adata_sub.n_vars
)


de_iso = sc.get.rank_genes_groups_df(adata_sub, group='3' and '5')
sig = (de_iso["pvals_adj"] < 0.05) & (np.abs(de_iso["logfoldchanges"]) > 0.25)

# subset to sig hits
de_sig = de_iso.loc[sig].copy()

de_sig.to_csv("isoform_DE_cluster3_vs5_significant.csv", index=False)
print(de_sig.head())

# Top 20 isoforms upregulated in cluster 3
#print(de_iso.sort_values('logfoldchanges', ascending=False).head(20))

plt.figure(figsize=(6,5))
plt.scatter(de_iso['logfoldchanges'], -np.log10(de_iso['pvals_adj']), s=10, alpha=0.6)
plt.axvline(0.25, linestyle='--', color='red', label='log2FC > 0.25')
plt.axvline(-0.25, linestyle='--', color='red')
plt.axhline(-np.log10(0.05), linestyle='--', color='blue', label='adjusted p value < 0.05')
plt.xlabel('log2FC (cluster 3 / 5)')
plt.ylabel('-log10(adjusted p value)')
plt.title('Isoform DE: Cluster 3 vs 5 (nuclei)')
plt.legend()
plt.tight_layout()
plt.show()