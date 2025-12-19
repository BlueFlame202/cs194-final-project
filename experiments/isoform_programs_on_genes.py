import scanpy as sc
import scvi
import pandas as pd

sc.settings.figdir = '/Users/aolhava/Desktop/Berkeley/2025 Fall/CS 294/Project/figures'

def main():
    # loading data and leiden clustering
    genes_combined = sc.read_h5ad('data/mrVI_integrated_sc_sn_gene-exp.h5ad')
    isoforms_combined = sc.read_h5ad('data/scanvi_integrated_sc_sn_isoform-exp.h5ad')
    sc.tl.leiden(isoforms_combined, resolution=0.5, flavor = "igraph")

    sc.pp.neighbors(isoforms_combined, use_rep='X_scVI')
    sc.tl.umap(isoforms_combined, min_dist=0.3)
    sc.pl.umap(isoforms_combined, color=['batch', 'leiden'], save = '_scvi_isoform_expression_3.png')

    # combining gene expression and isoform cluster then plotting it =
    genes_combined.obs['isoform_cluster'] = isoforms_combined.obs['leiden']

    sc.pp.neighbors(genes_combined, use_rep='X_mrVI')
    sc.tl.umap(genes_combined, min_dist=0.3)
    sc.pl.umap(genes_combined, color=['batch', 'isoform_cluster'], save = '_gene_expression_scanVI_isoform_cluster.png')

if __name__ == "__main__":
    main()