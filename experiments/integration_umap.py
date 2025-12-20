
import scanpy as sc

import os

################################################################
# DATA
# data source: # https://www.biorxiv.org/content/10.1101/2025.09.11.675724v1.full - PacBio long reads
# files paths 
CARMELLE_DATA_DIR = "/mnt/lareaulab/carmelle/longread_sc/lung"
AATH_DATA_DIR = "/Users/aathreyakadambi/Documents/school/berkeley/fa25/cs194/final_project/data"
AATH_LAMBDA_DATA_DIR = "/home/ubuntu/workspace/data"
ALEX_DATA_DIR = "/Users/aolhava/Desktop/Berkeley/2025 Fall/CS 294/Project/data"
DATA_DIR = ALEX_DATA_DIR

longbench_gene_integrated = f"{DATA_DIR}/mrVI_integrated_sc_sn_gene-exp.h5ad"
longbench_isoform_integrated = f"{DATA_DIR}/isoform_data_integrated_12.18.2025.h5ad"

################################################################

################################################################
# RESULTS
AATH_RESULTS_DIR = "/Users/aathreyakadambi/Documents/school/berkeley/fa25/cs194/final_project/results"
AATH_LAMBDA_RESULTS_DIR = "/home/ubuntu/workspace/results"
ALEX_RESULTS_DIR = "/Users/aolhava/Desktop/Berkeley/2025 Fall/CS 294/Project/cs194-final-project/results/cluster_analysis_longbench/alex"
RESULTS_DIR = ALEX_RESULTS_DIR

################################################################

sc.settings.figdir = RESULTS_DIR

genes = False
isoforms = True

# genes
if genes: 
    print("Gene UMAPS...")
    genes_combined = sc.read_h5ad(longbench_gene_integrated)
    found_genes_latent = False

    if "X_scVI" in genes_combined.obsm:
        print("X_scVI UMAP...")
        sc.pp.neighbors(genes_combined, use_rep="X_scVI")
        sc.tl.umap(genes_combined, min_dist=0.3)
        sc.pl.umap(genes_combined, color=['batch'], wspace=0.55, save = '_scvi_gene_expression.png')
        found_genes_latent = True

    if "X_scANVI" in genes_combined.obsm:
        print("X_scANVI UMAP...")
        sc.pp.neighbors(genes_combined, use_rep="X_scANVI")
        sc.tl.umap(genes_combined, min_dist=0.3)
        sc.pl.umap(genes_combined, color=['batch'], wspace=0.55, save = '_scanvi_gene_expression.png')
        found_genes_latent = True

    if "X_mrVI" in genes_combined.obsm:
        print("X_mrVI UMAP...")
        sc.pp.neighbors(genes_combined, use_rep="X_mrVI")
        sc.tl.umap(genes_combined, min_dist=0.3)
        sc.pl.umap(genes_combined, color=['batch'], wspace=0.55, save = '_mrvi_gene_expression.png')
        found_genes_latent = True

    if not found_genes_latent:
        print("Warning: X_scVI/X_scANVI/X_mrVI not found in genes_combined.obsm.")

# isoforms
if isoforms:
    print("Isoform UMAPS...")
    isoforms_combined = sc.read_h5ad(longbench_isoform_integrated)
    found_isoforms_latent = False

    if "X_scVI" in isoforms_combined.obsm:
        print("X_scVI UMAP...")
        sc.pp.neighbors(isoforms_combined, use_rep="X_scVI")
        sc.tl.umap(isoforms_combined, min_dist=0.3)
        sc.pl.umap(isoforms_combined, color=['batch'], wspace=0.55, save = '_scvi_isoform_expression.png')
        found_isoforms_latent = True

    if "X_scANVI" in isoforms_combined.obsm:
        print("X_scANVI UMAP...")
        sc.pp.neighbors(isoforms_combined, use_rep="X_scANVI")
        sc.tl.umap(isoforms_combined, min_dist=0.3)
        sc.pl.umap(isoforms_combined, color=['batch'], wspace=0.55, save = '_scanvi_isoform_expression.png')
        found_isoforms_latent = True

    if "X_mrANVI" in isoforms_combined.obsm:
        print("X_mrVI UMAP...")
        sc.pp.neighbors(isoforms_combined, use_rep="X_mrVI")
        sc.tl.umap(isoforms_combined, min_dist=0.3)
        sc.pl.umap(isoforms_combined, color=['batch'], wspace=0.55, save = '_mrvi_isoform_expression.png')
        found_isoforms_latent = True

    if not found_isoforms_latent:
        print("Warning: X_scVI/X_scANVI/X_mrVI not found in isoforms_combined.obsm.")