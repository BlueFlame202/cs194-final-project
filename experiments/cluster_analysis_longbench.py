
import scanpy as sc
import scvi

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import os

from sklearn.metrics import silhouette_score

################################################################
# DATA
# data source: # https://www.biorxiv.org/content/10.1101/2025.09.11.675724v1.full - PacBio long reads
# files paths 
CARMELLE_DATA_DIR = "/mnt/lareaulab/carmelle/longread_sc/lung"
AATH_DATA_DIR = "/Users/aathreyakadambi/Documents/school/berkeley/fa25/cs194/final_project/data"
AATH_LAMBDA_DATA_DIR = "/home/ubuntu/workspace/data"
DATA_DIR = AATH_LAMBDA_DATA_DIR

longbench_gene_integrated = f"{DATA_DIR}/processed/sc_sn_gene_integration_longbench/post_scvi_scanvi_integrated.h5ad"
longbench_isoform_integrated = f"{DATA_DIR}/processed/sc_sn_isoform_integration_longbench/post_scvi_scanvi_integrated.h5ad"
################################################################

################################################################
# RESULTS
AATH_RESULTS_DIR = "/Users/aathreyakadambi/Documents/school/berkeley/fa25/cs194/final_project/results"
AATH_LAMBDA_RESULTS_DIR = "/home/ubuntu/workspace/results"
RESULTS_DIR = AATH_LAMBDA_RESULTS_DIR

EXP_NAME = "cluster_analysis_longbench"
EXP_RESULTS_DIR = os.path.join(RESULTS_DIR, EXP_NAME)
EXP_PROC_DATA_DIR = DATA_DIR + "/processed/" + EXP_NAME
os.makedirs(EXP_RESULTS_DIR, exist_ok=True)
os.makedirs(EXP_PROC_DATA_DIR, exist_ok=True)
################################################################

genes_combined = sc.read_h5ad(longbench_gene_integrated)
isoforms_combined = sc.read_h5ad(longbench_isoform_integrated)

# begin ChatGPT
def leiden_silhouette(adata, latent_key, resolutions, prefix):
    """
    Run neighbors+Leiden at each resolution on adata using latent_key,
    compute silhouette scores, and save results to adata.obs.
    """
    # If neighbors have already been computed for another rep, clear them
    keys_to_remove = ['neighbors']
    for k in keys_to_remove:
        if k in adata.uns:
            # Removing old neighbors info; optional but avoids confusion
            del adata.uns[k]

    scores = []
    for res in resolutions:
        # build neighbors graph based on the latent space
        sc.pp.neighbors(adata, use_rep=latent_key)
        
        # run Leiden clustering
        leiden_key = f"{prefix}_leiden_res_{res:.2f}"
        sc.tl.leiden(adata, resolution=res, key_added=leiden_key)
        
        # compute silhouette score for clusters
        labels = adata.obs[leiden_key].astype(str)
        sil = silhouette_score(adata.obsm[latent_key], labels)
        scores.append(sil)
        print(f"{latent_key} @ res={res:.2f}: silhouette={sil:.3f}")

    # plot silhouette scores
    plt.figure(figsize=(5,4))
    plt.plot(resolutions, scores, marker="o")
    plt.xlabel("resolution")
    plt.ylabel("silhouette score")
    plt.title(f"Leiden silhouette ({latent_key})")
    plt.tight_layout()
    plt.savefig(f"{EXP_RESULTS_DIR}/{prefix}_leiden_silhouette.png")
    plt.close()

    return scores

# Gene latent
gene_resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]
gene_sils = leiden_silhouette(genes_combined, "X_scVI", gene_resolutions, "gene_scVI")

# scANVI latent
scanvi_resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]
scanvi_sils = leiden_silhouette(genes_combined, "X_scANVI", scanvi_resolutions, "scanvi")
# end ChatGPT