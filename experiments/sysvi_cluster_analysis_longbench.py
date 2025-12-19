
import scanpy as sc

import matplotlib.pyplot as plt

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

longbench_gene_integrated = f"{DATA_DIR}/processed/sysvi_sc_sn_gene_integration_longbench/post_sysvi_integrated.h5ad"
longbench_isoform_integrated = f"{DATA_DIR}/processed/sysvi_sc_sn_isoform_integration_longbench/post_sysvi_integrated.h5ad"
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
def leiden_silhouette(adata, latent_key, resolutions, prefix, plot_clusters_umap=None):
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

    # begin Cursor AI
    if plot_clusters_umap:
        # Define UMAP key for this latent space
        umap_key = f"{latent_key}_umap"
        # construct 2x3 grid
        fig, axs = plt.subplots(2, 3, figsize=(16, 10))
        axs = axs.flatten()
        n = len(resolutions)
        for i, res in enumerate(resolutions):
            leiden_key = f"{prefix}_leiden_res_{res:.2f}"
            # compute UMAP if not already present for this latent space
            if umap_key not in adata.obsm:
                sc.tl.umap(
                    adata, 
                    min_dist=0.3, 
                    spread=1.0, 
                    random_state=123, 
                    key_added=umap_key, 
                    neighbors_key=None
                )
            umap_coords = adata.obsm[umap_key]
            labels = adata.obs[leiden_key].astype(str)

            # Pick a colormap with enough distinguishable colors, or fallback
            unique_labels = sorted(labels.unique())
            n_clusters = len(unique_labels)
            cmap = plt.get_cmap('tab20' if n_clusters <= 20 else 'tab20b')

            for idx, group in enumerate(unique_labels):
                group_mask = labels == group
                axs[i].scatter(
                    umap_coords[group_mask, 0],
                    umap_coords[group_mask, 1],
                    label=f'Cluster {group}',
                    alpha=0.7,
                    s=8,
                    color=cmap(idx % cmap.N)
                )
            axs[i].set_title(f"Leiden (res={res:.2f})")
            axs[i].set_xlabel("UMAP1")
            axs[i].set_ylabel("UMAP2")
            axs[i].legend(markerscale=2, fontsize=8, loc='best', frameon=False)

        # Remove unused subplots if any
        for j in range(n, 6):
            fig.delaxes(axs[j])

        fig.suptitle(f"{prefix}: UMAP colored by Leiden clusters", fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])
        plt.savefig(f"{EXP_RESULTS_DIR}/{prefix}_leiden_umap_grid.png")
        plt.close()

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

# Generate UMAP embeddings for genes_combined and isoforms_combined

# For gene-level data
if "X_scVI" in genes_combined.obsm:
    print("Computing UMAP for genes_combined (X_sysVI)...")
    sc.pp.neighbors(genes_combined, use_rep="X_sysVI")
    sc.tl.umap(genes_combined, key_added="X_sysVI_umap")

# For isoform-level data
if "X_scVI" in isoforms_combined.obsm:
    print("Computing UMAP for isoforms_combined (X_scVI)...")
    sc.pp.neighbors(isoforms_combined, use_rep="X_sysVI")
    sc.tl.umap(isoforms_combined, key_added="X_sysVI_umap")


# Gene latent
resolutions = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
gene_sils = leiden_silhouette(genes_combined, "X_sysVI", resolutions, "gene_sysVI", plot_clusters_umap=True)

# Isoform latent
isoform_sils = leiden_silhouette(isoforms_combined, "X_sysVI", resolutions, "isoform_sysVI", plot_clusters_umap=True)
