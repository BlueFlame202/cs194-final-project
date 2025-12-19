
import scanpy as sc
import scvi

import matplotlib.pyplot as plt
import pandas as pd

import os

################################################################
# DATA
# data source: # https://www.biorxiv.org/content/10.1101/2025.09.11.675724v1.full - PacBio long reads
# files paths 
CARMELLE_DATA_DIR = "/mnt/lareaulab/carmelle/longread_sc/lung"
AATH_DATA_DIR = "/Users/aathreyakadambi/Documents/school/berkeley/fa25/cs194/final_project/data"
AATH_LAMBDA_DATA_DIR = "/home/ubuntu/workspace/data"
DATA_DIR = AATH_LAMBDA_DATA_DIR

# sn
sn_gene = f'{DATA_DIR}/longbench_data/sn/GSM9135509_pb_sn_genes_count_matrix.tsv.gz'
# sc
sc_gene = f'{DATA_DIR}/longbench_data/sc/GSM9135508_pb_sc_genes_count_matrix.tsv.gz'

# cell line data
sc_cell_lines_path = f'{DATA_DIR}/longbench_data/cell_line_labels/sc_cell_line_labels.csv'
sn_cell_lines_path = f'{DATA_DIR}/longbench_data/cell_line_labels/sn_cell_line_labels.csv'
################################################################

################################################################
# RESULTS
AATH_RESULTS_DIR = "/Users/aathreyakadambi/Documents/school/berkeley/fa25/cs194/final_project/results"
AATH_LAMBDA_RESULTS_DIR = "/home/ubuntu/workspace/results"
RESULTS_DIR = AATH_LAMBDA_RESULTS_DIR

EXP_NAME = "sc_sn_gene_integration_longbench"
EXP_RESULTS_DIR = os.path.join(RESULTS_DIR, EXP_NAME)
EXP_PROC_DATA_DIR = DATA_DIR + "/processed/" + EXP_NAME
os.makedirs(EXP_RESULTS_DIR, exist_ok=True)
os.makedirs(EXP_PROC_DATA_DIR, exist_ok=True)
################################################################

# ------------------------------
# Load data
# ------------------------------
print("Loading single nuclei and single cell gene expression matrices...")
adata_sn_gene = sc.read_csv(sn_gene).T
adata_sc_gene = sc.read_csv(sc_gene).T

# QC
print("Calculating QC metrics...")
sc.pp.calculate_qc_metrics(adata_sc_gene, inplace=True)
sc.pp.calculate_qc_metrics(adata_sn_gene, inplace=True)

adata_sn_gene.obs["batch"] = "single nuclei"
adata_sc_gene.obs["batch"] = "single cell"

adata_sn_gene.obs_names_make_unique()
adata_sc_gene.obs_names_make_unique()
adata_sn_gene.var_names_make_unique()
adata_sc_gene.var_names_make_unique()

print("Concatenating anndata objects...")
genes_combined = sc.concat(
    [adata_sc_gene, adata_sn_gene],
    join="inner",
    label="dataset",
    keys=["single cell", "single nuclei"]
)

# ------------------------------
# Cell line mapping
# ------------------------------
print("Loading cell line labels and mapping to cells...")

sc_cell_lines = pd.read_csv(sc_cell_lines_path)
sn_cell_lines = pd.read_csv(sn_cell_lines_path)

combined_cell_lines = pd.concat([sc_cell_lines, sn_cell_lines])
combined_cell_lines.drop_duplicates(subset=["cell_id"], inplace=True)
filtered_cell_lines = combined_cell_lines[combined_cell_lines["correlation"] > 0.4]

cell_line_map = filtered_cell_lines.set_index("cell_id")["assigned_cell_line"]
genes_combined.obs["cell_line"] = genes_combined.obs.index.map(cell_line_map)
genes_combined.obs["cell_line"] = genes_combined.obs["cell_line"].fillna("Unknown")
genes_combined.obs["cell_line"] = genes_combined.obs["cell_line"].astype("category")

# ------------------------------
# Filtering
# ------------------------------
print("Filtering cells and genes...")
sc.pp.filter_cells(genes_combined, min_genes=100)
sc.pp.filter_genes(genes_combined, min_cells=3)

# HVGs
print("Selecting highly variable genes (HVGs)...")
genes_combined.layers["counts"] = genes_combined.X
sc.pp.highly_variable_genes( # ask alex and carmelle about this one
    genes_combined, 
    layer="counts",
    flavor="seurat_v3",
    n_top_genes=3000,
    subset=True
)

# ------------------------------
# SCVI/SCANVI training
# ------------------------------
print("Setting up AnnData for SCVI/SCANVI...")
scvi.model.SCVI.setup_anndata(
    genes_combined,
    layer="counts",
    batch_key="batch"
)

# Train scVI
print("Training SCVI model...")
vae = scvi.model.SCVI(genes_combined, n_layers=2, n_latent=20)
vae.train(max_epochs=250, accelerator="gpu")
history_vae = vae.history
print(vae.history.keys())

# Plot and save scVI training loss
scvi_loss_plot_path = os.path.join(EXP_RESULTS_DIR, "scvi_training_loss.png")
plt.figure(figsize=(6,4))
plt.plot(history_vae["elbo_train"], label="train ELBO")
if "elbo_validation" in history_vae:
    plt.plot(history_vae["elbo_validation"], label="validation ELBO")
plt.xlabel("Epoch")
plt.ylabel("ELBO")
plt.title("scVI training loss")
plt.legend()
plt.tight_layout()
plt.savefig(scvi_loss_plot_path)
plt.close()

# Save scVI latent representation
print("Extracting SCVI latent representation...")
latent_scvi = vae.get_latent_representation()
genes_combined.obsm["X_scVI"] = latent_scvi  # save scVI embedding

# Train scANVI using pretrained scVI
print("Training SCANVI model...")
scanvi = scvi.model.SCANVI.from_scvi_model(
    vae,
    unlabeled_category="Unknown",
    labels_key="cell_line" # use as noisy proxy for cell type
)
scanvi.train(max_epochs=250, accelerator="gpu")
history_scanvi = scanvi.history

# Plot and save scANVI training loss
scanvi_loss_plot_path = os.path.join(EXP_RESULTS_DIR, "scanvi_training_loss.png")
plt.figure(figsize=(6,4))
plt.plot(history_scanvi["elbo_train"], label="train ELBO")
if "elbo_validation" in history_scanvi:
    plt.plot(history_scanvi["elbo_validation"], label="validation ELBO")
plt.xlabel("Epoch")
plt.ylabel("ELBO")
plt.title("scANVI training loss")
plt.legend()
plt.tight_layout()
plt.savefig(scanvi_loss_plot_path)
plt.close()

# Save scANVI latent representation
print("Extracting SCANVI latent representation...")
latent_scanvi = scanvi.get_latent_representation()
genes_combined.obsm["X_scANVI"] = latent_scanvi  # save scANVI embedding

# Save integrated AnnData
output_path = f"{EXP_PROC_DATA_DIR}/post_scvi_scanvi_integrated.h5ad"
print(f"Writing integrated AnnData to {output_path}")
genes_combined.write_h5ad(output_path)