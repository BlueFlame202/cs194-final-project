
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
sn_isoform = f'{DATA_DIR}/longbench_data/sn/GSM9135509_pb_sn_transcript_count_matrix.mtx.gz'
# sc
sc_isoform = f'{DATA_DIR}/longbench_data/sc/GSM9135508_pb_sc_transcript_count_matrix.mtx.gz'

sn_features_path = f'{DATA_DIR}/longbench_data/sn/GSM9135509_pb_sn_transcript_count_features.tsv.gz'
sn_cells_path = f'{DATA_DIR}/longbench_data/sn/GSM9135509_pb_sn_transcript_count_barcodes.tsv.gz'

sc_features_path = f'{DATA_DIR}/longbench_data/sc/GSM9135508_pb_sc_transcript_count_features.tsv.gz'
sc_cells_path = f'{DATA_DIR}/longbench_data/sc/GSM9135508_pb_sc_transcript_count_barcodes.tsv.gz'

# cell line data
sc_cell_lines_path = f'{DATA_DIR}/longbench_data/cell_line_labels/sc_cell_line_labels.csv'
sn_cell_lines_path = f'{DATA_DIR}/longbench_data/cell_line_labels/sn_cell_line_labels.csv'
################################################################

################################################################
# RESULTS
AATH_RESULTS_DIR = "/Users/aathreyakadambi/Documents/school/berkeley/fa25/cs194/final_project/results"
AATH_LAMBDA_RESULTS_DIR = "/home/ubuntu/workspace/results"
RESULTS_DIR = AATH_LAMBDA_RESULTS_DIR

EXP_NAME = "sc_sn_isoform_integration_longbench"
EXP_RESULTS_DIR = os.path.join(RESULTS_DIR, EXP_NAME)
EXP_PROC_DATA_DIR = DATA_DIR + "/processed/" + EXP_NAME
os.makedirs(EXP_RESULTS_DIR, exist_ok=True)
os.makedirs(EXP_PROC_DATA_DIR, exist_ok=True)
################################################################

# ------------------------------
# Load data
# ------------------------------
print("Loading single nuclei and single cell isoform expression matrices...")
adata_sn_isoform = sc.read_mtx(sn_isoform)
adata_sc_isoform = sc.read_mtx(sc_isoform)

print(adata_sc_isoform.X)
print(adata_sn_isoform.X)

sn_features = pd.read_csv(sn_features_path, header=None)[0].tolist()
sn_cells = pd.read_csv(sn_cells_path, header=None)[0].tolist()

sc_features = pd.read_csv(sc_features_path, header=None)[0].tolist()
sc_cells = pd.read_csv(sc_cells_path, header=None)[0].tolist()

adata_sn_isoform.var_names = sn_features
adata_sn_isoform.obs_names = sn_cells

adata_sc_isoform.var_names = sc_features
adata_sc_isoform.obs_names = sc_cells

# QC
print("Calculating QC metrics...")
sc.pp.calculate_qc_metrics(adata_sc_isoform, inplace=True)
sc.pp.calculate_qc_metrics(adata_sn_isoform, inplace=True)

adata_sn_isoform.obs["batch"] = "single nuclei"
adata_sc_isoform.obs["batch"] = "single cell"

adata_sn_isoform.obs_names_make_unique()
adata_sc_isoform.obs_names_make_unique()
adata_sn_isoform.var_names_make_unique()
adata_sc_isoform.var_names_make_unique()

print("Concatenating anndata objects...")
isoforms_combined = sc.concat(
    [adata_sc_isoform, adata_sn_isoform],
    join="inner",
    label="dataset",
    keys=["single cell", "single nuclei"]
)

isoforms_combined.var_names_make_unique()
isoforms_combined.obs_names_make_unique()

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
isoforms_combined.obs["cell_line"] = isoforms_combined.obs.index.map(cell_line_map)
isoforms_combined.obs["cell_line"] = isoforms_combined.obs["cell_line"].fillna("Unknown")
isoforms_combined.obs["cell_line"] = isoforms_combined.obs["cell_line"].astype("category")

# ------------------------------
# Filtering
# ------------------------------
print("Filtering cells and isoforms...")
sc.pp.filter_cells(isoforms_combined, min_genes=100)
sc.pp.filter_genes(isoforms_combined, min_cells=3)

# HVGs
print("Selecting highly variable isoforms (HVIs)...")
isoforms_combined.layers["counts"] = isoforms_combined.X
sc.pp.highly_variable_genes(
    isoforms_combined, 
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
    isoforms_combined,
    layer="counts",
    batch_key="batch"
)

# Train scVI
print("Training SCVI model...")
vae = scvi.model.SCVI(isoforms_combined, n_layers=2, n_latent=20)
vae.train(max_epochs=250, accelerator="gpu")
history_vae = vae.history

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
isoforms_combined.obsm["X_scVI"] = latent_scvi  # save scVI embedding

# Train scANVI with error handling
try:
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
    isoforms_combined.obsm["X_scANVI"] = latent_scanvi

except Exception as e:
    print("Warning: SCANVI training failed.")
    print(e)
    print("Proceeding with scVI latent space only.")

# ------------------------------
# Save integrated AnnData
# ------------------------------
output_path = f"{EXP_PROC_DATA_DIR}/post_scvi_scanvi_integrated.h5ad"
print(f"Writing integrated AnnData to {output_path}")
isoforms_combined.write_h5ad(output_path)