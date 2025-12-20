# Exploring splicing landscapes at the single-cell level using long reads

## Experiments

### Longbench LUAD Dataset

Simple exploration of the LUAD dataset is done in `notebooks/data_raw_viz.ipynb` along with all other notebooks aside from `notebooks/CRC_isoform_clusters_analysis.ipynb`.

Code for reproducing the single nuclei and single cell batches is in the `experiments/*integration_gene.py` and `experiments/*integration_isoform.py` files. Cluster analysis is in `experiments/cluster_analysis_longbench*.py files`, differential expression analysis is in `experiments/DE_anirudh.py`, and isoform programs analysis is in `experiments/isoform_programs_on_genes.py`. Isoform clusters analysis is shown in `notebooks/lung-cell-line-isoform-clusters-analysis.ipynb`.

### CRC Dataset

Analysis and visualization of the CRC dataset is done in `notebooks/crc_isoform_clusters_analysis.ipynb`.