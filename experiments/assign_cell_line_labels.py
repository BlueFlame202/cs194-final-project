import pandas as pd
import scanpy as sc
import numpy as np
import os
import glob
from scipy.spatial.distance import cdist

def load_bulk_data(bulk_dir):
    bulk_profiles = {}
    
    # Get all bulk files
    files = glob.glob(os.path.join(bulk_dir, "*_matrix.tsv.gz"))
    print(f"Found {len(files)} bulk files.")
    
    for f in files:
        # Extract cell line name from filename
        # Format: GSM9135496_pb_bulk_H1975_matrix.tsv.gz -> H1975
        basename = os.path.basename(f)
        parts = basename.split('_')
        # Assuming format GSM..._pb_bulk_[CELL_LINE]_matrix.tsv.gz
        # parts: [GSM..., pb, bulk, H1975, matrix.tsv.gz]
        # Careful with index if parts vary, but this seems consistent in the file list
        try:
            cell_line_idx = parts.index('bulk') + 1
            cell_line = parts[cell_line_idx]
        except (ValueError, IndexError):
            print(f"Warning: Could not parse cell line from {basename}. Using full name.")
            cell_line = basename.split('.')[0]
        
        print(f"Processing {cell_line} from {basename}...")
        
        # Read file
        # Columns: tname, len, num_reads
        df = pd.read_csv(f, sep='\t')
        
        # Extract ENSG ID from tname
        # tname format: ENST...|ENSG...|...
        # We need to handle potential format variations. 
        # Based on user check: ENST...|ENSG...|...
        def extract_ensg(tname):
            try:
                return tname.split('|')[1]
            except IndexError:
                return tname # Fallback
        
        df['gene_id'] = df['tname'].apply(extract_ensg)
        
        # Aggregate by gene_id (sum num_reads)
        gene_counts = df.groupby('gene_id')['num_reads'].sum()
        
        # Normalize to CPM
        total_counts = gene_counts.sum()
        if total_counts > 0:
            cpm = (gene_counts / total_counts) * 1e6
        else:
            cpm = gene_counts
        
        bulk_profiles[cell_line] = cpm
        
    # Combine into a single DataFrame
    bulk_df = pd.DataFrame(bulk_profiles)
    
    # Fill NaNs with 0
    bulk_df = bulk_df.fillna(0)
    
    # Log1p transform
    bulk_df = np.log1p(bulk_df)
    
    return bulk_df

def main():
    # Datasets to process
    datasets = [
        ("sc", "data/longbench_data/sc/GSM9135508_pb_sc_genes_count_matrix.tsv.gz"),
        ("sn", "data/longbench_data/sn/GSM9135509_pb_sn_genes_count_matrix.tsv.gz")
    ]
    
    bulk_dir = 'data/longbench_data/bulk'
    output_dir = 'outputs'
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Loading Bulk Data...")
    bulk_df = load_bulk_data(bulk_dir)
    print(f"Loaded bulk profiles for {len(bulk_df.columns)} cell lines.")

    for name, path in datasets:
        print(f"\nProcessing {name.upper()} data from {path}...")
        output_file = os.path.join(output_dir, f'{name}_cell_line_labels.csv')

        print(f"Loading {name.upper()} Data...")
        adata = sc.read_csv(path).transpose()
        print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes.")
        
        # Basic preprocessing
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        
        # Intersect genes
        common_genes = list(set(adata.var_names) & set(bulk_df.index))
        print(f"Number of common genes: {len(common_genes)}")
        
        if len(common_genes) == 0:
            print(f"Error: No common genes found for {name} data.")
            continue

        # Subset
        sc_subset = adata[:, common_genes]
        sc_data = sc_subset.X
        
        if hasattr(sc_data, "toarray"):
            sc_data = sc_data.toarray()
            
        bulk_data_subset = bulk_df.loc[common_genes].values.T # Shape: (n_cell_lines, n_genes)
        
        print("Calculating correlations...")
        dists = cdist(sc_data, bulk_data_subset, metric='correlation')
        corrs = 1 - dists
        
        # Assign labels
        best_idx = np.argmax(corrs, axis=1)
        max_corrs = np.max(corrs, axis=1)
        
        cell_lines = bulk_df.columns
        assigned_labels = cell_lines[best_idx]
        
        # Create result DataFrame
        results = pd.DataFrame({
            'cell_id': adata.obs_names,
            'assigned_cell_line': assigned_labels,
            'correlation': max_corrs
        })
        
        # Print stats
        print(f"\nAssignment Stats ({name.upper()}):")
        print(results['assigned_cell_line'].value_counts())
        
        # Save
        results.to_csv(output_file, index=False)
        print(f"Saved {name.upper()} assignments to {output_file}")

if __name__ == "__main__":
    main()
