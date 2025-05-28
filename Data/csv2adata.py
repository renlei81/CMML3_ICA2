import pandas as pd
import numpy as np
import scipy.sparse as sparse
import scanpy as sc

def make_index_unique(index, description="index"):
    if index.is_unique:
        return index.astype(str)
    print(f"⚠️ Warning: {description} not unique, deduplicating.")
    return pd.Index(pd.io.parsers.ParserBase({'names': index})._maybe_dedup_names(index.astype(str)))

def csv_to_adata_from_counts(file_path):
    print(f"\n📂 Loading file: {file_path}")
    df = pd.read_csv(file_path, index_col=0)
    print(f"🧾 Raw shape: {df.shape}")

    # 如果列数远大于行数，说明原始格式是 gene × cell，需要转置
    if df.shape[1] > df.shape[0] * 2:
        df = df.T
        print(f"🔁 Transposed matrix to shape: {df.shape} (cells × genes)")
    else:
        print(f"✅ Using matrix as-is: {df.shape} (cells × genes)")

    df.index = make_index_unique(df.index, "cell names")
    df.columns = make_index_unique(df.columns, "gene names")

    adata = sc.AnnData(X=sparse.csr_matrix(df.values, dtype=np.float32),
                       obs=pd.DataFrame(index=df.index),
                       var=pd.DataFrame(index=df.columns))
    print(f"✅ Created AnnData: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


import os
from csv2adata import csv_to_adata_from_counts

input_dir = "/root/Desktop/my_pan/workspace/Data/lung_csv"
output_dir = "/root/Desktop/my_pan/workspace/Data/lung_h5ad_output"
os.makedirs(output_dir, exist_ok=True)

for fname in os.listdir(input_dir):
    if not fname.endswith(".csv"):
        continue

    input_path = os.path.join(input_dir, fname)
    sample_name = fname.replace(".csv", "")
    output_path = os.path.join(output_dir, f"{sample_name}.h5ad")

    try:
        adata = csv_to_adata_from_counts(input_path)
        adata.write_h5ad(output_path)
        print(f"💾 Saved to {output_path}")
    except Exception as e:
        print(f"❌ Failed to process {fname}: {e}")
