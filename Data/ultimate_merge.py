import scanpy as sc
import os
from anndata import AnnData
from typing import List, Optional

def safe_merge_and_fix(
    h5ad_paths: List[str],
    batch_names: Optional[List[str]] = None,
    output_path: str = "fixed_merged_data.h5ad"
):
    """
    Merge multiple .h5ad files and fix orientation + batch + obs_names.
    """

    print(f"📦 Merging {len(h5ad_paths)} files...")

    adatas = []

    for i, path in enumerate(h5ad_paths):
        print(f"🔹 Loading: {path}")
        adata = sc.read_h5ad(path)

        # Try to guess if .X is genes x cells (wrong)
        if adata.shape[0] < adata.shape[1]:
            print("Detected likely transposed matrix, fixing...")
            adata = AnnData(
                X=adata.X.T,
                obs=adata.var.copy(),
                var=adata.obs.copy()
            )

        # Ensure unique names
        adata.obs_names_make_unique()
        adata.var_names_make_unique()

        # Add batch info
        if batch_names:
            batch_name = batch_names[i]
        else:
            # Use filename as batch
            batch_name = os.path.splitext(os.path.basename(path))[0]

        adata.obs["batch"] = batch_name
        adata.obs["batch"] = adata.obs["batch"].astype(str)

        adatas.append(adata)

    # Merge all
    print("🧬 Concatenating all AnnData objects...")
    adata_merged = adatas[0].concatenate(adatas[1:], batch_key="tmp", index_unique=None)

    # Drop extra tmp column added by concatenate
    if "tmp" in adata_merged.obs:
        adata_merged.obs.drop(columns=["tmp"], inplace=True)

    adata_merged.obs_names_make_unique()
    adata_merged.var_names_make_unique()

    print(f"💾 Saving merged object to: {output_path}")
    adata_merged.write_h5ad(output_path)
    print("✅ Merge & fix complete!")

    return adata_merged
