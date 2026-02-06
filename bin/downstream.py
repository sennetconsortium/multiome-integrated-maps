#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import json
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import os
import scanpy as sc
from muon import prot as pt

from plot_utils import new_plot


def updata_metadata(metadata_json, cell_count, file_size):
    with open(metadata_json, 'r') as json_file:
        metadata = json.load(json_file)
    metadata["Processed Total Cell Count"] = cell_count
    metadata["Processed File Size"] = file_size
    uuid = metadata["Data Product UUID"]
    with open(f"{uuid}.json", "w") as outfile:
        json.dump(metadata, outfile)



def main(mudata_raw: Path, tissue: str, metadata: Path):
    # Read the raw mudata
    expr = mu.read(str(mudata_raw))

    expr.obs['num_genes_rna'] = (expr['rna'].X > 0).sum(axis=1)
    mu.pp.filter_obs(expr, 'num_genes_rna', lambda x: x > 200)
    rna_expr = expr["rna"]
    rna_expr.X = rna_expr.layers["spliced"]
    print(rna_expr)
    atac_cbg_expr = expr["atac_cell_by_gene"]
    atac_cbb_expr = expr["atac_cell_by_bin"]
    print("Construct muon object with RNA and ATAC cell-by-gene data.")
    mdata_raw = mu.MuData({"rna": rna_expr, "atac_cbg": atac_cbg_expr})
    print(mdata_raw)
    original_obs = expr.obs.copy()
    uuid = expr.uns["uuid"]
    original_uns = expr.uns.copy()

    mdata_rna = mdata_raw.mod["rna"]
    mdata_ataccbg = mdata_raw.mod["atac_cbg"]

    ## RNA normalization
    print("Normalizing RNA data...")
    sc.pp.normalize_total(mdata_rna, target_sum=1e4)
    sc.pp.log1p(mdata_rna)

    # ATAC TFIDF normalization
    mu.atac.pp.tfidf(mdata_ataccbg, scale_factor=1e4)

    # find highly variable genes for ATAC
    sc.pp.highly_variable_genes(mdata_ataccbg, min_mean=0.02, max_mean=4, min_disp=0.5)
    print("Found", np.sum(mdata_ataccbg.var.highly_variable), "highly variable genes in ATAC-seq.")

    ## Downstream analysis for Atac
    print("Performing downstream analysis for ATAC...")
    sc.tl.pca(mdata_ataccbg)
    sc.pp.neighbors(mdata_ataccbg)

    for axis in [0, 1]:
        neighbor_counts = np.array((mdata_ataccbg.obsp["distances"] > 0).sum(axis=axis)).flatten()
        mdata_ataccbg.obs.loc[:, f"neighbor_counts_ax{axis}"] = neighbor_counts
        cells_with_neighbors = (neighbor_counts > 0).astype(int)
        mdata_ataccbg.obs.loc[:, f"has_neighbors_ax{axis}"] = cells_with_neighbors
    mdata_ataccbg.obs.loc[:, "has_neighbors"] = mdata_ataccbg.obs.loc[:, "has_neighbors_ax1"]
    mu.pp.filter_obs(mdata_ataccbg, "has_neighbors", lambda x: x > 0)

    sc.tl.leiden(mdata_ataccbg)
    sc.tl.umap(mdata_ataccbg)
    with new_plot():
        sc.pl.umap(mdata_ataccbg, color="leiden", legend_loc="on data")
        plt.savefig(f"{uuid}_leiden_cluster_atac.png")

    ## Downstream analysis for RNA
    print("Performing downstream analysis for RNA...")

    # find highly variable genes for RNA-seq
    sc.pp.highly_variable_genes(mdata_rna, min_mean=0.02, max_mean=4, min_disp=0.5)
    print("Found", np.sum(mdata_rna.var.highly_variable), "highly variable genes in RNA-seq.")

    # scaling the data
    mdata_rna.raw = mdata_rna
    print("Scaling the log-normalised counts to zero mean and unit variance...")
    sc.pp.scale(mdata_rna, max_value=10)

    # PCA and neighborhood graph
    print("Performing PCA...")
    sc.tl.pca(mdata_rna)
    print("Constructing neighborhood graph...")
    sc.pp.neighbors(mdata_rna)

    for axis in [0, 1]:
        neighbor_counts = np.array((mdata_rna.obsp["distances"] > 0).sum(axis=axis)).flatten()
        mdata_rna.obs.loc[:, f"neighbor_counts_ax{axis}"] = neighbor_counts
        cells_with_neighbors = (neighbor_counts > 0).astype(int)
        mdata_rna.obs.loc[:, f"has_neighbors_ax{axis}"] = cells_with_neighbors
    mdata_rna.obs.loc[:, "has_neighbors"] = mdata_rna.obs.loc[:, "has_neighbors_ax1"]
    mu.pp.filter_obs(mdata_rna, "has_neighbors", lambda x: x > 0)

    # Clustering
    print("Performing leiden clustering with the computed neighbourhood graph...")
    sc.tl.leiden(mdata_rna)
    sc.tl.umap(mdata_rna)
    with new_plot():
        sc.pl.umap(mdata_rna, color="leiden", legend_loc="on data")
        plt.savefig(f"{uuid}_leiden_cluster_rna.png")

    ## Multi-omics integration
    mdata_raw.update()
    mu.pp.intersect_obs(mdata_raw)
    
    # # Merge the original .obs back with the modified .obs
    # # original_obs_filtered = original_obs.loc[mdata_raw.obs.index]

    print(mdata_raw.obs_keys())
    print(mdata_raw)
    mdata_raw.write("multiome_normalized.h5mu")
    
    ## Multi-omics factor analysis
    mu.tl.mofa(mdata_raw, outfile="multiome_mofa.hdf5", n_factors=30)

    # Multiplex clustering
    sc.pp.neighbors(mdata_raw["rna"])
    sc.pp.neighbors(mdata_raw["atac_cbg"])
    print(mdata_raw)
    print(mdata_raw.obs_keys())
    sc.pp.neighbors(mdata_raw, use_rep="X_mofa", key_added="mofa")
    sc.tl.umap(mdata_raw, neighbors_key="mofa")
    sc.tl.leiden(mdata_raw, resolution=1.0, neighbors_key="mofa", key_added="leiden_wnn")

    with new_plot():
        sc.pl.umap(mdata_raw, color="leiden_wnn", legend_loc="on data")
        plt.savefig(f"{uuid}_leiden_cluster_combined.png")

    # Add the cell-by-bin data back for output
    mdata_raw.mod["atac_cbb"] = atac_cbb_expr
    mdata_raw_copy = mdata_raw.copy()
    print(mdata_raw_copy.obs_keys())
    for column in original_obs.columns:
        mdata_raw_copy.obs[column] = original_obs[column]
    columns_to_keep = ['hubmap_id', 'age', 'sex', 'height', 'weight', 'bmi', 'cause_of_death', 'race', 'barcode', 'dataset', 'cell_id', 'num_genes_rna', 'leiden_wnn', "tissue"]
    mdata_raw_copy.obs = mdata_raw_copy.obs[columns_to_keep]
    mdata_raw_copy.obs["cell_id"] = mdata_raw_copy.obs["cell_id"].astype(str)
    mdata_raw_copy.var["highly_variable"] = mdata_raw_copy.var["highly_variable"].astype(bool)
    print(mdata_raw_copy)
    print(mdata_raw_copy.obs_keys())
    print(mdata_raw.var["highly_variable"])
    for key in original_uns:
        mdata_raw_copy.uns[key] = original_uns[key]

    mdata_raw_copy.write(f"{tissue}_processed.h5mu")

    processed_cell_count = mdata_raw_copy.obs.shape[0]
    processed_file_size = os.path.getsize(f"{tissue}_processed.h5mu")
    updata_metadata(metadata, processed_cell_count, processed_file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--mudata_raw", type=Path)
    p.add_argument("--tissue", type=str)
    p.add_argument("--metadata_json", type=Path)
    args = p.parse_args()

    main(args.mudata_raw, args.tissue, args.metadata_json)
