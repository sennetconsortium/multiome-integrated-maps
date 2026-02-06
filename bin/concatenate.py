#!/usr/bin/env python3

import json
from argparse import ArgumentParser
from collections import defaultdict
from datetime import datetime
from os import fspath, walk, listdir
from pathlib import Path
from typing import Dict, Tuple

import anndata
import mudata as md
import muon as mu
import numpy as np
import os
import pandas as pd
import requests
import uuid
import yaml


def get_tissue_type(dataset: str) -> str:
    organ_dict = yaml.load(open("/opt/organ_types.yaml"), Loader=yaml.BaseLoader)
    organ_code = requests.get(
        f"https://entity.api.hubmapconsortium.org/dataset/{dataset}/organs/"
    )
    organ_name = organ_dict[organ_code]
    return organ_name.replace(" (Left)", "").replace(" (Right)", "")


def convert_tissue_code(tissue_code):
    with open("/opt/organ_types.yaml", 'r') as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
    tissue_name = data.get(tissue_code)['description']
    return tissue_name


def find_files(directory, patterns):
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            for pattern in patterns:
                if filepath.match(pattern):
                    return filepath
                

def find_file_pairs(directory):
    raw_mdata_patterns = ["mudata_raw.h5mu"]
    processed_mdata_patterns = ["secondary_analysis.h5mu"]
    raw_mdata_file = find_files(directory, raw_mdata_patterns)
    processed_mdata_file = find_files(directory, processed_mdata_patterns)
    return (raw_mdata_file, processed_mdata_file)


def make_unique_barcodes(mdata_file, tissue_type: str = None):
    data_set_dir = fspath(mdata_file.parent.stem)
    tissue_type = tissue_type if tissue_type else get_tissue_type(data_set_dir)
    
    # Load MuData object
    mdata = mu.read_h5mu(mdata_file)
    mdata_copy = mdata.copy()

    # Make index unique to this dataset
    mdata_copy.obs["barcode"] = mdata.obs.index
    mdata_copy.obs["barcode"] = mdata_copy.obs["barcode"].str.replace("BAM_data#", "", regex=False)
    mdata_copy.obs["dataset"] = data_set_dir
    mdata_copy.obs["tissue"] = tissue_type
    
    cell_ids_list = [
        "-".join([data_set_dir, barcode]) for barcode in mdata_copy.obs["barcode"]
    ]
    mdata_copy.obs["cell_id"] = pd.Series(
        cell_ids_list, index=mdata_copy.obs.index, dtype=str
    )
    
    # Iterate through each modality in mdata, accessing by key
    print(list(mdata_copy.mod.keys()))
    for key in mdata_copy.mod.keys():
        mod_data = mdata_copy[key]
        print(key)
        print(mod_data.obs_keys())

        mod_data.obs["barcode"] = mod_data.obs.index
        if mod_data.obs["barcode"].str.contains("BAM_data#").any():
            mod_data.obs["barcode"] = mod_data.obs["barcode"].str.replace("BAM_data#", "", regex=False)

        cell_ids_list = [
            f"{data_set_dir}-{barcode}" for barcode in mod_data.obs["barcode"]
        ]
        mod_data.obs["cell_id"] = pd.Series(cell_ids_list, index=mod_data.obs.index, dtype=str)
        mod_data.obs.set_index("cell_id", drop=True, inplace=True)
        mod_data.obs["tissue"] = tissue_type

        print(mod_data.obs)
        mdata_copy.mod[key] = mod_data

    mdata_copy.obs["tissue"] = tissue_type
    
    return mdata_copy


def concatenate_modalities(raw_mdatas, modality_keys):
    concatenated_anndata_dict = {}

    # Iterate over each modality key
    for key in modality_keys:
        # Extract all Anndata objects for this modality across all muon objects
        anndata_list = [mdata.mod[key] for mdata in raw_mdatas if key in mdata.mod]
        
        # Concatenate the Anndata objects for this modality
        concatenated_anndata_dict[key] = anndata.concat(anndata_list, join='outer')

    return concatenated_anndata_dict


def concatenate_obs(raw_mdatas):
    obs_list = [mdata.obs for mdata in raw_mdatas]
    concatenated_obs = pd.concat(obs_list, axis=0)
    return concatenated_obs


def concat_mudatas(concatenated_anndata_dict, mudata_obs):
    # Create a new mudata object from the concatenated anndata objects
    new_mudata = mu.MuData(concatenated_anndata_dict)
    new_mudata.obs = mudata_obs
    return new_mudata


def create_json(tissue, data_product_uuid, creation_time, uuids, hbmids, cell_count, file_size):
    bucket_url = f"https://hubmap-data-products.s3.amazonaws.com/{data_product_uuid}/"
    metadata = {
        "Data Product UUID": data_product_uuid,
        "Tissue": convert_tissue_code(tissue),
        "Assay": "10X Multiome",
        "URL": bucket_url + f"{tissue}.h5mu",
        "Creation Time": creation_time,
        "Dataset UUIDs": uuids,
        "Dataset HBMIDs": hbmids,
        "Raw Total Cell Count": cell_count,
        "Raw File Size": file_size
    }
    print("Writing metadata json")
    with open(f"{data_product_uuid}.json", "w") as outfile:
        json.dump(metadata, outfile)


def annotate_mudata(mdata, uuids_df):
    merged = uuids_df.merge(mdata.obs, left_on="uuid", right_on="dataset", how="inner")
    merged = merged.set_index(mdata.obs.index)
    merged = merged.drop(columns=["Unnamed: 0"])
    merged = merged.fillna(np.nan)
    merged["age"] = pd.to_numeric(merged["age"])
    return merged


def main(data_directory: Path, uuids_file: Path, tissue: str = None):
    output_file_name = f"{tissue}_raw" if tissue else "multiome"
    uuids_df = pd.read_csv(uuids_file, sep="\t", dtype=str)
    uuids_list = uuids_df["uuid"].to_list()
    hbmids_list = uuids_df["hubmap_id"].to_list()
    directories = [data_directory / Path(uuid) for uuid in uuids_df["uuid"]]
    # Load files
    file_pairs = [find_file_pairs(directory) for directory in directories if len(listdir(directory))>1]
    print("Annotating objects")
    raw_mdatas = [
        make_unique_barcodes(file_pair[0], tissue)
        for file_pair in file_pairs
    ]


    # processed_mdatas = [     
    #     load_mudata(file_pair[1])
    #     for file_pair in file_pairs
    # ]

    print("Concatenating objects")
    modality_keys = ["rna", "atac_cell_by_bin", "atac_cell_by_gene"]
    concatenated_anndata = concatenate_modalities(raw_mdatas, modality_keys)
    concat_obs = concatenate_obs(raw_mdatas)
    raw_mdata_concat = concat_mudatas(concatenated_anndata, concat_obs)
    raw_mdata_concat.obs = annotate_mudata(raw_mdata_concat, uuids_df)
    columns_to_keep = [
        "hubmap_id", "age", "sex", "height", "weight", "bmi", "cause_of_death", "race", "barcode", "dataset", "cell_id", "tissue"
    ]
    raw_mdata_concat.obs = raw_mdata_concat.obs[columns_to_keep]


    creation_time = str(datetime.now())
    data_product_uuid = str(uuid.uuid4())
    total_cell_count = raw_mdata_concat.obs.shape[0]
    raw_mdata_concat.uns["creation_data_time"] = creation_time
    raw_mdata_concat.uns["datasets"] = hbmids_list
    raw_mdata_concat.uns["uuid"] = data_product_uuid
    raw_mdata_concat.write(f"{output_file_name}.h5mu")
    print(raw_mdata_concat.obs)
    print(raw_mdata_concat.obs_keys())
    print(raw_mdata_concat.mod["atac_cell_by_bin"].obs)
    print(raw_mdata_concat.mod["atac_cell_by_gene"].obs)
    print(raw_mdata_concat.mod["rna"].obs)
    file_size = os.path.getsize(f"{output_file_name}.h5mu")
    create_json(tissue, data_product_uuid, creation_time, uuids_list, hbmids_list, total_cell_count, file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("data_directory", type=Path)
    p.add_argument("uuids_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.data_directory, args.uuids_file, args.tissue)