#!/usr/bin/env python3
import json
import os
from argparse import ArgumentParser
from pathlib import Path


def get_uuid(data_product_metadata):
    with open(data_product_metadata, 'r') as json_file:
        metadata = json.load(json_file)
    uuid = metadata["Data Product UUID"]
    return uuid


def upload_to_ec2(metadata_json, uuid, ssh_key, umap):
    os.system(f"scp -i {ssh_key} {metadata_json} main_user@ec2-44-213-71-141.compute-1.amazonaws.com:/pipeline_outputs/{uuid}.json")
    os.system(f"scp -i {ssh_key} {umap} main_user@ec2-44-213-71-141.compute-1.amazonaws.com:/pipeline_outputs/{umap.name}")


def main(metadata_json, umap, ssh_key):
    uuid = get_uuid(metadata_json)
    upload_to_ec2(metadata_json, uuid, ssh_key, umap)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--data_product_metadata", type=Path)
    p.add_argument("--umap", type=Path)
    p.add_argument("--ssh_key", type=Path)
    args = p.parse_args()

    main(args.data_product_metadata, args.umap, args.ssh_key)
