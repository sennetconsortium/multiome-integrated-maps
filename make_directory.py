#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import requests
import subprocess


def main(uuids_file: Path, tissue: str):
    uuids = pd.read_csv(uuids_file, sep="\t")["uuid"].dropna()
    h5ads_base_directory = Path(f"{tissue}_h5mus")
    h5ads_base_directory.mkdir(exist_ok=True)

    for uuid in uuids:
        h5ads_directory = h5ads_base_directory / uuid
        h5ads_directory.mkdir(parents=True, exist_ok=True)
        out_file = h5ads_directory / "mudata_raw.h5mu"
        url = f"https://assets.api.sennetconsortium.org/{uuid}/mudata_raw.h5mu"
        # Check if file exists and downlaod
        try:
            head = requests.head(url, allow_redirects=True)
            if head.status_code != 200:
                print(f"File for {uuid} does not exist on server (status {head.status_code}).")
                continue
        except requests.RequestException as e:
            print(f"Error checking {uuid}: {e}")
            continue
        try:
            print(f"Downloading {uuid} ...")
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(out_file, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            print(f"Downloaded {uuid} successfully.")
        except requests.RequestException as e:
            print(f"Failed to download {uuid}: {e}")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("uuids_file", type=Path)
    p.add_argument("tissue", type=str)

    args = p.parse_args()

    main(args.uuids_file, args.tissue)
