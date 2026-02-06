#!/usr/bin/env python3
from argparse import ArgumentParser

import pandas as pd
import requests
import yaml
import json

organ_types_yaml_file = "bin/organ_types.yaml"
organ_uberon_file = "bin/organs.json"


def get_uuids(organ_uberon: str, organism: str):
    params = {
        "status": "Published",
        "dataset_type": "10X Multiome",
    }

    if organ_uberon:
        params["origin_samples.organ"] = organ_uberon

    url = "https://search.api.sennetconsortium.org/param-search/datasets"
    response = requests.get(url, params=params)
    print("Request URL:", response.url)
    print("Response status code: ", response.status_code)

    # Handle a successful response
    if response.status_code == 200:
        return process_response(response, organism)

    # Handle 303 redirection
    elif response.status_code == 303:
        redirection_url = (
            response.text.strip()
        )  # Get redirection URL from response body
        print("Following redirection URL: ", redirection_url)

        # Make a request to the redirection URL
        redirection_response = requests.get(redirection_url)
        if redirection_response.status_code == 200:
            return process_response(redirection_response, organism)

    # Handle other error responses
    else:
        print(f"Error {response.status_code}: {response.text}")
        return [], [], []


def process_response(response, organism):
    data = response.json()
    items = data

    uuids = []
    sennet_ids = []
    donor_metadata_list = []

    for item in items:
        sources = item.get("sources", [])
        if not sources:
            continue
        source = sources[0]
        if source.get("source_type").lower() == organism.lower():
            uuids.append(item.get("uuid"))
            sennet_ids.append(item.get("sennet_id"))

            # Attempt to extract donor metadata
            metadata = item.get("sources")[0]
            donor_metadata_list.append(extract_donor_metadata(metadata))

    return uuids, sennet_ids, donor_metadata_list


def extract_donor_metadata(metadata):
    donor_info = {
        "age": None,
        "sex": None,
        "height": None,
        "weight": None,
        "bmi": None,
        "cause_of_death": None,
        "race": None,
        "social_history": None,
        "abo_blood_type": None,
        "mechanism_of_injury": None,
    }

    donor_metadata = metadata.get("mapped_metadata", {})

    for key in donor_metadata:
        if key == "abo_blood_group_system":
            donor_info["abo_blood_type"] = donor_metadata[key].get("value_display")
        elif key == "age":
            donor_info["age"] = donor_metadata[key].get("value_display")
        elif key == "body_mass_index":
            donor_info["bmi"] = donor_metadata[key].get("value_display")
        elif key == "cause_of_death":
            donor_info["cause_of_death"] = donor_metadata[key].get("value_display")
        elif key == "height":
            donor_info["height"] = donor_metadata[key].get("value_display")
        elif key == "mechanism_of_injury":
            donor_info["mechanism_of_injury"] = donor_metadata[key].get("value_display")
        elif key == "race":
            donor_info["race"] = donor_metadata[key].get("value_display")
        elif key == "sex":
            donor_info["sex"] = donor_metadata[key].get("value_display")
        elif key == "social_history":
            donor_info["social_history"] = donor_metadata[key].get("value_display")
        elif key == "weight":
            donor_info["weight"] = donor_metadata[key].get("value_display")

    return donor_info


def get_organ_uberon(organ_name):
    term = organ_name.lower().strip()
    with open(organ_uberon_file) as f:
        data = json.load(f)
    for entry in data:
        if entry.get("term", "").lower() == term:
            return entry["organ_uberon"]
    for entry in data:
        cat = entry.get("category")
        if cat and cat.get("term", "").lower() == term:
            return cat["organ_uberon"]
    return None


def main(tissue_type: str, organism:str):
    organ_dict = yaml.load(open(organ_types_yaml_file), Loader=yaml.BaseLoader)
    for key in organ_dict:
        organ_dict[key] = organ_dict[key]["description"]
    uberon_code = get_organ_uberon(tissue_type)
    uuids_list, sennet_ids_list, donor_metadata = get_uuids(uberon_code, organism)
    uuids_df = pd.DataFrame()
    uuids_df["uuid"] = pd.Series(uuids_list, dtype=str)
    uuids_df["sennet_id"] = pd.Series(sennet_ids_list, dtype=str)
    donor_metadata_df = pd.DataFrame(donor_metadata)
    result_df = pd.concat([uuids_df, donor_metadata_df], axis=1)
    key_for_tissue = [key for key, value in organ_dict.items() if value == tissue_type]
    if key_for_tissue:
        output_file_name = f"{key_for_tissue[0].lower()}.tsv"
    else:
        output_file_name = "rna.tsv"
    result_df['organism'] = organism
    result_df['tissue'] = tissue_type
    print(result_df)
    result_df.to_csv(output_file_name, sep="\t")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("tissue_type", type=str, nargs="?", help="Type of tissue (optional)")
    p.add_argument("organism", type=str)
    args = p.parse_args()

    main(args.tissue_type, args.organism)
