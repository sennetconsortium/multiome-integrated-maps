#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for concatenating sc-ATAC-seq datasets into MuData object

requirements:
  ScatterFeatureRequirement: {}

inputs: 
    data_directory:
        label: "Path to directory containing cell by gene and cell by bin files"
        type: Directory
    
    uuids_file:
        label: "Path to a file containing a list of uuids and other metadata for the dataset to be indexed"
        type: File
    
    tissue:
        label: "Two letter tissue type code"
        type: string?
      
    access_key_id:
        label: "AWS access key id"
        type: string
    
    secret_access_key:
        label: "AWS secret access key"
        type: string

outputs:
    mudata_raw:
        type: File
        outputSource: concatenate/mudata_raw
    
    muon_processed: 
        type: File
        outputSource: downstream/muon_processed
    
    joint_embedding:
        type: File
        outputSource: downstream/joint_embedding
    
    final_metadata_json:
        type: File
        outputSource: downstream/final_metadata_json

steps:

    - id: concatenate
      in: 
        - id: data_directory
          source: data_directory
        - id: uuids_file
          source: uuids_file
        - id: tissue
          source: tissue
    
      out:
        - mudata_raw
        - metadata_json
      run: steps/concatenate.cwl
      label: "Concatenates h5ad files in directory"

    - id: downstream
      in:
        - id: mudata_raw
          source: concatenate/mudata_raw
        - id: tissue
          source: tissue

        - id: metadata_json
          source: concatenate/metadata_json
      
      out:
        - muon_processed
        - mofa_out
        - joint_embedding
        - rna_embedding
        - atac_embedding
        - final_metadata_json
      run: steps/downstream.cwl


    - id: upload
      in: 
        - id: mudata_raw
          source: concatenate/mudata_raw
        - id: muon_processed
          source: downstream/muon_processed
        - id: final_metadata_json
          source: downstream/final_metadata_json
        - id: joint_embedding
          source: downstream/joint_embedding
        - id: access_key_id
          source: access_key_id
        - id: secret_access_key
          source: secret_access_key
    
      out:
        - finished_text
      run: steps/upload.cwl
      label: "Uploads the pipeline outputs to s3"
      