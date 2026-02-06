cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: hubmap/multi-data-products
baseCommand: /opt/upload.py

inputs:
    mudata_raw:
        type: File
        doc: The raw h5ad file
        inputBinding:
            position: 0
    
    muon_processed:
        type: File
        doc: processed data product
        inputBinding: 
            position: 1

    final_metadata_json:
        type: File
        doc: data product metadata json
        inputBinding:
            position: 2
    
    joint_embedding:
        type: File
        doc: combined umap (atac and rna)
        inputBinding:
            position: 3
    
    access_key_id:
        type: string?
        doc: AWS access key id
        inputBinding:
            position: 4
    
    secret_access_key:
        type: string
        doc: AWS secret access key
        inputBinding:
            position: 5

outputs: 
    finished_text:
        type: File
        outputBinding:
            glob: "*.txt"