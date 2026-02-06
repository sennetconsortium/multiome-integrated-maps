cwlVersion: v1.1
class: CommandLineTool
label: Downstream analysis for RNA and ATAC 
requirements:
  DockerRequirement:
      dockerPull: hubmap/multi-data-products:latest
baseCommand: /opt/downstream.py

inputs:
  mudata_raw:
    type: File
    inputBinding:
      position: 0
      prefix: "--mudata_raw"

  tissue:
    type: string
    inputBinding:
      position: 1
      prefix: "--tissue"
  
  metadata_json:
    type: File
    inputBinding:
      position: 2
      prefix: "--metadata_json"
outputs:
  muon_processed:
    type: File
    outputBinding:
      glob: "*_processed.h5mu"
  mofa_out:
    type: File
    outputBinding:
      glob: "multiome_mofa.hdf5"
  joint_embedding:
    type: File
    outputBinding:
      glob: "*leiden_cluster_combined.png"
  rna_embedding:
    type: File
    outputBinding:
      glob: "*leiden_cluster_rna.png"
  atac_embedding:
    type: File
    outputBinding:
      glob: "*leiden_cluster_atac.png"
  final_metadata_json:
    type: File
    outputBinding:
      glob: "*.json"
  