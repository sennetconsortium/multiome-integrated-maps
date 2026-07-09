#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: sennet/multi-maps:latest

inputs:
  muon_processed:
    type: File
    inputBinding:
      position: 0
  tissue:
    type: string
    inputBinding:
      position: 1
  metadata_json:
    type: File
    inputBinding:
      position: 2

outputs:
  annotated_mudata:
    type: File
    outputBinding:
      glob: '*_processed.h5mu'
  metadata_with_cell_types:
    type: File
    outputBinding:
      glob: "*.json"

baseCommand: ['python3', '/opt/pan_organ_azimuth.py']