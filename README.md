# read_based_mod_coPhasing


step 1:
get clairS model automatically from the bam file
`wf_snakemake/MultiTumorNormalPair/GET_CLAIR_MODELS/get_basecall_n_clairs_models.smk`
input: samples.yaml 
output: .csv with sample names and clairS models

step 2.
run read based pahsing with sample sheet generated in step 1
```
wf_snakemake/MultiTumorNormalPair/read_based_coPhasing.smk
```