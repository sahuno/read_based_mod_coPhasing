## USAGE:
# 1. Run specific rule
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/GET_CLAIR_MODELS/get_basecall_n_clairs_models.smk \
# results/summary/complete_summary.csv --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" --keep-going -np

# --config samples_file=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/GET_CLAIR_MODELS/config.yaml \

import os
import subprocess
import pandas as pd
import re
from pathlib import Path
from packaging.version import Version, InvalidVersion
import yaml

# Configuration
configfile: "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/GET_CLAIR_MODELS/config.yaml"

# Load samples from YAML
with open(config["samples_file"], 'r') as f:
    SAMPLES = yaml.safe_load(f)["SAMPLES"]

# Helper functions
def get_basecall_model(bam_path):
    """Extract basecall_model from BAM @RG header"""
    hdr = subprocess.check_output(["samtools", "view", "-H", bam_path]).decode()
    for line in hdr.splitlines():
        if line.startswith("@RG") and "DS:basecall_model=" in line:
            m = re.search(r"DS:basecall_model=([^ \t]+)", line)
            if m:
                return m.group(1)
    return None

def derive_clairs_model(
    bc_model,
    default_cancer_model="ont_r10_dorado_sup_5khz_ssrs",
    user_cancer_model=None,
    sup_4khz_model="ont_r10_dorado_sup_4khz",
    hac_5khz_model="ont_r10_dorado_hac_5khz",
    r9_sup_model="ont_r9_guppy"
):
    """
    Given a basecall_model string, decide which ClairS model to use:
    
    1. If bc_model contains both "r9.4.1" and "_sup":
       Return r9_sup_model.
    
    2. If bc_model contains both "r10.4.1" and "_sup":
       a) If version == 4.1.0, return sup_4khz_model.
       b) If version > 4.1.0:
          - If user_cancer_model is provided, return user_cancer_model.
          - Otherwise, return default_cancer_model.
       c) Otherwise (< 4.1.0), return None.
    
    3. If bc_model contains both "r10.4.1" and "_hac":
       a) If version == 4.2.0, return hac_5khz_model.
       b) Otherwise, return None.
    
    4. In all other cases, return None.
    """
    if not bc_model:
        return None

    # Split off the "@v<version>" suffix
    try:
        name, ver_str = bc_model.split("@v", 1)
    except ValueError:
        return None

    name = name.lower()

    # 1) Handle R9.4.1 SUP cases first
    if "r9.4.1" in name and "_sup" in name:
        return r9_sup_model
    
    # 2) Handle R10.4.1 "sup" cases
    if "r10.4.1" in name and "_sup" in name:
        try:
            version = Version(ver_str)
        except InvalidVersion:
            return None
        # Exact match v4.1.0 → use the 4 kHz SUP model
        if version == Version("4.1.0"):
            return sup_4khz_model
        # Strictly greater than v4.1.0 → prefer user-provided, else default
        if version > Version("4.1.0"):
            return user_cancer_model or default_cancer_model
        # Otherwise (< 4.1.0) → no match for SUP
        return None
    
    # 3) Handle R10.4.1 "hac" cases
    if "r10.4.1" in name and "_hac" in name:
        try:
            version = Version(ver_str)
        except InvalidVersion:
            return None
        # Exact match v4.2.0 → use the 5 kHz HAC model
        if version == Version("4.2.0"):
            return hac_5khz_model
        # All other versions → no match
        return None

    # 4) No matching flow cell / accuracy pattern
    return None

# Generate all sample-type combinations
def get_all_sample_files():
    """Get all BAM files from the samples YAML"""
    files = []
    for patient, sample_types in SAMPLES.items():
        for sample_type, samples in sample_types.items():
            for sample_id, bam_path in samples.items():
                files.append({
                    'patient': patient,
                    'sample_type': sample_type,
                    'sample_id': sample_id,
                    'bam_path': bam_path
                })
    return files

def validate_samples():
    """Validate sample structure and paths"""
    import os
    issues = []
    
    for sample in ALL_SAMPLES:
        # Check if BAM file exists
        if not os.path.exists(sample['bam_path']):
            issues.append(f"BAM file not found: {sample['bam_path']}")
        
        # Check sample_id consistency
        if sample['sample_type'] == 'TUMOR' and '_T' not in sample['sample_id']:
            issues.append(f"Warning: TUMOR sample {sample['sample_id']} doesn't contain '_T'")
        elif sample['sample_type'] == 'NORMAL' and '_N' not in sample['sample_id']:
            issues.append(f"Warning: NORMAL sample {sample['sample_id']} doesn't contain '_N'")
    
    if issues:
        print("=== VALIDATION ISSUES ===")
        for issue in issues:
            print(f"  - {issue}")
        print("========================")
    
    return len(issues) == 0

def get_bam_path(patient, sample_type, sample_id):
    """Get BAM path for given wildcards with better error handling"""
    matching_samples = [
        sample for sample in ALL_SAMPLES 
        if sample['patient'] == patient 
        and sample['sample_type'] == sample_type 
        and sample['sample_id'] == sample_id
    ]
    
    if not matching_samples:
        available_samples = [
            f"{s['patient']}/{s['sample_type']}/{s['sample_id']}" 
            for s in ALL_SAMPLES
        ]
        raise ValueError(
            f"No matching sample found for {patient}/{sample_type}/{sample_id}.\n"
            f"Available samples:\n" + "\n".join(f"  - {s}" for s in available_samples)
        )
    
    if len(matching_samples) > 1:
        raise ValueError(f"Multiple samples found for {patient}/{sample_type}/{sample_id}")
    
    return matching_samples[0]['bam_path']

ALL_SAMPLES = get_all_sample_files()

# Validate samples
validate_samples()

# Debug: Print sample information
print("=== DEBUG: Sample Information ===")
for sample in ALL_SAMPLES:
    print(f"Patient: {sample['patient']}, Type: {sample['sample_type']}, ID: {sample['sample_id']}")
print(f"Total samples: {len(ALL_SAMPLES)}")
print("=================================")

# Define output files
def get_output_files():
    """Define all output files based on samples"""
    outputs = []
    for sample in ALL_SAMPLES:
        # Basecall model extraction outputs
        outputs.append(f"results/basecall_models/{sample['patient']}/{sample['sample_type']}/{sample['sample_id']}_basecall_model.txt")
        # ClairS model assignment outputs
        outputs.append(f"results/clairs_models/{sample['patient']}/{sample['sample_type']}/{sample['sample_id']}_clairs_model.txt")
    
    # Summary outputs
    outputs.extend([
        "results/summary/basecall_models_summary.csv",
        "results/summary/clairs_models_summary.csv",
        "results/summary/complete_summary.csv"
    ])
    return outputs

# Main rule
rule all:
    input:
        get_output_files()

# Rule to extract basecall model from BAM file
rule extract_basecall_model:
    input:
        bam = lambda wildcards: get_bam_path(wildcards.patient, wildcards.sample_type, wildcards.sample_id)
    output:
        "results/basecall_models/{patient}/{sample_type}/{sample_id}_basecall_model.txt"
    singularity:
        config.get("samtools_container", "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0")
    shell:
        """
        python3 -c "
import subprocess
import re
import sys

def get_basecall_model(bam_path):
    hdr = subprocess.check_output(['samtools', 'view', '-H', bam_path]).decode()
    for line in hdr.splitlines():
        if line.startswith('@RG') and 'DS:basecall_model=' in line:
            m = re.search(r'DS:basecall_model=([^ \\t]+)', line)
            if m:
                return m.group(1)
    return 'Not_Found'

basecall_model = get_basecall_model('{input.bam}')
with open('{output}', 'w') as f:
    f.write(basecall_model)
"
        """

# Rule to derive ClairS model from basecall model
rule derive_clairs_model:
    input:
        "results/basecall_models/{patient}/{sample_type}/{sample_id}_basecall_model.txt"
    output:
        "results/clairs_models/{patient}/{sample_type}/{sample_id}_clairs_model.txt"
    params:
        default_cancer_model = config.get("default_cancer_model", "ont_r10_dorado_sup_5khz_ssrs"),
        user_cancer_model = config.get("user_cancer_model", None),
        sup_4khz_model = config.get("sup_4khz_model", "ont_r10_dorado_sup_4khz"),
        hac_5khz_model = config.get("hac_5khz_model", "ont_r10_dorado_hac_5khz"),
        r9_sup_model = config.get("r9_sup_model", "ont_r9_guppy")
    singularity:
        config.get("python_container", "docker://python:3.11-slim")
    shell:
        """
        pip install packaging && python3 -c "
from packaging.version import Version, InvalidVersion

def derive_clairs_model(
    bc_model,
    default_cancer_model='{params.default_cancer_model}',
    user_cancer_model={params.user_cancer_model},
    sup_4khz_model='{params.sup_4khz_model}',
    hac_5khz_model='{params.hac_5khz_model}',
    r9_sup_model='{params.r9_sup_model}'
):
    if not bc_model or bc_model == 'Not_Found':
        return 'No_Model_Available'

    try:
        name, ver_str = bc_model.split('@v', 1)
    except ValueError:
        return 'No_Model_Available'

    name = name.lower()

    # 1) Handle R9.4.1 SUP cases first (no version checking needed)
    if 'r9.4.1' in name and '_sup' in name:
        return r9_sup_model

    # 2) Handle R10.4.1 SUP cases
    if 'r10.4.1' in name and '_sup' in name:
        try:
            version = Version(ver_str)
        except InvalidVersion:
            return 'No_Model_Available'
        if version == Version('4.1.0'):
            return sup_4khz_model
        if version > Version('4.1.0'):
            return user_cancer_model or default_cancer_model
        return 'No_Model_Available'
    
    # 3) Handle R10.4.1 HAC cases
    if 'r10.4.1' in name and '_hac' in name:
        try:
            version = Version(ver_str)
        except InvalidVersion:
            return 'No_Model_Available'
        if version == Version('4.2.0'):
            return hac_5khz_model
        return 'No_Model_Available'

    return 'No_Model_Available'

with open('{input}', 'r') as f:
    basecall_model = f.read().strip()

clairs_model = derive_clairs_model(basecall_model)

with open('{output}', 'w') as f:
    f.write(clairs_model)
"
        """

# Rule to create basecall models summary
rule summarize_basecall_models:
    input:
        [f"results/basecall_models/{s['patient']}/{s['sample_type']}/{s['sample_id']}_basecall_model.txt" for s in ALL_SAMPLES]
    output:
        "results/summary/basecall_models_summary.csv"
    singularity:
        config.get("pandas_container", "docker://quay.io/biocontainers/pandas:1.5.3")
    run:
        import pandas as pd
        
        data = []
        for sample in ALL_SAMPLES:
            basecall_file = f"results/basecall_models/{sample['patient']}/{sample['sample_type']}/{sample['sample_id']}_basecall_model.txt"
            
            with open(basecall_file, 'r') as f:
                basecall_model = f.read().strip()
            
            data.append({
                'patient': sample['patient'],
                'sample_type': sample['sample_type'], 
                'sample_id': sample['sample_id'],
                'bam_path': sample['bam_path'],
                'basecall_model': basecall_model
            })
        
        df = pd.DataFrame(data)
        df.to_csv(output[0], index=False)

# Rule to create ClairS models summary
rule summarize_clairs_models:
    input:
        [f"results/clairs_models/{s['patient']}/{s['sample_type']}/{s['sample_id']}_clairs_model.txt" for s in ALL_SAMPLES]
    output:
        "results/summary/clairs_models_summary.csv"
    singularity:
        config.get("pandas_container", "docker://quay.io/biocontainers/pandas:1.5.3")
    run:
        import pandas as pd
        
        data = []
        for sample in ALL_SAMPLES:
            clairs_file = f"results/clairs_models/{sample['patient']}/{sample['sample_type']}/{sample['sample_id']}_clairs_model.txt"
            
            with open(clairs_file, 'r') as f:
                clairs_model = f.read().strip()
            
            data.append({
                'patient': sample['patient'],
                'sample_type': sample['sample_type'],
                'sample_id': sample['sample_id'],
                'bam_path': sample['bam_path'],
                'clairs_model': clairs_model
            })
        
        df = pd.DataFrame(data)
        df.to_csv(output[0], index=False)

# Rule to create complete summary
rule complete_summary:
    input:
        basecall_summary = "results/summary/basecall_models_summary.csv",
        clairs_summary = "results/summary/clairs_models_summary.csv"
    output:
        "results/summary/complete_summary.csv"
    singularity:
        config.get("pandas_container", "docker://quay.io/biocontainers/pandas:1.5.3")
    run:
        import pandas as pd
        
        basecall_df = pd.read_csv(input.basecall_summary)
        clairs_df = pd.read_csv(input.clairs_summary)
        
        # Merge the dataframes
        merge_cols = ['patient', 'sample_type', 'sample_id', 'bam_path']
        complete_df = pd.merge(basecall_df, clairs_df[merge_cols + ['clairs_model']], on=merge_cols)
        
        complete_df.to_csv(output[0], index=False)