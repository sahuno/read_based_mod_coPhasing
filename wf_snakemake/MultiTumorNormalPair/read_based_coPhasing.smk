import pandas as pd
import yaml
import os
import itertools

# snakemake --snakefile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/read_based_coPhasing.smk \
# --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/config/slurmMinimal \
# --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" --keep-going --rerun-incomplete -np


# ─── 0) LOAD CONFIGURATION ────────────────────────────────────────────────────
configfile: "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/config.yaml"

# Extract configuration values
SAMPLES_FILE = config["samples_file"]
SEED = config["misc"]["random_seed"]
OUTDIR = config["output_dir"]

# Container images
CLAIRS_IMG = config["containers"]["clairs"]
MODKIT_IMG = config["containers"]["modkit"]
ONT_TOOLS_IMG = config["containers"]["ont_tools"]
SIF_SNIFFLES = config["containers"]["sniffles"]
CONDAPREFIX = config["conda_prefix"]

# Reference files
REF = config["reference"]["genome"]
TANDEM_REPS = config["reference"]["tandem_repeats"]
CONTIG = config["reference"]["contig"]
GENOMESIZES = config["reference"]["genome_sizes"]

# Computational resources
THREADS = config["resources"]["threads"]["default"]
THREADS_SV = config["resources"]["threads"]["sv_calling"]

# Output directories
SV_OUTDIR = config["sv_output_dir"]

# Variant calling parameters
NORMAL_FILTER_EXPR = config["variant_calling"]["sv_filters"]["normal_filter"]
TUMOR_FILTER_EXPR = config["variant_calling"]["sv_filters"]["tumor_filter"]

# ─── LOAD AND PROCESS SAMPLES ─────────────────────────────────────────────────
# Read sample metadata
samples_df = pd.read_csv(SAMPLES_FILE)

# Create tumor-normal pairs for each patient
def create_tn_pairs(df):
    """Create all tumor-normal pairs for each patient"""
    pairs = []
    for patient in df['patient'].unique():
        patient_df = df[df['patient'] == patient]
        tumors = patient_df[patient_df['sample_type'] == 'TUMOR']
        normals = patient_df[patient_df['sample_type'] == 'NORMAL']
        
        # Create all combinations of tumor-normal pairs
        for _, tumor in tumors.iterrows():
            for _, normal in normals.iterrows():
                pair_id = f"{tumor['sample_id']}_vs_{normal['sample_id']}"
                pairs.append({
                    'pair_id': pair_id,
                    'patient': patient,
                    'tumor_id': tumor['sample_id'],
                    'normal_id': normal['sample_id'],
                    'tumor_bam': tumor['bam_path'],
                    'normal_bam': normal['bam_path'],
                    'basecall_model': tumor['basecall_model'],
                    'clairs_model': tumor['clairs_model']
                })
    return pd.DataFrame(pairs)

# Create tumor-normal pairs dataframe
tn_pairs_df = create_tn_pairs(samples_df)
PAIRS = list(tn_pairs_df['pair_id'])

# Get unique sample IDs for individual sample processing
ALL_SAMPLES = list(samples_df['sample_id'].unique())
TUMOR_SAMPLES = list(samples_df[samples_df['sample_type'] == 'TUMOR']['sample_id'])
NORMAL_SAMPLES = list(samples_df[samples_df['sample_type'] == 'NORMAL']['sample_id'])

# Create lookup dictionaries
sample_info = samples_df.set_index('sample_id').to_dict('index')
pair_info = tn_pairs_df.set_index('pair_id').to_dict('index')

# Methylation parameters
minCov = config["methylation"]["modkit"]["min_coverage"]
haplotypes = config["methylation"]["haplotypes"]
MODCODES = config["methylation"]["mod_codes"]
SAMPLING_FRAC = config["methylation"]["modkit"]["sampling_fraction"]
FILTER_THRESHOLD = config["methylation"]["modkit"]["filter_threshold"]

# Helper functions to get sample paths
def get_sample_bam(sample_id):
    return sample_info[sample_id]['bam_path']

def get_sample_type(sample_id):
    return sample_info[sample_id]['sample_type'].lower()

# ─── FINAL OUTPUTS ─────────────────────────────────────────────────────────────
rule all:
    input:
        # Clair3 outputs for tumor-normal pairs
        expand(f"{OUTDIR}/run_clair/{{pair}}/{{pair}}.snv.vcf.gz", pair=PAIRS),
        expand(f"{OUTDIR}/run_clair/{{pair}}/clair3_normal_germline_output.vcf.gz", pair=PAIRS),
        expand(f"{OUTDIR}/run_clair/{{pair}}/done.{{pair}}.txt", pair=PAIRS),

        # Filtered PASS VCFs
        expand(f"{OUTDIR}/filter_pass/tumor/{{pair}}/{{tumor_id}}.pass.vcf.gz", 
               zip, pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS]),
        expand(f"{OUTDIR}/filter_pass/normal/{{pair}}/{{normal_id}}.pass.vcf.gz", 
               zip, pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS]),

        # VCF validation (ensures variants exist before phasing)
        expand(f"{OUTDIR}/vcf_validation/tumor/{{pair}}/{{tumor_id}}.validated", 
               zip, pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS]),
        expand(f"{OUTDIR}/vcf_validation/normal/{{pair}}/{{normal_id}}.validated", 
               zip, pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS]),

        # Individual sample processing (longphase modcall)
        expand(f"{OUTDIR}/longphase_modcall/{{sample}}/modcall_{{sample}}.vcf", sample=ALL_SAMPLES),
        expand(f"{OUTDIR}/longphase_modcall/{{sample}}/done.{{sample}}.txt", sample=ALL_SAMPLES),

        # Pair-specific co-phasing
        expand(f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/tumor/{{tumor_id}}_mod.vcf",
               zip, pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS]),
        expand(f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/normal/{{normal_id}}_mod.vcf",
               zip, pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS]),

        # Haplotagged BAMs
        expand(f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}_haplotagged.bam",
               zip, pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS]),
        expand(f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/normal/{{normal_id}}_haplotagged.bam",
               zip, pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS]),

        # SV calling outputs
        expand(f"{SV_OUTDIR}/{{pair}}/normal/{{normal_id}}.normal.sniffles.filtered.vcf.gz",
               zip, pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS]),
        expand(f"{SV_OUTDIR}/{{pair}}/tumor/{{tumor_id}}.tumor.sniffles.filtered.vcf.gz",
               zip, pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS]),
        expand(f"{SV_OUTDIR}/{{pair}}/somatic/{{pair}}.somatic.filtered.sv.vcf.gz", pair=PAIRS),

        # Methylation outputs
        expand(f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}.bed",
               pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS], haplotype=haplotypes),
        expand(f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}.bed",
               pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS], haplotype=haplotypes),

        # DMR outputs
        expand(f"{OUTDIR}/modkit_dmr_pairCpGs/{{pair}}/HP{{haplotype}}/{{pair}}_HP{{haplotype}}_TumorNormal_dmr_diff.bed.gz",
               pair=PAIRS, haplotype=haplotypes),

# ─── 2) RUN Clair3 (tumor+normal pair) ────────────────────────────────────────
rule run_clair:
    input:
        tumor  = lambda wc: pair_info[wc.pair]['tumor_bam'],
        normal = lambda wc: pair_info[wc.pair]['normal_bam'],
    output:
        tumor_vcf       = f"{OUTDIR}/run_clair/{{pair}}/{{pair}}.snv.vcf.gz",
        normal_germ_vcf = f"{OUTDIR}/run_clair/{{pair}}/clair3_normal_germline_output.vcf.gz",
        doneflag        = f"{OUTDIR}/run_clair/{{pair}}/done.{{pair}}.txt",
    threads: config["resources"]["threads"]["clair3"]
    resources:
        mem_mb = config["resources"]["memory_mb"]["clair3"]
    singularity: CLAIRS_IMG
    params:
        ref           = REF,
        platform      = lambda wc: pair_info[wc.pair]['clairs_model'],
        ctg_names     = CONTIG if CONTIG else "",
        conda_prefix  = CONDAPREFIX,
        qual_threshold = config["variant_calling"]["clair3"]["quality_threshold"],
        enable_germline = "--enable_clair3_germline_output" if config["variant_calling"]["clair3"]["enable_germline_output"] else "",
    log:
        f"{config['logging']['log_dir']}/run_clair/{{pair}}/{{pair}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/run_clair/{wildcards.pair}
        
        # Define platforms that support indel calling
        INDEL_PLATFORMS="ont_r10_dorado_sup_5khz_ssrs ont_r10_dorado_sup_5khz_ss ont_r10_guppy ont_r10 ont hifi_revio hifi_sequel2 ont_r10_dorado_sup_5khz ont_r10_dorado_hac_5khz ont_r10_dorado_sup_4khz ont_r10_dorado_hac_5khz_liquid hifi_revio_ssrs hifi_revio_ss"
        
        # Check if platform supports indel calling
        ENABLE_INDEL=""
        if [[ " $INDEL_PLATFORMS " =~ " {params.platform} " ]]; then
            ENABLE_INDEL="--enable_indel_calling"
            echo "Platform {params.platform} supports indel calling" | tee -a {log}
        else
            echo "WARNING: Platform {params.platform} does not support indel calling" | tee -a {log}
        fi
        
        # Run Clair3
        /opt/bin/run_clairs \
          --tumor_bam_fn   {input.tumor} \
          --normal_bam_fn  {input.normal} \
          --output_prefix {wildcards.pair}.snv \
          --ref_fn         {params.ref} \
          --threads        {threads} \
          --qual           {params.qual_threshold} \
          --platform       '{params.platform}' \
          $ENABLE_INDEL \
          {params.enable_germline} \
          --output_dir     {OUTDIR}/run_clair/{wildcards.pair} \
          --conda_prefix   {params.conda_prefix} \
          --include_all_ctgs \
          $(if [ -n "{params.ctg_names}" ]; then echo "--ctg_name {params.ctg_names}"; fi) \
        2>&1 | tee -a {log}
        
        # Check if the command succeeded
        if [ ${{PIPESTATUS[0]}} -eq 0 ]; then
            # Check if output files exist
            if [ -f "{output.tumor_vcf}" ] && [ -f "{output.normal_germ_vcf}" ]; then
                touch {output.doneflag}
                echo "Clair3 completed successfully" >> {log}
            else
                # Handle case where Clair3 succeeded but produced no variants
                echo "Clair3 completed but no output files found. Creating empty VCFs..." >> {log}
                
                # Create temporary uncompressed VCF files
                TEMP_TUMOR="{OUTDIR}/run_clair/{wildcards.pair}/{wildcards.pair}.snv.vcf"
                TEMP_NORMAL="{OUTDIR}/run_clair/{wildcards.pair}/clair3_normal_germline_output.vcf"
                
                # Create empty VCF files with proper headers
                echo '##fileformat=VCFv4.2' > "$TEMP_TUMOR"
                echo '##source=ClairS' >> "$TEMP_TUMOR"
                echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL' >> "$TEMP_TUMOR"
                
                # Copy for normal VCF
                cp "$TEMP_TUMOR" "$TEMP_NORMAL"
                
                # Compress and index
                bgzip -f "$TEMP_TUMOR"
                bgzip -f "$TEMP_NORMAL"
                tabix -p vcf "{output.tumor_vcf}"
                tabix -p vcf "{output.normal_germ_vcf}"
                
                touch {output.doneflag}
                echo "Created empty VCF files" >> {log}
            fi
        else
            echo "Clair3 failed with exit code ${{PIPESTATUS[0]}}" >> {log}
            exit 1
        fi
        """

# ─── 3) FILTER PASS for tumor & normal ─────────────────────────────────────────
rule filter_pass_tumor:
    input:
        raw_vcf = f"{OUTDIR}/run_clair/{{pair}}/{{pair}}.snv.vcf.gz"
    output:
        pass_vcf = f"{OUTDIR}/filter_pass/tumor/{{pair}}/{{tumor_id}}.pass.vcf.gz",
        pass_tbi = f"{OUTDIR}/filter_pass/tumor/{{pair}}/{{tumor_id}}.pass.vcf.gz.tbi"
    singularity: MODKIT_IMG
    shell:
        r"""
        mkdir -p $(dirname {output.pass_vcf})
        bcftools view -f PASS {input.raw_vcf} -Oz -o {output.pass_vcf}
        tabix -p vcf {output.pass_vcf}
        """

rule filter_pass_normal:
    input:
        raw_vcf = f"{OUTDIR}/run_clair/{{pair}}/clair3_normal_germline_output.vcf.gz"
    output:
        pass_vcf = f"{OUTDIR}/filter_pass/normal/{{pair}}/{{normal_id}}.pass.vcf.gz",
        pass_tbi = f"{OUTDIR}/filter_pass/normal/{{pair}}/{{normal_id}}.pass.vcf.gz.tbi"
    singularity: MODKIT_IMG
    shell:
        r"""
        mkdir -p $(dirname {output.pass_vcf})
        bcftools view -f PASS {input.raw_vcf} -Oz -o {output.pass_vcf}
        tabix -p vcf {output.pass_vcf}
        """

# ─── 3b) VALIDATE VCFs have variants ─────────────────────────────────
rule validate_vcf_has_variants:
    input:
        vcf = lambda wc: f"{OUTDIR}/filter_pass/{wc.type}/{wc.pair}/{wc.sample_id}.pass.vcf.gz"
    output:
        validation = f"{OUTDIR}/vcf_validation/{{type}}/{{pair}}/{{sample_id}}.validated"
    log:
        f"{config['logging']['log_dir']}/vcf_validation/{{type}}/{{pair}}/{{sample_id}}.log"
    shell:
        r"""
        # Count variants (excluding header lines)
        VARIANT_COUNT=$(zcat {input.vcf} | grep -v "^#" | wc -l)
        
        echo "Validating {input.vcf}" | tee {log}
        echo "Found $VARIANT_COUNT variants" | tee -a {log}
        
        if [ "$VARIANT_COUNT" -eq 0 ]; then
            echo "ERROR: No PASS variants found in {input.vcf}" | tee -a {log}
            echo "This sample cannot be phased - no heterozygous SNPs available" | tee -a {log}
            echo "Consider using full BAM files instead of subsampled" | tee -a {log}
            exit 1
        else
            echo "PASS: Found $VARIANT_COUNT variants for phasing" | tee -a {log}
            echo "$VARIANT_COUNT" > {output.validation}
        fi
        """

# ─── 4) LONGPHASE MODCALL (per individual sample) ──────────────────────────────────
rule longphase_modcall:
    input:
        bam = lambda wc: get_sample_bam(wc.sample)
    output:
        out_modcall = f"{OUTDIR}/longphase_modcall/{{sample}}/modcall_{{sample}}.vcf",
        done_modcall= f"{OUTDIR}/longphase_modcall/{{sample}}/done.{{sample}}.txt"
    params:
        reference_genome   = REF,
        threads            = THREADS,
        out_modcall_suffix= f"{OUTDIR}/longphase_modcall/{{sample}}/modcall_{{sample}}"
    threads: THREADS
    singularity: ONT_TOOLS_IMG
    log:
        f"logs/longphase_modcall/{{sample}}/{{sample}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/longphase_modcall/{wildcards.sample}
        longphase modcall \
          -b {input.bam} \
          -t {params.threads} \
          -o {params.out_modcall_suffix} \
          -r {params.reference_genome} \
        && touch {output.done_modcall} 2> {log}
        """

# ─── 5) LONGPHASE COPHASE for tumor in pair ──────────────────────────────────
rule mod_SNV_coPhase_tumor:
    input:
        bamfile     = lambda wc: pair_info[wc.pair]['tumor_bam'],
        modcallfile = lambda wc: f"{OUTDIR}/longphase_modcall/{pair_info[wc.pair]['tumor_id']}/modcall_{pair_info[wc.pair]['tumor_id']}.vcf",
        snpFile     = f"{OUTDIR}/filter_pass/tumor/{{pair}}/{{tumor_id}}.pass.vcf.gz",
        validation  = f"{OUTDIR}/vcf_validation/tumor/{{pair}}/{{tumor_id}}.validated"  # Require validation
    output:
        co_phased_mod_vcf  = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/tumor/{{tumor_id}}_mod.vcf",
        co_phased_vcf      = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/tumor/{{tumor_id}}.vcf",
        done_mod_SNV_coPhase= f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/tumor/done.{{tumor_id}}.txt"
    params:
        reference_genome    = REF,
        threads             = THREADS,
        out_coPhase_prefix  = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/tumor/{{tumor_id}}",
        ont_flag            = "--ont" if config["phasing"]["longphase"]["use_ont_model"] else ""
    threads: THREADS
    singularity: ONT_TOOLS_IMG
    log:
        f"{config['logging']['log_dir']}/mod_SNV_coPhase/{{pair}}/tumor/{{tumor_id}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/mod_SNV_coPhase/{wildcards.pair}/tumor
        
        longphase phase \
          -s {input.snpFile} \
          --mod-file {input.modcallfile} \
          -b {input.bamfile} \
          -r {params.reference_genome} \
          -t {params.threads} \
          -o {params.out_coPhase_prefix} \
          {params.ont_flag} \
        && touch {output.done_mod_SNV_coPhase} 2>&1 | tee -a {log}
        """

# ─── 5b) LONGPHASE COPHASE for normal in pair ──────────────────────────────────
rule mod_SNV_coPhase_normal:
    input:
        bamfile     = lambda wc: pair_info[wc.pair]['normal_bam'],
        modcallfile = lambda wc: f"{OUTDIR}/longphase_modcall/{pair_info[wc.pair]['normal_id']}/modcall_{pair_info[wc.pair]['normal_id']}.vcf",
        snpFile     = f"{OUTDIR}/filter_pass/normal/{{pair}}/{{normal_id}}.pass.vcf.gz",
        validation  = f"{OUTDIR}/vcf_validation/normal/{{pair}}/{{normal_id}}.validated"  # Require validation
    output:
        co_phased_mod_vcf  = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/normal/{{normal_id}}_mod.vcf",
        co_phased_vcf      = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/normal/{{normal_id}}.vcf",
        done_mod_SNV_coPhase= f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/normal/done.{{normal_id}}.txt"
    params:
        reference_genome    = REF,
        threads             = THREADS,
        out_coPhase_prefix  = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/normal/{{normal_id}}",
        ont_flag            = "--ont" if config["phasing"]["longphase"]["use_ont_model"] else ""
    threads: THREADS
    singularity: ONT_TOOLS_IMG
    log:
        f"{config['logging']['log_dir']}/mod_SNV_coPhase/{{pair}}/normal/{{normal_id}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/mod_SNV_coPhase/{wildcards.pair}/normal
        
        longphase phase \
          -s {input.snpFile} \
          --mod-file {input.modcallfile} \
          -b {input.bamfile} \
          -r {params.reference_genome} \
          -t {params.threads} \
          -o {params.out_coPhase_prefix} \
          {params.ont_flag} \
        && touch {output.done_mod_SNV_coPhase} 2>&1 | tee -a {log}
        """

# ─── 6) LONGPHASE HAPLOTAG for tumor ─────────────────────────────────
rule mod_SNV_coPhased_haplotagged_tumor:
    input:
        bamfile           = lambda wc: pair_info[wc.pair]['tumor_bam'],
        co_phased_mod_vcf = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/tumor/{{tumor_id}}_mod.vcf",
        co_phased_vcf     = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/tumor/{{tumor_id}}.vcf"
    params:
        reference_genome       = REF,
        threads                = THREADS,
        out_haplotagged_prefix = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}_haplotagged"
    output:
        haplotagged_bam = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}_haplotagged.bam",
        haplotagged_bai = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}_haplotagged.bam.bai",
        done = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/done.{{tumor_id}}.txt"
    threads: THREADS
    singularity: ONT_TOOLS_IMG
    log:
        f"logs/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/mod_SNV_coPhased_haplotagged/{wildcards.pair}/tumor
        longphase haplotag \
          -s {input.co_phased_vcf} \
          --mod-file {input.co_phased_mod_vcf} \
          -b {input.bamfile} \
          -r {params.reference_genome} \
          -t {params.threads} \
          -o {params.out_haplotagged_prefix} \
        && touch {output.done} 2> {log}
        samtools index {output.haplotagged_bam}
        """

# ─── 6b) LONGPHASE HAPLOTAG for normal ─────────────────────────────────
rule mod_SNV_coPhased_haplotagged_normal:
    input:
        bamfile           = lambda wc: pair_info[wc.pair]['normal_bam'],
        co_phased_mod_vcf = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/normal/{{normal_id}}_mod.vcf",
        co_phased_vcf     = f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/normal/{{normal_id}}.vcf"
    params:
        reference_genome       = REF,
        threads                = THREADS,
        out_haplotagged_prefix = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/normal/{{normal_id}}_haplotagged"
    output:
        haplotagged_bam = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/normal/{{normal_id}}_haplotagged.bam",
        haplotagged_bai = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/normal/{{normal_id}}_haplotagged.bam.bai",
        done = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/normal/done.{{normal_id}}.txt"
    threads: THREADS
    singularity: ONT_TOOLS_IMG
    log:
        f"logs/mod_SNV_coPhased_haplotagged/{{pair}}/normal/{{normal_id}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/mod_SNV_coPhased_haplotagged/{wildcards.pair}/normal
        longphase haplotag \
          -s {input.co_phased_vcf} \
          --mod-file {input.co_phased_mod_vcf} \
          -b {input.bamfile} \
          -r {params.reference_genome} \
          -t {params.threads} \
          -o {params.out_haplotagged_prefix} \
        && touch {output.done} 2> {log}
        samtools index {output.haplotagged_bam}
        """

# ─── 7) STRUCTURAL VARIANT CALLING ──────────────────────────
rule normal_sv:
    input:
        bam = lambda wc: pair_info[wc.pair]['normal_bam']
    output:
        raw_vcf      = f"{SV_OUTDIR}/{{pair}}/normal/{{normal_id}}.normal.sniffles.vcf.gz",
        snf          = f"{SV_OUTDIR}/{{pair}}/normal/{{normal_id}}.normal.sniffles.snf",
        filtered_vcf = f"{SV_OUTDIR}/{{pair}}/normal/{{normal_id}}.normal.sniffles.filtered.vcf.gz"
    threads: THREADS_SV
    singularity: SIF_SNIFFLES
    log:
        f"logs/sv/{{pair}}/normal/{{normal_id}}.log"
    params:
        tandem      = TANDEM_REPS,
        ref         = REF,
        contig      = CONTIG,
        filter_expr = NORMAL_FILTER_EXPR
    shell:
        r"""
        mkdir -p {SV_OUTDIR}/{wildcards.pair}/normal
        sniffles --input {input.bam} --vcf {output.raw_vcf} \
                 --output-rnames \
                 --tandem-repeats {params.tandem} --reference {params.ref} --snf {output.snf} \
                 2> {log}
        bcftools view -f "{params.filter_expr}" {output.raw_vcf} -Oz -o {output.filtered_vcf} 2>> {log}
        tabix -f -p vcf {output.filtered_vcf} 2>> {log}
        """

rule tumor_sv:
    input:
        bam = lambda wc: pair_info[wc.pair]['tumor_bam']
    output:
        raw_vcf      = f"{SV_OUTDIR}/{{pair}}/tumor/{{tumor_id}}.tumor.sniffles.vcf.gz",
        snf          = f"{SV_OUTDIR}/{{pair}}/tumor/{{tumor_id}}.tumor.sniffles.snf",
        filtered_vcf = f"{SV_OUTDIR}/{{pair}}/tumor/{{tumor_id}}.tumor.sniffles.filtered.vcf.gz"
    threads: THREADS_SV
    singularity: SIF_SNIFFLES
    log:
        f"logs/sv/{{pair}}/tumor/{{tumor_id}}.log"
    params:
        tandem      = TANDEM_REPS,
        ref         = REF,
        contig      = CONTIG,
        filter_expr = TUMOR_FILTER_EXPR
    shell:
        r"""
        mkdir -p {SV_OUTDIR}/{wildcards.pair}/tumor
        sniffles --input {input.bam} --vcf {output.raw_vcf} \
                 --output-rnames \
                 --tandem-repeats {params.tandem} --reference {params.ref} --snf {output.snf} \
                 2> {log}
        bcftools view -f "{params.filter_expr}" {output.raw_vcf} -Oz -o {output.filtered_vcf} 2>> {log}
        tabix -f -p vcf {output.filtered_vcf} 2>> {log}
        """

# ─── 7b) SOMATIC SV (isec tumor vs normal) ────────────────────────────
rule somatic_sv:
    input:
        tumor  = lambda wc: f"{SV_OUTDIR}/{wc.pair}/tumor/{pair_info[wc.pair]['tumor_id']}.tumor.sniffles.filtered.vcf.gz",
        normal = lambda wc: f"{SV_OUTDIR}/{wc.pair}/normal/{pair_info[wc.pair]['normal_id']}.normal.sniffles.filtered.vcf.gz"
    output:
        sv_filtered = f"{SV_OUTDIR}/{{pair}}/somatic/{{pair}}.somatic.filtered.sv.vcf.gz"
    threads: THREADS_SV
    singularity: SIF_SNIFFLES
    log:
        f"logs/sv/{{pair}}/somatic.log"
    shell:
        r"""
        mkdir -p {SV_OUTDIR}/{wildcards.pair}/somatic/isec
        bcftools isec -c none \
            --output-type z \
            -p {SV_OUTDIR}/{wildcards.pair}/somatic/isec \
            --threads {threads} \
            {input.tumor} {input.normal} \
            2> {log}
        mv {SV_OUTDIR}/{wildcards.pair}/somatic/isec/0000.vcf.gz {output.sv_filtered}
        tabix -f -p vcf {output.sv_filtered} 2>> {log}
        """

# ─── METHYLATION PILEUP for tumor ────────────────────────────
rule modkit_pileupCpGsBed_tumor:
    input:
        haplotagged_bam = lambda wc: f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{wc.pair}/tumor/{pair_info[wc.pair]['tumor_id']}_haplotagged.bam"
    output:
        done = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/tumor/done.{{tumor_id}}.txt",
        mp_phased1_Bed = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/tumor/{{tumor_id}}_1.bed",
        mp_phased2_Bed = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/tumor/{{tumor_id}}_2.bed",
        mp_phasedUngrouped_Bed = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/tumor/{{tumor_id}}_ungrouped.bed",
    singularity: MODKIT_IMG
    params:
        ref = REF,
        mincov = config["methylation"]["modkit"]["min_coverage"],
        filter_threshold = config["methylation"]["modkit"]["filter_threshold"],
        sampling_frac = config["methylation"]["modkit"]["sampling_fraction"],
        combine_strands = "--combine-strands" if config["methylation"]["modkit"]["combine_strands"] else "",
        cpg_only = "--cpg" if config["methylation"]["modkit"]["cpg_only"] else "",
        partition_tag = config["phasing"]["longphase"]["partition_tag"],
        outdir = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/tumor/"
    threads: THREADS
    log:
        f"{config['logging']['log_dir']}/modkit_pileup/{{pair}}/tumor/{{tumor_id}}.log"
    shell:
        r"""
        mkdir -p {params.outdir}
        modkit pileup \
        --ref {params.ref} \
        {params.cpg_only} {params.combine_strands} \
        --filter-threshold {params.filter_threshold} \
        --sampling-frac {params.sampling_frac} \
        --threads {threads} \
        --prefix {wildcards.tumor_id} \
        --partition-tag {params.partition_tag} \
        {input.haplotagged_bam} {params.outdir} && touch {output.done} 2>&1 | tee -a {log}
        """

# ─── METHYLATION PILEUP for normal ────────────────────────────
rule modkit_pileupCpGsBed_normal:
    input:
        haplotagged_bam = lambda wc: f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{wc.pair}/normal/{pair_info[wc.pair]['normal_id']}_haplotagged.bam"
    output:
        done = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/normal/done.{{normal_id}}.txt",
        mp_phased1_Bed = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/normal/{{normal_id}}_1.bed",
        mp_phased2_Bed = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/normal/{{normal_id}}_2.bed",
        mp_phasedUngrouped_Bed = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/normal/{{normal_id}}_ungrouped.bed",
    singularity: MODKIT_IMG
    params:
        ref = REF,
        mincov = config["methylation"]["modkit"]["min_coverage"],
        filter_threshold = config["methylation"]["modkit"]["filter_threshold"],
        sampling_frac = config["methylation"]["modkit"]["sampling_fraction"],
        combine_strands = "--combine-strands" if config["methylation"]["modkit"]["combine_strands"] else "",
        cpg_only = "--cpg" if config["methylation"]["modkit"]["cpg_only"] else "",
        partition_tag = config["phasing"]["longphase"]["partition_tag"],
        outdir = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/normal/"
    threads: THREADS
    log:
        f"{config['logging']['log_dir']}/modkit_pileup/{{pair}}/normal/{{normal_id}}.log"
    shell:
        r"""
        mkdir -p {params.outdir}
        modkit pileup \
        --ref {params.ref} \
        {params.cpg_only} {params.combine_strands} \
        --filter-threshold {params.filter_threshold} \
        --sampling-frac {params.sampling_frac} \
        --threads {threads} \
        --prefix {wildcards.normal_id} \
        --partition-tag {params.partition_tag} \
        {input.haplotagged_bam} {params.outdir} && touch {output.done} 2>&1 | tee -a {log}
        """

# ─── DMR ANALYSIS ────────────────────────────
rule modkit_dmr_pair:
    input:
        tumor_bed = lambda wc: f"{OUTDIR}/modkit_pileupCpGsBed/{wc.pair}/tumor/{pair_info[wc.pair]['tumor_id']}_{wc.haplotype}.bed",
        normal_bed = lambda wc: f"{OUTDIR}/modkit_pileupCpGsBed/{wc.pair}/normal/{pair_info[wc.pair]['normal_id']}_{wc.haplotype}.bed"
    output:
        dmr_diff = f"{OUTDIR}/modkit_dmr_pairCpGs/{{pair}}/HP{{haplotype}}/{{pair}}_HP{{haplotype}}_TumorNormal_dmr_diff.bed.gz"
    singularity: MODKIT_IMG
    params:
        ref = REF,
        mincov = minCov
    threads: THREADS
    log:
        f"{config['logging']['log_dir']}/dmr/{{pair}}/HP{{haplotype}}/dmr.log"
    shell:
        r"""
        mkdir -p $(dirname {output.dmr_diff})
        
        # Run modkit dmr
        modkit dmr pair \
            -a {input.tumor_bed} \
            -b {input.normal_bed} \
            --ref {params.ref} \
            --base C \
            --threads {threads} \
            --log-filepath {output.dmr_diff}.log \
            -o {output.dmr_diff}.tmp 2>&1 | tee -a {log}
        
        # Sort and compress
        sort -k1,1 -k2,2n {output.dmr_diff}.tmp | bgzip -c > {output.dmr_diff}
        tabix -p bed {output.dmr_diff}
        rm -f {output.dmr_diff}.tmp
        """