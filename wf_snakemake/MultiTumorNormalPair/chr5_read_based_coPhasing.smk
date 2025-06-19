import pandas as pd
import yaml
import os
import itertools

# snakemake --snakefile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/chr5_read_based_coPhasing.smk \
# --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/config/slurmMinimalist \
# --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" --keep-going --rerun-incomplete -np

##TODO: Request longphase modcall team to add trailing `>` to output VCFs
# +        sed '/^##FORMAT=<ID=DP,/ s/$/>/' {output.out_modcall} \
# +          | bgzip -c > {output.out_modcallgz}

# ─── 0) LOAD CONFIGURATION ────────────────────────────────────────────────────
configfile: "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/chr5_config.yaml"

# ─── LOAD CHROMOSOMES FROM REFERENCE ─────────────────────────────────────────
def get_chromosomes_from_fai(fai_path):
    """Read chromosome names from FASTA index file"""
    chromosomes = []
    try:
        with open(fai_path, 'r') as f:
            for line in f:
                chrom = line.strip().split('\t')[0]
                chromosomes.append(chrom)
    except FileNotFoundError:
        print(f"Warning: FAI file not found at {fai_path}")
        print("Using default chromosome list")
        # Fallback to default chromosomes if FAI not found
        chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    return chromosomes

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


# Get chromosomes from reference FAI file
FAI_FILE = REF + ".fai"
CHROMOSOMES = get_chromosomes_from_fai(FAI_FILE)
print(f"Found {len(CHROMOSOMES)} chromosomes/contigs in reference")

# Optional: Filter to only major chromosomes if desired
# This filters to only chr1-22, chrX, chrY, chrM (adjust pattern as needed)
MAJOR_CHROMOSOMES = CHROMOSOMES
# MAJOR_CHROMOSOMES = [c for c in CHROMOSOMES if c.startswith('chr') and not '_' in c]
# print(f"Using {len(MAJOR_CHROMOSOMES)} major chromosomes for parallel processing")

TYPE = ["tumor", "normal"]

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
        *[
            expand(
                f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{t}/{{sid}}.phase.{suf}.vcf.gz",
                zip,
                pair=PAIRS,
                sid=[ pair_info[p][f"{t}_id"] for p in PAIRS ]
            )
            for t in TYPE
            for suf in ("snp","sv","mod")
        ],
        # the done flags
        *[
            expand(
                f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{t}/done.{{sid}}.txt",
                zip,
                pair=PAIRS,
                sid=[ pair_info[p][f"{t}_id"] for p in PAIRS ]
            )
            for t in TYPE
        ],
        *[
            expand(
                f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{t}/{{sid}}_haplotagged.{ext}",
                zip,
                pair=PAIRS,
                sid=[ pair_info[p][f"{t}_id"] for p in PAIRS ],
            )
            for t in TYPE
            for ext in ("bam", "bam.bai", "subsampled.bam", "subsampled.bam.bai")
        ],
        # done flags
        *[
            expand(
                f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{t}/done.{{sid}}.txt",
                zip,
                pair=PAIRS,
                sid=[ pair_info[p][f"{t}_id"] for p in PAIRS ],
            )
            for t in TYPE
        ],
        # SNV+SV+mod co‐phased haplotagged BAMs
        *[
            expand(
                f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{t}/{{sid}}_haplotagged.{ext}",
                zip,
                pair=PAIRS,
                sid=[ pair_info[p][f"{t}_id"] for p in PAIRS ],
            )
            for t in TYPE
            for ext in ("bam", "bam.bai", "subsampled.bam", "subsampled.bam.bai")
        ],
        # done flags
        *[
            expand(
                f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{t}/done.{{sid}}.txt",
                zip,
                pair=PAIRS,
                sid=[ pair_info[p][f"{t}_id"] for p in PAIRS ],
            )
            for t in TYPE
        ],
        expand(
            f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/tumor/{{tumor_id}}.bed.gz",
            zip,
            pair=PAIRS,
            tumor_id=[ pair_info[p]['tumor_id'] for p in PAIRS ]
        ),
        expand(
            f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/tumor/{{tumor_id}}.bed.gz.tbi",
            zip,
            pair=PAIRS,
            tumor_id=[ pair_info[p]['tumor_id'] for p in PAIRS ]
        ),
        expand(
            f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/tumor/done.{{tumor_id}}.txt",
            zip,
            pair=PAIRS,
            tumor_id=[ pair_info[p]['tumor_id'] for p in PAIRS ]
        ),

        # modkit_pileup_unphased_TN outputs for normal
        expand(
            f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/normal/{{normal_id}}.bed.gz",
            zip,
            pair=PAIRS,
            normal_id=[ pair_info[p]['normal_id'] for p in PAIRS ]
        ),
        expand(
            f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/normal/{{normal_id}}.bed.gz.tbi",
            zip,
            pair=PAIRS,
            normal_id=[ pair_info[p]['normal_id'] for p in PAIRS ]
        ),
        expand(
            f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/normal/done.{{normal_id}}.txt",
            zip,
            pair=PAIRS,
            normal_id=[ pair_info[p]['normal_id'] for p in PAIRS ]
        ),
## Pair-specific mod-snv-sv co-phasing
        # expand(f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/tumor/{{tumor_id}}_mod.vcf",
        #        zip, pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS]),
        # expand(f"{OUTDIR}/mod_SNV_coPhase/{{pair}}/normal/{{normal_id}}_mod.vcf",
        #        zip, pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS]),

        # Haplotagged BAMs
        expand(f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}_haplotagged.bam",
               zip, pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS]),
        expand(f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/normal/{{normal_id}}_haplotagged.bam",
               zip, pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS]),
        # Haplotagged BAMs; subsampled (if applicable)
        expand(f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}_haplotagged.subsampled.bam",
               zip, pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS]),
        expand(f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/normal/{{normal_id}}_haplotagged.subsampled.bam",
               zip, pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS]),
#subsampledHaplotagged_bam = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}_haplotagged.subsampled.bam",
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


        expand(f"{OUTDIR}/highQual_modkit_pileup_sortTabixBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}_sorted.bed",
               pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS], haplotype=haplotypes),
        expand(f"{OUTDIR}/highQual_modkit_pileup_sortTabixBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}_sorted.bed.gz",
               pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS], haplotype=haplotypes),

            # methylation sort tabix outputs;
        expand(f"{OUTDIR}/modkit_pileup_sortTabixBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}_sorted.bed",
               pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS], haplotype=haplotypes),
        expand(f"{OUTDIR}/modkit_pileup_sortTabixBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}_sorted.bed.gz",
               pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS], haplotype=haplotypes),

            # High-quality methylation outputs
        expand(f"{OUTDIR}/highQual_modkit_pileup_sortTabixBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}_sorted.bed",
               pair=PAIRS, tumor_id=[pair_info[p]['tumor_id'] for p in PAIRS], haplotype=haplotypes),
        expand(f"{OUTDIR}/highQual_modkit_pileup_sortTabixBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}_sorted.bed.gz",
               pair=PAIRS, normal_id=[pair_info[p]['normal_id'] for p in PAIRS], haplotype=haplotypes),

        # DMR outputs
        expand(f"{OUTDIR}/modkit_dmr_pairCpGs/{{pair}}/HP{{haplotype}}/{{pair}}_HP{{haplotype}}_TumorNormal_dmr_diff.bed.gz",
               pair=PAIRS, haplotype=haplotypes),
        expand(f"{OUTDIR}/modkit_dmr_pairCpGs/{{pair}}/HP{{haplotype}}/{{pair}}_HP{{haplotype}}_TumorNormal_raw_segmentation.bed.gz",
               pair=PAIRS, haplotype=haplotypes),
#        seg_file = f"{OUTDIR}/modkit_dmr_pairCpGs/{{pair}}/HP{{haplotype}}/{{pair}}_HP{{haplotype}}_TumorNormal_raw_segmentation.bed.gz"
##DMR unphased
        *[
            expand(
                f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}{s}",
                pair=PAIRS
            )
            for s in [
                "_unphased_raw_dmr.bed",
                "_unphased_raw_segmentation.bed",
                "_unphased_TN_dmr.bed.gz",
                "_unphased_TN_dmr.bed.gz.tbi",
                "_unphased_TN_dmr_segmentation.bed.gz",
                "_unphased_TN_dmr_segmentation.bed.gz.tbi",
                "_unphased_TN_dmr_diff.bed",
                "_unphased_TN_dmr_diff.bed.gz",
                "_unphased_TN_dmr_diff.bed.gz.tbi",
            ]
        ]








# ─── 2) RUN Clair3 (tumor+normal pair) ────────────────────────────────────────
# ─── 2a) RUN Clair3 by chromosome ────────────────────────────────────────

# ─── 2a) RUN Clair3 by chromosome ────────────────────────────────────────
rule run_clair_by_chr:
    input:
        tumor  = lambda wc: pair_info[wc.pair]['tumor_bam'],
        normal = lambda wc: pair_info[wc.pair]['normal_bam'],
    output:
        tumor_vcf       = temp(f"{OUTDIR}/run_clair/{{pair}}/by_chr/{{chr}}/{{pair}}.snv.vcf.gz"),
        normal_germ_vcf = temp(f"{OUTDIR}/run_clair/{{pair}}/by_chr/{{chr}}/clair3_normal_germline_output.vcf.gz"),
        doneflag        = temp(f"{OUTDIR}/run_clair/{{pair}}/by_chr/{{chr}}/done.txt"),
    threads: config["resources"]["threads"]["clair3"]
    resources:
        mem_mb = config["resources"]["memory_mb"]["clair3"] // 2  # Less memory needed per chromosome
    singularity: CLAIRS_IMG
    params:
        ref           = REF,
        platform      = lambda wc: pair_info[wc.pair]['clairs_model'],
        conda_prefix  = CONDAPREFIX,
        qual_threshold = config["variant_calling"]["clair3"]["quality_threshold"],
        enable_germline = "--enable_clair3_germline_output" if config["variant_calling"]["clair3"]["enable_germline_output"] else "",
    log:
        f"{config['logging']['log_dir']}/run_clair/{{pair}}/by_chr/{{chr}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/run_clair/{wildcards.pair}/by_chr/{wildcards.chr}
        
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
        
        # Run Clair3 on specific chromosome
        echo "Processing chromosome {wildcards.chr} for pair {wildcards.pair}" | tee -a {log}
        
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
          --output_dir     {OUTDIR}/run_clair/{wildcards.pair}/by_chr/{wildcards.chr} \
          --conda_prefix   {params.conda_prefix} \
          --ctg_name       {wildcards.chr} \
        2>&1 | tee -a {log}
        
        # Check if the command succeeded
        if [ ${{PIPESTATUS[0]}} -eq 0 ]; then
            # Check if output files exist
            if [ -f "{output.tumor_vcf}" ] && [ -f "{output.normal_germ_vcf}" ]; then
                touch {output.doneflag}
                echo "Clair3 completed successfully for {wildcards.chr}" >> {log}
            else
                # Handle case where Clair3 succeeded but produced no variants
                echo "Clair3 completed but no output files found for {wildcards.chr}. Creating empty VCFs..." >> {log}
                
                # Create temporary uncompressed VCF files
                TEMP_TUMOR="{OUTDIR}/run_clair/{wildcards.pair}/by_chr/{wildcards.chr}/{wildcards.pair}.snv.vcf"
                TEMP_NORMAL="{OUTDIR}/run_clair/{wildcards.pair}/by_chr/{wildcards.chr}/clair3_normal_germline_output.vcf"
                
                # Create empty VCF files with proper headers
                echo '##fileformat=VCFv4.2' > "$TEMP_TUMOR"
                echo '##source=ClairS' >> "$TEMP_TUMOR"
                echo '##contig=<ID={wildcards.chr}>' >> "$TEMP_TUMOR"
                echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL' >> "$TEMP_TUMOR"
                
                # Copy for normal VCF
                cp "$TEMP_TUMOR" "$TEMP_NORMAL"
                
                # Compress and index
                bgzip -f "$TEMP_TUMOR"
                bgzip -f "$TEMP_NORMAL"
                tabix -p vcf "{output.tumor_vcf}"
                tabix -p vcf "{output.normal_germ_vcf}"
                
                touch {output.doneflag}
                echo "Created empty VCF files for {wildcards.chr}" >> {log}
            fi
        else
            echo "Clair3 failed for {wildcards.chr} with exit code ${{PIPESTATUS[0]}}" >> {log}
            exit 1
        fi
        """

# ─── 2b) MERGE Clair3 chromosome results ────────────────────────────────────
rule merge_clair_results:
    input:
        tumor_vcfs = lambda wc: expand(f"{OUTDIR}/run_clair/{wc.pair}/by_chr/{{chr}}/{wc.pair}.snv.vcf.gz", 
                                      chr=MAJOR_CHROMOSOMES),
        normal_vcfs = lambda wc: expand(f"{OUTDIR}/run_clair/{wc.pair}/by_chr/{{chr}}/clair3_normal_germline_output.vcf.gz", 
                                       chr=MAJOR_CHROMOSOMES),
        doneflags = lambda wc: expand(f"{OUTDIR}/run_clair/{wc.pair}/by_chr/{{chr}}/done.txt", 
                                     chr=MAJOR_CHROMOSOMES)
    output:
        tumor_vcf       = f"{OUTDIR}/run_clair/{{pair}}/{{pair}}.snv.vcf.gz",
        normal_germ_vcf = f"{OUTDIR}/run_clair/{{pair}}/clair3_normal_germline_output.vcf.gz",
        doneflag        = f"{OUTDIR}/run_clair/{{pair}}/done.{{pair}}.txt",
    threads: 4
    singularity: ONT_TOOLS_IMG
    log:
        f"{config['logging']['log_dir']}/run_clair/{{pair}}/merge.log"
    params:
        chromosomes = MAJOR_CHROMOSOMES  # Pass chromosome list as parameter
    shell:
        r"""
        mkdir -p {OUTDIR}/run_clair/{wildcards.pair}
        
        echo "Starting merge of Clair3 results for pair {wildcards.pair}" | tee {log}
        echo "Chromosomes to merge: {params.chromosomes}" | tee -a {log}
        
        # Build lists of VCFs maintaining chromosome order from reference
        TUMOR_VCFS=""
        NORMAL_VCFS=""
        MISSING_CHRS=""
        
        for chr in {params.chromosomes}; do
            tumor_vcf_path="{OUTDIR}/run_clair/{wildcards.pair}/by_chr/${{chr}}/{wildcards.pair}.snv.vcf.gz"
            normal_vcf_path="{OUTDIR}/run_clair/{wildcards.pair}/by_chr/${{chr}}/clair3_normal_germline_output.vcf.gz"
            
            if [ -f "$tumor_vcf_path" ] && [ -f "$normal_vcf_path" ]; then
                TUMOR_VCFS="$TUMOR_VCFS $tumor_vcf_path"
                NORMAL_VCFS="$NORMAL_VCFS $normal_vcf_path"
                echo "Found VCFs for chromosome $chr" | tee -a {log}
            else
                MISSING_CHRS="$MISSING_CHRS $chr"
                echo "Warning: Missing VCF files for chromosome $chr" | tee -a {log}
            fi
        done
        
        # Check if we have any VCFs to merge
        if [ -z "$TUMOR_VCFS" ]; then
            echo "ERROR: No VCF files found to merge!" | tee -a {log}
            exit 1
        fi
        
        if [ -n "$MISSING_CHRS" ]; then
            echo "WARNING: Missing results for chromosomes:$MISSING_CHRS" | tee -a {log}
        fi
        
        # Merge tumor VCFs
        echo "Merging tumor VCFs..." | tee -a {log}
        bcftools concat $TUMOR_VCFS \
            --allow-overlaps \
            --output-type z \
            --output {output.tumor_vcf} \
            --threads {threads} 2>&1 | tee -a {log}
        
        if [ ${{PIPESTATUS[0]}} -ne 0 ]; then
            echo "ERROR: Failed to merge tumor VCFs" | tee -a {log}
            exit 1
        fi
        
        # Index tumor VCF
        tabix -p vcf {output.tumor_vcf}
        
        # Merge normal germline VCFs
        echo "Merging normal germline VCFs..." | tee -a {log}
        bcftools concat $NORMAL_VCFS \
            --allow-overlaps \
            --output-type z \
            --output {output.normal_germ_vcf} \
            --threads {threads} 2>&1 | tee -a {log}
        
        if [ ${{PIPESTATUS[0]}} -ne 0 ]; then
            echo "ERROR: Failed to merge normal VCFs" | tee -a {log}
            exit 1
        fi
        
        # Index normal VCF
        tabix -p vcf {output.normal_germ_vcf}
        
        # Verify outputs and report statistics
        if [ -f "{output.tumor_vcf}" ] && [ -f "{output.normal_germ_vcf}" ]; then
            # Count variants in final VCFs
            TUMOR_COUNT=$(bcftools view -H {output.tumor_vcf} | wc -l)
            NORMAL_COUNT=$(bcftools view -H {output.normal_germ_vcf} | wc -l)
            
            echo "====== Merge Summary ======" | tee -a {log}
            echo "Successfully merged $(echo $TUMOR_VCFS | wc -w) chromosome VCFs" | tee -a {log}
            echo "Total tumor variants: $TUMOR_COUNT" | tee -a {log}
            echo "Total normal germline variants: $NORMAL_COUNT" | tee -a {log}
            echo "==========================" | tee -a {log}
            
            touch {output.doneflag}
        else
            echo "ERROR: Output files not created properly" | tee -a {log}
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
    singularity: ONT_TOOLS_IMG
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
    singularity: ONT_TOOLS_IMG
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
# ─── LONGPHASE MODCALL by chromosome ──────────────────────────────────
rule longphase_modcall_by_chr:
    input:
        bam = lambda wc: get_sample_bam(wc.sample)
    output:
        out_modcall = temp(f"{OUTDIR}/longphase_modcall/{{sample}}/by_chr/modcall_{{sample}}_{{chr}}.vcf"),
        out_modcallgz = f"{OUTDIR}/longphase_modcall/{{sample}}/by_chr/modcall_{{sample}}_{{chr}}.vcf.gz",
        out_modcallTabix = f"{OUTDIR}/longphase_modcall/{{sample}}/by_chr/modcall_{{sample}}_{{chr}}.vcf.gz.tbi"
    params:
        reference_genome = REF,
        threads = THREADS,
        out_modcall_suffix = f"{OUTDIR}/longphase_modcall/{{sample}}/by_chr/modcall_{{sample}}_{{chr}}"
    threads: THREADS
    singularity: ONT_TOOLS_IMG
    log:
        f"logs/longphase_modcall/{{sample}}/{{chr}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/longphase_modcall/{wildcards.sample}/by_chr
        longphase modcall \
          -b {input.bam} \
          -t {params.threads} \
          -o {params.out_modcall_suffix} \
          -r {params.reference_genome} \
          --region {wildcards.chr} \
        2> {log}

        sed '/^##FORMAT=<ID=DP,/ s/$/>/' {output.out_modcall} | bgzip -c > {output.out_modcallgz}
        #bgzip -c {output.out_modcall} > {output.out_modcallgz}
        tabix -p vcf {output.out_modcallgz}
        """

rule merge_longphase_modcall:
    input:
        vcfs = lambda wc: expand(f"{OUTDIR}/longphase_modcall/{wc.sample}/by_chr/modcall_{wc.sample}_{{chr}}.vcf.gz", 
                                chr=MAJOR_CHROMOSOMES),
        vcfsIndx = lambda wc: expand(f"{OUTDIR}/longphase_modcall/{wc.sample}/by_chr/modcall_{wc.sample}_{{chr}}.vcf.gz.tbi", 
                                chr=MAJOR_CHROMOSOMES)
    output:
        out_modcall = f"{OUTDIR}/longphase_modcall/{{sample}}/modcall_{{sample}}.vcf",
        done_modcall = f"{OUTDIR}/longphase_modcall/{{sample}}/done.{{sample}}.txt"
    threads: 2
    singularity: ONT_TOOLS_IMG
    log:
        f"logs/longphase_modcall/{{sample}}/merge.log"
    shell:
        r"""
        # Concatenate all chromosome VCFs
        bcftools concat {input.vcfs} \
            --allow-overlaps \
            -o {output.out_modcall} 2>&1 | tee {log}
        
        # Verify output
        if [ -f "{output.out_modcall}" ]; then
            touch {output.done_modcall}
            echo "Successfully merged $(echo {input.vcfs} | wc -w) chromosome VCFs" | tee -a {log}
        else
            echo "ERROR: Merge failed" | tee -a {log}
            exit 1
        fi
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
        subsampledHaplotagged_bam = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}_haplotagged.subsampled.bam",
        subsampledHaplotagged_bai = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/tumor/{{tumor_id}}_haplotagged.subsampled.bam.bai",
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


        #DELETE AFTER TESTING; CHANGE `SN:5` values 
        samtools view -h {output.haplotagged_bam} \
        | awk 'BEGIN{{OFS="\t"}}
         /^@SQ/ {{ if ($2=="SN:5") print; next }}
         /^@/  {{ print;      next }}
                {{ print }}
        ' \
        | samtools view -b -o {output.subsampledHaplotagged_bam} -
        
        samtools index {output.subsampledHaplotagged_bam}
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
        subsampledHaplotagged_bam = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/normal/{{normal_id}}_haplotagged.subsampled.bam",
        subsampledHaplotagged_bai = f"{OUTDIR}/mod_SNV_coPhased_haplotagged/{{pair}}/normal/{{normal_id}}_haplotagged.subsampled.bam.bai",
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

        #DELETE AFTER TESTING; CHANGE `SN:5` values 
        samtools view -h {output.haplotagged_bam} \
        | awk 'BEGIN{{OFS="\t"}}
         /^@SQ/ {{ if ($2=="SN:5") print; next }}
         /^@/  {{ print;      next }}
                {{ print }}
        ' \
        | samtools view -b -o {output.subsampledHaplotagged_bam} -
        
        samtools index {output.subsampledHaplotagged_bam}
        """

# ─── 8) SNV-SV-MOD CO-PHASING ─────────────────────────────────────────────
rule longphase_SNV_SV_MOD_co_phase:
    """
    Co-phase SNPs, SVs and 5mC modifications in one LongPhase invocation.
    Produces three phased VCFs (snp, sv, mod) plus a done file.
    """
    input:
        bam     = lambda wc: pair_info[wc.pair][f"{wc.type}_bam"],
        snp_vcf = f"{OUTDIR}/filter_pass/{{type}}/{{pair}}/{{sample_id}}.pass.vcf.gz",
        sv_vcf  = f"{SV_OUTDIR}/{{pair}}/{{type}}/{{sample_id}}.{{type}}.sniffles.filtered.vcf.gz",
        mod_vcf = f"{OUTDIR}/longphase_modcall/{{sample_id}}/modcall_{{sample_id}}.vcf"
    output:
        tmp_phased_snp = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}.vcf",
        tmp_phased_sv  = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}_SV.vcf",
        tmp_phased_mod = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}_mod.vcf",
        # Temporary files to handle a bug in LongPhase
        #tmp2_phased_snp = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}.phase.snp.vcf",
        #tmp2_phased_sv  = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}.phase.sv.vcf",
        #tmp2_phased_mod = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}.phase.mod.vcf",
        # Final phased VCFs
        phased_snp = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}.phase.snp.vcf.gz",
        phased_sv  = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}.phase.sv.vcf.gz",
        phased_mod = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}.phase.mod.vcf.gz",
        done       = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/done.{{sample_id}}.txt"
    params:
        ref     = REF,
        threads = THREADS
    singularity: ONT_TOOLS_IMG
    threads: THREADS
    log:
        f"{config['logging']['log_dir']}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample_id}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/SNV_SV_MOD_co_phase/{wildcards.pair}/{wildcards.type}

        longphase phase \
          -s {input.snp_vcf} \
          --sv-file {input.sv_vcf} \
          --mod-file {input.mod_vcf} \
          -b {input.bam} \
          -r {params.ref} \
          -t {params.threads} \
          -o {OUTDIR}/SNV_SV_MOD_co_phase/{wildcards.pair}/{wildcards.type}/{wildcards.sample_id} \
          --ont \
        2>&1 | tee {log}

        # compress & index each output VCF
        bgzip -f {output.tmp_phased_snp} -o {output.phased_snp}
        bgzip -f {output.tmp_phased_sv} -o {output.phased_sv}
        bgzip -f {output.tmp_phased_mod} -o {output.phased_mod}

        tabix -p vcf {output.phased_snp}
        tabix -p vcf {output.phased_sv}
        tabix -p vcf {output.phased_mod}

        touch {output.done}
        """

        #UNTIL LONGPHASE FIXES THIS BUG, we need to rename the output files
        #mv {output.tmp_phased_snp} {output.tmp2_phased_snp}
        #mv {output.tmp_phased_sv} {output.tmp2_phased_sv}
        #mv {output.tmp_phased_mod} {output.tmp2_phased_mod} 

rule SNV_SV_MOD_haplotag:
    """
    Haplotag a co-phased SNV+SV+mod VCF trio in one go, for either tumor or normal.
    """
    input:
        bam         = lambda wc: pair_info[wc.pair][f"{wc.type}_bam"],
        phased_snp  = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample}}.phase.snp.vcf.gz",
        phased_sv   = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample}}.phase.sv.vcf.gz",
        phased_mod  = f"{OUTDIR}/SNV_SV_MOD_co_phase/{{pair}}/{{type}}/{{sample}}.phase.mod.vcf.gz"
    output:
        hap_bam     = f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{{type}}/{{sample}}_haplotagged.bam",
        hap_bai     = f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{{type}}/{{sample}}_haplotagged.bam.bai",
        sub_bam     = f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{{type}}/{{sample}}_haplotagged.subsampled.bam",
        sub_bai     = f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{{type}}/{{sample}}_haplotagged.subsampled.bam.bai",
        done        = f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{{type}}/done.{{sample}}.txt"
    params:
        out_prefix  = f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{{pair}}/{{type}}/{{sample}}_haplotagged"
    threads: THREADS
    singularity: ONT_TOOLS_IMG
    log:
        f"{config['logging']['log_dir']}/haplotag/{{pair}}/{{type}}/{{sample}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/{wildcards.pair}/{wildcards.type}

        longphase haplotag \
          -s {input.phased_snp} \
          --sv-file {input.phased_sv} \
          --mod-file {input.phased_mod} \
          -b {input.bam} \
          -r {REF} \
          -t {threads} \
          -o {params.out_prefix} \
        2>&1 | tee {log}

        # index the full haplotagged BAM
        samtools index {output.hap_bam}

        # subsample out only 5 (or whatever filtering you want)
        samtools view -h {output.hap_bam} \
        | awk 'BEGIN{{OFS="\t"}}
           /^@SQ/ {{ if ($2=="SN:5") print; next }}
           /^@/  {{ print;      next }}
                  {{ print }}
        ' \
        | samtools view -b -o {output.sub_bam}

        samtools index {output.sub_bam}

        touch {output.done}
        """

# ─── OPTIONAL: WhatsHap SNV phasing ──────────────────────────────────────
rule whatshap_phase:
    """
    Phase SNPs (and small indels) with WhatsHap.
    Produces a phased VCF and its index.
    """
    input:
        vcf = f"{OUTDIR}/filter_pass/{{type}}/{{pair}}/{{sample_id}}.pass.vcf.gz",
        bam = lambda wc: pair_info[wc.pair][f"{wc.type}_bam"]
    output:
        phased = f"{OUTDIR}/whatsHap_phase/{{pair}}/{{type}}/{{sample_id}}.whatsHap.phased.vcf.gz",
        tbi    = f"{OUTDIR}/whatsHap_phase/{{pair}}/{{type}}/{{sample_id}}.whatsHap.phased.vcf.gz.tbi"
    params:
        ref     = REF,
        threads = THREADS
    singularity: ONT_TOOLS_IMG
    threads: THREADS
    log:
        f"{config['logging']['log_dir']}/whatsHap/{{pair}}/{{type}}/{{sample_id}}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.phased})
        whatshap phase \
            --reference {params.ref} \
            --output {output.phased} \
            --output-read-phases {output.phased}.reads \
            {input.vcf} {input.bam} \
            --threads {threads} 2>&1 | tee {log}
        tabix -p vcf {output.phased}
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
    singularity: ONT_TOOLS_IMG
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
    singularity: ONT_TOOLS_IMG
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
    singularity: ONT_TOOLS_IMG
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

rule modkit_pileup_unphased_TN:
    """
    Run an unphased modkit pileup on tumor or normal haplotagged BAM,
    then sort, bgzip & tabix the output bed.
    """
    input:
        # input callable is allowed
        haplotagged_bam=lambda wc: (
            f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/"
            f"{wc.pair}/{wc.type}/"
            f"{wc.sample_id}_haplotagged.subsampled.bam"
        )
    output:
        # these must be static patterns with wildcards
        tmp_bed   = f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/{{type}}/{{sample_id}}.bed",
        bed   = f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/{{type}}/{{sample_id}}.bed.gz",
        tbi   = f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/{{type}}/{{sample_id}}.bed.gz.tbi",
        done  = f"{OUTDIR}/modkit_pileup_unphased_TN/{{pair}}/{{type}}/done.{{sample_id}}.txt"
    params:
        ref             = REF,
        filter_threshold= config["methylation"]["modkit"]["filter_threshold"],
        sampling_frac   = config["methylation"]["modkit"]["sampling_fraction"],
        combine_strands = "--combine-strands" if config["methylation"]["modkit"]["combine_strands"] else "",
        cpg_only        = "--cpg"             if config["methylation"]["modkit"]["cpg_only"] else "",
        # outdir can still be a lambda
        outdir          = lambda wc: f"{OUTDIR}/modkit_pileup_unphased_TN/{wc.pair}/{wc.type}/"
    threads: THREADS
    log:
        # log must also be a static pattern
        f"{config['logging']['log_dir']}/modkit_pileup_unphased_TN/{{pair}}/{{type}}/{{sample_id}}.log"
    singularity: MODKIT_IMG
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}

        # 1) run modkit pileup
        modkit pileup \
            --ref {params.ref} \
            {params.cpg_only} {params.combine_strands} \
            --filter-threshold {params.filter_threshold} \
            --sampling-frac {params.sampling_frac} \
            --threads {threads} \
            {input.haplotagged_bam} {output.tmp_bed} \
        2>&1 | tee {log}

        # 2) sort + compress
        sort -k1,1 -k2,2n \
            {output.tmp_bed} \
        | bgzip -c > {output.bed}

        # 3) index
        tabix -p bed {output.bed}

        # 4) done
        touch {output.done}
        """

rule modkit_dmr_unphased_TN:
    """
    Call DMRs comparing unphased normal vs tumor 5mC BED outputs,
    emit raw DMRs + segmentation, then compress/index both,
    and finally filter “different” DMRs and compress/index those.
    """
    input:
        normal_bed = lambda wc: f"{OUTDIR}/modkit_pileup_unphased_TN/{wc.pair}/normal/{pair_info[wc.pair]['normal_id']}.bed.gz",
        tumor_bed  = lambda wc: f"{OUTDIR}/modkit_pileup_unphased_TN/{wc.pair}/tumor/{pair_info[wc.pair]['tumor_id']}.bed.gz"
    output:
        # raw outputs
        raw_dmr      = f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}_unphased_raw_dmr.bed",
        raw_seg      = f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}_unphased_raw_segmentation.bed",
        # compressed/indexed raw
        dmr_gz       = f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}_unphased_TN_dmr.bed.gz",
        dmr_tbi      = f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}_unphased_TN_dmr.bed.gz.tbi",
        seg_gz       = f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}_unphased_TN_dmr_segmentation.bed.gz",
        seg_tbi      = f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}_unphased_TN_dmr_segmentation.bed.gz.tbi",
        # filtered “different” DMRs
        diff_bed     = f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}_unphased_TN_dmr_diff.bed",
        diff_gz      = f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}_unphased_TN_dmr_diff.bed.gz",
        diff_tbi     = f"{OUTDIR}/modkit_dmr_unphased_TN/{{pair}}/{{pair}}_unphased_TN_dmr_diff.bed.gz.tbi"
    params:
        ref            = REF,
        base           = "C",
        mincov         = minCov,
        min_dmr_sites  = 3
    threads: THREADS
    singularity: MODKIT_IMG
    log:
        f"{config['logging']['log_dir']}/dmr_unphased_TN/{{pair}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.raw_dmr})

        # 1) call raw DMR + segmentation
        modkit dmr pair \
          -a {input.normal_bed} \
          -b {input.tumor_bed} \
          --ref {params.ref} \
          --base {params.base} \
          --min-valid-coverage {params.mincov} \
          --threads {threads} \
          --segment {output.raw_seg} \
          -o {output.raw_dmr} \
        2>&1 | tee {log}

        # 2) compress & index raw DMR
        sort -k1,1 -k2,2n {output.raw_dmr} \
          | bgzip -c > {output.dmr_gz}
        tabix -p bed {output.dmr_gz}

        # 3) compress & index raw segmentation
        sort -k1,1 -k2,2n {output.raw_seg} \
          | bgzip -c > {output.seg_gz}
        tabix -p bed {output.seg_gz}

        # 4) filter “different” DMRs
        awk -F '\t' 'NR>1 && $4=="different" && $6>={params.min_dmr_sites} && ($14>=0.5||$14<=-0.5) && ($15*$16>0)' \
          {output.raw_seg} > {output.diff_bed}

        # 5) compress & index filtered DMRs
        sort -k1,1 -k2,2n {output.diff_bed} \
          | bgzip -c > {output.diff_gz}
        tabix -p bed {output.diff_gz}
        """

#--prefix {wildcards.sample_id} \
# ─── METHYLATION PILEUP for tumor ────────────────────────────
rule modkit_pileupCpGsBed:
    """
    Run modkit pileup for either tumor or normal haplotagged BAMs,
    writing three bed outputs plus a done flag.
    """
    input:
        haplotagged_bam = (
            f"{OUTDIR}/SNV_SV_MOD_coPhased_haplotagged/"
            f"{{pair}}/{{type}}/{{sample_id}}_haplotagged.subsampled.bam"
        )
    output:
        mp1   = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/{{type}}/{{sample_id}}_1.bed",
        mp2   = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/{{type}}/{{sample_id}}_2.bed",
        mpU   = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/{{type}}/{{sample_id}}_ungrouped.bed",
        done  = f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/{{type}}/done.{{sample_id}}.txt"
    params:
        ref             = REF,
        mincov          = config["methylation"]["modkit"]["min_coverage"],
        filter_threshold= config["methylation"]["modkit"]["filter_threshold"],
        sampling_frac   = config["methylation"]["modkit"]["sampling_fraction"],
        combine_strands = "--combine-strands" if config["methylation"]["modkit"]["combine_strands"] else "",
        cpg_only        = "--cpg"             if config["methylation"]["modkit"]["cpg_only"] else "",
        outdir          = lambda wc: f"{OUTDIR}/modkit_pileupCpGsBed/{wc.pair}/{wc.type}/",
        prefix = lambda wc: wc.sample_id
    threads: THREADS
    singularity: MODKIT_IMG
    log:
        f"{config['logging']['log_dir']}/modkit_pileupCpGsBed/{{pair}}/{{type}}/{{sample_id}}.log"
    shell:
        r"""
        mkdir -p {params.outdir}

        modkit pileup \
          --ref {params.ref} \
          {params.cpg_only} {params.combine_strands} \
          --filter-threshold {params.filter_threshold} \
          --sampling-frac {params.sampling_frac} \
          --threads {threads} \
          --partition-tag HP \
          --prefix {params.prefix} \
          {input.haplotagged_bam} \
          {params.outdir} \
        2>&1 | tee -a {log}

        touch {output.done}
        """
#          --prefix {wildcards.sample_id} \

rule modkit_pileup_sortTabixBedNormal:
    input:
        bed= f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}.bed"
    output:
        sorted_bed= f"{OUTDIR}/modkit_pileup_sortTabixBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}_sorted.bed",
        sorted_tabixed= f"{OUTDIR}/modkit_pileup_sortTabixBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}_sorted.bed.gz"
    singularity: MODKIT_IMG
    params:
        mincov= minCov
    threads: THREADS
    shell:
        """
        echo "Sorting and filtering {input.bed} for min coverage {params.mincov}"
        sort -k1,1 -k2,2n {input.bed} > {output.sorted_bed} 
        bgzip -c {output.sorted_bed} > {output.sorted_tabixed}
        tabix -p bed {output.sorted_tabixed}
        """
rule modkit_pileup_sortTabixBedTumor:
    input:
        bed= f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}.bed"
    output:
        sorted_bed= f"{OUTDIR}/modkit_pileup_sortTabixBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}_sorted.bed",
        sorted_tabixed= f"{OUTDIR}/modkit_pileup_sortTabixBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}_sorted.bed.gz"
    singularity: MODKIT_IMG
    params:
        mincov= minCov
    threads: THREADS
    shell:
        """
        echo "Sorting and filtering {input.bed} for min coverage {params.mincov}"
        sort -k1,1 -k2,2n {input.bed} > {output.sorted_bed} 
        bgzip -c {output.sorted_bed} > {output.sorted_tabixed}
        tabix -p bed {output.sorted_tabixed}
        """

rule highQual_modkit_pileup_sortTabixBedNormal:
    input:
        bed= f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}.bed"
    output:
        sorted_bed= f"{OUTDIR}/highQual_modkit_pileup_sortTabixBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}_sorted.bed",
        sorted_tabixed= f"{OUTDIR}/highQual_modkit_pileup_sortTabixBed/{{pair}}/normal/{{normal_id}}_{{haplotype}}_sorted.bed.gz"
    singularity: MODKIT_IMG
    params:
        mincov= minCov
    threads: THREADS
    shell:
        """
        echo "Sorting and filtering {input.bed} for min coverage {params.mincov}"
        sort -k1,1 -k2,2n {input.bed} | awk -v min_cov={params.mincov} '$5 >= min_cov' > {output.sorted_bed} 
        bgzip -c {output.sorted_bed} > {output.sorted_tabixed}
        tabix -p bed {output.sorted_tabixed}
        """

rule highQual_modkit_pileup_sortTabixBedTumor:
    input:
        bed= f"{OUTDIR}/modkit_pileupCpGsBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}.bed"
    output:
        sorted_bed= f"{OUTDIR}/highQual_modkit_pileup_sortTabixBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}_sorted.bed",
        sorted_tabixed= f"{OUTDIR}/highQual_modkit_pileup_sortTabixBed/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}_sorted.bed.gz"
    singularity: MODKIT_IMG
    params:
        mincov= minCov
    threads: THREADS
    shell:
        """
        echo "Sorting and filtering {input.bed} for min coverage {params.mincov}"
        sort -k1,1 -k2,2n {input.bed} | awk -v min_cov={params.mincov} '$5 >= min_cov' > {output.sorted_bed} 
        bgzip -c {output.sorted_bed} > {output.sorted_tabixed}
        tabix -p bed {output.sorted_tabixed}
        """

# rule highQual_modkit_pileup_5mC_5hmC_BigWigs:
#     input:
#         bed= f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_sortTabixBed/{{sample}}/{{sample}}_{{haplotypes}}_sorted.bed"
#     output:
#         bedmethyl= f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_BigWigs/{{sample}}/{{sample}}_{{haplotypes}}_{{mod}}_sorted.bw"
#     singularity: MODKIT_IMG
#     params:
#         ref= REF,
#         genomeSizes = GENOMESIZES
#     threads: THREADS
#     log:
#         "logs/highQual_modkit_pileup_5mC_5hmC_BigWigs/{sample}/{sample}_{haplotypes}_{mod}_sorted.log"
#     shell:
#         """
#         echo "Converting {input.bed} to bedmethyl format"
#         modkit bedmethyl tobigwig --mod-codes {wildcards.mod} \
#         --suppress-progress \
#         --nthreads {threads} \
#         --sizes {params.genomeSizes} \
#         --log-filepath {log} \
#           {input.bed} {output.bedmethyl}
#         """


# ─── DMR ANALYSIS ────────────────────────────
rule modkit_dmr_pair:
    input:
        tumor_bed = lambda wc: f"{OUTDIR}/modkit_pileup_sortTabixBed/{wc.pair}/tumor/{pair_info[wc.pair]['tumor_id']}_{wc.haplotype}_sorted.bed.gz",
        normal_bed = lambda wc: f"{OUTDIR}/modkit_pileup_sortTabixBed/{wc.pair}/normal/{pair_info[wc.pair]['normal_id']}_{wc.haplotype}_sorted.bed.gz"
    output:
        dmr_diff = f"{OUTDIR}/modkit_dmr_pairCpGs/{{pair}}/HP{{haplotype}}/{{pair}}_HP{{haplotype}}_TumorNormal_dmr_diff.bed.gz",
        seg_file = f"{OUTDIR}/modkit_dmr_pairCpGs/{{pair}}/HP{{haplotype}}/{{pair}}_HP{{haplotype}}_TumorNormal_raw_segmentation.bed.gz"

    singularity: MODKIT_IMG
    params:
        ref = REF,
        mincov = minCov,
        base           = "C"
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
            --base {params.base} \
            --min-valid-coverage {params.mincov} \
            --threads {threads} \
            --log-filepath {output.dmr_diff}.log \
            --segment {output.seg_file}.tmp \
            -o {output.dmr_diff}.tmp \
            2>&1 | tee -a {log}
        
        # Sort and compress
        sort -k1,1 -k2,2n {output.dmr_diff}.tmp | bgzip -c > {output.dmr_diff}
        tabix -p bed {output.dmr_diff}
        rm -f {output.dmr_diff}.tmp

        sort -k1,1 -k2,2n {output.seg_file}.tmp | bgzip -c > {output.seg_file}
        tabix -p bed {output.seg_file}
        rm -f {output.seg_file}.tmp
        """

        #         sorted_bed= f"{OUTDIR}/modkit_pileup_sortTabixBedTumor/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}_sorted.bed",
        # sorted_tabixed= f"{OUTDIR}/modkit_pileup_sortTabixBedTumor/{{pair}}/tumor/{{tumor_id}}_{{haplotype}}_sorted.bed.gz"