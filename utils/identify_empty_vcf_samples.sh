#!/bin/bash
# Identify samples with empty VCFs that should be excluded from phasing workflow

# sh ../utils/identify_empty_vcf_samples.sh results/
# filter_samples_csv.py
# python3 ../utils/filter_samples_csv.py /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/read_based_mod_coPhasing/wf_snakemake/MultiTumorNormalPair/GET_CLAIR_MODELS/models4SubSampled/results/summary/copy_complete_summary.csv samples_to_exclude.unique_pairs -o samples_clean.csv


echo "=============================================="
echo "Identifying Samples with No Variants"
echo "=============================================="

RESULTS_DIR="${1:-results}"
OUTPUT_FILE="${2:-samples_to_exclude.txt}"

if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: Results directory not found: $RESULTS_DIR"
    echo "Usage: $0 [results_directory] [output_file]"
    exit 1
fi

# Initialize lists
> "$OUTPUT_FILE"
> "${OUTPUT_FILE}.pairs"

EMPTY_SAMPLES=0
EMPTY_PAIRS=0
TOTAL_SAMPLES=0

echo "Scanning filtered VCFs..."
echo ""

# Check tumor samples
echo "=== TUMOR SAMPLES ==="
for vcf_dir in "$RESULTS_DIR"/filter_pass/tumor/*; do
    [ -d "$vcf_dir" ] || continue
    PAIR=$(basename "$vcf_dir")
    
    for vcf in "$vcf_dir"/*.pass.vcf.gz; do
        [ -f "$vcf" ] || continue
        
        TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
        SAMPLE=$(basename "$vcf" .pass.vcf.gz)
        VARIANT_COUNT=$(zcat "$vcf" 2>/dev/null | grep -v "^#" | wc -l)
        
        if [ "$VARIANT_COUNT" -eq 0 ]; then
            echo "❌ $SAMPLE (tumor) - 0 variants in pair $PAIR"
            echo "TUMOR:$SAMPLE:$PAIR" >> "$OUTPUT_FILE"
            echo "$PAIR" >> "${OUTPUT_FILE}.pairs"
            EMPTY_SAMPLES=$((EMPTY_SAMPLES + 1))
        else
            echo "✓ $SAMPLE (tumor) - $VARIANT_COUNT variants"
        fi
    done
done

# Check normal samples
echo ""
echo "=== NORMAL SAMPLES ==="
for vcf_dir in "$RESULTS_DIR"/filter_pass/normal/*; do
    [ -d "$vcf_dir" ] || continue
    PAIR=$(basename "$vcf_dir")
    
    for vcf in "$vcf_dir"/*.pass.vcf.gz; do
        [ -f "$vcf" ] || continue
        
        TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
        SAMPLE=$(basename "$vcf" .pass.vcf.gz)
        VARIANT_COUNT=$(zcat "$vcf" 2>/dev/null | grep -v "^#" | wc -l)
        
        if [ "$VARIANT_COUNT" -eq 0 ]; then
            echo "❌ $SAMPLE (normal) - 0 variants in pair $PAIR"
            echo "NORMAL:$SAMPLE:$PAIR" >> "$OUTPUT_FILE"
            echo "$PAIR" >> "${OUTPUT_FILE}.pairs"
            EMPTY_SAMPLES=$((EMPTY_SAMPLES + 1))
        else
            echo "✓ $SAMPLE (normal) - $VARIANT_COUNT variants"
        fi
    done
done

# Get unique pairs with issues
sort -u "${OUTPUT_FILE}.pairs" > "${OUTPUT_FILE}.unique_pairs"
EMPTY_PAIRS=$(wc -l < "${OUTPUT_FILE}.unique_pairs")

# Summary
echo ""
echo "=============================================="
echo "SUMMARY"
echo "=============================================="
echo "Total samples checked: $TOTAL_SAMPLES"
echo "Samples with NO variants: $EMPTY_SAMPLES"
echo "Pairs affected: $EMPTY_PAIRS"
echo ""
echo "Files created:"
echo "- $OUTPUT_FILE : List of samples with no variants"
echo "- ${OUTPUT_FILE}.unique_pairs : Unique pairs to exclude"

if [ "$EMPTY_SAMPLES" -gt 0 ]; then
    echo ""
    echo "⚠️  WARNING: Found $EMPTY_SAMPLES samples with no variants!"
    echo ""
    echo "These samples CANNOT be phased because they have no heterozygous SNPs."
    echo ""
    echo "OPTIONS:"
    echo "1. Exclude these pairs from the workflow"
    echo "2. Use full BAM files instead of subsampled"
    echo "3. Check if the correct reference genome was used"
    echo ""
    echo "To exclude these pairs, you can:"
    echo "- Remove them from your samples CSV file"
    echo "- Or add logic to skip them in the Snakefile"
    echo ""
    echo "Affected pairs:"
    cat "${OUTPUT_FILE}.unique_pairs" | head -10
    if [ "$EMPTY_PAIRS" -gt 10 ]; then
        echo "... and $((EMPTY_PAIRS - 10)) more"
    fi
else
    echo ""
    echo "✓ Good news: All samples have variants for phasing!"
fi

echo ""
echo "=============================================="