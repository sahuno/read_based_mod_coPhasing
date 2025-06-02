#!/usr/bin/env python3
"""
Filter samples CSV to remove pairs with no variants.
This helps create a clean dataset for the phasing workflow.
"""

import pandas as pd
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(
        description='Filter samples CSV to remove problematic pairs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # After identifying empty VCF samples:
  ./identify_empty_vcf_samples.sh results/
  
  # Then filter your samples:
  python3 filter_samples_csv.py samples.csv samples_to_exclude.unique_pairs -o samples_filtered.csv
  
  # Or manually specify pairs to exclude:
  python3 filter_samples_csv.py samples.csv --exclude-pairs PAIR1 PAIR2 PAIR3
        """
    )
    
    parser.add_argument('samples_csv', help='Original samples CSV file')
    parser.add_argument('exclude_file', nargs='?', 
                       help='File with pair IDs to exclude (one per line)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output filtered CSV file')
    parser.add_argument('--exclude-pairs', nargs='+',
                       help='Pair IDs to exclude (alternative to exclude_file)')
    parser.add_argument('--exclude-patients', nargs='+',
                       help='Patient IDs to exclude entirely')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be filtered without creating output')
    
    args = parser.parse_args()
    
    # Read samples
    try:
        df = pd.read_csv(args.samples_csv)
        original_count = len(df)
        print(f"Loaded {original_count} samples from {args.samples_csv}")
    except Exception as e:
        print(f"Error reading {args.samples_csv}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Determine pairs to exclude
    exclude_pairs = set()
    
    # From file
    if args.exclude_file and os.path.exists(args.exclude_file):
        with open(args.exclude_file, 'r') as f:
            file_pairs = [line.strip() for line in f if line.strip()]
            exclude_pairs.update(file_pairs)
            print(f"Loaded {len(file_pairs)} pairs to exclude from {args.exclude_file}")
    
    # From command line
    if args.exclude_pairs:
        exclude_pairs.update(args.exclude_pairs)
        print(f"Added {len(args.exclude_pairs)} pairs from command line")
    
    # Create pair IDs in dataframe for matching
    if exclude_pairs:
        # First, we need to identify which samples are involved in excluded pairs
        excluded_samples = set()
        
        for pair_id in exclude_pairs:
            # Extract sample IDs from pair ID (format: TUMOR_vs_NORMAL)
            if '_vs_' in pair_id:
                tumor_id, normal_id = pair_id.split('_vs_', 1)
                excluded_samples.add(tumor_id)
                excluded_samples.add(normal_id)
        
        print(f"\nExcluding samples involved in problematic pairs:")
        for sample in sorted(excluded_samples):
            print(f"  - {sample}")
        
        # Filter out these samples
        df_filtered = df[~df['sample_id'].isin(excluded_samples)]
    else:
        df_filtered = df.copy()
    
    # Exclude entire patients if requested
    if args.exclude_patients:
        print(f"\nExcluding entire patients: {', '.join(args.exclude_patients)}")
        df_filtered = df_filtered[~df_filtered['patient'].isin(args.exclude_patients)]
    
    # Summary
    filtered_count = len(df_filtered)
    removed_count = original_count - filtered_count
    
    print("\n" + "=" * 60)
    print("FILTERING SUMMARY")
    print("=" * 60)
    print(f"Original samples: {original_count}")
    print(f"Filtered samples: {filtered_count}")
    print(f"Removed samples: {removed_count}")
    
    # Check balance of tumor/normal per patient
    print("\nRemaining samples per patient:")
    for patient in sorted(df_filtered['patient'].unique()):
        patient_df = df_filtered[df_filtered['patient'] == patient]
        n_tumor = len(patient_df[patient_df['sample_type'] == 'TUMOR'])
        n_normal = len(patient_df[patient_df['sample_type'] == 'NORMAL'])
        print(f"  {patient}: {n_tumor} tumor, {n_normal} normal")
        
        if n_tumor == 0 or n_normal == 0:
            print(f"    ⚠️  WARNING: Patient {patient} has no {'tumor' if n_tumor == 0 else 'normal'} samples!")
    
    # Save or show results
    if args.dry_run:
        print("\n" + "=" * 60)
        print("DRY RUN - No output file created")
        print("=" * 60)
        print("\nFiltered dataset preview:")
        print(df_filtered.head())
    else:
        df_filtered.to_csv(args.output, index=False)
        print(f"\n✓ Filtered samples saved to: {args.output}")
        print(f"\nYou can now run the workflow with:")
        print(f"  snakemake --configfile config.yaml --config samples_file={args.output}")
    
    # Additional warnings
    patients_removed = set(df['patient'].unique()) - set(df_filtered['patient'].unique())
    if patients_removed:
        print(f"\n⚠️  WARNING: Entire patients removed: {', '.join(sorted(patients_removed))}")

if __name__ == "__main__":
    main()