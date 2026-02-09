#!/usr/bin/env python3
"""
Merge TRUST4 TCR/BCR repertoire reports across samples.
Calculate diversity metrics: Shannon entropy, Simpson's index, clonality.
"""

import argparse
import os
import glob
import math
import pandas as pd
from collections import Counter

def parse_args():
    parser = argparse.ArgumentParser(description='Merge TRUST4 TCR reports')
    parser.add_argument('--input-dir', default='.')
    parser.add_argument('--pattern', default='*_TRUST4_report.tsv')
    parser.add_argument('--output-summary', default='tcr_repertoire_summary.tsv')
    parser.add_argument('--output-diversity', default='tcr_diversity_metrics.tsv')
    parser.add_argument('--output-clonotypes', default='tcr_clonotype_tracking.tsv')
    return parser.parse_args()

def shannon_entropy(counts):
    """Calculate Shannon entropy from a list of counts."""
    total = sum(counts)
    if total == 0:
        return 0.0
    probs = [c / total for c in counts if c > 0]
    return -sum(p * math.log2(p) for p in probs)

def simpson_index(counts):
    """Calculate Simpson's diversity index."""
    total = sum(counts)
    if total <= 1:
        return 0.0
    return 1 - sum(c * (c - 1) for c in counts) / (total * (total - 1))

def clonality(counts):
    """Calculate clonality (1 - normalized Shannon entropy)."""
    n = len([c for c in counts if c > 0])
    if n <= 1:
        return 0.0
    max_entropy = math.log2(n)
    if max_entropy == 0:
        return 0.0
    return 1 - shannon_entropy(counts) / max_entropy

def main():
    args = parse_args()

    files = glob.glob(os.path.join(args.input_dir, args.pattern))
    if not files:
        print("No TRUST4 report files found")
        pd.DataFrame().to_csv(args.output_summary, sep='\t', index=False)
        return

    all_reports = []
    diversity_metrics = []
    all_clonotypes = []

    for f in files:
        sample_id = os.path.basename(f).replace('_TRUST4_report.tsv', '')
        try:
            df = pd.read_csv(f, sep='\t')
            df['sample_id'] = sample_id
            all_reports.append(df)

            # Calculate diversity for each chain
            for chain in ['TRA', 'TRB', 'TRG', 'TRD', 'IGH', 'IGK', 'IGL']:
                chain_data = df[df.iloc[:, 0].str.contains(chain, na=False)] if len(df) > 0 else pd.DataFrame()
                if len(chain_data) > 0 and 'count' in df.columns:
                    counts = chain_data['count'].tolist()
                    diversity_metrics.append({
                        'sample_id': sample_id,
                        'chain': chain,
                        'n_clonotypes': len(counts),
                        'total_reads': sum(counts),
                        'shannon_entropy': round(shannon_entropy(counts), 4),
                        'simpson_index': round(simpson_index(counts), 4),
                        'clonality': round(clonality(counts), 4),
                    })

            # Top clonotypes
            if 'CDR3aa' in df.columns:
                top = df.nlargest(20, 'count') if 'count' in df.columns else df.head(20)
                top['sample_id'] = sample_id
                all_clonotypes.append(top)

        except Exception as e:
            print(f"Warning: Error processing {f}: {e}")

    # Summary
    summary_df = pd.DataFrame(diversity_metrics)
    summary_df.to_csv(args.output_summary, sep='\t', index=False)

    # Diversity
    diversity_df = pd.DataFrame(diversity_metrics)
    diversity_df.to_csv(args.output_diversity, sep='\t', index=False)

    # Clonotypes
    if all_clonotypes:
        clonotype_df = pd.concat(all_clonotypes, ignore_index=True)
        clonotype_df.to_csv(args.output_clonotypes, sep='\t', index=False)

    print(f"Processed {len(files)} samples")
    print(f"Total diversity entries: {len(diversity_metrics)}")

if __name__ == '__main__':
    main()
