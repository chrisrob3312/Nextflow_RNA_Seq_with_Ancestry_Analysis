#!/usr/bin/env python3
"""
Merge and filter fusion calls from FusionCatcher and Arriba.
Identify high-confidence fusions supported by both callers.
"""

import argparse
import os
import pandas as pd
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Merge fusion results from multiple callers')
    parser.add_argument('--fusioncatcher-dir', default='.')
    parser.add_argument('--arriba-dir', default='.')
    parser.add_argument('--min-reads', type=int, default=3)
    parser.add_argument('--output-merged', default='merged_fusions.tsv')
    parser.add_argument('--output-hc', default='high_confidence_fusions.tsv')
    parser.add_argument('--output-summary', default='fusion_summary.tsv')
    return parser.parse_args()

def load_fusioncatcher(directory):
    """Load FusionCatcher results."""
    fusions = []
    for f in os.listdir(directory):
        if f.endswith('final-list_candidate-fusion-genes.txt'):
            sample_id = f.split('_fusioncatcher')[0]
            try:
                df = pd.read_csv(os.path.join(directory, f), sep='\t')
                df['sample_id'] = sample_id
                df['caller'] = 'fusioncatcher'
                # Standardize columns
                if 'Gene_1_symbol(5end_fusion_partner)' in df.columns:
                    df['gene5'] = df['Gene_1_symbol(5end_fusion_partner)']
                    df['gene3'] = df['Gene_2_symbol(3end_fusion_partner)']
                    df['fusion_id'] = df['gene5'] + '--' + df['gene3']
                fusions.append(df)
            except Exception as e:
                print(f"Warning: Could not parse {f}: {e}")
    return pd.concat(fusions) if fusions else pd.DataFrame()

def load_arriba(directory):
    """Load Arriba results."""
    fusions = []
    for f in os.listdir(directory):
        if f.endswith('.arriba.fusions.tsv'):
            sample_id = f.replace('.arriba.fusions.tsv', '')
            try:
                df = pd.read_csv(os.path.join(directory, f), sep='\t')
                df['sample_id'] = sample_id
                df['caller'] = 'arriba'
                if '#gene1' in df.columns:
                    df['gene5'] = df['#gene1']
                    df['gene3'] = df['gene2']
                    df['fusion_id'] = df['gene5'] + '--' + df['gene3']
                fusions.append(df)
            except Exception as e:
                print(f"Warning: Could not parse {f}: {e}")
    return pd.concat(fusions) if fusions else pd.DataFrame()

def main():
    args = parse_args()

    fc_fusions = load_fusioncatcher(args.fusioncatcher_dir)
    arriba_fusions = load_arriba(args.arriba_dir)

    # Merge all fusions
    all_fusions = pd.concat([fc_fusions, arriba_fusions], ignore_index=True)

    if len(all_fusions) == 0:
        pd.DataFrame(columns=['fusion_id', 'sample_id', 'gene5', 'gene3', 'callers']).to_csv(
            args.output_merged, sep='\t', index=False)
        pd.DataFrame(columns=['fusion_id', 'sample_id', 'gene5', 'gene3', 'callers']).to_csv(
            args.output_hc, sep='\t', index=False)
        return

    # Standardize fusion IDs (alphabetical gene pair)
    all_fusions['fusion_pair'] = all_fusions.apply(
        lambda row: '--'.join(sorted([str(row.get('gene5', '')), str(row.get('gene3', ''))])),
        axis=1
    )

    # Identify high-confidence: called by both tools in same sample
    caller_counts = all_fusions.groupby(['sample_id', 'fusion_pair'])['caller'].nunique().reset_index()
    caller_counts.columns = ['sample_id', 'fusion_pair', 'n_callers']

    hc = caller_counts[caller_counts['n_callers'] >= 2]
    hc_fusions = all_fusions.merge(hc[['sample_id', 'fusion_pair']], on=['sample_id', 'fusion_pair'])

    # Create summary
    summary = all_fusions.groupby('fusion_pair').agg(
        n_samples=('sample_id', 'nunique'),
        n_callers=('caller', 'nunique'),
        samples=('sample_id', lambda x: ','.join(x.unique())),
        callers=('caller', lambda x: ','.join(x.unique()))
    ).reset_index().sort_values('n_samples', ascending=False)

    # Save outputs
    cols_to_save = ['fusion_id', 'fusion_pair', 'sample_id', 'gene5', 'gene3', 'caller']
    cols_to_save = [c for c in cols_to_save if c in all_fusions.columns]

    all_fusions[cols_to_save].to_csv(args.output_merged, sep='\t', index=False)
    if len(hc_fusions) > 0:
        hc_fusions[cols_to_save].to_csv(args.output_hc, sep='\t', index=False)
    else:
        pd.DataFrame(columns=cols_to_save).to_csv(args.output_hc, sep='\t', index=False)
    summary.to_csv(args.output_summary, sep='\t', index=False)

    print(f"Total fusions: {len(all_fusions)}")
    print(f"High-confidence (multi-caller): {len(hc_fusions)}")
    print(f"Unique fusion pairs: {len(summary)}")

if __name__ == '__main__':
    main()
