#!/usr/bin/env python3
"""
Merge and filter neoantigen predictions from pVACseq and NeoFuse.
Calculate neoantigen burden per sample.
"""

import argparse
import os
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Filter and merge neoantigen predictions')
    parser.add_argument('--pvacseq-dir', default='.')
    parser.add_argument('--neofuse-dir', default='.')
    parser.add_argument('--binding-threshold', type=float, default=500)
    parser.add_argument('--output-merged', default='merged_neoantigens.tsv')
    parser.add_argument('--output-summary', default='neoantigen_summary.tsv')
    parser.add_argument('--output-burden', default='neoantigen_burden.tsv')
    return parser.parse_args()

def load_pvacseq_results(directory):
    """Load pVACseq filtered results."""
    results = []
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.endswith('.filtered.tsv'):
                filepath = os.path.join(root, f)
                try:
                    df = pd.read_csv(filepath, sep='\t')
                    df['source'] = 'pvacseq'
                    # Determine MHC class from path
                    if 'MHC_Class_I' in root:
                        df['mhc_class'] = 'I'
                    elif 'MHC_Class_II' in root:
                        df['mhc_class'] = 'II'
                    else:
                        df['mhc_class'] = 'unknown'
                    results.append(df)
                except Exception as e:
                    print(f"Warning: Could not parse {filepath}: {e}")
    return pd.concat(results) if results else pd.DataFrame()

def load_neofuse_results(directory):
    """Load NeoFuse results."""
    results = []
    for root, dirs, files in os.walk(directory):
        for f in files:
            if 'neoantigen' in f.lower() and f.endswith('.tsv'):
                filepath = os.path.join(root, f)
                try:
                    df = pd.read_csv(filepath, sep='\t')
                    df['source'] = 'neofuse'
                    results.append(df)
                except Exception as e:
                    print(f"Warning: Could not parse {filepath}: {e}")
    return pd.concat(results) if results else pd.DataFrame()

def main():
    args = parse_args()

    pvac = load_pvacseq_results(args.pvacseq_dir)
    neofuse = load_neofuse_results(args.neofuse_dir)

    all_neoantigens = pd.concat([pvac, neofuse], ignore_index=True)

    if len(all_neoantigens) == 0:
        pd.DataFrame().to_csv(args.output_merged, sep='\t', index=False)
        pd.DataFrame(columns=['sample_id', 'total_neoantigens', 'strong_binders', 'weak_binders']).to_csv(
            args.output_burden, sep='\t', index=False)
        return

    # Save merged
    all_neoantigens.to_csv(args.output_merged, sep='\t', index=False)

    # Calculate burden per sample
    sample_col = 'Sample' if 'Sample' in all_neoantigens.columns else 'sample_id'
    binding_col = next((c for c in all_neoantigens.columns if 'binding' in c.lower() or 'ic50' in c.lower()), None)

    if sample_col in all_neoantigens.columns:
        burden = all_neoantigens.groupby(sample_col).agg(
            total_neoantigens=('source', 'count')
        ).reset_index()

        if binding_col:
            strong = all_neoantigens[all_neoantigens[binding_col] <= 50]
            weak = all_neoantigens[(all_neoantigens[binding_col] > 50) & (all_neoantigens[binding_col] <= args.binding_threshold)]
            burden['strong_binders'] = all_neoantigens[all_neoantigens[binding_col] <= 50].groupby(sample_col).size().reindex(burden[sample_col]).fillna(0).values
            burden['weak_binders'] = weak.groupby(sample_col).size().reindex(burden[sample_col]).fillna(0).values

        burden.to_csv(args.output_burden, sep='\t', index=False)

    # Summary
    summary = pd.DataFrame({
        'metric': ['total_neoantigens', 'from_pvacseq', 'from_neofuse', 'n_samples'],
        'value': [len(all_neoantigens),
                  len(pvac) if len(pvac) > 0 else 0,
                  len(neofuse) if len(neofuse) > 0 else 0,
                  all_neoantigens[sample_col].nunique() if sample_col in all_neoantigens.columns else 0]
    })
    summary.to_csv(args.output_summary, sep='\t', index=False)

    print(f"Total neoantigens: {len(all_neoantigens)}")

if __name__ == '__main__':
    main()
