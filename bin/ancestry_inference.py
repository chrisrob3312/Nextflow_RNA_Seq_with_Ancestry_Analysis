#!/usr/bin/env python3
"""
Ancestry Inference from GRAF-pop SNP allele counts.
Uses supervised classification against reference panel to infer:
  - Continuous ancestry proportions (African, European, Amerindigenous, East Asian, South Asian)
  - Categorical GRAF ancestry assignments (EUR, AFR_AM, LA1, LA2, EAS, SAS)
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
import warnings
warnings.filterwarnings('ignore')

def parse_args():
    parser = argparse.ArgumentParser(description='Ancestry inference from GRAF SNP allele counts')
    parser.add_argument('--allele-counts', nargs='+', required=True, help='Allele count TSV files')
    parser.add_argument('--reference-panel', required=True, help='Reference panel VCF/TSV')
    parser.add_argument('--reference-labels', required=True, help='Population labels for reference')
    parser.add_argument('--categories', default='EUR,AFR_AM,LA1,LA2,EAS,SAS',
                       help='Comma-separated ancestry categories')
    parser.add_argument('--output-prefix', default='ancestry', help='Output file prefix')
    parser.add_argument('--plot-dir', default='ancestry_plots', help='Plot output directory')
    return parser.parse_args()

def load_allele_counts(files):
    """Load allele count files and build genotype matrix."""
    all_data = {}
    for f in files:
        sample_id = os.path.basename(f).replace('.graf_allele_counts.tsv', '')
        try:
            df = pd.read_csv(f, sep='\t', header=None,
                           names=['CHROM', 'POS', 'REF', 'ALT', 'DP', 'AD'])
            # Calculate allele frequency at each SNP
            afs = {}
            for _, row in df.iterrows():
                snp_id = f"{row['CHROM']}:{row['POS']}"
                try:
                    ads = str(row['AD']).split(',')
                    ref_count = int(ads[0])
                    alt_count = int(ads[1]) if len(ads) > 1 else 0
                    total = ref_count + alt_count
                    af = alt_count / total if total > 0 else np.nan
                    afs[snp_id] = af
                except (ValueError, IndexError):
                    afs[snp_id] = np.nan
            all_data[sample_id] = afs
        except Exception as e:
            print(f"Warning: Could not load {f}: {e}", file=sys.stderr)

    return pd.DataFrame(all_data).T

def infer_ancestry(genotype_matrix, ref_panel, ref_labels, categories):
    """
    Infer ancestry using PCA projection and supervised classification.
    Returns proportions (continuous) and categories (discrete).
    """
    categories_list = [c.strip() for c in categories.split(',')]

    # Align SNPs between study and reference
    common_snps = genotype_matrix.columns.intersection(ref_panel.columns)
    if len(common_snps) < 100:
        print(f"Warning: Only {len(common_snps)} common SNPs. Results may be unreliable.",
              file=sys.stderr)

    study_data = genotype_matrix[common_snps].fillna(0)
    ref_data = ref_panel[common_snps].fillna(0)

    # Combine for PCA
    combined = pd.concat([ref_data, study_data])
    scaler = StandardScaler()
    scaled = scaler.fit_transform(combined)

    # PCA
    n_components = min(20, scaled.shape[1], scaled.shape[0] - 1)
    pca = PCA(n_components=n_components)
    pca_coords = pca.fit_transform(scaled)

    ref_pca = pca_coords[:len(ref_data)]
    study_pca = pca_coords[len(ref_data):]

    # Supervised classification for categories
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(ref_pca, ref_labels['population'].values)

    predicted_categories = clf.predict(study_pca)
    predicted_probs = clf.predict_proba(study_pca)

    # Build proportions from class probabilities
    class_names = clf.classes_
    proportions = pd.DataFrame(predicted_probs, columns=class_names,
                                index=study_data.index)

    # Map to standard ancestry proportion columns
    ancestry_mapping = {
        'EUR': 'pct_european',
        'AFR_AM': 'pct_african',
        'AFR': 'pct_african',
        'LA1': 'pct_amerindigenous',
        'LA2': 'pct_amerindigenous',
        'AMR': 'pct_amerindigenous',
        'EAS': 'pct_east_asian',
        'SAS': 'pct_south_asian',
    }

    prop_df = pd.DataFrame(index=study_data.index)
    for std_col in ['pct_european', 'pct_african', 'pct_amerindigenous',
                    'pct_east_asian', 'pct_south_asian']:
        matching = [c for c, v in ancestry_mapping.items() if v == std_col and c in class_names]
        if matching:
            prop_df[std_col] = proportions[matching].sum(axis=1)
        else:
            prop_df[std_col] = 0.0

    # Categories
    cat_df = pd.DataFrame({
        'sample_id': study_data.index,
        'graf_category': predicted_categories,
        'confidence': predicted_probs.max(axis=1)
    })

    # PCA coordinates
    pca_df = pd.DataFrame(study_pca[:, :min(5, n_components)],
                           columns=[f'PC{i+1}' for i in range(min(5, n_components))],
                           index=study_data.index)
    pca_df.insert(0, 'sample_id', study_data.index)

    return prop_df, cat_df, pca_df

def create_plots(pca_df, cat_df, prop_df, plot_dir):
    """Generate ancestry visualization plots."""
    os.makedirs(plot_dir, exist_ok=True)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # PCA plot colored by category
        fig, ax = plt.subplots(figsize=(10, 8))
        merged = pca_df.merge(cat_df, on='sample_id')
        for cat in merged['graf_category'].unique():
            mask = merged['graf_category'] == cat
            ax.scatter(merged.loc[mask, 'PC1'], merged.loc[mask, 'PC2'],
                      label=cat, alpha=0.7, s=60)
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_title('Ancestry PCA - GRAF Categories')
        ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'ancestry_pca.png'), dpi=150)
        plt.close()

        # Ancestry proportion barplot
        fig, ax = plt.subplots(figsize=(14, 6))
        prop_df_sorted = prop_df.sort_values('pct_european')
        prop_df_sorted.plot(kind='bar', stacked=True, ax=ax, width=0.8)
        ax.set_ylabel('Ancestry Proportion')
        ax.set_title('Ancestry Proportions per Sample')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'ancestry_proportions.png'), dpi=150)
        plt.close()

        # Triangle plot (ternary)
        fig, ax = plt.subplots(figsize=(8, 8))
        if 'pct_european' in prop_df.columns and 'pct_african' in prop_df.columns:
            ax.scatter(prop_df['pct_european'], prop_df['pct_african'],
                      c=prop_df.get('pct_amerindigenous', 0), cmap='viridis',
                      alpha=0.7, s=60)
            ax.set_xlabel('European Proportion')
            ax.set_ylabel('African Proportion')
            ax.set_title('Ancestry Triangle Plot')
            plt.colorbar(ax.collections[0], label='Amerindigenous Proportion')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'ancestry_triangle.png'), dpi=150)
        plt.close()

    except ImportError:
        print("matplotlib not available. Skipping plots.", file=sys.stderr)

def main():
    args = parse_args()

    print("Loading allele counts...")
    genotype_matrix = load_allele_counts(args.allele_counts)
    print(f"Loaded {genotype_matrix.shape[0]} samples x {genotype_matrix.shape[1]} SNPs")

    print("Loading reference panel...")
    ref_panel = pd.read_csv(args.reference_panel, sep='\t', index_col=0)
    ref_labels = pd.read_csv(args.reference_labels, sep='\t')

    print("Inferring ancestry...")
    prop_df, cat_df, pca_df = infer_ancestry(
        genotype_matrix, ref_panel, ref_labels, args.categories
    )

    # Save outputs
    prop_df.insert(0, 'sample_id', prop_df.index)
    prop_df.to_csv(f'{args.output_prefix}_proportions.tsv', sep='\t', index=False)
    cat_df.to_csv(f'{args.output_prefix}_categories.tsv', sep='\t', index=False)
    pca_df.to_csv(f'{args.output_prefix}_pca.tsv', sep='\t', index=False)

    print("Creating plots...")
    create_plots(pca_df, cat_df, prop_df, args.plot_dir)

    print("Ancestry inference complete.")
    print(f"  Proportions: {args.output_prefix}_proportions.tsv")
    print(f"  Categories:  {args.output_prefix}_categories.tsv")
    print(f"  PCA coords:  {args.output_prefix}_pca.tsv")

if __name__ == '__main__':
    main()
