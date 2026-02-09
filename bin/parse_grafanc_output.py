#!/usr/bin/env python3
"""
Parse GRAF-anc output into standardized ancestry proportions and categories.

GRAF-anc outputs per individual:
  - SNP count used for inference
  - GD1, GD2, GD3: continental genetic distances (to EUR, AFR, EAS centroids)
  - Subcontinental scores: EA1-4, AF1-3, EU1-3, SA1-2, IC1-3
  - Pe, Pf, Pa: European, African, East Asian ancestry proportions (barycentric)
  - AncGroupID: 8 continental + 38 subcontinental assignments

GRAF-anc categories (continental):
  100=AFR, 200=MEN, 300=EUR, 400=SAS, 500=EAS, 600=AMR, 700=OCN, 800=MIX

This script maps GRAF-anc output to the pipeline's standard format:
  - pct_european, pct_african, pct_amerindigenous, pct_east_asian, pct_south_asian
  - graf_category: EUR, AFR_AM, LA1, LA2, EAS, SAS (user-facing labels)
"""

import argparse
import pandas as pd
import sys

# GRAF-anc AncGroupID to pipeline category mapping
GRAFANC_CONTINENTAL_MAP = {
    100: 'AFR', 101: 'AFR', 102: 'AFR', 103: 'AFR', 104: 'AFR',
    105: 'AFR', 106: 'AFR', 107: 'AFR_AM', 108: 'AFR',
    200: 'MEN', 201: 'MEN', 202: 'MEN', 203: 'MEN',
    300: 'EUR', 301: 'EUR', 302: 'EUR', 303: 'EUR', 304: 'EUR',
    305: 'EUR', 306: 'EUR', 307: 'EUR', 308: 'EUR',
    400: 'SAS', 401: 'SAS', 402: 'SAS', 403: 'SAS', 404: 'SAS', 405: 'SAS',
    500: 'EAS', 501: 'EAS', 502: 'EAS', 503: 'EAS', 504: 'EAS',
    505: 'EAS', 506: 'EAS', 507: 'EAS', 508: 'EAS', 509: 'EAS',
    510: 'EAS', 511: 'EAS',
    600: 'AMR', 601: 'LA1', 602: 'LA2', 603: 'AMR',
    700: 'OCN',
    800: 'MIX',
}

def parse_args():
    parser = argparse.ArgumentParser(description='Parse GRAF-anc output')
    parser.add_argument('--input', required=True, help='GRAF-anc results file')
    parser.add_argument('--output', required=True, help='Standardized output TSV')
    return parser.parse_args()

def main():
    args = parse_args()

    try:
        df = pd.read_csv(args.input, sep='\t', comment='#')
    except Exception as e:
        print(f"Error reading GRAF-anc output: {e}", file=sys.stderr)
        # Create empty output
        pd.DataFrame(columns=['sample_id', 'pct_european', 'pct_african',
                              'pct_amerindigenous', 'pct_east_asian', 'pct_south_asian',
                              'graf_category', 'grafanc_subcontinental',
                              'n_snps_used', 'confidence']).to_csv(args.output, sep='\t', index=False)
        return

    results = []
    for _, row in df.iterrows():
        sample_id = str(row.iloc[0]).strip()

        # GRAF-anc Pe, Pf, Pa are European, African, East Asian barycentric proportions
        pe = float(row.get('Pe', row.get('P(European)', 0)))
        pf = float(row.get('Pf', row.get('P(African)', 0)))
        pa = float(row.get('Pa', row.get('P(E.Asian)', 0)))

        # Remaining proportion approximated for South Asian and Amerindigenous
        # GRAF-anc barycentric coords sum to ~1 for three continental references
        # South Asian and Amerindigenous are inferred from subcontinental scores
        remaining = max(0, 1.0 - pe - pf - pa)

        # Get subcontinental assignment
        anc_group_id = int(row.get('AncGroupID', row.get('Anc_Group', 800)))
        graf_cat = GRAFANC_CONTINENTAL_MAP.get(anc_group_id, 'MIX')
        subcont = str(row.get('Anc_Subgroup', anc_group_id))

        # Estimate 5-way proportions from 3-way + subcontinental info
        if graf_cat == 'SAS':
            pct_south_asian = remaining + pa * 0.3
            pct_east_asian = pa * 0.7
            pct_amerindigenous = 0.0
        elif graf_cat in ('AMR', 'LA1', 'LA2'):
            pct_amerindigenous = remaining + pa * 0.3
            pct_east_asian = pa * 0.3
            pct_south_asian = 0.0
        else:
            pct_south_asian = remaining * 0.5
            pct_east_asian = pa
            pct_amerindigenous = remaining * 0.5

        n_snps = int(row.get('Num_Snps', row.get('SNP_Count', 0)))

        results.append({
            'sample_id': sample_id,
            'pct_european': round(pe, 4),
            'pct_african': round(pf, 4),
            'pct_amerindigenous': round(pct_amerindigenous, 4),
            'pct_east_asian': round(pct_east_asian, 4),
            'pct_south_asian': round(pct_south_asian, 4),
            'graf_category': graf_cat,
            'grafanc_subcontinental': subcont,
            'grafanc_group_id': anc_group_id,
            'n_snps_used': n_snps,
            'GD1': row.get('GD1', ''),
            'GD2': row.get('GD2', ''),
            'GD3': row.get('GD3', ''),
        })

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"Parsed {len(results)} samples from GRAF-anc output")

    # Summary statistics
    if len(out_df) > 0:
        print(f"\nAncestry category distribution:")
        print(out_df['graf_category'].value_counts().to_string())
        print(f"\nMean SNPs used: {out_df['n_snps_used'].mean():.0f}")

if __name__ == '__main__':
    main()
