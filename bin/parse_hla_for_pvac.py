#!/usr/bin/env python3
"""
Parse arcasHLA genotype JSON to extract HLA alleles in pVACseq format.
Output: comma-separated HLA alleles (e.g., HLA-A*02:01,HLA-A*03:01,HLA-B*07:02,...)
"""

import json
import sys

def parse_arcashla_json(json_file):
    with open(json_file) as f:
        data = json.load(f)

    alleles = []
    for gene, typed_alleles in data.items():
        if isinstance(typed_alleles, list):
            for allele in typed_alleles:
                # arcasHLA format: A*02:01:01 -> HLA-A*02:01
                parts = allele.split(':')
                if len(parts) >= 2:
                    short_allele = ':'.join(parts[:2])
                    if not short_allele.startswith('HLA-'):
                        short_allele = f"HLA-{short_allele}"
                    alleles.append(short_allele)

    return ','.join(alleles)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: parse_hla_for_pvac.py <arcashla_genotype.json>", file=sys.stderr)
        sys.exit(1)

    alleles = parse_arcashla_json(sys.argv[1])
    print(alleles)
