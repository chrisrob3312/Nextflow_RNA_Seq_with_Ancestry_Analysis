#!/usr/bin/env python3
"""
Validate pipeline samplesheet CSV.
Checks: required columns, file existence, valid values, duplicate sample IDs.
"""

import argparse
import csv
import os
import sys

REQUIRED_COLUMNS = ['sample_id', 'bam', 'bai']
OPTIONAL_COLUMNS = [
    'batch', 'sex', 'age', 'tumor_purity', 'cytomolecular_subgroup',
    'relapse_status', 'adi_quartile', 'timepoint', 'disease_stage',
    'blast_percentage'
]
VALID_SEX = {'M', 'F', 'NA', ''}
VALID_RELAPSE = {'relapse', 'no_relapse', 'NA', ''}
VALID_TIMEPOINT = {'diagnostic', 'relapse', 'NA', ''}

def validate(samplesheet_path, check_files=True):
    errors = []
    warnings = []
    sample_ids = set()

    with open(samplesheet_path) as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames

        # Check required columns
        for col in REQUIRED_COLUMNS:
            if col not in headers:
                errors.append(f"Missing required column: {col}")

        if errors:
            return errors, warnings

        for i, row in enumerate(reader, start=2):
            sid = row['sample_id']

            # Check duplicate sample IDs
            if sid in sample_ids:
                errors.append(f"Row {i}: Duplicate sample_id '{sid}'")
            sample_ids.add(sid)

            # Check BAM/BAI files exist
            if check_files:
                for col in ['bam', 'bai']:
                    if not os.path.exists(row[col]):
                        errors.append(f"Row {i} ({sid}): File not found: {row[col]}")

            # Validate sex
            if 'sex' in row and row['sex'] not in VALID_SEX:
                warnings.append(f"Row {i} ({sid}): Unexpected sex value: {row['sex']}")

            # Validate numeric fields
            for num_col in ['age', 'tumor_purity', 'blast_percentage']:
                if num_col in row and row[num_col] and row[num_col] != 'NA':
                    try:
                        val = float(row[num_col])
                        if num_col == 'tumor_purity' and not (0 <= val <= 1):
                            warnings.append(f"Row {i} ({sid}): tumor_purity should be 0-1, got {val}")
                        if num_col == 'blast_percentage' and not (0 <= val <= 100):
                            warnings.append(f"Row {i} ({sid}): blast_percentage should be 0-100, got {val}")
                    except ValueError:
                        errors.append(f"Row {i} ({sid}): Non-numeric value in {num_col}: {row[num_col]}")

            # Validate relapse_status
            if 'relapse_status' in row and row['relapse_status'] not in VALID_RELAPSE:
                warnings.append(f"Row {i} ({sid}): Unexpected relapse_status: {row['relapse_status']}")

            # Validate timepoint
            if 'timepoint' in row and row['timepoint'] not in VALID_TIMEPOINT:
                warnings.append(f"Row {i} ({sid}): Unexpected timepoint: {row['timepoint']}")

    return errors, warnings

def main():
    parser = argparse.ArgumentParser(description='Validate pipeline samplesheet')
    parser.add_argument('samplesheet', help='Path to samplesheet CSV')
    parser.add_argument('--no-check-files', action='store_true', help='Skip file existence checks')
    args = parser.parse_args()

    errors, warnings = validate(args.samplesheet, check_files=not args.no_check_files)

    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)

    if errors:
        print(f"\nValidation FAILED with {len(errors)} errors and {len(warnings)} warnings.",
              file=sys.stderr)
        sys.exit(1)
    else:
        print(f"Validation PASSED ({len(warnings)} warnings).", file=sys.stderr)

if __name__ == '__main__':
    main()
