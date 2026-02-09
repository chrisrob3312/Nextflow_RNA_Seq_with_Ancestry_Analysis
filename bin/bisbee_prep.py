#!/usr/bin/env python3
"""
Bisbee Prep: Convert SplAdder HDF5 count files to Bisbee-compatible TSV format.

Reads SplAdder merge_graphs_*.counts.hdf5 files and extracts inclusion/exclusion
junction counts per event per sample. Computes PSI (Percent Spliced In) values
and outputs one TSV file per event type for downstream Bisbee analysis.

Supported event types:
  - exon_skip: cassette exon skipping
  - intron_retention: intron retention
  - alt_3prime: alternative 3' splice site
  - alt_5prime: alternative 5' splice site
  - mult_exon_skip: multiple exon skipping
"""

import argparse
import os
import sys
import glob
import warnings

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EVENT_TYPES = [
    "exon_skip",
    "intron_retention",
    "alt_3prime",
    "alt_5prime",
    "mult_exon_skip",
]

# Mapping from SplAdder HDF5 dataset names to inclusion/exclusion count arrays.
# SplAdder stores counts in datasets named like:
#   event_counts  -> shape (n_events, n_samples, n_positions)
# The positions index encodes inclusion vs exclusion junctions depending on event type.
# For most event types the layout is:
#   inclusion junctions = positions 0,1  (flanking the included region)
#   exclusion junctions = position 2     (the skipping junction)
# Intron retention uses a different layout:
#   inclusion (intron retained) = position 0
#   exclusion (intron spliced)  = position 1
INCLUSION_EXCLUSION_MAP = {
    "exon_skip":       {"inc_idx": [0, 1], "exc_idx": [2]},
    "intron_retention": {"inc_idx": [0],    "exc_idx": [1]},
    "alt_3prime":      {"inc_idx": [0],    "exc_idx": [1]},
    "alt_5prime":      {"inc_idx": [0],    "exc_idx": [1]},
    "mult_exon_skip":  {"inc_idx": [0, 1], "exc_idx": [2]},
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare SplAdder HDF5 count files for Bisbee analysis"
    )
    parser.add_argument(
        "--spladder-dir",
        required=True,
        help="Directory containing SplAdder merge_graphs_*.counts.hdf5 files",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for per-event-type TSV files",
    )
    parser.add_argument(
        "--min-total-count",
        type=int,
        default=10,
        help="Minimum total junction count (inclusion + exclusion) to report an event-sample pair (default: 10)",
    )
    return parser.parse_args()


def compute_psi(inclusion, exclusion):
    """Compute PSI = inclusion / (inclusion + exclusion), returning NaN when denominator is 0."""
    total = inclusion + exclusion
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        psi = np.where(total > 0, inclusion / total, np.nan)
    return psi


def extract_event_ids_from_hdf5(h5file, event_type):
    """
    Extract event identifiers from the HDF5 file.

    SplAdder stores event coordinates in datasets such as 'event_pos' or
    'chrms'/'gene_names'. We build an event_id string from available fields.
    """
    event_ids = []

    # Try to read gene names and chromosome info
    gene_names = None
    chrms = None
    event_pos = None

    if "gene_names" in h5file:
        gene_names = h5file["gene_names"][:]
        if hasattr(gene_names[0], "decode"):
            gene_names = np.array([g.decode("utf-8") for g in gene_names])

    if "chrms" in h5file:
        chrms = h5file["chrms"][:]
        if hasattr(chrms[0], "decode"):
            chrms = np.array([c.decode("utf-8") for c in chrms])

    if "event_pos" in h5file:
        event_pos = h5file["event_pos"][:]

    n_events = h5file["event_counts"].shape[0]

    for i in range(n_events):
        parts = [event_type]
        if chrms is not None and i < len(chrms):
            parts.append(str(chrms[i]))
        if event_pos is not None and i < len(event_pos):
            coords = event_pos[i]
            # Flatten coordinates to a colon-separated string
            parts.append(":".join(str(int(c)) for c in coords if c > 0))
        if not event_pos is not None:
            parts.append(str(i))
        event_ids.append("_".join(parts))

    genes = []
    for i in range(n_events):
        if gene_names is not None and i < len(gene_names):
            genes.append(str(gene_names[i]))
        else:
            genes.append("unknown")

    return event_ids, genes


def extract_sample_names(h5file):
    """
    Extract sample identifiers from the HDF5 file.

    SplAdder typically stores them in a 'samples' or 'strains' dataset.
    """
    for key in ("samples", "strains", "sample_names"):
        if key in h5file:
            names = h5file[key][:]
            if hasattr(names[0], "decode"):
                names = [s.decode("utf-8") for s in names]
            else:
                names = [str(s) for s in names]
            return names

    # Fallback: generate numeric sample names
    n_samples = h5file["event_counts"].shape[1]
    return [f"sample_{i}" for i in range(n_samples)]


def process_hdf5_file(h5_path, event_type, min_total_count):
    """
    Read a single SplAdder HDF5 counts file and return a list of row dicts.

    Each row dict has keys: event_id, gene, sample_id, inclusion_count,
    exclusion_count, psi.
    """
    try:
        import h5py
    except ImportError:
        print("ERROR: h5py is required to read SplAdder HDF5 files. "
              "Install with: pip install h5py", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(h5_path):
        print(f"WARNING: HDF5 file not found: {h5_path}", file=sys.stderr)
        return []

    rows = []

    try:
        with h5py.File(h5_path, "r") as h5:
            if "event_counts" not in h5:
                print(f"WARNING: No 'event_counts' dataset in {h5_path}", file=sys.stderr)
                return []

            counts = h5["event_counts"][:]  # shape: (n_events, n_samples, n_positions)

            if counts.ndim != 3:
                print(
                    f"WARNING: Unexpected event_counts shape {counts.shape} in {h5_path}",
                    file=sys.stderr,
                )
                return []

            n_events, n_samples, n_positions = counts.shape

            if n_events == 0 or n_samples == 0:
                print(f"WARNING: Empty counts array in {h5_path}", file=sys.stderr)
                return []

            event_ids, gene_names = extract_event_ids_from_hdf5(h5, event_type)
            sample_names = extract_sample_names(h5)

            # Determine inclusion/exclusion index mapping
            ie_map = INCLUSION_EXCLUSION_MAP.get(event_type)
            if ie_map is None:
                print(f"WARNING: Unknown event type '{event_type}'", file=sys.stderr)
                return []

            # Validate that index positions exist in the data
            max_idx = max(max(ie_map["inc_idx"]), max(ie_map["exc_idx"]))
            if max_idx >= n_positions:
                # Fall back to simpler two-column layout
                inc_idx = [0]
                exc_idx = [min(1, n_positions - 1)]
            else:
                inc_idx = ie_map["inc_idx"]
                exc_idx = ie_map["exc_idx"]

            # Sum across inclusion and exclusion junction positions
            inclusion = counts[:, :, inc_idx].sum(axis=2).astype(np.float64)
            exclusion = counts[:, :, exc_idx].sum(axis=2).astype(np.float64)

            psi_matrix = compute_psi(inclusion, exclusion)

            # Build output rows, filtering by minimum total count
            for ev_i in range(n_events):
                for s_j in range(n_samples):
                    inc_val = int(inclusion[ev_i, s_j])
                    exc_val = int(exclusion[ev_i, s_j])
                    total = inc_val + exc_val

                    if total < min_total_count:
                        continue

                    rows.append({
                        "event_id": event_ids[ev_i],
                        "gene": gene_names[ev_i],
                        "sample_id": sample_names[s_j],
                        "inclusion_count": inc_val,
                        "exclusion_count": exc_val,
                        "psi": round(psi_matrix[ev_i, s_j], 6)
                        if not np.isnan(psi_matrix[ev_i, s_j])
                        else "NA",
                    })

    except Exception as e:
        print(f"ERROR: Failed to process {h5_path}: {e}", file=sys.stderr)
        return []

    return rows


def write_tsv(rows, output_path):
    """Write a list of row dicts to a TSV file."""
    header = ["event_id", "gene", "sample_id", "inclusion_count", "exclusion_count", "psi"]

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[col]) for col in header) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not os.path.isdir(args.spladder_dir):
        print(f"ERROR: SplAdder directory does not exist: {args.spladder_dir}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # Discover HDF5 count files
    h5_files = glob.glob(os.path.join(args.spladder_dir, "merge_graphs_*.counts.hdf5"))

    if not h5_files:
        print(f"WARNING: No merge_graphs_*.counts.hdf5 files found in {args.spladder_dir}",
              file=sys.stderr)
        # Write empty output files so downstream processes do not fail
        for event_type in EVENT_TYPES:
            write_tsv([], os.path.join(args.output_dir, f"{event_type}_counts.tsv"))
        print("Created empty output files for all event types.")
        return

    print(f"Found {len(h5_files)} HDF5 count file(s) in {args.spladder_dir}")

    files_processed = 0

    for event_type in EVENT_TYPES:
        # SplAdder naming convention: merge_graphs_{event_type}.counts.hdf5
        pattern = f"merge_graphs_{event_type}.counts.hdf5"
        matching = [f for f in h5_files if os.path.basename(f) == pattern]

        if not matching:
            print(f"  No HDF5 file found for event type '{event_type}', writing empty TSV.")
            write_tsv([], os.path.join(args.output_dir, f"{event_type}_counts.tsv"))
            continue

        h5_path = matching[0]
        print(f"  Processing {event_type}: {os.path.basename(h5_path)}")

        rows = process_hdf5_file(h5_path, event_type, args.min_total_count)
        output_path = os.path.join(args.output_dir, f"{event_type}_counts.tsv")
        write_tsv(rows, output_path)

        n_events = len(set(r["event_id"] for r in rows)) if rows else 0
        n_samples = len(set(r["sample_id"] for r in rows)) if rows else 0
        print(f"    -> {len(rows)} rows ({n_events} events x {n_samples} samples)")
        files_processed += 1

    print(f"\nBisbee prep complete. Processed {files_processed} event type(s).")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
