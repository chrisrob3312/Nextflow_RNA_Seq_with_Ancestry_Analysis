#!/usr/bin/env python3
"""
Bisbee Prot: Predict protein-level effects of differential splicing events.

Takes significant differential splicing results and a reference GTF + FASTA to
classify the protein impact of each splicing change:
  - Frameshift
  - In-frame insertion or deletion
  - NMD (nonsense-mediated decay) target
  - Domain disruption
  - Truncation

Input:
  - Differential splicing TSV from bisbee_diff.R
  - Reference GTF annotation
  - Reference genome FASTA

Output:
  - protein_effects.tsv with per-event protein impact classifications
"""

import argparse
import os
import sys
import re
import warnings


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Standard codon table (excluding stop codons)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}

# NMD rule: premature stop > 50 nt upstream of last exon-exon junction
NMD_THRESHOLD_NT = 50

# Effect type constants
EFFECT_FRAMESHIFT = "frameshift"
EFFECT_INFRAME_INSERTION = "in-frame_insertion"
EFFECT_INFRAME_DELETION = "in-frame_deletion"
EFFECT_NMD_TARGET = "NMD_target"
EFFECT_TRUNCATION = "truncation"
EFFECT_DOMAIN_DISRUPTION = "domain_disruption"
EFFECT_UNKNOWN = "unknown"
EFFECT_NO_CDS = "no_CDS_overlap"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Predict protein effects of differential splicing events"
    )
    parser.add_argument(
        "--diff-results",
        required=True,
        help="Differential splicing results TSV from bisbee_diff.R",
    )
    parser.add_argument(
        "--gtf",
        required=True,
        help="Reference GTF annotation file",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Reference genome FASTA file",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for protein effect predictions",
    )
    parser.add_argument(
        "--padj-threshold",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold for selecting significant events (default: 0.05)",
    )
    parser.add_argument(
        "--dpsi-threshold",
        type=float,
        default=0.1,
        help="Minimum absolute delta-PSI for significance (default: 0.1)",
    )
    parser.add_argument(
        "--domain-bed",
        default=None,
        help="Optional BED file of known protein domains (columns: chrom, start, end, domain_name, gene)",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

def parse_gtf_gene_structures(gtf_path):
    """
    Parse a GTF file to extract CDS and exon structures per gene/transcript.

    Returns:
        dict mapping gene_name -> list of transcript dicts, each with:
            - transcript_id: str
            - chrom: str
            - strand: str
            - exons: list of (start, end) tuples (0-based)
            - cds: list of (start, end) tuples (0-based)
    """
    genes = {}

    if not os.path.isfile(gtf_path):
        print(f"ERROR: GTF file not found: {gtf_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing GTF: {gtf_path}")

    with open(gtf_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type not in ("exon", "CDS"):
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])        # End is already exclusive in 0-based
            strand = fields[6]

            # Parse attributes
            attrs = fields[8]
            gene_name = _extract_attribute(attrs, "gene_name")
            transcript_id = _extract_attribute(attrs, "transcript_id")

            if gene_name is None or transcript_id is None:
                continue

            if gene_name not in genes:
                genes[gene_name] = {}

            if transcript_id not in genes[gene_name]:
                genes[gene_name][transcript_id] = {
                    "transcript_id": transcript_id,
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                    "cds": [],
                }

            tx = genes[gene_name][transcript_id]

            if feature_type == "exon":
                tx["exons"].append((start, end))
            elif feature_type == "CDS":
                tx["cds"].append((start, end))

    # Sort exons and CDS by start position
    for gene_name in genes:
        for tx_id in genes[gene_name]:
            tx = genes[gene_name][tx_id]
            tx["exons"].sort()
            tx["cds"].sort()

    n_genes = len(genes)
    n_tx = sum(len(txs) for txs in genes.values())
    print(f"  Parsed {n_genes} genes, {n_tx} transcripts")

    return genes


def _extract_attribute(attr_string, key):
    """Extract a value for a given key from a GTF attribute string."""
    # Pattern: key "value"; or key "value"
    pattern = rf'{key}\s+"([^"]+)"'
    match = re.search(pattern, attr_string)
    if match:
        return match.group(1)
    return None


# ---------------------------------------------------------------------------
# FASTA reading
# ---------------------------------------------------------------------------

def load_fasta_index(fasta_path):
    """
    Load a FASTA file via pysam if available, else build a simple dict.

    Returns a callable that takes (chrom, start, end) and returns sequence.
    """
    # Try pysam first (handles .fai indexed FASTA efficiently)
    try:
        import pysam
        fa = pysam.FastaFile(fasta_path)

        def fetch_seq(chrom, start, end):
            try:
                return fa.fetch(chrom, start, end).upper()
            except (KeyError, ValueError):
                return None

        return fetch_seq
    except ImportError:
        pass

    # Fallback: load entire FASTA into memory (for small genomes or testing)
    print("WARNING: pysam not available; loading FASTA into memory (may be slow)", file=sys.stderr)

    sequences = {}
    current_chrom = None
    current_seq = []

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_chrom is not None:
                    sequences[current_chrom] = "".join(current_seq)
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_chrom is not None:
            sequences[current_chrom] = "".join(current_seq)

    def fetch_seq(chrom, start, end):
        seq = sequences.get(chrom)
        if seq is None:
            return None
        return seq[start:end]

    return fetch_seq


# ---------------------------------------------------------------------------
# Domain loading
# ---------------------------------------------------------------------------

def load_domains(domain_bed_path):
    """
    Load protein domain BED file.

    Returns dict mapping gene_name -> list of (chrom, start, end, domain_name).
    """
    if domain_bed_path is None or not os.path.isfile(domain_bed_path):
        return {}

    domains = {}
    with open(domain_bed_path, "r") as fh:
        for line in fh:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            domain_name = fields[3]
            gene_name = fields[4] if len(fields) > 4 else "unknown"

            if gene_name not in domains:
                domains[gene_name] = []
            domains[gene_name].append((chrom, start, end, domain_name))

    return domains


# ---------------------------------------------------------------------------
# Event coordinate parsing
# ---------------------------------------------------------------------------

def parse_event_coordinates(event_id):
    """
    Parse a Bisbee event ID to extract chromosome and coordinates.

    Event IDs from bisbee_prep.py have the format:
        {event_type}_{chrom}_{coord1}:{coord2}:...

    Returns:
        (event_type, chrom, coords_list) or (event_type, None, []) if unparseable.
    """
    parts = event_id.split("_", 2)
    if len(parts) < 2:
        return (event_id, None, [])

    # Event type may be multi-word (e.g., exon_skip, alt_3prime, mult_exon_skip)
    # Try to match known event types
    for n_underscores in [2, 1]:
        candidate = "_".join(parts[:n_underscores + 1]) if len(parts) > n_underscores else None
        if candidate and candidate in (
            "exon_skip", "intron_retention", "alt_3prime", "alt_5prime", "mult_exon_skip"
        ):
            remainder = "_".join(parts[n_underscores + 1:])
            break
    else:
        # Fallback: first part is event type
        candidate = parts[0]
        remainder = "_".join(parts[1:])

    # Parse chrom and coordinates from remainder
    # Expected: chr1_12345:23456:34567 or chr1_12345:23456
    rem_parts = remainder.split("_", 1)
    chrom = rem_parts[0] if rem_parts else None

    coords = []
    if len(rem_parts) > 1:
        coord_str = rem_parts[1]
        try:
            coords = [int(c) for c in coord_str.split(":") if c]
        except ValueError:
            coords = []

    return (candidate, chrom, coords)


# ---------------------------------------------------------------------------
# Protein effect prediction
# ---------------------------------------------------------------------------

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def translate_sequence(dna_seq):
    """Translate a DNA sequence to amino acids. Returns (protein_seq, has_stop)."""
    protein = []
    has_stop = False
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i + 3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            has_stop = True
            break
        protein.append(aa)
    return "".join(protein), has_stop


def predict_effect_for_event(event_id, gene, event_type, gene_structures,
                              fetch_seq, domains, chrom, coords):
    """
    Predict the protein-level effect of a single splicing event.

    Returns a dict with: effect_type, affected_domains, nmd_prediction, protein_change
    """
    result = {
        "event_id": event_id,
        "gene": gene,
        "effect_type": EFFECT_UNKNOWN,
        "affected_domains": "none",
        "nmd_prediction": "NA",
        "protein_change": "NA",
    }

    # Look up gene structure
    if gene not in gene_structures or not gene_structures[gene]:
        result["effect_type"] = EFFECT_NO_CDS
        return result

    # Use the longest CDS transcript as representative
    transcripts = gene_structures[gene]
    best_tx = None
    best_cds_len = 0
    for tx_id, tx in transcripts.items():
        cds_len = sum(e - s for s, e in tx["cds"])
        if cds_len > best_cds_len:
            best_cds_len = cds_len
            best_tx = tx

    if best_tx is None or not best_tx["cds"]:
        result["effect_type"] = EFFECT_NO_CDS
        return result

    tx_chrom = best_tx["chrom"]
    tx_strand = best_tx["strand"]
    exons = best_tx["exons"]
    cds_regions = best_tx["cds"]

    # If we have coordinates, compute the affected region length
    if len(coords) < 2:
        result["effect_type"] = EFFECT_UNKNOWN
        result["protein_change"] = "insufficient_coordinates"
        return result

    # Determine the size of the alternatively spliced region
    # For exon_skip: the skipped exon is defined by coords
    # For intron_retention: the retained intron region
    # For alt_3/5prime: the alternative region between splice sites
    if event_type in ("exon_skip", "mult_exon_skip"):
        # Spliced region = the exon(s) being skipped
        # coords typically: [exon_start, exon_end, ...]
        if len(coords) >= 2:
            spliced_region_length = coords[1] - coords[0]
        else:
            spliced_region_length = 0
    elif event_type == "intron_retention":
        if len(coords) >= 2:
            spliced_region_length = coords[1] - coords[0]
        else:
            spliced_region_length = 0
    elif event_type in ("alt_3prime", "alt_5prime"):
        if len(coords) >= 2:
            spliced_region_length = abs(coords[1] - coords[0])
        else:
            spliced_region_length = 0
    else:
        spliced_region_length = 0

    # Check if the event overlaps CDS
    event_start = min(coords[:2]) if len(coords) >= 2 else 0
    event_end = max(coords[:2]) if len(coords) >= 2 else 0
    overlaps_cds = any(
        s < event_end and e > event_start for s, e in cds_regions
    )

    if not overlaps_cds:
        result["effect_type"] = EFFECT_NO_CDS
        result["protein_change"] = "UTR_or_non-coding"
        return result

    # Determine frameshift vs in-frame based on region length mod 3
    if spliced_region_length % 3 != 0:
        result["effect_type"] = EFFECT_FRAMESHIFT

        # Frameshift events are likely NMD targets: check if a premature stop
        # would occur more than NMD_THRESHOLD_NT upstream of the last EJC
        # Simplified heuristic: frameshifts in internal exons are likely NMD targets
        if len(exons) > 1:
            last_ejc_pos = exons[-2][1] if tx_strand == "+" else exons[1][0]
            # If event is sufficiently upstream of last EJC, predict NMD
            if tx_strand == "+":
                dist_to_last_ejc = last_ejc_pos - event_end
            else:
                dist_to_last_ejc = event_start - exons[1][0]

            if dist_to_last_ejc > NMD_THRESHOLD_NT:
                result["nmd_prediction"] = "likely_NMD"
            else:
                result["nmd_prediction"] = "may_escape_NMD"
        else:
            result["nmd_prediction"] = "single_exon_gene"

        # Attempt to compute the actual protein change
        if fetch_seq is not None and chrom is not None:
            seq = fetch_seq(tx_chrom, event_start, event_end)
            if seq is not None:
                if tx_strand == "-":
                    seq = reverse_complement(seq)
                aa_len = spliced_region_length // 3
                result["protein_change"] = f"frameshift_affecting_{aa_len}+_codons"

    else:
        # In-frame change
        aa_change = spliced_region_length // 3

        if event_type in ("exon_skip", "mult_exon_skip"):
            result["effect_type"] = EFFECT_INFRAME_DELETION
            result["protein_change"] = f"deletion_of_{aa_change}_aa"
        elif event_type == "intron_retention":
            # Intron retention: check for in-frame stop codons
            if fetch_seq is not None and chrom is not None:
                intron_seq = fetch_seq(tx_chrom, event_start, event_end)
                if intron_seq is not None:
                    if tx_strand == "-":
                        intron_seq = reverse_complement(intron_seq)
                    _, has_stop = translate_sequence(intron_seq)
                    if has_stop:
                        result["effect_type"] = EFFECT_TRUNCATION
                        result["protein_change"] = f"premature_stop_in_retained_intron"
                        result["nmd_prediction"] = "likely_NMD"
                    else:
                        result["effect_type"] = EFFECT_INFRAME_INSERTION
                        result["protein_change"] = f"insertion_of_{aa_change}_aa"
                else:
                    result["effect_type"] = EFFECT_INFRAME_INSERTION
                    result["protein_change"] = f"insertion_of_{aa_change}_aa"
            else:
                result["effect_type"] = EFFECT_INFRAME_INSERTION
                result["protein_change"] = f"insertion_of_{aa_change}_aa"
        elif event_type in ("alt_3prime", "alt_5prime"):
            if spliced_region_length > 0:
                result["effect_type"] = EFFECT_INFRAME_INSERTION
                result["protein_change"] = f"alt_region_{aa_change}_aa"
            else:
                result["effect_type"] = EFFECT_INFRAME_DELETION
                result["protein_change"] = f"alt_region_{aa_change}_aa"
        else:
            result["effect_type"] = EFFECT_INFRAME_DELETION if aa_change > 0 else EFFECT_UNKNOWN
            result["protein_change"] = f"in-frame_change_{aa_change}_aa"

        result["nmd_prediction"] = "unlikely_NMD" if result["nmd_prediction"] == "NA" else result["nmd_prediction"]

    # Check domain disruption
    affected = check_domain_disruption(gene, event_start, event_end, chrom, domains)
    if affected:
        if result["effect_type"] in (EFFECT_INFRAME_DELETION, EFFECT_INFRAME_INSERTION):
            result["effect_type"] = EFFECT_DOMAIN_DISRUPTION
        result["affected_domains"] = ";".join(affected)

    return result


def check_domain_disruption(gene, event_start, event_end, chrom, domains):
    """Check if the splicing event overlaps any known protein domains."""
    if not domains or gene not in domains:
        return []

    affected = []
    for d_chrom, d_start, d_end, d_name in domains[gene]:
        if chrom and d_chrom != chrom:
            continue
        # Check overlap
        if d_start < event_end and d_end > event_start:
            affected.append(d_name)
    return affected


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # Validate inputs
    if not os.path.isfile(args.diff_results):
        print(f"ERROR: Differential results file not found: {args.diff_results}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.gtf):
        print(f"ERROR: GTF file not found: {args.gtf}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.fasta):
        print(f"ERROR: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # ---- Load differential splicing results --------------------------------
    try:
        import pandas as pd
        diff_df = pd.read_csv(args.diff_results, sep="\t")
    except ImportError:
        # Fallback without pandas
        diff_df = _read_tsv_simple(args.diff_results)
    except Exception as e:
        print(f"ERROR: Could not read {args.diff_results}: {e}", file=sys.stderr)
        sys.exit(1)

    if isinstance(diff_df, dict):
        # Simple dict-of-lists format from fallback reader
        n_rows = len(diff_df.get("event_id", []))
    else:
        n_rows = len(diff_df)

    if n_rows == 0:
        print("WARNING: No events in differential splicing results.")
        _write_empty_output(args.output_dir)
        return

    print(f"Loaded {n_rows} events from differential splicing results")

    # Filter to significant events
    if isinstance(diff_df, dict):
        sig_indices = []
        for i in range(n_rows):
            try:
                padj = float(diff_df["padj"][i]) if diff_df["padj"][i] not in ("NA", "NaN", "") else 1.0
                dpsi = abs(float(diff_df["deltaPSI"][i])) if diff_df["deltaPSI"][i] not in ("NA", "NaN", "") else 0.0
            except (ValueError, KeyError):
                continue
            if padj < args.padj_threshold and dpsi >= args.dpsi_threshold:
                sig_indices.append(i)

        sig_events = {k: [v[i] for i in sig_indices] for k, v in diff_df.items()}
        n_sig = len(sig_indices)
    else:
        # pandas DataFrame
        diff_df["padj"] = pd.to_numeric(diff_df["padj"], errors="coerce")
        diff_df["deltaPSI"] = pd.to_numeric(diff_df["deltaPSI"], errors="coerce")
        mask = (diff_df["padj"] < args.padj_threshold) & (diff_df["deltaPSI"].abs() >= args.dpsi_threshold)
        sig_events = diff_df[mask].reset_index(drop=True)
        n_sig = len(sig_events)

    if n_sig == 0:
        print("No significant events passed thresholds.")
        _write_empty_output(args.output_dir)
        return

    print(f"Processing {n_sig} significant events for protein effects")

    # ---- Load gene structures from GTF -------------------------------------
    gene_structures = parse_gtf_gene_structures(args.gtf)

    # ---- Load FASTA --------------------------------------------------------
    print(f"Loading FASTA: {args.fasta}")
    fetch_seq = load_fasta_index(args.fasta)

    # ---- Load domains (optional) -------------------------------------------
    domains = load_domains(args.domain_bed)
    if domains:
        print(f"Loaded domains for {len(domains)} genes")

    # ---- Predict protein effects -------------------------------------------
    print("Predicting protein effects...")
    results = []

    if isinstance(sig_events, dict):
        for i in range(n_sig):
            event_id = sig_events["event_id"][i]
            gene = sig_events["gene"][i]
            event_type = sig_events.get("event_type", ["unknown"] * n_sig)[i]

            event_type_parsed, chrom, coords = parse_event_coordinates(event_id)
            if event_type == "unknown" and event_type_parsed:
                event_type = event_type_parsed

            result = predict_effect_for_event(
                event_id, gene, event_type, gene_structures,
                fetch_seq, domains, chrom, coords
            )
            results.append(result)
    else:
        for _, row in sig_events.iterrows():
            event_id = str(row["event_id"])
            gene = str(row["gene"])
            event_type = str(row.get("event_type", "unknown"))

            event_type_parsed, chrom, coords = parse_event_coordinates(event_id)
            if event_type == "unknown" and event_type_parsed:
                event_type = event_type_parsed

            result = predict_effect_for_event(
                event_id, gene, event_type, gene_structures,
                fetch_seq, domains, chrom, coords
            )
            results.append(result)

    # ---- Write output ------------------------------------------------------
    output_path = os.path.join(args.output_dir, "protein_effects.tsv")
    header = ["event_id", "gene", "effect_type", "affected_domains", "nmd_prediction", "protein_change"]

    with open(output_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in results:
            fh.write("\t".join(str(r[col]) for col in header) + "\n")

    # Print summary
    effect_counts = {}
    for r in results:
        et = r["effect_type"]
        effect_counts[et] = effect_counts.get(et, 0) + 1

    print(f"\nProtein effect predictions written to: {output_path}")
    print(f"Total events processed: {len(results)}")
    print("Effect type summary:")
    for et, count in sorted(effect_counts.items(), key=lambda x: -x[1]):
        print(f"  {et}: {count}")

    nmd_count = sum(1 for r in results if r["nmd_prediction"] == "likely_NMD")
    domain_count = sum(1 for r in results if r["affected_domains"] != "none")
    print(f"Predicted NMD targets: {nmd_count}")
    print(f"Events affecting known domains: {domain_count}")


def _read_tsv_simple(path):
    """Simple TSV reader returning dict of lists (fallback when pandas unavailable)."""
    data = {}
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for col in header:
            data[col] = []
        for line in fh:
            if line.strip() == "":
                continue
            fields = line.rstrip("\n").split("\t")
            for i, col in enumerate(header):
                data[col].append(fields[i] if i < len(fields) else "")
    return data


def _write_empty_output(output_dir):
    """Write empty protein_effects.tsv output."""
    header = ["event_id", "gene", "effect_type", "affected_domains", "nmd_prediction", "protein_change"]
    output_path = os.path.join(output_dir, "protein_effects.tsv")
    with open(output_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
    print(f"Empty output written to: {output_path}")


if __name__ == "__main__":
    main()
