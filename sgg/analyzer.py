"""Core analysis module - sequence parsing, genotype identification, pathogenicity, trees."""

import csv
import io
import json
import os
import shutil
import subprocess
import zlib
from typing import ClassVar

from Bio import Phylo

from sgg.palette import PALETTE

# ============================================================================
# External tool lookup (mafft, FastTree, iqtree2 - all provided as AUR/pacman deps)
# ============================================================================


def _require(cmd_names: list[str], package_hint: str) -> str:
    for name in cmd_names:
        found = shutil.which(name)
        if found:
            return found
    raise SystemExit(
        f"Could not find {'/'.join(cmd_names)} on PATH. Install it, e.g. `pacman -S {package_hint}`."
    )


def get_mafft_cmd() -> str:
    return _require(["mafft"], "mafft")


def get_fasttree_cmd() -> str:
    return _require(["FastTree", "fasttree"], "fasttree")


def get_iqtree_cmd() -> str:
    return _require(["iqtree2", "iqtree"], "iqtree")


# ============================================================================
# Entry / gene configuration
# ============================================================================


def load_gene_config(gene_path: str) -> dict:
    config = {}
    config_path = os.path.join(gene_path, "_config.json")
    if os.path.exists(config_path):
        with open(config_path) as f:
            raw = json.load(f)
        if "cleavage_start" in raw:
            config["cleavage_start"] = raw["cleavage_start"]
        if "genotype_pattern" in raw:
            config["genotype_pattern"] = raw["genotype_pattern"]

    motifs_files = [f for f in os.listdir(gene_path) if f.endswith("_motifs.csv")]
    if motifs_files:
        motifs_by_type = {}
        with open(os.path.join(gene_path, motifs_files[0]), newline="") as f:
            for row in csv.DictReader(f):
                motif = row["motif"].strip().upper()
                label = row["label"].strip()
                type_name = row["type"].strip().lower()
                if type_name not in motifs_by_type:
                    motifs_by_type[type_name] = {}
                motifs_by_type[type_name][motif] = label
        if motifs_by_type:
            config["motifs_by_type"] = motifs_by_type
    return config


def load_entry_config(entry_path: str) -> dict:
    genes = [
        f for f in os.listdir(entry_path) if os.path.isdir(os.path.join(entry_path, f))
    ]

    if genes:
        return {
            "multi": True,
            "genes": {
                gene: load_gene_config(os.path.join(entry_path, gene))
                for gene in sorted(genes)
            },
        }
    return load_gene_config(entry_path)


_FASTA_EXTS = {".fasta", ".fas", ".fa", ".txt"}


def load_all_references(path: str):
    """Load reference FASTA(s) from an entry directory.

    Returns (references, files_count, total_count, errors). For a multi-gene entry,
    `references` is {gene: {header: seq}}; for mono-gene it's {header: seq} directly.
    """
    subdirs = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
    if subdirs:
        genes_ref = {}
        genes_count = {}
        all_errors = []
        for gene in sorted(subdirs):
            gene_path = os.path.join(path, gene)
            combined = {}
            errors = []
            fas_files = [
                f
                for f in os.listdir(gene_path)
                if os.path.splitext(f)[1].lower() in _FASTA_EXTS
            ]
            for filename in fas_files:
                try:
                    seqs = FASTAParser.parse_file(os.path.join(gene_path, filename))
                    if seqs:
                        combined.update(seqs)
                except (OSError, ValueError, UnicodeDecodeError) as e:
                    errors.append(f"{gene}/{filename}: {e}")
            genes_ref[gene] = combined
            genes_count[gene] = len(combined)
            all_errors.extend(errors)
        return genes_ref, len(subdirs), sum(genes_count.values()), all_errors

    combined_sequences = {}
    errors = []
    fas_files = [
        f for f in os.listdir(path) if os.path.splitext(f)[1].lower() in _FASTA_EXTS
    ]
    for filename in fas_files:
        try:
            seqs = FASTAParser.parse_file(os.path.join(path, filename))
            if seqs:
                combined_sequences.update(seqs)
        except (OSError, ValueError, UnicodeDecodeError) as e:
            errors.append(f"{filename}: {e}")
    return combined_sequences, len(fas_files), len(combined_sequences), errors


# ============================================================================
# FASTA parsing
# ============================================================================


class FASTAParser:
    @staticmethod
    def parse_file(filepath: str) -> dict[str, str]:
        sequences = {}
        current_header = None
        current_seq = []

        try:
            with open(filepath, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith(">"):
                        if current_header:
                            sequences[current_header] = "".join(current_seq)
                        current_header = line[1:]
                        current_seq = []
                    else:
                        current_seq.append(line.upper())

                if current_header:
                    sequences[current_header] = "".join(current_seq)

        except (OSError, ValueError, UnicodeDecodeError) as e:
            print(f"Error parsing FASTA file: {e}")
            return {}

        return sequences

    @staticmethod
    def parse_text(fasta_text: str) -> dict[str, str]:
        sequences = {}
        lines = fasta_text.strip().split("\n")

        current_header = None
        current_seq = []

        for line in lines:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(
                    "".join(c for c in line.upper() if c in "ATCGNRYSWKMBDHV")
                )

        if current_header:
            sequences[current_header] = "".join(current_seq)

        return sequences


# ============================================================================
# Similarity
# ============================================================================


class SequenceSimilarity:
    @staticmethod
    def hamming_distance(seq1: str, seq2: str) -> float:
        if len(seq1) != len(seq2):
            min_len = min(len(seq1), len(seq2))
            seq1 = seq1[:min_len]
            seq2 = seq2[:min_len]

        if len(seq1) == 0:
            return 0.0

        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return round((matches / len(seq1)) * 100, 2)

    @staticmethod
    def pairwise_similarity(seq1: str, seq2: str) -> float:
        from Bio.Align import PairwiseAligner

        if len(seq1) == 0 and len(seq2) == 0:
            return 100.0
        if len(seq1) == 0 or len(seq2) == 0:
            return 0.0

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = 0
        aligner.open_gap_score = -1
        aligner.extend_gap_score = 0

        score = aligner.score(seq1, seq2)
        max_score = max(len(seq1), len(seq2))
        return round((score / max_score) * 100, 2)


# ============================================================================
# Pathogenicity / cleavage site analysis
# ============================================================================


class CleavageSiteAnalyzer:
    CLEAVAGE_START = 333  # 0-indexed
    CLEAVAGE_LENGTH = 24
    WINDOW = 30

    CODON_TABLE: ClassVar[dict[str, str]] = {
        "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
        "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
        "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
        "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
        "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
        "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
        "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
        "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
        "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
        "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
        "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
        "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
        "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
        "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
        "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
    }
    # https://gist.github.com/juanfal/09d7fb53bd367742127e17284b9c47bf

    @staticmethod
    def translate_to_protein(dna_sequence: str) -> str:
        protein = []
        for i in range(0, len(dna_sequence) - 2, 3):
            codon = dna_sequence[i : i + 3].upper()
            if "N" in codon or len(codon) < 3:
                protein.append("X")
            else:
                protein.append(CleavageSiteAnalyzer.CODON_TABLE.get(codon, "X"))
        return "".join(protein)

    @staticmethod
    def motif_matches(motif: str, sequence: str) -> bool:
        if len(motif) > len(sequence):
            return False
        for i in range(len(sequence) - len(motif) + 1):
            match = True
            for j, aa in enumerate(motif):
                if aa != "X" and sequence[i + j] != aa:
                    match = False
                    break
            if match:
                return True
        return False

    @staticmethod
    def analyze(
        sequence: str,
        cleavage_start: int,
        motifs_by_type: dict | None = None,
    ):
        motifs_by_type = motifs_by_type or {}

        WINDOW = 29 * 3  # tolerance for indels

        result = {
            "cleavage_region_found": False,
            "cleavage_nucleotides": None,
            "cleavage_protein": None,
            "pathogenicity": "Undetermined",
            "motif_type": "No known motif found in cleavage region",
            "motif_category": None,
        }

        result_plus_one = result.copy()
        result_minus_one = result.copy()

        if len(sequence) < cleavage_start - WINDOW:
            return result, result_plus_one, result_minus_one

        region_start = max(0, cleavage_start - WINDOW)
        region_end = min(len(sequence), cleavage_start + WINDOW)
        cleavage_region_nuc = sequence[region_start:region_end]

        region_start_plus_one = max(0, cleavage_start + 1 - WINDOW)
        region_end_plus_one = min(len(sequence), cleavage_start + 1 + WINDOW)
        cleavage_region_nuc_plus_one = sequence[region_start_plus_one:region_end_plus_one]

        region_start_minus_one = max(0, cleavage_start - 1 - WINDOW)
        region_end_minus_one = min(len(sequence), cleavage_start - 1 + WINDOW)
        cleavage_region_nuc_minus_one = sequence[region_start_minus_one:region_end_minus_one]

        cleavage_region_prot = CleavageSiteAnalyzer.translate_to_protein(cleavage_region_nuc)
        cleavage_region_prot_plus_one = CleavageSiteAnalyzer.translate_to_protein(
            cleavage_region_nuc_plus_one
        )
        cleavage_region_prot_minus_one = CleavageSiteAnalyzer.translate_to_protein(
            cleavage_region_nuc_minus_one
        )

        result["cleavage_region_found"] = True
        result["cleavage_nucleotides"] = cleavage_region_nuc
        result["cleavage_protein"] = cleavage_region_prot

        result_plus_one["cleavage_region_found"] = True
        result_plus_one["cleavage_nucleotides"] = cleavage_region_nuc_plus_one
        result_plus_one["cleavage_protein"] = cleavage_region_prot_plus_one

        result_minus_one["cleavage_region_found"] = True
        result_minus_one["cleavage_nucleotides"] = cleavage_region_nuc_minus_one
        result_minus_one["cleavage_protein"] = cleavage_region_prot_minus_one

        for type_name, motifs in motifs_by_type.items():
            for frames, result_dict in [
                (cleavage_region_prot, result),
                (cleavage_region_prot_plus_one, result_plus_one),
                (cleavage_region_prot_minus_one, result_minus_one),
            ]:
                for motif, label in motifs.items():
                    if CleavageSiteAnalyzer.motif_matches(motif, frames):
                        result_dict["pathogenicity"] = type_name
                        result_dict["motif_type"] = motif
                        result_dict["motif_category"] = label
                        return result, result_plus_one, result_minus_one

        return result, result_plus_one, result_minus_one


# ============================================================================
# Genotype identification
# ============================================================================


class GenotypeIdentifier:
    def __init__(self, references: dict[str, str]):
        self.references = references

    def identify(self, input_sequence: str, method="hamming", top_n=3) -> list[tuple]:
        if method == "hamming":
            similarity_func = SequenceSimilarity.hamming_distance
        else:
            similarity_func = SequenceSimilarity.pairwise_similarity

        genotype_scores: dict[str, list[float]] = {}
        genotype_best_match: dict[str, tuple[float, str]] = {}

        for header, ref_sequence in self.references.items():
            similarity = similarity_func(input_sequence, ref_sequence)

            parts = header.split("|")
            genotype = parts[2] if len(parts) >= 3 else "Unknown"

            if genotype not in genotype_scores:
                genotype_scores[genotype] = []
                genotype_best_match[genotype] = (0, header)

            genotype_scores[genotype].append(similarity)

            if similarity > genotype_best_match[genotype][0]:
                genotype_best_match[genotype] = (similarity, header)

        results = []
        for genotype, scores in genotype_scores.items():
            avg_score = sum(scores) / len(scores)
            best_score, best_header = genotype_best_match[genotype]
            results.append(
                (genotype, round(avg_score, 2), len(scores), best_header, round(best_score, 2))
            )

        results.sort(key=lambda x: x[1], reverse=True)
        return results[:top_n]


# ============================================================================
# Orchestration
# ============================================================================


def analyze_sequence(
    input_fasta: str,
    reference_sequences: dict[str, str],
    top_matches: int = 3,
    similarity_method: str = "hamming",
    pathogenicity_config: dict | None = None,
) -> dict:
    parsed_input = FASTAParser.parse_text(input_fasta)

    if not parsed_input:
        return {"error": "Could not parse input FASTA"}

    input_header = next(iter(parsed_input.keys()))
    input_sequence = parsed_input[input_header]

    identifier = GenotypeIdentifier(reference_sequences)
    genotype_matches = identifier.identify(
        input_sequence, method=similarity_method, top_n=top_matches
    )

    cfg = pathogenicity_config or {}
    pathogenicity_configured = "cleavage_start" in cfg

    if pathogenicity_configured:
        cleavage_main, cleavage_plus_one, cleavage_minus_one = CleavageSiteAnalyzer.analyze(
            input_sequence,
            cleavage_start=cfg["cleavage_start"],
            motifs_by_type=cfg.get("motifs_by_type"),
        )
    else:
        _empty = {
            "cleavage_region_found": False,
            "cleavage_nucleotides": None,
            "cleavage_protein": None,
            "pathogenicity": "Not configured",
            "motif_type": None,
            "motif_category": None,
        }
        cleavage_main = cleavage_plus_one = cleavage_minus_one = _empty.copy()

    return {
        "input_header": input_header,
        "sequence_length": len(input_sequence),
        "genotype_matches": genotype_matches,
        "pathogenicity_configured": pathogenicity_configured,
        "cleavage_main": cleavage_main,
        "cleavage_plus_one": cleavage_plus_one,
        "cleavage_minus_one": cleavage_minus_one,
        "error": None,
    }


def unpack_top_match(top_match):
    return {
        "genotype": top_match[0],
        "avg_similarity": top_match[1],
        "sample_count": top_match[2],
        "best_header": top_match[3],
        "best_score": top_match[4],
    }


def tree_to_newick(tree) -> str:
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue()


def build_tree_fasttree(aln_file):
    result = subprocess.run(
        [get_fasttree_cmd(), "-nt", "-gtr", aln_file],
        capture_output=True,
        text=True,
        check=False,
    )
    newick_str = result.stdout
    if not newick_str or not newick_str.strip():
        return None
    tree = Phylo.read(io.StringIO(newick_str), "newick")
    tree.root_at_midpoint()
    return tree


def build_tree_iqtree2(aln_file):
    subprocess.run(
        [get_iqtree_cmd(), "-s", aln_file, "-m", "MFP", "-bb", "1000", "-nt", "AUTO", "-redo"],
        capture_output=True,
        text=True,
        check=False,
    )
    treefile = aln_file + ".treefile"
    if not os.path.exists(treefile):
        return None
    tree = Phylo.read(treefile, "newick")
    tree.root_at_midpoint()
    # IQ-TREE2 expresses ultrafast bootstrap as 0-100, FastTree as 0-1 - normalize here
    for clade in tree.find_clades():
        if clade.confidence is not None:
            clade.confidence = clade.confidence / 100
    return tree


def align_sequences_mafft(input_fasta_path, output_fasta_path):
    with open(output_fasta_path, "w") as out:
        subprocess.run(
            [get_mafft_cmd(), "--localpair", input_fasta_path],
            stdout=out,
            stderr=subprocess.PIPE,
            check=False,
        )


def find_closest_neighbours(query_sequence, reference_sequences, n=20):
    similarities = []
    for header, seq in reference_sequences.items():
        score = SequenceSimilarity.hamming_distance(query_sequence, seq)
        similarities.append((header, seq, score))
    similarities.sort(key=lambda x: x[2], reverse=True)
    return similarities[:n]


def clean_sequence(seq):
    valid = set("ATCGNatcgnRYSWKMBDHVryswkmbdhv")
    return "".join(c for c in seq if c in valid)


def write_temp_fasta(query_header, query_sequence, neighbours, output_path):
    with open(output_path, "w") as f:
        f.write(f">QUERY_{query_header}\n{clean_sequence(query_sequence)}\n")
        f.writelines(f">{header}\n{clean_sequence(seq)}\n" for header, seq, score in neighbours)


def get_color(name):
    if not name:
        return "#888888"
    parts = name.split("|")
    if len(parts) < 3:
        return "#888888"
    genotype = parts[2]
    if genotype in ("?", "UNKNOWN", "") or genotype.startswith("UNCL"):
        return "#888888"
    idx = zlib.crc32(genotype.encode()) % len(PALETTE)
    return PALETTE[idx]
