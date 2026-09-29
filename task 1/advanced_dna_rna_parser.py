"""
CS6120 

CHARUL - 002535330

Assignment 1 - FASTA DNA/RNA Parser

Parse a FASTA-like file and classify each record as DNA, RNA or Invalid.
"""

import os
import re
import sys
import unittest
from collections import Counter


VALID = set("ACGTUNRYSWKMBDHV")   # IUPAC nucleotide codes
HEADER = re.compile(r"^>")        # header lines start with '>'
WS = re.compile(r"\s+")           # whitespace to strip from sequence lines
 
CODONS = {
    "DNA": {"start": {"ATG"}, "stop": {"TAA", "TAG", "TGA"}, "core": "ACGT"},
    "RNA": {"start": {"AUG"}, "stop": {"UAA", "UAG", "UGA"}, "core": "ACGU"},
}


# Part 1: Parsing + classification

def parse_fasta(text):
    """Return (header, merged_sequence) for each record.
    A record's sequence is every non-header line until the next '>' or EOF,
    with all whitespace/newlines removed and joined into one string.
    """
    header, parts = None, []
    for line in text.splitlines():
        if HEADER.match(line):
            if header is not None:
                yield header, "".join(parts)
            header, parts = line.rstrip(), []
        elif header is not None:
            parts.append(WS.sub("", line))
    if header is not None:
        yield header, "".join(parts)


def classify(seq):
    """DNA = has T, no U. RNA = has U, no T. Otherwise Invalid.
 
    Invalid also covers: empty sequence, any character outside the IUPAC set,
    or both T and U present. A sequence with only A/C/G (+ambiguous) and no
    T/U is ambiguous by definition; defaulted to DNA since ACGT is the DNA
    alphabet.
    """
    s = seq.upper()
    if not s or any(c not in VALID for c in s):
        return "Invalid"
    has_t, has_u = "T" in s, "U" in s
    if has_t and has_u:
        return "Invalid"
    if has_u:
        return "RNA"
    return "DNA"


# Part 2: Find open reading frames (ORFs) in each valid DNA/RNA sequence.

def find_orfs(seq, kind):
    """Return list of (start_index, stop_index, subsequence).
 
    Iterative approach:
      1. Every position where the start codon occurs is a candidate start.
      2. From each start, step forward by 3 (in-frame) looking for a stop.
      3. First in-frame stop closes the ORF. A codon holding an ambiguous
         letter simply won't equal any start/stop, so it is skipped as
         non-matching and scanning continues (per the ATN clarification).
    """
    starts = CODONS[kind]["start"]
    stops = CODONS[kind]["stop"]
    n = len(seq)
    orfs = []
    for i in range(n - 2):
        if seq[i:i + 3] in starts:
            j = i + 3
            while j + 3 <= n:
                if seq[j:j + 3] in stops:
                    orfs.append((i, j, seq[i:j + 3]))
                    break
                j += 3
    return orfs


# Part 3: Nucleotide frequencies + summary stats

def counts(seq, kind):
    """Return (ordered core-count dict, ambiguous count)."""
    c = Counter(seq)
    core = {base: c.get(base, 0) for base in CODONS[kind]["core"]}
    return core, len(seq) - sum(core.values())


# Drivers

def report_part1(records):
    print("Part 1: Parsing & Classification")
    for header, seq, kind in records:
        print(header)
        print(seq)
        print(f"Type: {kind}\n")
 
 
def report_part2(records):
    print("Part 2: Open Reading Frames")
    for header, seq, kind in records:
        print(header)
        print(f"Type: {kind}")
        if kind == "Invalid":
            print("Skipped (Invalid sequence)\n")
            continue
        orfs = find_orfs(seq, kind)
        print(f"Found {len(orfs)} ORF(s):")
        for start, stop, sub in orfs:
            print(f"Start at {start}, Stop at {stop}, ORF: {sub}")
        print()
 
 
def report_part3(records):
    print("Part 3: Statistics")
    dna_lengths, rna_lengths, invalid = [], [], 0
    for header, seq, kind in records:
        print(header)
        print(f"Type: {kind}")
        print(f"Length: {len(seq)}")
        if kind == "Invalid":
            invalid += 1
            print("Skipped (Invalid sequence)\n")
            continue
        core, ambiguous = counts(seq, kind)
        pairs = ", ".join(f"{b}={n}" for b, n in core.items())
        print(f"Counts: {pairs}, ambiguous={ambiguous}")
        print(f"Number of ORFs: {len(find_orfs(seq, kind))}\n")
        (dna_lengths if kind == "DNA" else rna_lengths).append(len(seq))
 
    mean = lambda xs: round(sum(xs) / len(xs), 2) if xs else "N/A"
    
    print(f"Valid DNA sequences: {len(dna_lengths)}")
    print(f"Mean DNA length: {mean(dna_lengths)}")
    print(f"Valid RNA sequences: {len(rna_lengths)}")
    print(f"Mean RNA length: {mean(rna_lengths)}")
    print(f"Invalid sequences: {invalid}")
 
 
REPORTS = {"1": report_part1, "2": report_part2, "3": report_part3}

def load(path):
    """Read a FASTA file into a list of (header, seq, type)."""
    with open(path, encoding="utf-8") as f:
        return [(h, s, classify(s)) for h, s in parse_fasta(f.read())]
 
 
def run_file(path, part):
    if not os.path.exists(path):
        print(f"[missing] {path} not found; skipping.\n")
        return
    records = load(path)
    if part == "all":
        for i, key in enumerate(("1", "2", "3")):
            if i:
                print()
            REPORTS[key](records)
    else:
        REPORTS[part](records)

 
# unit tests
 
class Tests(unittest.TestCase):
    def test_parse_merges_multiline_and_strips_whitespace(self):
        text = ">a\nAC GT\nAC\n>b\nTTTT\n"
        self.assertEqual(list(parse_fasta(text)),
                         [(">a", "ACGTAC"), (">b", "TTTT")])
 
    def test_classify(self):
        self.assertEqual(classify("ACGTN"), "DNA")
        self.assertEqual(classify("ACGUN"), "RNA")
        self.assertEqual(classify("ACGTU"), "Invalid")   # both T and U
        self.assertEqual(classify("ACGTX"), "Invalid")   # bad char
        self.assertEqual(classify(""), "Invalid")        # empty
        self.assertEqual(classify("ACGN"), "DNA")        # no T/U -> DNA default
 
    def test_find_orfs_basic_and_multiple(self):
        self.assertEqual(find_orfs("ATGACGTAA", "DNA"),
                         [(0, 6, "ATGACGTAA")])
        self.assertEqual(len(find_orfs("ATGAAATAGATGCCCTGA", "DNA")), 2)
 
    def test_find_orfs_ambiguous_codon_is_skipped_not_matched(self):
        # ambiguous middle codon (ANG) is non-matching, scan continues to TAA
        self.assertEqual(find_orfs("ATGANGTAA", "DNA"),
                         [(0, 6, "ATGANGTAA")])
        # ambiguous start (ATN) never matches -> no ORF
        self.assertEqual(find_orfs("ATNCCCTAA", "DNA"), [])
        # no in-frame stop -> no ORF
        self.assertEqual(find_orfs("ATGCCCCCC", "DNA"), [])
 
    def test_find_orfs_rna(self):
        self.assertEqual(find_orfs("AUGCCCUAA", "RNA"),
                         [(0, 6, "AUGCCCUAA")])
 
    def test_counts_core_and_ambiguous(self):
        core, amb = counts("ATGACGTAANNRY", "DNA")
        self.assertEqual(core, {"A": 4, "C": 1, "G": 2, "T": 2})
        self.assertEqual(amb, 4)
 
 
def run_tests():
    suite = unittest.TestLoader().loadTestsFromTestCase(Tests)
    unittest.TextTestRunner(verbosity=2).run(suite)
 
 
# main
 
def main(argv):
    if len(argv) >= 2 and argv[1] == "test":
        run_tests()
        return
 
    if len(argv) == 1:
        # Default: each provided test file with its matching part.
        defaults = [
            ("data/test_part1_in.fasta", "1"),
            ("data/test_part2_in.fasta", "2"),
            ("data/test_part3_in.fasta", "3"),
        ]
        for i, (path, part) in enumerate(defaults):
            if i:
                print("\n" + "=" * 60 + "\n")
            print(f"# File: {path}")
            run_file(path, part)
        return
 
    path = argv[1]
    part = argv[2] if len(argv) > 2 else "all"
    if part not in REPORTS and part != "all":
        sys.exit(f"Unknown part '{part}'. Use 1, 2, 3, or all.")
    run_file(path, part)
 
 
if __name__ == "__main__":
    main(sys.argv)