# Building a Parser for DNA/RNA Sequences

## Name: . CHARUL

## NUID: 002535330

Parses FASTA-like files, classifies each sequence as DNA / RNA / Invalid,
locates open reading frames (ORFs) and reports nucleotide frequencies plus
run-level summary statistics.

## How to run

Place the three test files in a `data/` folder next to the script:

```
data/test_part1_in.fasta
data/test_part2_in.fasta
data/test_part3_in.fasta
```

Then:

```
python advanced_dna_rna_parser.py
```

With no arguments the script runs each test file with its matching part:
Part 1 on `test_part1_in.fasta`, Part 2 on `test_part2_in.fasta`, Part 3 on `test_part3_in.fasta`.

Other modes:

```
python advanced_dna_rna_parser.py <file>          # run all 3 parts on one file
python advanced_dna_rna_parser.py <file> 1        # Part 1 only
python advanced_dna_rna_parser.py <file> 2        # Part 2 only
python advanced_dna_rna_parser.py <file> 3        # Part 3 only
python advanced_dna_rna_parser.py test            # run the unit tests
```

No third-party dependencies; standard library only (Python 3).

## Assumptions

### Nucleotide validity and classification

- Valid characters are the IUPAC nucleotide codes: `A C G T U N R Y S W K M B D H V`.
- DNA = contains `T` and no `U`; RNA = contains `U` and no `T`.
- Invalid = an empty sequence, any character outside the IUPAC set or a sequence containing both `T` and `U`.
- A sequence with only `A/C/G` (plus ambiguous letters) and no `T` or `U` is
  ambiguous by definition; it defaults to DNA, since ACGT is the DNA alphabet.

### Ambiguous letters

- Ambiguous letters (`N`, `R`, `Y`, ...) are valid for classification but are
  never treated as a concrete base. In ORF detection they are handled for
  free by exact codon matching: a codon such as `ANG` or `ATN` cannot equal
  any start or stop codon, so it is treated as non-matching. A mid-ORF
  ambiguous codon is skipped and scanning continues; an ambiguous start
  simply is not recognized as a start.
- In Part 3 counts, ambiguous letters are excluded from the core `A/C/G/T`
  (or `A/C/G/U`) counts and reported together as a single `ambiguous` figure.

### ORF logic

- The iterative approach is used: every position where the start codon occurs
  is a candidate start (starts are not restricted to global frame 1) and from
  each start the sequence is read in steps of 3 until the first in-frame stop
  codon, which closes the ORF.
- Start/stop indices are the first-base index of each codon; the reported ORF
  subsequence includes the stop codon.
- Length/mean statistics use only valid sequences; invalid sequences are
  counted separately and excluded from the DNA/RNA length means. A mean over
  zero sequences is reported as `N/A` rather than dividing by zero.

## Regex usage

Regex is used in the parsing layer, kept deliberately small so the logic stays
readable:

- `^>` detects header lines (records begin with `>`).
- `\s+` strips all whitespace/newlines from sequence lines before they are
  joined, which is what merges a multiline record into one contiguous string.

Codon and base validity are handled with exact set membership rather than
regex, because start/stop codons and the IUPAC alphabet are small fixed sets;
set lookups are clearer and avoid subtle regex edge cases with ambiguous
letters. A pure-regex ORF alternative would be, for DNA,
`ATG(?:[ACGT]{3})*?(?:TAA|TAG|TGA)`, but the iterative scan makes the in-frame
stepping and ambiguous-codon handling explicit.

## Unit tests

`python advanced_dna_rna_parser.py test` runs an embedded `unittest` suite
covering multiline parsing/whitespace stripping, DNA/RNA/Invalid
classification (including the both-`T`-and-`U` and no-`T`/`U` cases), ORF
detection (basic, multiple, RNA, ambiguous-skip, ambiguous-start, no-stop),
and nucleotide counting.

## Files

- `advanced_dna_rna_parser.py` - the single submission script (all three parts + tests)
- `README.md` - this file
- `sample_output.txt` - captured console output
- `data/` - input FASTA files
