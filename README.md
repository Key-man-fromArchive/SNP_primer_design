# Custom ASQ Marker Design Script

## Overview

This script is designed to create custom Allele-Specific Quantification (ASQ) markers by designing primers specific to single nucleotide polymorphisms (SNPs) in a given DNA sequence. It includes functions to introduce mismatches at specific positions in the primers to enhance specificity, generate reverse complements of sequences, design forward and reverse primers using the Primer3 library, and optimize PCR conditions for the designed primers.

## Features

- **Primer Design with SNP Specificity**: Generates primers that target specific SNP variants in a DNA sequence.
- **Mismatch Introduction**: Introduces mismatches at the antepenultimate position to increase specificity.
- **Reverse Complement Calculation**: Generates the reverse complement of sequences when necessary.
- **Primer3 Integration**: Utilizes the Primer3 library to calculate melting temperatures (Tm) and optimize primer designs.
- **PCR Condition Optimization**: Ensures primers are designed within specified melting temperature and GC content ranges.

## Requirements

- Python 3.x
- Primer3-py library

## Installation

To install the required library, use:

```bash
pip install primer3-py
```

## Usage

### Command Line

Run the script with the input DNA sequence file:

```bash
python script.py <input_file>
```

### Input File

The input file should be in FASTA format, containing the DNA sequence with SNP notation (e.g., `[A/T]`).

### Example Input

```text
>example_sequence
ATGCGT[A/T]GCTAGCTAGCTAG
```

## Functions

### introduce_mismatch(base)

Introduces a mismatch at the antepenultimate position based on predefined rules.

- **Args**: 
  - `base` (str): The original base at the antepenultimate position.
- **Returns**: 
  - `str`: The mismatched base.

### reverse_complement(seq)

Generates the reverse complement of a given DNA sequence.

- **Args**: 
  - `seq` (str): The input DNA sequence.
- **Returns**: 
  - `str`: The reverse complement of the input sequence.

### design_fw_primers(sequence, snp_pos, snp_variants, header, prefix='', fw_min_tm=55, fw_max_tm=60, dna_conc=50, na_conc=50, mg_conc=2.5, dntp_conc=0.2, max_amplicon_size=300)

Designs forward primers for the given sequence with SNP at `snp_pos` using Primer3.

- **Args**: 
  - `sequence` (str): The input DNA sequence.
  - `snp_pos` (int): The position of the SNP in the sequence (0-based index).
  - `snp_variants` (list): List of SNP variants.
  - `header` (str): The header from the FASTA file to be included in primer names.
  - Additional parameters for primer design.
- **Returns**: 
  - `dict`: A dictionary containing the designed forward primers.

### design_rev_primer(sequence, snp_pos, header, prefix='', rv_min_tm=61, rv_max_tm=68, dna_conc=50, na_conc=50, mg_conc=2.5, dntp_conc=0.2, max_amplicon_size=300)

Designs a reverse primer for the given sequence using Primer3.

- **Args**: 
  - `sequence` (str): The input DNA sequence.
  - `snp_pos` (int): The position of the SNP in the sequence (0-based index).
  - `header` (str): The header from the FASTA file to be included in primer names.
  - Additional parameters for primer design.
- **Returns**: 
  - `dict`: A dictionary containing the designed reverse primer.

### design_primers(sequence, snp_pos, snp_variants, header, is_reversed=False, fw_min_tm=55, fw_max_tm=60, rv_min_tm=61, rv_max_tm=68)

Designs both forward and reverse primers for the given sequence with SNP at `snp_pos` using Primer3.

- **Args**: 
  - `sequence` (str): The input DNA sequence.
  - `snp_pos` (int): The position of the SNP in the sequence (0-based index).
  - `snp_variants` (list): List of SNP variants.
  - `header` (str): The header from the FASTA file to be included in primer names.
  - `is_reversed` (bool): Whether the input sequence is reverse complemented.
  - Additional parameters for primer design.
- **Returns**: 
  - `dict`: A dictionary containing the designed primers.

## Output

The script outputs designed primers, including their sequences, melting temperatures, GC content, and positions. It also calculates the amplicon sizes for each forward primer.

### Example Output

```text
Designed Primers:
>AsqF_example_sequence_A_ 59.3'C
ATGCGTCCGCTAG
>Rev_primer_example_sequence_ 62.1'C
GCTAGCTAGCTAG
```

## License

This project is licensed under the MIT License.

## Contact

For any questions or issues, please contact [Kibeom park] at [ceo@invirutsech.com].

---
