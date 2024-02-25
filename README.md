This repository contains a utility called **fast_seqs** for working with fastq sequences.

This Python module provides a set of tools for working with biological sequences, including DNA, RNA, and protein sequences.

## Installation

To install fast_seqs tools you need to clone the git repository using the following command:
```bash
git clone git@github.com:Alisa411/fast_seqs.git \
cd fast_seqs
```
## Features
- Filter FASTQ: filter_fastq function allows filtering FASTQ files based on specified criteria such as GC content, sequence length, and quality score threshold.
- Biological Sequence Abstractions: Abstract base class BiologicalSequence defines common operations for biological sequences, such as DNA, RNA, or protein sequences.
- Nucleic Acid Sequences: NucleicAcidSequence class represents a nucleic acid sequence (e.g., DNA or RNA) and provides functionality for complementation and GC content calculation.
- DNA Sequences: DNASequence class extends NucleicAcidSequence and provides additional functionality for transcription into RNA.
- RNA Sequences: RNASequence class extends NucleicAcidSequence and defines the RNA-specific alphabet.
- Amino Acid Sequences: AminoAcidSequence class represents an amino acid sequence and provides functionality for translation into RNA, accounting for degeneracy in the genetic code.

## Example usage:

```python
from biological_sequence_toolbox import filter_fastq, DNASequence, RNASequence, AminoAcidSequence

# Example usage of filter_fastq function
filter_fastq("input.fastq", output_filename="filtered_output")

# Example usage of DNASequence
dna_seq = DNASequence("ATGC")
print("Complement:", dna_seq.complement())

# Example usage of RNASequence
rna_seq = RNASequence("AUGC")
print("GC Content:", rna_seq.gc_content())

# Example usage of AminoAcidSequence
amino_acid_seq = AminoAcidSequence("MPL")
rna_seq = amino_acid_seq.translate_to_rna()
print("Translated RNA Sequence:", rna_seq)
```

#### Running Protein Sequence Manipulation

To use these protein sequence manipulation tools, you can call the `main_protein_tools` function. This function accepts variable arguments, where the first n-1 arguments should be protein sequences, and the last argument should be a string specifying the action to be performed.

- Arguments:
  - `*args` (str): Variable number of arguments. The first n-1 arguments should be protein sequences, and the last argument should be a string specifying the action to be performed.

- Supported Actions:
  - "get_pI": Calculate isoelectric points for each amino acid in the sequence.
  - "calculate_aa_freq": Calculate the frequency of each amino acid in a protein sequence.
  - "translate_protein_rna": Translate an amino acid sequence to RNA.
  - "convert_to_3L_code": Convert one-letter amino acid sequence to three-letter coding.
  - "protein_mass": Calculate the molecular weight of the protein sequence.

Example usage:

```python
# Calculate the isoelectric points for each amino acid in the sequence
result = main_protein_tools("ACDE", "get_pI")
print(result)

# Calculate the frequency of each amino acid in the sequence
result = main_protein_tools("ACDE", "calculate_aa_freq")
print(result)

# Translate an amino acid sequence to RNA
result = main_protein_tools("ACDE", "translate_protein_rna")
print(result)

# Convert a one-letter amino acid sequence to three-letter coding
result = main_protein_tools("ACDE", "convert_to_3L_code")
print(result)

# Calculate the molecular weight of the protein sequence
result = main_protein_tools("ACDE", "protein_mass")
print(result)
```

## Contacts 
Alisa Fedorenko
contanct me via e-mail: afedorenko00@gmail.com
