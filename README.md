This repository contains a utility called **fast_seqs** for working with fastq sequences.

_fast_seqs is a bionformatic tool for beggining your steps in bioinformatical analysis. It works with nucleic acid sequences, protein and fastq sequences. Read the **Usage** section to get to know abilities of this tool.

## Installation

To install fast_seqs tools you need to clone the git repository using the following command:
```bash
git clone git@github.com:Alisa411/fast_seqs.git \
cd fast_seqs
```

## Usage

### Working with Nucleic Acids

#### Nucleic Acid Manipulation

The script allows you to perform various operations on nucleic acid sequences. It includes the following functions:

- **Transcribe**
  - Description: Transcribes a DNA sequence to an RNA sequence.
  - Arguments:
    - `sequence` (str): A DNA sequence.
  - Returns:
    - `str`: The transcribed RNA sequence.

- **Reverse**
  - Description: Reverses a given sequence.
  - Arguments:
    - `sequence` (str): The input sequence.
  - Returns:
    - `str`: The reversed sequence.

- **Complement**
  - Description: Finds the complement of a DNA or RNA sequence.
  - Arguments:
    - `sequence` (str): A DNA or RNA sequence.
  - Returns:
    - `str`: The complement sequence.

- **Reverse Complement**
  - Description: Finds the reverse complement of a DNA or RNA sequence.
  - Arguments:
    - `sequence` (str): A DNA or RNA sequence.
  - Returns:
    - `str`: The reverse complement sequence.

#### Running Sequence Manipulation

To use these nucleic acid manipulation tools, you can call the `main_dna_rna_tools` function. This function accepts variable arguments, where the first n-1 arguments should be DNA/RNA sequences, and the last argument should be a string specifying the action to be performed.

- Arguments:
  - `*args` (str): Variable number of arguments. The first n-1 arguments should be DNA/RNA sequences, and the last argument should be a string specifying the action to be performed.

- Supported Actions:
  - "transcribe": Transcribe a DNA sequence to an RNA sequence.
  - "reverse": Reverse a sequence.
  - "complement": Find the complement of a DNA or RNA sequence.
  - "reverse_complement": Find the reverse complement of a DNA or RNA sequence.

- Returns:
  - `str` or `list`: The result of the specified action on the input sequences.

Example usage:

```python
# Transcribe a DNA sequence to RNA
result = main_dna_rna_tools("ATCG", "transcribe")
print(result)

# Reverse a sequence
result = main_dna_rna_tools("ATCG", "reverse")
print(result)

# Find the complement of a sequence
result = main_dna_rna_tools("ATCG", "complement")
print(result)

# Find the reverse complement of a sequence
result = main_dna_rna_tools("ATCG", "reverse_complement")
print(result)
```


### Working with Protein Sequences

#### Protein Sequence Manipulation

The script allows you to perform various operations on protein sequences. It includes the following functions:

- **Get pI (Isoelectric Point) for Each Amino Acid**
  - Description: Calculates the isoelectric point value for each amino acid individually.
  - Arguments:
    - `sequence` (str): The protein sequence for which to calculate isoelectric points.
    - `pI_values` (dict): Acid dissociation constants for each amino acid (optional).
  - Returns:
    - `str`: A string containing:
      - The original sequence.
      - A list of tuple pairs of amino acid and corresponding isoelectric point.

- **Calculate Amino Acid Frequency**
  - Description: Calculates the frequency of each amino acid in a protein sequence.
  - Arguments:
    - `sequences` (str): A protein sequence or sequences.
  - Returns:
    - `dict`: A dictionary with the frequency of each amino acid.

- **Convert One-Letter Protein Sequence to Three-Letter Coding**
  - Description: Converts a one-letter amino acid sequence to a three-letter amino acid coding.
  - Arguments:
    - `seq` (str): A sequence of amino acids.
  - Returns:
    - `str`: The same sequence in three-letter coding.

- **Calculate Protein Mass**
  - Description: Calculates the molecular weight of a protein sequence using monoisotopic masses.
  - Arguments:
    - `seq` (str): A sequence of amino acids.
  - Returns:
    - `float`: The molecular weight.

- **Translate Protein to RNA**
  - Description: Translates an amino acid sequence to RNA, using random codons for amino acids with multiple codons.
  - Arguments:
    - `seq` (str): A sequence of amino acids.
  - Returns:
    - `str`: The RNA sequence.

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


### Working with Fastq Sequences

#### Fastq Sequence Analysis

The script allows you to analyze Fastq sequences based on specified criteria. It includes the following functions:

- **Calculate GC Content**
  - Description: Calculates the GC content in a Fastq sequence as a percentage.
  - Arguments:
    - `sequence` (str): The Fastq sequence.
  - Returns:
    - `float`: GC content as a percentage.

- **Calculate Sequence Length**
  - Description: Calculates the length of a Fastq sequence.
  - Arguments:
    - `sequence` (str): The Fastq sequence.
  - Returns:
    - `int`: The length of the sequence as a number.

- **Calculate Mean Encoding Offset for Quality Scores**
  - Description: Calculates the mean encoding offset for quality scores in the input string.
  - Arguments:
    - `quality_string` (str): A string of quality scores.
  - Returns:
    - `float`: The mean encoding offset.

#### Processing Fastq Sequences

To process Fastq sequences based on specified criteria and obtain filtered results, you can call the `main_fastq_tools` function. This function accepts a dictionary of Fastq sequences and several optional parameters for filtering.

- Arguments:
  - `seqs` (dict): A dictionary of Fastq sequences.
  - `gc_bounds` (tuple or float): GC content filter bounds (default is (0, 100)).
  - `length_bounds` (tuple or float): Sequence length filter bounds (default is (0, 2^32)).
  - `quality_threshold` (int): Quality score threshold (default is 0).

- Returns:
  - `dict`: Filtered Fastq sequences as a dictionary.

Example usage:

```python
# Example Fastq sequences
EXAMPLE_FASTQ = {
    '@SRX079804:1:SRR292678:1:1101:21885:21885': (
        'ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
        'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD')
}

# Process Fastq sequences with filtering criteria
filtered_sequences = main_fastq_tools(
    seqs=EXAMPLE_FASTQ,
    gc_bounds=(40, 60),  # GC content from 40% to 60%
    length_bounds=(10, 100),  # Sequence length from 10 to 100
    quality_threshold=30,  # Quality score threshold
)
print(filtered_sequences)
```

## Troubleshooting

1. **Invalid Sequence Format:**
   - **Issue:** If you encounter errors related to invalid sequence formats, ensure that your input sequences match the expected format (e.g., DNA, RNA, or protein sequences) for the respective functions.
   - **Solution:** Check your input sequences and ensure they are correctly formatted for the function you are using.

2. **Unsupported Action:**
   - **Issue:** If you receive a "No such action" error, it means that you provided an unsupported action when using the main functions for DNA/RNA or protein sequences.
   - **Solution:** Review the supported actions listed in the README and make sure you specify one of these actions when calling the main functions.

3. **Sequence Length Errors:**
   - **Issue:** If you encounter issues related to sequence length, such as sequences being too short or too long, ensure that your input sequences meet the specified length criteria.
   - **Solution:** Adjust the length bounds (if applicable) or provide sequences that meet the required length criteria.

4. **Quality Score Errors (Fastq Sequences):**
   - **Issue:** When working with Fastq sequences, if you face problems related to quality scores, check the quality scores in your Fastq sequences and verify that they are correctly encoded.
   - **Solution:** Ensure that the quality scores are encoded following the expected encoding offset (default is 33) and that they match the sequence length.

5. **Undefined Variables:**
   - **Issue:** If you encounter errors related to undefined variables (e.g., missing imports), ensure that you have imported the required modules and dictionaries (e.g., `dna_rna_dict` for DNA/RNA scripts and `protein_dict` for protein scripts).
   - **Solution:** Confirm that you have correctly imported the necessary modules and dictionaries at the beginning of your script.

6. **Function Parameter Errors:**
   - **Issue:** If you encounter errors related to function parameters, double-check that you are providing the correct arguments and data types to the functions.
   - **Solution:** Review the function signatures and documentation to ensure you are passing the expected arguments and data types.

7. **Mixed U and T in DNA Sequences:**
   - **Issue:** DNA sequences should not contain both 'U' and 'T'. If you encounter this error, check your DNA sequence data.
   - **Solution:** Ensure that your DNA sequences contain either 'U' (uracil) for RNA or 'T' (thymine) for DNA, but not both.

8. **Dictionary Errors:**
   - **Issue:** If you encounter errors related to missing dictionary values (e.g., amino acid pI values), make sure that the required dictionaries (e.g., `AA_pI`, `AA_MONOISOTOPIC_MASS_DICT`) are correctly defined and imported.
   - **Solution:** Verify the contents of your dictionaries and ensure they match the expected format.

By addressing these common issues, you can troubleshoot and resolve potential problems when using the provided scripts for DNA/RNA, protein, or Fastq sequence analysis.

## Contacts 
Alisa Fedorenko
contanct me via e-mail: afedorenko00@gmail.com
