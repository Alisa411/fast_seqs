from data_processing_scripts.dna_rna_tools import main_dna_rna_tools
from data_processing_scripts.das_protein_tools import main_protein_tools
from data_processing_scripts.fastq_tools import main_fastq_tools


# Call the main_dna_rna_tools function with the necessary arguments.
result = main_dna_rna_tools("ATcg", "reverse")
print(result)

# Call the main_protein_tools function with the necessary arguments.
result = main_protein_tools("ACDE", "protein_mass")
print(result)

# Call the main_fastq_tools function with the necessary arguments.
EXAMPLE_FASTQ = {
    '@SRX079804:1:SRR292678:1:1101:21885:21885': (
        'ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
        'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
}
filtered_sequences = main_fastq_tools(
    seqs=EXAMPLE_FASTQ,
    gc_bounds=(0, 80),  # GC content от 40% до 60%
    length_bounds=(10, 100),  # Длина последовательности от 10 до 100
    quality_threshold=0  # Порог качества
)

print(filtered_sequences)
