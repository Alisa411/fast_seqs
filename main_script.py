from data_processing_scripts.dna_rna_tools import transcribe, reverse, complement, reverse_complement
from data_processing_scripts.das_protein_tools import get_pI, calculate_aa_freq, translate_protein_rna, convert_to_3L_code, protein_mass
from data_processing_scripts.fastq_script import main_fastq_tools, parse_file, save_filtered_fastq
import data_processing_scripts.dna_rna_dict as drd
import data_processing_scripts.protein_dict as prd
import os


def main_dna_rna_tools(*args: str):
    """
    Run various DNA/RNA sequence manipulation tools.

    Args:
        *args: Variable number of arguments. The first n-1 arguments should be DNA/RNA sequences,
                    and the last argument should be a string specifying the action to be performed.

    Returns:
        str or list: The result of the specified action on the input sequences.

    """

    action = args[-1].lower()
    sequences = args[:-1]
    results = []

    action_list = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    for sequence in sequences:
        if all(base in drd.DNA_LETTERS or base in drd.RNA_LETTERS for base in sequence):
            if 'U' in sequence and 'T' in sequence:
                raise ValueError("Invalid sequence: Contains both U and T")
            elif action == 'transcribe' and not set(sequence).issubset(drd.DNA_LETTERS):
                raise ValueError("Invalid DNA sequence for transcription")
            result = action_list[action](sequence)
        else:
            raise ValueError(f"Invalid sequence for procedure: {action}")
        results.append(result)

    return results if len(results) > 1 else results[0]


def main_protein_tools(*args: str):
    """
    Main function to perform various actions on protein sequences.

    Args:
    - *args: Variable number of arguments. The first n-1 arguments should be protein sequences,
             and the last argument should be a string specifying the action to be performed.

    Returns:
    - The result of the specified action on the input protein sequences.

    Raises:
    - ValueError: If the specified action is not supported or if there is an error in the number of sequences.
                  Also raised if the input sequences are not valid protein sequences.

    Supported Actions:
    - "get_pI": Calculate isoelectric points for each amino acid in the sequence.
    - "calculate_aa_freq": Calculate the frequency of each amino acid in a protein sequence.
    - "translate_protein_rna": Translate amino acid sequence to RNA, using random codons for each amino acid.
    - "convert_to_3L_code": Convert one-letter amino acid sequence to three-letter coding.
    - "protein_mass": Calculate the molecular weight of the protein sequence.
    """

    action = args[-1]
    sequences = args[:-1]
    action_list = {
        "get_pI": get_pI,
        "calculate_aa_freq": calculate_aa_freq,
        "translate_protein_rna": translate_protein_rna,
        "convert_to_3L_code": convert_to_3L_code,
        "protein_mass": protein_mass,
    }

    if action not in action_list:
        raise ValueError(f"No such action: {action}")

    for sequence in sequences:
        if not set(sequence).issubset(prd.AA_LETTERS):
            raise ValueError(f"The sequence is not protein sequence: {sequence}")

    result = action_list[action](*sequences)

    return result

def main(input_path: str, output_filename: str = None):
    """
    Process a FASTQ file, apply filtering, and save the filtered sequences to a new FASTQ file.

    Args:
        input_path (str): The path to the input FASTQ file.
        output_filename (str, optional): The name for the output FASTQ file. If not provided, the input file name is used.

    Returns:
        None: This function does not return a value but saves the filtered sequences to a new FASTQ file.
    """


    parsed_data = parse_file(input_path)
    filtered_seqs = main_fastq_tools(parsed_data)
    if output_filename is None:
        output_filename = input_path.split("/")[-1]
    else:
        output_filename = output_filename + ".fastq"

    os.makedirs("./fastq_filtrator_resuls", exist_ok=True)
    output_path = "./fastq_filtrator_resuls/" + output_filename
    save_filtered_fastq(output_path, filtered_seqs)