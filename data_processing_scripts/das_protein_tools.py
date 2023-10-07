# importing necessary modules
import protein_dict as prd
from random import choice


# Function to determine is the sequence is a protein or not
def is_protein(seq: str) -> bool:
    """
    This function checks if the sequence is a protein or not

    Arguments:
        seq (str): A sequence of aminoacids

    Output:
        returns True or False
    """
    unique_chars = set(seq)
    aminoacids = set(prd.AA_MONOISOTOPIC_MASS_DICT.keys())
    return unique_chars.issubset(aminoacids)


# Function to get pI for each aa
def get_pI(
    sequence: str,
    pI_values: dict = None,
) -> str:
    """
    Gives isoelectric point value for each aminoacid individually

    Args:
    - sequence (str): sequence for which to calculate isoelectric point
    - pI_values (dict): acid dissociation constants for each aminoacid

    Return:
    - str: string, which contains:
            - an original sequence,
            - list of tuple pairs of aminoacid and corresponding isoelectric point,
    """

    if pI_values is None:
        # Default pKa_values if not provided
        pI_values = prd.AA_pI

    aminoacid_pIs = []

    # Calculate pI for each amino acid in the sequence while preserving case
    analysed_aa = []
    for aa in sequence:
        aa_upper = aa.upper()
        if aa_upper not in analysed_aa:
            if aa_upper in pI_values:
                pI = pI_values[aa_upper]
                analysed_aa.append(aa_upper)
                if aa.isupper():
                    aminoacid_pIs.append((aa_upper, pI))
                else:
                    aminoacid_pIs.append((aa, pI))
        else:
            continue

    return f"Sequence: {sequence}. Isoelectric point of each aminoacid: {aminoacid_pIs}"


# Function to calculate frequency of unique aminoacid in the sequence
def calculate_aa_freq(sequences: str) -> dict:
    """
    Calculates the frequency of each amino acid in a protein sequence or sequences.

    :param sequences: protein sequence or sequences
    :type sequences: str or list of str
    :return: dictionary with the frequency of each amino acid
    :rtype: dict
    """

    # Creating a dictionary with aminoacid frequencies:
    amino_acid_frequency = {}

    for amino_acid in sequences:
        # If the aminoacid has been already in:
        if amino_acid in amino_acid_frequency:
            amino_acid_frequency[amino_acid] += 1
        # If the aminoacid hasn't been already in:
        else:
            amino_acid_frequency[amino_acid] = 1

    return amino_acid_frequency


# Function to convert one-letter protein sequence to three-letter protein sequence
def convert_to_3L_code(seq: str) -> str:
    """
    This function takes one letter aminoacids sequence and convert's it to three leter coding

    Arguments:
        seq (str): A sequence of aminoacids

    Output:
        same sequence but in three-letter coding
    """
    seq = seq.upper()
    sequence = "".join(prd.AA_ONE_TO_THREE_LETTER[aa] for aa in seq)
    return sequence[:-1]


# Function to calculate protein mass
def protein_mass(seq: str) -> float:
    """
    This function takes aminoacids sequence and counts it's summary molecular weight using monoisotopic masses

    Arguments:
        seq (str): A sequence of aminoacids

    Output:
        returns molecular weight
    """
    seq = seq.upper()
    mass = sum(prd.AA_MONOISOTOPIC_MASS_DICT[aa] for aa in seq)
    return mass


# Function to translate Protein to RNA
def translate_protein_rna(seq: str) -> str:
    """
    This function takes  aminoacid sequence and translates in to the RNA.
    As most of the aminoacids are coded with several different codons,
    this function will take a random codon of the set for such aminoacids.

    Arguments:
        seq (str): A sequence of RNA molecule

    Output:
        returns sequence of aminoacids
    """
    seq = seq.upper()
    rna = ""
    for aa in seq:
        codon = choice(prd.AA_CODON_DICT[aa])
        rna += codon
    return rna


def main(*args: str):
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
    action_list = {
        "get_pI": get_pI,
        "calculate_aa_freq": calculate_aa_freq,
        "translate_protein_rna": translate_protein_rna,
        "convert_to_3L_code": convert_to_3L_code,
        "protein_mass": protein_mass,
    }

    if action not in action_list:
        raise ValueError(f"No such action: {action}")

    for sequence in args[:-1]:
        if not set(sequence).issubset(prd.AA_LETTERS):
            raise ValueError(f"The sequence is not protein sequence: {sequence}")

    result = action_list[action](*args[:-1])

    return result
