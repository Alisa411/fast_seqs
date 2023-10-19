# Import dna_rna_dict.py containing dictionaries for working with dna and rna sequences
import dna_rna_dict as drd


def transcribe(sequence):
    """
    Transcribe a DNA sequence to an RNA sequence.

    Args:
        sequence (str): A DNA sequence.

    Returns:
        str: The transcribed RNA sequence.
    """
    if set(sequence).issubset(drd.DNA_LETTERS):
        return ''.join(drd.TRANSCRIBE_DICT[base] if base in drd.TRANSCRIBE_DICT else base for base in sequence)
    else:
        raise ValueError("Invalid DNA sequence for transcription")


def reverse(sequence):
    """
    Reverse a sequence.

    Args:
        sequence (str): The input sequence.

    Returns:
        str: The reversed sequence.
    """
    return sequence[::-1]


def complement(sequence):
    """
    Find the complement of a DNA or RNA sequence.

    Args:
        sequence (str): A DNA or RNA sequence.

    Returns:
        str: The complement sequence.
    """
    return ''.join(drd.COMPLEMENT_DICT[base] if base in drd.COMPLEMENT_DICT else base for base in sequence)


def reverse_complement(sequence):
    """
    Find the reverse complement of a DNA or RNA sequence.

    Args:
        sequence (str): A DNA or RNA sequence.

    Returns:
        str: The reverse complement sequence.
    """
    return reverse(complement(sequence))
