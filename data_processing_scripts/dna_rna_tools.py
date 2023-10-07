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
