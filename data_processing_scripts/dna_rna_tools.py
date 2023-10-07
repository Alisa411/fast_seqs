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
    return ''.join(drd.TRANSCRIBE_DICT[base] if base in drd.TRANSCRIBE_DICT else base for base in sequence)


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


def run_dna_rna_tools(*arguments):
    """
    Run various DNA/RNA sequence manipulation tools.

    Args:
        *arguments: Variable number of arguments. The first n-1 arguments should be DNA/RNA sequences,
                    and the last argument should be a string specifying the action to be performed.

    Returns:
        str or list: The result of the specified action on the input sequences.

    """

    action = arguments[-1].lower()
    sequences = arguments[:-1]
    results = []

    for sequence in sequences:
        if all(base in 'ATCGUatcgu' for base in sequence):
            if action == 'transcribe':
                result = transcribe(sequence)
            elif action == 'reverse':
                result = reverse(sequence)
            elif action == 'complement':
                result = complement(sequence)
            elif action == 'reverse_complement':
                result = reverse_complement(sequence)
        else:
            result = f"Invalid procedure: {action}"
    else:
        result = "Invalid dna/rna sequence"

        results.append(result)

        return results if len(results) > 1 else results[0]
