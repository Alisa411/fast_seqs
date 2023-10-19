# Import dna_rna_dict.py containing dictionaries for working with dna and rna sequences
import dna_rna_dict as drd
import os

# GC content frequency function

def gc_content(sequence: str) -> float:
    """
    Calculates the GC content in a fastq sequence.

    :param sequence: fastq sequence
    :type sequence: str
    :return: GC content as a percentage
    :rtype: float
    """
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len(sequence)
    gc_content1 = round((gc_count / total_count) * 100, 2)
    return gc_content1


# fastq sequence length function


def seq_length(sequence: str) -> int:
    """
    Calculates the length of fastq sequence.

    :param sequence: fastq sequence
    :type sequence: str
    :return: the length as a number
    :rtype: int
    """
    return len(sequence)


# the mean encoding offset for quality scores function


def mean_encoding_offset(quality_string: str) -> float:
    """
    Calculate the mean encoding offset for quality scores in the input string.

    :param quality_string: A string of quality scores
    :type quality_string: str
    :return: The mean encoding offset
    :rtype: float
    """
    total_offset = 0
    for char in quality_string:
        offset = ord(char) - 33  # Assuming 33 as the default encoding offset
        total_offset += offset
    mean_offset = total_offset / len(quality_string)
    return mean_offset


def main_fastq_tools(seqs, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0):
    """
    Process fastq sequences based on specified criteria and return filtered results.

    :param seqs: Dictionary of fastq sequences
    :type seqs: dict
    :param gc_bounds: GC content filter bounds (default is (0, 100))
    :type gc_bounds: tuple or float
    :param length_bounds: Sequence length filter bounds (default is (0, 2**32))
    :type length_bounds: tuple or float
    :param quality_threshold: Quality score threshold (default is 0)
    :type quality_threshold: int
    :return: Filtered fastq sequences
    :rtype: dict
    """
    filtered_seqs = {}

    # Checking if it is tuples or not
    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)

    for seq_name, (sequence, quality_string) in seqs.items():
        # Checking if it is a fastq sequence
        if not all(letter in drd.DNA_LETTERS for letter in sequence):
            print(f"Skipping non-fastq sequence: {seq_name}")
            continue

        # Calculate GC content
        gc = gc_content(sequence)

        # Calculate sequence length
        seq_len = seq_length(sequence)

        # Calculate mean encoding offset for quality scores
        mean_offset = mean_encoding_offset(quality_string)

        # Check if sequence meets criteria
        if gc_bounds[0] <= gc <= gc_bounds[1] and \
                length_bounds[0] <= seq_len <= length_bounds[1] and \
                mean_offset >= quality_threshold:
            filtered_seqs[seq_name] = (sequence, quality_string) # back the dictionary with filtered sequences

    return filtered_seqs


# Here should be a python script


def parse_file(filename):
    with open(filename, 'r') as f:
        lines = f.read().split('\n')

        data = {}

        i = 0
        # print(len(lines))
        # print(*lines, sep="\nx\n")
        lines = lines[:-1]
        while i < len(lines):
            line = lines[i].strip()
            # Parse 4 lines per 1 sequence
            # The first line should contain '@' sign
            if line.startswith('@'):
                sequence_id = line.split(' ')[0]  #SEQ_ID except paired-end or not and index info
                sequence = lines[i + 1].strip()  # SEQ_FASTA
                quality = lines[i + 3].strip()  # QUALITY

                data[sequence_id] = (sequence, quality)

                i += 4  # iteration on 4 lines per 1 sequence


    return data


def save_filtered_fastq(filename, filtered_data):
    """
    Save filtered FASTQ sequences to a new FASTQ file.

    Args:
        filename (str): The name of the output FASTQ file.
        filtered_data (dict): A dictionary of filtered FASTQ sequences, where keys are sequence IDs,
            and values are tuples of sequence and quality score strings.

    Returns:
        None: This function does not return a value but saves the filtered sequences to the specified FASTQ file.
    """

    with open(filename, 'w') as f:
        for sequence_id, (sequence, quality) in filtered_data.items():
            f.write(f"{sequence_id}\n{sequence}\n+\n{quality}\n")




