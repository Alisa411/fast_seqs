#GC content frequency function

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
    gc_content = round((gc_count / total_count) * 100, 2)
    return gc_content

#fastq sequence length function

def seq_length(sequence: str) -> int:
    """
    Calculates the length of fastq sequence.

    :param sequence: fastq sequence
    :type sequence: str
    :return: the length as a number
    :rtype: int
    """
    total_count = len(sequence)
    return total_count

#the mean encoding offset for quality scores function
def mean_encoding_offset(quality_string):
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





