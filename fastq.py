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
    gc_content1 = round((gc_count / total_count) * 100, 2)
    return gc_content1

#fastq sequence length function


def seq_length(sequence: str) -> int:
    """
    Calculates the length of fastq sequence.

    :param sequence: fastq sequence
    :type sequence: str
    :return: the length as a number
    :rtype: int
    """
    return len(sequence)


#the mean encoding offset for quality scores function


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


fastq_letters = set("ATGCatgc")


def main(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
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

    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)


    for seq_name, (sequence, quality_string) in seqs.items():
        # Checking if it is a fastq sequence
        if not all(letter in fastq_letters for letter in sequence):
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
            filtered_seqs[seq_name] = sequence

    print(filtered_seqs)

    return filtered_seqs


# Example usage:
# seqs = {
#     # 'name' : ('sequence', 'quality')
#     '@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
#     '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
#     '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
#     '@SRX079804:1:SRR292678:1:1101:47176:47176': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'),
#     '@SRX079804:1:SRR292678:1:1101:149302:149302': ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA', '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A'),
#     '@SRX079804:1:SRR292678:1:1101:170868:170868': ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT', 'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB'),
#     '@SRX079804:1:SRR292678:1:1101:171075:171075': ('CATTATAGTAATACGGAAGATGACTTGCTGTTATCATTACAGCTCCATCGCATGAATAATTCTCTAATATAGTTGTCAT', 'HGHHHHGFHHHHFHHEHHHHFGEHFGFGGGHHEEGHHEEHBHHFGDDECEGGGEFGF<FGGIIGEBGDFFFGFFGGFGF'),
#     '@SRX079804:1:SRR292678:1:1101:175500:175500': ('GACGCCGTGGCTGCACTATTTGAGGCACCTGTCCTCGAAGGGAAGTTCATCTCGACGCGTGTCACTATGACATGAATG', 'GGGGGFFCFEEEFFDGFBGGGA5DG@5DDCBDDE=GFADDFF5BE49<<<BDD?CE<A<8:59;@C.C9CECBAC=DE'),
#     '@SRX079804:1:SRR292678:1:1101:190136:190136': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT', 'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE'),
#     '@SRX079804:1:SRR292678:1:1101:190845:190845': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC', 'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE'),
#     '@SRX079804:1:SRR292678:1:1101:198993:198993': ('AGTTATTTATGCATCATTCTCATGTATGAGCCAACAAGATAGTACAAGTTTTATTGCTATGAGTTCAGTACAACA', '<<<=;@B??@<>@><48876EADEG6B<A@*;398@.=BB<7:>.BB@.?+98204<:<>@?A=@EFEFFFEEFB'),
#     '@SRX079804:1:SRR292678:1:1101:204480:204480': ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG', '<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')
#     }
#
# main(seqs, gc_bounds=80, length_bounds=100, quality_threshold=10)



