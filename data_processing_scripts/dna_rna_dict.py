# DNA set of letters
DNA_LETTERS = set("ATGCatgc")


# RNA set of letters
RNA_LETTERS = set("AUGCaugc")


# Dictionary: keys - nucleotide letter of dna; values - complement nucleotide letter of rna
TRANSCRIBE_DICT = {
    'T': 'U',
    't': 'u'
}


# Dictionary: keys - nucleotide letter of dna; values - complement nucleotide letter of dna
COMPLEMENT_DICT = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'a': 't',
    't': 'a',
    'c': 'g',
    'g': 'c'}
