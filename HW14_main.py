from abc import ABC, abstractmethod

class BiologicalSequence(ABC):
    """
        Abstract base class representing a biological sequence.

        Defines common operations for biological sequences, such as DNA, RNA, or protein sequences.

        Attributes:
            sequence (str): The sequence of characters representing the biological sequence.

        Methods:
            __len__(): Returns the length of the biological sequence.
            __getitem__(index): Gets the character at the specified index or returns a subsequence.
            __str__(): Returns a string representation of the biological sequence.
            is_valid_alphabet(): Checks if all characters in the sequence belong to a valid alphabet.
        """
    def __init__(self, sequence):
        self.sequence = sequence

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

    @abstractmethod
    def is_valid_alphabet(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
        Represents a nucleic acid sequence, such as DNA or RNA.

        Inherits from BiologicalSequence and adds functionality specific to nucleic acid sequences.

        Attributes:
            sequence (str): The sequence of characters representing the nucleic acid sequence.

        Methods:
            complement(): Returns the complement of the nucleic acid sequence.
            gc_content(as_percentage=False): Calculates the GC content of the nucleic acid sequence.

        Overrides:
            is_valid_alphabet(): Checks if all characters in the sequence belong to a valid nucleic acid alphabet.
        """

    COMPLEMENT_DICT = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c'}

class NucleicAcidSequence():
    def

class DNASequence(NucleicAcidSequence):
    DNA_LETTERS = set("ATGCatgc")

class RNASequence(NucleicAcidSequence):
    RNA_LETTERS = set("AUGCaugc")