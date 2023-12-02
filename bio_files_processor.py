import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> str:
    """
    Convert a multi-line FASTA file into a one-line FASTA file format.

    Sometimes we deal with broken fasta files with sequences split across multiple lines.
    This function concatenates these lines into one-line sequences.

    Args:
        input_fasta (str): Input multi-line FASTA file.
        output_fasta (str, optional): Output one-line FASTA file name.
                                     If not provided, a default name will be generated.

    Returns:
        str: A message indicating the status of the operation, including the name of the output file.
    """

    sequences = {}
    sequence_name = ""
    sequence = ""

    with open(input_fasta, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if sequence_name:
                    sequences[sequence_name] = sequence
                sequence_name = line
                sequence = ""
            else:
                sequence += line

        if sequence_name:
            sequences[sequence_name] = sequence

    # Determine the output one-line FASTA file path
    if output_fasta is None:
        input_filename = os.path.splitext(os.path.basename(input_fasta))[0]
        output_fasta = f"{input_filename}_oneline.fasta"
    # Write down format to output without initial format
    elif not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'

    # Write the sequences in one-line format to the output file
    with open(output_fasta, "w") as file:
        for sequence_name, sequence in sequences.items():
            file.write(f"{sequence_name}\n{sequence}\n")

    return f"Successfully converted to one-line FASTA: {output_fasta}"


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: str = None) -> str:
    """
    Shifts the starting position of a sequence within a one-line FASTA file.

    Args:
        input_fasta (str): Input one-line FASTA file.
        shift (int): The amount to shift the starting position of the sequence.
                     A positive shift moves the sequence to the right, while a negative shift moves it to the left.
        output_fasta (str): Output one-line FASTA file name.

    Returns:
        str: A message indicating the status of the operation, including the name of the output file.
    """
    sequences = ""
    sequence_name = ""
    sequence = ""

    with open(input_fasta, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_name:
                    sequences += f"{sequence_name}\n{sequence}\n"
                sequence_name = line
                sequence = ""
            else:
                sequence += line

        if sequence_name:
            sequences += f"{sequence_name}\n{sequence}\n"

    if not sequences or '>' not in sequences:
        return "Error: No valid sequence found in the input file"

    header_end = sequences.find('\n') if '\n' in sequences else len(sequences)
    header = sequences[:header_end]

    sequence = sequences[header_end:].replace('\n', '')

    if len(sequence) == 0:
        return "Error: Empty sequence found in the input file"

    shifted_sequence = sequence[shift % len(sequence):] + sequence[:shift % len(sequence)]

    if output_fasta is None:
        output_fasta = f"{os.path.splitext(input_fasta)[0]}_shifted.fasta"
    elif not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'

    with open(output_fasta, "w") as outfile:
        outfile.write(f"{header}\n{shifted_sequence}\n")

    return f"Successfully shifted start position by {shift} to create: {output_fasta}"

