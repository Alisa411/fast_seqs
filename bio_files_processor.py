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
    elif not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'

    # Write the sequences in one-line format to the output file
    with open(output_fasta, "w") as file:
        for sequence_name, sequence in sequences.items():
            file.write(f"{sequence_name}\n{sequence}\n")

    return f"Successfully converted to one-line FASTA: {output_fasta}"
