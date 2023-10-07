

def transcribe(sequence):
    return ''.join(TRANSCRIBE_DICT[base] if base in TRANSCRIBE_DICT else base for base in sequence)


def reverse(sequence):
    return sequence[::-1]


def complement(sequence):
    return ''.join(COMPLEMENT_DICT[base] if base in COMPLEMENT_DICT else base for base in sequence)


def reverse_complement(sequence):
    return reverse(complement(sequence))


def run_dna_rna_tools(*argument):
    action = argument[-1].lower()
    sequences = argument[:-1]
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
