import pandas as pd
from Bio import SeqIO

def calculate_gc(sequence):
    """Calculates GC content percentage for a given DNA sequence."""
    sequence = sequence.upper()
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    total_bases = len(sequence)
    if total_bases == 0:
        return 0.0
    return ((g_count + c_count) / total_bases) * 100.0

def sliding_window_cub(fasta_path, window_size=5000, step_size=1000):
    """
    Scans a FASTA genome using a sliding window.
    Returns a DataFrame with window start, end, and GC content.
    """
    record = next(SeqIO.parse(fasta_path, "fasta"))
    sequence = str(record.seq)
    seq_len = len(sequence)

    windows = []
    for start in range(0, seq_len - window_size + 1, step_size):
        end = start + window_size
        subseq = sequence[start:end]
        gc = calculate_gc(subseq)
        windows.append({
            "chrom": record.id,
            "start": start,
            "end": end,
            "gc_content": gc
        })

    return pd.DataFrame(windows)
