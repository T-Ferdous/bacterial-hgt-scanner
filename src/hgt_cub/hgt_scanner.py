import pandas as pd
import numpy as np

def scan_hgt_islands(df, z_threshold=2.0):
    """
    Identifies candidate HGT regions based on GC content Z-score.
    """
    mean_gc = df["gc_content"].mean()
    std_gc = df["gc_content"].std()

    if std_gc == 0:
        df["z_score"] = 0.0
    else:
        df["z_score"] = (df["gc_content"] - mean_gc) / std_gc

    # Filter regions exceeding the Z-score threshold (absolute deviation)
    candidate_islands = df[df["z_score"].abs() >= z_threshold].copy()
    return candidate_islands

def export_to_bed(df, output_path):
    """
    Exports a DataFrame of HGT regions to a standard BED file format.
    """
    bed_df = df[["chrom", "start", "end", "gc_content", "z_score"]].copy()
    bed_df.to_csv(output_path, sep="\t", header=False, index=False)
