import pytest
import pandas as pd
from src.hgt_cub.cub_filter import calculate_gc, sliding_window_cub
from src.hgt_cub.hgt_scanner import scan_hgt_islands, export_to_bed

def test_calculate_gc():
    assert calculate_gc("ATGC") == 50.0
    assert calculate_gc("GGCC") == 100.0
    assert calculate_gc("AATT") == 0.0

def test_sliding_window_cub(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">seq1\n" + "ATGC" * 3000 + "\n")
    
    df = sliding_window_cub(str(fasta_file), window_size=5000, step_size=1000)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 8
    assert "gc_content" in df.columns

def test_scan_hgt_islands():
    # 9 normal windows around 50% GC, 1 strong outlier at 80% GC
    gc_values = [50.0, 50.1, 49.9, 50.2, 49.8, 50.0, 50.1, 49.9, 50.0, 80.0]
    data = {
        "chrom": ["seq1"] * len(gc_values),
        "start": [i * 1000 for i in range(len(gc_values))],
        "end": [(i * 1000) + 5000 for i in range(len(gc_values))],
        "gc_content": gc_values
    }
    df = pd.DataFrame(data)
    islands = scan_hgt_islands(df, z_threshold=2.0)
    
    assert isinstance(islands, pd.DataFrame)
    assert len(islands) >= 1
    assert 80.0 in islands["gc_content"].values

def test_export_to_bed(tmp_path):
    data = {
        "chrom": ["seq1"],
        "start": [1000],
        "end": [6000],
        "gc_content": [80.0],
        "z_score": [2.5]
    }
    df = pd.DataFrame(data)
    out_file = tmp_path / "out.bed"
    export_to_bed(df, str(out_file))
    assert out_file.exists()
