import os
import matplotlib.pyplot as plt
import seaborn as sns
from src.hgt_cub.cub_filter import sliding_window_cub
from src.hgt_cub.hgt_scanner import scan_hgt_islands

def plot_genome_hgt(fasta_path, output_img="docs/genome_gc_hgt_track.png"):
    print("Extracting sliding window GC content...")
    df = sliding_window_cub(fasta_path, window_size=5000, step_size=1000)
    islands = scan_hgt_islands(df, z_threshold=2.0)

    df["start_mb"] = df["start"] / 1e6
    mean_gc = df["gc_content"].mean()

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(12, 5), dpi=300)

    ax.plot(df["start_mb"], df["gc_content"], color="#1f77b4", linewidth=1.2, label="GC Content (%)")
    ax.axhline(mean_gc, color="black", linestyle="--", linewidth=1, label=f"Mean GC ({mean_gc:.1f}%)")

    if not islands.empty:
        for idx, row in islands.iterrows():
            ax.axvspan(row["start"] / 1e6, row["end"] / 1e6, color="#e74c3c", alpha=0.35)
        ax.axvspan(0, 0, color="#e74c3c", alpha=0.35, label=f"HGT Candidate (|Z| >= 2.0, n={len(islands)})")

    ax.set_title("Bacterial Genome GC-Content Profile & Predicted HGT Islands", fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel("Genomic Position (Mb)", fontsize=11)
    ax.set_ylabel("GC Content (%)", fontsize=11)
    ax.set_ylim(df["gc_content"].min() - 2, df["gc_content"].max() + 2)
    ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.9)

    plt.tight_layout()

    os.makedirs(os.path.dirname(output_img), exist_ok=True)
    plt.savefig(output_img, bbox_inches="tight")
    plt.close()
    print(f"Plot saved successfully to {output_img}")

if __name__ == "__main__":
    fasta = "data/raw/E_coli_CFT073.fasta"
    if os.path.exists(fasta):
        plot_genome_hgt(fasta)
    else:
        print(f"Error: {fasta} not found. Please run main.py first.")
