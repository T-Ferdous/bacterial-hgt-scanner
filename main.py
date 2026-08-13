import os
import urllib.request
from src.hgt_cub.cub_filter import sliding_window_cub
from src.hgt_cub.hgt_scanner import scan_hgt_islands, export_to_bed

def main():
    # File paths
    fasta_path = "data/raw/E_coli_CFT073.fasta"
    bed_path = "data/processed/predicted_hgt_islands.bed"

    # Ensure directories exist
    os.makedirs("data/raw", exist_ok=True)
    os.makedirs("data/processed", exist_ok=True)

    # Download E. coli CFT073 genome if not present
    if not os.path.exists(fasta_path):
        print("Downloading E. coli CFT073 genome from NCBI...")
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_004431.1&rettype=fasta&retmode=text"
        urllib.request.urlretrieve(url, fasta_path)
        print("Download complete.")

    print("Running 5 kb sliding window CUB filter...")
    df_windows = sliding_window_cub(fasta_path, window_size=5000, step_size=1000)
    print(f"Processed {len(df_windows)} windows.")

    print("Scanning for HGT candidate islands (|Z| >= 2.0)...")
    df_islands = scan_hgt_islands(df_windows, z_threshold=2.0)
    print(f"Detected {len(df_islands)} candidate HGT windows.")

    print(f"Exporting results to {bed_path}...")
    export_to_bed(df_islands, bed_path)
    print("Pipeline execution complete!")

if __name__ == "__main__":
    main()
