import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def generate_density_plot(csv_path, output_image, organism_name, genome_size_mb):
    if not os.path.exists(csv_path):
        print(f"[-] Error: {csv_path} not found. Check your results directory.")
        return

    # Load the hotspot coordinates
    df = pd.read_csv(csv_path)
    
    # Convert base pairs to Megabases (Mb) for clean visualization
    df['Start_Mb'] = df['Start_BP'] / 1e6

    # Set up clean styling
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 5))

    # Plot a kernel density estimate curve combined with a histogram
    sns.histplot(data=df, x='Start_Mb', bins=80, kde=True, color='#1f77b4', edgecolor='black', alpha=0.6)
    
    # Set titles and labels suitable for a thesis defense
    plt.title(f'Chromosomal Distribution Frequency of Predicted Genomic Islands\nLineage: {organism_name}', 
              fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Chromosomal Position (Megabases)', fontsize=12)
    plt.ylabel('Detected Hotspot Frequency (Windows)', fontsize=12)
    
    # Set X-axis limit to match the actual genome size
    plt.xlim(0, genome_size_mb)
    
    plt.tight_layout()
    plt.savefig(output_image, dpi=300)
    print(f"[+] Chromosomal density plot successfully saved to: {output_image}")

if __name__ == "__main__":
    os.makedirs('results', exist_ok=True)
    
    # Pointing to the exact file verified by your ls command: PA14_HGT_hotspots.csv
    generate_density_plot(
        csv_path='results/PA14_HGT_hotspots.csv',
        output_image='results/PA14_island_density.png',
        organism_name='Pseudomonas aeruginosa PA14',
        genome_size_mb=6.53
    )
