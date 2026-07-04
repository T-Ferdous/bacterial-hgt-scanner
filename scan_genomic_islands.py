import os

def calculate_gc_content(sequence):
    if not sequence:
        return 0
    g_count = sequence.count('G') + sequence.count('g')
    c_count = sequence.count('C') + sequence.count('c')
    return (g_count + c_count) / len(sequence) * 100

def scan_genome(fasta_path, window_size=5000, step_size=2500):
    print(f"=== Scanning Genome: {fasta_path} ===")
    
    # Read the FASTA sequence
    header = ""
    sequence_fragments = []
    
    if not os.path.exists(fasta_path):
        print(f"ERROR: File {fasta_path} not found.")
        return
        
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if not header:
                    header = line.strip()
            else:
                sequence_fragments.append(line.strip())
                
    full_sequence = "".join(sequence_fragments)
    genome_length = len(full_sequence)
    
    # Calculate overall host baseline GC
    baseline_gc = calculate_gc_content(full_sequence)
    print(f"Genome Length: {genome_length} bp")
    print(f"Host Baseline GC Content: {baseline_gc:.2f}%")
    print("-" * 50)
    print(f"{'Window Coordinates':<25} | {'Window GC%':<12} | {'Divergence Status'}")
    print("-" * 50)
    
    # Sliding window scanning
    island_count = 0
    for i in range(0, genome_length - window_size, step_size):
        window = full_sequence[i:i+window_size]
        window_gc = calculate_gc_content(window)
        
        # Flag significant downward GC deviations (typical for mobile genomic islands)
        if window_gc < (baseline_gc - 6.0):
            start_coord = i
            end_coord = i + window_size
            print(f"{start_coord:,} - {end_coord:,} bp | {window_gc:.2f}%     | [!] Potential Low-GC Island Detected")
            island_count += 1
            
    print("-" * 50)
    print(f"Scan complete. Found {island_count} suspicious low-GC horizontal integration hotspots.")

# Run the scanner on your PA14 reference dataset
# (Adjust path if your local filename varies slightly)
scan_genome("extra_strains/bac1.fna") 
