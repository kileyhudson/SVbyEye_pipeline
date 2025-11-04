#!/usr/bin/env python3
"""
Create synthetic test assembly FASTA files for testing the pipeline
"""

import random
import os

def generate_random_sequence(length, gc_content=0.5):
    """Generate a random DNA sequence with specified GC content"""
    bases_gc = ['G', 'C']
    bases_at = ['A', 'T']

    num_gc = int(length * gc_content)
    num_at = length - num_gc

    sequence = (
        [random.choice(bases_gc) for _ in range(num_gc)] +
        [random.choice(bases_at) for _ in range(num_at)]
    )
    random.shuffle(sequence)
    return ''.join(sequence)

def generate_similar_sequence(template, similarity=0.95):
    """Generate a sequence similar to template with specified similarity"""
    length = len(template)
    num_changes = int(length * (1 - similarity))

    # Start with template
    seq_list = list(template)

    # Make random changes
    positions = random.sample(range(length), num_changes)
    for pos in positions:
        current = seq_list[pos]
        # Pick a different base
        bases = ['A', 'T', 'G', 'C']
        bases.remove(current)
        seq_list[pos] = random.choice(bases)

    return ''.join(seq_list)

def write_fasta(filename, sequences):
    """Write sequences to FASTA file"""
    with open(filename, 'w') as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            # Write in 60 character lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def main():
    """Create test assembly files"""

    random.seed(42)  # For reproducibility

    os.makedirs("data/assemblies", exist_ok=True)

    print("Creating test assembly files...")

    # TEST001_hap1 - contig1
    print("  TEST001_hap1...")
    contig1_length = 60000
    contig1_seq = generate_random_sequence(contig1_length)

    # Create SD regions that are similar to each other
    # SD1: 1000-3000, SD2: 10000-12000
    sd1_template = contig1_seq[1000:3000]
    sd2_similar = generate_similar_sequence(sd1_template, similarity=0.95)
    contig1_seq = contig1_seq[:10000] + sd2_similar + contig1_seq[12000:]

    # SD1: 20000-22500, SD2: 50000-52500
    sd1_template2 = contig1_seq[20000:22500]
    sd2_similar2 = generate_similar_sequence(sd1_template2, similarity=0.92)
    contig1_seq = contig1_seq[:50000] + sd2_similar2 + contig1_seq[52500:]

    write_fasta("data/assemblies/TEST001_hap1.fasta", [
        ("contig1", contig1_seq)
    ])

    # TEST001_hap2 - contig2
    print("  TEST001_hap2...")
    contig2_length = 60000
    contig2_seq = generate_random_sequence(contig2_length)

    # SD1: 5000-7500, SD2: 30000-32500
    sd1_template = contig2_seq[5000:7500]
    sd2_similar = generate_similar_sequence(sd1_template, similarity=0.94)
    contig2_seq = contig2_seq[:30000] + sd2_similar + contig2_seq[32500:]

    write_fasta("data/assemblies/TEST001_hap2.fasta", [
        ("contig2", contig2_seq)
    ])

    # TEST002_hap1 - contig3
    print("  TEST002_hap1...")
    contig3_length = 60000
    contig3_seq = generate_random_sequence(contig3_length)

    # SD1: 15000-18000, SD2: 45000-48000
    sd1_template = contig3_seq[15000:18000]
    sd2_similar = generate_similar_sequence(sd1_template, similarity=0.90)
    contig3_seq = contig3_seq[:45000] + sd2_similar + contig3_seq[48000:]

    write_fasta("data/assemblies/TEST002_hap1.fasta", [
        ("contig3", contig3_seq)
    ])

    # TEST002_hap2 - contig4
    print("  TEST002_hap2...")
    contig4_length = 60000
    contig4_seq = generate_random_sequence(contig4_length)

    # SD1: 8000-11000, SD2: 25000-28000
    sd1_template = contig4_seq[8000:11000]
    sd2_similar = generate_similar_sequence(sd1_template, similarity=0.93)
    contig4_seq = contig4_seq[:25000] + sd2_similar + contig4_seq[28000:]

    write_fasta("data/assemblies/TEST002_hap2.fasta", [
        ("contig4", contig4_seq)
    ])

    print("Done! Created 4 test assembly files in data/assemblies/")
    print("\nAssembly files:")
    for f in os.listdir("data/assemblies"):
        filepath = os.path.join("data/assemblies", f)
        size = os.path.getsize(filepath)
        print(f"  {f}: {size} bytes")

if __name__ == "__main__":
    main()
