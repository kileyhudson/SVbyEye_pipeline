#!/usr/bin/env python3
"""Generate synthetic assemblies used by the automated tests."""

import os
import random
from typing import Iterable, Tuple


def generate_random_sequence(length: int, gc_content: float = 0.5) -> str:
    """Construct a random DNA sequence with an approximate GC content."""

    bases_gc = ["G", "C"]
    bases_at = ["A", "T"]

    num_gc = int(length * gc_content)
    num_at = max(length - num_gc, 0)

    sequence = [random.choice(bases_gc) for _ in range(num_gc)]
    sequence.extend(random.choice(bases_at) for _ in range(num_at))
    random.shuffle(sequence)
    return "".join(sequence)


def generate_similar_sequence(template: str, similarity: float = 0.95) -> str:
    """Perturb ``template`` to achieve the requested similarity."""

    if not template:
        return template

    length = len(template)
    num_changes = int(length * (1 - similarity))

    seq_list = list(template)
    positions = random.sample(range(length), num_changes)
    for pos in positions:
        current = seq_list[pos]
        bases = ["A", "T", "G", "C"]
        bases.remove(current)
        seq_list[pos] = random.choice(bases)

    return "".join(seq_list)


def write_fasta(filename: str, sequences: Iterable[Tuple[str, str]]) -> None:
    """Write one or more sequences to ``filename`` in FASTA format."""

    with open(filename, "w", encoding="utf-8") as handle:
        for name, seq in sequences:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")


def main() -> None:
    """Create reproducible synthetic assemblies for regression testing."""

    random.seed(42)
    output_dir = "data/assemblies"
    os.makedirs(output_dir, exist_ok=True)

    print("Generating synthetic assemblies for integration tests")

    # TEST001_hap1 – contig1
    print("  Building TEST001_hap1")
    contig1_seq = generate_random_sequence(60000)
    sd1_template = contig1_seq[1000:3000]
    contig1_seq = contig1_seq[:10000] + generate_similar_sequence(sd1_template, 0.95) + contig1_seq[12000:]
    sd1_template2 = contig1_seq[20000:22500]
    contig1_seq = contig1_seq[:50000] + generate_similar_sequence(sd1_template2, 0.92) + contig1_seq[52500:]
    write_fasta(os.path.join(output_dir, "TEST001_hap1.fasta"), [("contig1", contig1_seq)])

    # TEST001_hap2 – contig2
    print("  Building TEST001_hap2")
    contig2_seq = generate_random_sequence(60000)
    template = contig2_seq[5000:7500]
    contig2_seq = contig2_seq[:30000] + generate_similar_sequence(template, 0.94) + contig2_seq[32500:]
    write_fasta(os.path.join(output_dir, "TEST001_hap2.fasta"), [("contig2", contig2_seq)])

    # TEST002_hap1 – contig3
    print("  Building TEST002_hap1")
    contig3_seq = generate_random_sequence(60000)
    template = contig3_seq[15000:18000]
    contig3_seq = contig3_seq[:45000] + generate_similar_sequence(template, 0.90) + contig3_seq[48000:]
    write_fasta(os.path.join(output_dir, "TEST002_hap1.fasta"), [("contig3", contig3_seq)])

    # TEST002_hap2 – contig4
    print("  Building TEST002_hap2")
    contig4_seq = generate_random_sequence(60000)
    template = contig4_seq[8000:11000]
    contig4_seq = contig4_seq[:25000] + generate_similar_sequence(template, 0.93) + contig4_seq[28000:]
    write_fasta(os.path.join(output_dir, "TEST002_hap2.fasta"), [("contig4", contig4_seq)])

    print("Synthetic assemblies written to data/assemblies")
    print("\nAssembly inventory:")
    for filename in sorted(os.listdir(output_dir)):
        filepath = os.path.join(output_dir, filename)
        size = os.path.getsize(filepath)
        print(f"  {filename}: {size} bytes")


if __name__ == "__main__":
    main()
