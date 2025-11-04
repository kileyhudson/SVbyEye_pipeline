# SVbyEye Pipeline for Segmental Duplications

## What This Pipeline Does

This pipeline visualizes segmental duplication (SD) pairs using SVbyEye. It:

1. Takes your SD calls (TSV file with coordinates)
2. Filters to interesting SDs (those with genes, or top N by size)
3. Extracts sequences for each SD pair
4. Aligns SD1 vs SD2 using minimap2
5. Creates beautiful visualizations showing sequence similarity/differences

## Quick Start

### 1. Setup your input files

Put your files in the right places:
- SD calls: `data/sd_calls.tsv` (your table with SD coordinates)
- Assembly FASTA: `data/assemblies/` (your haplotype assemblies)
- Reference genome: `data/reference/` (T2T reference, optional)

### 2. Edit the config file

Open `config.yaml` and set:
- Path to your SD calls
- Path to your assembly files
- How many SDs to visualize (start with 20!)
- Filtering criteria

### 3. Run the pipeline

```bash
# Dry run to see what will happen (doesn't actually run anything)
snakemake -n

# Run with 4 cores
snakemake --cores 4

# Run specific parts
snakemake --cores 4 filter_sds  # Just filter SDs
snakemake --cores 4 all_plots   # Make all visualizations
```

### 4. Check your results

- Filtered SDs: `results/filtered_sds.tsv`
- Alignments: `results/alignments/`
- Plots: `results/plots/`
- Summary report: `results/summary_report.html`

## Understanding the Output

Each SD pair gets:
- A PAF alignment file (shows how sequences align)
- A PNG/PDF visualization showing:
  - **Blue lines** = forward strand matches
  - **Red/orange lines** = reverse strand matches
  - **Breaks/gaps** = insertions/deletions between the SDs

## Common Issues

**"I don't have assembly FASTAs, just BAM files"**
- You'll need to extract sequences from BAM or get FASTA assemblies

**"4115 SDs is too many!"**
- Start with `max_sds: 20` in config
- Filter by `priority: genes` to get interesting ones first

**"I don't know what coordinates to use"**
- If you have T2T coordinates, use those (easier to share)
- If only contig coordinates, use those (need assembly FASTA)

**"Snakemake is confusing"**
- Think of it like a recipe book
- Each "rule" is a recipe step
- Snakemake reads the recipes and follows them in order

## Learning Snakemake

Key concepts:
- **rule**: A step in the pipeline (like "filter SDs" or "make plot")
- **input**: What files this step needs
- **output**: What files this step creates
- **shell**: The command to run
- **wildcard**: A placeholder (like {sample} or {sd_id}) that gets filled in

Example rule:
```python
rule align_sequences:
    input:
        seq1 = "sequences/{sd_id}_SD1.fa",
        seq2 = "sequences/{sd_id}_SD2.fa"
    output:
        paf = "alignments/{sd_id}.paf"
    shell:
        "minimap2 -x asm20 -c --eqx {input.seq1} {input.seq2} > {output.paf}"
```

This says: "To make an alignment PAF, take two sequences and run minimap2"

## File Structure

```
SVbyEye_pipeline/
├── Snakefile              # Main workflow (the recipe book)
├── config.yaml            # Settings (what you edit)
├── scripts/
│   ├── filter_sds.py      # Filter SD calls
│   ├── extract_sequences.py  # Get sequences from assemblies
│   └── visualize_sd.R     # Make SVbyEye plots
├── data/                  # Your input files go here
│   ├── sd_calls.tsv
│   └── assemblies/
└── results/               # Output appears here
    ├── filtered_sds.tsv
    ├── sequences/
    ├── alignments/
    └── plots/
```

## Next Steps

1. Start small (20 SDs)
2. Check the plots make sense
3. Adjust filters if needed
4. Scale up to more SDs
5. Add gene annotations to plots (optional)

Don't worry! The pipeline does the hard work. You just need to set up the config file and run it.
