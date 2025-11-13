# SVbyEye Pipeline

The SVbyEye pipeline is a Snakemake workflow that turns a curated list of segmental duplication (SD) pairs and pre-computed PAF alignments into plots and an optional HTML report using the [SVbyEye](https://github.com/daewoooo/SVbyEye) R package. It is designed to sit downstream of your own SD discovery process and assumes that the alignments you want to visualise already exist.

## Repository Contents

- `Snakefile` – defines the end-to-end workflow.
- `config.yaml` – configuration template describing inputs, plotting options, and resource hints.
- `scripts/prepare_sd_table.py` – assigns stable SD identifiers and derived metadata.
- `scripts/prepare_alignment_manifest.py` – validates that every SD pair has an accompanying PAF file and records the path.
- `scripts/visualize_sd.R` – renders one PNG plot per SD using SVbyEye.
- `scripts/generate_report.py` – assembles an HTML summary that embeds the plots.
- `scripts/create_test_data.py` – produces synthetic PAF files that match the bundled example SD callset.
- `environment.yml` – conda environment specification with Snakemake, Python, R, and the libraries required by SVbyEye.
- `install_svbyeye.R` – helper script that installs SVbyEye and its R dependencies from CRAN/Bioconductor/GitHub.

## Requirements

- Linux or macOS shell environment.
- Conda or Mamba for environment management (recommended for reproducing the workflow environment).
- Ability to install Snakemake ≥7.32, Python 3.11, R 4.4, and the SVbyEye package (all provided through `environment.yml`).
- Pre-computed PAF files generated with minimap2 (or another aligner that emits a `cg` tag).

## Environment Setup

1. Clone this repository and move into the checkout.
2. Create the analysis environment (using mamba for speed if available):
   ```bash
   mamba env create -f environment.yml  # or: conda env create -f environment.yml
   conda activate svbyeye_pipeline
   ```
3. Install SVbyEye and its dependencies:
   ```bash
   Rscript install_svbyeye.R
   ```
   The script detects missing CRAN/Bioconductor packages, installs them if needed, and then retrieves SVbyEye from GitHub.

The `setup.sh` helper mirrors these steps and additionally generates `activate_pipeline.sh`. If you choose to run it, update the embedded R command to `Rscript install_svbyeye.R` first—the script still references the historical filename `install_svbyeye_github.R`.

## Input Preparation

1. **Segmental duplication table (`sd_calls`)**
   - Provide a tab-delimited file (default: `data/sd_calls.tsv`).
   - Required columns: `Sample`, `Haplotype`, `SD1_original`, `SD2_original` (formatted as `chrom:start-end`).
   - Optional columns that unlock richer summaries: `Gene_names`, `SD1_T2T`, `SD2_T2T`, `SD1_has_exon`, `SD2_has_exon`.
2. **PAF alignments**
   - Generate one minimap2 (or equivalent) PAF file per SD pair using parameters that emit the `cg` tag (e.g. `-c --eqx`).
   - Place the files in a directory of your choosing (default: `data/alignments`).
   - Update `config.yaml` → `alignments.filename_pattern` so the workflow can find each file. Supported placeholders are `{sd_id}`, `{sample}`, `{haplotype}`, and `{index}` (zero-based index padded to five digits).

For testing or demos, run `python scripts/create_test_data.py` after populating `data/sd_calls.tsv`; it will fabricate simple PAF records that match each SD pair.

## Configuring the Workflow

Edit `config.yaml` to point at your inputs and tune plotting behaviour:

- `sd_calls`: path to the SD table described above.
- `alignments`: `dir` and `filename_pattern` controlling where PAF files are discovered.
- `svbyeye`: plotting knobs consumed by `scripts/visualize_sd.R` (colour scheme, bin size, indel highlighting, figure dimensions, etc.).
- `add_gene_annotations`: when `true`, the plotting script will log a reminder if gene annotations are requested (actual overlays are not yet implemented).
- `output_dir`: directory that will receive all generated artefacts (default `results/`).
- `generate_report`: toggle HTML report creation; disable if you only want the tables and plots.
- `threads` / `memory`: optional Snakemake resource hints for cluster execution.

Run the built-in sanity check before launching the workflow:
```bash
snakemake test_config
```
This rule verifies that the configured SD table exists and echoes the resolved alignment directory/pattern.

## Running the Pipeline

1. Perform a dry run to ensure every rule resolves correctly:
   ```bash
   snakemake -n
   ```
2. Execute the workflow once satisfied with the configuration:
   ```bash
   snakemake --cores 4
   ```
   Adjust the core count to match your environment.

The default `all` target produces the SD metadata table and, when enabled, the HTML summary report. Invoke `snakemake all_plots` if you only want the individual PNG figures.

## Outputs

All deliverables land under `output_dir` (default `results/`):

- `sd_table.tsv` – enriched SD table with `SD_ID`, per-interval lengths, and the average size column.
- `sd_summary.txt` – concise text summary listing every SD identifier and T2T availability.
- `alignment_manifest.tsv` – mapping between each `SD_ID` and the validated PAF path.
- `plots/<SD_ID>.png` – one SVbyEye plot per SD (PNG format).
- `plots_complete.txt` – marker file emitted when all plots are generated.
- `summary_report.html` – optional report combining metadata and plots (only when `generate_report: true`).
- `logs/` – Snakemake log files keyed by rule.

Use `snakemake clean` to remove every generated file, or `snakemake clean_intermediate` to keep the final tables/plots while deleting the manifest and completion marker.

## Workflow Details

The Snakefile defines the following rule order:

1. `prepare_sd_table` – reads the SD callset, assigns sequential `SD_ID`s (`SD_00000`, `SD_00001`, …), calculates per-interval sizes, and writes both the enriched table and a text summary.
2. `prepare_alignment_manifest` (checkpoint) – checks that a PAF exists for each SD using the configured directory/pattern. Missing files trigger an error before plotting begins.
3. `visualize_sd` – renders PNG plots through SVbyEye. Empty PAFs produce an explicit placeholder image so downstream steps still complete.
4. `collect_plots` – waits for every plot and records the list into `plots_complete.txt`.
5. `generate_report` – builds `summary_report.html` using the SD metadata and the generated PNG files.

## Continuous Integration

The GitHub Actions workflow at `.github/workflows/test-pipeline.yml` demonstrates a full automated run: it creates the conda environment, installs SVbyEye, replaces the default SD table with bundled test data, fabricates alignments with `scripts/create_test_data.py`, and executes the workflow end-to-end. Use it as a reference for adapting the pipeline to your infrastructure.

## Getting Help

If you run into issues with the SVbyEye plotting functions themselves, consult the upstream [SVbyEye documentation](https://github.com/daewoooo/SVbyEye). For questions about this workflow (configuration, Snakemake rules, or automation), open an issue in this repository.
