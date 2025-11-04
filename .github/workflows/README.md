# GitHub Actions Workflows

This directory contains CI/CD workflows for testing the SVbyEye pipeline.

## Workflows

### 1. `test-pipeline.yml` - Full Pipeline Test

**Triggers:**
- Push to `main` or `claude/*` branches
- Pull requests to `main`
- Manual trigger (workflow_dispatch)

**What it does:**
1. Sets up Python, R, minimap2, samtools
2. Installs SVbyEye from GitHub
3. Tests each script individually:
   - SD filtering
   - Sequence extraction
   - Minimap2 alignment
   - SVbyEye visualization
4. Runs the full Snakemake pipeline
5. Uploads results as artifacts

**Artifacts generated:**
- `pipeline-results`: Full pipeline outputs (plots, reports, logs)
- `test-plots`: Individual test visualizations

**Runtime:** ~10-15 minutes

### 2. `quick-test.yml` - Quick Validation

**Triggers:**
- Push to `claude/*` branches
- Pull requests

**What it does:**
1. Validates YAML syntax
2. Checks Python script syntax
3. Checks R script syntax
4. Lints Snakefile
5. Verifies test data exists

**Runtime:** ~2-3 minutes

## Manual Trigger

To manually run the full pipeline test:

1. Go to the **Actions** tab on GitHub
2. Select **Test SVbyEye Pipeline**
3. Click **Run workflow**
4. Select your branch
5. Click **Run workflow**

## Viewing Results

After a workflow runs:

1. Go to the **Actions** tab
2. Click on the workflow run
3. Scroll down to **Artifacts**
4. Download `pipeline-results` or `test-plots`
5. Extract and view the plots!

## Troubleshooting

If the workflow fails:

1. Click on the failed workflow run
2. Expand the failed step to see error messages
3. Check the logs for details

Common issues:
- **SVbyEye installation fails**: Check if the GitHub repo is accessible
- **Sequence extraction fails**: Check test data integrity
- **Snakemake errors**: Check Snakefile syntax
