# VIROnator

A tool for sensitive detection and quantitative reporting of viral presence in human samples sequenced by WGS.

This workflow works modularly to process and align sequencing datasets. It consists of sequential modules whose parameters and execution are fully controlled by the user via the configuration file (`config/ssc_config.yaml`):

1. **Standalone Utility:** Reference Sanitization & BED Generation (`rename_fasta_contigs.sh`)
2. **Module 1:** Unmapped Extraction Module (`ssc_unmapped.job`)
3. **Module 2:** Viral DB Alignment Module (`ssc_alignment.job`)
4. **Module 3:** Reporting Module (`ssc_reporting.job`)

---

## Standalone Utility: Reference Sanitization <small>(rename_fasta_contigs.sh)</small>

When adopting a new vector or plasmid reference dataset, raw FASTA header names often contain special characters (Greek letters like `α`, `β`, `γ`, trademarks `™`, `®`, symbols `+`, `%`, `#`, commas, slashes, and spaces). These special characters cause parsing failures or silent truncation in alignment tools like `bwa mem` and `samtools`.

`scripts/rename_fasta_contigs.sh` is an optional **standalone utility script** to preprocess and sanitize FASTA contig header names. **It is not integrated into the automated Snakemake workflow and should be run manually prior to reference indexing when adding new reference files.**

#### Usage Command
```bash
# 1. Standard plasmid / vector sanitization:
./scripts/rename_fasta_contigs.sh SnapGene.fa SnapGene_plasmids_modified

# 2. NCBI viral reference sanitization (using -ncbi flag to keep clean accession IDs):
./scripts/rename_fasta_contigs.sh --ncbi HumanViral_Reference_02-07-2022.fa HumanViral_Reference_02-07-2022_modified
```

---

## Workflow Execution & Modules

For detailed documentation, see [docs/ssc_extraction.md](docs/ssc_extraction.md).

### Step-by-Step Setup

### Phase 1: Prepare your Workspace
Before configuring any files, run these steps in your terminal:
1. **Clone the repository:** Clone the `VIROnator` repository to your environment.
2. **Move into the repository folder:**
   ```bash
   cd VIROnator
   ```
3. **Copy your sample list file here:** Copy the file containing the list of sample IDs you want to process directly into this cloned `VIROnator` folder.
   * **What is in this file?** It is a simple text file listing one sample ID per line. The very first line (header) of the file **must** be the word `SAMPLE`. For example:
     ```text
     SAMPLE
     sample_id_1
     sample_id_2
     ```

### Phase 2: Configure your Settings
Open `config/ssc_config.yaml` in a text editor. Configure these settings before running the pipeline:

> [!IMPORTANT]
> This pipeline is designed to run on Google Cloud Platform. Ensure that your GCS parameters (such as `data_bucket`, `output_bucket`, and project details) are properly configured to point to your GCP cloud buckets.

##### 1. Common / Joint Variables
* **`data_bucket`**: GCS source bucket containing raw CRAM files.
* **`output_bucket`**: GCS destination bucket for outputs.
* **`work_dir`**: Path to repository folder on server.
* **`samples_list`**: Name of sample ID list file.

##### 2. Reporting Module Variables
* **`reporting_module`**: Switch (`"on"` or `"off"`) to enable/disable consolidated reporting.
* **`sample_metadata_file`**: Sample metadata TSV file (`SSC_sample_metadata.tsv`).
* **`viral_rename_map_file`**: Name mapping TSV (`HumanViral_Reference_02-07-2022_modified.rename_map.tsv`).
* **`viral_bed_file`**: Target BED regions file (`HumanViral_Reference_02-07-2022_modified.bed`).

### Phase 3: Execute the Workflow

First, run Snakemake on the head node to compile job and resource configuration files:
```bash
snakemake --cores 1
```

#### Module 1: Unmapped Extraction Module (`ssc_unmapped.job`)
```bash
batchRun -multibatch samples_p2_base -config config/batch_jobexec_unmapped.config -non-spot config/ssc_unmapped.job -investigator MDJ -pau 0
```

#### Module 2: Viral DB Alignment Module (`ssc_alignment.job`)
```bash
batchRun -multibatch samples_p2_base -config config/batch_jobexec_vironator.config -non-spot config/ssc_alignment.job -investigator MDJ -pau 0
```

#### Module 3: Reporting Module (`ssc_reporting.job`)
* **Output Path:** `/mnt/disks/staff/SSC_hg38_reports/phase2/[project]/<sample_id>/<sample_id>_viral_report.tsv`
* **Batch Command:**
  ```bash
  batchRun -multibatch samples_p2_base -config config/batch_jobexec_reporting.config -non-spot config/ssc_reporting.job -investigator MDJ -pau 0
  ```

#### Module 4: Coverage Calculation Module (`ssc_coverage.job`)
* **Output Path:** `/mnt/disks/staff/SSC_hg38_coverage/phase2/[project]/<sample_id>_coverage.tsv`
* **Batch Command:**
  ```bash
  batchRun -multibatch samples_p2_base -config config/batch_jobexec_coverage.config -non-spot config/ssc_coverage.job -investigator MDJ -pau 0
  ```

---

## Merging Cloud Outputs into Master Reports

### 1. Merging All Sample Viral Reports
Once all parallel cloud batch jobs finish, combine all individual sample reports into one master file named with the phase and project:

```bash
PHASE="phase2"
PROJECT="base"  # Or sub-project name e.g. Project_CCDG_...

(
  echo -e "Sample_ID\tVirus_Accession\tVirus_Length\tVirus_Mapped_Reads\tNormalized_Coverage\tPhysical_Coverage\tHuman_Genome_Size\tSample_Read_Depth\tViral_Copy_Number\tVirus_Name_Sanitized\tSpecimen\tPhase\tProject\tSource_File";
  gsutil -m cat "gs://ml-phi-staff-m277455-p-rsa-us-central1-p-a3d4/SSC_hg38_reports/${PHASE}/${PROJECT}/**/*.tsv" | awk 'FNR==1 && NR!=1{next} {print}'
) > "/mnt/disks/staff/SSC_hg38_reports/${PHASE}/cohort_master_viral_report_${PHASE}_${PROJECT}.tsv"
```

### 2. Merging All Cloud Sample Coverages into Master Coverage File
To combine all cloud-calculated sample coverage outputs into a master file named with the phase and project:

```bash
PHASE="phase2"
PROJECT="base"  # Or sub-project name e.g. Project_CCDG_...

(
  echo -e "sample\tcoverage\tspecimen\tphase\tproject\tcram_url";
  gsutil -m cat "gs://ml-phi-staff-m277455-p-rsa-us-central1-p-a3d4/SSC_hg38_coverage/${PHASE}/${PROJECT}/**/*.tsv" | awk 'FNR==1 && NR!=1{next} {print}'
) > "/mnt/disks/staff/SSC_hg38_coverage/${PHASE}/cohort_master_coverage_${PHASE}_${PROJECT}.tsv"
```

### 3. Generating Cohort Statistical Summaries (`stats_module`)
To generate a 9-column cohort summary snapshot (`cohort_stats_summary.tsv` & `cohort_stats_summary.md`) directly in your `VIROnator` root directory:

1. Set `stats_module: "on"` and specify `master_report_file` in `config/ssc_config.yaml`:
   ```yaml
   stats_module: "on"
   master_report_file: "master_all_cohorts_viral_report_final.tsv"
   ```
2. Execute Snakemake:
   ```bash
   snakemake --cores 1
   ```

### 4. Generating Distribution TIFF Plots & Virus Statistics (`distributions_module`)
To generate publication-quality 300 DPI **TIFF distribution figures** (`plots/`) and the consolidated **`virus_stats_summary.tsv`** (`stats/`):

1. Configure Section 7 in `config/ssc_config.yaml`:
   ```yaml
   distributions_module: "on"
   master_report_file: "master_all_cohorts_viral_report_final.tsv"
   distributions_script: "generate_distributions.py"
   target_strategies:
     - "exogeneSR_viral_clean_filtered.sorted.flags.cram"
     - "exogeneSR_viral_raw_filtered.sorted.flags.cram"
   ```
2. Execute Snakemake:
   ```bash
   snakemake --cores 1
   ```
3. Outputs:
   - **`stats/virus_stats_summary.tsv`**: Per-virus prevalence, positivity ratios, and total assigned reads (un-signed percentage values).
   - **`plots/`**: 300 DPI publication-grade TIFF figures (`.tif`) with 2-panel layout, median reference lines, and metric callout badges.

#### 14-Column Master Report Schema:
1. `Sample_ID`
2. `Virus_Accession` (`None` if zero reads detected)
3. `Virus_Length` (`0` if zero reads detected)
4. `Virus_Mapped_Reads` (`0` if zero reads detected)
5. `Normalized_Coverage` (Rounded to 6 decimals: `0.000000`)
6. `Physical_Coverage` (Percentage formatted without `%` sign: `00.00`)
7. `Human_Genome_Size`
8. `Sample_Read_Depth`
9. `Viral_Copy_Number` (Rounded to 6 decimals: `0.000000`)
10. `Virus_Name_Sanitized` (`None` if zero reads detected)
11. `Specimen`
12. `Phase`
13. `Project` (or `base` if project is empty)
14. `Source_File` (CRAM source file name)
