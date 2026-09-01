# VIROnator

A tool for sensitive detection and quantitative reporting of viral presence in human samples sequenced by WGS.

This workflow works modularly to process and align sequencing datasets. It consists of sequential modules whose parameters and execution are fully controlled by the user via the configuration file (`config/ssc_config.yaml`):

1. **Standalone Utility:** Reference Sanitization & BED Generation (`rename_fasta_contigs.sh`)
2. **Module 1:** Coverage Calculation Module (`ssc_coverage.job`)
3. **Module 2:** Unmapped Extraction Module (`ssc_unmapped.job`)
4. **Module 3:** Viral DB Alignment Module (`ssc_alignment.job`)
5. **Module 4:** Consolidated Reporting Module (`ssc_reporting.job`)
6. **Module 5:** Cohort Stats Summary Module (`generate_stats.py`)
7. **Module 6:** Distributions & Plots Module (`generate_distributions.py`)
8. **Module 7:** SAM Flag Comparison Module (`ssc_flag_comparison.job`)
9. **Module 8:** NCBI RefSeq Refinement Module (`ssc_refinement.job`)

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
* **`data_dir`**: Shared local mount directory for reading raw input files (`/mnt/disks/lab`).
* **`data_subdir`**: Middle subdirectory path inside `data_dir` (e.g. `"tertiary/SSC_hg38/WGS"`, or `""` for flat datasets).
* **`input_keyword`**: Optional middle keyword following sample ID dot (e.g. `"final"`, `"sorted"`, or `""` for none).
* **`input_suffix`**: Target file extension following keyword (e.g. `"cram"`, `"bam"`).
* **`output_dir`**: Shared local mount directory for output files (`/mnt/disks/staff`).
* **`work_dir`**: Active repository root directory on head node (default `"."`).
* **`data_bucket`**: Source GCS bucket containing raw sequence files.
* **`output_bucket`**: Destination GCS bucket for output data and cluster execution logs.
* **`gcp_account`**: GCP Service Account email address for cluster execution.
* **`pi`**: Principal Investigator initials/tag (e.g. `"MDJ"`).
* **`pau`**: Project Allocation Unit code (e.g. `0`).
* **`task_name`**: Cluster job log identifier prefix tag.
* **`dataset`**: Dataset cohort prefix label (e.g. `"SSC"`).
* **`phase`**: Dataset phase identifier (e.g. `"2"`, or `"none"` / `""` for flat datasets).
* **`project`**: Sub-project folder identifier (e.g. `"base"`, or `"none"` / `""` for flat datasets).
* **`samples_list`**: Name of text file listing target sample IDs (header must be `SAMPLE`).
* **`ref_dir`**: Directory path where reference genomes and database files are stored (`/mnt/disks/staff/refs`).
* **`scripts_dir`**: Repository folder path containing helper scripts (`/mnt/disks/staff/scripts`).

##### 2. Module Execution Switches
To execute specific pipeline stages, turn `"on"` or `"off"` the switches for the features and modules you want to execute in `config/ssc_config.yaml` (e.g. `coverage_module: "on"`, `unmapped_extraction: "on"`, `viral_db_alignment: "on"`, `reporting_module: "on"`, `stats_module: "on"`, `distributions_module: "on"`, `flag_comparison_module: "on"`, or `refinement_module: "on"`).

### Phase 3: Execute the Workflow

First, run Snakemake on the head node to compile job and resource configuration files:
```bash
snakemake --cores 1
```

#### Module 1: Coverage Calculation Module (`ssc_coverage.job`)
```bash
batchRun -multibatch samples_list.txt -config config/batch_jobexec_coverage.config -non-spot config/ssc_coverage.job -investigator MDJ -pau 0
```
Combine all cloud-calculated sample coverage outputs into a master file named as `${DATASET}_master_coverage.tsv`. This file will be used for `sample_metadata_file` in Module 4 - Reporting Module.

#### Module 2: Unmapped Extraction Module (`ssc_unmapped.job`)
```bash
batchRun -multibatch samples_list.txt -config config/batch_jobexec_unmapped.config -non-spot config/ssc_unmapped.job -investigator MDJ -pau 0
```

#### Module 3: Viral DB Alignment Module (`ssc_alignment.job`)
```bash
batchRun -multibatch samples_list.txt -config config/batch_jobexec_vironator.config -non-spot config/ssc_alignment.job -investigator MDJ -pau 0
```

#### Module 4: Consolidated Reporting Module (`ssc_reporting.job`)
```bash
batchRun -multibatch samples_list.txt -config config/batch_jobexec_reporting.config -non-spot config/ssc_reporting.job -investigator MDJ -pau 0
```
- **Generated Per-Sample File**: `<sample_id>_viral_report.tsv`
- **Combined Master Output File**: Once all parallel cloud batch jobs finish, combine all individual sample reports into one master file named with the phase and project: `master_all_cohorts_viral_report_final.tsv` (used as input for Module 5 Stats).

#### Module 5: Cohort Stats Summary Module (`stats_module`)
```bash
snakemake cohort_stats_summary.tsv --cores 1
```
- **Generated Summary Files**:
  - `cohort_stats_summary.tsv` (Tab-delimited cohort statistics snapshot)
  - `cohort_stats_summary.md` (Markdown-formatted summary report)
  - `<stats_out_dirname>/virus_stats_summary_<PHASE>_<PROJECT>.tsv` (Detailed per-virus summary table containing prevalence, sample counts, read counts, mean mapped reads, and preliminary classifications: `sporadic_noise`, `systematic_noise`, `virome`, `infection`).
  - `<stats_out_dirname>/virus_stats_summary_<PHASE>_<PROJECT>_OVERALL.tsv` (High-level executive dataset summary table reporting total positive samples, cohort size, overall positivity rate, total mapped reads, unique viruses detected, top prevalent virus, mean reads per positive sample, and overall classification).

#### Module 6: Distributions & Plots Module (`distributions_module`)
```bash
snakemake generate_distributions --cores 1
```
- **Generated Summary Files & Figures**:
  - `<stats_out_dirname>/virus_stats_summary_<PHASE>_<PROJECT>.tsv` (13-column per-virus summary table)
  - `<stats_out_dirname>/virus_stats_summary_<PHASE>_<PROJECT>_OVERALL.tsv` (11-column cohort overview table)
  - `<plots_out_dirname>/<VIRUS_NAME>.tif` (2-panel 300 DPI publication TIFF figures)

#### Module 7: SAM Flag Comparison Module (`ssc_flag_comparison.job`)
```bash
batchRun -multibatch samples_all_cohorts -config config/batch_jobexec_flag_comparison.config -non-spot config/ssc_flag_comparison.job -investigator MDJ -pau 0
```

#### Module 8: NCBI RefSeq Refinement Module (`ssc_refinement.job`)
```bash
# 1. Build taxonomy index TSV:
snakemake config/db_metadata/viral_reference_taxonomy_index.tsv --cores 1

# 2. Download RefSeq complete genomes & build date-stamped Bowtie2 combined reference:
snakemake config/ncbi_download.completed --cores 1

# 3. Compile refinement jobexec cluster files:
snakemake config/ssc_refinement.job --cores 1

# 4. Submit parallel batch run across target refinement samples:
batchRun -multibatch config/refinement_samples.tsv -config config/batch_jobexec_refinement.config -non-spot config/ssc_refinement.job -investigator MDJ -pau 0

# 5. Merge all cohort outputs into master report (run on head node after cluster jobs finish):
snakemake cohort_refinement_master.tsv --cores 1
```

### 7. Systematic SAM Flag Comparison Module (`flag_comparison_module`)
Evaluates the clean strategy across 3 SAM flag filtering commands to isolate the impact of aligner `-f 2` flags vs manual bitwise flags:

> [!NOTE]
> This module is not estimating viral abundance. It is estimating high-confidence, reference-discriminating read-pair support. That is useful, however it will undercount viruses when related references share homologous sequence. This is for the reference-unique evidence. We report high-specificity reference-discriminating read support.

1. Enable Module 7 in `config/ssc_config.yaml`:
   ```yaml
   flag_comparison_module: "on"
   cohort_scope: "combined_all"
   ```
2. Compile job and resource config files on head node:
   ```bash
   snakemake --cores 1
   ```
3. Submit parallel batch run across all cohort samples (single batch across all phases & projects):
   ```bash
   batchRun -multibatch samples_all_cohorts -config config/batch_jobexec_flag_comparison.config -non-spot config/ssc_flag_comparison.job -investigator MDJ -pau 0
   ```
4. Merge all cohort outputs into single master TSV report:
   ```bash
   snakemake cohort_flag_comparison_master.tsv --cores 1
   ```

#### 7-Column Flag Comparison Master TSV Schema (`cohort_flag_comparison_master.tsv`):
1. `phase` — Cohort phase (e.g. `phase1`)
2. `project` — Cohort project identifier (e.g. `base`)
3. `sample` — Sample ID
4. `strategy` — Strategy evaluated (`clean`)
5. `original_command_flags_reads` — Read pair count using original command WITH `-f 2`
6. `original_command_flags_modified_reads` — Read pair count using original command WITHOUT `-f 2`
7. `manual_flags_reads` — Read pair count using manual `gawk` bitwise flag verification

### 8. NCBI RefSeq Refinement Module (`refinement_module`)
Re-aligns candidate viral reads against complete RefSeq viral genomes for detected species groups using Bowtie2 multimapping (`-k 10`) and a Hierarchical Two-Tier Classification Cascade (`Reference_Unique` vs `Species_Supportive`):

1. Configure Module 8 in `config/ssc_config.yaml`:
   ```yaml
   refinement_module: "on"
   make_taxonomy_index_stage: "on"
   ncbi_download_stage: "on"
   build_combined_ref_stage: "on"
   refinement_alignment_stage: "on"
   ```
2. Build 14-column taxonomy index TSV:
   ```bash
   snakemake config/db_metadata/viral_reference_taxonomy_index.tsv --cores 1
   ```
3. Download complete RefSeq genomes & build Bowtie2 combined index:
   ```bash
   snakemake config/ncbi_download.completed --cores 1
   ```
4. Compile jobexec refinement cluster files:
   ```bash
   snakemake config/ssc_refinement.job --cores 1
   ```
5. Submit parallel batch run across target refinement samples:
   ```bash
   batchRun -multibatch config/refinement_samples.tsv -config config/batch_jobexec_refinement.config -non-spot config/ssc_refinement.job -investigator MDJ -pau 0
   ```

> [!IMPORTANT]
> **Post-Cluster Consolidation**: `batchRun` executes sample jobs in parallel across cluster nodes, writing per-sample TSVs under `/mnt/disks/staff/SSC_hg38_refinement/`. Once all parallel cluster jobs have finished, run Step 6 separately on the head node to aggregate all results into `cohort_refinement_master.tsv`:

6. Consolidate cohort output into master report (run on head node after cluster jobs finish):
   ```bash
   snakemake cohort_refinement_master.tsv --cores 1
   ```

#### 11-Column Refinement Master TSV Schema (`cohort_refinement_master.tsv`):
1. `Sample_ID` — Sample accession identifier
2. `Phase` — Cohort phase tag (e.g. `phase1`)
3. `Project` — Cohort project tag (or `base` if empty)
4. `Virus_Accession` — Best-supported RefSeq Accession ID
5. `Species_TaxID` — NCBI Species Taxonomy ID
6. `Species_Name` — Taxonomic Species Scientific Name
7. `Reference_Unique_Reads` — Read pairs uniquely supporting this specific reference
8. `Species_Supportive_Reads` — Read pairs supporting this species group (group-collapsed)
9. `Total_Refined_Reads` — Total read pairs (`Reference_Unique_Reads` + `Species_Supportive_Reads`)
10. `Refined_Copy_Number` — Reference-discriminating copy-number lower bound
11. `Classification_Status` — Call status (`CONFIRMED_UNIQUE` or `SPECIES_COLLAPSED`)
