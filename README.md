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
   - **`stats/virus_stats_summary_<PHASE>_<PROJECT>.tsv`**: Detailed per-virus summary table containing prevalence, sample counts, read counts, mean mapped reads, and preliminary classifications (`sporadic_noise`, `systematic_noise`, `virome`, `infection`).
   - **`stats/virus_stats_summary_<PHASE>_<PROJECT>_OVERALL.tsv`**: High-level executive dataset summary table reporting total positive samples, cohort size, overall positivity rate, total mapped reads, unique viruses detected, top prevalent virus, mean reads per positive sample, and overall classification.
   - **`plots/`**: 300 DPI publication-grade TIFF figures (`.tif`) with 2-panel layout and horizontal viral prevalence bar plots.

#### Distributions Module TSV Output Schemas:

##### A. Per-Virus Summary Table (`virus_stats_summary_<PHASE>_<PROJECT>.tsv`) — 13 Columns:
1. `Phase` — Cohort phase (e.g. `all_cohorts` or `phase1`)
2. `Project` — Cohort project identifier (e.g. `combined` or `base`)
3. `Strategy` — Alignment strategy label (e.g. `clean_flags`)
4. `Virus_Accession` — RefSeq accession ID of the viral species (e.g. `NC_001829`)
5. `Virus_Name_Sanitized` — Full viral species name
6. `Positive_Samples_Count` — Number of positive samples containing this virus
7. `Total_Viral_Positive_Samples` — Total number of positive samples in the cohort (with $\ge 1$ viral hit)
8. `Total_Cohort_Samples` — Total number of samples in the entire cohort ($N = 9,040$)
9. `Pct_Of_Viral_Positive_Samples` — % of virus-positive samples containing this species
10. `Pct_Of_Total_Cohort_Samples` — % of total cohort samples containing this species (Prevalence)
11. `Total_Reads_Assigned` — Total mapped reads assigned to this virus across the cohort
12. `Mean_Mapped_Reads_Per_Positive_Sample` — Average mapped reads per positive sample ($\frac{\text{Total Reads}}{\text{Positive Samples}}$)
13. `Preliminary_Classification` — Automated classification (`sporadic_noise`, `systematic_noise`, `virome`, `infection`)

##### B. Executive Dataset Overview Table (`virus_stats_summary_<PHASE>_<PROJECT>_OVERALL.tsv`) — 11 Columns:
1. `Phase` — Cohort phase (e.g. `all_cohorts`)
2. `Project` — Cohort project identifier (e.g. `combined`)
3. `Strategy` — Alignment strategy label (e.g. `clean_flags`)
4. `Total_Viral_Positive_Samples` — Total samples testing positive for any virus
5. `Total_Cohort_Samples` — Total sample count in the dataset ($N = 9,040$)
6. `Cohort_Positivity_Pct` — Overall viral positivity rate (% of total cohort)
7. `Total_Mapped_Viral_Reads` — Total mapped viral reads across all detected viruses
8. `Unique_Viruses_Detected` — Total number of distinct viral species identified
9. `Top_Prevalent_Virus` — Most prevalent viral species in the dataset
10. `Mean_Mapped_Reads_Per_Positive_Sample` — Cohort-wide average mapped reads per positive sample
11. `Preliminary_Classification` — Overall dataset-wide classification

#### 14-Column Master Report Schema:
1. `Sample_ID`
2. `Virus_Accession` (`None` if zero reads detected)
3. `Virus_Length` (`0` if zero reads detected)
4. `Virus_Mapped_Reads` (`0` if zero reads detected)
5. `Normalized_Coverage` (Rounded to 6 decimals: `0.000000`)
6. `Physical_Coverage` (Percentage formatted without `%` sign: `00.00`)
7. `Human_Genome_Size`
8. `Sample_Read_Depth`
9. `Viral_Copy_Number` (Reference-discriminating copy-number lower bound, rounded to 6 decimals: `0.000000`)
10. `Virus_Name_Sanitized` (`None` if zero reads detected)
11. `Specimen`
12. `Phase`
14. `Source_File` (CRAM source file name)

### 8. Systematic SAM Flag Comparison Module (`flag_comparison_module`)
Evaluates the clean strategy across 3 SAM flag filtering commands to isolate the impact of aligner `-f 2` flags vs manual bitwise flags:

> [!NOTE]
> This module is not estimating viral abundance. It is estimating high-confidence, reference-discriminating read-pair support. That is useful, however it will undercount viruses when related references share homologous sequence. This is for the reference-unique evidence. We report high-specificity reference-discriminating read support.

1. Enable Module 8 in `config/ssc_config.yaml`:
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

### 9. NCBI RefSeq Refinement Module (`refinement_module`)
Re-aligns candidate viral reads against complete RefSeq viral genomes for detected species groups using Bowtie2 multimapping (`-k 10`) and a Hierarchical Two-Tier Classification Cascade (`Reference_Unique` vs `Species_Supportive`):

1. Configure Section 9 in `config/ssc_config.yaml`:
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
