# SSC Dataset Extraction & Alignment Workflow

This tool works modularly to process and align sequencing datasets. It consists of sequential modules whose parameters and execution are fully controlled by the user via the configuration file (`config/ssc_config.yaml`):

1. **Unmapped Extraction Module:** Extracts unmapped human reads from CRAM files. This preprocesses the reads for subsequent viral analysis.
2. **Viral DB Alignment Module (VIROnator):** Maps the preprocessed unmapped reads to the viral database, filters out decoy/RNA sequences, and runs final read filtering based on targeted BED regions.
3. **Reporting Module:** Consolidates metrics across output CRAM files (`clean`, `clean_flags`, `raw`, `raw_flags`), evaluates physical coverage, merges sample metadata, and produces unified 14-column TSV reports.

## Prerequisites
- Snakemake
- Cluster modules: `samtools`, `bwa`, `python`, `jobexec/2.0.1`

## Storage Recommendations: GCS vs. Local Disk

1. **References Directory:** All reference genomes and database files (human, viral, decoy, rna) must be stored directly in the same flat directory on the GCS mount (specifically `/mnt/disks/staff/refs/`) without any nested subfolders.
2. **Access:** Because the batch jobs run on separate ephemeral cloud VM instances, they must be stored on the GCS bucket (`/mnt/disks/staff/refs/`) so the batch execution nodes can access them.

---

## Standalone Utility: FASTA Header Cleaner <small>(rename_fasta_contigs.sh)</small>

When adopting a new vector or plasmid reference dataset, raw FASTA header names often contain special characters (Greek letters like `α`, `β`, `γ`, trademarks `™`, `®`, symbols `+`, `%`, `#`, commas, slashes, and spaces). These special characters cause parsing failures or silent truncation in alignment tools like `bwa mem` and `samtools`.

`scripts/rename_fasta_contigs.sh` is an optional **standalone utility script** to preprocess and sanitize FASTA contig header names. **It is not integrated into the automated Snakemake workflow and should be run manually prior to reference indexing when adding new reference files.**

### How Headers Are Sanitized
1. **Whitespace Replacement:** All spaces, tabs, and non-breaking whitespace characters are immediately converted into underscores `_`.
2. **Greek Transliteration:** Converts Greek characters to English words (`α` $\rightarrow$ `alpha`, `β` $\rightarrow$ `beta`, `λ` $\rightarrow$ `lambda`, `μ` $\rightarrow$ `mu`, `φ` $\rightarrow$ `phi`).
3. **Trademark Removal:** Strips `™`, `®`, `©`, `(TM)`, `(R)`.
4. **Symbol Translation:** Converts mathematical/strand symbols (`(+)` $\rightarrow$ `_plus_`, `(-)` $\rightarrow$ `_minus_`, `+` $\rightarrow$ `_plus_`, `&` $\rightarrow$ `_and_`).
5. **Special Character Conversion:** Any remaining non-alphanumeric character is converted to `_`, and consecutive delimiters are collapsed.
6. **Collision Deduplication:** If two distinct raw headers clean down to the exact same string, numeric suffixes (`_2`, `_3`) are automatically appended to guarantee that **every contig ID is unique**.

7. **NCBI Accession Mode (`-n` / `--ncbi`):** When processing NCBI-formatted FASTA headers (e.g. `>NC_001829 Adeno-associated virus - 4, complete genome.`), passing `--ncbi` strips everything after the first whitespace to keep the accession ID (`NC_001829`) as the contig ID while preserving the original raw header string in `rename_map.tsv`.

### Standalone Execution Syntax

Command options (`-p` / `--prefix` and `-n` / `--ncbi`) **must be placed before the input FASTA file**:

```bash
./scripts/rename_fasta_contigs.sh [OPTIONS] <INPUT_FASTA> <OUTPUT_PREFIX>
```

#### Available Options
* `-p <PREFIX>` / `--prefix <PREFIX>`: Prepends a prefix string to all contig IDs (e.g., `-p "mm39_"` for mouse contigs to prevent collisions with human contigs).
* `-n` / `--ncbi`: Activates NCBI mode, stripping descriptions after the first whitespace to leave accession IDs (e.g., `NC_001829`), while recording full raw headers in `<OUTPUT_PREFIX>.rename_map.tsv`.

#### Usage Examples

```bash
# 1. Standard plasmid / vector sanitization:
./scripts/rename_fasta_contigs.sh SnapGene.fa SnapGene_plasmids_modified

# 2. Mouse reference sanitization (prepends "mm39_" prefix):
./scripts/rename_fasta_contigs.sh -p "mm39_" mm39_raw.fa mm39_ms_modified

# 3. NCBI viral reference sanitization (strips descriptions after accession ID):
./scripts/rename_fasta_contigs.sh --ncbi HumanViral_Reference_02-07-2022.fa HumanViral_Reference_02-07-2022_modified

# 4. Combining prefix and NCBI options together:
./scripts/rename_fasta_contigs.sh -p "mm39_" --ncbi input_ncbi_raw.fa output_modified
```

### Generated Output Files & Storage Location

Executing the sanitizer script creates five output files in the **exact directory specified by the `OUTPUT_PREFIX` argument** (typically your shared references directory `/mnt/disks/staff/refs/`):

1. **`<OUTPUT_PREFIX>.renamed.fa`** $\rightarrow$ Cleaned FASTA file with sanitized, BWA-compatible headers (used for reference building and indexing).
2. **`<OUTPUT_PREFIX>.rename_map.tsv`** $\rightarrow$ Lookup map table containing 2 columns: `final_clean_id` and `original_header`.
3. **`<OUTPUT_PREFIX>.name_collisions.tsv`** $\rightarrow$ Log of name collisions and suffix resolution.
4. **`<OUTPUT_PREFIX>.contigs.txt`** $\rightarrow$ Clean 1-column list of sanitized contig IDs, stored directly alongside the FASTA file in `/mnt/disks/staff/refs/`.
5. **`<OUTPUT_PREFIX>.bed`** $\rightarrow$ Clean 3-column BED file (`contig_id\t0\tlength`) covering the full span of each contig in the reference, stored directly in `/mnt/disks/staff/refs/`.

---

## Setup and Execution

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
* **`data_bucket`**: The name of the GCS source bucket where your raw CRAM files are stored.
* **`data_dir`**: The local folder mount path on the server for `data_bucket` (where raw inputs are read).
* **`output_bucket`**: The name of the GCS destination bucket where outputs/logs are deposited. (Note: this is the name of the cloud bucket itself, whereas `output_dir` is the local directory path on the server where this bucket is mounted).
* **`output_dir`**: The local folder mount path on the server for `output_bucket` (where outputs will be written).
* **`work_dir`**: The path to your cloned repository folder on the server (which acts as your active working directory).
* **`scripts_dir`**: The folder inside the repository containing custom Python helper scripts.
* **`db_metadata_dir`**: The folder inside the repository containing BED and contig database files.
* **`samples_list`**: The name of the sample ID list file copied into the repository folder (must have its first line/header written as `SAMPLE`).
* **`dataset` / `phase` / `project`**: Dataset run metadata:
  - **`dataset`**: The prefix/label of the dataset cohort (e.g., `"SSC"`).
  - **`phase`**: The dataset phase number (e.g., `"2"`).
  - **`project`**: An optional sub-project folder identifier.
* **`placeholder_file`**: An empty text file (e.g. `test.txt`) used for creating the file structure/hierarchy on Google Cloud Storage (GCS). Because GCS is an object store, directories do not exist on the cloud unless they contain at least one file.

##### 2. Extraction Variables (Configure ONLY if running the Unmapped Extraction module)
* **`unmapped_extraction`**: Switch (`"on"` or `"off"`) to enable/disable extraction.
* **`ref_genome`**: Full reference genome path.

##### 3. Alignment Variables (Configure if running the Alignment module OR both Extraction + Alignment)
* **`viral_db_alignment`**: Switch (`"on"` or `"off"`) to enable/disable viral database alignment.
* **`align_and_filter`**: Switch (`"on"` or `"off"`) to run full alignments and downstream streaming read filters.
* **`align_only`**: Switch (`"on"` or `"off"`) to run full BWA alignments against the foreign and combined reference databases, skipping downstream streaming read filters.
* **`dual_mode`**: Switch (`"on"` or `"off"`). When set to `"on"`, outputs BOTH Mode 1 (`exogeneSR_viral_clean`) with plasmid/mouse subtracted AND Mode 2 (`exogeneSR_viral_raw`) with all viral hits preserved.
* **Count Statistics Output Files**:
  - `plasmid_pairs.count`: 3-column TSV (`plasmid_pair_count\ttotal_pairs\tplasmid_fraction`)
  - `mouse_pairs.count`: 3-column TSV (`mouse_pair_count\ttotal_pairs\tmouse_fraction`)
  - `unmapped_pairs.count`: 3-column TSV (`unmapped_pair_count\ttotal_pairs\tunmapped_fraction`)

##### 4. Reporting Module Variables (Configure if running Module 3 Reporting)
* **`reporting_module`**: Switch (`"on"` or `"off"`) to enable/disable consolidated reporting.
* **`sample_metadata_file`**: Sample metadata TSV file (`SSC_sample_metadata.tsv`) with columns `sample`, `coverage`, `specimen`, `phase`, `project`.
* **`viral_rename_map_file`**: Name mapping TSV (`HumanViral_Reference_02-07-2022_modified.rename_map.tsv`), generated by `rename_fasta_contigs.sh` (sanitized using the `-ncbi` argument).
* **`viral_bed_file`**: Target BED regions file (`HumanViral_Reference_02-07-2022_modified.bed`), generated automatically during the reference sanitization script.

### Phase 3: Execute the Workflow

Before submitting a job, run Snakemake on the head node to generate the compiled job and dedicated resource configuration files:
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
```bash
batchRun -multibatch samples_p2_base -config config/batch_jobexec_reporting.config -non-spot config/ssc_reporting.job -investigator MDJ -pau 0
```

#### Module 4: SAM Flag Comparison Module (`ssc_flag_comparison.job`)
```bash
batchRun -multibatch samples_all_cohorts -config config/batch_jobexec_flag_comparison.config -non-spot config/ssc_flag_comparison.job -investigator MDJ -pau 0
```

---

## SAM Flag Comparison Module Overview <small>(ssc_flag_comparison.job)</small>

### Overview
Evaluates the clean strategy across 3 SAM flag filtering commands to systematically quantify the difference between aligner `-f 2` flags vs manual bitwise flags.

> [!NOTE]
> This module is not estimating viral abundance. It is estimating high-confidence, reference-discriminating read-pair support. That is useful, however it will undercount viruses when related references share homologous sequence. This is for the reference-unique evidence. We report high-specificity reference-discriminating read support.

### Automated Cohort Consolidation:
Once all batch jobs complete, merge all per-sample TSV files across all phases and projects into a single master TSV report:

```bash
snakemake cohort_flag_comparison_master.tsv --cores 1
```

### Outputs (`/mnt/disks/staff/SSC_hg38_reports/phase2/[project]/<sample_id>/`):

When executed in parallel across multiple samples via `batchRun`, each worker node generates an individual, isolated report file inside its sample folder:
`<sample_id>_viral_report.tsv`

### Merging All Sample Reports into a Single Cohort Master Report
Once all parallel batch jobs complete, combine all individual sample reports from their sample subfolders into a single unified cohort report for downstream statistical analysis:

```bash
# Navigate to your cohort reports folder:
cd /mnt/disks/staff/SSC_hg38_reports/phase2/[project]

# Combine header from the first report + all data rows (skipping duplicate header lines) across sample subfolders:
awk 'FNR==1 && NR!=1{next} {print}' */*_viral_report.tsv > cohort_master_viral_report.tsv
```

## Coverage Calculation Module Overview <small>(ssc_coverage.job)</small>

### Overview
Calculates actual mean human genome read depth across WGS samples directly from original human CRAM files on Google Cloud Storage.

### Execution Command:
```bash
# 1. Compile coverage job and resource configuration files:
snakemake --cores 1

# 2. Execute parallel batch run across cohort samples:
batchRun -multibatch samples_p2_base -config config/batch_jobexec_coverage.config -non-spot config/ssc_coverage.job -investigator MDJ -pau 0
```

### Merging Cloud Coverage Outputs into Master Coverage TSV File:
Once all cloud batch jobs complete, combine all individual sample coverage files directly from Google Cloud Storage into a single master coverage file named with the phase and project:

```bash
PHASE="phase2"
PROJECT="base"  # Or sub-project name e.g. Project_CCDG_...

(
  echo -e "sample\tcoverage\tspecimen\tphase\tproject\tcram_url";
  gsutil -m cat "gs://ml-phi-staff-m277455-p-rsa-us-central1-p-a3d4/SSC_hg38_coverage/${PHASE}/${PROJECT}/**/*.tsv" | awk 'FNR==1 && NR!=1{next} {print}'
) > "/mnt/disks/staff/SSC_hg38_coverage/${PHASE}/cohort_master_coverage_${PHASE}_${PROJECT}.tsv"
```

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
13. `Project` (or `base` if project is empty)
14. `Source_File` (CRAM source file name)
