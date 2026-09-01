# SSC Dataset Extraction & Alignment Workflow

This tool works modularly to process and align sequencing datasets. It consists of sequential modules whose parameters and execution are fully controlled by the user via the configuration file (`config/ssc_config.yaml`):

1. **Module 1: Coverage Calculation Module:** Calculates per-sample human genome sequencing read depth and coverage metrics.
2. **Module 2: Unmapped Extraction Module:** Extracts unmapped human reads from CRAM files to preprocess candidate reads for viral detection.
3. **Module 3: Viral DB Alignment Module (VIROnator):** Maps candidate reads to the viral database, filters out decoy/RNA sequences, and runs targeted BED filtering.
4. **Module 4: Consolidated Reporting Module:** Consolidates metrics across output CRAM files (`clean`, `clean_flags`, `raw`), evaluates physical coverage, merges sample metadata, and produces unified 14-column TSV reports.
5. **Module 5: Cohort Stats Summary Module:** Generates cohort-wide statistical snapshots (`cohort_stats_summary.tsv` and `cohort_stats_summary.md`).
6. **Module 6: Distributions & Plots Module:** Produces publication-ready 2-panel 300 DPI TIFF figures and per-virus classification tables.
7. **Module 7: SAM Flag Comparison Module:** Evaluates clean strategy alignment metrics across SAM flag filtering variants.
8. **Module 8: NCBI RefSeq Refinement Module:** Re-aligns candidate reads against complete RefSeq viral genomes using Bowtie2 multimapping (`-k 10`) and a Two-Tier Classification Cascade (`Reference_Unique` vs `Species_Supportive`).

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

## Standalone Utility: Taxonomy Index Builder <small>(make_taxonomy_index.py)</small>

`scripts/make_taxonomy_index.py` is a standalone helper script that pre-processes NCBI taxonomy dump files (`nodes.dmp`, `names.dmp`) and accession-to-taxid mapping files (`nucl_gb.accession2taxid`) into a unified, 14-column viral taxonomy lookup index (`viral_reference_taxonomy_index.tsv`).

### Downloading NCBI Taxonomy Source Files
Before executing the script, download the latest NCBI taxonomy dump and accession mapping files:

```bash
mkdir -p /tmp/ncbi_taxdump && cd /tmp/ncbi_taxdump

# Download NCBI taxonomy nodes and names:
wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz nodes.dmp names.dmp

# Download GenBank nucleotide accession-to-taxid mapping:
wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz
```

### Execution Command
```bash
python3 scripts/make_taxonomy_index.py \
    --accession-taxid /tmp/ncbi_taxdump/nucl_gb.accession2taxid \
    --nodes /tmp/ncbi_taxdump/nodes.dmp \
    --names /tmp/ncbi_taxdump/names.dmp \
    --output config/db_metadata/viral_reference_taxonomy_index.tsv
```

### Generated Output 14-Column Index Schema (`viral_reference_taxonomy_index.tsv`)
1. `accession` — RefSeq / GenBank accession ID (e.g. `NC_001829`)
2. `accession_version` — Accession version (e.g. `NC_001829.1`)
3. `accession_taxid` — NCBI TaxID assigned to accession
4. `accession_name` — Scientific name assigned to accession
5. `species_taxid` — TaxID of the species rank
6. `species_name` — Species scientific name
7. `genus_taxid` — TaxID of the genus rank
8. `genus_name` — Genus scientific name
9. `family_taxid` — TaxID of the family rank
10. `family_name` — Family scientific name
11. `order_taxid` — TaxID of the order rank
12. `order_name` — Order scientific name
13. `realm_taxid` — TaxID of the realm rank
14. `realm_name` — Realm scientific name

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

##### Module 1 Variables (Coverage Calculation Module)
* **`coverage_module`**: Switch (`"on"` or `"off"`) to enable/disable coverage calculation.
* **`coverage_out_dirname`**: Directory for output coverage TSVs.

##### Module 2 Variables (Unmapped Extraction Module)
* **`unmapped_extraction`**: Switch (`"on"` or `"off"`) to enable/disable extraction.
* **`ref_genome`**: Full reference genome path.

##### Module 3 Variables (Viral DB Alignment Module)
* **`viral_db_alignment`**: Switch (`"on"` or `"off"`) to enable/disable viral database alignment.
* **`align_and_filter`**: Switch (`"on"` or `"off"`) to run full alignments and downstream streaming read filters.
* **`align_only`**: Switch (`"on"` or `"off"`) to run full BWA alignments against the foreign and combined reference databases, skipping downstream streaming read filters.
* **`dual_mode`**: Switch (`"on"` or `"off"`). When set to `"on"`, outputs BOTH Mode 1 (`exogeneSR_viral_clean`) with plasmid/mouse subtracted AND Mode 2 (`exogeneSR_viral_raw`) with all viral hits preserved.

##### Module 4 Variables (Consolidated Reporting Module)
* **`reporting_module`**: Switch (`"on"` or `"off"`) to enable/disable consolidated reporting.
* **`sample_metadata_file`**: Sample metadata TSV file (`SSC_sample_metadata.tsv`) with columns `sample`, `coverage`, `specimen`, `phase`, `project`.
* **`viral_rename_map_file`**: Name mapping TSV (`HumanViral_Reference_02-07-2022_modified.rename_map.tsv`), generated by `rename_fasta_contigs.sh` (sanitized using the `-ncbi` argument).
* **`viral_bed_file`**: Target BED regions file (`HumanViral_Reference_02-07-2022_modified.bed`), generated automatically during the reference sanitization script.

### Phase 3: Execute the Workflow

Before submitting a job, run Snakemake on the head node to generate the compiled job and dedicated resource configuration files:
```bash
snakemake --cores 1
```

#### Module 1: Coverage Calculation Module (`ssc_coverage.job`)
```bash
batchRun -multibatch samples_list.txt -config config/batch_jobexec_coverage.config -non-spot config/ssc_coverage.job -investigator MDJ -pau 0
```

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

#### Module 5 & 6: Cohort Stats & Distributions Modules (Local Execution)
```bash
snakemake cohort_stats_summary.tsv --cores 1
snakemake generate_distributions --cores 1
```

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

---

## NCBI RefSeq Refinement Module Overview <small>(ssc_refinement.job)</small>

### Overview & Purpose
The Refinement Module re-aligns candidate viral reads (from the `clean` strategy) against an expanded reference index containing **all complete RefSeq viral genomes** for detected viral groups. It uses Bowtie2 multimapping (`-k 10`) and a **Hierarchical Two-Tier Classification Cascade** (`Reference_Unique` vs `Species_Supportive`).

> [!NOTE]
> This module is not estimating viral abundance. It is estimating high-confidence, reference-discriminating read-pair support. That is useful, however it will undercount viruses when related references share homologous sequence. This is for the reference-unique evidence. We report high-specificity reference-discriminating read support.

---

### The 4 Submodules / Stages

#### Stage 1: Taxonomy Index Construction (`make_taxonomy_index_stage`)
* **Purpose**: Parses NCBI taxdump files (`nodes.dmp`, `names.dmp`) and accession-to-taxid mapping files (`nucl_gb.accession2taxid`) cached in `/mnt/disks/staff/refs/ncbi_taxdump/`, cross-referencing your database FASTA (`HumanViral_Reference_02-07-2022.fa`) entries up to the highest rank (**Realm**).
* **Snakemake Command**:
  ```bash
  snakemake config/db_metadata/viral_reference_taxonomy_index.tsv --cores 1
  ```
* **Output File**: `config/db_metadata/viral_reference_taxonomy_index.tsv` (14-column taxonomy index TSV).

#### Stage 2: NCBI Complete RefSeq Genome Fetcher (`ncbi_download_stage`)
* **Purpose**: Extracts unique viral accessions from Column 2 (`Virus_Accession`) of `master_all_cohorts_viral_report_final.tsv`, looks up their `species_taxid` (or `genus_taxid`), downloads **all complete nucleotide RefSeq genomes** under detected taxIDs from NCBI Entrez into `/mnt/disks/staff/refs/ncbi_download/`, and runs contig header sanitization (`rename_fasta_contigs.sh --ncbi`).
* **Snakemake Command**:
  ```bash
  snakemake config/ncbi_download.completed --cores 1
  ```
* **Outputs**:
  - `/mnt/disks/staff/refs/ncbi_download/refinement_species_<TAXID>.renamed.fa` (Sanitized RefSeq FASTA files)
  - `/mnt/disks/staff/SSC_hg38_refinement/logs/refinement_ncbi_fetch.log` (Audit fetch log)
  - `config/ncbi_download.completed` (Stage completion token)

#### Stage 3: Date-Stamped Combined Reference & Bowtie2 Index Builder (`build_combined_ref_stage`)
* **Purpose**: Concatenates Human (GRCh38), Mouse (`mm39`), Plasmids (`SnapGene`), and all downloaded RefSeq viral FASTAs into a single date-stamped reference and builds `samtools faidx` (`.fai`) and `bowtie2-build` (`.bt2` / `.bt2l`) binary indices.
* **Snakemake Command**:
  ```bash
  snakemake config/ssc_refinement.job --cores 1
  ```
* **Outputs (`/mnt/disks/staff/refs/`)**:
  - `combined_human_mouse_plasmids_refinement_viral_2026-08-29.fa` (Combined FASTA)
  - `combined_human_mouse_plasmids_refinement_viral_2026-08-29.fa.fai` (FASTA Index)
  - `combined_human_mouse_plasmids_refinement_viral_2026-08-29.fa.1.bt2` through `.4.bt2` (or `.1.bt2l` through `.4.bt2l`) (Bowtie2 Binary Indices)

#### Stage 4: Multimapping Alignment & Two-Tier Classification Cascade (`refinement_alignment_stage`)
* **Purpose**: Executes Bowtie2 multimapping alignment (`-k 10 --very-sensitive`) per sample, followed by the **Hierarchical Two-Tier Classification Cascade**:
  - **Tier 1 (Strict Primary Unique Check)**: MAPQ $\ge 20$, $F\,12 == 0$, $F\,2048 == 0$, CIGAR $150M$/$151M$, $150 \le |TLEN| \le 1500$, unique accession hit $\rightarrow$ `REFERENCE_UNIQUE`. (Strictly locked in Tier 1! Does NOT inspect alternative alignments or species-supportive assignments).
  - **Tier 2 (Species-Supportive Cascade Check)**: Evaluates remaining reads to test if target accessions share the same `species_taxid` $\rightarrow$ `SPECIES_SUPPORTIVE`.
* **Batch Execution Command**:
  ```bash
  batchRun -multibatch config/refinement_samples.tsv -config config/batch_jobexec_refinement.config -non-spot config/ssc_refinement.job -investigator MDJ -pau 0
  ```

> [!IMPORTANT]
> **Post-Cluster Consolidation**: Because `batchRun` dispatches sample alignments to compute nodes in parallel, each job writes its own individual TSV report under `/mnt/disks/staff/SSC_hg38_refinement/phase<PHASE>/[PROJECT]/<SAMPLE_ID>/`.
> Once all parallel cluster jobs have completed, you **must run the consolidation command separately on the head node** to aggregate all per-sample reports into `cohort_refinement_master.tsv`:

* **Master Consolidation Command (Run on Head Node After Cluster Jobs Finish)**:
  ```bash
  snakemake cohort_refinement_master.tsv --cores 1
  ```

---

### Output Directory Structure & File Contents

#### 1. Per-Sample Refinement Output Directory:
`Path: /mnt/disks/staff/SSC_hg38_refinement/phase<PHASE>/[PROJECT]/<SAMPLE_ID>/`

* **`exogeneSR_viral_refinement.cram` & `.crai`**:
  - **Content**: Compressed, indexed Bowtie2 multimapping CRAM file containing all aligned candidate reads against the date-stamped combined reference.
* **`bowtie2_refinement.log`**:
  - **Content**: Bowtie2 alignment summary log detailing overall alignment rate, mapped pair count, and multimapping distribution.
* **`<PHASE>_<PROJECT>_<SAMPLE_ID>_refined_report.tsv`**:
  - **Content**: 11-column per-sample TSV report summarizing refined read counts and lower-bound copy numbers per viral species.

#### 2. Master Cohort Output TSV File (`cohort_refinement_master.tsv`):
`Path: /mnt/disks/staff/VIROnator/cohort_refinement_master.tsv`

#### 11-Column Schema:
1. `Sample_ID` — Sample accession identifier (e.g., `SS0012978`)
2. `Phase` — Cohort phase tag (e.g., `phase1`)
3. `Project` — Cohort project tag (or `base` if empty)
4. `Virus_Accession` — Best-supported RefSeq Accession ID (e.g., `NC_001829.1`)
5. `Species_TaxID` — NCBI Species Taxonomy ID (e.g., `10509`)
6. `Species_Name` — Taxonomic Species Scientific Name (e.g., `Human gammaherpesvirus 4`)
7. `Reference_Unique_Reads` — Read pairs uniquely supporting this specific reference (Tier 1 strict count)
8. `Species_Supportive_Reads` — Read pairs supporting this species group (Tier 2 group-collapsed count)
9. `Total_Refined_Reads` — Total refined read pairs (`Reference_Unique_Reads` + `Species_Supportive_Reads`)
10. `Refined_Copy_Number` — Reference-discriminating copy-number lower bound
11. `Classification_Status` — Call status (`CONFIRMED_UNIQUE` or `SPECIES_COLLAPSED`)

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
