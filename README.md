# VIROnator

A tool for sensitive detection of viral presence in human samples sequenced by WGS.

This workflow works modularly to process and align sequencing datasets. It consists of sequential modules whose parameters and execution are fully controlled by the user via the configuration file (`config/ssc_config.yaml`):

1. **Unmapped Extraction Module:** Extracts unmapped human reads from CRAM files. This preprocesses the reads for subsequent viral analysis.
2. **Viral DB Alignment Module (VIROnator):** Maps the preprocessed unmapped reads to the viral database, filters out decoy/RNA sequences, and runs final read filtering based on targeted BED regions.
3. **Reporting Module:** Consolidates metrics across output CRAM files (`clean`, `clean_flags`, `raw`, `raw_flags`), evaluates physical coverage, merges sample metadata, and produces unified 14-column TSV reports.

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

### Usage Instructions

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
* **`contamination_quantification`**: Submodule switch (`"on"` or `"off"`). When set to `"on"`, high-confidence plasmid (vector) and mouse contigs are blacklisted from the clean viral CRAM and quantified as separate pair count statistics (`plasmid_pairs.count`, `mouse_pairs.count`). When `"off"`, all viral-hitting reads are passed directly down the pipeline.
* **`dual_mode`**: Switch (`"on"` or `"off"`). When set to `"on"`, outputs BOTH Mode 1 (`exogeneSR_viral_clean`) with plasmid/mouse subtracted AND Mode 2 (`exogeneSR_viral_raw`) with all viral hits preserved.
* **References & Metadata (`ref_dir` / `/mnt/disks/staff/refs/`)**: Path to the reference directory on worker VMs. **All reference FASTA files, contig ID lists (`.contigs.txt`), and target region BED files (`.bed`) are stored directly in `ref_dir`**.
  - `ref_human_full`, `ref_human_no_ebv`, `ref_viral_base`, `ref_plasmids`, `ref_mouse`, `ref_decoy`, `ref_rna`
  - `viral_contigs_file`, `viral_bed_file`, `plasmid_contigs_file`, `mouse_contigs_file`
* **Tool paths (`bwa_bin`)**: Path to BWA binary.
* **`align_and_filter`**: Switch (`"on"` or `"off"`) to run full alignments and downstream streaming read filters.
* **`align_only`**: Switch (`"on"` or `"off"`) to run full BWA alignments against the foreign and combined reference databases, skipping downstream streaming read filters.
* **Count Statistics Output Files**:
  - `plasmid_pairs.count`: 3-column TSV (`plasmid_pair_count\ttotal_pairs\tplasmid_fraction`)
  - `mouse_pairs.count`: 3-column TSV (`mouse_pair_count\ttotal_pairs\tmouse_fraction`)
  - `unmapped_pairs.count`: 3-column TSV (`unmapped_pair_count\ttotal_pairs\tunmapped_fraction`)

##### 4. Reporting Module Variables (Configure if running Module 3 Reporting)
* **`reporting_module`**: Switch (`"on"` or `"off"`) to enable/disable consolidated reporting.
* **`sample_metadata_file`**: Sample metadata TSV file (`sample_metadata.tsv`) with columns `sample`, `coverage`, `specimen`, `phase`, `project`.
* **`viral_rename_map_file`**: Name mapping TSV (`HumanViral_Reference_02-07-2022_modified.rename_map.tsv`), generated by `rename_fasta_contigs.sh` (sanitized using the `-ncbi` argument).
* **`viral_bed_file`**: Target BED regions file (`HumanViral_Reference_02-07-2022_modified.bed`), generated automatically during the reference sanitization script.

#### Directory & Folder Hierarchy Schema
To help you understand how files are organized and how to configure the path variables, here is the structure of the repository, shared references, and output directories:

```text
VIROnator/                                    <-- Cloned Repository folder (work_dir)
├── Snakefile                                 <-- Main Snakemake workflow executor
├── test.txt                                  <-- Placeholder file for GCS directory structure creation
├── samples_p2_base                           <-- Your sample ID list file
│
├── config/                                   <-- Configuration & Job Templates folder
│   ├── ssc_config.yaml                       <-- Main pipeline configuration file
│   ├── batch_jobexec_resources.config.template
│   ├── ssc_unmapped.job.template             <-- Phase 1 unmapped extraction template
│   └── ssc_alignment.job.template            <-- Phase 2 viral alignment template
│
├── rules/                                    <-- Snakemake Logic & Compilation Rules
│   └── ssc.smk                               <-- Pipeline compilation rules & sync logic
│
└── scripts/                                  <-- Python & Bash helper scripts (scripts_dir)
    ├── init_ref.sh                           <-- Reference initialization & combined FASTA builder
    ├── fq_to_bam_py3.py                      <-- FASTQ/BAM read pair extraction helper
    ├── aln_match_filter_py3.py               <-- Decoy & RNA alignment match filter
    ├── filter_reads.py                       <-- Streaming alignment candidate filter
    └── rename_fasta_contigs.sh              <-- Standalone FASTA header cleaner & BED generator

/mnt/disks/staff/refs/                        <-- Shared Reference Mount Directory (ref_dir)
├── *.renamed.fa                              <-- Cleaned reference FASTA files
├── *.contigs.txt                             <-- 1-column contig ID list files
├── *.bed                                     <-- 3-column target region BED files
├── *.rename_map.tsv                          <-- FASTA header sanitization lookup maps
└── *.name_collisions.tsv                     <-- Log of duplicate header suffix resolutions

/mnt/disks/staff/                             <-- Output GCS Mount Directory (output_dir)
├── SSC_hg38_unmapped/                        <-- Phase 1 unmapped extraction outputs
└── SSC_hg38_vironator/                       <-- Phase 2 viral alignment & count outputs
    └── phase2/<sample_id>/
        ├── exogeneSR_viral_clean.sorted.cram  <-- Mode 1: Clean viral CRAM (plasmid/mouse subtracted)
        ├── exogeneSR_viral_clean.sorted.flags.cram
        ├── exogeneSR_viral_raw.sorted.cram    <-- Mode 2: Raw viral CRAM (unsubtracted, if dual_mode=on)
        ├── exogeneSR_viral_raw.sorted.flags.cram
        ├── plasmid_hits.sorted.cram          <-- Aligned plasmid/vector CRAM
        ├── mouse_hits.sorted.cram            <-- Aligned mouse contamination CRAM
        ├── plasmid_pairs.count               <-- 3-column plasmid count TSV
        ├── mouse_pairs.count                 <-- 3-column mouse count TSV
        └── unmapped_pairs.count              <-- 3-column unmapped count TSV
```

### Phase 3: Execute the Workflow

#### Option A: Running Module 1 (Unmapped Extraction ONLY)
If you only want to extract unmapped human reads:
1. Open `config/ssc_config.yaml`:
   - Set `unmapped_extraction: "on"`
   - Set `viral_db_alignment: "off"`
   - Set `reporting_module: "off"`
2. Run Snakemake on the head node:
   ```bash
   snakemake --cores 1
   ```
3. Submit the batch execution job to the cloud:
   ```bash
   batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_unmapped.job -investigator MDJ -pau 0
   ```

#### Option B: Running Module 2 (Viral DB Alignment ONLY)
If you already extracted unmapped reads and want to align them to the viral database:
1. Open `config/ssc_config.yaml`:
   - Set `unmapped_extraction: "off"`
   - Set `viral_db_alignment: "on"`
   - Set `reporting_module: "off"`
2. Run Snakemake on the head node:
   ```bash
   snakemake --cores 1
   ```
3. Submit the batch execution job to the cloud:
   ```bash
   batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_alignment.job -investigator MDJ -pau 0
   ```

#### Option C: Running Module 3 (Reporting Module ONLY)
If you already ran viral alignment and want to generate consolidated reports:
1. Open `config/ssc_config.yaml`:
   - Set `unmapped_extraction: "off"`
   - Set `viral_db_alignment: "off"`
   - Set `reporting_module: "on"`
2. Run Snakemake on the head node:
   ```bash
   snakemake --cores 1
   ```
3. Submit the batch execution job to the cloud:
   ```bash
   batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_reporting.job -investigator MDJ -pau 0
   ```
