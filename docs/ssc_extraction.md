# SSC Dataset Extraction & Alignment Workflow

This tool works modularly to process and align sequencing datasets. It consists of two sequential modules whose parameters and execution are fully controlled by the user via the configuration file (`config/ssc_config.yaml`):

1. **Unmapped Extraction Module:** Extracts unmapped human reads from CRAM files. This preprocesses the reads for subsequent viral analysis.
2. **Viral DB Alignment Module (VIROnator):** Maps the preprocessed unmapped reads to the viral database, filters out decoy/RNA sequences, and runs final read filtering based on targeted BED regions.

## Prerequisites
- Snakemake
- Cluster modules: `samtools`, `bwa`, `python`, `jobexec/2.0.1`

## Storage Recommendations: GCS vs. Local Disk

1. **References Directory:** All reference genomes and database files (human, viral, decoy, rna) must be stored directly in the same flat directory on the GCS mount (specifically `/mnt/disks/staff/refs/`) without any nested subfolders.
2. **Access:** Because the batch jobs run on separate ephemeral cloud VM instances, they must be stored on the GCS bucket (`/mnt/disks/staff/refs/`) so the batch execution nodes can access them.

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
* **`filter_only`**: Switch (`"on"` or `"off"`) to bypass alignment steps and run final streaming read filtering on pre-aligned CRAMs/BAMs.
* **Count Statistics Output Files**:
  - `plasmid_pairs.count`: 3-column TSV (`plasmid_pair_count\ttotal_pairs\tplasmid_fraction`)
  - `mouse_pairs.count`: 3-column TSV (`mouse_pair_count\ttotal_pairs\tmouse_fraction`)
  - `unmapped_pairs.count`: 3-column TSV (`unmapped_pair_count\ttotal_pairs\tunmapped_fraction`)

#### Directory & Folder Hierarchy Schema
To help you understand how files are organized and how to configure the path variables, here is the structure of the cloned repository folder and mounted reference directory:

```text
VIROnator/                                    <-- Cloned Repository folder (This is your work_dir)
├── Snakefile                                 <-- Snakemake workflow executor blueprint
├── test.txt                                  <-- Empty placeholder file (placeholder_file)
├── samples_p2_base                           <-- Your sample ID list file (Copied here in Phase 1)
│
├── config/                                   <-- Configuration templates folder
│   ├── ssc_config.yaml                       <-- Main configuration file (Your settings go here)
│   ├── batch_jobexec_resources.config.template
│   ├── ssc_unmapped.job.template
│   └── ssc_alignment.job.template
│
├── config/db_metadata/                       <-- Database region coordinate folder (db_metadata_dir)
│   └── targeted_contigs.txt                  <-- Legacy target chromosome lists (if used)
│
└── scripts/                                  <-- Python & Bash helper scripts folder (scripts_dir)
    ├── init_ref.sh                           <-- Script to initialize and index reference genomes
    ├── fq_to_bam_py3.py                      <-- Script to convert FASTQ reads back to BAM format
    ├── aln_match_filter_py3.py               <-- Script to map and match alignments
    ├── filter_reads.py                       <-- Script to run final target read filtering
    └── rename_fasta_contigs.sh              <-- Standalone script to clean and rename FASTA headers

/mnt/disks/staff/refs/                        <-- Shared Reference Mount Directory (ref_dir)
├── *.renamed.fa                              <-- Cleaned reference FASTA files
├── *.contigs.txt                             <-- 1-column contig ID list files
├── *.bed                                     <-- 3-column target region BED files
└── *.rename_map.tsv                          <-- FASTA header sanitization lookup tables
```

* **`VIROnator/` (Repository Root / `work_dir`):** The primary folder where Snakemake is executed. You copy your sample list text file (e.g. `samples_p2_base`) directly here.
* **`config/`:** Holds Snakemake compilation rules and batch job templates.
* **`/mnt/disks/staff/refs/` (`ref_dir`):** Holds all reference FASTA files, `.contigs.txt` lists, and `.bed` coordinate files used by `ssc_alignment.job`.
* **`scripts/` (`scripts_dir`):** Holds the custom helper Python and Bash scripts that Snakemake compiles and runs on the cloud batch nodes.

---

## Alignment Module Script Overview (`ssc_alignment.job`)

The alignment job script (`config/ssc_alignment.job.template` compiled into `config/ssc_alignment.job`) runs a 4-step pipeline for viral candidate detection and filtering:

1. **Step 1: Alignment to Viral Database**
   Extracts candidate reads from the unmapped CRAM file using `samtools fastq` and aligns them single-end to the target reference with `bwa mem`. Candidate mapped reads are filtered against the target contig list (`viral_contigs_file`).
2. **Step 2: Low-Complexity & Candidate Read ID Extraction**
   Extracts unique read IDs from the mapped candidates and prepares the candidate whitelist (`viral_reads_se.keep`).
3. **Step 3: Mate Retrieval & Paired FASTQ Reconstruction**
   Executes `scripts/fq_to_bam_py3.py` against the unmapped CRAM file to retrieve both mates for every candidate read ID and reconstructs paired FASTQ files (`viral_1.fq`, `viral_2.fq`).
4. **Step 4: Paired-End Alignment & Streaming Read Filtering**
   Aligns read pairs to the combined human + viral reference using `bwa mem -Y`. Outputs are streamed through decoy/RNA reference filters and `scripts/filter_reads.py` to produce the final filtered CRAM (`exogeneSR_viral_filtered.sorted.cram`).

---

## Standalone Utility: FASTA Header Cleaner (`rename_fasta_contigs.sh`)

When adopting a new vector or plasmid reference dataset, raw FASTA header names often contain special characters (Greek letters like $\alpha, \beta, \gamma$, trademarks $\text{™}, \text{®}$, symbols $+$, $\%$, $\#$, commas, slashes, and spaces). These special characters cause parsing failures or silent truncation in alignment tools like `bwa mem` and `samtools`.

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

> [!IMPORTANT]
> **Configuration Requirement:** In `config/ssc_config.yaml`, you **must specify the sanitized contig names and sanitized FASTA paths**. BWA and samtools operate on the sanitized headers; providing unsanitized names in `ssc_config.yaml` will prevent candidate contig filtering from matching.

---

## Reference Database Naming & Automated Merging Workflow

To keep reference genomes and metadata structured and reproducible across datasets, specify the individual component references in `config/ssc_config.yaml`:

```yaml
# In config/ssc_config.yaml (relative to ref_dir /mnt/disks/staff/refs/)
ref_human_full: "GRCh38_full_analysis_set_plus_decoy_hla_modified.fa"
ref_human_no_ebv: "GRCh38_full_analysis_set_plus_decoy_hla_no_ebv_modified.fa"
ref_viral_base: "HumanViral_Reference_02-07-2022.fa"
ref_plasmids: "SnapGene_plasmids_modified.fa"
ref_mouse: "mm39_ms_modified.fa"
ref_vir_cont: "HumanViral_Reference_02-07-2022_SnapGene_plasmids_modified_mm39_ms_modified.fa"
ref_combined_default: "combined_GRCh38_full_analysis_set_plus_decoy_hla_no_ebv_modified_HumanViral_Reference_02-07-2022_SnapGene_plasmids_modified_mm39_ms_modified.fa"
ref_decoy: "GRCh38_Decoy_Reference_4-30-2019.fa"
ref_rna: "GRCh38_RNA_Reference_4-30-2019.fa"
```

### Automated Combination & Smart Skip Logic
Inside the alignment job script (`ssc_alignment.job`), reference merging and indexing are dynamically managed:

1. **Foreign Database (`REF_VIR_CONT`):** If `REF_VIR_CONT` (Viral Base + Plasmids + Mouse) and its indices (`.ann`, `.fai`) already exist, **no merging or indexing is performed**. Otherwise, it automatically concatenates the component FASTAs and generates the BWA and samtools indices.
2. **Combined Reference Database (`REF_COMBINED`):** If `REF_COMBINED` (Human No-EBV + `REF_VIR_CONT`) and its indices (`.ann`, `.fai`) already exist, **no merging or indexing is performed**. Otherwise, it constructs the combined FASTA and indices using clean component naming (`combined_<human_basename>_<viral_basename>.fa`).

### File Naming Summary
* **Raw / Unstripped References:** `<name>_raw.fa` (e.g., `SnapGene_plasmids_raw.fa`)
* **Sanitized / Cleaned References:** `<name>_modified.fa` or `<name>_cleaned.renamed.fa`
* **Sanitized Contig List Files:** `<name>_modified.txt` or `<name>.contigs.txt`
* **Foreign DB Reference (`REF_VIR_CONT`):** `<viral_base>_<plasmids>_<mouse>.fa`
* **Combined DB Reference (`REF_COMBINED`):** `combined_<human_basename>_<viral_basename>.fa`

### Phase 3: Compile and Execute

1. **Initialize:** Run Snakemake locally to compile the configurations:
   ```bash
   snakemake --cores 1
   ```

2. **Load Environment:** Load the required environment modules:
   ```bash
   module load samtools bwa python jobexec/2.0.1
   ```

4. **Submit Batch Jobs:**
   Based on the module you set to `"on"` in Step 1, run the corresponding command:
   - **For Unmapped Extraction:**
     ```bash
     batchRun -multibatch <samples_file> -config config/batch_jobexec_resources.config -non-spot config/ssc_unmapped.job -investigator <pi_id> -pau <pau_id>
     ```
   - **For Viral Alignment:**
     ```bash
     batchRun -multibatch <samples_file> -config config/batch_jobexec_resources.config -non-spot config/ssc_alignment.job -investigator <pi_id> -pau <pau_id>
     ```

---

### 💡 Running Test Runs / Dry Runs
Before running a full production pipeline, it is recommended to run a test run to verify the setup. To isolate your test outputs from production, change these four variables in `config/ssc_config.yaml`:

1. **GCS Output Folders:**
   * `unmapped_out_dirname: "SSC_hg38_unmapped_test"`
   * `vironator_out_dirname: "SSC_hg38_vironator_test"`
2. **Jobexec State Folders:**
   * `unmapped_jobexec_dirname: "jobexec_unmapped_test"`
   * `vironator_jobexec_dirname: "jobexec_vironator_test"`

*Remember to run `snakemake --cores 1` (or `snakemake --cores 1 -f` to force recompile) after changing the YAML file to compile these new paths into your batch job scripts.*

---

### Sample List Path Resolution
* **Relative Paths:** If you specify a relative path for `samples_list` (e.g., `"samples_p2_base"`), Snakemake dynamically appends it to the `work_dir` configured in the YAML file (e.g. `work_dir/samples_p2_base`).
* **Absolute Paths:** If your samples list is stored outside the working directory, you must provide an absolute path starting with `/` (e.g. `/home/user/working/samples_p2_base`) to bypass automatic appending.
