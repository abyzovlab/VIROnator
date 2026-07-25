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
* **References (`ref_human_no_ebv`, `ref_human_full`, `ref_viral`, etc.)**: Path names relative to `ref_dir`.
* **Metadata files (`viral_contigs_file`, `viral_bed_file`)**: File names stored in `db_metadata_dir`.
* **Tool paths (`bwa_bin`)**: Path to BWA binary.
* **`alignment_mode`**: Select `"align_and_filter"` to run full alignments, or `"filter_only"` to run final read filtering on existing BAMs.

#### Directory & Folder Hierarchy Schema
To help you understand how files are organized and how to configure the path variables, here is the structure of the cloned repository folder on your server:

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
│   └── ssc_align.job.template
│
├── config/db_metadata/                       <-- Database region coordinate folder (db_metadata_dir)
│   ├── targeted_contigs.txt                  <-- List of target viral chromosomes
│   └── targeted_regions.bed                  <-- Coordinates of target viral regions
│
└── scripts/                                  <-- Python helper scripts folder (scripts_dir)
    ├── init_ref.sh                           <-- Script to initialize and index reference genomes
    ├── fq_to_bam_py3.py                      <-- Script to convert FASTQ reads back to BAM format
    ├── aln_match_filter_py3.py               <-- Script to map and match alignments
    └── filter_reads.py                       <-- Script to run final target read filtering
```

* **`VIROnator/` (Repository Root / `work_dir`):** The primary folder where Snakemake is executed. You copy your sample list text file (e.g. `samples_p2_base`) directly here.
* **`config/`:** Holds Snakemake compilation rules and batch job templates.
* **`config/db_metadata/` (`db_metadata_dir`):** Holds the specific coordinates and chromosome lists used during the alignment module to filter reads down to targets.
* **`scripts/` (`scripts_dir`):** Holds the custom helper Python and Bash scripts that Snakemake compiles and runs on the cloud batch nodes.

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
     batchRun -multibatch <samples_file> -config config/batch_jobexec_resources.config -non-spot config/ssc_align.job -investigator <pi_id> -pau <pau_id>
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
