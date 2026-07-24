# VIROnator

A tool for sensitive detection of viral presence in human samples sequenced by WGS.

This workflow works modularly to process and align sequencing datasets. It consists of two sequential modules whose parameters and execution are fully controlled by the user via the configuration file (`config/ssc_config.yaml`):

1. **Unmapped Extraction Module:** Extracts unmapped human reads from CRAM files. This preprocesses the reads for subsequent viral analysis.
2. **Viral DB Alignment Module (VIROnator):** Maps the preprocessed unmapped reads to the viral database, filters out decoy/RNA sequences, and runs final read filtering based on targeted BED regions.

### Usage Instructions

For detailed documentation, see [docs/ssc_extraction.md](docs/ssc_extraction.md).

### Step-by-Step Setup

#### Phase 1: Prepare your Workspace
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

#### Phase 2: Configure your Settings
Open `config/ssc_config.yaml` in a text editor. Configure these settings before running the pipeline:

> [!IMPORTANT]
> This pipeline is designed to run on Google Cloud Platform. Ensure that your GCS parameters (such as `data_bucket`, `output_bucket`, and project details) are properly configured to point to your GCP cloud buckets.

##### 1. Common / Joint Variables
* **`data_dir` / `output_dir`**: local folder mount paths. `data_dir` is the path where raw input datasets (BAM/CRAM files) are stored/mounted. `output_dir` is the path where the preprocessed unmapped reads and all final output files will be written.
* **`placeholder_file`**: An empty text file (e.g. `test.txt`) used for creating the file structure/hierarchy on Google Cloud Storage (GCS). Because GCS is an object store, directories do not exist on the cloud unless they contain at least one file.
* **`data_bucket`**: The name of the source GCS bucket where input CRAM files are stored.
* **`output_bucket`**: The name of the destination GCS bucket where outputs/logs are deposited. (Note: this is the name of the cloud bucket itself, whereas `output_dir` is the local directory path on the server where this bucket is mounted).
* **`module` / `phase` / `project`**: Dataset run metadata.
* **`samples_list`**: Sample ID list path (the file must have its first line/header written as `SAMPLE` and be copied into the repository folder as specified in Phase 1).
* **`work_dir`**: The absolute path to your cloned repository folder on your machine (e.g., `/home/user/working/VIROnator`).
* **`ref_dir`**: Base directory for all reference genomes (GCS mount; all reference files must reside directly in this directory with no subfolders).
* **`scripts_dir`**: Folder in the repo containing custom python scripts (e.g. `scripts/`).
* **`db_metadata_dir`**: Folder in the repo containing BED and contig lists (e.g. `config/db_metadata/`).

##### 2. Extraction Variables (Configure ONLY if running the Unmapped Extraction module)
* **`unmapped_extraction`**: Switch (`"on"` or `"off"`) to enable/disable extraction.
* **`ref_genome`**: Full reference genome path.

##### 3. Alignment Variables (Configure if running the Alignment module OR both Extraction + Alignment)
* **`viral_db_alignment`**: Switch (`"on"` or `"off"`) to enable/disable viral database alignment.
* **References (`ref_human_no_ebv`, `ref_human_full`, `ref_viral`, etc.)**: Path names relative to `ref_dir`.
* **Metadata files (`viral_contigs_file`, `viral_bed_file`)**: File names stored in `db_metadata_dir`.
* **Tool paths (`bwa_bin`)**: Path to BWA binary.
* **`alignment_mode`**: Select `"align_and_filter"` to run full alignments, or `"filter_only"` to run final read filtering on existing BAMs.

#### Phase 3: Compile and Execute

1. **Run Snakemake:** Execute Snakemake locally to compile the configurations:
   ```bash
   snakemake --cores 1
   ```

2. **Load Modules:** Load the execution modules on the head node:
   ```bash
   module load samtools bwa python jobexec/2.0.1
   ```

#### 4. Submit the Batch Jobs
Based on the module you set to `"on"` in Step 1, run the corresponding command:
* **For Unmapped Extraction:**
  ```bash
  batchRun -multibatch <samples_file> -config config/batch_jobexec_resources.config -non-spot config/ssc_unmapped.job -investigator <pi_id> -pau <pau_id>
  ```
* **For Viral Alignment:**
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

### ⚠️ Sample List Path Resolution
* **Relative Paths:** If you specify a relative path for `samples_list` (e.g., `"samples_p2_base"`), Snakemake dynamically appends it to the `work_dir` configured in the YAML file (e.g. `work_dir/samples_p2_base`).
* **Absolute Paths:** If your samples list is stored outside the working directory, you must provide an absolute path starting with `/` (e.g. `/home/user/working/samples_p2_base`) to bypass automatic appending.
