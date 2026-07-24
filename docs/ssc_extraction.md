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
Open `config/ssc_config.yaml` in a text editor. Here are the key variables you need to configure and what they mean in plain language:

> [!IMPORTANT]
> This pipeline is designed to run on Google Cloud Platform. Ensure that your GCS parameters (such as `data_bucket`, `output_bucket`, and project details) are properly configured to point to your GCP cloud buckets.

#### Global Settings:
* **`data_bucket`**: The name of your Google Cloud storage bucket where the raw input dataset (CRAM files) is stored.
* **`output_bucket`**: The name of your Google Cloud storage bucket where you want all output files, results, and logs to be saved.
* **`work_dir`**: The absolute path to your cloned repository folder on your machine (e.g., `/home/user/working/VIROnator`).
* **`samples_list`**: The name of the sample list file you copied into the repository in Phase 1 (e.g., `samples_p2_base`).

#### Module Toggles (ON / OFF):
You specify which parts of the pipeline to run by setting these switches to `"on"` or `"off"`:
* **`unmapped_extraction`**: Set to `"on"` to extract human unmapped reads from CRAM files. Set to `"off"` to skip this step.
* **`viral_db_alignment`**: Set to `"on"` to align the preprocessed unmapped reads to the viral database. Set to `"off"` to skip this step.
* *Toggling these determines which batch job files Snakemake compiles and which command you run in Step 4.*

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

### ⚠️ Sample List Path Resolution
* **Relative Paths:** If you specify a relative path for `samples_list` (e.g., `"samples_p2_base"`), Snakemake dynamically appends it to the `work_dir` configured in the YAML file (e.g. `work_dir/samples_p2_base`).
* **Absolute Paths:** If your samples list is stored outside the working directory, you must provide an absolute path starting with `/` (e.g. `/home/user/working/samples_p2_base`) to bypass automatic appending.
