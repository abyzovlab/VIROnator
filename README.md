# VIROnator

A tool for sensitive detection of viral presence in human samples sequenced by WGS.

## SSC Workflow Modules

This workflow contains two main modules:
1. **Unmapped Extraction Module:** Extracts unmapped human reads from CRAM files.
2. **Viral DB Alignment Module (VIROnator):** Maps unmapped reads to the viral database, filters out decoys/RNA, and runs final read filtering based on targeted BED regions.

### Configuration Variables

Configure these settings in the file `config/ssc_config.yaml` before running the pipeline:

#### 1. Common / Joint Variables
* **`lab_dir` / `staff_dir`**: Directory mount paths.
* **`placeholder_file`**: GCS folder placeholder (`test.txt`).
* **`lab_bucket`**: Source GCS bucket where input CRAM files are stored.
* **`staff_bucket`**: Destination GCS bucket where outputs/logs are deposited.
* **`module` / `phase` / `project`**: Dataset run metadata.
* **`samples_list`**: Sample ID list path (the file must have its first line/header written as `SAMPLE`).
* **`ref_dir`**: Base directory for all reference genomes (GCS mount; all reference files must reside directly in this directory with no subfolders).
* **`scripts_dir`**: Folder in the repo containing custom python scripts (e.g. `scripts/`).
* **`db_metadata_dir`**: Folder in the repo containing BED and contig lists (e.g. `config/db_metadata/`).

#### 2. Extraction Variables
* **`unmapped_extraction`**: Switch (`"on"` or `"off"`) to enable/disable extraction.
* **`ref_genome`**: Full reference genome path.

#### 3. Alignment Variables
* **`viral_db_alignment`**: Switch (`"on"` or `"off"`) to enable/disable viral database alignment.
* **References (`ref_human_no_ebv`, `ref_human_full`, `ref_viral`, etc.)**: Path names relative to `ref_dir`.
* **Metadata files (`viral_contigs_file`, `viral_bed_file`)**: File names stored in `db_metadata_dir`.
* **Tool paths (`bwa_bin`)**: Path to BWA binary.
* **`alignment_mode`**: Select `"align_and_filter"` to run full alignments, or `"filter_only"` to run final read filtering on existing BAMs.

---

### Usage Instructions

For detailed documentation, see [docs/ssc_extraction.md](docs/ssc_extraction.md).

#### 1. Configure the Run
Open `config/ssc_config.yaml` to set up run parameters and toggle modules:
> [!IMPORTANT]
> This pipeline is designed to run on Google Cloud Platform. Ensure that your GCS parameters (such as `lab_bucket`, `staff_bucket`, and project details) are properly configured to point to your GCP cloud buckets.

* **Module Switches:** Indicate with `"on"` or `"off"` which module you want to run:
  - `unmapped_extraction`: Toggles the human unmapped reads extraction module.
  - `viral_db_alignment`: Toggles the viral database alignment module.
  *Toggling these determines which batch job files are compiled by Snakemake and which command you run in Step 4.*
* **Repository & Working Directory:** It is recommended to set `work_dir` in the YAML to the path of your cloned repository, and place your sample list file (e.g., `samples_p2_base`) directly inside the cloned repository directory to keep paths simple.
* **Sample List Format:** The sample list file must have its first line/header written as `SAMPLE`.

#### 2. Run Snakemake
Execute Snakemake locally to compile the configurations:
```bash
snakemake --cores 1
```

#### 3. Load Modules
Load the execution modules on the head node:
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
