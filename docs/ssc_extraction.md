# SSC Dataset Extraction & Alignment Workflow

This workflow contains two main modules:
1. **Unmapped Extraction Module:** Extracts unmapped human reads from CRAM files.
2. **Viral DB Alignment Module (VIROnator):** Maps unmapped reads to the viral database, filters out decoys/RNA, and runs final read filtering based on targeted BED regions.

## Prerequisites
- Snakemake
- Cluster modules: `samtools`, `bwa`, `python`, `jobexec/2.0.1`

## Storage Recommendations: GCS vs. Local Disk

1. **References Directory:** All reference genomes and database files (human, viral, decoy, rna) must be stored directly in the same flat directory on the GCS mount (specifically `/mnt/disks/staff/refs/`) without any nested subfolders.
2. **Access:** Because the batch jobs run on separate ephemeral cloud VM instances, they must be stored on the GCS bucket (`/mnt/disks/staff/refs/`) so the batch execution nodes can access them.

---

## Setup and Execution

1. **Configure:** Open `config/ssc_config.yaml` to set up run parameters and toggle modules:
   * **`work_dir`**: The path to your local repository folder (where the sample lists are stored).
   * **Sample List format**: The sample list file (e.g., `samples_p2_clean_base`) **must** have its first line/header written as `SAMPLE`.
   * **Module Switches**:
     - To run extraction only: set `unmapped_extraction: "on"` and `viral_db_alignment: "off"`.
     - To run alignment only: set `unmapped_extraction: "off"` and `viral_db_alignment: "on"`.
     - To run both: set both to `"on"`.

2. **Initialize:** Run Snakemake locally to create directories and compile job files:
   ```bash
   snakemake --cores 1
   ```

3. **Load Environment:** Load the required environment modules:
   ```bash
   module load samtools bwa python jobexec/2.0.1
   ```

4. **Submit Batch Jobs:**
   - **For Unmapped Extraction:**
     ```bash
     batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_unmapped.job
     ```
   - **For Viral Alignment:**
     ```bash
     batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_align.job
     ```
