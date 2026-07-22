# SSC Unmapped Read Extraction Module

This module prepares output directories on cloud mounts and dynamically configures the `ssc.job` and `batch_jobexec_resources.config` files inside the `config/` directory.

## Prerequisites
- Snakemake
- samtools (in container environment)

## Setup and Execution

1. **Configure:** Open `config/ssc_config.yaml` and set your target variables:
   ```yaml
   # Target Run Metadata
   phase: "4"
   project: "Project_REI_12386_B01_GRM_WGS.cram.2017-06-02"
   
   # Storage Buckets
   lab_bucket: "ml-phi-abyzov-alexej-m124423-rsa-us-central1-p-d3c4" # Source data bucket
   staff_bucket: "ml-phi-staff-m277455-p-rsa-us-central1-p-a3d4"     # Output destination bucket
   ```

2. **Initialize:** Run Snakemake to create folders and generate your configuration files:
   ```bash
   snakemake --cores 1
   ```

3. **Load Environment:** Load the required cluster/batch modules before running the batch submission:
   ```bash
   module load samtools jobexec/2.0.1
   ```

4. **Submit:** Launch parallel cloud extraction:
   ```bash
   batchRun -multibatch samples_p4_clean_base -config config/batch_jobexec_resources.config -non-spot config/ssc.job
   ```
