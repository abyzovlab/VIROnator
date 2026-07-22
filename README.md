# VIROnator

A tool for sensitive detection of viral presence in human samples sequenced by WGS.

## SSC Unmapped Read Extraction Module

This module prepares output directories on cloud mounts and configures the batch submission template (`config/ssc.job`) and resource file (`config/batch_jobexec_resources.config`) for a selected project and phase.

### Configuration Variables

Set these variables in the configuration file `config/ssc_config.yaml` before running the pipeline:

* **`phase`**: The target phase of the dataset (e.g., `"4"`).
* **`project`**: The target project directory containing raw samples (e.g., `"Project_REI_12386_B01_GRM_WGS.cram.2017-06-02"`).
* **`lab_bucket`**: The source Google Cloud Storage (GCS) bucket where raw CRAM files are stored (mounts to `/mnt/disks/lab`).
* **`staff_bucket`**: The destination GCS bucket where unmapped output CRAMs are deposited (mounts to `/mnt/disks/staff`).

### Usage Instructions

For detailed documentation, see [docs/ssc_extraction.md](docs/ssc_extraction.md).

#### 1. Configure the Run
Open the configuration file `config/ssc_config.yaml` and set your target phase, project, and GCS buckets.

#### 2. Run Snakemake
Execute Snakemake locally to initialize the cloud directories and generate both configured config files:
```bash
snakemake --cores 1
```

#### 3. Load Modules
Load the execution modules on the cluster:
```bash
module load samtools jobexec/2.0.1
```

#### 4. Submit the Batch Job
Submit the parallel extraction jobs to the cloud:
```bash
batchRun -multibatch samples_p4_clean_base -config config/batch_jobexec_resources.config -non-spot config/ssc.job
```
