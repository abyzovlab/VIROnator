# SSC Dataset Extraction, Alignment & Reporting Workflow

A modular pipeline for sensitive detection and quantitative reporting of viral presence in human samples sequenced by WGS.

The pipeline consists of three independent, sequential modules and one standalone reference preparation utility:

- **Standalone Utility:** Reference Sanitization & BED Generation (`rename_fasta_contigs.sh`)
- **Module 1:** Unmapped Extraction Module
- **Module 2:** Viral DB Alignment Module (VIROnator)
- **Module 3:** Reporting Module

---

## Prerequisites
- Snakemake
- Environment modules / tools: `samtools`, `bwa`, `python` (>= 3.8)

## Storage Recommendations: GCS vs. Local Disk

1. **References Directory:** All reference genomes, database files (human, viral, decoy, rna), BED files, and metadata files must be stored directly in the same flat directory on the GCS mount (`/mnt/disks/staff/refs/`) without any nested subfolders.
2. **Access:** Because batch execution nodes run on ephemeral cloud VMs, files must be stored on the GCS bucket (`/mnt/disks/staff/refs/`) so worker VMs can access them directly.

---

## Standalone Utility: Reference Sanitization & BED Generation

Raw FASTA headers (especially for custom plasmid or viral databases) often contain special characters (Greek letters `α`, `β`, trademarks `™`, spaces, commas, slashes) that break alignment tools (`bwa`, `samtools`).

`scripts/rename_fasta_contigs.sh` is a **standalone utility script** run manually **prior to running any workflow module**.

### Usage Command:
```bash
# For viral database (using -ncbi flag to keep clean accession IDs):
bash scripts/rename_fasta_contigs.sh -i /mnt/disks/staff/refs/HumanViral_Reference_02-07-2022.fa -ncbi

# For plasmid / custom databases:
bash scripts/rename_fasta_contigs.sh -i /mnt/disks/staff/refs/SnapGene_plasmids.fa
```

### Outputs Generated in `/mnt/disks/staff/refs/`:
1. `*.renamed.fa` $\rightarrow$ Cleaned, sanitized FASTA sequence file.
2. `*.contigs.txt` $\rightarrow$ 1-column list of sanitized contig IDs.
3. `*.bed` $\rightarrow$ 3-column BED file (`contig_id\t0\tlength`) used for interval filtering and genome lengths.
4. `*.rename_map.tsv` $\rightarrow$ 2-column TSV (`final_clean_id\toriginal_header`) used by the Reporting Module to map accession IDs back to descriptive virus names.

---

## Setup and Workflow Execution

### Phase 1: Workspace Setup
1. **Clone repository:**
   ```bash
   cd VIROnator
   ```
2. **Copy sample list file:**
   Copy your text file containing sample IDs into the `VIROnator` folder. The first line must be the header `SAMPLE`:
   ```text
   SAMPLE
   SSC01120
   SSC00505
   ```

### Phase 2: Pipeline Configuration (`config/ssc_config.yaml`)

Edit `config/ssc_config.yaml` to toggle modules on/off and configure paths:

#### 1. Common / Joint Variables
* `data_bucket`: GCS source bucket containing raw CRAM files (mounts to `/mnt/disks/lab`).
* `output_bucket`: GCS destination bucket for outputs (mounts to `/mnt/disks/staff`).
* `work_dir`: Path to repository folder on server.
* `ref_dir`: Shared reference directory (`/mnt/disks/staff/refs/`).
* `scripts_dir`: Shared scripts directory (`/mnt/disks/staff/scripts/`).
* `samples_list`: Name of sample ID list file.
* `dataset` / `phase` / `project`: Dataset phase and optional project subfolder identifier.

---

## Module 1: Unmapped Extraction Module

### Overview
Extracts unmapped human reads from raw sample CRAM files to preprocess reads for downstream viral analysis.

### Configuration (`config/ssc_config.yaml`):
```yaml
unmapped_extraction: "on"
ref_cram_decoder: "GRCh38_full_analysis_set_plus_decoy_hla.fa"
```

### Execution Script:
`config/ssc_unmapped.job.template` $\rightarrow$ compiled into `config/ssc_unmapped.job`.

### Outputs:
`/mnt/disks/staff/SSC_hg38_unmapped/phase2/[project]/<sample_id>.unmapped.cram`

---

## Module 2: Viral DB Alignment Module (VIROnator)

### Overview
Maps extracted unmapped reads against the viral reference database, quantifies plasmid and mouse contamination, filters decoy and RNA alignments, and streams filtered read pairs into CRAM format.

### Configuration (`config/ssc_config.yaml`):
```yaml
viral_db_alignment: "on"
contamination_quantification: "on"
dual_mode: "on"
```

### Script Workflow (`config/ssc_alignment.job`):
1. **Candidate Selection:** Single-end BWA alignment to foreign database (`REF_VIR_CONT`).
2. **Contamination Blacklisting:** Blacklists high-confidence plasmid (`plasmid_contigs_file`) and mouse (`mouse_contigs_file`) reads.
3. **Mate Retrieval:** Reconstructs paired-end FASTQ files using `scripts/fq_to_bam_py3.py`.
4. **Combined Alignment & Decoy/RNA Filtering:** Aligns to combined reference, subtracts decoy (`ref_decoy`) and RNA (`ref_rna`) hits, and converts to sorted CRAM format.

### Outputs (`/mnt/disks/staff/SSC_hg38_vironator/phase2/[project]/<sample_id>/`):
* `exogeneSR_viral_clean.sorted.cram` (+ `.crai`) $\rightarrow$ Mode 1: Clean viral CRAM (plasmid/mouse subtracted)
* `exogeneSR_viral_clean.sorted.flags.cram` (+ `.crai`) $\rightarrow$ Mode 1: Clean viral CRAM (strict SAM flag filtering)
* `exogeneSR_viral_raw.sorted.cram` (+ `.crai`) $\rightarrow$ Mode 2: Raw viral CRAM (unsubtracted, if `dual_mode: "on"`)
* `exogeneSR_viral_raw.sorted.flags.cram` (+ `.crai`) $\rightarrow$ Mode 2: Raw viral CRAM (strict SAM flag filtering)
* `plasmid_hits.sorted.cram` (+ `.crai`) $\rightarrow$ Aligned plasmid read pairs
* `mouse_hits.sorted.cram` (+ `.crai`) $\rightarrow$ Aligned mouse read pairs
* `plasmid_pairs.count`, `mouse_pairs.count`, `unmapped_pairs.count` $\rightarrow$ 3-column count TSVs

---

## Module 3: Reporting Module

### Overview
Evaluates all 4 viral output CRAM files (`clean`, `clean_flags`, `raw`, `raw_flags`) per sample, calculates breadth of physical coverage via `samtools mpileup`, matches sample sequencing metadata, and outputs a consolidated 14-column master report.

### Configuration (`config/ssc_config.yaml`):
```yaml
reporting_module: "on"
sample_metadata_file: "sample_metadata.tsv"
viral_rename_map_file: "HumanViral_Reference_02-07-2022_modified.rename_map.tsv" # Sanitized via -ncbi
viral_bed_file: "HumanViral_Reference_02-07-2022_modified.bed" # Sanitization BED file
report_script: "generate_report.py"
```

### Metadata Input Schema (`/mnt/disks/staff/refs/sample_metadata.tsv`):
```tsv
sample	coverage	specimen	phase	project
SSC01120	30.5	Blood_PBMC	2	base
SSC00505	28.2	Saliva	2	base
```

### Outputs (`/mnt/disks/staff/SSC_hg38_reports/phase2/[project]/`):
`<sample_id>_viral_report.tsv`

#### 14-Column Report Format:
1. `Phase`
2. `Project` (or `base` if project is empty)
3. `Source_File` (CRAM source file name)
4. `Sample_ID`
5. `Virus_Accession`
6. `Virus_Length`
7. `Virus_Mapped_Reads`
8. `Normalized_Coverage`
9. `Physical_Coverage` (Percentage formatted without `%` sign, e.g. `100.00`)
10. `Human_Genome_Size`
11. `Sample_Read_Depth`
12. `Viral_Copy_Number`
13. `Specimen`
14. `Virus_Name_Sanitized` (Last column)

---

## Directory & Folder Hierarchy Schema

```text
VIROnator/                                    <-- Cloned Repository folder (work_dir)
├── Snakefile                                 <-- Main Snakemake workflow executor
├── samples_p2_base                           <-- Sample ID list file
├── config/
│   ├── ssc_config.yaml                       <-- Main configuration file
│   ├── batch_jobexec_resources.config.template
│   ├── ssc_unmapped.job.template             <-- Module 1 job template
│   ├── ssc_alignment.job.template            <-- Module 2 job template
│   └── ssc_reporting.job.template            <-- Module 3 job template
├── rules/
│   └── ssc.smk                               <-- Snakemake rules
└── scripts/
    ├── generate_report.py                    <-- Module 3 Python reporting script
    ├── init_ref.sh                           <-- Reference initialization helper
    ├── fq_to_bam_py3.py                      <-- Read pair extraction helper
    ├── aln_match_filter_py3.py               <-- Decoy & RNA match filter
    ├── filter_reads.py                       <-- Candidate streaming filter
    └── rename_fasta_contigs.sh              <-- Standalone reference sanitization utility

/mnt/disks/staff/refs/                        <-- Shared Reference Directory (ref_dir)
├── *.renamed.fa                              <-- Cleaned reference FASTA files
├── *.contigs.txt                             <-- 1-column contig ID list files
├── *.bed                                     <-- 3-column target region BED files
├── *.rename_map.tsv                          <-- FASTA header sanitization lookup maps
└── sample_metadata.tsv                       <-- Sample sequencing depth & specimen metadata

/mnt/disks/staff/                             <-- Output GCS Mount Directory (output_dir)
├── SSC_hg38_unmapped/                        <-- Module 1 unmapped CRAM outputs
├── SSC_hg38_vironator/                       <-- Module 2 viral alignment CRAM outputs
└── SSC_hg38_reports/                         <-- Module 3 consolidated TSV reports
```
