# SSC Dataset Extraction, Alignment & Reporting Workflow

This tool works modularly to process and align sequencing datasets. It consists of three sequential modules and one standalone reference preparation utility:

1. **Standalone Utility:** Reference Sanitization & BED Generation (`rename_fasta_contigs.sh`)
2. **Module 1:** Unmapped Extraction Module (`ssc_unmapped.job`)
3. **Module 2:** Viral DB Alignment Module (`ssc_alignment.job`)
4. **Module 3:** Reporting Module (`ssc_reporting.job`)

---

## Prerequisites & Storage Setup

- **Snakemake**
- **Cluster modules:** `samtools`, `bwa`, `python`, `jobexec/2.0.1`

### Storage Recommendations: GCS vs. Local Disk
1. **References Directory:** All reference genomes and database files (human, viral, decoy, rna, bed, metadata) must be stored directly in the same flat directory on the GCS mount (`/mnt/disks/staff/refs/`) without any nested subfolders.
2. **Access:** Because the batch jobs run on separate ephemeral cloud VM instances, they must be stored on the GCS bucket (`/mnt/disks/staff/refs/`) so worker VMs can access them.

---

# Standalone Utility: Reference Sanitization (`rename_fasta_contigs.sh`)

When adopting a new vector or plasmid reference dataset, raw FASTA header names often contain special characters (Greek letters like `α`, `β`, `γ`, trademarks `™`, `®`, symbols `+`, `%`, `#`, commas, slashes, and spaces). These special characters cause parsing failures or silent truncation in alignment tools like `bwa mem` and `samtools`.

`scripts/rename_fasta_contigs.sh` is a **standalone utility script** to preprocess and sanitize FASTA contig header names. **It is run manually prior to reference indexing when adding new reference files.**

### How Headers Are Sanitized
1. **Whitespace Replacement:** All spaces, tabs, and non-breaking whitespace characters are immediately converted into underscores `_`.
2. **Greek Transliteration:** Converts Greek characters to English words (`α` $\rightarrow$ `alpha`, `β` $\rightarrow$ `beta`, `λ` $\rightarrow$ `lambda`, `μ` $\rightarrow$ `mu`, `φ` $\rightarrow$ `phi`).
3. **Trademark Removal:** Strips `™`, `®`, `©`, `(TM)`, `(R)`.
4. **Symbol Translation:** Converts mathematical/strand symbols (`(+)` $\rightarrow$ `_plus_`, `(-)` $\rightarrow$ `_minus_`, `+` $\rightarrow$ `_plus_`, `&` $\rightarrow$ `_and_`).
5. **Special Character Conversion:** Any remaining non-alphanumeric character is converted to `_`, and consecutive delimiters are collapsed.
6. **Collision Deduplication:** If two distinct raw headers clean down to the exact same string, numeric suffixes (`_2`, `_3`) are automatically appended to guarantee that **every contig ID is unique**.
7. **NCBI Accession Mode (`-n` / `--ncbi`):** When processing NCBI-formatted FASTA headers (e.g. `>NC_001829 Adeno-associated virus - 4, complete genome.`), passing `--ncbi` strips everything after the first whitespace to keep the accession ID (`NC_001829`) as the contig ID while preserving the original raw header string in `rename_map.tsv`.

### Execution Syntax

Command options (`-p` / `--prefix` and `-n` / `--ncbi`) **must be placed before the input FASTA file**:

```bash
./scripts/rename_fasta_contigs.sh [OPTIONS] <INPUT_FASTA> <OUTPUT_PREFIX>
```

#### Usage Examples
```bash
# 1. Standard plasmid / vector sanitization:
./scripts/rename_fasta_contigs.sh SnapGene.fa SnapGene_plasmids_modified

# 2. Mouse reference sanitization (prepends "mm39_" prefix):
./scripts/rename_fasta_contigs.sh -p "mm39_" mm39_raw.fa mm39_ms_modified

# 3. NCBI viral reference sanitization (strips descriptions after accession ID):
./scripts/rename_fasta_contigs.sh --ncbi HumanViral_Reference_02-07-2022.fa HumanViral_Reference_02-07-2022_modified
```

### Generated Output Files in `/mnt/disks/staff/refs/`:
1. **`<OUTPUT_PREFIX>.renamed.fa`** $\rightarrow$ Cleaned FASTA file with sanitized, BWA-compatible headers.
2. **`<OUTPUT_PREFIX>.rename_map.tsv`** $\rightarrow$ Lookup map table containing 2 columns: `final_clean_id` and `original_header`.
3. **`<OUTPUT_PREFIX>.contigs.txt`** $\rightarrow$ Clean 1-column list of sanitized contig IDs.
4. **`<OUTPUT_PREFIX>.bed`** $\rightarrow$ Clean 3-column BED file (`contig_id\t0\tlength`).

---

# Module 1: Unmapped Extraction Module (`ssc_unmapped.job`)

### Overview
Extracts unmapped human reads from raw sample CRAM files to preprocess reads for downstream viral analysis.

### Configuration (`config/ssc_config.yaml`):
```yaml
unmapped_extraction: "on"
ref_cram_decoder: "GRCh38_full_analysis_set_plus_decoy_hla.fa"
```

### Batch Execution Command:
```bash
batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_unmapped.job -investigator MDJ -pau 0
```

### Outputs:
`/mnt/disks/staff/SSC_hg38_unmapped/phase2/[project]/<sample_id>.unmapped.cram`

---

# Module 2: Viral DB Alignment Module (`ssc_alignment.job`)

### Overview
Maps extracted unmapped reads against the viral reference database, quantifies plasmid and mouse contamination, filters decoy/RNA alignments, and streams filtered read pairs into CRAM format.

### Configuration (`config/ssc_config.yaml`):
```yaml
viral_db_alignment: "on"
contamination_quantification: "on"
dual_mode: "on"
```

### Script Workflow:
1. **Step 1: Alignment to Viral Database** $\rightarrow$ Single-end BWA alignment to foreign database (`REF_VIR_CONT`).
2. **Step 2: Read Filtering & Contamination Submodule** $\rightarrow$ Blacklists plasmid (`plasmid_contigs_file`) and mouse (`mouse_contigs_file`) reads.
3. **Step 3: Mate Retrieval & FASTQ Reconstruction** $\rightarrow$ Reconstructs paired FASTQs using `scripts/fq_to_bam_py3.py`.
4. **Step 4: Paired Alignment & Streaming Read Filtering** $\rightarrow$ Aligns to combined reference, subtracts decoy (`ref_decoy`) and RNA (`ref_rna`) hits, and converts to sorted CRAM format.

### Batch Execution Command:
```bash
batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_alignment.job -investigator MDJ -pau 0
```

### Outputs (`/mnt/disks/staff/SSC_hg38_vironator/phase2/[project]/<sample_id>/`):
* `exogeneSR_viral_clean.sorted.cram` (+ `.crai`) $\rightarrow$ Mode 1: Clean viral CRAM (plasmid/mouse subtracted)
* `exogeneSR_viral_clean.sorted.flags.cram` (+ `.crai`) $\rightarrow$ Mode 1: Clean viral CRAM (strict SAM flags)
* `exogeneSR_viral_raw.sorted.cram` (+ `.crai`) $\rightarrow$ Mode 2: Raw viral CRAM (unsubtracted, if `dual_mode: "on"`)
* `exogeneSR_viral_raw.sorted.flags.cram` (+ `.crai`) $\rightarrow$ Mode 2: Raw viral CRAM (strict SAM flags)
* `plasmid_hits.sorted.cram`, `mouse_hits.sorted.cram` $\rightarrow$ Contamination CRAMs
* `plasmid_pairs.count`, `mouse_pairs.count`, `unmapped_pairs.count` $\rightarrow$ 3-column count TSVs

---

# Module 3: Reporting Module (`ssc_reporting.job`)

### Overview
Evaluates all 4 viral output CRAM files (`clean`, `clean_flags`, `raw`, `raw_flags`) per sample, calculates physical breadth of coverage via `samtools mpileup`, matches sample metadata (`sample_metadata.tsv`), and outputs a consolidated 14-column master report.

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

### Batch Execution Command:
```bash
batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_reporting.job -investigator MDJ -pau 0
```

### Outputs (`/mnt/disks/staff/SSC_hg38_reports/phase2/[project]/`):
`<sample_id>_viral_report.tsv`

#### 14-Column Master Report Schema:
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
├── test.txt                                  <-- Placeholder file for GCS directory structure creation
├── samples_p2_base                           <-- Your sample ID list file
│
├── config/                                   <-- Configuration & Job Templates folder
│   ├── ssc_config.yaml                       <-- Main pipeline configuration file
│   ├── batch_jobexec_resources.config.template
│   ├── ssc_unmapped.job.template             <-- Module 1 job template
│   ├── ssc_alignment.job.template            <-- Module 2 job template
│   └── ssc_reporting.job.template            <-- Module 3 job template
│
├── rules/                                    <-- Snakemake Logic & Compilation Rules
│   └── ssc.smk                               <-- Pipeline compilation rules & sync logic
│
└── scripts/                                  <-- Python & Bash helper scripts (scripts_dir)
    ├── generate_report.py                    <-- Module 3 Python reporting script
    ├── init_ref.sh                           <-- Reference initialization helper
    ├── fq_to_bam_py3.py                      <-- FASTQ/BAM read pair extraction helper
    ├── aln_match_filter_py3.py               <-- Decoy & RNA alignment match filter
    ├── filter_reads.py                       <-- Streaming alignment candidate filter
    └── rename_fasta_contigs.sh              <-- Standalone FASTA header cleaner & BED generator

/mnt/disks/staff/refs/                        <-- Shared Reference Mount Directory (ref_dir)
├── *.renamed.fa                              <-- Cleaned reference FASTA files
├── *.contigs.txt                             <-- 1-column contig ID list files
├── *.bed                                     <-- 3-column target region BED files
├── *.rename_map.tsv                          <-- FASTA header sanitization lookup maps
└── sample_metadata.tsv                       <-- Sample sequencing depth & specimen metadata

/mnt/disks/staff/                             <-- Output GCS Mount Directory (output_dir)
├── SSC_hg38_unmapped/                        <-- Module 1 unmapped extraction outputs
├── SSC_hg38_vironator/                       <-- Module 2 viral alignment outputs
└── SSC_hg38_reports/                         <-- Module 3 consolidated TSV reports
```
