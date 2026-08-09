# VIROnator

A tool for sensitive detection and quantitative reporting of viral presence in human samples sequenced by WGS.

This workflow works modularly to process and align sequencing datasets. It consists of three sequential modules and one standalone reference preparation utility:

- **Standalone Utility:** Reference Sanitization & BED Generation (`rename_fasta_contigs.sh`)
- **Module 1:** Unmapped Extraction Module
- **Module 2:** Viral DB Alignment Module (VIROnator)
- **Module 3:** Reporting Module

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

---

## Workflow Execution & Modules

### Step-by-Step Setup

1. **Clone repository:**
   ```bash
   cd VIROnator
   ```
2. **Copy sample list file:**
   Copy your sample ID list file into the `VIROnator` folder. The first line must be `SAMPLE`:
   ```text
   SAMPLE
   SSC01120
   SSC00505
   ```
3. **Configure Settings:**
   Edit `config/ssc_config.yaml` to toggle module switches on/off.

---

### Module 1: Unmapped Extraction Module
* Extracts unmapped human reads from CRAM files to preprocess reads for viral analysis.
* **Toggle:** `unmapped_extraction: "on"` in `config/ssc_config.yaml`.
* **Output Path:** `/mnt/disks/staff/SSC_hg38_unmapped/phase2/[project]/<sample_id>.unmapped.cram`

---

### Module 2: Viral DB Alignment Module (VIROnator)
* Maps unmapped reads to the viral database, blacklists plasmid/mouse contamination, filters decoy/RNA sequences, and streams filtered read pairs into CRAM format.
* **Toggle:** `viral_db_alignment: "on"`, `contamination_quantification: "on"`, `dual_mode: "on"`.
* **Output Path:** `/mnt/disks/staff/SSC_hg38_vironator/phase2/[project]/<sample_id>/`
  - `exogeneSR_viral_clean.sorted.cram` (+ `.crai`)
  - `exogeneSR_viral_clean.sorted.flags.cram` (+ `.crai`)
  - `exogeneSR_viral_raw.sorted.cram` (+ `.crai`)
  - `exogeneSR_viral_raw.sorted.flags.cram` (+ `.crai`)
  - `plasmid_hits.sorted.cram`, `mouse_hits.sorted.cram`
  - `plasmid_pairs.count`, `mouse_pairs.count`, `unmapped_pairs.count`

---

### Module 3: Reporting Module
* Consolidates metrics across all 4 output CRAM files (`clean`, `clean_flags`, `raw`, `raw_flags`), calculates physical breadth of coverage via `samtools mpileup`, matches sample metadata (`sample_metadata.tsv`), and produces unified 14-column TSV reports.
* **Toggle:** `reporting_module: "on"` in `config/ssc_config.yaml`.
* **Output Path:** `/mnt/disks/staff/SSC_hg38_reports/phase2/[project]/<sample_id>_viral_report.tsv`

#### 14-Column Master Report Schema:
1. `Phase`
2. `Project` (or `base` if project is empty)
3. `Source_File`
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

### Running Batch Execution Jobs

```bash
# Module 1 Batch Execution:
batchRun -multibatch samples_p2_base_test -config config/batch_jobexec_resources.config -non-spot config/ssc_unmapped.job -investigator MDJ -pau 0

# Module 2 Batch Execution:
batchRun -multibatch samples_p2_base_test -config config/batch_jobexec_resources.config -non-spot config/ssc_alignment.job -investigator MDJ -pau 0

# Module 3 Batch Execution:
batchRun -multibatch samples_p2_base_test -config config/batch_jobexec_resources.config -non-spot config/ssc_reporting.job -investigator MDJ -pau 0
```

For full details, see [docs/ssc_extraction.md](docs/ssc_extraction.md).
