# VIROnator

A tool for sensitive detection and quantitative reporting of viral presence in human samples sequenced by WGS.

This workflow works modularly to process and align sequencing datasets. It consists of three sequential modules and one standalone reference preparation utility:

1. **Standalone Utility:** Reference Sanitization & BED Generation (`rename_fasta_contigs.sh`)
2. **Module 1:** Unmapped Extraction Module (`ssc_unmapped.job`)
3. **Module 2:** Viral DB Alignment Module (`ssc_alignment.job`)
4. **Module 3:** Reporting Module (`ssc_reporting.job`)

---

# Standalone Utility: Reference Sanitization (`rename_fasta_contigs.sh`)

When adopting a new vector or plasmid reference dataset, raw FASTA header names often contain special characters (Greek letters like `α`, `β`, `γ`, trademarks `™`, `®`, symbols `+`, `%`, `#`, commas, slashes, and spaces). These special characters cause parsing failures or silent truncation in alignment tools like `bwa mem` and `samtools`.

`scripts/rename_fasta_contigs.sh` is an optional **standalone utility script** to preprocess and sanitize FASTA contig header names. **It is run manually prior to reference indexing when adding new reference files.**

### Usage Command:
```bash
# 1. Standard plasmid / vector sanitization:
./scripts/rename_fasta_contigs.sh SnapGene.fa SnapGene_plasmids_modified

# 2. NCBI viral reference sanitization (using -ncbi flag to keep clean accession IDs):
./scripts/rename_fasta_contigs.sh --ncbi HumanViral_Reference_02-07-2022.fa HumanViral_Reference_02-07-2022_modified
```

---

# Module 1: Unmapped Extraction Module (`ssc_unmapped.job`)

Extracts unmapped human reads from raw sample CRAM files to preprocess reads for downstream viral analysis.

* **Configuration:** Set `unmapped_extraction: "on"` in `config/ssc_config.yaml`.
* **Output Path:** `/mnt/disks/staff/SSC_hg38_unmapped/phase2/[project]/<sample_id>.unmapped.cram`
* **Batch Command:**
  ```bash
  batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_unmapped.job -investigator MDJ -pau 0
  ```

---

# Module 2: Viral DB Alignment Module (`ssc_alignment.job`)

Maps extracted unmapped reads against the viral reference database, quantifies plasmid and mouse contamination, filters decoy/RNA alignments, and streams filtered read pairs into CRAM format.

* **Configuration:** Set `viral_db_alignment: "on"`, `contamination_quantification: "on"`, `dual_mode: "on"` in `config/ssc_config.yaml`.
* **Output Path:** `/mnt/disks/staff/SSC_hg38_vironator/phase2/[project]/<sample_id>/`
  - `exogeneSR_viral_clean.sorted.cram` (+ `.crai`)
  - `exogeneSR_viral_clean.sorted.flags.cram` (+ `.crai`)
  - `exogeneSR_viral_raw.sorted.cram` (+ `.crai`)
  - `exogeneSR_viral_raw.sorted.flags.cram` (+ `.crai`)
  - `plasmid_hits.sorted.cram`, `mouse_hits.sorted.cram`
  - `plasmid_pairs.count`, `mouse_pairs.count`, `unmapped_pairs.count`
* **Batch Command:**
  ```bash
  batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_alignment.job -investigator MDJ -pau 0
  ```

---

# Module 3: Reporting Module (`ssc_reporting.job`)

Consolidates metrics across all 4 output CRAM files (`clean`, `clean_flags`, `raw`, `raw_flags`), calculates physical breadth of coverage via `samtools mpileup`, matches sample metadata (`sample_metadata.tsv`), and produces unified 14-column TSV reports.

* **Configuration:** Set `reporting_module: "on"` in `config/ssc_config.yaml`.
* **Output Path:** `/mnt/disks/staff/SSC_hg38_reports/phase2/[project]/<sample_id>_viral_report.tsv`
* **Batch Command:**
  ```bash
  batchRun -multibatch samples_p2_base -config config/batch_jobexec_resources.config -non-spot config/ssc_reporting.job -investigator MDJ -pau 0
  ```

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

## Detailed Instructions

For full documentation and setup details, see [docs/ssc_extraction.md](docs/ssc_extraction.md).
