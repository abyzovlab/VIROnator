# VIROnator

A tool for sensitive detection and quantitative reporting of viral presence in human samples sequenced by WGS.

This workflow works modularly to process and align sequencing datasets. It consists of sequential modules whose parameters and execution are fully controlled by the user via the configuration file (`config/ssc_config.yaml`):

1. **Standalone Utility:** Reference Sanitization & BED Generation (`rename_fasta_contigs.sh`)
2. **Module 1:** Coverage Calculation Module (`ssc_coverage.job`)
3. **Module 2:** Unmapped Extraction Module (`ssc_unmapped.job`)
4. **Module 3:** Alignment Module (`ssc_alignment.job`)
5. **Module 4:** Reporting Module (`ssc_reporting.job`)
6. **Module 5:** Cohort Stats Summary Module (`generate_stats.py`)
7. **Module 6:** Distributions & Plots Module (`generate_distributions.py`)
8. **Module 7:** SAM Flag Comparison Module (`ssc_flag_comparison.job`)
9. **Module 8:** NCBI RefSeq Refinement Module (`ssc_refinement.job`)
10. **Module 9:** Heatmaps Module (`generate_heatmap.py`)

---

## Standalone Utility: Reference Sanitization <small>(rename_fasta_contigs.sh)</small>

When adopting a new vector or plasmid reference dataset, raw FASTA header names often contain special characters (Greek letters like `α`, `β`, `γ`, trademarks `™`, `®`, symbols `+`, `%`, `#`, commas, slashes, and spaces). These special characters cause parsing failures or silent truncation in alignment tools like `bwa mem` and `samtools`.

`scripts/rename_fasta_contigs.sh` is a **standalone utility script** to preprocess and sanitize FASTA contig header names. **It is run manually prior to reference indexing when preparing reference databases.**

> [!IMPORTANT]
> All sanitized output files (FASTA `.renamed.fa`, rename map `.rename_map.tsv`, contigs list `.renamed.contigs.txt`, and BED file `.renamed.bed`) **must be stored directly in `/mnt/disks/staff/refs/`** so that all downstream pipeline modules can access them properly.

> [!NOTE]
> The `.renamed` extension tag is automatically appended to generated files by default (e.g. `.renamed.fa`, `.renamed.contigs.txt`, `.renamed.bed`). Therefore, the `OUTPUT_PREFIX` command argument should **not** contain the word `"renamed"`, but simply the base reference name (e.g. `HumanViral_Reference_02-07-2022_modified` or `mm39_modified`).

#### Sanitization Commands & Examples:

1. **NCBI Viral Reference Database Sanitization**:
   Use `--ncbi` to extract clean accession IDs and avoid overly long contig names coming from full record descriptions:
   ```bash
   scripts/rename_fasta_contigs.sh --ncbi HumanViral_Reference_02-07-2022.fa HumanViral_Reference_02-07-2022_modified
   ```
   **Output Files Generated**:
   - `HumanViral_Reference_02-07-2022_modified.renamed.fa` (Sanitized reference FASTA)
   - `HumanViral_Reference_02-07-2022_modified.rename_map.tsv` (Accession header mapping TSV)
   - `HumanViral_Reference_02-07-2022_modified.renamed.contigs.txt` (Contigs ID list)
   - `HumanViral_Reference_02-07-2022_modified.renamed.bed` (Target BED regions file)

2. **Mouse Genome Reference Sanitization**:
   Use `-p "mm39_"` to prefix contig names and differentiate mouse contigs from human/viral contigs:
   ```bash
   scripts/rename_fasta_contigs.sh -p "mm39_" mm39.fa mm39_modified
   ```

3. **Plasmid Reference Database Sanitization**:
   ```bash
   scripts/rename_fasta_contigs.sh SnapGene.fa SnapGene_modified
   ```

---

### Standalone Utility: Viral Name Shortener <small>(shorten_viral_names.py)</small>

`scripts/shorten_viral_names.py` is a **standalone utility script** that processes a `.renamed_map.tsv` file (created during reference sanitization) and applies rule-based text replacements to abbreviate and clean viral descriptive names for cleaner figure labels in stats and heatmaps.

It supports both local file paths and Google Cloud Storage (`gs://`) URIs.

#### Shortening Rules Applied:
- **Comma Stripping**: Automatically removes all commas (`,`) from headers first.
- **Removals**: Strips genomic descriptors such as ` genomic RNA complete genome strain:`, ` complete genome strain:`, ` genomic RNA complete`, ` genomic DNA partial`, ` genomic DNA`, ` genomic RNA`, ` complete virion genome`, ` virion genome`, ` complete genome`, `sense strand`, and space-dash-space combinations (` - `).
- **Replacements**: Standardizes viral family/species names (e.g. `Torque teno mini virus` $\rightarrow$ `TTV mini`, `Torque teno virus` $\rightarrow$ `TTV`, `Adeno-associated virus` $\rightarrow$ `AAV`, `Simian adenovirus` $\rightarrow$ `Simian Ad`, `Human papillomavirus` / `Human papillomavirus virus` $\rightarrow$ `HPV`, `Human immunodeficiency virus` $\rightarrow$ `HIV`, `Human alphaherpesvirus 1-3` $\rightarrow$ `HHV-1`-`HHV-3`, `Human gammaherpesvirus 4` $\rightarrow$ `EBV`, `Human betaherpesvirus 5` $\rightarrow$ `CMV`, `Human betaherpesvirus 6A/6B/7` $\rightarrow$ `HHV-6A`/`HHV-6B`/`HHV-7`, `Molluscum contagiosum virus` $\rightarrow$ `MCV`, `Porcine endogenous retrovirus` $\rightarrow$ `PERV`).

#### Usage Example:
```bash
# Process map file stored on Google Cloud Storage (automatically saves to gs://<output_bucket>/refs/<map_name>_short.tsv)
python3 scripts/shorten_viral_names.py gs://<output_bucket>/refs/HumanViral_Reference_02-07-2022_modified.renamed_map.tsv
```

---

## Workflow Execution & Modules

For detailed documentation, see [docs/ssc_extraction.md](docs/ssc_extraction.md).

### Step-by-Step Setup

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
Open `config/ssc_config.yaml` in a text editor. Configure these settings before running the pipeline:

> [!IMPORTANT]
> This pipeline is designed to run on Google Cloud Platform. Ensure that your GCS parameters (such as `data_bucket`, `output_bucket`, and project details) are properly configured to point to your GCP cloud buckets.

#### 1. Common / Joint Variables
* **`data_dir`**: Shared local mount directory for reading raw input files (`/mnt/disks/lab`).
* **`data_subdir`**: Middle subdirectory path inside `data_dir` (e.g. `"tertiary/SSC_hg38/WGS"`, or `""` for flat datasets).
* **`input_keyword`**: Optional middle keyword following sample ID dot (e.g. `"final"`, `"sorted"`, or `""` for none).
* **`input_suffix`**: Target file extension following keyword (e.g. `"cram"`, `"bam"`).
* **`output_dir`**: Shared local mount directory for output files (`/mnt/disks/staff`).
* **`work_dir`**: Active repository root directory on head node (default `"."`).
* **`data_bucket`**: Source GCS bucket containing raw sequence files.
* **`output_bucket`**: Destination GCS bucket for output data and cluster execution logs.
* **`gcp_account`**: GCP Service Account email address for cluster execution.
* **`pi`**: Principal Investigator initials/tag (e.g. `"MDJ"`).
* **`pau`**: Project Allocation Unit code (e.g. `0`).
* **`task_name`**: Cluster job log identifier prefix tag.
* **`dataset`**: Dataset cohort prefix label (e.g. `"SSC"`).
* **`genome_build`**: Human reference genome assembly build label (e.g. `"hg38"`).
* **`phase`**: Dataset phase identifier (e.g. `"2"`, or `"none"` / `""` for flat datasets).
* **`project`**: Sub-project folder identifier (e.g. `"base"`, or `"none"` / `""` for flat datasets).
* **`samples_list`**: Name of text file listing target sample IDs (header must be `SAMPLE`).
* **`ref_dir`**: Directory path where reference genomes and database files are stored (`/mnt/disks/staff/refs`).
* **`scripts_dir`**: Repository folder path containing helper scripts (`/mnt/disks/staff/scripts`).

#### 4. Module 4: Reporting Module Sub-Switches & Execution

Module 4 handles viral quantification, taxonomy annotation, and master report generation. It contains **two independent submodules** controlled via sub-switches in `config/ssc_config.yaml`:

```yaml
reporting_module: "off"                # Master switch for Module 4 (default: off)

# Submodule 4A: Database Taxonomy Index Builder
make_taxonomy_index_submodule: "off"    # "on" or "off" (default: off)
viral_database_file: "HumanViral_Reference_02-07-2022.fa" # Database FASTA file (in ref_dir)
taxonomy_index_file: "{db_name}_taxonomy_index.tsv"        # Generated 14-column index TSV

# Submodule 4B: Master Sample Report Generator
generate_report_submodule: "off"       # "on" or "off" (default: off)
sample_metadata_file: "{dataset}_{genome_build}_master_coverage.tsv"
viral_rename_map_file: "{db_name}_map.tsv"                 # Generated by rename_fasta_contigs.sh
viral_bed_file: "{db_name}.bed"                            # Generated during reference sanitization
```

##### Submodule 4A Execution (Taxonomy Index Construction):
When `make_taxonomy_index_submodule: "on"` is set, running Snakemake compiles the 14-column taxonomy index for the whole viral database directly on the head node. **This submodule does not require a sample list file.**
```bash
snakemake --cores 1
```
*Comment / Note*: This step is to be created **only once** per reference database release as it indexes all accession IDs for the entire database used.

##### Submodule 4B Execution (Master Sample Reporting):
When `generate_report_submodule: "on"` is set, Snakemake compiles `ssc_reporting.job` and checks the `samples_list`. You can submit the sample batch execution to the cluster using:
```bash
batchRun -multibatch <SAMPLE_LIST> -config config/batch_jobexec_reporting.config -non-spot config/ssc_reporting.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>
```

#### Module 0: Human Genome Alignment Module (`ssc_human_alignment.job`)
Aligns paired-end FASTQ reads directly to the GRCh38 human reference genome using multithreaded BWA-MEM and streams output into coordinate-sorted, indexed CRAM files (`.cram` and `.cram.crai`).
```bash
# 1. Compile human alignment jobexec cluster files:
snakemake config/batch_jobexec_human_alignment.config config/ssc_human_alignment.job --cores 1

# 2. Submit parallel batch run across samples:
batchRun -multibatch <SAMPLE_LIST> -config config/batch_jobexec_human_alignment.config -non-spot config/ssc_human_alignment.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>
# Or via jobexec:
jobexec batch config/batch_jobexec_human_alignment.config config/ssc_human_alignment.job <SAMPLE_LIST>
```

#### Module 1: Coverage Calculation Module (`ssc_coverage.job`)

```bash
batchRun -multibatch <SAMPLE_LIST> -config config/batch_jobexec_coverage.config -non-spot config/ssc_coverage.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>
```
Combine all cloud-calculated sample coverage outputs into a master file named as `{dataset}_{genome_build}_master_coverage.tsv` and place it directly into your reference metadata directory at `/mnt/disks/staff/refs/`. Later specify this filename as `sample_metadata_file` in `config/ssc_config.yaml` to run **Module 4 (Reporting Module)** to calculate sample mean read depth and normalized viral copy numbers.

#### Module 2: Unmapped Extraction Module (`ssc_unmapped.job`)
```bash
batchRun -multibatch <SAMPLE_LIST> -config config/batch_jobexec_unmapped.config -non-spot config/ssc_unmapped.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>
```

#### Module 3: Alignment Module (`ssc_alignment.job`)
```bash
batchRun -multibatch <SAMPLE_LIST> -config config/batch_jobexec_vironator.config -non-spot config/ssc_alignment.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>
```

#### Module 4: Reporting Module (`ssc_reporting.job`)
```bash
batchRun -multibatch <SAMPLE_LIST> -config config/batch_jobexec_reporting.config -non-spot config/ssc_reporting.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>
```
- **Generated Per-Sample File**: `<sample_id>_viral_report.tsv`
- **Combined Master Output File**: Once all parallel cloud batch jobs finish, combine all individual sample reports into one master file named `{dataset}_{genome_build}_master_report.tsv` and place it in `/mnt/disks/staff/refs/` (used as input for Module 5 Stats, Module 6 Distributions, and Module 8 Refinement).

#### Module 5: Cohort Stats Summary Module (`stats_module`)
```bash
snakemake generate_stats --cores 1
```
Module 5 generates basic tabular statistics on the viruses that are found and summary statistics plots. Run this command on the head node after compiling your configuration. Snakemake reads your master report at `/mnt/disks/staff/refs/{dataset}_{genome_build}_master_report.tsv` and deposits summary files directly into `/mnt/disks/staff/{dataset}_{genome_build}_stats/`.
- **Generated Summary Files & Figures**:
  - `/mnt/disks/staff/<stats_out_dirname>/{dataset}_{genome_build}_stats_summary.tsv` (Tab-delimited 9-column cohort statistics snapshot)
  - `/mnt/disks/staff/<stats_out_dirname>/{dataset}_{genome_build}_stats_summary.md` (Markdown-formatted executive summary report)
  - `/mnt/disks/staff/<stats_out_dirname>/{dataset}_{genome_build}_virus_stats_{phase}_{project}_VIRUSES_percentages.tiff` (Horizontal bar plot showing positive sample percentages per virus)
  - `/mnt/disks/staff/<stats_out_dirname>/{dataset}_{genome_build}_virus_stats_{phase}_{project}_VIRUSES_reads.tiff` (Horizontal bar plot showing total assigned reads per virus)
  - `/mnt/disks/staff/<stats_out_dirname>/{dataset}_{genome_build}_virus_stats_{phase}_{project}_VIRUSES_classes.tiff` (Vertical bar plot showing occurrence counts across 4 preliminary classification categories, ordered highest to lowest)
  - `/mnt/disks/staff/<stats_out_dirname>/{dataset}_{genome_build}_virus_stats_{phase}_{project}_reads.tiff` (2-panel dataset-wide overall read count distribution figure)
  - `/mnt/disks/staff/<stats_out_dirname>/{dataset}_{genome_build}_virus_stats_{phase}_{project}_copy_number.tiff` (2-panel dataset-wide overall copy number distribution figure)
  - `/mnt/disks/staff/<stats_out_dirname>/{dataset}_{genome_build}_heatmap_read_counts.tiff` (Publication-quality 300 DPI hierarchical clustering heatmap for viral read counts)
  - `/mnt/disks/staff/<stats_out_dirname>/{dataset}_{genome_build}_heatmap_copy_number.tiff` (Publication-quality 300 DPI hierarchical clustering heatmap for viral copy number values)

#### Module 6: Distributions & Plots Module (`distributions_module`)
```bash
snakemake generate_distributions --cores 1
```
- **Generated Summary Files & Figures**:
  - `<stats_out_dirname>/{dataset}_{genome_build}_virus_stats_{phase}_{project}_VIRUSES.tsv` (Detailed per-virus summary table containing prevalence, sample counts, read counts, mean mapped reads, and preliminary classifications: `sporadic_noise`, `systematic_noise`, `virome`, `infection`).
  - `<stats_out_dirname>/{dataset}_{genome_build}_virus_stats_{phase}_{project}_OVERALL.tsv` (High-level executive dataset summary table reporting total positive samples, cohort size, overall positivity rate, total mapped reads, unique viruses detected, top prevalent virus, mean reads per positive sample, and overall classification).
  - `<plots_out_dirname>/VIRUS_<PHASE>_<PROJECT>_<STRATEGY>_<VIRUS_ACC>_<VIRUS_NAME>.tiff` (2-panel 300 DPI publication TIFF figure per detected virus species).
  - `<plots_out_dirname>/<PHASE>_<PROJECT>_<STRATEGY>_VIRAL_POSITIVITY_RATES.tiff` (Cohort-wide horizontal bar plot showing viral prevalence % across all detected viruses for the entire phase/project).
  - `<plots_out_dirname>/<PHASE>_<PROJECT>_<STRATEGY>_OVERALL_SUMMARY.tiff` (Executive 2-panel overall dataset summary figure for the entire phase/project).

#### Module 7: SAM Flag Comparison Module (`ssc_flag_comparison.job`)
```bash
batchRun -multibatch <SAMPLE_LIST> -config config/batch_jobexec_flag_comparison.config -non-spot config/ssc_flag_comparison.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>
```

#### Module 8: NCBI RefSeq Refinement Module (`ssc_refinement.job`)
```bash
# 1. Build taxonomy index TSV:
snakemake config/db_metadata/viral_reference_taxonomy_index.tsv --cores 1

# 2. Download RefSeq complete genomes & build date-stamped Bowtie2 combined reference:
snakemake config/ncbi_download.completed --cores 1

# 3. Compile refinement jobexec cluster files:
snakemake config/ssc_refinement.job --cores 1

# 4. Submit parallel batch run across target refinement samples:
batchRun -multibatch config/refinement_samples.tsv -config config/batch_jobexec_refinement.config -non-spot config/ssc_refinement.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>

# 5. Merge all cohort outputs into master report (run on head node after cluster jobs finish):
snakemake cohort_refinement_master.tsv --cores 1
```

#### Module 10: Discordant Reads & Viral Integration Breakpoints Module (`ssc_discordant.job`)
```bash
# 1. Compile discordant jobexec cluster files:
snakemake config/ssc_discordant.job --cores 1

# 2. Submit parallel batch run across samples:
jobexec batch config/batch_jobexec_resources.config config/ssc_discordant.job samples_list.txt
# Or with module-specific config:
jobexec batch config/batch_jobexec_discordant.config config/ssc_discordant.job samples_list.txt

# 3. Consolidate all per-sample reports into master cohort report (run on head node after cluster jobs finish):
snakemake cohort_breakpoints_master.tsv --cores 1
```

### 7. Systematic SAM Flag Comparison Module (`flag_comparison_module`)
Evaluates the clean strategy across 3 SAM flag filtering commands to isolate the impact of aligner `-f 2` flags vs manual bitwise flags:

> [!NOTE]
> **`clean_flags` Strategy Rationale:**
> The `clean_flags` strategy mitigates the overly restrictive read reduction applied by the `f2` (bitwise flag 2 / proper pair) filter in the alternative "no-flag" option. Relying strictly on proper pair alignment statistics can inadvertently discard valid candidate reads, particularly when alignment statistics or pairing metrics are unaligned or not entirely solid. Retaining reads via `clean_flags` ensures a more comprehensive and robust viral read capture.

> [!NOTE]
> This module is not estimating viral abundance. It is estimating high-confidence, reference-discriminating read-pair support. That is useful, however it will undercount viruses when related references share homologous sequence. This is for the reference-unique evidence. We report high-specificity reference-discriminating read support.

1. Enable Module 7 in `config/ssc_config.yaml`:
   ```yaml
   flag_comparison_module: "on"
   cohort_scope: "combined_all"
   ```
2. Compile job and resource config files on head node:
   ```bash
   snakemake --cores 1
   ```
3. Submit parallel batch run across all cohort samples (single batch across all phases & projects):
   ```bash
   batchRun -multibatch <SAMPLE_LIST> -config config/batch_jobexec_flag_comparison.config -non-spot config/ssc_flag_comparison.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>
   ```
4. Merge all cohort outputs into single master TSV report:
   ```bash
   snakemake cohort_flag_comparison_master.tsv --cores 1
   ```

#### 7-Column Flag Comparison Master TSV Schema (`cohort_flag_comparison_master.tsv`):
1. `phase` — Cohort phase (e.g. `phase1`)
2. `project` — Cohort project identifier (e.g. `base`)
3. `sample` — Sample ID
4. `strategy` — Strategy evaluated (`clean`)
5. `original_command_flags_reads` — Read pair count using original command WITH `-f 2`
6. `original_command_flags_modified_reads` — Read pair count using original command WITHOUT `-f 2`
7. `manual_flags_reads` — Read pair count using manual `gawk` bitwise flag verification

### 8. NCBI RefSeq Refinement Module (`refinement_module`)
Re-aligns candidate viral reads against complete RefSeq viral genomes for detected species groups using Bowtie2 multimapping (`-k 10`) and a Hierarchical Two-Tier Classification Cascade (`Reference_Unique` vs `Species_Supportive`):

1. Configure Module 8 in `config/ssc_config.yaml`:
   ```yaml
   refinement_module: "on"
   make_taxonomy_index_stage: "on"
   ncbi_download_stage: "on"
   build_combined_ref_stage: "on"
   refinement_alignment_stage: "on"
   ```
2. Build 14-column taxonomy index TSV:
   ```bash
   snakemake config/db_metadata/viral_reference_taxonomy_index.tsv --cores 1
   ```
3. Download complete RefSeq genomes & build Bowtie2 combined index:
   ```bash
   snakemake config/ncbi_download.completed --cores 1
   ```
4. Compile jobexec refinement cluster files:
   ```bash
   snakemake config/ssc_refinement.job --cores 1
   ```
5. Submit parallel batch run across target refinement samples:
   ```bash
   batchRun -multibatch config/refinement_samples.tsv -config config/batch_jobexec_refinement.config -non-spot config/ssc_refinement.job -investigator <INVESTIGATOR_TAG> -pau <PAU_CODE>
   ```

> [!IMPORTANT]
> **Post-Cluster Consolidation**: `batchRun` executes sample jobs in parallel across cluster nodes, writing per-sample TSVs under `/mnt/disks/staff/SSC_hg38_refinement/`. Once all parallel cluster jobs have finished, run Step 6 separately on the head node to aggregate all results into `cohort_refinement_master.tsv`:

6. Consolidate cohort output into master report (run on head node after cluster jobs finish):
   ```bash
   snakemake cohort_refinement_master.tsv --cores 1
   ```

#### 11-Column Refinement Master TSV Schema (`cohort_refinement_master.tsv`):
1. `Sample_ID` — Sample accession identifier
2. `Phase` — Cohort phase tag (e.g. `phase1`)
3. `Project` — Cohort project tag (or `base` if empty)
4. `Virus_Accession` — Best-supported RefSeq Accession ID
5. `Species_TaxID` — NCBI Species Taxonomy ID
6. `Species_Name` — Taxonomic Species Scientific Name
7. `Reference_Unique_Reads` — Read pairs uniquely supporting this specific reference
8. `Species_Supportive_Reads` — Read pairs supporting this species group (group-collapsed)
9. `Total_Refined_Reads` — Total read pairs (`Reference_Unique_Reads` + `Species_Supportive_Reads`)
10. `Refined_Copy_Number` — Reference-discriminating copy-number lower bound
11. `Classification_Status` — Call status (`CONFIRMED_UNIQUE` or `SPECIES_COLLAPSED`)

### 10. Discordant Reads & Viral Integration Breakpoints Module (`discordant_module`)
Extracts chimeric discordant read pairs and split reads from unfiltered alignment CRAM files (`exogeneSR_viral_clean.cram` or `exogeneSR_viral_raw.cram`), discovers host-viral genomic integration sites using a rightward sliding window clustering algorithm (Seed / Anchor approach), estimates single-nucleotide or interval integration boundaries, assesses call confidence, annotates nearest human genomic genes (GENCODE), and generates dual genomic arc visualizations.

1. Configure Module 10 in `config/ssc_config.yaml`:
   ```yaml
   discordant_module: "on"
   discordant_input_strategy: "exogeneSR_viral_clean.cram"
   integration_cluster_window: 500
   integration_min_support_reads: 2
   plot_integration_plots: "on"
   plot_min_support_reads: 3
   plot_top_loci_count: 5
   ```
2. Compile job and resource config files:
   ```bash
   snakemake config/ssc_discordant.job --cores 1
   ```
3. Submit parallel batch run across cohort samples:
   ```bash
   jobexec batch config/batch_jobexec_resources.config config/ssc_discordant.job samples_list.txt
   ```
4. Consolidate cohort output into master report (run on head node after cluster jobs finish):
   ```bash
   snakemake cohort_breakpoints_master.tsv --cores 1
   ```

#### Algorithmic Workflow & Methodology:
1. **Extraction of Discordant Pairs & Split Reads**:
   - Inspects the unfiltered CRAM file where inter-chromosomal chimeric pairs and soft-clipped reads are retained.
   - **Discordant Pairs**: Exactly one mate aligns to a human chromosome (`chr1`-`chr22`, `chrX`, `chrY`, `chrM`) and one mate aligns to a viral contig with `MAPQ >= 20` (unmapped, secondary, supplementary, and duplicate reads excluded).
   - **Split / Soft-Clipped Reads**: Primary alignment on human or viral contig with `MAPQ >= 10`, soft-clipped by $\ge 25$ bp, containing a cross-species supplementary alignment (`SA` tag) with `MAPQ >= 10`.
2. **Rightward Sliding Window Clustering (Seed / Anchor Approach)**:
   - Discordant pairs are sorted by human coordinate and grouped strictly by:
     - `(human_chr, viral_accession, structural_orientation)`
   - **Structural Orientation Classes**:
     - **Colinear**: $(+, -)$ and $(-, +)$ read pairs (opposite strands).
     - **Inverted**: $(+, +)$ and $(-,-)$ read pairs (same strands).
   - The first discordant pair in an orientation group is designated as the cluster **seed / anchor**.
   - As subsequent candidate pairs are scanned, the human window extends rightward from the anchor boundary (`human_max_end + window_size`, default 500 bp), while requiring the viral distance between reads to be $\le window\_size$. Split reads matching the genomic interval are merged into the cluster.
3. **Breakpoint Estimation**:
   - **Single-Nucleotide Precision ($N_{split} \ge 1$)**: When split reads are present, the exact junction point is calculated directly from the read coordinate where sequence transitions from human to viral genome.
   - **Interval Estimation ($N_{split} = 0$)**: When supported only by discordant pairs, the breakpoint is estimated from the median human coordinate and median viral coordinate of supporting reads.
4. **Call Confidence Rating**:
   - **High Confidence**: Supported by $\ge 1$ split read AND $\ge 3$ total supporting reads ($N_{split} \ge 1$ and $N_{total} \ge 3$).
   - **Medium Confidence**: ($N_{split} \ge 1$ and $N_{total} < 3$) OR ($N_{split} = 0$ and $N_{discordant} \ge 3$).
   - **Low Confidence**: $N_{split} = 0$ and $N_{discordant} < 3$.
5. **Human Gene Annotation**:
   - Automatically utilizes GENCODE v38 gene models (downloaded to `/mnt/disks/staff/refs/` if absent).
   - Identifies if the breakpoint falls within an intragenic element (gene name, biotype, strand) or is intergenic, recording the nearest upstream (`-X bp`) and downstream (`+Y bp`) genes.

#### Output Files & Locations:
- **Per-Sample Directory**: `gs://<output_bucket>/<vironator_out_dirname>/[phase]/[project]/<sample>/`
  - `exogeneSR_discordant_pairs.cram` & `.crai`: Filtered discordant pairs CRAM.
  - `exogeneSR_discordant_pairs.tsv`: Table of extracted discordant read pairs.
  - `exogeneSR_softclipped_split.cram` & `.crai`: Filtered split reads CRAM.
  - `exogeneSR_softclipped_split.tsv`: Table of extracted split/soft-clipped reads.
  - `<sample>_integration_breakpoints.tsv`: Locus-level summary of integration breakpoints.
  - `<sample>_locus_rank<N>_<chr>_<virus>_arc.tiff` (300 DPI publication TIFF) & `.png`: Dual genomic arc and coverage locus plot (for loci with $\ge 3$ supporting reads).
- **Consolidated Cohort Master Report**:
  - `cohort_breakpoints_master.tsv`: Generated on the head node after all parallel batch jobs complete via `snakemake cohort_breakpoints_master.tsv --cores 1`.

#### 20-Column Breakpoints TSV Schema (`<sample>_integration_breakpoints.tsv` and `cohort_breakpoints_master.tsv`):
1. `sample_id` — Sample identifier (e.g. `SS0012978`)
2. `phase` — Cohort phase tag (e.g. `phase1` or `unknown`)
3. `project` — Cohort project tag (e.g. `base`)
4. `human_chr` — Human host chromosome (e.g. `chr8`)
5. `human_coord_window` — Spanned human genomic interval (e.g. `chr8:128748200-128748950`)
6. `human_breakpoint_est` — Exact split junction or median human coordinate (e.g. `128748575`)
7. `viral_accession` — NCBI Viral Reference Accession ID (e.g. `NC_001357.1`)
8. `viral_name_sanitized` — Human-readable viral species name (e.g. `Human papillomavirus type 18`)
9. `viral_coord_window` — Spanned viral genomic interval (e.g. `NC_001357.1:3100-3450`)
10. `viral_breakpoint_est` — Exact split junction or median viral coordinate (e.g. `3275`)
11. `orientation` — Structural orientation (`colinear` or `inverted`)
12. `discordant_read_count` — Supporting discordant read pair count
13. `split_read_count` — Supporting chimeric split read count
14. `total_evidence_reads` — Total supporting reads (`discordant_read_count` + `split_read_count`)
15. `confidence` — Categorical confidence tier (`High`, `Medium`, or `Low`)
16. `intragenic_element` — Intragenic gene symbol if breakpoint falls within a gene body (or `NA`)
17. `upstream_element` — Nearest upstream gene symbol (e.g. `PVT1`)
18. `upstream_distance_bp` — Distance in bp to nearest upstream gene (e.g. `-1420`)
19. `downstream_element` — Nearest downstream gene symbol (e.g. `MYC`)
20. `downstream_distance_bp` — Distance in bp to nearest downstream gene (e.g. `+5230`)
