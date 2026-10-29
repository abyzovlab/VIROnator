import os

configfile: "config/ssc_config.yaml"

include: "rules/ssc.smk"

# Validate samples_list file only if batch job modules requiring samples are active
batch_modules_active = any([
    config.get("human_alignment_module", "off") == "on",
    config.get("unmapped_extraction", "off") == "on",
    config.get("viral_db_alignment", "off") == "on",
    (config.get("reporting_module", "off") == "on" and config.get("generate_report_submodule", "off") == "on"),
    config.get("coverage_module", "off") == "on",
    config.get("refinement_module", "off") == "on",
    config.get("discordant_module", "off") == "on"
])

if batch_modules_active:
    samples_list_path = config.get("samples_list", "")
    if samples_list_path and not os.path.isabs(samples_list_path):
        if os.path.exists(samples_list_path):
            samples_list_path = os.path.abspath(samples_list_path)
        else:
            work_dir = config.get("work_dir", ".")
            samples_list_path = os.path.join(work_dir, samples_list_path)

    if not samples_list_path or not os.path.exists(samples_list_path):
        raise FileNotFoundError(
            f"Samples list file not found at: {samples_list_path}. Please place it in your VIROnator directory."
        )

project_part = f"{config['project']}/" if config["project"] else ""

# Build targets dynamically based on config switches
targets = []

# Always target the resource configurations compilation
targets.extend([
    "config/batch_jobexec_resources.config",
    "config/batch_jobexec_unmapped.config",
    "config/batch_jobexec_vironator.config",
    "config/batch_jobexec_reporting.config",
    "config/batch_jobexec_coverage.config",
    "config/batch_jobexec_flag_comparison.config",
    "config/batch_jobexec_refinement.config",
    "config/batch_jobexec_discordant.config",
    "config/batch_jobexec_human_alignment.config",
    "config/ssc_coverage.job"
])

# Module 0: Human Genome Alignment (FASTQ -> CRAM via BWA-MEM)
if config.get("human_alignment_module", "off") == "on":
    targets.append("config/ssc_human_alignment.job")

# Module 1: Coverage calculation
if config.get("coverage_module", "off") == "on":
    targets.append("config/coverage_dir.created")

# Module 1: Unmapped extraction
if config.get("unmapped_extraction", "on") == "on":
    targets.append("config/unmapped_dir.created")
    targets.append("config/ssc_unmapped.job")

# Module 2: Viral DB alignment (VIROnator)
if config.get("viral_db_alignment", "off") == "on":
    targets.append("config/vironator_dir.created")
    targets.append("config/ssc_alignment.job")

if config.get("reporting_module", "off") == "on":
    if config.get("make_taxonomy_index_submodule", "off") == "on":
        tax_idx = os.path.join("config/db_metadata", f"{config.get('db_name', 'HumanViral_Reference_02-07-2022')}_taxonomy_index.tsv")
        targets.append(tax_idx)
    if config.get("generate_report_submodule", "off") == "on":
        targets.append("config/reports_dir.created")
        targets.append("config/ssc_reporting.job")

# Module 5: Stats Summary & Plots
if config.get("stats_module", "off") == "on":
    _ds = config.get("dataset", "DATASET")
    _gb = config.get("genome_build", config.get("build", "hg38"))
    targets.append(f"{_ds}_{_gb}_stats_summary.tsv")
    targets.append(f"{_ds}_{_gb}_stats_summary.md")

# Module 5: Distributions Plots & Virus Stats
if config.get("distributions_module", "off") == "on":
    targets.append("config/distributions.done")

# Module 7: Flag Comparison & IGV Snapshot Generation Module
if config.get("flag_comparison_module", "off") == "on":
    targets.append("config/flag_difference_crams.done")

# Module 8: Refinement Module
if config.get("refinement_module", "off") == "on":
    targets.append("config/ssc_refinement.job")

# Module 10: Discordant Reads & Viral Integration Breakpoints Module
if config.get("discordant_module", "off") == "on":
    targets.append("config/ssc_discordant.job")

rule all:
    input:
        targets


