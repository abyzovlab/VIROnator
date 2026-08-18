import os

configfile: "config/ssc_config.yaml"

include: "rules/ssc.smk"

# Validate that the samples list file exists in current execution directory or work_dir
samples_list_path = config["samples_list"]
if not os.path.isabs(samples_list_path):
    if os.path.exists(samples_list_path):
        samples_list_path = os.path.abspath(samples_list_path)
    else:
        work_dir = config.get("work_dir", ".")
        samples_list_path = os.path.join(work_dir, samples_list_path)

if not os.path.exists(samples_list_path):
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
    "config/ssc_coverage.job"
])

# Module 0: Coverage calculation
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

# Module 3: Reporting
if config.get("reporting_module", "off") == "on":
    targets.append("config/reports_dir.created")
    targets.append("config/ssc_reporting.job")

rule all:
    input:
        targets
