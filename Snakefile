import os

configfile: "config/ssc_config.yaml"

include: "rules/ssc.smk"

# Validate that the samples list file exists in the working directory
work_dir = config.get("work_dir", ".")
samples_list_path = config["samples_list"]
if not os.path.isabs(samples_list_path):
    samples_list_path = os.path.join(work_dir, samples_list_path)

if not os.path.exists(samples_list_path):
    raise FileNotFoundError(
        f"Samples list file not found at: {samples_list_path}. Please place it in your working directory."
    )

# Build targets dynamically based on config switches
targets = []

# Always target the resource configurations compilation
targets.append("config/batch_jobexec_resources.config")

# Module 1: Unmapped extraction
if config.get("unmapped_extraction", "on") == "on":
    targets.append(
        os.path.join(
            config["staff_dir"],
            f"SSC_hg38_unmapped/phase{config['phase']}/{config['project']}/test.txt"
        )
    )
    targets.append("config/ssc_unmapped.job")

# Module 2: Viral DB alignment (VIROnator)
if config.get("viral_db_alignment", "on") == "on":
    targets.append(
        os.path.join(
            config["staff_dir"],
            f"SSC_hg38_vironator/phase{config['phase']}/{config['project']}/test.txt"
        )
    )
    targets.append("config/ssc_align.job")

rule all:
    input:
        targets
