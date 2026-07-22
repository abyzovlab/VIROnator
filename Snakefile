import os

configfile: "config/ssc_config.yaml"

include: "rules/ssc.smk"

rule all:
    input:
        # Target 1: The GCS output directory placeholder
        os.path.join(
            config["staff_dir"],
            f"SSC_hg38_unmapped/phase{config['phase']}/{config['project']}/test.txt"
        ),
        # Target 2: The dynamically configured job file inside config/
        "config/ssc.job",
        # Target 3: The dynamically configured resource config inside config/
        "config/batch_jobexec_resources.config"
