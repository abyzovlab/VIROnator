# Logic to prepare output directories and dynamically configure the job and resource files

import os

rule create_output_directory:
    """
    Creates the target output directory and copies the placeholder file
    into it to ensure directory persistence on GCS.
    """
    input:
        placeholder=config["placeholder_file"]
    output:
        placeholder=os.path.join(config["staff_dir"], "SSC_hg38_unmapped/phase{phase}/{project}/test.txt")
    shell:
        """
        out_dir=$(dirname {output.placeholder})
        mkdir -p "$out_dir"
        cp {input.placeholder} {output.placeholder}
        """

rule generate_job_file:
    """
    Generates the final ssc.job file from ssc.job.template inside the config directory
    by substituting variables defined in the Snakemake configuration.
    """
    input:
        template="config/ssc.job.template"
    output:
        job="config/ssc.job"
    run:
        with open(input.template, "r") as f:
            content = f.read()
        
        # Safely replace placeholders without affecting bash variable references (e.g. ${PHASE})
        formatted_content = (
            content.replace("{phase}", str(config["phase"]))
            .replace("{project}", str(config["project"]))
            .replace("{ref_genome}", str(config["ref_genome"]))
            .replace("{lab_dir}", str(config["lab_dir"]))
            .replace("{staff_dir}", str(config["staff_dir"]))
            .replace("{staff_bucket}", str(config["staff_bucket"]))
        )
        
        with open(output.job, "w") as f:
            f.write(formatted_content)

rule generate_resources_config:
    """
    Generates the final batch_jobexec_resources.config from its template inside the config directory
    by substituting the GCS bucket variables defined in the Snakemake configuration.
    """
    input:
        template="config/batch_jobexec_resources.config.template"
    output:
        config="config/batch_jobexec_resources.config"
    run:
        with open(input.template, "r") as f:
            content = f.read()
        
        formatted_content = (
            content.replace("{staff_bucket}", str(config["staff_bucket"]))
            .replace("{lab_bucket}", str(config["lab_bucket"]))
        )
        
        with open(output.config, "w") as f:
            f.write(formatted_content)
