# Logic to prepare output directories and dynamically configure the job and resource files

import os

project_part = f"{config['project']}/" if config["project"] else ""

rule create_output_directory:
    """
    Initializes the target unmapped output directory on GCS using gsutil.
    """
    input:
        placeholder=config["placeholder_file"]
    output:
        token="config/unmapped_dir.created"
    shell:
        """
        gsutil cp {input.placeholder} gs://{config[staff_bucket]}/{config[unmapped_out_dirname]}/phase{config[phase]}/{project_part}test.txt
        touch {output.token}
        """

rule generate_job_file:
    """
    Generates the final ssc_unmapped.job file from ssc_unmapped.job.template inside the config directory
    by substituting variables defined in the Snakemake configuration.
    """
    input:
        template="config/ssc_unmapped.job.template"
    output:
        job="config/ssc_unmapped.job"
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
            .replace("{unmapped_out_dirname}", str(config["unmapped_out_dirname"]))
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

rule create_vironator_directory:
    """
    Initializes the target vironator output directory on GCS using gsutil.
    """
    input:
        placeholder=config["placeholder_file"]
    output:
        token="config/vironator_dir.created"
    shell:
        """
        gsutil cp {input.placeholder} gs://{config[staff_bucket]}/{config[vironator_out_dirname]}/phase{config[phase]}/{project_part}test.txt
        touch {output.token}
        """

rule generate_align_job_file:
    """
    Generates the final ssc_align.job file from ssc_align.job.template inside the config directory
    by substituting variables defined in the Snakemake configuration.
    """
    input:
        template="config/ssc_align.job.template"
    output:
        job="config/ssc_align.job"
    run:
        with open(input.template, "r") as f:
            content = f.read()
        
        # Safely replace placeholders without affecting bash variable references (e.g. ${PHASE})
        formatted_content = (
            content.replace("{phase}", str(config["phase"]))
            .replace("{project}", str(config["project"]))
            .replace("{staff_bucket}", str(config["staff_bucket"]))
            .replace("{staff_dir}", str(config["staff_dir"]))
            .replace("{ref_human_full}", os.path.join(config["ref_dir"], config["ref_human_full"]))
            .replace("{ref_human_no_ebv}", os.path.join(config["ref_dir"], config["ref_human_no_ebv"]))
            .replace("{ref_viral}", os.path.join(config["ref_dir"], config["ref_viral"]))
            .replace("{ref_combined_default}", os.path.join(config["ref_dir"], config["ref_combined_default"]))
            .replace("{ref_decoy}", os.path.join(config["ref_dir"], config["ref_decoy"]))
            .replace("{ref_rna}", os.path.join(config["ref_dir"], config["ref_rna"]))
            .replace("{viral_contigs_path}", os.path.join(config["db_metadata_dir"], config["viral_contigs_file"]))
            .replace("{viral_bed_path}", os.path.join(config["db_metadata_dir"], config["viral_bed_file"]))
            .replace("{bwa_bin}", str(config["bwa_bin"]))
            .replace("{exogene_python_dir}", str(config["scripts_dir"]))
            .replace("{init_ref_script_path}", os.path.join(config["scripts_dir"], config["init_ref_script"]))
            .replace("{combined_refs_dir_path}", os.path.join(config["ref_dir"], config["combined_refs_dir"]))
            .replace("{alignment_mode}", str(config["alignment_mode"]))
            .replace("{unmapped_out_dirname}", str(config["unmapped_out_dirname"]))
            .replace("{vironator_out_dirname}", str(config["vironator_out_dirname"]))
        )
        
        with open(output.job, "w") as f:
            f.write(formatted_content)
