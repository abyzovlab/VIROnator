# Logic to prepare output directories and dynamically configure the job and resource files

import os

project_part = f"{config['project']}/" if config["project"] else ""

rule create_output_directory:
    """
    Initializes the target unmapped output directory on GCS using gsutil.
    """
    input:
        placeholder=config["placeholder_file"],
        config_file="config/ssc_config.yaml"
    output:
        token="config/unmapped_dir.created"
    shell:
        """
        gsutil cp {input.placeholder} gs://{config[output_bucket]}/{config[unmapped_out_dirname]}/phase{config[phase]}/{project_part}test.txt
        touch {output.token}
        """

rule generate_job_file:
    """
    Generates the final ssc_unmapped.job file from ssc_unmapped.job.template inside the config directory
    by substituting variables defined in the Snakemake configuration.
    """
    input:
        template="config/ssc_unmapped.job.template",
        config_file="config/ssc_config.yaml"
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
            .replace("{data_dir}", str(config["data_dir"]))
            .replace("{output_dir}", str(config["output_dir"]))
            .replace("{output_bucket}", str(config["output_bucket"]))
            .replace("{unmapped_out_dirname}", str(config["unmapped_out_dirname"]))
            .replace("{unmapped_jobexec_dirname}", str(config["unmapped_jobexec_dirname"]))
        )
        
        with open(output.job, "w") as f:
            f.write(formatted_content)

rule generate_resources_config:
    """
    Generates the final batch_jobexec_resources.config from its template inside the config directory
    by substituting the GCS bucket variables defined in the Snakemake configuration.
    """
    input:
        template="config/batch_jobexec_resources.config.template",
        config_file="config/ssc_config.yaml"
    output:
        resources_config="config/batch_jobexec_resources.config"
    run:
        with open(input.template, "r") as f:
            content = f.read()
        
        formatted_content = (
            content.replace("{output_bucket}", str(config["output_bucket"]))
            .replace("{data_bucket}", str(config["data_bucket"]))
            .replace("{gcp_account}", str(config["gcp_account"]))
            .replace("{pi}", str(config["pi"]))
            .replace("{pau}", str(config["pau"]))
            .replace("{task_name}", str(config["task_name"]))
            .replace("{unmapped_jobexec_dirname}", str(config["unmapped_jobexec_dirname"]))
        )
        
        with open(output.resources_config, "w") as f:
            f.write(formatted_content)

rule create_vironator_directory:
    """
    Initializes the target vironator output directory on GCS using gsutil.
    """
    input:
        placeholder=config["placeholder_file"],
        config_file="config/ssc_config.yaml"
    output:
        token="config/vironator_dir.created"
    shell:
        """
        gsutil cp {input.placeholder} gs://{config[output_bucket]}/{config[vironator_out_dirname]}/phase{config[phase]}/{project_part}test.txt
        touch {output.token}
        """

rule generate_align_job_file:
    """
    Generates the final ssc_align.job file from ssc_alignment.job.template inside the config directory
    by substituting variables defined in the Snakemake configuration.
    """
    input:
        template="config/ssc_alignment.job.template",
        config_file="config/ssc_config.yaml"
    output:
        job="config/ssc_align.job"
    run:
        import shutil
        # Auto-sync repository scripts and db_metadata to shared output_dir mount (/mnt/disks/staff)
        staff_scripts = os.path.join(config["output_dir"], "scripts")
        staff_metadata = os.path.join(config["output_dir"], "db_metadata")
        os.makedirs(staff_scripts, exist_ok=True)
        os.makedirs(staff_metadata, exist_ok=True)
        if os.path.exists("scripts"):
            shutil.copytree("scripts", staff_scripts, dirs_exist_ok=True)
        if os.path.exists("config/db_metadata"):
            shutil.copytree("config/db_metadata", staff_metadata, dirs_exist_ok=True)

        with open(input.template, "r") as f:
            content = f.read()
        
        # Safely replace placeholders without affecting bash variable references (e.g. ${PHASE})
        formatted_content = (
            content.replace("{phase}", str(config["phase"]))
            .replace("{project}", str(config["project"]))
            .replace("{output_bucket}", str(config["output_bucket"]))
            .replace("{output_dir}", str(config["output_dir"]))
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
            .replace("{unmapped_jobexec_dirname}", str(config["unmapped_jobexec_dirname"]))
            .replace("{vironator_jobexec_dirname}", str(config["vironator_jobexec_dirname"]))
        )
        
        with open(output.job, "w") as f:
            f.write(formatted_content)
