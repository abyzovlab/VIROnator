# Logic to prepare output directories and dynamically configure the job and resource files

import os

project_part = f"{config['project']}/" if config["project"] else ""

rule create_output_directory:
    """
    Initializes the target unmapped output directory on GCS using gsutil if not already initialized.
    """
    input:
        placeholder=config["placeholder_file"],
        config_file="config/ssc_config.yaml"
    output:
        token="config/unmapped_dir.created"
    shell:
        """
        TARGET_PATH="gs://{config[output_bucket]}/{config[unmapped_out_dirname]}/phase{config[phase]}/{project_part}test.txt"
        if ! gsutil -q stat "$TARGET_PATH"; then
            gsutil cp {input.placeholder} "$TARGET_PATH"
        fi
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
            .replace("{ref_cram_decoder}", os.path.join(config["ref_dir"], config.get("ref_cram_decoder", "GRCh38_full_analysis_set_plus_decoy_hla.fa")))
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
    Generates module-specific batch_jobexec resource config files from template.
    """
    input:
        template="config/batch_jobexec_resources.config.template",
        config_file="config/ssc_config.yaml"
    output:
        resources_config="config/batch_jobexec_resources.config",
        unmapped_config="config/batch_jobexec_unmapped.config",
        vironator_config="config/batch_jobexec_vironator.config",
        reporting_config="config/batch_jobexec_reporting.config",
        coverage_config="config/batch_jobexec_coverage.config"
    run:
        with open(input.template, "r") as f:
            template_content = f.read()
        
        base_sub = (
            template_content.replace("{output_bucket}", str(config["output_bucket"]))
            .replace("{data_bucket}", str(config["data_bucket"]))
            .replace("{gcp_account}", str(config["gcp_account"]))
            .replace("{pi}", str(config["pi"]))
            .replace("{pau}", str(config["pau"]))
            .replace("{task_name}", str(config["task_name"]))
        )
        
        # 1. Unmapped config
        unmapped_content = base_sub.replace("{jobexec_dirname}", str(config.get("unmapped_jobexec_dirname", "jobexec_unmapped")))
        with open(output.unmapped_config, "w") as f:
            f.write(unmapped_content)

        # 2. VIROnator config
        vironator_content = base_sub.replace("{jobexec_dirname}", str(config.get("vironator_jobexec_dirname", "jobexec_vironator")))
        with open(output.vironator_config, "w") as f:
            f.write(vironator_content)

        # 3. Reporting config
        reporting_content = base_sub.replace("{jobexec_dirname}", str(config.get("reporting_jobexec_dirname", "jobexec_reporting")))
        with open(output.reporting_config, "w") as f:
            f.write(reporting_content)

        # 4. Coverage config
        coverage_content = base_sub.replace("{jobexec_dirname}", "jobexec_coverage")
        with open(output.coverage_config, "w") as f:
            f.write(coverage_content)

        # Active default config based on active module switch
        if config.get("coverage_module", "off") == "on":
            active_content = coverage_content
        elif config.get("reporting_module", "off") == "on":
            active_content = reporting_content
        elif config.get("viral_db_alignment", "off") == "on":
            active_content = vironator_content
        else:
            active_content = unmapped_content
            
        with open(output.resources_config, "w") as f:
            f.write(active_content)

rule create_vironator_directory:
    """
    Initializes the target vironator output directory on GCS using gsutil if not already initialized.
    """
    input:
        placeholder=config["placeholder_file"],
        config_file="config/ssc_config.yaml"
    output:
        token="config/vironator_dir.created"
    shell:
        """
        TARGET_PATH="gs://{config[output_bucket]}/{config[vironator_out_dirname]}/phase{config[phase]}/{project_part}test.txt"
        if ! gsutil -q stat "$TARGET_PATH"; then
            gsutil cp {input.placeholder} "$TARGET_PATH"
        fi
        touch {output.token}
        """

rule generate_align_job_file:
    """
    Generates the final ssc_alignment.job file from ssc_alignment.job.template inside the config directory
    by substituting variables defined in the Snakemake configuration.
    """
    input:
        template="config/ssc_alignment.job.template",
        config_file="config/ssc_config.yaml"
    output:
        job="config/ssc_alignment.job"
    run:
        import shutil, subprocess
        # Auto-sync repository scripts and db_metadata to shared output_dir mount (/mnt/disks/staff)
        staff_scripts = os.path.join(config["output_dir"], "scripts")
        staff_metadata = os.path.join(config["output_dir"], "db_metadata")
        try:
            os.makedirs(staff_scripts, exist_ok=True)
            os.makedirs(staff_metadata, exist_ok=True)
            if os.path.exists("scripts"):
                shutil.copytree("scripts", staff_scripts, dirs_exist_ok=True)
            if os.path.exists("config/db_metadata"):
                shutil.copytree("config/db_metadata", staff_metadata, dirs_exist_ok=True)
        except Exception:
            # Fallback to gsutil if local /mnt/disks mount is restricted on head node
            bucket = config.get("output_bucket")
            if bucket:
                subprocess.run(f"gsutil -q cp -r scripts/* gs://{bucket}/scripts/ 2>/dev/null", shell=True)
                subprocess.run(f"gsutil -q cp -r config/db_metadata/* gs://{bucket}/db_metadata/ 2>/dev/null", shell=True)

        with open(input.template, "r") as f:
            content = f.read()
        
        # Safely replace placeholders without affecting bash variable references (e.g. ${PHASE})
        formatted_content = (
            content.replace("{phase}", str(config["phase"]))
            .replace("{project}", str(config["project"]))
            .replace("{output_bucket}", str(config["output_bucket"]))
            .replace("{output_dir}", str(config["output_dir"]))
            .replace("{ref_cram_decoder}", os.path.join(config["ref_dir"], config.get("ref_cram_decoder", "GRCh38_full_analysis_set_plus_decoy_hla.fa")))
            .replace("{ref_human_full}", os.path.join(config["ref_dir"], config["ref_human_full"]))
            .replace("{ref_human_no_ebv}", os.path.join(config["ref_dir"], config["ref_human_no_ebv"]))
            .replace("{ref_viral_base}", os.path.join(config["ref_dir"], config.get("ref_viral_base", "HumanViral_Reference_02-07-2022.fa")))
            .replace("{ref_plasmids}", os.path.join(config["ref_dir"], config.get("ref_plasmids", "SnapGene_plasmids_modified.fa")))
            .replace("{ref_mouse}", os.path.join(config["ref_dir"], config.get("ref_mouse", "mm39_ms_modified.fa")))
            .replace("{ref_vir_cont}", os.path.join(config["ref_dir"], config.get("ref_vir_cont", config.get("ref_viral", "HumanViral_Reference_02-07-2022_SnapGene_plasmids_modified_mm39_ms_modified.fa"))))
            .replace("{ref_combined_default}", os.path.join(config["ref_dir"], config["ref_combined_default"]))
            .replace("{ref_decoy}", os.path.join(config["ref_dir"], config["ref_decoy"]))
            .replace("{ref_rna}", os.path.join(config["ref_dir"], config["ref_rna"]))
            .replace("{viral_contigs_path}", os.path.join(config["ref_dir"], config["viral_contigs_file"]))
            .replace("{viral_bed_path}", os.path.join(config["ref_dir"], config["viral_bed_file"]))
            .replace("{plasmid_contigs_path}", os.path.join(config["ref_dir"], config.get("plasmid_contigs_file", "SnapGene_modified.contigs.txt")))
            .replace("{mouse_contigs_path}", os.path.join(config["ref_dir"], config.get("mouse_contigs_file", "mm39_modified.contigs.txt")))
            .replace("{contamination_quantification}", str(config.get("contamination_quantification", "off")))
            .replace("{bwa_bin}", str(config["bwa_bin"]))
            .replace("{python_bin}", str(config.get("python_bin", "python3")))
            .replace("{exogene_python_dir}", str(config["scripts_dir"]))
            .replace("{init_ref_script_path}", os.path.join(config["scripts_dir"], config["init_ref_script"]))
            .replace("{combined_refs_dir_path}", os.path.normpath(os.path.join(config["ref_dir"], config.get("combined_refs_dir", ""))))
            .replace("{align_and_filter}", str(config.get("align_and_filter", "on")))
            .replace("{align_only}", str(config.get("align_only", "off")))
            .replace("{filter_only}", str(config.get("filter_only", "off")))
            .replace("{dual_mode}", str(config.get("dual_mode", "on")))
            .replace("{unmapped_out_dirname}", str(config["unmapped_out_dirname"]))
            .replace("{vironator_out_dirname}", str(config["vironator_out_dirname"]))
            .replace("{unmapped_jobexec_dirname}", str(config["unmapped_jobexec_dirname"]))
            .replace("{vironator_jobexec_dirname}", str(config["vironator_jobexec_dirname"]))
        )
        
        with open(output.job, "w") as f:
            f.write(formatted_content)

rule create_reports_directory:
    """
    Initializes the target reports output directory on GCS using gsutil if not already initialized.
    """
    input:
        placeholder=config["placeholder_file"],
        config_file="config/ssc_config.yaml"
    output:
        token="config/reports_dir.created"
    shell:
        """
        TARGET_PATH="gs://{config[output_bucket]}/{config[reports_out_dirname]}/phase{config[phase]}/{project_part}test.txt"
        if ! gsutil -q stat "$TARGET_PATH"; then
            gsutil cp {input.placeholder} "$TARGET_PATH"
        fi
        touch {output.token}
        """

rule create_coverage_directory:
    """
    Initializes the target coverage output directory on GCS using gsutil if not already initialized.
    """
    input:
        placeholder=config["placeholder_file"],
        config_file="config/ssc_config.yaml"
    output:
        token="config/coverage_dir.created"
    shell:
        """
        TARGET_PATH="gs://{config[output_bucket]}/{config[coverage_out_dirname]}/phase{config[phase]}/{project_part}test.txt"
        if ! gsutil -q stat "$TARGET_PATH"; then
            gsutil cp {input.placeholder} "$TARGET_PATH"
        fi
        touch {output.token}
        """

rule generate_reporting_job_file:
    """
    Generates the final ssc_reporting.job file from ssc_reporting.job.template inside the config directory
    by substituting variables defined in the Snakemake configuration.
    """
    input:
        template="config/ssc_reporting.job.template",
        config_file="config/ssc_config.yaml"
    output:
        job="config/ssc_reporting.job"
    run:
        import shutil, subprocess
        # Auto-sync repository scripts and db_metadata to shared output_dir mount (/mnt/disks/staff)
        staff_scripts = os.path.join(config["output_dir"], "scripts")
        staff_metadata = os.path.join(config["output_dir"], "db_metadata")
        try:
            os.makedirs(staff_scripts, exist_ok=True)
            os.makedirs(staff_metadata, exist_ok=True)
            if os.path.exists("scripts"):
                shutil.copytree("scripts", staff_scripts, dirs_exist_ok=True)
            if os.path.exists("config/db_metadata"):
                shutil.copytree("config/db_metadata", staff_metadata, dirs_exist_ok=True)
        except Exception:
            bucket = config.get("output_bucket")
            if bucket:
                subprocess.run(f"gsutil -q cp -r scripts/* gs://{bucket}/scripts/ 2>/dev/null", shell=True)
                subprocess.run(f"gsutil -q cp -r config/db_metadata/* gs://{bucket}/db_metadata/ 2>/dev/null", shell=True)

        with open(input.template, "r") as f:
            content = f.read()
        
        formatted_content = (
            content.replace("{phase}", str(config["phase"]))
            .replace("{project}", str(config["project"]))
            .replace("{output_bucket}", str(config["output_bucket"]))
            .replace("{output_dir}", str(config["output_dir"]))
            .replace("{ref_human_full}", os.path.join(config["ref_dir"], config["ref_human_full"]))
            .replace("{ref_human_no_ebv}", os.path.join(config["ref_dir"], config["ref_human_no_ebv"]))
            .replace("{ref_vir_cont}", os.path.join(config["ref_dir"], config.get("ref_vir_cont", config.get("ref_viral", "HumanViral_Reference_02-07-2022_SnapGene_plasmids_modified_mm39_ms_modified.fa"))))
            .replace("{viral_rename_map_path}", os.path.join(config["ref_dir"], config.get("viral_rename_map_file", "HumanViral_Reference_02-07-2022_modified.rename_map.tsv")))
            .replace("{viral_bed_path}", os.path.join(config["ref_dir"], config["viral_bed_file"]))
            .replace("{sample_metadata_path}", os.path.join(config["ref_dir"], config.get("sample_metadata_file", "SSC_sample_metadata.tsv")))
            .replace("{bwa_bin}", str(config["bwa_bin"]))
            .replace("{python_bin}", str(config.get("python_bin", "python3")))
            .replace("{report_script_path}", os.path.join(config["scripts_dir"], config.get("report_script", "generate_report.py")))
            .replace("{combined_refs_dir_path}", os.path.normpath(os.path.join(config["ref_dir"], config.get("combined_refs_dir", ""))))
            .replace("{vironator_out_dirname}", str(config["vironator_out_dirname"]))
            .replace("{reports_out_dirname}", str(config["reports_out_dirname"]))
            .replace("{vironator_jobexec_dirname}", str(config["vironator_jobexec_dirname"]))
            .replace("{reporting_jobexec_dirname}", str(config.get("reporting_jobexec_dirname", "jobexec_reporting")))
        )
        
        with open(output.job, "w") as f:
            f.write(formatted_content)

rule generate_coverage_job_file:
    """
    Generates ssc_coverage.job from ssc_coverage.job.template.
    """
    input:
        template="config/ssc_coverage.job.template",
        config_file="config/ssc_config.yaml"
    output:
        job="config/ssc_coverage.job"
    run:
        with open(input.template, "r") as f:
            content = f.read()
        
        formatted_content = (
            content.replace("{phase}", str(config["phase"]))
            .replace("{project}", str(config["project"]))
            .replace("{output_bucket}", str(config["output_bucket"]))
            .replace("{output_dir}", str(config["output_dir"]))
            .replace("{data_dir}", str(config["data_dir"]))
            .replace("{sample_metadata_path}", os.path.join(config["ref_dir"], config.get("sample_metadata_file", "SSC_sample_metadata.tsv")))
        )
        
        with open(output.job, "w") as f:
            f.write(formatted_content)

rule generate_stats:
    """
    Generates cohort_stats_summary.tsv and cohort_stats_summary.md from master report.
    """
    input:
        script=os.path.join(config["scripts_dir"], config.get("stats_script", "generate_stats.py")),
        master_report=config.get("master_report_file", "cohort_master_viral_report_phase2_base.tsv"),
        config_file="config/ssc_config.yaml"
    output:
        tsv="cohort_stats_summary.tsv",
        md="cohort_stats_summary.md"
    shell:
        """
        python3 {input.script} --input-report {input.master_report} --out-dir .
        """
