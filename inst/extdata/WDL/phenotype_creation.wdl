version 1.0

task phecode_generation {
    
    input {
        File health_data
        File sex_info
        File? control_exclusions
    }

    Int cores = 16
    String phecode_save_file = "phecodes.RDS"
    String range_ID_save_file = "range_ID.RDS"

    command <<<
       Rscript `Rscript -e 'cat(system.file("extdata/scripts/phenotype_generation","03a_phecode_generation.R", package = "DeepPheWAS"))'` \
            --health_data ~{health_data} \
            --sex_info ~{sex_info} \
            --N_cores ~{cores} \
            --phecode_save_file ~{phecode_save_file} \
            --range_ID_save_file ~{range_ID_save_file}
    >>>

    output {
        File phecode_file = phecode_save_file
        File? range_ID_file = range_ID_save_file
        Array[File]+ out = glob("*.RDS")
    }

	runtime {
		cpu: cores
		memory: "120 GB"
	}
}

task data_field_phenotypes {

    input {
        File min_data
        File? phewas_manifest
    }

    String phenotype_save_file = "data_field_phenotypes.RDS"

    command <<<
       Rscript `Rscript -e 'cat(system.file("extdata/scripts/phenotype_generation","03b_data_field_phenotypes.R", package = "DeepPheWAS"))'` \
            --min_data ~{min_data} \
            --phenotype_save_file ~{phenotype_save_file} \
            ~{"--PheWAS_manifest_overide " + phewas_manifest}
    >>>

    output {
        File out = phenotype_save_file
    }

    runtime {
        cpu: 1
        memory: "200 GB"
    }
}

task creating_concepts {

    input {
        File GPP
        File health_data
        File? phewas_manifest
        File? code_list
    }

    String concept_save_file = "concepts.RDS"
    String all_dates_save_file = "all_dates.RDS"

    command <<<
		Rscript `Rscript -e 'cat(system.file("extdata/scripts/phenotype_generation","03c_creating_concepts.R", package = "DeepPheWAS"))'` \
			--GPP ~{GPP} \
			--concept_save_file ~{concept_save_file} \
			--all_dates_save_file ~{all_dates_save_file} \
			--health_data ~{health_data} \
			~{"--PheWAS_manifest_overide " + phewas_manifest } \
			~{"--code_list_folder_override " + code_list}
    >>>

    output {
        File concept_file = concept_save_file
        File all_dates_file = all_dates_save_file
        Array[File]+ out = glob("*.RDS")
    }

    runtime {
        cpu: 1
        memory: "120 GB"
    }
}

task primary_care_quantitative_phenotypes {

    input {
        File GPC
        File DOB
        File? phewas_manifest
        File? code_list
    }

    String phenotype_save_file = "PQP.RDS"
    Int cores = 16

    command <<<
       Rscript `Rscript -e 'cat(system.file("extdata/scripts/phenotype_generation","03d_primary_care_quantitative_phenotypes.R", package = "DeepPheWAS"))'` \
        --GPC ~{GPC} \
        --DOB ~{DOB} \
        --phenotype_save_file ~{phenotype_save_file} \
        ~{"--PheWAS_manifest_overide " + phewas_manifest} \
        ~{"--code_list_overide " + code_list} \
        --N_cores ~{cores}
    >>>

    output {
        File out = phenotype_save_file
    }

    runtime {
        cpu: cores
        memory: "120 GB"
    }
}

task formula_phenotypes {

    input {
        File min_data
        File? data_field_phenotypes
        File? sex_info
        File? phewas_manifest
    }

    String phenotype_save_file = "formula_phenotypes.RDS"

    command <<<
       Rscript `Rscript -e 'cat(system.file("extdata/scripts/phenotype_generation","04_formula_phenotypes.R", package = "DeepPheWAS"))'` \
            --min_data ~{min_data} \
            --phenotype_save_file ~{phenotype_save_file} \
            ~{"--data_field_phenotypes " + data_field_phenotypes} \
            ~{"--sex_info " + sex_info} \
            ~{"--PheWAS_manifest_overide " + phewas_manifest}
    >>>

    output {
        File out = phenotype_save_file
    }

    runtime {
        cpu: 1
        memory: "120 GB"
    }
}

task composite_phenotypes {

    input {
        File? phecode_file
        File? range_ID_file
        File data_field_file
        File? concept_file
        File? all_dates_file
        File? PQP_file
        File? formula_file
        File? control_populations
        File? composite_phenotype_map_overide
    }

    String phenotype_save_file = "composite_phenotypes.RDS"

    command <<<
       Rscript `Rscript -e 'cat(system.file("extdata/scripts/phenotype_generation","05_composite_phenotypes.R", package = "DeepPheWAS"))'` \
       --phenotype_files ~{data_field_file}~{"," + phecode_file}~{"," + range_ID_file}~{"," + concept_file}~{"," + all_dates_file}~{"," + PQP_file}~{"," + formula_file} \
       --phenotype_save_file composite_phenotypes.RDS \
       ~{"--control_populations " + control_populations} \
       ~{"--composite_phenotype_map_overide " + composite_phenotype_map_overide}
    >>>

    output {
        File out = phenotype_save_file
    }

    runtime {
        cpu: 1
        memory: "120 GB"
    }
}
