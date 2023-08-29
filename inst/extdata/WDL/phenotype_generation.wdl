version 1.0

import "data_wrangling.wdl"
import "phenotype_creation.wdl"

workflow phenotype_generation {

	input {
        File? phewas_manifest = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/PheWAS_manifest.csv"
		File? code_list
	}

	call data_wrangling.minimum_data {
	}
	
	call data_wrangling.data_preparation {
		input:
			min_data = minimum_data.out,
	}

	call phenotype_creation.phecode_generation {
		input:
			health_data = data_preparation.health_records,
			sex_info = data_preparation.combined_sex,
			control_exclusions = data_preparation.control_exclusions
	}

    call phenotype_creation.data_field_phenotypes {
        input:
            min_data = minimum_data.out,
            phewas_manifest = phewas_manifest
	}

    call phenotype_creation.creating_concepts {
        input:
            GPP = data_preparation.GP_P,
            health_data = data_preparation.health_records,
            phewas_manifest = phewas_manifest,
            code_list = code_list
    }

    call phenotype_creation.primary_care_quantitative_phenotypes {
        input:
            GPC = data_preparation.GP_C,
            DOB = data_preparation.DOB,
            phewas_manifest = phewas_manifest,
            code_list = code_list
    }

    call phenotype_creation.formula_phenotypes {
        input:
            min_data = minimum_data.out,
            data_field_phenotypes = data_field_phenotypes.out,
            sex_info = data_preparation.combined_sex,
            phewas_manifest = phewas_manifest
    }

    call phenotype_creation.composite_phenotypes {
        input:
            phecode_file = phecode_generation.phecode_file,
			range_ID_file = phecode_generation.range_ID_file,
			concept_file = creating_concepts.concept_file,
			all_dates_file = creating_concepts.all_dates_file,
			data_field_file = data_field_phenotypes.out,
			PQP_file = primary_care_quantitative_phenotypes.out,
			formula_file = formula_phenotypes.out,
			control_populations = data_preparation.control_populations
	}

	output {
		File min_data = minimum_data.out
		Array[File]+ out = data_preparation.out
		Array[File]+ phecode = phecode_generation.out
		File data_field_pheno = data_field_phenotypes.out
		Array[File]+ concept = creating_concepts.out
		File pqp = primary_care_quantitative_phenotypes.out
		File formula = formula_phenotypes.out
		File composite = composite_phenotypes.out
	}
}
