version 1.0

task minimum_data {

	input {
		Array[File]+ files = ["dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/ukb_test_tab.csv"]
		File? exclusions = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/exclusions.csv"
		String save_loc = "."
	}

	Int cores = 16
	String min_data = save_loc + "/minimum_tab_data.gz"

	command <<<
		[ -d ~{save_loc} ] || mkdir ~{save_loc}

		Rscript `Rscript -e 'cat(system.file("extdata/scripts/phenotype_generation","01_minimum_data.R", package = "DeepPheWAS"))'` \
			--data_files ~{sep="," files} \
			~{"--exclusions " + exclusions} \
			--save_loc ~{min_data} \
			--N_cores ~{cores}
	>>>

	output {
		File out = min_data
	}

	runtime {
		cpu: cores
		memory: "120 GB"
	}
}

task data_preparation {

	input {
		String save_loc = "."
		File min_data
		File GPC = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/gp_clinical.tsv"
		File GPP = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/gp_scripts.tsv"
		File hesin_diag = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/hesin_diag.tsv"
		File HESIN = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/hesin.tsv"
		File hesin_opr = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/hesin_oper.tsv"
		File death_cause = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/death_cause.tsv"
		File death = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/death.tsv"
		File? king_coef = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/ukb43027_rel_s488264.dat"
	}

	command <<<
		[ -d ~{save_loc} ] || mkdir ~{save_loc}

		Rscript `Rscript -e 'cat(system.file("extdata/scripts/phenotype_generation","02_data_preparation.R", package = "DeepPheWAS"))'` \
			--save_location ~{save_loc} \
			--min_data ~{min_data} \
			--GPC ~{GPC} \
			--GPP ~{GPP} \
			--hesin_diag ~{hesin_diag} \
			--HESIN ~{HESIN} \
			--hesin_oper ~{hesin_opr} \
			--death_cause ~{death_cause} \
			--death ~{death} \
			--king_coef ~{king_coef}
	>>>

	output {
		File health_records = "health_records.txt.gz"
		File combined_sex = "combined_sex"
		File? control_exclusions = "control_exclusions"
		File? related_callrate = "related_callrate"
        File GP_P = "GP_P_edit.txt.gz"
		File GP_C = "GP_C_edit.txt.gz"
		File? GP_ID = "GP_C_ID.txt.gz"
		File? control_populations = "control_populations"
		File DOB = "DOB"
		Array[File] out = glob(save_loc + "/*")
	}

	runtime {
		cpu: 1
		memory: "120 GB"
	}
}
