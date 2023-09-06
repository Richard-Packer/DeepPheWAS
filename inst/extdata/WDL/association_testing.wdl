version 1.0

workflow association_testing {

	input {
		Boolean prepare_phenotypes = false
		File? phewas_manifest = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/PheWAS_manifest.csv"
		File snp_list
		String analysis_name
		File? covariates = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/covariates_plink"
		File? phenotype_inclusion_file
	}

	if (prepare_phenotypes) {
		call phenotype_preparation {
			input:
				phewas_manifest = phewas_manifest
		}

	}
	
	call extracting_snps {
		input:
			snp_list = snp_list,
			variant_save_name = analysis_name
	}

	call PLINK_association_testing {
		input:
			pgen = extracting_snps.pgen,
			psam = extracting_snps.psam,
			pvar = extracting_snps.pvar,
			covariates = covariates,
			analysis_name = analysis_name,
			phewas_manifest = phewas_manifest,
			phenotype_inclusion_file = phenotype_inclusion_file
	}
	
	call tables_graphs {
		input:
			results_file = PLINK_association_testing.results,
        	analysis_name = analysis_name,
        	snp_list = snp_list,
			phewas_manifest = phewas_manifest
	}

	output {
		Array[File?] prep_pheno_out = if prepare_phenotypes then [ phenotype_preparation.phenotype, phenotype_preparation.stats ] else []
		File assoc_out = PLINK_association_testing.results
		Array[File?] tables_graphs_out = tables_graphs.out
	}
}

task phenotype_preparation {

	input {
		Array[File]+ phenotype_files = [ "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/composite_phenotypes.RDS", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/concepts.RDS", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/all_dates.RDS", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/range_ID.RDS", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/phecodes.RDS", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/PQP.RDS", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/formula_phenotypes.RDS", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/data_field_phenotypes.RDS" ]
		File groupings = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/ancestry_panUKB"
		String phenotype_filtered_save_name = "IVNT_RR_N_filtered_neuro_PanUKB_ancestry"
		Boolean relate_remove = true
		File? kinship_file = "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/related_callrate"
		Boolean IVNT = true
		String stats_save = "IVNT_RR_N_filtered_neuro_PanUKB_ancestry_stats"
		File? phewas_manifest
	}

	command <<<
		SCRIPT=`Rscript -e 'cat(system.file("extdata/scripts/association_testing","01_phenotype_preparation.R", package = "DeepPheWAS"))'` && \
		echo $SCRIPT && \
		Rscript $SCRIPT \ 
			--phenotype_filtered_save_name ~{phenotype_filtered_save_name} \
			--phenotype_files ~{sep="," phenotype_files} \
			~{"--groupings " + groupings } \
			~{true="--relate_remove" false="" relate_remove} \
			~{"--kinship_file " + kinship_file} \
			~{true="--IVNT" false="" IVNT} \
			--stats_save ~{stats_save} \
			~{"--PheWAS_manifest_overide " + phewas_manifest} \
	>>>

	output {
		File phenotype = phenotype_filtered_save_name
		File stats = stats_save
	}

	runtime {
		cpu: 1
		memory: "200 GB"
	}
}

task extracting_snps {

	input {
		Array[File]+ bgens = ["dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c10_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c6_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c17_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c9_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c14_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c21_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c4_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cXY_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c11_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c20_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c19_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c7_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c12_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c22_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c18_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c5_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c1_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c15_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c8_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c2_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c16_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cX_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c3_b0_v3.bgen","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c13_b0_v3.bgen"]
		Array[File]+ bgis = ["dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c10_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c6_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c17_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c9_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c14_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c21_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c4_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cXY_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c11_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c20_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c19_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c7_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c12_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c22_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c18_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c5_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c1_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c15_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c8_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c2_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c16_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cX_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c3_b0_v3.bgen.bgi","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c13_b0_v3.bgen.bgi"]
		Array[File]+ sample_files = ["dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c10_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c6_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c2_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c11_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c8_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c15_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c4_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c12_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c5_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c13_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cXY_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c17_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c18_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c14_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c7_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c19_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c1_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c9_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c16_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cX_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c21_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c22_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c3_b0_v3.sample","dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c20_b0_v3.sample"]
		File snp_list
		String variant_save_name
	}

	command <<<
		echo ~{sep="," bgis} &&
		paste <(tr ',' '\n' <<< "~{sep=',' bgens}" | sort) <(tr ',', '\n' <<< "~{sep=',' sample_files}" | sort) | gawk 'BEGIN { OFS = ","; print "chromosome,genetic_file_location,psam_fam_sample_file_location" } { sub(".bgen", "", $1); print gensub(".+ukb22828_c(.+)_b0_v3", "\\1", "g", $1), $1, $2 }' > genetic_file_guide.csv &&
		Rscript `Rscript -e 'cat(system.file("extdata/scripts/association_testing","02_extracting_snps.R", package = "DeepPheWAS"))'` \
		--SNP_list ~{snp_list} \
		--genetic_file_guide genetic_file_guide.csv \
		--analysis_folder . \
		--bgen_input \
		--ref_bgen ref-first \
		--variant_save_name ~{variant_save_name}
	>>>

	output {
		File pgen = variant_save_name + ".pgen"
		File psam = variant_save_name + ".psam"
		File pvar = variant_save_name + ".pvar"
	}

	runtime {
		cpu: 1
		memory: "24 GB"
	}
}	

task PLINK_association_testing {

	input {
		Array[File]+ phenotypes = [ "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/phenotypes/AFR_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/phenotypes/AMR_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/phenotypes/CSA_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/phenotypes/EAS_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/phenotypes/EUR_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz", "dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/phenotypes/MID_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz" ]
		File pgen
		File psam
		File pvar
		File? covariates
		String analysis_name
		File? phewas_manifest
		File? phenotype_inclusion_file
	}

	Int cores = 8

	command <<<
		echo ~{psam} ~{pvar} &&
		`Rscript -e 'cat(system.file("extdata/scripts/association_testing","03a_PLINK_association_testing.R", package = "DeepPheWAS"))'` \
		--analysis_folder ./ \
		--phenotype_files ~{sep="," phenotypes} \
		--variants_for_association ~{pgen} \
		~{"--PheWAS_manifest_overide " + phewas_manifest} \
		~{"--covariate " + covariates} \
		--analysis_name ~{analysis_name} \
		--plink_exe "plink2 --threads ~{cores}" \
		--save_plink_tables \
		~{"--phenotype_inclusion_file " + phenotype_inclusion_file}
	>>>

	output {
		File results = "association_results/" + analysis_name + "_association_results_list.gz"
	}

	runtime {
		cpu: cores
		memory: "200 GB"
	}
}

task tables_graphs {

	input {
		File results_file
		String analysis_name
		File snp_list
		File? phewas_manifest
	}

	command <<<
		`Rscript -e 'cat(system.file("extdata/scripts/association_testing","05_tables_graphs.R", package = "DeepPheWAS"))'` \
			--results_file ~{results_file} \
			--analysis_name ~{analysis_name} \
			--plink_results \
			--SNP_list ~{snp_list} \
			~{"--PheWAS_manifest_overide " + phewas_manifest} \
			--per_group_name_graph \
			--save_table_per_snp \
			--save_table_per_group_name \
			--sex_split \
			--save_folder "./"
	>>>

	output {
		Array[File]+ out = glob("*")
	}

	runtime {
		cores: 1
		memory: "64 GB"
	}
}
