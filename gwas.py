gwas_do.py

# el VCF con los genotipos
ts_vcf = config.VCF_TS
# un fichero con los tamaños de los cromosomas
genome_fai_path = TOMATO_GENOME_40_FAI

# Creación de los ficheros de entrada del plink
ts_snps_bfiles_res = plink.create_plink_bfiles(
            vcf_path=ts_vcf, out_base_path=ts_snps_bfiles_base_path
        )

# SNPs con poco LD para hacer el PCA        
plink.create_ld_indep_variants_file(
            bfiles_base_path=ts_snps_bfiles_res["bfiles_base_path"],
            out_base_path=pca_base_path,
            pruned_vars_list_path=pruned_vars_list_path,
            variant_filters=plink.VariantFilters(max_missing_rate=pca_max_missing_rate),
            window_size=ld_pruning_window_size,
            step_size=ld_pruning_step_size,
        )
        
# PCA para poder descontar el efecto de la estructura

pca_variant_filters = plink.VariantFilters(
            max_missing_rate=pca_max_missing_rate,
            lists_of_vars_to_keep_paths=[pruned_vars_list_path],
        )
pca_res = plink.do_pca(
            ts_snps_bfiles_base_path,
            out_base_path=pca_base_path,
            variant_filters=pca_variant_filters,
        )
        
# El GWAS usando el plink
gwas_kwargs = {
            "bfiles_base_path": gwas_bfiles_res["bfiles_base_path"],
            "phenotype_dframe": quant_traits_data,
            "out_base_name": base_name,
            "out_dir": out_dir,
            "traits": traits_to_use,
            "genome_fai_path": genome_fai_path,
            "qualitative": True,
        }
if structure_method == "pca":
	gwas_kwargs["covars_path"] = pca_res["eigenvec_path"]
do_gwas_analysis(**gwas_kwargs)


# Lectura de los ficheros de resultados
gwas_results.py

