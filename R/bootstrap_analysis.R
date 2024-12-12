

# run the bootstrap analysis on data to get distribution of SBS attribution

trinuc_list_raw <- vector(mode = "list",length = bootstrap_num)
trinuc_list_biological <- vector(mode = "list",length = bootstrap_num)
selection_list <- vector(mode = "list",length = bootstrap_num)
attribution_list <- vector(mode = "list",length = bootstrap_num)

pbapply::pboptions(type = "none") # runs faster 

# cesa <- gene_mutation_rates(cesa = cesa, covariates = ces.refset.hg38$covariates$GBM)

signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "GBM", treatment_naive = TRUE)

cesa <- gene_mutation_rates(cesa = cesa, covariates = ces.refset.hg38$covariates$GBM)


for(boot_ind in 1:bootstrap_num){
  
  cesa <- trinuc_mutation_rates(cesa = cesa,
                                signature_set = "COSMIC_v3.2",
                                signature_exclusions = signature_exclusions,
                                bootstrap_mutations = T,cores = 4,
                                sig_averaging_threshold = 1
  )
  
  
  trinuc_list_raw[[boot_ind]] <- cesa$mutational_signatures$raw_attributions
  trinuc_list_biological[[boot_ind]] <- cesa$mutational_signatures$biological_weights
  # 
  cesa <- ces_variant(cesa = cesa, cores = 4,variants = cesa$variants)
  # 
  selection_list[[boot_ind]] <- cesa$selection$selection.1
  # 
  attribution_list[[boot_ind]] <- mutational_signature_effects(cesa = cesa,effects = cesa$selection$selection.1)
  
  if(boot_ind < bootstrap_num){cesa <- clear_trinuc_rates_and_signatures(cesa)}
  cesa <- clear_effect_output(cesa)
  
  print(boot_ind)
  
}


