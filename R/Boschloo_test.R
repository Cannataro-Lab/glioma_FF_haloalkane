# Boschloo's exact test script 


# filtering step from previous analyses e.g. https://www.synapse.org/#!Synapse:syn31770233/wiki/617920


# input: all_data_filtered_maf from main_analysis.R 


# build listed vectors

# all_data_filtered_maf <- all_data_filtered_maf |>  
#   filter(variant_type == "snv")

VAF_list <- vector(mode = "list",length = nrow(all_data_filtered_maf))


for(rowind in 1:nrow(all_data_filtered_maf)){
  
  VAF_list[[rowind]] <- c(all_data_filtered_maf$`AltDepth(Tumor)`[rowind],
                          all_data_filtered_maf$`RefDepth(Tumor)`[rowind],
                          all_data_filtered_maf$`AltDepth(Normal)`[rowind],
                          all_data_filtered_maf$`RefDepth(Normal)`[rowind]) 
  
}


boschloo_tester <- function(input_vec,
                            test="two.sided"
){
  
  

  b_test_cen_p <- boschloo(x1 = input_vec[1], n1 = input_vec[1] + input_vec[2],
                           x2 = input_vec[3], n2 = input_vec[3] + input_vec[4],
                           alternative = test)

  b_test_minlike <- boschloo(x1 = input_vec[1], n1 = input_vec[1] + input_vec[2],
                             x2 = input_vec[3], n2 = input_vec[3] + input_vec[4],
                             alternative = test#,
                             #tsmethod = "minlike"
                             )
  
  
  # if(input_index %% 1 == 0) system(paste("echo 'now processing:",input_index,"'"))
  
  
  
  return(list(
    btest_cen_p = b_test_cen_p$p.value,
    btest_min_p = b_test_minlike$p.value
  ))
  
  
  
} 


pvals <- parallel::mclapply(X = VAF_list, FUN = boschloo_tester,mc.cores = cores_to_use)





saveRDS(pvals, file = "output_data/pvals.rds")


