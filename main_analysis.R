
# main analysis script 





# packages ---------------------------------------------------------------------
library(cancereffectsizeR)
library(tidyverse)
library(plotly)
library(ggrepel)
library(data.table)
library(patchwork)

# for Boschloo test 
library(exact2x2)
library(parallel)


# load in raw project data -----------------------------------------------------
source("R/load_project_data.R")


# data filtering ---------------------------------------------------------------


# our previous filtering https://www.synapse.org/#!Synapse:syn31770233/wiki/617920
# used the "common_variant" call from vcf2maf
# https://github.com/mskcc/vcf2maf/blob/main/docs/vep_maf_readme.txt
# which is "allele frequency across at least one gnomAD subpopulation is >0.04%"
#

all_data_filtered <- all_data |>
  filter(`Max. Sub-Population Frequency` < 0.0004) |> 
  filter(`Max. Population Frequency` < 0.0004)


# load into cesR to annotate 
all_data_filtered_maf <- preload_maf(maf = all_data_filtered,
                                     refset = "ces.refset.hg38",
                                     sample_col = "tumor_sample",
                                     ref_col = "Ref",
                                     chr_col = "Chrom",
                                     start_col = "Position",
                                     tumor_allele_col = "Alt",keep_extra_columns = T)

# filter out repetitive regions
all_data_filtered_maf <- all_data_filtered_maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]



# +Boschloo test --------------------------------------------------------------- 

# only want variants in the next test that have allele frequencies 
# (NA introduced to some variants in preload_maf() when snv doublets re-mapped to dbs)
 
 

all_data_filtered_maf <- all_data_filtered_maf |>  
  filter(!is.na(`RefDepth(Tumor)`))


# do not want potential dbs variants in analysis

all_data_filtered_maf <- all_data_filtered_maf |> 
  filter(is.na(problem))


# This takes a long time, so running once then commenting out. 
# Rerun the script if packages are updated or input changes 


# cores_to_use <- 8
# source("R/Boschloo_test.R") # saves pvals.rds to output_data/

pvals <- readRDS("output_data/pvals.rds")

pvals_min <- sapply(pvals, function(x) x$btest_min_p)

all_data_filtered_maf <- all_data_filtered_maf |> mutate(pval_min = pvals_min)


saveRDS(all_data_filtered_maf,file = "output_data/all_data_filtered_maf_afterBoschloo.rds")



# load in data with Boschloo P values saved in the above script 
all_data_filtered_maf <- readRDS("output_data/all_data_filtered_maf_afterBoschloo.rds")


all_data_filtered_maf <- all_data_filtered_maf |> 
  filter(pval_min <= 0.05)


# save final data from this analysis post-filtering 
saveRDS(all_data_filtered_maf, file = "output_data/all_data_filtered_maf.rds")

# save final data after cleaning for easy access
write_csv(x = all_data_filtered_maf, file = "output_data/sequencing_data_filtered_maf.csv")


# for easy access to filtered data
# all_data_filtered_maf <- read_csv(file = "output_data/sequencing_data_filtered_maf.csv")

# Bootstrap analysis -----------------------------------------------------------


# initiate cancereffectsizeR analysis to facilitate signature deconvolution
cesa <- CESAnalysis(refset = "ces.refset.hg38")

cesa <- cesa |> 
  load_maf(maf = all_data_filtered_maf, maf_name = "all_data_filtered") 


bootstrap_num <- 1000
set.seed(1234) # reproducible analyses

# run bootstrap loop 
source("R/bootstrap_analysis.R")


# save for subsequent analyses
saveRDS(trinuc_list_raw, file = "output_data/bootstrap/trinuc_list_raw_pval_min.rds")
saveRDS(object = trinuc_list_biological, file = "output_data/bootstrap/trinuc_list_biological_pval_min.rds")
saveRDS(selection_list,file = "output_data/bootstrap/selection_list_pval_min.rds")
saveRDS(attribution_list, file = "output_data/bootstrap/attribution_list_pval_min.rds")



# load in bootstrap analyses-----
trinuc_raw <- readRDS(file = "output_data/bootstrap/trinuc_list_raw_pval_min.rds")
biological_raw <- readRDS(file = "output_data/bootstrap/trinuc_list_biological_pval_min.rds")
selection_raw <- readRDS(file = "output_data/bootstrap/selection_list_pval_min.rds")
attributions_raw <- readRDS(file = "output_data/bootstrap/attribution_list_pval_min.rds")




# run cancereffectsizeR analysis from multiple data sources for 
# gene-level significance values  
# cesR single run -------------------------------------------------------------- 

# restarts cesa object using all_data_filtered_maf and 
# tcga / GLASS data which is loaded in within the following script: 

source("R/cesa_singlerun.R")
save_cesa(cesa = cesa, file = "output_data/cesa_singlerun.rds")

# saves output, can reload here: 

cesa <- load_cesa(file = "output_data/cesa_singlerun.rds")




# Figure: SBS42 v. FF years v. occupation highlights ---------------------------

source("R/figure_SBS42_v_FFyear_v_occupation.R")

job_fig_median




# Figure: SBS42 v. FF years v. variant highlights ------------------------------


source("R/figure_SBS42_v_FFyear_v_variants.R")

variant_fig_median







# figure build -----------------------------------------------------------------

# job_fig + variant_fig + plot_annotation(tag_levels = 'A')
# 
# ggsave(filename = "figures/SBS42_v_ffyears_JobAndVariants.png",
#        width = 14,height = 6)


fig1 <- job_fig_median + variant_fig_median + plot_annotation(tag_levels = 'A')

ggsave(plot = fig1,filename = "figures/SBS42_v_ffyears_JobAndVariants_median.png",
       width = 14,height = 5)
ggsave(plot = fig1,filename = "figures/SBS42_v_ffyears_JobAndVariants_median.pdf",
       width = 14,height = 5)
ggsave(plot = fig1,filename = "figures/SBS42_v_ffyears_JobAndVariants_median.EPS",
       width = 14,height = 5)

# 


# statistics -------------------------------------------------------------------



# + number of samples w/ positive mean SBS 42 signature?

trinuc_raw_means_cov_alljob_interest |> 
  filter(median_value > 0) |> 
  nrow() 

nrow(trinuc_raw_means_cov_alljob_interest)

nrow(trinuc_raw_means_cov_alljob_interest)



# + correlation with FF years? ----- 


FF_no_occupation_highlighted <- trinuc_raw_means_cov_alljob_interest |> 
  filter(FF_years>0) |>  # filter for firefighters
  filter(is.na(occupation)) # filter for samples without occupation highlighted


FF_no_occupation_highlighted |> 
  ggplot(aes(x=FF_years, y=median_value)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  coord_cartesian(xlim = c(0,40),y=c(0,15))

FF_no_occupation_highlighted_lm <- lm(data = FF_no_occupation_highlighted, median_value ~ FF_years)

FF_no_occupation_highlighted_lm_summary <- summary(FF_no_occupation_highlighted_lm)

FF_no_occupation_highlighted_lm_summary




# output SBS data --------------------------------------------------------------

write.csv(x = trinuc_raw_means_cov_alljob_interest,file = "output_data/cov_data_w_SBS42.csv")






# Tumor mutation burden ------ 

TMB_all <- all_data_filtered_maf |> 
  count(Unique_Patient_Identifier,ff) |> 
  arrange(ff)

write_csv(x = TMB_all, file = "output_data/TMB_all.csv")

summary(TMB_all |> filter(ff==T) |> pull(n))
summary(TMB_all |> filter(ff==F) |> pull(n))




# signature summary ------ 

median_tmb_among_bootstrap <- trinuc_raw |> 
  data.table::rbindlist() |> 
  pivot_longer(cols = -Unique_Patient_Identifier,
               names_to = "sbs_signature",
               values_to = "attributed_TMB") |> 
  group_by(Unique_Patient_Identifier,sbs_signature) |>
  summarize(median_attr = median(attributed_TMB)) 
  


signatures_detectable <- setdiff(colnames(cesa$mutational_signatures$raw_attributions),signature_exclusions)
signatures_detectable <- signatures_detectable[grep(pattern = "SBS",x = signatures_detectable)] # keep just signature names


median_tmb_among_bootstrap_detect <- median_tmb_among_bootstrap |> 
  left_join(TMB_all) |>
  filter(sbs_signature %in% signatures_detectable) |> 
  dplyr::rename(TMB = n)


mean_median_tmbattr <- median_tmb_among_bootstrap_detect |> 
  mutate(median_boot_percent = median_attr / TMB) |> 
  group_by(sbs_signature,ff) |> 
  summarize(mean_median_boot_percent = mean(median_boot_percent)) |> 
  arrange(desc(mean_median_boot_percent))

write_csv(x = mean_median_tmbattr |> arrange(ff), file = "output_data/mean_median_attr_SBS.csv")


mean_median_tmbattr |> 
  filter(ff==T) |> 
  print(n = Inf)
  # filter(mean_median_boot_percent>0)


mean_median_tmbattr |> 
  filter(ff==F) |> 
  print(n=Inf)
  # filter(mean_median_boot_percent>0)


