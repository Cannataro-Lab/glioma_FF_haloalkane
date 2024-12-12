


## TCGA

# Loading in TCGA GBM data downloaded from the NCI GDC July 5 2023


# obtained via 
# cancereffectsizeR::get_TCGA_project_MAF(project = "GBM",filename = "input_data/tcga_gbm.maf",exclude_TCGA_nonprimary = T)

tcga_gbm <- data.table::fread(input = "input_data/tcga_gbm.maf")





# tcga_gbm |> 
#   count(Unique_Patient_Identifier)





## GLASS

# Glass data downloaded from https://www.synapse.org/#!Synapse:syn17038081/wiki/585622 , 2022-05-31 release 
  

glass_data_passgeno <- fread(input = "input_data/synapse_GLASS_variants_passgeno.csv")

glass_data_anno <- fread(input = "input_data/synapse_GLASS_variants_anno")

# glass_data_passgeno |>  
#   count(aliquot_barcode)

# GLSS-19-0267-

# GLSS_ex <- glass_data_passgeno |> 
#   filter(stringr::str_detect("GLSS-19-0267-", string = aliquot_barcode))

GLASS_wxs_primary <- glass_data_passgeno |> 
  # filter(stringr::str_detect("-WXS-", string = aliquot_barcode)) |> 
  filter(stringr::str_detect("GLSS-", string = aliquot_barcode)) |>
  filter(stringr::str_detect("-TP-", string = aliquot_barcode)) |> 
  filter(ssm2_pass_call == "t")

# verify one sample per patient 

# GLASS_wxs_primary |> 
#   mutate(substr = stringr::str_sub(aliquot_barcode,1,12)) |> 
#   count(aliquot_barcode, substr)
# 
# GLASS_wxs_primary |>  
#   count(aliquot_barcode)

# yes

# # use same chromosome names 
# glass_data_anno$chrom |> unique() %in% GLASS_wxs_primary$chrom |>  unique()
# GLASS_wxs_primary$chrom |>  unique() %in% glass_data_anno$chrom |> unique()

glass_data_anno_match <- glass_data_anno |>  
  mutate(matcher = paste(chrom,start,end,alt)) |>  
  select(ref, matcher)

# now with reference allele 
GLASS_wxs_primary <- GLASS_wxs_primary |> 
  mutate(matcher = paste(chrom,start,end,alt)) |>  
  left_join(glass_data_anno_match, by = "matcher")

glass_clin_surgeries <- fread(input = "input_data/clinical_surgeries.csv")

glass_clin_surgeries <- glass_clin_surgeries |> 
  select(sample_barcode, histology)

GLASS_wxs_primary <- GLASS_wxs_primary |> 
  mutate(sample_barcode = stringr::str_sub(aliquot_barcode,1,15))

GLASS_wxs_primary <- GLASS_wxs_primary |> 
  left_join(glass_clin_surgeries, by = "sample_barcode")


# filter by glioblastoma 

GLASS_wxs_primary <- GLASS_wxs_primary |>  
  filter(histology == "Glioblastoma")







# GLASS_wxs_primary |> 
#   count(sample_barcode) |> 
#   ggplot(aes(x=n)) + 
#   geom_histogram()



# Initiate cancereffectsizeR analysis


# 
# 
# ff_maf <- preload_maf(maf = all_data_ff, 
#                       refset = "ces.refset.hg38",
#                       sample_col = "tumor_sample",
#                       ref_col = "Ref",
#                       chr_col = "Chrom",
#                       start_col = "Position",
#                       tumor_allele_col = "Alt")
# 
# non_ff_maf <- preload_maf(maf = all_data_nonFF, 
#                           refset = "ces.refset.hg38",
#                           sample_col = "tumor_sample",
#                           ref_col = "Ref",
#                           chr_col = "Chrom",
#                           start_col = "Position",
#                           tumor_allele_col = "Alt")


tcga_maf <- preload_maf(maf = tcga_gbm, 
                        refset = "ces.refset.hg38")


GLASS_wxs_primary$chrom[which(GLASS_wxs_primary$chrom == 23)] <- "chrX"


# split GLASS into whole exome and whole genome
GLASS_wgs_primary <- GLASS_wxs_primary |> 
  filter(stringr::str_detect(string = aliquot_barcode, pattern = "-WGS-"))
GLASS_wxs_primary <- GLASS_wxs_primary |> 
  filter(stringr::str_detect(string = aliquot_barcode, pattern = "-WXS-"))


# remove tumors from wxs that are also in wgs
GLASS_wxs_primary <- GLASS_wxs_primary[-which(GLASS_wxs_primary$sample_barcode %in% GLASS_wgs_primary$sample_barcode),]


glass_wxs_maf <- preload_maf(maf = GLASS_wxs_primary, 
                             refset = "ces.refset.hg38",
                             sample_col = "sample_barcode",
                             ref_col = "ref",
                             chr_col = "chrom",
                             start_col = "start",
                             tumor_allele_col = "alt",
                             chain_file = "input_data/chain_file/hg19ToHg38.over.chain")

glass_wgs_maf <- preload_maf(maf = GLASS_wgs_primary, 
                             refset = "ces.refset.hg38",
                             sample_col = "sample_barcode",
                             ref_col = "ref",
                             chr_col = "chrom",
                             start_col = "start",
                             tumor_allele_col = "alt",
                             chain_file = "input_data/chain_file/hg19ToHg38.over.chain")


# saveRDS(object = glass_wgs_maf,file = "output_data/glass_wgs_maf.rds")



cesa <- CESAnalysis(refset = "ces.refset.hg38")


cesa <- load_maf(cesa, maf = glass_wgs_maf, maf_name = "GLASS_wgs",coverage = "genome")

cesa <- cesa |> 
  load_maf(maf = glass_wxs_maf,maf_name = "GLASS_wxs")


cesa <- cesa |> 
  load_maf(maf = tcga_maf, maf_name = "TCGA")

# cesa <- cesa |> 
#   load_maf(maf = non_ff_maf, maf_name = "non_FF") 
# 
# cesa <- cesa |> 
#   load_maf(maf = ff_maf, maf_name = "FF")

cesa <- load_maf(cesa = cesa, maf = all_data_filtered_maf, maf_name = "FF_and_non_FF")

# save_cesa(cesa = cesa, file = "output_data/cesa_objects/cesa_after_loadmaf.rds")







# Cancer effect size, single run on all data 






# cesa <- load_cesa(file = "output_data/cesa_objects/cesa_after_loadmaf.rds")

cesa <- gene_mutation_rates(cesa = cesa, covariates = ces.refset.hg38$covariates$GBM)

signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "GBM", treatment_naive = TRUE)

cesa <- trinuc_mutation_rates(cesa = cesa,
                              signature_set = "COSMIC_v3.2",
                              signature_exclusions = signature_exclusions,
                              cores = 4
)


cesa <- ces_variant(cesa = cesa, cores = 4)





