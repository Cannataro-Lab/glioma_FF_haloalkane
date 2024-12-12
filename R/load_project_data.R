

# data from R/original_data_anonymize.R


cov_data <- read_csv("input_data/cov_data.csv")

all_data <- read_csv("input_data/all_data.csv")


# 
# 
# 
# 
# # directories from collaboration
# # to be updated once data is published online with manuscript
# 
# non_ff_data_directory <- "~/Library/CloudStorage/Box-Box/Firefighter_Study/WES Reports_18 non-frefighter gliomas_Dr. Knight YCGA 10-25-23/"
# 
# ff_data_directory <- "~/Library/CloudStorage/Box-Box/Firefighter_Study/Whole Exome Tumor & Normal Sequencing WES Data/"
# 
# cov_file <- "~/Library/CloudStorage/Box-Box/Covariate files/RAP_covariates_55_for_WES_NTA_2023-08-16.csv"
# 
# 
# 
# 
# # load in non FF data ---- 
# 
# # all of the files 
# xl_files <- dir(path = non_ff_data_directory,
#                 pattern = "report.xlsx",full.names = T)
# 
# # load in just somatic data
# all_somatic_report <- purrr::map(xl_files,readxl::read_excel,sheet = "Somatic")
# 
# # purrr::map(all_somatic_report,nrow)
# # all_somatic_report
# 
# 
# tumor_names <- dir(path = non_ff_data_directory,
#                    pattern = "report.xlsx")
# 
# tumor_names <- stringr::str_remove(string = tumor_names,pattern = ".report.xlsx")
# 
# names(all_somatic_report) <- tumor_names
# # 
# 
# all_data_nonFF <- dplyr::bind_rows(all_somatic_report,.id = "tumor_sample")
# 
# 
# 
# # load in FF data ---- 
# 
# 
# 
# xl_files <- dir(path = ff_data_directory,
#                 pattern = "report.xlsx",full.names = T)
# 
# 
# all_somatic_report <- purrr::map(xl_files,readxl::read_excel,sheet = "Somatic")
# 
# # purrr::map(all_somatic_report,nrow)
# # all_somatic_report
# 
# tumor_names <- dir(path = ff_data_directory,
#                    pattern = "report.xlsx")
# 
# tumor_names <- stringr::str_remove(string = tumor_names,pattern = ".report.xlsx")
# 
# names(all_somatic_report) <- tumor_names
# 
# 
# all_data_ff <- dplyr::bind_rows(all_somatic_report,.id = "tumor_sample")
# 
# 
# 
# all_data <- bind_rows(all_data_nonFF |> mutate(ff="F"),
#                       all_data_ff |> mutate(ff="T"))
# 
# 
# 
# 
# cov_data <- fread(input = cov_file)


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




# all_data <- readRDS("troubleshoot/all_data.rds")




