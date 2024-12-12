
Files you need to add to `/input_data` to run the full analysis. 


- Sequencing data and covariate data associated with this manuscript (to be published along with manuscript)
- Synapse GLASS data, accessed from  https://www.synapse.org/#!Synapse:syn17038081/wiki/585622 , 2022-05-31 release 
- TCGA GBM data, may be stored with `cancereffectsizeR::get_TCGA_project_MAF(project = "GBM",filename = "input_data/tcga_gbm.maf",exclude_TCGA_nonprimary = T)`
- COSMIC cancer gene census, can be obtained by downloading from https://cancer.sanger.ac.uk/census 
- hg19toHg38 chain file, obtained via https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/ 