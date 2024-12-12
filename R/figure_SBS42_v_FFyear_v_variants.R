
# figure: variants of interest determine 




dndscv_results <- cesa$dNdScv_results$rate_grp_1 |> 
  arrange(qallsubs_cv)

dndscv_sig_genes <- dndscv_results |> 
  filter(qallsubs_cv <= 0.1) |> 
  pull(gene_name)


attribution_df <- data.table::rbindlist(lapply(attributions_raw, function(x) x$mutational_sources$source_probabilities),idcol = "n_boot",fill = T)

# just project data 
attribution_df <- attribution_df |> 
  filter(Unique_Patient_Identifier %in% trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier)

project_maf_sig_genes <- cesa$maf |> 
  filter(Unique_Patient_Identifier %in% trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier) |> 
  filter(top_gene %in% dndscv_sig_genes)

attribution_df_long <- attribution_df |> 
  pivot_longer(cols = starts_with("SBS")) |> 
  group_by(variant_id, Unique_Patient_Identifier, name) |> 
  summarize(mean_val = mean(value),
            median_val = median(value)) |> 
  ungroup()

# TP53_R273H_ENSP00000269305.4
# PIK3CA_E365K_ENSP00000263967.3
# which(project_maf_sig_genes$variant_id == "TP53_R273H_ENSP00000269305.4")
# 
# # EGFR_L861Q_ENSP00000275493.2
# 
# which(attribution_df_long$variant_id == "EGFR_L861Q_ENSP00000275493.2")
# 
# which(cesa$maf$variant_id == "EGFR_L861Q_ENSP00000275493.2")
# cesa$maf[which(cesa$maf$top_consequence == "EGFR L861Q"),]

attr_sig_genes <- attribution_df_long |> 
  filter(name == "SBS42") |> 
  left_join(cesa$variants) |> 
  filter(gene %in% dndscv_sig_genes)



attr_sig_genes <- attr_sig_genes |> 
  arrange(desc(mean_val)) |> 
  filter(median_val > .3) |> 
  mutate(label = paste0(variant_name,": ", round(median_val,3)))

attr_sig_genes <- attr_sig_genes |> 
  select(Unique_Patient_Identifier, label)



# attr_sig_genes 




trinuc_raw_means_cov_alljob_interest$variant_label <- NA
trinuc_raw_means_cov_alljob_interest$variant_inclusion <- NA



tumors_to_annot <- unique(attr_sig_genes$Unique_Patient_Identifier)

for(tumor_ind in 1:length(tumors_to_annot)){
  
  this_tumor <- tumors_to_annot[tumor_ind]
  
  if(length(which(attr_sig_genes$Unique_Patient_Identifier == this_tumor)) == 1){
    this_label <- attr_sig_genes$label[which(attr_sig_genes$Unique_Patient_Identifier == this_tumor)]
    trinuc_raw_means_cov_alljob_interest$variant_label[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == this_tumor)] <- this_label
    trinuc_raw_means_cov_alljob_interest$variant_inclusion[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == this_tumor)] <- "Significantly mutated gene"
  }
  
  if(length(which(attr_sig_genes$Unique_Patient_Identifier == this_tumor)) > 1){
    
    trinuc_raw_means_cov_alljob_interest$variant_label[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == this_tumor)] <- paste(attr_sig_genes$label[which(attr_sig_genes$Unique_Patient_Identifier == this_tumor)],collapse = "; ")
    trinuc_raw_means_cov_alljob_interest$variant_inclusion[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == this_tumor)] <- "Significantly mutated gene"
    
  }
}




# if( length(unique(attr_sig_genes$Unique_Patient_Identifier)) == nrow(attr_sig_genes)){
#   for(row_ind in 1:nrow(attr_sig_genes)){
#     
#     trinuc_raw_means_cov_alljob_interest$variant_label[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == attr_sig_genes$Unique_Patient_Identifier[row_ind])] <- attr_sig_genes$label[row_ind]
#     
#   }
# }else{
#   message("More than one significantly mutated gene in a sample!")
# }


variant_fig <- trinuc_raw_means_cov_alljob_interest |> 
  ggplot(aes(x=FF_years, y=mean_value)) + 
  geom_point() + 
  geom_text_repel(aes(label = variant_label),
                  min.segment.length = unit(0, 'lines'),
                  box.padding = 0.5,nudge_x = 1) + 
  labs(y="Mean number of mutations attributable to SBS42", 
       x= "Firefighter years",color="Firefighter status") +
  theme_bw()


cosmic_genes <- fread("input_data/Census_allTue May 7 18 41 10 2024.tsv")

# AGS_12 ----- 

AGS_12_attr <- attribution_df_long |> 
  filter(Unique_Patient_Identifier == "AGS_12") |> 
  filter(name == "SBS42") |> 
  arrange(desc(mean_val)) |> 
  left_join(cesa$selection$selection.1) 

# NOTCH1 implicated in glioma e.g. https://pubmed.ncbi.nlm.nih.gov/33076453/ 

AGS_12_attr <- AGS_12_attr |> 
  select(variant_id, mean_val,median_val) |> 
  left_join(cesa$variants) |> 
  mutate(label = paste0(variant_name,": ", round(median_val,3))) |>
  filter(gene == "NOTCH1") |> 
  left_join(cesa$selection$selection.1,by="variant_id") |> 
  left_join(cosmic_genes, by=c("gene.x"="Gene Symbol")) |> 
  relocate(variant_id, mean_val, selection_intensity,Tier, Hallmark,`Tumour Types(Somatic)`)

AGS_12_attr

this_label <- AGS_12_attr$label
trinuc_raw_means_cov_alljob_interest$variant_label[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_12")] <- this_label
trinuc_raw_means_cov_alljob_interest$variant_inclusion[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_12")] <- "COSMIC tier 1 gene"




# AGS_29 ---- 
AGS_29_attr <- attribution_df_long |> 
  filter(Unique_Patient_Identifier == "AGS_29") |> 
  filter(name == "SBS42") |> 
  left_join(cesa$variants) |> 
  arrange(desc(median_val)) |> 
  left_join(cesa$selection$selection.1,by="variant_id") |> 
  left_join(cosmic_genes, by=c("gene.x"="Gene Symbol")) |> 
  relocate(variant_id, mean_val, selection_intensity,Tier, Hallmark,`Tumour Types(Somatic)`)


# https://www.nature.com/articles/s41467-018-08087-9
# Following EGFR stimulation, CIC repression is relieved, allowing for the expression of target genes. The best-characterized CIC targets in mammalian cells are the oncogenic transcription factors ETV1, ETV4, and ETV55, which mediate cell proliferation, motility and invasion downstream of Ras6.


AGS_29_attr <- AGS_29_attr |> 
  filter(gene.x == "ETV1") |> 
  left_join(cesa$variants, by="variant_id") |> 
  mutate(label = paste0(variant_name,": ", round(median_val,3)))

this_label <- AGS_29_attr$label
trinuc_raw_means_cov_alljob_interest$variant_label[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_29")] <- this_label
trinuc_raw_means_cov_alljob_interest$variant_inclusion[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_29")] <- "COSMIC tier 1 gene"


# AGS_43 ----- 

AGS_43_attr <- attribution_df_long |> 
  filter(Unique_Patient_Identifier == "AGS_43") |> 
  filter(name == "SBS42") |> 
  left_join(cesa$variants) |> 
  arrange(desc(median_val)) |> 
  left_join(cesa$selection$selection.1,by="variant_id") |> 
  left_join(cosmic_genes, by=c("gene.x"="Gene Symbol")) |> 
  relocate(variant_id, mean_val, selection_intensity,Tier, Hallmark,`Tumour Types(Somatic)`)

# ROS1
# https://pubmed.ncbi.nlm.nih.gov/38629293/
# ROS1 facilitates the progression of various malignancies via self-mutations or rearrangements.
# https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=ROS1

AGS_43_attr <- AGS_43_attr |> 
  filter(gene.x == "ROS1") |> 
  left_join(cesa$variants, by= "variant_id") |> 
  mutate(label = paste0(variant_name,": ", round(median_val,3)))


this_label <- AGS_43_attr$label
trinuc_raw_means_cov_alljob_interest$variant_label[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_43")] <- this_label
trinuc_raw_means_cov_alljob_interest$variant_inclusion[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_43")] <- "COSMIC tier 1 gene"






# AGS_17

AGS_17_attr <- attribution_df_long |> 
  filter(Unique_Patient_Identifier == "AGS_17") |> 
  filter(name == "SBS42") |> 
  left_join(cesa$variants) |> 
  arrange(desc(mean_val)) |> 
  left_join(cesa$selection$selection.1,by="variant_id") |> 
  left_join(cosmic_genes, by=c("gene.x"="Gene Symbol")) |> 
  relocate(variant_id, mean_val, selection_intensity,Tier, Hallmark,`Tumour Types(Somatic)`)


# already highlighting dndscv result. 



# AGS_45 ----- 

AGS_45_attr <- attribution_df_long |> 
  filter(Unique_Patient_Identifier == "AGS_45") |> 
  filter(name == "SBS42") |> 
  left_join(cesa$variants) |> 
  arrange(desc(mean_val)) |> 
  left_join(cesa$selection$selection.1,by="variant_id") |> 
  left_join(cosmic_genes, by=c("gene.x"="Gene Symbol")) |> 
  relocate(variant_id, mean_val, selection_intensity,Tier, Hallmark,`Tumour Types(Somatic)`)


# https://cancer.sanger.ac.uk/cosmic/census-page/NCOA2
# NCOA2 labeled as oncogene and TSG... knock down decreases division of cell lines etc. 

AGS_45_attr <- AGS_45_attr |> 
  filter(gene.x == "NCOA2") |> 
  left_join(cesa$variants, by= "variant_id") |>
  mutate(label = paste0(variant_name,": ", round(median_val,3)))


this_label <- AGS_45_attr$label
trinuc_raw_means_cov_alljob_interest$variant_label[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_45")] <- this_label
trinuc_raw_means_cov_alljob_interest$variant_inclusion[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_45")] <- "COSMIC tier 1 gene"



# AGS_17 also has 
# CTNNB1_D32N_ENSP00000344456.5



# AGS_31 ----- 


AGS_31_attr <- attribution_df_long |> 
  filter(Unique_Patient_Identifier == "AGS_31") |> 
  filter(name == "SBS42") |> 
  left_join(cesa$variants) |> 
  arrange(desc(mean_val)) |> 
  left_join(cesa$selection$selection.1,by="variant_id") |> 
  left_join(cosmic_genes, by=c("gene.x"="Gene Symbol")) |> 
  relocate(variant_id, mean_val, selection_intensity,Tier, Hallmark,`Tumour Types(Somatic)`)

AGS_31_attr |> 
  filter(!is.na(selection_intensity) | !is.na(Tier))

# ZMYM3
# COSMIC tier 2 gene https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=ZMYM3 

# AGS_31_attr <- AGS_31_attr |> 
#   filter(gene.x == "ZMYM3") |> 
#   left_join(cesa$variants, by= "variant_id") |>
#   mutate(label = paste0(variant_name,": ", round(mean_val,3)))
# 
# 
# this_label <- AGS_31_attr$label
# trinuc_raw_means_cov_alljob_interest$variant_label[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_31")] <- this_label
# 




# 
# attr_project_sig <- attr_project |> 
#   pivot_longer(cols = starts_with("SBS")) |> 
#   group_by(variant_id, Unique_Patient_Identifier, name) |> 
#   summarize(mean_val = mean(value)) |> 
#   ungroup() |> 
#   filter(name == "SBS42") |> 
#   left_join(cesa$variants) |> 
#   filter(gene %in% dndscv_sig_genes)

# 
# 
# attr_project_sig |> 
#   filter(Unique_Patient_Identifier == "AGS_17") |> # CHECK HERE for best option
#   filter(gene == "PIK3CA" & variant_type == "aac") |> 
#   mutate(label = paste0(variant_name,": ", round(mean_val,3))) |> 
#   pull(label) ->
#   label_AGS_17
# 
# trinuc_raw_means_cov_alljob_interest$variant_label[which(trinuc_raw_means_cov_alljob_interest$Unique_Patient_Identifier == "AGS_17")] <- 
#   label_AGS_17
# 

# 



trinuc_raw_means_cov_alljob_interest |> 
  ggplot(aes(x=FF_years, y=mean_value)) + 
  geom_point() + 
  geom_text_repel(aes(label = Unique_Patient_Identifier),
                  min.segment.length = unit(0, 'lines'),box.padding = 0.4,nudge_x = 1) + 
  labs(y="Mean number of mutations attributable to SBS42", 
       x= "Firefighter years",color="Variant inclusion criteria") +
  theme_bw()


trinuc_raw_means_cov_alljob_interest$variant_inclusion <- factor(
  trinuc_raw_means_cov_alljob_interest$variant_inclusion,
  levels = c("Significantly mutated gene","COSMIC tier 1 gene","COSMIC tier 2 gene"))


# dxcode	Final diagnosis per study pathologist (or UCSF pathologist)	
# 01	gbm
# 02	anap astro
# 03	astro gr 2
# 04	anap oligo
# 05	oligo gr 2
# 06	anap oligoastro
# 07	oligoastro gr 2
# 08	anap ependymoma
# 09	ependymoma
# 10	jpa
# 11	medullo
# 12	other
# 13	mixed
# 14	ganglio
# 15	oligo NOS
# 16	ependy NOS
# 17	astro NOS


trinuc_raw_means_cov_alljob_interest$dxcode_fct <- 
  fct_recode(as.character(trinuc_raw_means_cov_alljob_interest$dxcode), 
             `Glioblastoma multiforme` = "1",
             `Astrocytoma grade 3` = "2",
             `Astrocytoma grade 2` = "3",
             `Oligodendroglioma grade 3` = "4",
             `Oligodendroglioma grade 2` = "5",
             `Oligoastrocytoma grade 2` = "7",
             `Other` = "12")

trinuc_raw_means_cov_alljob_interest$dxcode_fct <- 
  fct_relevel(trinuc_raw_means_cov_alljob_interest$dxcode_fct,
              "Glioblastoma multiforme",
              "Astrocytoma grade 3",
              "Astrocytoma grade 2",
              "Oligodendroglioma grade 3",
              "Oligodendroglioma grade 2",
              "Oligoastrocytoma grade 2",
              "Other")

# variant_fig <- trinuc_raw_means_cov_alljob_interest |> 
#   ggplot(aes(x=FF_years, y=mean_value)) + 
#   # geom_label(aes(color=variant_inclusion),alpha=0,size = 0,label="") + 
#   geom_point(aes(shape=dxcode_fct),size=4,alpha=0.5) +
#   # geom_jitter(aes(shape=dxcode_fct),size=3,alpha=0.5,height = 0) + 
#   geom_text_repel(aes(label = variant_label,color=variant_inclusion),
#                   min.segment.length = unit(0, 'lines'),
#                   box.padding = 0.6,nudge_x = 2,
#                   key_glyph = "rect",fontface = "bold",size=5) + 
#   labs(y="Mean number of mutations attributable to SBS42", 
#        x= "Firefighter years",color="Variant inclusion criteria",shape = "Dx code") +
#   theme_bw() + 
#   scale_color_discrete(na.translate=F) + 
#   scale_shape_manual(values = c(16, 17, 15, 3,7,8,18))


variant_labels <- trinuc_raw_means_cov_alljob_interest |> 
  select(Unique_Patient_Identifier,FF_years,median_value,variant_label,variant_inclusion) |> 
  filter(!is.na(variant_label)) |> 
  mutate(new_x_text = case_when(
    Unique_Patient_Identifier == "AGS_08" ~ FF_years + 13,
    Unique_Patient_Identifier == "AGS_17" ~ FF_years + 9,
    Unique_Patient_Identifier == "AGS_12" ~ FF_years + 7,
    Unique_Patient_Identifier == "AGS_09" ~ FF_years + 7.5,
    Unique_Patient_Identifier == "AGS_43" ~ FF_years + 5,
    Unique_Patient_Identifier == "AGS_29" ~ FF_years + 9,
    Unique_Patient_Identifier == "AGS_45" ~ FF_years - 9,
    Unique_Patient_Identifier == "AGS_14" ~ FF_years - 8,
    Unique_Patient_Identifier == "AGS_52" ~ FF_years + 8,
    Unique_Patient_Identifier == "AGS_23" ~ FF_years + 13.75,
    TRUE ~ FF_years)) |> 
  mutate(new_y_text = case_when(
    Unique_Patient_Identifier == "AGS_08" ~ median_value + 0.75,
    Unique_Patient_Identifier == "AGS_43" ~ median_value -1.5,
    Unique_Patient_Identifier == "AGS_29" ~ median_value + 0,
    Unique_Patient_Identifier == "AGS_12" ~ median_value + 1,
    Unique_Patient_Identifier == "AGS_09" ~ median_value + 0.6,
    Unique_Patient_Identifier == "AGS_14" ~ median_value - 1,
    TRUE ~ median_value)) |> 
  mutate(new_x_seg = case_when(
    Unique_Patient_Identifier == "AGS_08" ~ FF_years + 1,
    Unique_Patient_Identifier == "AGS_17" ~ FF_years + 2,
    Unique_Patient_Identifier == "AGS_12" ~ FF_years + 3,
    Unique_Patient_Identifier == "AGS_09" ~ FF_years + 1,
    Unique_Patient_Identifier == "AGS_29" ~ FF_years + 3,
    Unique_Patient_Identifier == "AGS_45" ~ FF_years - 1.5 ,
    Unique_Patient_Identifier == "AGS_14" ~ FF_years - 1 ,
    Unique_Patient_Identifier == "AGS_52" ~ FF_years + 1 ,
    Unique_Patient_Identifier == "AGS_23" ~ FF_years + 1 ,
    TRUE ~ FF_years)) |> 
  mutate(new_y_seg = case_when(
    Unique_Patient_Identifier == "AGS_08" ~ median_value + 0.5,
    Unique_Patient_Identifier == "AGS_43" ~ median_value -.8,
    Unique_Patient_Identifier == "AGS_29" ~ median_value + 0,
    Unique_Patient_Identifier == "AGS_12" ~ median_value + 0.5,
    Unique_Patient_Identifier == "AGS_09" ~ median_value + .5,
    Unique_Patient_Identifier == "AGS_14" ~ median_value - .5,
    TRUE ~ median_value))


variant_fig_median <- trinuc_raw_means_cov_alljob_interest |> 
  ggplot(aes(x=FF_years, y=median_value)) + 
  # geom_label(aes(color=variant_inclusion),alpha=0,size = 0,label="") + 
  geom_text(data = variant_labels, 
            aes(label = variant_label,
                color=variant_inclusion,x=new_x_text,y=new_y_text
            ),fontface = "bold",size=5,key_glyph = "rect") +
  
  geom_segment(data=variant_labels, 
               aes(x= FF_years, xend = new_x_seg,y=median_value,yend=new_y_seg,
                   color=variant_inclusion),key_glyph = "rect") + 
  geom_point(aes(shape=dxcode_fct),size=4,alpha=0.5) +
  # geom_jitter(aes(shape=dxcode_fct),size=3,alpha=0.5,height = 0) + 
  # geom_text_repel(aes(label = variant_label,color=variant_inclusion),
  #                 min.segment.length = unit(0, 'lines'),
  #                 box.padding = .2,nudge_x = .8,nudge_y = -1,
  #                 key_glyph = "rect",fontface = "bold",size=5) + 
  

  
  labs(y="Median number of mutations attributable to SBS42", 
       x= "Firefighter years",color="Variant inclusion criteria",shape = "Histologic diagnosis") +
  theme_bw() + 
  scale_color_manual(na.translate=F,values=c("#1b9e77","#d95f02")) + 
  scale_shape_manual(values = c(16, 17, 15, 3,7,8,18))

# variant_fig





