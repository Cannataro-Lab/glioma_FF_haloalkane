
# fig 1 recreate -----

trinuc_raw_df <- data.table::rbindlist(trinuc_raw,idcol = "boot_n")


mean_42_cesr <- trinuc_raw_df |> 
  select(-boot_n) |> 
  pivot_longer(-Unique_Patient_Identifier) |> 
  filter(name == "SBS42") |> 
  group_by(Unique_Patient_Identifier) |> 
  summarize(mean_value = mean(value),
            median_value = median(value),
            q2.5 = quantile(value,0.025),
            q97.5 = quantile(value,0.975))



FF_samples <- all_data_filtered_maf |> 
  filter(ff==T) |> 
  pull(Unique_Patient_Identifier) |> 
  unique()

nonFF_samples <- all_data_filtered_maf |> 
  filter(ff==F) |> 
  pull(Unique_Patient_Identifier) |> 
  unique()

mean_42_cesr <- mean_42_cesr |> 
    mutate(ff = case_when(
      Unique_Patient_Identifier %in%  FF_samples ~ "FF",
      Unique_Patient_Identifier %in% nonFF_samples ~ "nonFF"))

# mean_42_cesr <- mean_42_cesr |> 
#   mutate(ff = case_when(
#     Unique_Patient_Identifier %in%  all_data_ff$tumor_sample ~ "FF", 
#     Unique_Patient_Identifier %in% all_data_nonFF$tumor_sample ~ "nonFF"))

mean_42_cesr |> 
  left_join(cov_data, by = c("Unique_Patient_Identifier" = "AGSID")) |> 
  ungroup() -> 
  trinuc_raw_means_cov


trinuc_raw_means_cov_alljob <- trinuc_raw_means_cov |> 
  mutate(FF_years = case_when(
    is.na(Yrs_Firefighting_Exposure) ~ 0, 
    !is.na(Yrs_Firefighting_Exposure) ~ Yrs_Firefighting_Exposure)) |> 
  unite(all_jobs, starts_with(c("JOB","ACTIV")),sep = "\n",na.rm = T) |> 
  mutate(all_jobs = paste(Unique_Patient_Identifier, all_jobs))




trinuc_raw_means_cov_alljob_interest <- trinuc_raw_means_cov_alljob |> 
  mutate(ff_cases_lab = case_when(ff == "FF" ~ "Firefighter",
                                  ff == "nonFF" ~ "Not firefighter")) |> 
  mutate(occupation = case_when(
    Unique_Patient_Identifier == "AGS_43" ~ "Truck/Car Mechanic,\n Machine Maintenance,\nincluding diesel repair & maintenance",
    Unique_Patient_Identifier == "AGS_09" ~ "Farmer,\nMarine Engineer,\nSpray insecticide,\nInspect/repair ships, \nMilitary petroleum transport",
    Unique_Patient_Identifier == "AGS_17" ~ "Painted cars",
    Unique_Patient_Identifier == "AGS_12" ~ "Marine Engineer,\nBoilermaker/Pipefitter,\nMechanic "
    ))


# trinuc_raw_means_cov_alljob_interest_cesr <- trinuc_raw_means_cov_alljob_interest |> left_join(mean_42_cesr,by = "Unique_Patient_Identifier")


library(ggrepel)

job_fig <- trinuc_raw_means_cov_alljob_interest |> 
  ggplot(aes(x=FF_years, y=mean_value)) + 
  # geom_point(aes(shape = ff_cases_lab), size=2) + 
  geom_point(size=4,alpha=0.5) + 
  geom_text_repel(aes(label = occupation),
                  min.segment.length = unit(0, 'lines'),
                  box.padding = 0.4,nudge_x = 1,
                  fontface = "bold",size = 5) + 
  labs(y="Mean number of mutations attributable to SBS42", 
       x= "Firefighter years",color="Firefighter status",
       shape = "Firefighter status") +
  theme_bw() 


trinuc_raw_means_cov_alljob_interest |> select(Unique_Patient_Identifier,occupation) |> 
  filter(!is.na(occupation))


job_fig_median <- trinuc_raw_means_cov_alljob_interest |> 
  ggplot(aes(x=FF_years, y=median_value)) + 
  # geom_point(aes(shape = ff_cases_lab), size=2) + 
  geom_point(size=4,alpha=0.5) + 
  # geom_text_repel(aes(label = occupation),
  #                 min.segment.length = unit(0, 'lines'),
  #                 box.padding = 1,nudge_x = 1.2,
  #                 fontface = "bold",size = 5) + 
  geom_text(data = trinuc_raw_means_cov_alljob_interest |> 
              filter(Unique_Patient_Identifier == "AGS_17"),
            aes(label=occupation,x = FF_years +1),
            size=5,hjust = 0,fontface = "bold") + 
  geom_segment(data = trinuc_raw_means_cov_alljob_interest |> 
                 filter(Unique_Patient_Identifier == "AGS_17"),
               aes(x = FF_years +1,xend = FF_years),size=1) + 
  
  geom_text(data = trinuc_raw_means_cov_alljob_interest |> 
              filter(Unique_Patient_Identifier == "AGS_12"),
            aes(label=occupation,x = FF_years + 9.5, y=median_value +1),
            size=5,hjust = 0.5,fontface = "bold") + 
  geom_segment(data = trinuc_raw_means_cov_alljob_interest |> 
                 filter(Unique_Patient_Identifier == "AGS_12"),
               aes(x = FF_years +2,xend = FF_years,
                   y =median_value +1, yend=median_value),size=1) +
  
  
  geom_text(data = trinuc_raw_means_cov_alljob_interest |> 
              filter(Unique_Patient_Identifier == "AGS_09"),
            aes(label=occupation,x = FF_years + 12, y=median_value - 3),
            size=5,hjust = 0.5,fontface = "bold") + 
  geom_segment(data = trinuc_raw_means_cov_alljob_interest |> 
                 filter(Unique_Patient_Identifier == "AGS_09"),
               aes(x = FF_years + 6,xend = FF_years,
                   y =median_value - 3, yend=median_value),size=1) +
  
  
  geom_text(data = trinuc_raw_means_cov_alljob_interest |> 
              filter(Unique_Patient_Identifier == "AGS_43"),
            aes(label=occupation,x = FF_years + 10, y=median_value - 7),
            size=5,hjust = 0.5,fontface = "bold") + 
  geom_segment(data = trinuc_raw_means_cov_alljob_interest |> 
                 filter(Unique_Patient_Identifier == "AGS_43"),
               aes(x = FF_years + 2.5,xend = FF_years,
                   y =median_value - 6, yend=median_value),size=1) +
  
  labs(y="Median number of mutations attributable to SBS42", 
       x= "Firefighter years",color="Firefighter status",
       shape = "Firefighter status") +
  theme_bw() 

job_fig_median

# AGS40183

# trinuc_raw_df |> 
#   select(-boot_n) |> 
#   pivot_longer(-Unique_Patient_Identifier) |> 
#   filter(Unique_Patient_Identifier == "AGS40183") |> 
#   group_by(name) |> 
#   summarize(meanval = mean(value)) |> 
#   filter(meanval > 0)