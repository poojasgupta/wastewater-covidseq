# R-code associated with manuscript titled 
# "Wastewater genomic surveillance captures early detection of Omicron in Utah"
# Load required R-packages

library(dplyr); library(tidyverse); library(ggplot2)
library(ggpubr); library(ggthemes); library(gridExtra);
library(lubridate); library(tibbletime)
library(scales); library(viridis)
library(patchwork) ; library(broom.mixed)


setwd("/Users/pgupta/Library/CloudStorage/GoogleDrive-pgupta@utah.gov/My Drive/ww_manuscript/Mspectrum_analysis")

# Read data and format collection date
lin_abund_data <- read.csv('2022-06-14_WW_feyja_varaints_SC2_lineage_abundance_ediDates.csv')
lin_abund_data$collection_date <- ymd(lin_abund_data$collection_date)


# Filter the data for the dates in the study
lin_abund_data_filt <- lin_abund_data %>% filter(between(collection_date, as.Date('2021-11-02'), as.Date('2022-03-29'))) 

# Create a 'Parent_lineage' column to group sub-lineages and rename as WHO lineages
lin_abund_data_filt2 <- lin_abund_data_filt %>% 
  mutate(Parent_lineage = case_when(
    lineages == "A" ~ "A",
    lineages == "B" ~ "B",
    lineages == "B.1.1.7" ~ "Alpha",
    lineages == "B.1.351" ~ "Beta",
    lineages == "P.1" ~ "Gamma",
    lineages == "P.2" ~ "Zeta",
    lineages == "B.1.617.1" ~ "Kappa",
    lineages == "B.1.617.2" ~ "Delta",
    lineages == "B.1.427" ~ "Epsilon",
    lineages == "B.1.429" ~ "Epsilon",
    lineages == "B.1.525" ~ "Eta",
    lineages == "B.1.526" ~ "Iota",
    lineages == "B.1.621" ~ "Mu",
    lineages == "B.1.1.529" ~ "Omicron",
    startsWith(lineages, "BA.1") ~ "Omicron (BA.1)",
    startsWith(lineages, "BA.2") ~ "Omicron (BA.2)",
    startsWith(lineages, "BA.3") ~ "Omicron (BA.3)",
    startsWith(lineages, "BA.4") ~ "Omicron (BA.4)",
    startsWith(lineages, "BA.5") ~ "Omicron (BA.5)",
    startsWith(lineages, "Q") ~ "Alpha",
    startsWith(lineages, "AY.") ~ "Delta (AY)",
    startsWith(lineages, "B.1.") ~ "B.1",
    startsWith(lineages, "AC") ~ "B.1.1.405",
    startsWith(lineages, "AD") ~ "B.1.1.315",
    startsWith(lineages, "AE") ~ "B.1.1.306",
    startsWith(lineages, "AF") ~ "B.1.1.305",
    startsWith(lineages, "AH") ~ "B.1.1.241",
    startsWith(lineages, "AK") ~ "B.1.1.232",
    startsWith(lineages, "AS") ~ "B.1.1.317",
    startsWith(lineages, "AT") ~ "B.1.1",
    startsWith(lineages, "AZ") ~ "B.1.1.318",
    startsWith(lineages, "A.1.") ~ "A.1",
    startsWith(lineages, "A.2.") ~ "A.2",
    startsWith(lineages, "A.4.") ~ "A.4",
    startsWith(lineages, "A.5.") ~ "A.5",
    startsWith(lineages, "A.6.") ~ "A.6",
    startsWith(lineages, "A.7.") ~ "A.7",
    startsWith(lineages, "A.") ~ "A.",
    startsWith(lineages, "C.") ~ "B.1.1.1",
    startsWith(lineages, "D") ~ "B.1.1.25",
    startsWith(lineages, "L") ~ "B.1.1.10",
    startsWith(lineages, "M") ~ "B.1.1.216",
    startsWith(lineages, "N") ~ "B.1.1.33",
    startsWith(lineages, "P") ~ "B.1.1.28",
    startsWith(lineages, "R") ~ "B.1.1.316",
    startsWith(lineages, "B.1") ~ "B.1",
    startsWith(lineages, "B.2") ~ "B.2",
    startsWith(lineages, "X") ~ "Recombinant",
    startsWith(lineages, "proposed") ~ "Unassigned",
    startsWith(lineages, "miscRecombOrContam") ~ "Unassigned",
    TRUE ~ "NA"))

write.csv(lin_abund_data, file = "2022-06-14_WW_feyja_varaints_SC2_lineage_abundance_ediDates_addParentLin_20230320.csv")
# This however did not work as expected so after this step, the file was modified externally in excel
# to correct the grouping of specifically B.1.1 and B.1 lineages
# and the file was uploaded back

# Read data and format collection date column again
new_lin_abund_data <- read.csv('2022-06-14_WW_feyja_varaints_SC2_lineage_abundance_ediDates_addParentLin_corrected20230322.csv')
new_lin_abund_data <- transform(new_lin_abund_data,collection_date=as.Date(as.character(collection_date),"%Y-%m-%d"))

# Re-filtering the dates again, just in case and remove unnecessary columns
new_lin_abund_data_filter <- new_lin_abund_data %>% 
  filter(between(collection_date, as.Date('2021-11-02'), as.Date('2022-03-29'))) %>% 
  select(-lin_grp, -idx_name)

##-----------------------------------------------------------------------
# Creating second level of grouping lineages by 
# creating another column 'Parent_lineage_grp'
# This df is used to make the detailed lineage plot in the Supplementary
##-----------------------------------------------------------------------

new_lin_abund_data_filter <- new_lin_abund_data_filter %>% 
   mutate(Parent_lineage_grp = case_when(
    Parent_lineage == "A" ~ "A",
    Parent_lineage == "B" ~ "B",
    Parent_lineage == "Delta" ~ "Delta",
    Parent_lineage ==  "Alpha" ~ "Other WHO lineages",
    Parent_lineage == "Beta" ~ "Other WHO lineages",
    Parent_lineage == "Gamma" ~ "Other WHO lineages",
    Parent_lineage == "Theta" ~ "Other WHO lineages",
    Parent_lineage == "Zeta" ~ "Other WHO lineages",
    Parent_lineage == "Kappa" ~ "Other WHO lineages",
    Parent_lineage == "B.1.617.2" ~ "Delta",
    Parent_lineage == "Epsilon" ~ "Other WHO lineages",
    Parent_lineage == "Eta" ~ "Other WHO lineages",
    Parent_lineage == "Iota" ~ "Other WHO lineages",
    Parent_lineage == "Mu" ~ "Other WHO lineages",
    Parent_lineage == "B.1" ~ "B.1*",
    Parent_lineage == "B.1.1" ~ "B.1.1*",
    Parent_lineage == "B.1.1.25" ~ "B.1.1*",
    Parent_lineage == "B.1.1.33" ~ "B.1.1*",
    Parent_lineage == "B.1.1.7" ~ "B.1.1*",
    Parent_lineage == "B.1.1.10" ~ "B.1.1*",
    Parent_lineage == "B.1.1.1" ~ "B.1.1.1*",
    Parent_lineage =="Recombinant" ~ "Recombinant",
     Parent_lineage =="miscRecombOrContam" ~ "Unassigned",
    Parent_lineage == "other" ~ "Unassigned",
    startsWith(Parent_lineage, "Omicron") ~ "Omicron",
     TRUE ~ "NA"))

# Summarize the data
new_lin_abund_data_filter_summ <- new_lin_abund_data_filter %>% filter(lineages != "None") %>%
  group_by(collection_date) %>%  count(Parent_lineage_grp)

# Get color-blind friendly colors 
cbPalette <- c("#8dd3c7", "#fb9a99","#bc80bd","#fccde5" ,"#CC69A2", "#80b1d3", "#fdb462",  
 "#ffffb3", "#d9d9d9", "#999999",  "#ccebc5", "#ffed6f")

# Detailed lineage plot as in Supplementary
ww_p1 <- ggplot(new_lin_abund_data_filter_summ , aes(x=collection_date, y=n, fill= Parent_lineage_grp))+ 
  geom_bar(position="fill", stat="identity", color="light grey", linewidth=0.05) + 
  labs(y = "Proportion of SARS-CoV-2 lineages", x = "Sample collection date")  + theme_pubr() + 
   scale_fill_manual(values = darker_cbPalette ) + 
  theme(legend.title=element_blank(), legend.text=element_text(size=10),
    axis.text=element_text(size=10),  
    axis.title=element_text(size=11), plot.title = element_text(face="bold"),
        panel.grid.major.x = element_line(color = "dark grey", linewidth = 0.2))+ 
  scale_x_date(date_labels = "%d %b %Y",   date_breaks = "2 weeks")+ rotate_x_text(angle = 45, hjust = 0.8, vjust = 0.8) 

ggsave(filename = "lin_prop_ww_plot_detailed_2023-03-22.tiff", plot = ww_p1, width = 6.8, height = 4.8, dpi = 600, units = "in")

##-----------------------------------------------------
# Creating third level of grouping lineages by creating 
# another column 'Category'
# This is used for creating lineage plots in the main text
##-----------------------------------------------------

new_lin_abund_data_filter_cat <- new_lin_abund_data_filter %>% 
   mutate(Category = case_when(
    Parent_lineage == "Delta" ~ "Delta",
    startsWith(Parent_lineage, "Omicron") ~ "Omicron",
    startsWith(Parent_lineage, "Other") ~ "Other",
    TRUE ~ "Other"))
write.csv(new_lin_abund_data_filter_cat , file = "2022-06-14_WW_feyja_varaints_SC2_lineage_abundance_ediDates_addParentLin_corrected_filter_date_final.csv")

# Summarize the data
new_lin_abund_data_filter_cat_summ <- new_lin_abund_data_filter_cat %>%  
  group_by(collection_date) %>%  count(Parent_lineage, Category)

# Make all lineage bar plot for Delta and Omicron
plotA <- ggplot(new_lin_abund_data_filter_cat_summ , aes(x=collection_date, y=n, fill= Category))+ 
  geom_bar(position="fill", stat="identity", color="light grey", size=0.05) + 
  labs(title = "A", y = "Proportion of SARS-CoV-2 lineages", x="")  + theme_pubr() + 
  scale_fill_viridis_d(option="viridis", alpha=0.8)+
  theme(legend.title=element_blank(), legend.text=element_text(size=8), legend.position = "top",
    axis.text=element_text(size=7),  axis.title=element_text(size=9),
    plot.title=element_text(size=10), legend.key.size = unit(0.4, 'cm'),
    panel.grid.major.x = element_line(color = "dark grey", linewidth = 0.2), 
  plot.margin = unit(c(0, 0.5, 0, 0), "cm"))+
  scale_x_date(date_labels = "%d %b %Y",   date_breaks = "2 weeks")+ rotate_x_text(angle = 45, hjust = 0.8, vjust = 0.8)

##-----------------------------------------------------
# Omicron specific plot
# Filter out only Omicron lineages
new_lin_abund_data_filter_cat_summ2 <- new_lin_abund_data_filter_cat_summ %>% 
filter(str_detect(Category, 'Omicron'))
##-----------------------------------------------------

# Make bar plot for Omicron
plotC <- ggplot(new_lin_abund_data_filter_cat_summ2 , aes(x=collection_date, y=n, fill= Parent_lineage)) +
  geom_bar(position="fill", stat="identity", color="light grey", size=0.05) + 
  labs(title = "C",y = "Proportion of Omicron lineages", x="")  + theme_pubr() + 
  theme(legend.title=element_blank(), legend.text=element_text(size=8), legend.position = "top",
    axis.text=element_text(size=7),  axis.title=element_text(size=9),
    plot.title=element_text(size=10), legend.key.size = unit(0.4, 'cm'),
    panel.grid.major.x = element_line(color = "dark grey", linewidth = 0.2),
    plot.margin = unit(c(0, 0.5, 0, 0), "cm"))+
  scale_fill_manual( values = c("Omicron" = "#9FDA3AFF", "Omicron BA.1" = "#277F8EFF", 
    "Omicron BA.2" = "#4AC16DFF", "Omicron BA.3" = "#1FA187FF" ), 
    labels=c("Omicron", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.3)"))+
  rotate_x_text(angle = 45, hjust = 0.8, vjust = 0.8) + scale_x_date(date_labels = "%d %b-%Y", date_breaks = "2 weeks")

# Merge with clinical plots from 'clinical_plots.R'
plot_top <- ggarrange(plotA, plotB, common.legend = TRUE, align = "h")
plot_bottom <- ggarrange(plotC, plotD, common.legend = TRUE, align = "h")

ww_clinical_p1 <- plot_top + plot_bottom + plot_layout(ncol=1)

# Formatting layout for mspectrum
ww_clinical_p2 <- patchwork::patchworkGrob(ww_clinical_p1)
ww_clinical_p3 <- gridExtra::grid.arrange(ww_clinical_p2, bottom = "Sample collection date" )

ww_clinical_p3 <- gridExtra::grid.arrange(ww_clinical_p2, bottom = grid::textGrob('Sample collection date', vjust = -1.8, gp=grid::gpar(fontsize=9)))

# Save plot 
ggsave(filename = "lin_prop_ww_clinical_plot_2023-03-22.tiff", plot = ww_clinical_p3 , 
width = 6.875, height = 6, dpi = 600, units = "in") #seems to be more useful

#########################################################################
# Omicron bubble plot showing lineage abundance data
#########################################################################

#Reading in the lineage data after adding a 'Compress lineage' column externally
lin_data_edi <- read.csv('2022-06-14_WW_feyja_varaints_SC2_lineage_abundance_ediDates_addParentLin_corrected_filter_date_final_forMS.csv')

# Format collection date
lin_data_edi$Collection_date <- ymd(lin_data_edi$Collection_date)

# Filter data for Omicron lineages and samples collected in December and January
jan.dec.subset <- lin_data_edi %>% filter(Parent_lineage_grp =='Omicron') %>%
  filter(between(Collection_date, as.Date('2021-12-03'), as.Date("2022-01-28")))

jan.dec.subset$Lineage_abundance <- round(jan.dec.subset$Lineage_abundance, 2)
  
# Plot bubble plot
set.seed(101441)
omicron_bubblep_jd_edi <- jan.dec.subset %>%
  ggplot(aes(x=Collection_date, y=Compress_lineages, size=Lineage_abundance, fill=WW_site_shrtnm)) +
  geom_point(alpha=0.5, shape=21, color="black", position = position_jitter()) +
  scale_size_area(
    max_size = 5, breaks=c(0.2,0.4,0.6,0.8),
    labels=c("0.20","0.40","0.60","0.80"),
    name = "Lineage abundance", guide_legend(),
    limits = c(0.2, 1), oob = scales::squish) +
  labs(y = "Omicron lineages", x = "Sample collection date") + theme_pubr(border=TRUE) +
  scale_fill_viridis_d(option="magma", 
    guide = guide_legend(override.aes = list(size=4,alpha=0.6), title="Wastewater sampling sites", 
    ncol=3, title.position = "top"))+
  theme(axis.title=element_text(size=11, face = "bold"), axis.text=element_text(size=8),
        legend.text=element_text(size=8), legend.key.size = unit(0.4, 'cm'), 
        legend.title=element_text(size=10), legend.position = "bottom",
        panel.grid.major.x = element_line(color = "light grey", linewidth = 0.2),
  plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm"))+
  rotate_x_text(angle = 45, hjust = 0.8, vjust = 0.8) + scale_x_date(date_labels = "%d %b-%Y",   date_breaks = "4 days")+
  theme(legend.direction = "vertical", legend.box = "horizontal")

#Save plot
ggsave(filename = "Omicron_abundance_jan_dec_2023-03-23.tiff", plot= omicron_bubblep_jd_edi, 
  width = 6.875, height = 8, dpi = 600, units = "in") 