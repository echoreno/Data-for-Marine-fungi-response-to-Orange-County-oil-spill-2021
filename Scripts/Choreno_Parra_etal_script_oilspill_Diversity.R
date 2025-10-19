#!/usr/bin/env Rscript

# R script
# Analysis of fungal ITS2 data for marine fungi after the Orange County oil spill, 2021.

# This script contains R code for the analysis of alpha (richness) and beta diversity, trend analysis,
# comparison of relative abundances of genera, and analysis of indicator OTUs

# R version: 4.3.0

# --Reproducibility note-- 
# This script illustrates the analysis workflow used in the study.
# During code cleaning from an earlier version, some details changed 
# (e.g., object naming, sample IDs, ordering, and random seeds).
# Therefore, running this cleaned script may not reproduce the original
# results exactly (e.g., identical OTUs, identical ordination coordinates, 
# or the exact p-values), though the qualitative conclusions are unchanged.


library(ggrepel)
library(BiocManager)
library(phyloseq)
library(vegan) 
library(tidyverse)
library(ggpubr)
library(devtools)
library(randomcoloR)
library(microbiome)
library(pairwiseAdonis)
library(reshape2)
library(RColorBrewer)
library(lme4)
library(ggplot2)
library(nlme)
library(MuMIn)
library(forcats)
library(ggalluvial)
library(ggplotify)
library(gtable)
library(nlme)
library(readxl)
library(tibble)
library(trend)
library(indicspecies)


# ================================================================================================================ #
#### 1. Import tables ####

# NOTE: This script uses OTU tables and taxonomy assignation obtained with the AMPtk pipeline, after filtering only Fungi.
# This information is contained in "Spreadsheet_S2.xlsx". 

# To obtain the original OTU tables and taxonomy assignation,
# run the AMPtk pipeline (in "Choreno_Parra_etal_script_oilspill_AMPtk.sh") 
# on raw sequences (available from NCBI’s SRA under BioProject ID PRJNA1344880).  
# Then import the "outputs" as a phyloseq object and filter only OTUs that are Fungi.


# Import metadata

raw_metadata <- read_excel("Spreadsheet_S1.xlsx", skip = 1) %>%
  slice_head(n = -1)

metadata_tbl <- raw_metadata %>%
  select(sample_name = "Sample name",
    sample_point = "Sample point",
    oiling_degree = "Oiling degree",
    time_days = "Time (days)")

metadata_tbl <- metadata_tbl %>%
  mutate(sample_name = as.character(sample_name),
    sample_point = as.factor(sample_point),
    oiling_degree = as.factor(oiling_degree),
    time_days = as.numeric(time_days),
    time_days_categoric = factor(case_when(
      time_days == 18 ~ "18d",
      time_days == 34 ~ "34d",
      time_days == 49 ~ "49d",
      time_days == 66 ~ "66d"),
      levels = c("18d", "34d", "49d", "66d"),
      ordered = TRUE))

# Import OTU table and taxonomy

# Spreadsheet_S2 results from filtering for "Fungi" in the original OTU table obtained from the AMPtk pipeline

otu_table_taxonomy <- read_excel("Spreadsheet_S2.xlsx")
otu_table_taxonomy <- otu_table_taxonomy %>%
  rename(OTUID = "OUT ID",
         growthForm = "Growth form")

taxonomy_tbl <- otu_table_taxonomy %>%
  select(OTUID, Kingdom, Phylum, Class, Order, Family, Genus, Species, growthForm)

otu_tbl <- otu_table_taxonomy %>% select(OTUID,
    "1_T1_1L", "2_T1_1L", "3_T1_2L", "4_T1_2L", "5_T1_3L", "6_T1_3L", "7_T1_4L", "8_T1_4L", "9_T1_5VL", "10_T1_5VL",
    "1_T2_1L", "2_T2_1L", "3_T2_2L", "4_T2_2L", "5_T2_3L", "6_T2_3L", "7_T2_4L", "8_T2_4L", "9_T2_5VL", "10_T2_5VL",
    "1_T3_1L", "2_T3_1L", "3_T3_2L", "4_T3_2L", "5_T3_3L", "6_T3_3L", "7_T3_4L", "8_T3_4L", "9_T3_5VL", "10_T3_5VL",
    "1_T4_1L", "2_T4_1L", "3_T4_2L", "4_T4_2L", "5_T4_3L", "6_T4_3L", "7_T4_4L", "8_T4_4L", "9_T4_5VL", "10_T4_5VL")

otu_tbl <- otu_tbl %>%
  column_to_rownames("OTUID") %>%
  as.matrix()




# ================================================================================================================ #
#### 2. Convert to phyloseq ####


# Phyloseq OTU table
otu_ps <- otu_table(otu_tbl, taxa_are_rows = TRUE)

# Phyloseq taxonomy table
tax_mat <- taxonomy_tbl %>%
  column_to_rownames("OTUID") %>%
  as.matrix()
tax_ps <- tax_table(tax_mat)

# Phyloseq metadata
metadata_ps <- metadata_tbl %>%
  column_to_rownames("sample_name") %>%
  sample_data()

# Merge phyloseq tables
fungi_physeq <- phyloseq(otu_ps, tax_ps, metadata_ps)
fungi_physeq




# ================================================================================================================ #
#### 3. Filtering marine fungi ####

marine_physeq <- subset_taxa(fungi_physeq, 
                             growthForm %in% c("Microfungus","Yeast",NA,"Facultative Yeast","NULL","Dimorphic Yeast",
                                               "Dimorphic-Facultative Yeast","Dimorphic-Microfungus","Dimorphic",
                                               "Facultative Yeast-Yeast","Dimorphic-Yeast","Dimorphic Yeast-Yeast",
                                               "Thallus","Chytrid-Microfungus","Dimorphic-Facultative Yeast-Microfungus-Tremelloid",
                                               "Dimorphic-Dimorphic Yeast","Dimorphic-Microfungus; Facultative Yeast (Tedersoo et al. 2014)",
                                               "Microfungus-Yeast","Facultative Yeast-Microfungus","Microfungus-Thallus"))




# ================================================================================================================ #
#### 4. Rarefying  ####


# curve of raw reads
otu.rare = otu_table(marine_physeq)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)
otu.rarecurve = rarecurve(otu.rare, step = 10000, label = T, cex = 0.5)
otu.rarecurve

# Decide reads for rarefaction
sort(sample_sums(marine_physeq))

# Final rarefication:
marine_physeq_rar <- rarefy_even_depth(marine_physeq, sample.size = 1185,    
                                       rngseed = 1744, replace = TRUE, 
                                       trimOTUs = TRUE, verbose = TRUE)

marine_physeq_rar #650 = 0.061, #1040 = 0.059, #1744 = 0.062

# curve of rarefied reads
otu.rare = otu_table(marine_physeq_rar)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)
otu.rarecurve = rarecurve(otu.rare, step = 10000, label = T)

marine_physeq_rar # FINAL TABLE FOR ANALYSES




# ================================================================================================================ #
#### 5. Community composition ####

# This part of the code averages counts per block, ranks abundances at each block
# and calculates Bray-Curtis distances


##### 5.1 Get average and rank #####

# Create blocks in metadata (sample_point × time_days)
meta <- as.data.frame(sample_data(marine_physeq_rar))
meta$sample_name <- rownames(meta)
meta$block <- paste(meta$sample_point, meta$time_days, sep = "_")

# OTU table to data frame (samples as rows)
otu_df <- as.data.frame(t(otu_table(marine_physeq_rar))) %>%
  rownames_to_column("sample_name")

# Merge metadata and average subsamples per block
otu_block_avg <- left_join(otu_df, meta[, c("sample_name", "block")], by = "sample_name") %>%
  select(-sample_name) %>%
  group_by(block) %>%
  summarise(across(everything(), ~ mean(as.numeric(.), na.rm = TRUE)))

otu_block_avg <- otu_block_avg %>%
  column_to_rownames("block")

# Transpose back (taxa as rows, blocks as columns)
otu_block_avg <- as.data.frame(t(otu_block_avg))

# Rank OTUs within each block
otu_block_rank <- apply(otu_block_avg, 2, rank, ties.method = "average")

# Build metadata for blocks
meta_block <- meta %>%
  group_by(block) %>%
  slice(1) %>%
  ungroup() %>%
  as.data.frame() %>%
  column_to_rownames("block")

# Reuse taxonomy and build phyloseq object
tax_mat <- as(tax_table(marine_physeq_rar), "matrix")[rownames(otu_block_rank), , drop = FALSE]
marine_physeq_avg_ranked <- phyloseq(otu_table(otu_block_rank, taxa_are_rows = TRUE),
                          sample_data(meta_block), 
                          tax_table(tax_mat))

# Calculate Bray–Curtis distances and ordination
marine_bray_avg_ranked <- phyloseq::distance(marine_physeq_avg_ranked, method = "bray", weighted = FALSE)


##### 5.2.NMDS #####

set.seed(600)
marine_nmds_bray_avg_ranked = ordinate(marine_physeq_avg_ranked, method="NMDS", distance=marine_bray_avg_ranked)

# Color palette for plots
oil_palette <- colorRampPalette(c("gray40", rev(brewer.pal(4, "Blues")),"#1dd3b0"))(4)


##### 5.3 Figure 2A #####

Figure2A <- plot_ordination(marine_physeq_avg_ranked,marine_nmds_bray_avg_ranked) +
  geom_point(aes(fill = time_days_categoric, shape = oiling_degree),
             size = 3.3, stroke = 0.3, color = "black") + 
  scale_fill_manual(values = oil_palette, name = "Time") +
  scale_shape_manual(values = c("Light" = 21, "Very_light" = 24),
                     labels = c("Light" = "Light", "Very_light" = "Very light"),
                     name = "Oiling degree") +
  labs(x = "NMDS 1", y = "NMDS 2") + 
  annotate("text", x = Inf, y = Inf,
           label = paste("Stress =", formatC(marine_nmds_bray_avg_ranked$stress, digits = 2, format = "f")),
           hjust = 3.7, vjust = 2, size = 3, color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3.3)),
         shape = guide_legend(override.aes = list(fill = NA, size = 3.3))) +
  guides(color = guide_legend(override.aes = list(size = 3.3)),
         shape = guide_legend(override.aes = list(size = 3.3)))

plot(Figure2A)


##### 5.4 Permanova#####

meta_block$time_days_categoric <- as.factor(meta_block$time_days_categoric)
meta_block$sample_point <- as.factor(meta_block$sample_point)
meta_block$oiling_degree <- as.factor(meta_block$oiling_degree)

# Run model   

set.seed(238)
permanova_results <- adonis2(
  marine_bray_avg_ranked ~  oiling_degree + time_days_categoric + oiling_degree * time_days_categoric,
  data = meta_block, strata = NULL, permutations = 999, by = "term")
permanova_results 

# Pairwise comparisons for time days
pairwise.adonis2(marine_bray_avg_ranked ~ time_days_categoric,
                 data = meta_block, strata = NULL, perm = 999)




# ================================================================================================================ #
#### 6. Trend analysis ####

# Uses marine_physeq_avg_ranked (block-averaged and ranked) and marine_bray_avg_ranked (Bray–Curtis distances)


##### 6.1. Prepare data #####

# Melt distance matrix (upper triangle only)
bray_matrix <- as.matrix(marine_bray_avg_ranked)
bray_long <- as.data.frame(as.table(bray_matrix)) %>%
  transmute(
    Sample1    = as.character(Var1),
    Sample2    = as.character(Var2),
    BrayCurtis = as.numeric(Freq)) %>%
  filter(Sample1 < Sample2)   # keep the upper triangle (unique pairs)

# Get metadata for both samples
meta_trend <- as(sample_data(marine_physeq_avg_ranked), "data.frame") %>%
  rownames_to_column("block")

bray_long <- bray_long %>%
  left_join(meta_trend,  by = c("Sample1" = "block")) %>%
  rename(Oiling1 = oiling_degree, TimeCat1 = time_days_categoric) %>%
  left_join(meta_trend,  by = c("Sample2" = "block")) %>%
  rename(Oiling2 = oiling_degree, TimeCat2 = time_days_categoric)

# Keep Light vs Very_light pairs at the same timepoint
bray_pairs <- bray_long %>%
  filter(((Oiling1 == "Very_light" & Oiling2 == "Light") |
            (Oiling1 == "Light"      & Oiling2 == "Very_light"))) %>%
  filter(TimeCat1 == TimeCat2)

# Summarize per timepoint
bray_summary <- bray_pairs %>%
  group_by(TimeCat1) %>%
  summarise(mean_bray = mean(BrayCurtis),
    sd_bray   = sd(BrayCurtis),
    se_bray   = sd_bray / sqrt(n()),
    n_pairs   = n(),
    .groups   = "drop") %>%
  arrange(factor(TimeCat1, levels = c("18d","34d","49d","66d")))


##### 6.2. Trend analysis #####

# Mann–Kendall trend test on the means
mk_results <- mk.test(bray_summary$mean_bray)
tau_value <- round(mk_results$estimates[3], 3)
p_value   <- round(mk_results$p.value, 3)


##### 6.3. Figure 2B ##### 

Figure2B <- ggplot(bray_summary, aes(x = as.numeric(gsub("d","", TimeCat1)), y = mean_bray)) +
  geom_jitter(data = bray_pairs,
              aes(x = as.numeric(gsub("d","", TimeCat1)),
                  y = BrayCurtis, fill = TimeCat1, color = TimeCat1),
              width = 1.5, height = 0, size = 1.5, alpha = 0.7,
              shape = 23, stroke = 0.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_bray - se_bray, ymax = mean_bray + se_bray),
                width = 1, color = "black") +
  geom_line(color = "black", size = 0.2) +
  geom_point(aes(fill = TimeCat1),
             size = 3, shape = 23, color = "black", stroke = 0.5) +
  scale_fill_manual(values = oil_palette) +
  scale_color_manual(values = oil_palette) +
  scale_x_continuous(breaks = c(18,34,49,66)) +
  scale_y_continuous(limits = c(0.0, 0.15)) +
  annotate("text", x = Inf, y = Inf,
           label = bquote(italic("Tau") ~ "=" ~ .(tau_value) * "," ~ italic("P") ~ "=" ~ .(p_value)),
           color = "black", size = 3, hjust = 1.05, vjust = 1.5) +
  labs(x = "Time after oil spill (days)",
       y = "Mean dissimilarity\nbetween oiling degrees",
       fill = "Time") +
  guides(color = "none",
         fill = guide_legend(override.aes = list(shape = 21, color = NA))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"),
        legend.position = "right",
        legend.margin = margin(0,0,0,0))

plot(Figure2B)




#### 7. Richness #####

# Obtain alpha diversity (per sample) and join metadata

# Alpha diversity with sample IDs as column
alpha_df <- estimate_richness(marine_physeq_rar, measures = "Observed")
alpha_df$sample_name <- sample_names(marine_physeq_rar)

# Metadata with sample IDs as column
meta_alpha <- as(sample_data(marine_physeq_rar), "data.frame") %>%
  rownames_to_column("sample_name")

alpha_meta <- left_join(alpha_df, meta_alpha, by = "sample_name")

alpha_meta$sample_point  <- as.factor(alpha_meta$sample_point)
alpha_meta$oiling_degree <- factor(alpha_meta$oiling_degree, levels = c("Light","Very_light"))
alpha_meta$time_days     <- as.numeric(alpha_meta$time_days)

##### 7.1. Mixed-effects model #####

# This model accounts for non-independence among subsamples and repeated measures

marine_model_richness = lme(Observed ~ time_days + oiling_degree + time_days * oiling_degree, 
                            random = ~1|sample_point/time_days,
                            method = "REML", data = alpha_meta)

summary(marine_model_richness) 


##### 7.2. Figure 3 #####

Figure3 <- ggplot(alpha_meta,
                  aes(x = time_days, y = Observed)) +
  geom_point(aes(fill = time_days_categoric,shape=oiling_degree), size = 3.3, alpha = 1, na.rm = T)+
  geom_point(aes(color = time_days_categoric,shape = oiling_degree), size = 3.3, alpha = 1, na.rm = T)+
  geom_point(aes(shape = oiling_degree), size = 3.3, alpha = 1, na.rm = T,color="black")+
  scale_color_manual(values = oil_palette) + 
  scale_fill_manual(values = oil_palette) + 
  scale_shape_manual(values = c("Light" = 21, "Very_light" = 24),
                     labels = c("Light" = "Light", "Very_light" = "Very light"))+
  scale_x_continuous(breaks=seq(0,140,20))+
  labs(x="Time after oil spill (days)", y="Richness", fill="Time", color="Time", shape="Oiling degree")+
  guides(color = "none", fill = "none") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=11,face="bold"), 
        axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size=9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"),
        legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt"),
        legend.box.margin = margin(t = -5, b = 5, l = 0, r = 0))

plot(Figure3)




#### 8. Relative abundance of genera #####


# Relative abundance for each OTU at each sample

Genus_abundance <- marine_physeq_rar %>%   
  tax_glom(taxrank = "Genus", NArm = FALSE) %>%
  transform_sample_counts(function(x) x / sum(x) * 100) %>%
  psmelt() %>%
  mutate(Genus = ifelse(is.na(Genus), "Unassigned", Genus)) %>% 
  select(Sample, Genus, Abundance, oiling_degree, time_days_categoric, sample_point)

# Make sure each samples sum 100%
Genus_abundance %>%
  group_by(Sample) %>%
  summarize(total_abundance = sum(Abundance)) %>%
  arrange(desc(abs(total_abundance - 1)))


##### 8.1. Per oiling degree #####
# Mean relative abundance of each genus per oiling degree

# Sum relative abundances of all the OTUs of the same Genus in each sample for each oiling 
Genus_mean_o <- Genus_abundance %>%  
  group_by(Genus, Sample, oiling_degree) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")%>% 
  group_by(oiling_degree)

# Summarizing
Genus_mean_o <- Genus_mean_o %>%
  group_by(Genus, oiling_degree) %>%
  summarize(mean_abundance = mean(Abundance),
            sd_abundance   = sd(Abundance, na.rm = TRUE),
            n              = n(),                           
            se_abundance   = sd(Abundance, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>% 
  group_by(oiling_degree) 

#Making sure abundances of each oiling level sum 100%
Genus_mean_o %>%
  group_by(oiling_degree) %>%
  summarise(total_abundance = sum(mean_abundance)) %>%
  arrange(oiling_degree)  

# Table for mean relative abundance of each genus per oiling degree: 
Genus_mean_o


##### 8.2. Per time #####
# Mean relative abundance of each genus per time

# sum relative abundances of all the OTUs of the same Genus in each sample for each oiling × time
Genus_mean_t <- Genus_abundance %>%  
  group_by(Genus, Sample, time_days_categoric) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop") %>% 
  group_by(time_days_categoric)

# Summarizing
Genus_mean_t <- Genus_mean_t %>%
  group_by(Genus, time_days_categoric) %>%
  summarize(mean_abundance = mean(Abundance),
            sd_abundance   = sd(Abundance, na.rm = TRUE),
            n              = n(),                               
            se_abundance   = sd(Abundance, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  group_by(time_days_categoric)

# Making sure abundances within each time sum 100% 
Genus_mean_t %>%
  group_by(time_days_categoric) %>%
  summarise(total_abundance = sum(mean_abundance)) %>%
  arrange(time_days_categoric)

# Table for mean relative abundance of each genus per time:  
Genus_mean_t


##### 8.3. Per oiling x time #####
# Mean relative abundance of each genus per oiling degree x time

# Sum relative abundances of all the OTUs of the same Genus in each sample for each oiling × time
Genus_mean_ot <- Genus_abundance %>%  
  group_by(Genus, Sample, oiling_degree, time_days_categoric) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop") %>% 
  group_by(oiling_degree, time_days_categoric)

# Summarizing
Genus_mean_ot <- Genus_mean_ot %>%
  group_by(Genus, oiling_degree, time_days_categoric) %>%
  summarize(mean_abundance = mean(Abundance),
            sd_abundance   = sd(Abundance, na.rm = TRUE),
            n              = n(),                               
            se_abundance   = sd(Abundance, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>% 
  group_by(oiling_degree, time_days_categoric)

# Making sure abundances within each oiling × time sum 100% 
Genus_mean_ot %>%
  group_by(oiling_degree, time_days_categoric) %>%
  summarise(total_abundance = sum(mean_abundance)) %>%
  arrange(oiling_degree, time_days_categoric)

# Table for mean relative abundance of each genus per oiling degree x time:  
Genus_mean_ot


##### 8.4. Mixed effects models #####

# Use Genus_abundance
length(unique(Genus_abundance$Genus)) # Number of genera

# Ranking abundances in long format

# Collapse to one row per Sample x Genus
Genus_abundance_collapsed <- Genus_abundance %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

# Merge metadata to abundances
sample_meta <- Genus_abundance %>%
  distinct(Sample, oiling_degree, time_days_categoric, sample_point)

Genus_abundance_collapsed <- Genus_abundance_collapsed %>%
  left_join(sample_meta, by = "Sample")

# Rank within each subsample (lowest = 1, highest = max)
Genus_ranked_abundance <- Genus_abundance_collapsed %>%
  group_by(Sample) %>%
  mutate(Abundance_rank = rank(Abundance, ties.method = "average")) %>%
  ungroup() %>%
  mutate(sample_point  = as.factor(sample_point),
         time_days_categoric = as.factor(time_days_categoric),
         oiling_degree = as.factor(oiling_degree),
         Genus         = as.factor(Genus))


###### 8.4.1. Model Genus x time x oiling ######

# Fit the model
model_ranked_Genus <- nlme::lme(
  Abundance_rank ~ Genus * oiling_degree * time_days_categoric,
  random = ~1 | sample_point/time_days_categoric,
  data   = Genus_ranked_abundance,
  method = "REML")

## 4) ANOVA table
anova(model_ranked_Genus)


###### 8.4.2. Pairwise comparisons ######

# Among time within each genus 
pw_genus_t <- emmeans::emmeans(model_ranked_Genus, pairwise ~ time_days_categoric | Genus, adjust = "tukey")
pw_genus_t <- pairs(pw_genus_t, adjust = "tukey")
pw_genus_t <- as.data.frame(pw_genus_t)
pw_genus_t_sig <- pw_genus_t %>% filter(p.value < 0.05) # Get significant pw comparisons
pw_genus_t_sig

# Among oiling within each genus
pw_genus_o <- emmeans::emmeans(model_ranked_Genus, pairwise ~ oiling_degree | Genus, adjust = "tukey")
pw_genus_o <- pairs(pw_genus_o, adjust = "tukey")
pw_genus_o <- as.data.frame(pw_genus_o)
pw_genus_o_sig <- pw_genus_o %>% filter(p.value < 0.05)
pw_genus_o_sig

# Interaction: oiling differences at specific times (per Genus)
pw_genus_o_by_t <- emmeans::emmeans(model_ranked_Genus, 
                                    pairwise ~ oiling_degree | Genus * time_days_categoric, adjust = "tukey")
pw_genus_o_by_t <- pairs(pw_genus_o_by_t, adjust = "tukey")
pw_genus_o_by_t <- as.data.frame(pw_genus_o_by_t)
pw_genus_o_by_t_sig <- pw_genus_o_by_t %>% filter(p.value < 0.05)
pw_genus_o_by_t_sig


##### 8.5. Figure 4 #####

# Select the genera and set factor orders
Genus_sel <- Genus_mean_ot %>%
  filter(Genus %in% c("Mycosphaerella","Papiliotrema", "Rhodotorula", "Malassezia")) %>%
  mutate(time_days_categoric     = factor(time_days_categoric, levels = c("18d","34d","49d","66d")),
    oiling_degree = factor(oiling_degree, levels = c("Very_light","Light")),
    Genus         = factor(Genus, levels = c("Mycosphaerella","Papiliotrema", "Rhodotorula",
                                             "Malassezia")))

# Establish y-axis limits
g_max15 <- c("Mycosphaerella")
g_max5  <- c("Papiliotrema","Rhodotorula","Malassezia")

time_lvls <- levels(Genus_sel$time_days_categoric)
if (is.null(time_lvls)) time_lvls <- unique(Genus_sel$time_days_categoric)

limits_df <- tibble::tibble(Genus = c(g_max15, g_max5),
                            ymax  = c(rep(15, length(g_max15)), rep(4.5, length(g_max5)))) %>%
  expand_grid(time_days_categoric = time_lvls) %>%
  mutate(Genus     = factor(Genus,     levels = levels(Genus_sel$Genus)),
         time_days_categoric = factor(time_days_categoric, levels = time_lvls))

pd <- position_dodge(width = 0.25)

Figure4 <- ggplot(Genus_sel, aes(x = time_days_categoric, y = mean_abundance,
                                 group = oiling_degree, shape = oiling_degree, linetype = oiling_degree)) +
  geom_blank(data = distinct(limits_df, Genus, time_days_categoric) %>% 
      mutate(y0 = 0), inherit.aes = FALSE, aes(x = time_days_categoric, y = y0, Genus = Genus)) +
  geom_blank(data = limits_df, inherit.aes = FALSE, aes(x = time_days_categoric, y = ymax, Genus = Genus)) +
  geom_line(color = "black", size = 0.2) +
  geom_errorbar(aes(x = time_days_categoric, group = oiling_degree,
                    ymin = mean_abundance - se_abundance,
                    ymax = mean_abundance + se_abundance),
                position = pd, width = 0.2, color = "black",
                linetype = "solid", size = 0.2, inherit.aes = FALSE) +
  geom_point(aes(fill = time_days_categoric), position = pd,
             size = 2.5, stroke = 0.5, color = "black") +
  scale_shape_manual(values = c("Light" = 21, "Very_light" = 24),
                     labels = c("Light" = "Light", "Very_light" = "Very light"),
                     name   = "Oiling degree") +
  scale_linetype_manual(values = c("Light" = "dashed", "Very_light" = "solid"),
                        labels = c("Light" = "Light", "Very_light" = "Very light"),
                        name   = "Oiling degree") +
  scale_fill_manual(values = oil_palette, guide = "none") +
  facet_wrap(~ Genus, scales = "free_y") +
  labs(x = "Time after oil spill (days)", y = "Relative abundance (%)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title   = element_text(size = 11, face = "bold"),
    axis.text.y  = element_text(size = 9),
    axis.text.x  = element_text(size = 9),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 11, face = "bold"),
    strip.text.x = element_text(face = "bold.italic", size = 11, color = "black"),
    strip.background = element_blank(),
    legend.key.size  = unit(0.6, "cm"),
    legend.position  = "top",
    legend.margin    = margin(t = -5, b = 5),
    legend.box.margin= margin(b = -13))

plot(Figure4)




#### 9. Indicator species ####

# Extract abundance data from the phyloseq object and transpose it
otu_df_indic <- as.data.frame(t(otu_table(marine_physeq_rar)))

# Extract metadata 
meta_indic <- as(sample_data(marine_physeq_rar), "data.frame")

# Create variable using oiling degree + time = sample type
meta_indic <- meta_indic %>% mutate(sample_type = paste(oiling_degree, time_days_categoric, sep = "_"))

# Add the column to the transposed abundance data
otu_df_indic$sample_type <- meta_indic$sample_type
otu_df_indic$sample_type

# Create a data frame suitable for multipatt
multipatt_data <- otu_df_indic %>%
  select(-sample_type)

##### 9.1. Multipatt analysis #####

set.seed(500)
indicator_otus <- multipatt(multipatt_data, otu_df_indic$sample_type, control = how(nperm=999), duleg=T)
summary(indicator_otus, indvalcomp=T, alpha=0.05)


##### 9.2. Identity of indicator OTUs ####

# Extract per-OTU results
sig_tbl <- indicator_otus$sign %>%
  as.data.frame() %>%
  rownames_to_column("OTU")

# Remove "s." prefix from the sample_type columns
colnames(sig_tbl) <- sub("^s\\.", "", colnames(sig_tbl))

# Keep only significant OTUs (alpha = 0.05)
sig_tbl_sig <- sig_tbl %>%
  filter(p.value <= 0.05) %>%
  arrange(p.value)

# Pick the group columns
group_cols <- c("Light_18d","Light_34d","Light_49d","Light_66d",
                "Very_light_18d","Very_light_34d","Very_light_49d","Very_light_66d")

# Create table of indicator OTUs in long format
indic_long <- sig_tbl_sig %>%
  select(OTU, all_of(group_cols), p.value) %>%
  pivot_longer(cols = all_of(group_cols),
               names_to = "sample_type",
               values_to = "is_indicator") %>%
  filter(is_indicator == 1) %>%
  select(OTU, sample_type, p.value) %>%
  separate(sample_type,
    into = c("oiling_degree", "time_days_categoric"),
    sep = "_(?=\\d)", remove = FALSE) %>%
  arrange(p.value, OTU, sample_type) %>%
  rename(OTUID = OTU)

# Get indicator OTUs from phyloseq object
otu_ids_to_keep <- indic_long$OTUID 
indicator_otus_identity <- prune_taxa(otu_ids_to_keep, marine_physeq_rar)

# Extract taxonomy of indicator OTUs 
tax_indicator_otus <- as.data.frame(tax_table(indicator_otus_identity))

# Add OTUID column
tax_indicator_otus$OTUID <- rownames(tax_indicator_otus)

# Merge table of indicator OTUs with taxonomy
indic_long$OTUID <- as.character(indic_long$OTUID)
indicator_otus_table <- merge(indic_long, tax_indicator_otus, by = "OTUID", all.x = TRUE)

# Final taxonomy table of indicators per oiling x time
indicator_otus_table




#### 10. FIGURES ARTICLE ####


##### Figure 2 #####

Figure2 <-ggarrange(Figure2A,Figure2B,
          nrow=1,ncol=2,common.legend = T,
          labels=c("A","B"), font.label = list(size = 10))

ggsave("Figure2.tiff", plot = Figure2, units = "in", 
       height = 3, width = 7, dpi=600)


##### Figure 3 #####

ggsave("Figure3.tiff", plot = Figure3, units = "in", 
       height = 3, width = 3.3, dpi=600)


##### Figure 4 #####

ggsave("Figure4.tiff", plot = Figure4, units = "in", 
       height = 4.5, width = 6, dpi = 600)
