# Analyses of shotgun metagenomic sequencing data for Tanaka et al. 
# Statistician: Matthew Kelly, MD, MPH

set.seed(1234)
library(phyloseq)
library(tidyverse)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

remove(list=ls())
# Need to set working directory
setwd("____________") 

phy.gvhd <- readRDS("phy.gvhd.rds")
metadata_gvhd <- data.frame(sample_data(phy.gvhd))

phy.gvhd <- transform_sample_counts(phy.gvhd, function(OTU) round(OTU))
sample_sums(phy.gvhd)
# 72 samples from 36 participants
nsamples(phy.gvhd)
ntaxa(phy.gvhd)
# 284 taxa
phy.genera <- tax_glom(phy.gvhd, taxrank = 'Genus')
ntaxa(phy.genera)
# 116 genera
phy.orders <- tax_glom(phy.gvhd, taxrank = 'Order')
ntaxa(phy.orders)
# 27 orders
phy.phyla <- tax_glom(phy.gvhd, taxrank = 'Phylum')
ntaxa(phy.phyla)
# 6 phyla
remove(phy.genera, phy.phyla)

######################################################################################################################################
# PATIENT CHARACTERISTICS

patients <- metadata_gvhd[which(metadata_gvhd$time_group=="pre"),]
cefepime <- patients[which(patients$treatment=="cefepime"),]
summary(cefepime$age)
table(cefepime$sex)
prop.table(table(cefepime$sex))
table(cefepime$race)
prop.table(table(cefepime$race))
table(cefepime$diagnosis)
prop.table(table(cefepime$diagnosis))
table(cefepime$donor_hsct)
prop.table(table(cefepime$donor_hsct))
table(cefepime$source_hsct)
prop.table(table(cefepime$source_hsct))
table(cefepime$hla)
prop.table(table(cefepime$hla))
table(cefepime$gvhd_prophy, cefepime$hsct_prep)
prop.table(table(cefepime$gvhd_prophy, cefepime$hsct_prep))

anaerobic <- patients[which(patients$treatment=="anaerobic"),]
summary(anaerobic$age)
table(anaerobic$sex)
prop.table(table(anaerobic$sex))
table(anaerobic$race)
prop.table(table(anaerobic$race))
table(anaerobic$diagnosis)
prop.table(table(anaerobic$diagnosis))
table(anaerobic$donor_hsct)
prop.table(table(anaerobic$donor_hsct))
table(anaerobic$source_hsct)
prop.table(table(anaerobic$source_hsct))
table(anaerobic$hla)
prop.table(table(anaerobic$hla))
table(anaerobic$gvhd_prophy, anaerobic$hsct_prep)
prop.table(table(anaerobic$gvhd_prophy, anaerobic$hsct_prep))

######################################################################################################################################
# ANALYSES OF ALPHA DIVERSITY

# Sample Reads
sum(sample_sums(phy.gvhd))
mean(sample_sums(phy.gvhd))

# Look at diversity of samples by antibiotic category
metadata_gvhd$treatment[metadata_gvhd$treatment=="cefepime"] <- "Cefepime"
metadata_gvhd$treatment[metadata_gvhd$treatment=="anaerobic"] <- "Anaerobic Antibiotics"
metadata_gvhd$treatment <- factor(metadata_gvhd$treatment, levels=c("Cefepime", "Anaerobic Antibiotics"))
metadata_gvhd$time_group <- factor(metadata_gvhd$time_group, levels=c("pre", "post"))

sdi_abx <- ggplot(metadata_gvhd, aes(x=time_group, y=Shannon)) + 
  geom_boxplot(aes(group=time_group), color="black", fill="grey50") +
  theme(axis.title.y = element_text(size=13), axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=13), axis.text.x = element_text(size=10),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text=element_text(size=14), legend.position="none") + ylab("Shannon Diversity Index") + xlab("Sample Timing") +
  scale_x_discrete(labels=c("Before Antibiotics", "During Antibiotics")) +
  facet_wrap(~treatment) 

chao_abx <- ggplot(metadata_gvhd, aes(x=time_group, y=Chao1)) + 
  geom_boxplot(aes(group=time_group), color="black", fill="grey50") +
  theme(axis.title.y = element_text(size=13), axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=13), axis.text.x = element_text(size=10),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text=element_text(size=14), legend.position="none") + ylab("Chao1 Index") + xlab("Sample Timing") + ylim(0, 60) +
  scale_x_discrete(labels=c("Before Antibiotics", "During Antibiotics")) +
  facet_wrap(~treatment) 

png(file="Results/Alpha Diversity by Antibiotics.png", 
    width = 12, height = 6, units = 'in', res = 300)
plot_grid(sdi_abx, chao_abx, labels="auto", label_size=16)
dev.off()
remove(sdi_abx, chao_abx)

# Did alpha diversity decline with cefepime treatment? 

cefepime <- subset(metadata_gvhd, treatment=="Cefepime")
tapply(cefepime$Shannon, cefepime$time_group, summary)
tapply(cefepime$Chao1, cefepime$time_group, summary)

cef_pre_sdi <- subset(cefepime, time_group=="pre", Shannon, drop = TRUE)
cef_post_sdi <- subset(cefepime, time_group=="post", Shannon, drop = TRUE)
cef_pre_chao <- subset(cefepime, time_group=="pre", Chao1, drop = TRUE)
cef_post_chao <- subset(cefepime, time_group=="post", Chao1, drop = TRUE)
wilcox.test(cef_pre_sdi, cef_post_sdi, paired = TRUE, alternative = "two.sided")
wilcox.test(cef_pre_chao, cef_post_chao, paired = TRUE, alternative = "two.sided")
remove(cef_pre_sdi, cef_pre_chao, cef_post_sdi, cef_post_chao)
# No change in SDI or Chao1 with cefepime treatment

# Did alpha diversity decline with anaerobic treatment? 

anaerobic <- subset(metadata_gvhd, treatment=="Anaerobic Antibiotics")
tapply(anaerobic$Shannon, anaerobic$time_group, summary)
tapply(anaerobic$Chao1, anaerobic$time_group, summary)

ana_pre_sdi <- subset(anaerobic, time_group=="pre", Shannon, drop = TRUE)
ana_post_sdi <- subset(anaerobic, time_group=="post", Shannon, drop = TRUE)
ana_pre_chao <- subset(anaerobic, time_group=="pre", Chao1, drop = TRUE)
ana_post_chao <- subset(anaerobic, time_group=="post", Chao1, drop = TRUE)
wilcox.test(ana_pre_sdi, ana_post_sdi, paired = TRUE, alternative = "two.sided")
wilcox.test(ana_pre_chao, ana_post_chao, paired = TRUE, alternative = "two.sided")
remove(ana_pre_sdi, ana_pre_chao, ana_post_sdi, ana_post_chao)
# Loss of diversity/richness as measured by SDI and Chao1 in anaerobic group

tapply(cefepime$abx_day, cefepime$time_group, summary)
tapply(anaerobic$abx_day, anaerobic$time_group, summary)
cefepime_days <- subset(cefepime, antibiotic_group=="cefepime_post", abx_day, drop = TRUE)
anaerobic_days <- subset(anaerobic, antibiotic_group=="anaerobic_post", abx_day, drop = TRUE)
wilcox.test(cefepime_days, anaerobic_days, alternative = "two.sided")
remove(cefepime_days, anaerobic_days)

######################################################################################################################################
# COMPOSITION PLOTS

# Transform to relative abundances
phy.gvhd <- transform_sample_counts(phy.gvhd, function(Abundance) Abundance/sum(Abundance))
sample_sums(phy.gvhd)  
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling
gvhd_df <- psmelt(phy.gvhd)

# Create dataframes with overall relative abundances of phyla and genera
gvhd_df$Phylum <- as.character(gvhd_df$Phylum)
phyla_abundances <- aggregate(gvhd_df$Abundance, by=list(Phylum=gvhd_df$Phylum), FUN=sum)
phyla_abundances$x <- (phyla_abundances$x)/(nsamples(phy))
phyla_abundances <- rename(phyla_abundances, c("x"="phyla_Ab"))
sum(phyla_abundances$phyla_Ab)    # Should sum to 1 (sum of relative abundances of phyla)
nrow(phyla_abundances)            # Corresponds to # of unique phyla
gvhd_df$Genus <- as.character(gvhd_df$Genus)
otu_abundances <- aggregate(gvhd_df$Abundance, by=list(Phylum=gvhd_df$Phylum, Genus=gvhd_df$Genus,
                                                       OTU=gvhd_df$OTU), FUN=mean)
genus_abundances <- aggregate(otu_abundances$x, by=list(Phylum=otu_abundances$Phylum,
                                                        Genus=otu_abundances$Genus), FUN=sum)
genus_abundances <- rename(genus_abundances, c("x"="genus_Ab"))
sum(genus_abundances$genus_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(genus_abundances)            # Corresponds to # of unique genera
abundances <- merge(genus_abundances, phyla_abundances, by="Phylum")

genus_abundances <- arrange(genus_abundances, genus_Ab, decreasing=TRUE)  
TOPGenera <- unique(genus_abundances$Genus[1:17])
genus_df <- genus_abundances[genus_abundances$Genus %in% TOPGenera,]
genus_df$Genus <- factor(genus_df$Genus, levels = genus_df$Genus[order(-genus_df$genus_Ab)])

otu_abundances <- rename(otu_abundances, c("x"="otu_Ab"))
sum(otu_abundances$otu_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(otu_abundances)            # Corresponds to # of unique genera

otu_abundances <- arrange(otu_abundances, otu_Ab, decreasing=TRUE)  
TOPOTU <- unique(otu_abundances$OTU[1:30])

m1 <- subset(gvhd_df, treatment=="cefepime" & time_group=="pre")
m1$gvhd_abundances <- (m1$Abundance)/(length(unique(m1$study_id))) 
m2 <- subset(gvhd_df, treatment=="cefepime" & time_group=="post")
m2$gvhd_abundances <- (m2$Abundance)/(length(unique(m2$study_id))) 
m3 <- subset(gvhd_df, treatment=="anaerobic" & time_group=="pre")
m3$gvhd_abundances <- (m3$Abundance)/(length(unique(m3$study_id))) 
m4 <- subset(gvhd_df, treatment=="anaerobic" & time_group=="post")
m4$gvhd_abundances <- (m4$Abundance)/(length(unique(m4$study_id))) 

m2b <- merge(m1, m2, all.x=TRUE, all.y=TRUE)
m3b <- merge(m2b, m3, all.x=TRUE, all.y=TRUE)
gvhd_final <- merge(m3b, m4, all.x=TRUE, all.y=TRUE)
remove(genus_abundances, genus_df, otu_abundances, phyla_abundances, abundances, m1, m2, m3, m4, m2b, m3b)

palette1 <- c("indianred4", "peru", "navajowhite3", "midnightblue", "tomato2", "mediumturquoise", "gray1", "thistle", 
              "dodgerblue3", "coral3", "goldenrod4", "mistyrose4", "mediumorchid4", 
              "darkslateblue", "darkolivegreen4", "black", "goldenrod2",  "wheat4", "darkblue", "firebrick", "gray70", "maroon", 
              "steelblue4", "skyblue2", "khaki2", "green4", "#6A3D9A")

gvhd_final$time_group <- factor(gvhd_final$time_group, levels=c("pre", "post"))
gvhd_final$treatment[gvhd_final$treatment=="cefepime"] <- "Cefepime"
gvhd_final$treatment[gvhd_final$treatment=="anaerobic"] <- "Anaerobic Antibiotics"
gvhd_final$treatment <- factor(gvhd_final$treatment, levels=c("Cefepime", "Anaerobic Antibiotics"))
gvhd_final$Order <- as.character(gvhd_final$Order)
gvhd_final$Order[gvhd_final$Order=="Firmicutes_unclassified"] <- "Unclassified Firmicutes"
gvhd_final$Order[gvhd_final$Order=="Proteobacteria_unclassified"] <- "Unclassified Proteobacteria"

# First look at the effect of anaerobic antibiotics vs. cefepime on the gut microbiota

plot_composition <- ggplot(arrange(gvhd_final, Order), aes(x=time_group, y=gvhd_abundances, fill=Order)) +
  geom_bar(stat="identity", position="stack") +
  theme(legend.text=element_text(size=16), legend.title=element_text(size=18), 
        axis.title.y = element_text(angle=90, size=18), 
        axis.text.y = element_text(size=18, colour="black"), axis.title.x = element_text(size=18), 
        axis.text.x = element_text(size=18, colour="black"), plot.title = element_text(size=18, hjust=1.4), 
        strip.text.x = element_text(size = 20), strip.background=element_rect(fill="white")) + 
  scale_fill_manual(values=palette1) + 
  guides(fill = guide_legend(ncol=1, byrow=TRUE, title="Orders")) + 
  scale_x_discrete(labels=c("Before Antibiotics", "During Antibiotics")) +
  ylab("Relative Abundance") + xlab("Sample Timing") +
  facet_wrap(~treatment)

png(file="Results/Composition Bar Plot by Antibiotics.png", 
    width = 16, height = 11, units = 'in', res = 600)
print(plot_composition)
dev.off()

######################################################################################################################################
# CEFEPIME TREATMENT

library(vegan)
phy.cefepime <- subset_samples(phy.gvhd, treatment=="cefepime")
adonis(distance(phy.cefepime, method="bray", strata=study_id) ~ time_group,
       data = cefepime)
# No difference in overall microbiome composition before and after cefepime treatment

# Compare composition at the ORDER level
phy.cefepime.order <- tax_glom(phy.cefepime, taxrank = 'Order')
ntaxa(phy.cefepime.order)
phy.cefepime.order <- transform_sample_counts(phy.cefepime.order, function(Abundance) Abundance/sum(Abundance))
sample_sums(phy.cefepime.order)  
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling
cefepime_order <- psmelt(phy.cefepime.order)
cefepime_order <- cefepime_order[,-1]
remove(phy.cefepime.order)
table(cefepime_order$Order)
# Confirm that there are 27 orders

# Limit analyses to orders present in at least 25% of samples
table(cefepime_order$Order)
cefepime_select <- subset(cefepime_order, Abundance>0)
table(cefepime_select$Order)
# Focus only on orders present in at least 13 samples

# Test for significant change in abundance of ORDERS with cefepime treatment

cefepime_Acidaminococcales <- subset(cefepime_order, Order=="Acidaminococcales")
cefepime_Acidaminococcales <- cefepime_Acidaminococcales[order(cefepime_Acidaminococcales$study_id),]
cefepime_Acidaminococcales <- cefepime_Acidaminococcales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_Acidaminococcales <- subset(cefepime_Acidaminococcales, time_group=="pre")
wilcox.test(cefepime_Acidaminococcales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_Acidaminococcales)

# Actinomycetales: P=0.40
# Bacteroidales: P=0.90
# Bifidobacteriales: P=0.18
# Burkholderiales: P=0.06
# Clostridiales: P=0.10
# Desulfovibrionales: P=0.09
# Eggerthellales: P=0.10
# Enterobacterales: P=0.01* 
# Erysipelotrichales: P=0.52
# Lactobacillales: P=0.04*
# Veillonellales: P=0.70

# THE FOLLOWING ORDERS CHANGED IN ABUNDANCE IN PAIRED SAMPLES WITH CEFEPIME: 
# Enterobacteriales (declined), Lactobacillales (declined)

######################################################################################################################################
# ANAEROBIC TREATMENT

phy.anaerobic <- subset_samples(phy.gvhd, treatment=="anaerobic")
adonis(distance(phy.anaerobic, method="bray", strata=study_id) ~ time_group,
       data = anaerobic)
# Markedly different overall microbiome composition after treatment with anaerobic antibiotics

# Compare composition at the Order level 
phy.anaerobic.order <- tax_glom(phy.anaerobic, taxrank = 'Order')
ntaxa(phy.anaerobic.order)
phy.anaerobic.order <- transform_sample_counts(phy.anaerobic.order, function(Abundance) Abundance/sum(Abundance))
sample_sums(phy.anaerobic.order)  
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling
anaerobic_order <- psmelt(phy.anaerobic.order)
anaerobic_order <- anaerobic_order[,-1]

# Limit analyses to orders present in at least 25% of samples
table(anaerobic_order$Order)
anaerobic_select <- subset(anaerobic_order, Abundance>0)
table(anaerobic_select$Order)
# Focus only on orders present in at least 5 samples

# Test for significant change in abundance of ORDERS with anaerobic antibiotics
anaerobic_Veillonellales <- subset(anaerobic_order, Order=="Veillonellales")
anaerobic_Veillonellales <- anaerobic_Veillonellales[order(anaerobic_Veillonellales$study_id),]
anaerobic_Veillonellales <- anaerobic_Veillonellales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_Veillonellales <- subset(anaerobic_Veillonellales, time_group=="pre")
wilcox.test(anaerobic_Veillonellales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_Veillonellales)

# Actinomycetales: P=0.83
# Bacteroidales: P=0.08
# Bifidobacteriales: P=0.01*
# Burkholderiales: P=0.67
# Clostridiales: P=0.002*
# Eggerthellales: P=0.59
# Enterobacterales: P>0.99
# Erysipelotrichales: P=0.32
# Lactobacillales: P=0.01*
# Veillonellales: P=0.94

# THE FOLLOWING ORDERS CHANGED IN ABUNDANCE IN PAIRED SAMPLES WITH ANAEROBIC ANTIBIOTICS: 
# Bifidobacteriales (decreased), Clostridiales (decreased), Lactobacillales (increased)

######################################################################################################################################
# TEST GENUS-LEVEL DIFFERENCES IN CLOSTRIDIALES AND BIFIDOBACTERIALES

# Compare composition at the GENUS level for Orders that differed after antibiotic exposure
phy.anaerobic.genus <- tax_glom(phy.anaerobic, taxrank = 'Genus')

phy.anaerobic.genus <- transform_sample_counts(phy.anaerobic.genus, function(Abundance) Abundance/sum(Abundance))
sample_sums(phy.anaerobic.genus)  
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling
anaerobic_genus <- psmelt(phy.anaerobic.genus)
anaerobic_genus <- anaerobic_genus[,-1]

# CLOSTRIDIALES
anaerobic_clostridiales_genus <- subset(anaerobic_genus, Order=="Clostridiales")
tapply(anaerobic_clostridiales_genus$Abundance, anaerobic_clostridiales_genus$Genus, FUN=mean)
# Focus on genera with mean relative abundance >0.005 (0.5%)

genus_test <- subset(anaerobic_clostridiales_genus, Genus=="Blautia")
genus_test <- genus_test[order(genus_test$study_id),]
tapply(genus_test$Abundance, genus_test$time_group, summary)
genus_test <- genus_test %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
genus_test <- subset(genus_test, time_group=="pre")
wilcox.test(genus_test$D_Abundance, mu = 0, alternative = "two.sided")
remove(genus_test)
# Blautia: P=0.01

genus_test <- subset(anaerobic_clostridiales_genus, Genus=="Clostridium")
genus_test <- genus_test[order(genus_test$study_id),]
tapply(genus_test$Abundance, genus_test$time_group, summary)
genus_test <- genus_test %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
genus_test <- subset(genus_test, time_group=="pre")
wilcox.test(genus_test$D_Abundance, mu = 0, alternative = "two.sided")
remove(genus_test)
# Clostridium: P=0.02

genus_test <- subset(anaerobic_clostridiales_genus, Genus=="Faecalibacterium")
genus_test <- genus_test[order(genus_test$study_id),]
tapply(genus_test$Abundance, genus_test$time_group, summary)
genus_test <- genus_test %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
genus_test <- subset(genus_test, time_group=="pre")
wilcox.test(genus_test$D_Abundance, mu = 0, alternative = "two.sided")
remove(genus_test)
# Faecalibacterium: P=0.18

genus_test <- subset(anaerobic_clostridiales_genus, Genus=="Lachnoclostridium")
genus_test <- genus_test[order(genus_test$study_id),]
tapply(genus_test$Abundance, genus_test$time_group, summary)
genus_test <- genus_test %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
genus_test <- subset(genus_test, time_group=="pre")
wilcox.test(genus_test$D_Abundance, mu = 0, alternative = "two.sided")
remove(genus_test)
# Lachnoclostridium: P=0.01

genus_test <- subset(anaerobic_clostridiales_genus, Genus=="Roseburia")
genus_test <- genus_test[order(genus_test$study_id),]
tapply(genus_test$Abundance, genus_test$time_group, summary)
genus_test <- genus_test %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
genus_test <- subset(genus_test, time_group=="pre")
wilcox.test(genus_test$D_Abundance, mu = 0, alternative = "two.sided")
remove(genus_test)
# Roseburia: P=0.06

genus_test <- subset(anaerobic_clostridiales_genus, Genus=="Ruthenibacterium")
genus_test <- genus_test[order(genus_test$study_id),]
tapply(genus_test$Abundance, genus_test$time_group, summary)
genus_test <- genus_test %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
genus_test <- subset(genus_test, time_group=="pre")
wilcox.test(genus_test$D_Abundance, mu = 0, alternative = "two.sided")
remove(genus_test)
# Ruthenibacterium: P=0.06

# BIFIDOBACTERIALES
anaerobic_bifidobacteriales_genus <- subset(anaerobic_genus, Order=="Bifidobacteriales")
tapply(anaerobic_bifidobacteriales_genus$Abundance, anaerobic_bifidobacteriales_genus$Genus, FUN=sum)
# Bifidobacterium

anaerobic_bifido <- subset(anaerobic_bifidobacteriales_genus, Genus=="Bifidobacterium")
anaerobic_bifido <- anaerobic_bifido[order(anaerobic_bifido$study_id),]
tapply(anaerobic_bifido$Abundance, anaerobic_bifido$time_group, summary)
anaerobic_bifido <- anaerobic_bifido %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_bifido <- subset(anaerobic_bifido, time_group=="pre")
wilcox.test(anaerobic_bifido$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_bifido, anaerobic_bifidobacteriales_genus)
# Bifidobacterium: P=0.02

######################################################################################################################################
# CHANGE IN MICROBIOME COMPOSITION IN CHILDREN WHO SUBSEQUENTLY DEVELOPED ACUTE GVHD OF THE GUT/LIVER

library(vegan)
phy.gvhd_y <- subset_samples(phy.gvhd, agvhd_gutliver=="1")
gvhd <- subset(metadata_gvhd, agvhd_gutliver=="1")
adonis(distance(phy.gvhd_y, method="bray", strata=study_id) ~ time_group,
       data = gvhd)
nsamples(phy.gvhd_y)
remove(gvhd)

# Compare composition at the Order level 
phy.gvhd.order <- tax_glom(phy.gvhd_y, taxrank = 'Order')
ntaxa(phy.gvhd.order)
phy.gvhd.order <- transform_sample_counts(phy.gvhd.order, function(Abundance) Abundance/sum(Abundance))
sample_sums(phy.gvhd.order)  
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling
gvhd_order <- psmelt(phy.gvhd.order)
gvhd_order <- gvhd_order[,-1]

gvhd_clostridiales <- subset(gvhd_order, Order=="Clostridiales")
gvhd_clostridiales <- gvhd_clostridiales[order(gvhd_clostridiales$study_id),]
gvhd_clostridiales <- gvhd_clostridiales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
gvhd_clostridiales <- subset(gvhd_clostridiales, time_group=="pre")
wilcox.test(gvhd_clostridiales$D_Abundance, mu = 0, alternative = "two.sided")
remove(gvhd_clostridiales)
# Clostridiales: P=0.06

gvhd_bifidobacteriales <- subset(gvhd_order, Order=="Bifidobacteriales")
gvhd_bifidobacteriales <- gvhd_bifidobacteriales[order(gvhd_bifidobacteriales$study_id),]
gvhd_bifidobacteriales <- gvhd_bifidobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
gvhd_bifidobacteriales <- subset(gvhd_bifidobacteriales, time_group=="pre")
wilcox.test(gvhd_bifidobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(gvhd_bifidobacteriales)
# Bifidobacteriales: P=0.006

remove(phy.gvhd.order, gvhd_order)

######################################################################################################################################
# HEATMAP

library(gtools)

phy.order <-  tax_glom(phy.gvhd, "Order")
taxa_names(phy.order) <-  as.character(tax_table(phy.order)[,c("Order")])

# Focus only on orders that are present in >=25% of samples
all_order <- psmelt(phy.order)
all_order <- all_order[,-1]

# Limit analyses to orders present in at least 25% of samples
table(all_order$Order)
all_select <- subset(all_order, Abundance>0)
table(all_select$Order)
# Focus only on orders present in at least 18 samples

phy.order <- subset_taxa(phy.order, Order!="Acidaminococcales" & Order!="Bacillales" & Order!="Campylobacterales" & 
                           Order!="Coriobacteriales" & Order!="Corynebacteriales" & Order!="Desulfovibrionales" & 
                           Order!="Firmicutes_unclassified" & Order!="Fusobacteriales" & Order!="Micrococcales" & 
                           Order!="Pasteurellales" & Order!="Propionibacteriales" & Order!="Proteobacteria_unclassified" &
                           Order!="Pseudomonadales" & Order!="Selenomonadales" & Order!="Tissierellales" & 
                           Order!="Verrucomicrobiales" & Order!="Xanthomonadales")

c_phy_pre <- data.frame(otu_table(subset_samples(phy.order, treatment == "cefepime" & time_group == "pre")))
names(c_phy_pre) <- paste0("Patient ", seq(1,26))
c_phy_pre <- c_phy_pre[,mixedsort(names(c_phy_pre))]

c_phy_post <-data.frame(otu_table( subset_samples(phy.order, treatment == "cefepime" & time_group == "post")))
names(c_phy_post) <- paste0("Patient ", seq(1,26))
c_phy_post <- c_phy_post[,mixedsort(names(c_phy_post))]

# Percent change
c_change <- 100*(c_phy_post - c_phy_pre)/(c_phy_pre)

# Anaerobic
a_phy_pre <- data.frame(otu_table(subset_samples(phy.order, treatment == "anaerobic" & time_group == "pre")))
names(a_phy_pre) <- paste0("Patient ", seq(27,28))
a_phy_pre <- a_phy_pre[,mixedsort(names(a_phy_pre))]

a_phy_post <-data.frame(otu_table( subset_samples(phy.order, treatment == "anaerobic" & time_group == "post")))
names(a_phy_post) <- paste0("Patient ", seq(27,36))
a_phy_post <- a_phy_post[,mixedsort(names(a_phy_post))]

a_change <- 100*(a_phy_post - a_phy_pre)/(a_phy_pre)

library("ComplexHeatmap")

library("viridis")
library("circlize")
jet.colors <-  col_fn <- colorRamp2(100000*c(-0.001,-0.0005,0,0.005,0.01), c("darkblue", "cornflowerblue", "white", "lightcoral", "red4"))

mat1 <- as.matrix(c_change)
ordered_names <- rownames(mat1)[order(rownames(mat1))]

mat1 <- mat1[match(ordered_names,rownames(mat1)),]
mat2 <- as.matrix(a_change)
mat2 <- mat2[match(ordered_names,rownames(mat2)),]

mat <-  cbind(mat1,mat2)

ht1 <- Heatmap(mat1,name = "% Change in Relative Abundance", 
               cluster_rows = FALSE,
               row_names_max_width = max_text_width(rownames(mat1), gp = gpar(fontsize = 12)),       
               col = jet.colors, cluster_columns = F,
               heatmap_legend_param = list(at = 100000*c(-0.001,-0.0005,0,0.005,0.01)),
               row_names_side = c("left"), column_title = "Cefepime", show_heatmap_legend = F) +
  Heatmap(mat2,name = "post",
          col = jet.colors,
          cluster_columns = F,
          heatmap_legend_param = list(at = 100000*c(-0.001,-0.0005,0,0.005,0.01)),
          show_row_names = F, column_title = "Anaerobic Antibiotics", show_heatmap_legend = F) 

png(file="Results/Heatmap by Antibiotics.png", width = 16, height = 8, units = 'in', res = 600)
draw(ht1)
decorate_heatmap_body("% Change in Relative Abundance", {
  
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  
})
decorate_heatmap_body("post", {
  
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  
})
dev.off()

png(file="Results/Heatmap Legend.png", width = 8, height = 8, units = 'in', res = 600)
col_fun <- colorRamp2(100000*c(-0.001,-0.0005,0,0.005,0.01), c("darkblue", "cornflowerblue", "white", "lightcoral", "red4"))
legend = Legend(col_fun = col_fun, legend_height=unit(4,"cm"), title="% Change in Relative Abundance", title_position="topcenter",
                legend_width=unit(5,"cm"), at = c(-100, 0, 100, 250, 500, 750, 1000), labels_rot = 45, grid_height = unit(0.75,"cm"),
                direction="horizontal")
draw(legend)
dev.off()

remove(a_change, a_phy_post, a_phy_pre, c_change, c_phy_post, c_phy_pre, ht1, legend, mat, mat1, mat2, 
       ordered_names, col_fn, col_fun, jet.colors)

######################################################################################################################################
# TESTING FOR CHANGES IN ABUNDANCE OF BUTYRATE GENES

# Butyrate genes DID NOT decline in abundance with CEFEPIME
cef_butyrate <- subset(metadata_gvhd, treatment=="Cefepime")
tapply(cef_butyrate$butyrate_genes, cef_butyrate$time_group, summary)
cef_pre_butyrate <- subset(cef_butyrate, antibiotic_group=="cefepime_pre", butyrate_genes, drop = TRUE)
cef_post_butyrate <- subset(cef_butyrate, antibiotic_group=="cefepime_post", butyrate_genes, drop = TRUE)
wilcox.test(cef_pre_butyrate, cef_post_butyrate, paired = TRUE, alternative = "two.sided")
remove(cef_butyrate, cef_pre_butyrate, cef_post_butyrate)
# P=0.22

# Butyrate genes declined in abundance with ANAEROBIC ANTIBIOTICS 
ana_butyrate <- subset(metadata_gvhd, treatment=="Anaerobic Antibiotics")
tapply(ana_butyrate$butyrate_genes, ana_butyrate$time_group, summary)
ana_pre_butyrate <- subset(ana_butyrate, antibiotic_group=="anaerobic_pre", butyrate_genes, drop = TRUE)
ana_post_butyrate <- subset(ana_butyrate, antibiotic_group=="anaerobic_post", butyrate_genes, drop = TRUE)
wilcox.test(ana_pre_butyrate, ana_post_butyrate, paired = TRUE, alternative = "two.sided")
remove(ana_butyrate, ana_pre_butyrate, ana_post_butyrate)
# P=0.01

# Butyrate genes were higher in pre- samples from patients who subsequently developed acute GVHD of gut/liver 
pre_butyrate <- subset(metadata_gvhd, time_group=="pre")
pre_gvhd_butyrate <- subset(pre_butyrate, agvhd_gutliver=="1", butyrate_genes, drop = TRUE)
pre_nogvhd_butyrate <- subset(pre_butyrate, agvhd_gutliver=="0", butyrate_genes, drop = TRUE)
summary(pre_gvhd_butyrate)
summary(pre_nogvhd_butyrate)
wilcox.test(pre_gvhd_butyrate, pre_nogvhd_butyrate, alternative = "two.sided")
remove(pre_butyrate, pre_gvhd_butyrate, pre_nogvhd_butyrate)
# P=1.1E-5

# Butyrate genes declined in abundance among patients with acute GVHD of gut/liver 
gvhd_butyrate <- subset(metadata_gvhd, agvhd_gutliver=="1")
gvhd_pre_butyrate <- subset(gvhd_butyrate, antibiotic_group=="cefepime_pre" | antibiotic_group=="anaerobic_pre", butyrate_genes, drop = TRUE)
gvhd_post_butyrate <- subset(gvhd_butyrate, antibiotic_group=="cefepime_post" | antibiotic_group=="anaerobic_post", butyrate_genes, drop = TRUE)
summary(gvhd_pre_butyrate)
summary(gvhd_post_butyrate)
wilcox.test(gvhd_pre_butyrate, gvhd_post_butyrate, paired = TRUE, alternative = "two.sided")
remove(gvhd_butyrate, gvhd_pre_butyrate, gvhd_post_butyrate)
# P=0.0001

# Butyrate genes DID NOT decline in abundance among patients without acute GVHD of gut/liver
no_gvhd_butyrate <- subset(metadata_gvhd, agvhd_gutliver=="0")
no_gvhd_pre_butyrate <- subset(no_gvhd_butyrate, antibiotic_group=="cefepime_pre" | antibiotic_group=="anaerobic_pre", butyrate_genes, drop = TRUE)
no_gvhd_post_butyrate <- subset(no_gvhd_butyrate, antibiotic_group=="cefepime_post" | antibiotic_group=="anaerobic_post", butyrate_genes, drop = TRUE)
summary(no_gvhd_pre_butyrate)
summary(no_gvhd_post_butyrate)
wilcox.test(no_gvhd_pre_butyrate, no_gvhd_post_butyrate, paired = TRUE, alternative = "two.sided")
remove(no_gvhd_butyrate, no_gvhd_pre_butyrate, no_gvhd_post_butyrate)
# P=0.30

metadata_gvhd$treatment <- factor(metadata_gvhd$treatment, levels=c("Cefepime", "Anaerobic Antibiotics"))
levels(metadata_gvhd$agvhd_gutliver) <- c(levels(metadata_gvhd$agvhd_gutliver), "No")
levels(metadata_gvhd$agvhd_gutliver) <- c(levels(metadata_gvhd$agvhd_gutliver), "Yes")
metadata_gvhd$agvhd_gutliver[metadata_gvhd$agvhd_gutliver == "0"] <- "No"
metadata_gvhd$agvhd_gutliver[metadata_gvhd$agvhd_gutliver == "1"] <- "Yes"

butyrate <- ggplot(metadata_gvhd, aes(x=time_group, y=butyrate_genes)) + 
  geom_point(aes(group=time_group)) + geom_line(aes(group=study_id, linetype = agvhd_gutliver)) +
  theme(axis.title.y = element_text(size=10), axis.text.y = element_text(size=9), legend.text=element_text(size=8), legend.title=element_text(size=9), 
        axis.title.x = element_text(size=10), axis.text.x = element_text(size=9),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text.y=element_text(size=10), legend.position="right") + ylab("Butyrate Gene Relative Abundance (RPKM)") + xlab("Sample Timing") +
  scale_x_discrete(labels=c("Before Antibiotics", "During Antibiotics")) + scale_color_manual(values = c("black", "red")) + ylim(0,1150) +
  guides(shape=guide_legend("Acute Gut/Liver GVHD"), linetype=guide_legend("Acute Gut/Liver GVHD")) +
  facet_wrap(~treatment) 

png(file="Results/Butyrate Genes by Antibiotics.png", width = 8, height = 3.5, units = 'in', res = 600)
plot(butyrate)
dev.off()
remove(butyrate)

# COMPARE CLOSTRIDIALES RELATIVE ABUNDANCE BY ANTIBIOTICS AND GVHD STATUS
phy.gvhd.order <- tax_glom(phy.gvhd, taxrank = 'Order')
gvhd_order <- psmelt(phy.gvhd.order)
clostridiales_df <- subset(gvhd_order, Order=="Clostridiales")
clostridiales_df <- clostridiales_df[,-1]

clostridiales_df$time_group <- factor(clostridiales_df$time_group, levels=c("pre", "post"))
clostridiales_df$treatment[clostridiales_df$treatment=="cefepime"] <- "Cefepime"
clostridiales_df$treatment[clostridiales_df$treatment=="anaerobic"] <- "Anaerobic Antibiotics"
clostridiales_df$treatment <- factor(clostridiales_df$treatment, levels=c("Cefepime", "Anaerobic Antibiotics"))
clostridiales_df$agvhd_gutliver <- as.factor(clostridiales_df$agvhd_gutliver)
levels(clostridiales_df$agvhd_gutliver) <- c(levels(clostridiales_df$agvhd_gutliver), "No")
levels(clostridiales_df$agvhd_gutliver) <- c(levels(clostridiales_df$agvhd_gutliver), "Yes")
clostridiales_df$agvhd_gutliver[clostridiales_df$agvhd_gutliver == "0"] <- "No"
clostridiales_df$agvhd_gutliver[clostridiales_df$agvhd_gutliver == "1"] <- "Yes"

clostridiales <- ggplot(clostridiales_df, aes(x=time_group, y=Abundance)) + 
  geom_point(aes(group=time_group)) + geom_line(aes(group=study_id, linetype = agvhd_gutliver)) +
  theme(axis.title.y = element_text(size=10), axis.text.y = element_text(size=9), legend.text=element_text(size=8), legend.title=element_text(size=9), 
        axis.title.x = element_text(size=10), axis.text.x = element_text(size=9),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text.y=element_text(size=10), legend.position="right") + ylab("Clostridiales Relative Abundance") + xlab("Sample Timing") +
  scale_x_discrete(labels=c("Before Antibiotics", "During Antibiotics")) + scale_color_manual(values = c("black", "red")) + ylim(0,1) +
  guides(shape=guide_legend("Acute Gut/Liver GVHD"), linetype=guide_legend("Acute Gut/Liver GVHD")) +
  facet_wrap(~treatment) 

png(file="Results/Clostridiales Abundance by Antibiotics.png", width = 8, height = 3.5, units = 'in', res = 600)
plot(clostridiales)
dev.off()
remove(clostridiales, clostridiales_df)