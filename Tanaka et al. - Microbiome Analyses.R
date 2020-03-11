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
setwd("C:/Users/msk37/Google Drive/Research/SAS_BMT/UMass_GVHD") 

phy.gvhd <- readRDS("phy.gvhd.rds")
metadata_gvhd <- data.frame(sample_data(phy.gvhd))

phy.gvhd <- transform_sample_counts(phy.gvhd, function(OTU) round(OTU))
sample_sums(phy.gvhd)
nsamples(phy.gvhd)
ntaxa(phy.gvhd)
# 229 taxa
phy.genera <- tax_glom(phy.gvhd, taxrank = 'Genus')
ntaxa(phy.genera)
# 79 genera
phy.orders <- tax_glom(phy.gvhd, taxrank = 'Order')
ntaxa(phy.orders)
# 18 orders
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
table(cefepime$hsct_source)
prop.table(table(cefepime$hsct_source))
table(cefepime$hla)
prop.table(table(cefepime$hla))

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
table(anaerobic$hsct_source)
prop.table(table(anaerobic$hsct_source))
table(anaerobic$hla)
prop.table(table(anaerobic$hla))

######################################################################################################################################
# ANALYSES OF ALPHA DIVERSITY

# Sample Reads
sum(sample_sums(phy.gvhd))
# total number of sequencing reads = 338917349
mean(sample_sums(phy.gvhd))
# mean number of sequencing reads = 6052096

# Look at diversity of samples by antibiotic category
metadata_gvhd$treatment[metadata_gvhd$treatment=="cefepime"] <- "Cefepime"
metadata_gvhd$treatment[metadata_gvhd$treatment=="anaerobic"] <- "Anaerobic Antibiotics"
metadata_gvhd$treatment <- factor(metadata_gvhd$treatment, levels=c("Cefepime", "Anaerobic Antibiotics"))
metadata_gvhd$time_group <- factor(metadata_gvhd$time_group, levels=c("pre", "post"))

sdi_abx <- ggplot(metadata_gvhd, aes(x=time_group, y=Shannon)) + 
  geom_boxplot(aes(group=time_group), color="black", fill="grey50") +
  theme(axis.title.y = element_text(size=13), axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=13), axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text=element_text(size=14), legend.position="none") + ylab("Shannon Diversity Index") + xlab("Sample Timing") +
        scale_x_discrete(labels=c("Before", "During")) +
        facet_wrap(~treatment) 

chao_abx <- ggplot(metadata_gvhd, aes(x=time_group, y=Chao1)) + 
  geom_boxplot(aes(group=time_group), color="black", fill="grey50") +
  theme(axis.title.y = element_text(size=13), axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=13), axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text=element_text(size=14), legend.position="none") + ylab("Chao1 Index") + xlab("Sample Timing") + ylim(0, 60) +
        scale_x_discrete(labels=c("Before", "During")) +
        facet_wrap(~treatment) 

png(file="Results/Alpha Diversity by Antibiotics.png", 
    width = 10, height = 5, units = 'in', res = 300)
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
# Time on antibiotics at time of sample obtained while on antibiotics tends be longer in patients on anaerobic antibiotics

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

palette1 <- c("indianred4", "tomato2", "midnightblue", "navajowhite3", "peru", "mediumturquoise", "thistle", "gray1",
                      "dodgerblue3", "coral3", "goldenrod4", "mistyrose4", "mediumorchid4", 
                     "darkslateblue", "darkolivegreen4", "goldenrod2", "black",  "wheat4")

gvhd_final$time_group <- factor(gvhd_final$time_group, levels=c("pre", "post"))
gvhd_final$treatment[gvhd_final$treatment=="cefepime"] <- "Cefepime"
gvhd_final$treatment[gvhd_final$treatment=="anaerobic"] <- "Anaerobic Antibiotics"
gvhd_final$treatment <- factor(gvhd_final$treatment, levels=c("Cefepime", "Anaerobic Antibiotics"))

# First look at the effect of anaerobic antibiotics vs. cefepime on the gut microbiota

plot_composition <- ggplot(arrange(gvhd_final, Order), aes(x=time_group, y=gvhd_abundances, fill=Order)) +
  geom_bar(stat="identity", position="stack") +
  theme(legend.text=element_text(size=16), legend.title=element_text(size=18), 
        axis.title.y = element_text(angle=90, size=18), 
        axis.text.y = element_text(size=16, colour="black"), axis.title.x = element_text(size=18), 
        axis.text.x = element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust=1.4), 
        strip.text.x = element_text(size = 16), strip.background=element_rect(fill="white")) + 
  scale_fill_manual(values=palette1) + 
  guides(fill = guide_legend(ncol=1, byrow=TRUE, title="Orders")) + 
  scale_x_discrete(labels=c("Before", "During")) +
  ylab("Relative Abundance") + xlab("Sample Timing") +
  facet_wrap(~treatment)

png(file="Results/Composition Bar Plot by Antibiotics.png", 
    width = 13, height = 8, units = 'in', res = 600)
print(plot_composition)
dev.off()

# Now look at changes in composition by presence of gut/liver GVHD

gvhd_df <- subset(gvhd_df, !is.na(agvhd_gutliver))
m1 <- subset(gvhd_df, agvhd_gutliver=="0" & time_group=="pre")
m1$gvhd_abundances <- (m1$Abundance)/(length(unique(m1$study_id))) 
m2 <- subset(gvhd_df, agvhd_gutliver=="0" & time_group=="post")
m2$gvhd_abundances <- (m2$Abundance)/(length(unique(m2$study_id))) 
m3 <- subset(gvhd_df, agvhd_gutliver=="1" & time_group=="pre")
m3$gvhd_abundances <- (m3$Abundance)/(length(unique(m3$study_id))) 
m4 <- subset(gvhd_df, agvhd_gutliver=="1" & time_group=="post")
m4$gvhd_abundances <- (m4$Abundance)/(length(unique(m4$study_id))) 

m2b <- merge(m1, m2, all.x=TRUE, all.y=TRUE)
m3b <- merge(m2b, m3, all.x=TRUE, all.y=TRUE)
gvhd_final <- merge(m3b, m4, all.x=TRUE, all.y=TRUE)
remove(m1, m2, m3, m4, m2b, m3b)

gvhd_final$time_group <- factor(gvhd_final$time_group, levels=c("pre", "post"))
gvhd_final$agvhd_gutliver[gvhd_final$agvhd_gutliver=="1"] <- "Acute Gut/Liver GVHD"
gvhd_final$agvhd_gutliver[gvhd_final$agvhd_gutliver=="0"] <- "No Acute Gut/Liver GVHD"
gvhd_final$agvhd_gutliver <- factor(gvhd_final$agvhd_gutliver, levels=c("No Acute Gut/Liver GVHD", "Acute Gut/Liver GVHD"))
gvhd_final$treatment <- factor(gvhd_final$treatment, levels=c("Cefepime", "Anaerobic Antibiotics"))

plot_composition <- ggplot(arrange(gvhd_final, Order), aes(x=time_group, y=gvhd_abundances, fill=Order)) +
   geom_bar(stat="identity", position="stack") +
   theme(legend.text=element_text(size=16), legend.title=element_text(size=18), 
         axis.title.y = element_text(angle=90, size=18), 
         axis.text.y = element_text(size=16, colour="black"), axis.title.x = element_text(size=18), 
         axis.text.x = element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust=1.4), 
         strip.text.x = element_text(size = 16), strip.background=element_rect(fill="white")) + 
   scale_fill_manual(values=palette1) + 
   guides(fill = guide_legend(ncol=1, byrow=TRUE, title="Orders")) + 
   scale_x_discrete(labels=c("Before", "During")) +
   ylab("Relative Abundance") + xlab("Sample Timing") +
   facet_wrap(~agvhd_gutliver)
 
png(file="Results/Composition Bar Plot by GVHD Status.png", 
     width = 13, height = 8, units = 'in', res = 600)
print(plot_composition)
dev.off()
remove(plot_composition, gvhd_final, gvhd_df)

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
# Confirm that there are 18 orders

# Test for significant change in abundance of ORDERS with cefepime treatment

cefepime_actinomycetales <- subset(cefepime_order, Order=="Actinomycetales")
cefepime_actinomycetales <- cefepime_actinomycetales[order(cefepime_actinomycetales$study_id),]
cefepime_actinomycetales <- cefepime_actinomycetales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_actinomycetales <- subset(cefepime_actinomycetales, time_group=="pre")
wilcox.test(cefepime_actinomycetales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_actinomycetales)
# P=0.36

cefepime_bacillales <- subset(cefepime_order, Order=="Bacillales")
cefepime_bacillales <- cefepime_bacillales[order(cefepime_bacillales$study_id),]
cefepime_bacillales <- cefepime_bacillales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_bacillales <- subset(cefepime_bacillales, time_group=="pre")
wilcox.test(cefepime_bacillales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_bacillales)
# P>0.99

cefepime_bacteroidales <- subset(cefepime_order, Order=="Bacteroidales")
cefepime_bacteroidales <- cefepime_bacteroidales[order(cefepime_bacteroidales$study_id),]
cefepime_bacteroidales <- cefepime_bacteroidales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_bacteroidales <- subset(cefepime_bacteroidales, time_group=="pre")
wilcox.test(cefepime_bacteroidales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_bacteroidales)
# P=0.92

cefepime_bifidobacteriales <- subset(cefepime_order, Order=="Bifidobacteriales")
cefepime_bifidobacteriales <- cefepime_bifidobacteriales[order(cefepime_bifidobacteriales$study_id),]
cefepime_bifidobacteriales <- cefepime_bifidobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_bifidobacteriales <- subset(cefepime_bifidobacteriales, time_group=="pre")
wilcox.test(cefepime_bifidobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_bifidobacteriales)
# P=0.54

cefepime_burkholderiales <- subset(cefepime_order, Order=="Burkholderiales")
cefepime_burkholderiales <- cefepime_burkholderiales[order(cefepime_burkholderiales$study_id),]
cefepime_burkholderiales <- cefepime_burkholderiales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_burkholderiales <- subset(cefepime_burkholderiales, time_group=="pre")
wilcox.test(cefepime_burkholderiales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_burkholderiales)
# P=0.48

cefepime_campylobacterales <- subset(cefepime_order, Order=="Campylobacterales")
remove(cefepime_campylobacterales)
# Only one sample with reads assigned to Campylobacterales in cefipime group

cefepime_clostridiales <- subset(cefepime_order, Order=="Clostridiales")
cefepime_clostridiales <- cefepime_clostridiales[order(cefepime_clostridiales$study_id),]
cefepime_clostridiales <- cefepime_clostridiales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_clostridiales <- subset(cefepime_clostridiales, time_group=="pre")
wilcox.test(cefepime_clostridiales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_clostridiales)
# P=0.06

cefepime_coriobacteriales <- subset(cefepime_order, Order=="Coriobacteriales")
cefepime_coriobacteriales <- cefepime_coriobacteriales[order(cefepime_coriobacteriales$study_id),]
cefepime_coriobacteriales <- cefepime_coriobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_coriobacteriales <- subset(cefepime_coriobacteriales, time_group=="pre")
wilcox.test(cefepime_coriobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_coriobacteriales)
# P=0.11

cefepime_desulfovibrionales <- subset(cefepime_order, Order=="Desulfovibrionales")
cefepime_desulfovibrionales <- cefepime_desulfovibrionales[order(cefepime_desulfovibrionales$study_id),]
cefepime_desulfovibrionales <- cefepime_desulfovibrionales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_desulfovibrionales <- subset(cefepime_desulfovibrionales, time_group=="pre")
wilcox.test(cefepime_desulfovibrionales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_desulfovibrionales)
# P=0.36

cefepime_enterobacteriales <- subset(cefepime_order, Order=="Enterobacteriales")
cefepime_enterobacteriales <- cefepime_enterobacteriales[order(cefepime_enterobacteriales$study_id),]
cefepime_enterobacteriales <- cefepime_enterobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_enterobacteriales <- subset(cefepime_enterobacteriales, time_group=="pre")
wilcox.test(cefepime_enterobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_enterobacteriales)
# P=0.17

cefepime_erysipelotrichales <- subset(cefepime_order, Order=="Erysipelotrichales")
cefepime_erysipelotrichales <- cefepime_erysipelotrichales[order(cefepime_erysipelotrichales$study_id),]
cefepime_erysipelotrichales <- cefepime_erysipelotrichales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_erysipelotrichales <- subset(cefepime_erysipelotrichales, time_group=="pre")
wilcox.test(cefepime_erysipelotrichales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_erysipelotrichales)
# P=0.36

cefepime_fusobacteriales <- subset(cefepime_order, Order=="Fusobacteriales")
cefepime_fusobacteriales <- cefepime_fusobacteriales[order(cefepime_fusobacteriales$study_id),]
cefepime_fusobacteriales <- cefepime_fusobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_fusobacteriales <- subset(cefepime_fusobacteriales, time_group=="pre")
wilcox.test(cefepime_fusobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_fusobacteriales)
# P>0.99

cefepime_lactobacillales <- subset(cefepime_order, Order=="Lactobacillales")
cefepime_lactobacillales <- cefepime_lactobacillales[order(cefepime_lactobacillales$study_id),]
cefepime_lactobacillales <- cefepime_lactobacillales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_lactobacillales <- subset(cefepime_lactobacillales, time_group=="pre")
wilcox.test(cefepime_lactobacillales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_lactobacillales)
# P=0.60

cefepime_pasteurellales <- subset(cefepime_order, Order=="Pasteurellales")
cefepime_pasteurellales <- cefepime_pasteurellales[order(cefepime_pasteurellales$study_id),]
cefepime_pasteurellales <- cefepime_pasteurellales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_pasteurellales <- subset(cefepime_pasteurellales, time_group=="pre")
wilcox.test(cefepime_pasteurellales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_pasteurellales)
# P=0.37

cefepime_pseudomonadales <- subset(cefepime_order, Order=="Pseudomonadales")
cefepime_pseudomonadales <- cefepime_pseudomonadales[order(cefepime_pseudomonadales$study_id),]
cefepime_pseudomonadales <- cefepime_pseudomonadales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_pseudomonadales <- subset(cefepime_pseudomonadales, time_group=="pre")
wilcox.test(cefepime_pseudomonadales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_pseudomonadales)
# P>0.99

cefepime_selenomonadales <- subset(cefepime_order, Order=="Selenomonadales")
cefepime_selenomonadales <- cefepime_selenomonadales[order(cefepime_selenomonadales$study_id),]
cefepime_selenomonadales <- cefepime_selenomonadales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_selenomonadales <- subset(cefepime_selenomonadales, time_group=="pre")
wilcox.test(cefepime_selenomonadales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_selenomonadales)
# P=0.75

cefepime_verrucomicrobiales <- subset(cefepime_order, Order=="Verrucomicrobiales")
cefepime_verrucomicrobiales <- cefepime_verrucomicrobiales[order(cefepime_verrucomicrobiales$study_id),]
cefepime_verrucomicrobiales <- cefepime_verrucomicrobiales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_verrucomicrobiales <- subset(cefepime_verrucomicrobiales, time_group=="pre")
wilcox.test(cefepime_verrucomicrobiales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_verrucomicrobiales)
# P=0.09

cefepime_xanthomonadales <- subset(cefepime_order, Order=="Xanthomonadales")
cefepime_xanthomonadales <- cefepime_xanthomonadales[order(cefepime_xanthomonadales$study_id),]
cefepime_xanthomonadales <- cefepime_xanthomonadales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
cefepime_xanthomonadales <- subset(cefepime_xanthomonadales, time_group=="pre")
wilcox.test(cefepime_xanthomonadales$D_Abundance, mu = 0, alternative = "two.sided")
remove(cefepime_xanthomonadales)
# P>0.99

# THERE ARE NO ORDER-LEVEL DIFFERENCES IN ABUNDANCE BETWEEN THE PRE- AND POST-CEFEPIME GROUPS

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

# Test for significant change in abundance of ORDERS with cefepime treatment

anaerobic_actinomycetales <- subset(anaerobic_order, Order=="Actinomycetales")
anaerobic_actinomycetales <- anaerobic_actinomycetales[order(anaerobic_actinomycetales$study_id),]
anaerobic_actinomycetales <- anaerobic_actinomycetales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_actinomycetales <- subset(anaerobic_actinomycetales, time_group=="pre")
wilcox.test(anaerobic_actinomycetales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_actinomycetales)
# P=0.35

anaerobic_bacillales <- subset(anaerobic_order, Order=="Bacillales")
remove(anaerobic_bacillales)
# no reads assigned to Bacillales in anaerobic patients

anaerobic_bacteroidales <- subset(anaerobic_order, Order=="Bacteroidales")
anaerobic_bacteroidales <- anaerobic_bacteroidales[order(anaerobic_bacteroidales$study_id),]
anaerobic_bacteroidales <- anaerobic_bacteroidales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_bacteroidales <- subset(anaerobic_bacteroidales, time_group=="pre")
wilcox.test(anaerobic_bacteroidales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_bacteroidales)
# P=0.04

anaerobic_bifidobacteriales <- subset(anaerobic_order, Order=="Bifidobacteriales")
anaerobic_bifidobacteriales <- anaerobic_bifidobacteriales[order(anaerobic_bifidobacteriales$study_id),]
anaerobic_bifidobacteriales <- anaerobic_bifidobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_bifidobacteriales <- subset(anaerobic_bifidobacteriales, time_group=="pre")
wilcox.test(anaerobic_bifidobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_bifidobacteriales)
# P=0.02

anaerobic_burkholderiales <- subset(anaerobic_order, Order=="Burkholderiales")
anaerobic_burkholderiales <- anaerobic_burkholderiales[order(anaerobic_burkholderiales$study_id),]
anaerobic_burkholderiales <- anaerobic_burkholderiales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_burkholderiales <- subset(anaerobic_burkholderiales, time_group=="pre")
wilcox.test(anaerobic_burkholderiales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_burkholderiales)
# P=0.58

anaerobic_campylobacterales <- subset(anaerobic_order, Order=="Campylobacterales")
anaerobic_campylobacterales <- anaerobic_campylobacterales[order(anaerobic_campylobacterales$study_id),]
anaerobic_campylobacterales <- anaerobic_campylobacterales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_campylobacterales <- subset(anaerobic_campylobacterales, time_group=="pre")
wilcox.test(anaerobic_campylobacterales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_campylobacterales)
# P>0.99

anaerobic_clostridiales <- subset(anaerobic_order, Order=="Clostridiales")
anaerobic_clostridiales <- anaerobic_clostridiales[order(anaerobic_clostridiales$study_id),]
anaerobic_clostridiales <- anaerobic_clostridiales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_clostridiales <- subset(anaerobic_clostridiales, time_group=="pre")
wilcox.test(anaerobic_clostridiales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_clostridiales)
# P=0.008

anaerobic_coriobacteriales <- subset(anaerobic_order, Order=="Coriobacteriales")
anaerobic_coriobacteriales <- anaerobic_coriobacteriales[order(anaerobic_coriobacteriales$study_id),]
anaerobic_coriobacteriales <- anaerobic_coriobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_coriobacteriales <- subset(anaerobic_coriobacteriales, time_group=="pre")
wilcox.test(anaerobic_coriobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_coriobacteriales)
# P=0.83

anaerobic_desulfovibrionales <- subset(anaerobic_order, Order=="Desulfovibrionales")
anaerobic_desulfovibrionales <- anaerobic_desulfovibrionales[order(anaerobic_desulfovibrionales$study_id),]
anaerobic_desulfovibrionales <- anaerobic_desulfovibrionales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_desulfovibrionales <- subset(anaerobic_desulfovibrionales, time_group=="pre")
wilcox.test(anaerobic_desulfovibrionales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_desulfovibrionales)
# P>0.99

anaerobic_enterobacteriales <- subset(anaerobic_order, Order=="Enterobacteriales")
anaerobic_enterobacteriales <- anaerobic_enterobacteriales[order(anaerobic_enterobacteriales$study_id),]
anaerobic_enterobacteriales <- anaerobic_enterobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_enterobacteriales <- subset(anaerobic_enterobacteriales, time_group=="pre")
wilcox.test(anaerobic_enterobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_enterobacteriales)
# P=0.15

anaerobic_erysipelotrichales <- subset(anaerobic_order, Order=="Erysipelotrichales")
anaerobic_erysipelotrichales <- anaerobic_erysipelotrichales[order(anaerobic_erysipelotrichales$study_id),]
anaerobic_erysipelotrichales <- anaerobic_erysipelotrichales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_erysipelotrichales <- subset(anaerobic_erysipelotrichales, time_group=="pre")
wilcox.test(anaerobic_erysipelotrichales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_erysipelotrichales)
# P=0.64

anaerobic_fusobacteriales <- subset(anaerobic_order, Order=="Fusobacteriales")
anaerobic_fusobacteriales <- anaerobic_fusobacteriales[order(anaerobic_fusobacteriales$study_id),]
anaerobic_fusobacteriales <- anaerobic_fusobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_fusobacteriales <- subset(anaerobic_fusobacteriales, time_group=="pre")
wilcox.test(anaerobic_fusobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_fusobacteriales)
# P>0.99

anaerobic_lactobacillales <- subset(anaerobic_order, Order=="Lactobacillales")
anaerobic_lactobacillales <- anaerobic_lactobacillales[order(anaerobic_lactobacillales$study_id),]
anaerobic_lactobacillales <- anaerobic_lactobacillales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_lactobacillales <- subset(anaerobic_lactobacillales, time_group=="pre")
wilcox.test(anaerobic_lactobacillales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_lactobacillales)
# P=0.25

anaerobic_pasteurellales <- subset(anaerobic_order, Order=="Pasteurellales")
anaerobic_pasteurellales <- anaerobic_pasteurellales[order(anaerobic_pasteurellales$study_id),]
anaerobic_pasteurellales <- anaerobic_pasteurellales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_pasteurellales <- subset(anaerobic_pasteurellales, time_group=="pre")
wilcox.test(anaerobic_pasteurellales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_pasteurellales)
# P=0.37

anaerobic_pseudomonadales <- subset(anaerobic_order, Order=="Pseudomonadales")
anaerobic_pseudomonadales <- anaerobic_pseudomonadales[order(anaerobic_pseudomonadales$study_id),]
anaerobic_pseudomonadales <- anaerobic_pseudomonadales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_pseudomonadales <- subset(anaerobic_pseudomonadales, time_group=="pre")
wilcox.test(anaerobic_pseudomonadales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_pseudomonadales)
# P=0.37

anaerobic_selenomonadales <- subset(anaerobic_order, Order=="Selenomonadales")
anaerobic_selenomonadales <- anaerobic_selenomonadales[order(anaerobic_selenomonadales$study_id),]
anaerobic_selenomonadales <- anaerobic_selenomonadales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_selenomonadales <- subset(anaerobic_selenomonadales, time_group=="pre")
wilcox.test(anaerobic_selenomonadales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_selenomonadales)
# P=0.64

anaerobic_verrucomicrobiales <- subset(anaerobic_order, Order=="Verrucomicrobiales")
anaerobic_verrucomicrobiales <- anaerobic_verrucomicrobiales[order(anaerobic_verrucomicrobiales$study_id),]
anaerobic_verrucomicrobiales <- anaerobic_verrucomicrobiales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_verrucomicrobiales <- subset(anaerobic_verrucomicrobiales, time_group=="pre")
wilcox.test(anaerobic_verrucomicrobiales$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_verrucomicrobiales)
# P>0.99

anaerobic_xanthomonadales <- subset(anaerobic_order, Order=="Xanthomonadales")
remove(anaerobic_xanthomonadales)
# Only one sample with reads assigned to Xanthomonadales in anaerobic group

# THE FOLLOWING ORDERS DECLINED IN ABUNDANCE IN PAIRED SAMPLES WITH ANAEROBIC ANTIBIOTICS: BACTEROIDALES, BIFIDOBACTERIALES, CLOSTRIDIALES

######################################################################################################################################
# TEST GENUS-LEVEL DIFFERENCES IN ORDERS THAT WERE DIFFERENTIALLY ABUNDANT BEFORE/DURING ANAEROBIC ANTIBIOTICS

# Compare composition at the GENUS level for Orders that differed after antibiotic exposure
phy.anaerobic.genus <- tax_glom(phy.anaerobic, taxrank = 'Genus')

phy.anaerobic.genus <- transform_sample_counts(phy.anaerobic.genus, function(Abundance) Abundance/sum(Abundance))
sample_sums(phy.anaerobic.genus)  
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling
anaerobic_genus <- psmelt(phy.anaerobic.genus)
anaerobic_genus <- anaerobic_genus[,-1]

# BACTEROIDALES
anaerobic_bacteroidales_genus <- subset(anaerobic_genus, Order=="Bacteroidales")
tapply(anaerobic_bacteroidales_genus$Abundance, anaerobic_bacteroidales_genus$Genus, FUN=sum)
# Prevotella, Bacteroides, Parabacteroides, Alistipes

anaerobic_prevotella <- subset(anaerobic_bacteroidales_genus, Genus=="Prevotella")
anaerobic_prevotella <- anaerobic_prevotella[order(anaerobic_prevotella$study_id),]
tapply(anaerobic_prevotella$Abundance, anaerobic_prevotella$time_group, summary)
anaerobic_prevotella <- anaerobic_prevotella %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_prevotella <- subset(anaerobic_prevotella, time_group=="pre")
wilcox.test(anaerobic_prevotella$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_prevotella)
# Prevotella: P>0.99

anaerobic_bacteroides <- subset(anaerobic_bacteroidales_genus, Genus=="Bacteroides")
anaerobic_bacteroides <- anaerobic_bacteroides[order(anaerobic_bacteroides$study_id),]
tapply(anaerobic_bacteroides$Abundance, anaerobic_bacteroides$time_group, summary)
anaerobic_bacteroides <- anaerobic_bacteroides %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_bacteroides <- subset(anaerobic_bacteroides, time_group=="pre")
wilcox.test(anaerobic_bacteroides$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_bacteroides)
# Bacteroides: P=0.20

anaerobic_parabacteroides <- subset(anaerobic_bacteroidales_genus, Genus=="Parabacteroides")
anaerobic_parabacteroides <- anaerobic_parabacteroides[order(anaerobic_parabacteroides$study_id),]
tapply(anaerobic_parabacteroides$Abundance, anaerobic_parabacteroides$time_group, summary)
anaerobic_parabacteroides <- anaerobic_parabacteroides %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_parabacteroides <- subset(anaerobic_parabacteroides, time_group=="pre")
wilcox.test(anaerobic_parabacteroides$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_parabacteroides)
# Parabacteroides: P=0.06

anaerobic_alistipes <- subset(anaerobic_bacteroidales_genus, Genus=="Alistipes")
anaerobic_alistipes <- anaerobic_alistipes[order(anaerobic_alistipes$study_id),]
tapply(anaerobic_alistipes$Abundance, anaerobic_alistipes$time_group, summary)
anaerobic_alistipes <- anaerobic_alistipes %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_alistipes <- subset(anaerobic_alistipes, time_group=="pre")
wilcox.test(anaerobic_alistipes$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_alistipes)
# Alistipes: P=0.79

remove(phy.anaerobic.genus, phy.anaerobic.order, anaerobic_bacteroidales_genus)

# CLOSTRIDIALES
anaerobic_clostridiales_genus <- subset(anaerobic_genus, Order=="Clostridiales")
tapply(anaerobic_clostridiales_genus$Abundance, anaerobic_clostridiales_genus$Genus, FUN=sum)
# Blautia, Clostridium, Subdoligranulum, Ruminococcus, Peptostreptococcae_noname
# Dorea, Faecalibacterium, Eubacterium, Lachnospiraceae_noname, Oscillibacter

anaerobic_blautia <- subset(anaerobic_clostridiales_genus, Genus=="Blautia")
anaerobic_blautia <- anaerobic_blautia[order(anaerobic_blautia$study_id),]
tapply(anaerobic_blautia$Abundance, anaerobic_blautia$time_group, summary)
anaerobic_blautia <- anaerobic_blautia %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_blautia <- subset(anaerobic_blautia, time_group=="pre")
wilcox.test(anaerobic_blautia$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_blautia)
# Blautia: P=0.04

anaerobic_clostridium <- subset(anaerobic_clostridiales_genus, Genus=="Clostridium")
anaerobic_clostridium <- anaerobic_clostridium[order(anaerobic_clostridium$study_id),]
tapply(anaerobic_clostridium$Abundance, anaerobic_clostridium$time_group, summary)
anaerobic_clostridium <- anaerobic_clostridium %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_clostridium <- subset(anaerobic_clostridium, time_group=="pre")
wilcox.test(anaerobic_clostridium$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_clostridium)
# Clostridium: P=0.008

anaerobic_subdoli <- subset(anaerobic_clostridiales_genus, Genus=="Subdoligranulum")
anaerobic_subdoli <- anaerobic_subdoli[order(anaerobic_subdoli$study_id),]
tapply(anaerobic_subdoli$Abundance, anaerobic_subdoli$time_group, summary)
anaerobic_subdoli <- anaerobic_subdoli %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_subdoli <- subset(anaerobic_subdoli, time_group=="pre")
wilcox.test(anaerobic_subdoli$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_subdoli)
# Subdoligranulum: P=0.04

anaerobic_rumino <- subset(anaerobic_clostridiales_genus, Genus=="Ruminococcus")
anaerobic_rumino <- anaerobic_rumino[order(anaerobic_rumino$study_id),]
tapply(anaerobic_rumino$Abundance, anaerobic_rumino$time_group, summary)
anaerobic_rumino <- anaerobic_rumino %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_rumino <- subset(anaerobic_rumino, time_group=="pre")
wilcox.test(anaerobic_rumino$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_rumino)
# Ruminococcus: P>0.99

anaerobic_pepto <- subset(anaerobic_clostridiales_genus, Genus=="Peptostreptococcaceae_noname")
anaerobic_pepto <- anaerobic_pepto[order(anaerobic_pepto$study_id),]
tapply(anaerobic_pepto$Abundance, anaerobic_pepto$time_group, summary)
anaerobic_pepto <- anaerobic_pepto %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_pepto <- subset(anaerobic_pepto, time_group=="pre")
wilcox.test(anaerobic_pepto$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_pepto)
# Unclassified Peptostreptococcaceae: P=0.04

anaerobic_dorea <- subset(anaerobic_clostridiales_genus, Genus=="Dorea")
anaerobic_dorea <- anaerobic_dorea[order(anaerobic_dorea$study_id),]
tapply(anaerobic_dorea$Abundance, anaerobic_dorea$time_group, summary)
anaerobic_dorea <- anaerobic_dorea %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_dorea <- subset(anaerobic_dorea, time_group=="pre")
wilcox.test(anaerobic_dorea$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_dorea)
# Dorea: P=0.06

anaerobic_faec <- subset(anaerobic_clostridiales_genus, Genus=="Faecalibacterium")
anaerobic_faec <- anaerobic_faec[order(anaerobic_faec$study_id),]
tapply(anaerobic_faec$Abundance, anaerobic_faec$time_group, summary)
anaerobic_faec <- anaerobic_faec %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_faec <- subset(anaerobic_faec, time_group=="pre")
wilcox.test(anaerobic_faec$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_faec)
# Faecalibacterium: P=0.18

anaerobic_eub <- subset(anaerobic_clostridiales_genus, Genus=="Eubacterium")
anaerobic_eub <- anaerobic_eub[order(anaerobic_eub$study_id),]
tapply(anaerobic_eub$Abundance, anaerobic_eub$time_group, summary)
anaerobic_eub <- anaerobic_eub %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_eub <- subset(anaerobic_eub, time_group=="pre")
wilcox.test(anaerobic_eub$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_eub)
# Eubacterium: P=0.10

anaerobic_lachno <- subset(anaerobic_clostridiales_genus, Genus=="Lachnospiraceae_noname")
anaerobic_lachno <- anaerobic_lachno[order(anaerobic_lachno$study_id),]
tapply(anaerobic_lachno$Abundance, anaerobic_lachno$time_group, summary)
anaerobic_lachno <- anaerobic_lachno %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_lachno <- subset(anaerobic_lachno, time_group=="pre")
wilcox.test(anaerobic_lachno$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_lachno)
# Unclassified Lachnospiraceae: P=0.04

anaerobic_osc <- subset(anaerobic_clostridiales_genus, Genus=="Oscillibacter")
anaerobic_osc <- anaerobic_osc[order(anaerobic_osc$study_id),]
tapply(anaerobic_osc$Abundance, anaerobic_osc$time_group, summary)
anaerobic_osc <- anaerobic_osc %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
anaerobic_osc <- subset(anaerobic_osc, time_group=="pre")
wilcox.test(anaerobic_osc$D_Abundance, mu = 0, alternative = "two.sided")
remove(anaerobic_osc, anaerobic_clostridiales_genus)
# Oscillibacter: P=0.06

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
# Bifidobacterium: P=0.04

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
# Clostridiales: P=0.10

gvhd_bifidobacteriales <- subset(gvhd_order, Order=="Bifidobacteriales")
gvhd_bifidobacteriales <- gvhd_bifidobacteriales[order(gvhd_bifidobacteriales$study_id),]
gvhd_bifidobacteriales <- gvhd_bifidobacteriales %>%
  group_by(study_id) %>%
  mutate(D_Abundance = Abundance[time_group == "post"] - Abundance)
gvhd_bifidobacteriales <- subset(gvhd_bifidobacteriales, time_group=="pre")
wilcox.test(gvhd_bifidobacteriales$D_Abundance, mu = 0, alternative = "two.sided")
remove(gvhd_bifidobacteriales)
# Bifidobacteriales: P=0.01

remove(phy.gvhd.order, gvhd_order)

######################################################################################################################################
# HEATMAP

library(gtools)

phy.order <-  tax_glom(phy.gvhd, "Order")
taxa_names(phy.order) <-  as.character(tax_table(phy.order)[,c("Order")])

# Focus only on orders that are present in >=25% of samples
phy.order <- subset_taxa(phy.order, Order!="Bacillales" & Order!="Campylobacterales" & Order!="Fusobacteriales" & Order!="Pasteurellales" &
                           Order!="Pseudomonadales" & Order!="Verrucomicrobiales" & Order!="Xanthomonadales")


c_phy_pre <- data.frame(otu_table(subset_samples(phy.order, treatment == "cefepime" & time_group == "pre")))
names(c_phy_pre) <- paste0("Patient ", seq(1,20))
c_phy_pre <- c_phy_pre[,mixedsort(names(c_phy_pre))]

c_phy_post <-data.frame(otu_table( subset_samples(phy.order, treatment == "cefepime" & time_group == "post")))
names(c_phy_post) <- paste0("Patient ", seq(1,20))
c_phy_post <- c_phy_post[,mixedsort(names(c_phy_post))]

# Percent change
c_change <- 100*(c_phy_post - c_phy_pre)/(c_phy_pre)

# Anaerobic
a_phy_pre <- data.frame(otu_table(subset_samples(phy.order, treatment == "anaerobic" & time_group == "pre")))
names(a_phy_pre) <- paste0("Patient ", seq(21,28))
a_phy_pre <- a_phy_pre[,mixedsort(names(a_phy_pre))]

a_phy_post <-data.frame(otu_table( subset_samples(phy.order, treatment == "anaerobic" & time_group == "post")))
names(a_phy_post) <- paste0("Patient ", seq(21,28))
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

png(file="Results/Heatmap by Antibiotics.png", width = 13, height = 8, units = 'in', res = 600)
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

# Butyrate genes declined in abundance with ANAEROBIC ANTIBIOTICS 
ana_butyrate <- subset(metadata_gvhd, treatment=="Anaerobic Antibiotics")
tapply(ana_butyrate$butyrate_genes, ana_butyrate$time_group, summary)
ana_pre_butyrate <- subset(ana_butyrate, antibiotic_group=="anaerobic_pre", butyrate_genes, drop = TRUE)
ana_post_butyrate <- subset(ana_butyrate, antibiotic_group=="anaerobic_post", butyrate_genes, drop = TRUE)
wilcox.test(ana_pre_butyrate, ana_post_butyrate, paired = TRUE, alternative = "two.sided")
remove(ana_butyrate, ana_pre_butyrate, ana_post_butyrate)

# Butyrate genes were higher in pre- samples from patients who subsequently developed acute GVHD of gut/liver 
pre_butyrate <- subset(metadata_gvhd, time_group=="pre")
pre_gvhd_butyrate <- subset(pre_butyrate, agvhd_gutliver=="1", butyrate_genes, drop = TRUE)
pre_nogvhd_butyrate <- subset(pre_butyrate, agvhd_gutliver=="0", butyrate_genes, drop = TRUE)
summary(pre_gvhd_butyrate)
summary(pre_nogvhd_butyrate)
wilcox.test(pre_gvhd_butyrate, pre_nogvhd_butyrate, alternative = "two.sided")
remove(pre_butyrate, pre_gvhd_butyrate, pre_nogvhd_butyrate)

# Butyrate genes declined in abundance among patients with acute GVHD of gut/liver 
gvhd_butyrate <- subset(metadata_gvhd, agvhd_gutliver=="1")
gvhd_pre_butyrate <- subset(gvhd_butyrate, antibiotic_group=="cefepime_pre" | antibiotic_group=="anaerobic_pre", butyrate_genes, drop = TRUE)
gvhd_post_butyrate <- subset(gvhd_butyrate, antibiotic_group=="cefepime_post" | antibiotic_group=="anaerobic_post", butyrate_genes, drop = TRUE)
summary(gvhd_pre_butyrate)
summary(gvhd_post_butyrate)
wilcox.test(gvhd_pre_butyrate, gvhd_post_butyrate, paired = TRUE, alternative = "two.sided")
remove(gvhd_butyrate, gvhd_pre_butyrate, gvhd_post_butyrate)

# Butyrate genes DID NOT decline in abundance among patients without acute GVHD of gut/liver
no_gvhd_butyrate <- subset(metadata_gvhd, agvhd_gutliver=="0")
no_gvhd_pre_butyrate <- subset(no_gvhd_butyrate, antibiotic_group=="cefepime_pre" | antibiotic_group=="anaerobic_pre", butyrate_genes, drop = TRUE)
no_gvhd_post_butyrate <- subset(no_gvhd_butyrate, antibiotic_group=="cefepime_post" | antibiotic_group=="anaerobic_post", butyrate_genes, drop = TRUE)
summary(no_gvhd_pre_butyrate)
summary(no_gvhd_post_butyrate)
wilcox.test(no_gvhd_pre_butyrate, no_gvhd_post_butyrate, paired = TRUE, alternative = "two.sided")
remove(no_gvhd_butyrate, no_gvhd_pre_butyrate, no_gvhd_post_butyrate)

metadata_gvhd$treatment <- factor(metadata_gvhd$treatment, levels=c("Cefepime", "Anaerobic Antibiotics"))
levels(metadata_gvhd$agvhd_gutliver) <- c(levels(metadata_gvhd$agvhd_gutliver), "No")
levels(metadata_gvhd$agvhd_gutliver) <- c(levels(metadata_gvhd$agvhd_gutliver), "Yes")
metadata_gvhd$agvhd_gutliver[metadata_gvhd$agvhd_gutliver == "0"] <- "No"
metadata_gvhd$agvhd_gutliver[metadata_gvhd$agvhd_gutliver == "1"] <- "Yes"

butyrate <- ggplot(metadata_gvhd, aes(x=time_group, y=butyrate_genes)) + 
  geom_point(aes(group=time_group, color=agvhd_gutliver)) + geom_line(aes(group=study_id, color=agvhd_gutliver)) +
  theme(axis.title.y = element_text(size=10), axis.text.y = element_text(size=9), legend.text=element_text(size=8), legend.title=element_text(size=9), 
        axis.title.x = element_text(size=10), axis.text.x = element_text(size=9),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text.y=element_text(size=10), legend.position="right") + ylab("Butyrate Gene Relative Abundance (RPKM)") + xlab("Sample Timing") +
  scale_x_discrete(labels=c("Before", "During")) + scale_color_manual(values = c("black", "red")) +
  guides(color=guide_legend("Acute Gut/Liver GVHD")) +
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
  geom_point(aes(group=time_group, color=agvhd_gutliver)) + geom_line(aes(group=study_id, color=agvhd_gutliver)) +
  theme(axis.title.y = element_text(size=10), axis.text.y = element_text(size=9), legend.text=element_text(size=8), legend.title=element_text(size=9), 
        axis.title.x = element_text(size=10), axis.text.x = element_text(size=9),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text.y=element_text(size=10), legend.position="right") + ylab("Clostridiales Relative Abundance") + xlab("Sample Timing") +
  scale_x_discrete(labels=c("Before", "During")) + scale_color_manual(values = c("black", "red")) +
  guides(color=guide_legend("Acute Gut/Liver GVHD")) +
  facet_wrap(~treatment) 

png(file="Results/Clostridiales Abundance by Antibiotics.png", width = 8, height = 3.5, units = 'in', res = 600)
plot(clostridiales)
dev.off()
remove(clostridiales, clostridiales_df)