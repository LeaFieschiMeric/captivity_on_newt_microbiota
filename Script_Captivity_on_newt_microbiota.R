##################################################################################
##   Script for the statistical analysis of the data presented in the article   ##
##    "Strong restructuration of skin microbiota in amphibians under captive    ##
##          management protocols challenges their conservation ex-situ".        ##                              ##
##################################################################################

# Data initially processed through DADA2 (https://benjjneb.github.io/dada2/tutorial.html)
# aligned with SILVA trained on EMPO (animal surface; https://earthmicrobiome.org/protocols-and-standards/empo/)
# and preprocessed using phyloseq (https://cran.r-project.org/web/packages/vegan/vegan.pdf).


### ------------------------------- Libraries loading + data prep ----------------

library(ggplot2)
library(phyloseq)
library(dplyr)
library(ggpubr)
library(btools)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(gratia)
library(mgcv)
library(magrittr)
library(data.table)
library(ape)
library(rbiom)
library(dendextend)
library(vegan)
library(ade4)
library(forcats)
library(dunn.test)
library(pairwiseAdonis)
library(MiscMetabar)
library(ggVennDiagram)
library(MicrobiotaProcess)
library(microbiomeMarker)
library(speedyseq)
library(multcomp)

load("Phyloseq_newts.RData")                       
load("Phyloseq_environmental.RData")

# Wild samples (n~40) -> variation between species at D0
initial_ps = subset_samples(filtered_ps, (Month_captivity == 0))
initial_ps = filter_taxa(initial_ps, function(x) sum(x) > 0, TRUE)
initial_inhib_ps = subset_taxa(initial_ps, Bd_inhibition == "inhibitory")

# Month1 samples (n~40)
M1_ps = subset_samples(filtered_ps, (Month_captivity == 1))
M1_ps = filter_taxa(M1_ps, function(x) sum(x) > 0, TRUE)

# Wild + Month1 samples (n~80) -> transfer into captivity 
M0M1_ps = subset_samples(filtered_ps, (Month_captivity == 0 | Month_captivity == 1))
M0M1_ps = filter_taxa(M0M1_ps, function(x) sum(x) > 0, TRUE) 

# Wild + Month1 samples PER SPECIES (n~40) -> transfer into captivity (replace with "Pal" for other species)
M0M1_alp_ps = subset_samples(M0M1_ps, (Species == "Alp"))
M0M1_alp_ps = filter_taxa(M0M1_alp_ps, function(x) sum(x) > 0, TRUE) 
M0M1_inhib_alp_ps = subset_taxa(M0M1_alp_ps, Bd_inhibition == "inhibitory")

# All samples PER SPECIES (n~220) -> variation throughout 10 months of experiment (replace with "Pal" for other species)
alp_ps = subset_samples(filtered_ps, (Species == "Alp")) 
alp_ps = filter_taxa(alp_ps, function(x) sum(x) > 0, TRUE)

# Environmental samples -> water poured in the aquaria (replace with "Cork" for other contamination)
ctrl_water_ps = subset_samples(raref_envtal_ps, (Type == "Water")) 
ctrl_water_ps = filter_taxa(ctrl_water_ps, function(x) sum(x) > 0, TRUE) 



### ------------------------------- Venn diagramm on total ASVs / pop ------------

# FIG. S2: Total phylotypes initially shared between species
vennlist <- get_vennlist(initial_ps, factorNames="Species")
ggVennDiagram(vennlist, label = "count", category.names = c("Alpine","Palmate"))
vennlist <- get_vennlist(initial_inhib_ps, factorNames="Species")
ggVennDiagram(vennlist, label = "count", category.names = c("Alpine","Palmate"))

# FIG. S4: Total phylotypes retained over the transfer into captivity, per species (replace by "Pal" in "M0M1_alp_ps" for other species)
vennlist <- get_vennlist(M0M1_alp_ps, factorNames="Month_captivity")
ggVennDiagram(vennlist, label = "count", category.names = c("Wild","Captive (1 Month)"))
vennlist <- get_vennlist(M0M1_inhib_alp_ps, factorNames="Month_captivity")
ggVennDiagram(vennlist, label = "count", category.names = c("Wild","Captive (1 Month)"))

# FIG. S5: Total phylotypes introduced by environmental contamination (replace by "Cork" and overwintering samples for other contamination)
vennlist <- get_vennlist(initial_ps, factorNames="Species")
wild <- list(c(vennlist$Alp,vennlist$Pal))
vennlist <- get_vennlist(M1_ps, factorNames="Phase")
M1 <- list(c(vennlist))
vennlist <- get_vennlist(ctrl_water_ps, factorNames="Species")
ctrl <- list(c(vennlist))
merged_list <- c(wild, M1, ctrl)
ggVennDiagram(merged_list, label = "count",
                     category.names = c("Wild","M1", "Water"))



### ------------------------------- Relative abundance plots ---------------------

# FIG. 1: Relative abundance per month and per species, at the phylum level (replace by "Pal" in "alp_ps" for other species)
ps_phylum = tax_glom(alp_ps,taxrank = "Phylum")
ps_phylum = ps_phylum %>% transform_sample_counts(function(x) {x/sum(x)} )
ps_phylum = psmelt(ps_phylum)
ps_barplot = ps_phylum %>%
  group_by(Phylum, Month_captivity) %>%
  summarise(nb = sum(Abundance))
ps_barplot$Phylum = as.character(ps_barplot$Phylum)
ps_barplot$Phylum[ps_barplot$nb < 0.005] = "Phylum < 0.5% abund."
ggplot(ps_barplot, aes(x = Month_captivity, y = nb, fill = Phylum, na.rm = TRUE)) +
  geom_bar(stat="identity",na.rm = TRUE, position="fill") +
  ylab("Relative Abundance \n") + xlab("Month") + 
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

# Addition of environmental samples
envtal_ps_phylum = tax_glom(raref_envtal_ps,taxrank = "Phylum")
envtal_ps_phylum = envtal_ps_phylum %>% transform_sample_counts(function(x) {x/sum(x)} )
envtal_ps_phylum = psmelt(envtal_ps_phylum)
envtal_ps_barplot = envtal_ps_phylum %>%
  group_by(Phylum, Type) %>%
  summarise(nb = sum(Abundance))
envtal_ps_barplot$Phylum = as.character(envtal_ps_barplot$Phylum)
envtal_ps_barplot$Phylum[envtal_ps_barplot$nb < 0.005] = "Phylum < 0.5% abund."
ggplot(envtal_ps_barplot %>% filter(Type != 'Spray' & Type != 'Aquaria') %>%
                mutate(Type = factor(Type, levels = c('Water', 'Cork'))),
              aes(x = Type, y = nb, fill = Phylum, na.rm = TRUE)) +
  geom_bar(stat="identity",na.rm = TRUE, position="fill") +
  ylab("Relative Abundance \n")

# FIG. 3: Relative abundance of Bd-inhibitory phylotypes (replace by "Pal" in "alp_ps" for other species)
taxo_alp <- data.frame(as(tax_table(alp_ps), "matrix"))
temp_alp_ps <- alp_ps %>% transmute_tax_table(Phylum_activity = Phylum_activity, 
                                              Bd_inhib_phylum = Bd_inhib_phylum)
taxo_alp <- data.frame(as(tax_table(temp_alp_ps), "matrix"))
ps_phylum = tax_glom(temp_alp_ps, taxrank = "Bd_inhib_phylum")
ps_phylum = ps_phylum %>% transform_sample_counts(function(x) {x/sum(x)} )
ps_phylum = psmelt(ps_phylum)
ps_barplot = ps_phylum %>%
  group_by(Bd_inhib_phylum, Month_captivity) %>%
  summarise(nb = sum(Abundance))
ps_barplot$Bd_inhib_phylum = factor(ps_barplot$Bd_inhib_phylum, levels = c('enhancing', 'not inhibitory', 'unknown', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes', 'Proteobacteria'))
ggplot(ps_barplot, aes(x = Month_captivity, y = nb, fill = Bd_inhib_phylum, na.rm = TRUE)) +
  geom_bar(stat="identity",na.rm = TRUE, position="fill") +
  ylab("Relative Abundance \n") + xlab("Status") +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))



### ------------------------------- Alpha diversity ------------------------------

alpha = estimate_richness(filtered_ps, split = TRUE, measures = NULL) 
alpha = rownames_to_column(alpha, var = "Sample_ID") 
metadata <- read.csv("metadata.csv", sep=";")
metadata <- metadata %>% mutate(Sample_ID = X)
alpha = alpha %>% inner_join(metadata, by = "Sample_ID")
lapply(alpha, class)
alpha <- alpha %>% dplyr::select(-X, -is.neg) %>%
  mutate(Sample_ID = factor(Sample_ID),
         Individual = factor(Individual),
         Sex = factor(Sex),
         Species = factor(Species),
         Month_captivity = factor(Month_captivity),
         Status = factor(Status),
         Phase = factor(Phase),
         SVL_cm = as.numeric(SVL_cm),
         Period = factor(Period),
         Aquarium_ID = factor(Aquarium_ID))

# FIG. 2: Monthly variation in alpha-diversity (replace with "y = Shannon" for the other figure)
ggplot(alpha, aes(x = Month_captivity, y = Chao1, fill = Species)) +
  geom_boxplot() +
  labs(y = 'Chao1 index', x = 'Month in captivity')

# Wild samples - differences in Chao1 between species (replace with "Shannon" for the other alpha-div measure)
mod1 = aov(Chao1 ~ Species * Sex + SVL_cm, data = alpha %>% filter(Month_captivity == "0")) 
shapiro.test(residuals(object = mod1))
bartlett.test(Chao1 ~ Species, data = alpha %>% filter(Month_captivity == "0")) 
bartlett.test(Chao1 ~ Sex, data = alpha %>% filter(Month_captivity == "0")) 
summary(mod1)
mod1 = aov(Chao1 ~ Species + Sex + SVL_cm, data = alpha %>% filter(Month_captivity == "0")) #rerun without the unsignificant interaction
shapiro.test(residuals(object = mod1))
bartlett.test(Chao1 ~ Species, data = alpha %>% filter(Month_captivity == "0")) 
bartlett.test(Chao1 ~ Sex, data = alpha %>% filter(Month_captivity == "0")) 
summary(mod1)

# Wild + Month1 samples - differences in Chao1 over transfer into captivity (replace with "Shannon" for the other alpha-div measure)
library(lme4)
mod2 = lmer(Chao1 ~ Species * Status + (1|Individual), data = alpha %>% filter(Month_captivity == "0" | Month_captivity == "1"))
shapiro.test(residuals(object = mod2)) 
mod2 = lmer(log(Chao1) ~ Species * Status + (1|Individual), data = alpha %>% filter(Month_captivity == "0" | Month_captivity == "1")) #log-transformation for normality of residuals
shapiro.test(residuals(object = mod2))
anova(mod2)
mod2 = lmer(log(Chao1) ~ Species + Status + (1|Individual), data = alpha %>% filter(Month_captivity == "0" | Month_captivity == "1")) #rerun without the unsignificant interaction
anova(mod2)

# All samples - differences in Chao1 over phase-shifts (replace with "Shannon" for the other alpha-div measure)
mod3 = lmer(Chao1 ~ Species * Period + (1|Individual)+ (1|Aquarium_ID) + (1|Period/Month_captivity), data = alpha)
mod3 = lmer(log(Chao1) ~ Species * Period + (1|Individual)+ (1|Aquarium_ID) + (1|Period/Month_captivity), data = alpha) #simplify to allow convergence
mod3 = lmer(Chao1 ~ Species * Period + (1|Individual)+ (1|Aquarium_ID) + (1|Month_captivity), data = alpha)             #simplify to allow convergence
shapiro.test(residuals(object = mod3))  
mod3 = lmer(log(Chao1) ~ Species * Period + (1|Individual)+ (1|Aquarium_ID) + (1|Month_captivity), data = alpha)  #log-transformation for normality of residuals
shapiro.test(residuals(object = mod3))  
anova(mod3)
mod3 = lmer(log(Chao1) ~ Species + Period + (1|Individual)+ (1|Aquarium_ID) + (1|Month_captivity), data = alpha)  #rerun without the unsignificant interaction
anova(mod3)
emmeans(mod3, pairwise ~ Period) # estimated marginal means and contrasts (replace with "Period | Species" for Shannon)



### ----------------------------- Beta diversity ---------------------------------

pick_new_outgroup = function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape")
  treeDT =
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  new.outgroup = treeDT[which.max(length)]$id
  return(new.outgroup) }
my.tree = phy_tree(filtered_ps)
out.group = pick_new_outgroup(my.tree)
out.group # value = db19a3ea000384e41633535b013c0743
new.tree = ape::root(my.tree, outgroup=out.group, resolve.root=TRUE)
phy_tree(filtered_ps) = new.tree

# Wild samples - differences in beta-div between species
rbiom_weighted = rbiom::unifrac(otu_table(initial_ps), weighted=TRUE, tree=new.tree)
ini_meta <- data.frame(sample_data(initial_ps))
adonis2(rbiom_weighted ~ Species*Sex + SVL_cm, data = ini_meta, permutations = 9999)
adonis2(rbiom_weighted ~ Species + Sex + SVL_cm, data = ini_meta, permutations = 9999) #rerun without the unsignificant interaction
disp = betadisper(rbiom_weighted, ini_meta$Species)
permutest(disp, pairwise = TRUE, permutations = 9999)   
disp = betadisper(rbiom_weighted, ini_meta$Sex)
permutest(disp, pairwise = TRUE, permutations = 9999)   

# FIG. S3: Beta diversity of wild samples
ord = ordinate(initial_ps, method = "PCoA", distance = rbiom_weighted)
plot_ordination(initial_ps, ord, type='samples', color="Species", shape="Sex", justDF = F) +
  stat_ellipse(aes(lty=factor(Sex))) +
  labs(color = "Species", shape = "Sex") 

# Wild + Month1 samples - change in beta-div over the transfer into captivity
rbiom_weighted = rbiom::unifrac(otu_table(M0M1_ps), weighted=TRUE, tree=new.tree)
tr_meta <- data.frame(sample_data(M0M1_ps))
adonis2(rbiom_weighted ~ Species * Status + Individual, data = tr_meta, permutations = 9999)
disp = betadisper(rbiom_weighted, tr_meta$Species)
permutest(disp, pairwise = TRUE, permutations = 9999)   
disp = betadisper(rbiom_weighted, tr_meta$Status)
permutest(disp, pairwise = TRUE, permutations = 9999)   

# FIG. 4A: Beta diversity of wild and captive (M1) samples
ord = ordinate(M0M1_ps, method = "PCoA", distance = rbiom_weighted)
plot_ordination(M0M1_ps, ord, type='samples', color="Species", shape="Status", justDF = F) +
  stat_ellipse(aes(lty=factor(Status))) +
  labs(color = "Species", shape = "Status")

# All samples - change in beta-div over phase-shifts
filtered_ps <- filtered_ps %>% mutate_sample_data(Period = case_when(Period == "Hibernation" ~ "Overwintering", TRUE ~ Period)) # first, split Wild vs. captive aquatic samples
sample_data(filtered_ps)$Phase <- factor(sample_data(filtered_ps)$Period, levels = c("Wild", "Aquatic1", "Overwintering", "Aquatic2"))
rbiom_weighted = rbiom::unifrac(otu_table(filtered_ps), weighted=TRUE, tree=new.tree)
beta_meta <- data.frame(sample_data(filtered_ps))
adonis2(rbiom_weighted ~ Species * Phase + Individual, data = beta_meta, permutations = 9999)
disp = betadisper(rbiom_weighted, beta_meta$Species)
permutest(disp, pairwise = TRUE, permutations = 9999)   
disp = betadisper(rbiom_weighted, beta_meta$Phase)
permutest(disp, pairwise = TRUE, permutations = 9999)   

# Species-specific change in beta-div over phase-shifts (Table S6; replace with "Pal" for the other species)
sample_data(alp_ps)$Period <- factor(sample_data(alp_ps)$Period, levels = c("Wild", "Aquatic1", "Hibernation", "Aquatic2")) 
rbiom_alp = rbiom::unifrac(otu_table(alp_ps),  weighted=TRUE, tree=new.tree)
alp_meta <- data.frame(sample_data(alp_ps))
pairwise.adonis2(rbiom_alp ~ Period, data = alp_meta, permutations = 9999)

# FIG. 4B: Average beta diversity of individual samples per phase
ord = ordinate(filtered_ps, method = "PCoA", distance = rbiom_weighted)
plot_ordination(filtered_ps, ord, type='samples', color="Species", shape="Phase", justDF = F) +
  stat_ellipse(aes(lty=factor(Phase))) +
  labs(color = "Species", shape = "Phase")



### ----------------- Libraries + datasets for differential abundance analysis ---

library(DESeq2)
load("Phyloseq_newts_unrarefied.RData") 
load("Phyloseq_environmental_unrarefied.RData")

# Wild samples (n~40) -> variation between species at D0
unrar_initial_ps = subset_samples(unrar_ps2, (Month_captivity == 0))
unrar_initial_ps = filter_taxa(unrar_initial_ps, function(x) sum(x) > 0, TRUE) 

# Wild + Month1 samples PER SPECIES (n~40) -> transfer into captivity
unrar_M0M1_Alp_ps = subset_samples(unrar_ps2, (Month_captivity == 0 & Species == "Alp" |
                                                 Month_captivity == 1 & Species == "Alp")) #replace with "Pal"
unrar_M0M1_Alp_ps = filter_taxa(unrar_M0M1_Alp_ps, function(x) sum(x) > 0, TRUE) 



### ----------------------------- Differential abundance analysis ----------------

library(Rmisc)

# Wild samples - differentially abundant phylotypes between species in the wild
deseq_env <- phyloseq_to_deseq2(unrar_initial_ps, ~ Species) 
deseq_env2 = DESeq(deseq_env, test= "Wald", fitType="parametric",sfType="poscounts") 
resultsNames(deseq_env2)                                                             
res <- results(deseq_env2, name = "Species_Pal_vs_Alp", alpha = 0.01, pAdjustMethod="BH")
res = res[order(res$padj, na.last = NA), ]     
sigtab01 = res[(res$padj < 0.05), ]           
sigtab01 = cbind(as(sigtab01, "data.frame"), as(tax_table(unrar_initial_ps)[rownames(sigtab01), ], "matrix"))
sigtab01 = sigtab01[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class","Order", "Family", "Genus", "Species")]
possigtab01 = sigtab01[(sigtab01$log2FoldChange >= 0) & (sigtab01$padj <= 0.05),]
negsigtab01 = sigtab01[(sigtab01$log2FoldChange <= 0) & (sigtab01$padj <= 0.05),]
possigtab01[["change"]] <- (rep(c("Increase"),dim(possigtab01)[1]))
negsigtab01[["change"]] <- (rep(c("Decrease"),dim(negsigtab01)[1]))
contrasts = rbind(possigtab01,negsigtab01)
#write.csv(contrasts, file = "Deseq2_contrasts_Initial_by_Species.csv") 

# Wild + Month1 samples - differentially abundant phylotypes over transfer into captivity (replace with "Pal" for the other species)
deseq_env <- phyloseq_to_deseq2(unrar_M0M1_Alp_ps, ~ Month)                          
deseq_env2 = DESeq(deseq_env, test= "Wald", fitType="parametric",sfType="poscounts")
resultsNames(deseq_env2) 
res <- results(deseq_env2, name = "Month_M1_vs_M0", alpha = 0.01, pAdjustMethod="BH")
res = res[order(res$padj, na.last = NA), ]    
sigtab01 = res[(res$padj < 0.05), ]           
summary(sigtab01)                            
sigtab01 = cbind(as(sigtab01, "data.frame"), as(tax_table(unrar_M0M1_Alp_ps)[rownames(sigtab01), ], "matrix"))
sigtab01 = sigtab01[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class","Order", "Family", "Genus", "Species", "Bd_inhibition")]
possigtab01 = sigtab01[(sigtab01$log2FoldChange >= 0) & (sigtab01$padj <= 0.05),]
negsigtab01 = sigtab01[(sigtab01$log2FoldChange <= 0) & (sigtab01$padj <= 0.05),]
possigtab01[["change"]] <- (rep(c("Increase"),dim(possigtab01)[1]))
negsigtab01[["change"]] <- (rep(c("Decrease"),dim(negsigtab01)[1]))
contrasts = rbind(possigtab01,negsigtab01)
contrasts
#write.csv(contrasts, file = "Deseq2_contrasts_M0M1_Alp.csv")  

#FIG.5A: Differencially abundant phylotypes over the transfer from the wild into captivity
contrasts_alp <- read.csv("Deseq2_contrasts_M0M1_Alp.csv")
contrasts_alp <- contrasts_alp %>% dplyr::rename("ASV" = "X")
sub <- subset(otu_table(unrar_M0M1_Alp_ps), rownames(otu_table(unrar_M0M1_Alp_ps)) %in% contrasts_alp$ASV)
deseqphyseq <- merge_phyloseq(sub, tax_table(unrar_M0M1_Alp_ps), sample_data(unrar_M0M1_Alp_ps), phy_tree(unrar_M0M1_Alp_ps))
deseqphyseq1 = transform_sample_counts(deseqphyseq, log10)
dp.melt = psmelt(deseqphyseq1)
dp.melt$Abundance = as.numeric(dp.melt$Abundance)
dp.melt[dp.melt == -Inf] <- 0  
x = tapply(contrasts_alp$log2FoldChange, contrasts_alp$Phylum, function(x) max(x)) 
x = sort(x, TRUE)
contrasts_alp$Phylum = factor(as.character(contrasts_alp$Phylum), levels = names(x))
x = tapply(contrasts_alp$log2FoldChange, contrasts_alp$Genus, function(x) max(x))  
x = sort(x, TRUE)
contrasts_alp$Genus = factor(as.character(contrasts_alp$Genus), levels = names(x))
contrasts_alp <- contrasts_alp %>% mutate(Bd_inhibition = case_when(Bd_inhibition == "inhibitory" ~ "inhibitory",
                                                                    TRUE ~ "not inhibitory"))
contrasts_pal <- read.csv("Deseq2_contrasts_M0M1_Pal.csv")
contrasts_pal <- contrasts_pal %>% dplyr::rename("ASV" = "X")
sub <- subset(otu_table(unrar_M0M1_Pal_ps), rownames(otu_table(unrar_M0M1_Pal_ps)) %in% contrasts_pal$ASV)
deseqphyseq <- merge_phyloseq(sub, tax_table(unrar_M0M1_Pal_ps), sample_data(unrar_M0M1_Pal_ps), phy_tree(unrar_M0M1_Pal_ps))
deseqphyseq1 = transform_sample_counts(deseqphyseq, log10)
dp.melt = psmelt(deseqphyseq1)
dp.melt$Abundance = as.numeric(dp.melt$Abundance)
dp.melt[dp.melt == -Inf] <- 0  
x = tapply(contrasts_pal$log2FoldChange, contrasts_pal$Phylum, function(x) max(x)) 
x = sort(x, TRUE)
contrasts_pal$Phylum = factor(as.character(contrasts_pal$Phylum), levels = names(x))
x = tapply(contrasts_pal$log2FoldChange, contrasts_pal$Genus, function(x) max(x))  
x = sort(x, TRUE)
contrasts_pal$Genus = factor(as.character(contrasts_pal$Genus), levels = names(x))
contrasts_pal <- contrasts_pal %>% mutate(Bd_inhibition = case_when(Bd_inhibition == "inhibitory" ~ "inhibitory",
                                                                    TRUE ~ "not inhibitory"))
contrasts_alp <- contrasts_alp %>% mutate(Species = "Alpine newts")
contrasts_pal <- contrasts_pal %>% mutate(Species = "Palmate newts")
contrasts <- full_join(contrasts_alp,contrasts_pal)
ggplot(contrasts %>% arrange(Bd_inhibition, decreasing = TRUE),
       aes(x = Phylum, y = log2FoldChange, color = Bd_inhibition)) +
  geom_point(aes(shape = Species, fill = Bd_inhibition), size = 3, 
             position = position_jitter()) +
  scale_shape_manual(values = c(18, 16)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
  scale_x_discrete(limits = c("Actinobacteriota", "Bacteroidota", "Campilobacterota", "Cyanobacteria",
    "Desulfobacterota", "Firmicutes", "Fusobacteriota", "Gemmatimonadota", "Myxococcota",
     "Patescibacteria", "Planctomycetota", "Proteobacteria", "Verrucomicrobiota"))



#... Apply the same code for the rest of the differential abundance analysis, using the following datasets:

# 1st aquatic phase + Overwintering samples PER SPECIES (n~140) -> variation over phase-shifts
unrar_P1H_Alp_ps = subset_samples(unrar_ps2, (Species == "Alp" & Status != "Wild" & Phase == "Aquatic_1" |
                                                Species == "Alp" & Phase == "Hibernation")) 
unrar_P1H_Alp_ps = filter_taxa(unrar_P1H_Alp_ps, function(x) sum(x) > 0, TRUE) 

# Overwintering + 2nd aquatic phase samples PER SPECIES (n~120) -> variation over phase-shifts
unrar_HP2_Alp_ps = subset_samples(unrar_ps2, (Species == "Alp" & Phase == "Hibernation" |
                                                Species == "Alp" & Phase == "Aquatic_2")) 
unrar_HP2_Alp_ps = filter_taxa(unrar_HP2_Alp_ps, function(x) sum(x) > 0, TRUE) 

# Wild + Month10 samples PER SPECIES (n~40) -> long-term effect of captivity
unrar_1stlast_Alp_ps = subset_samples(unrar_ps2, (Species == "Alp" & Month_captivity == "0" |
                                                    Species == "Alp" & Month_captivity == "10")) 
unrar_1stlast_Alp_ps = filter_taxa(unrar_1stlast_Alp_ps, function(x) sum(x) > 0, TRUE) 
