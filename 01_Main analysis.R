#### Main analysis & drawing figures ####

# List of contents
  # 1. Setting
  # 2. Input (& pre-processing)
  # 3. PCoA & Heatmap
  # 4. Alpha-diversity  & human read fraction
  # 5. F. nucleatum-based classfication
  # 6. IBD/non-IBD marker
  # 7. Microbial dissimilarity (Intra- or Inter-personal)
  # 8.  F. nucleatum-Prior/Posterior analysis
  # 9. Logistic regression models in server
  # 10. Modeling results and application
  # 11. K-mean clustering of microbes
  # 12. Replicate analysis
  # 13. Probiotic ingestoin & Bifidobacterium

#### 1. Setting ####

## Requirements
options(java.parameters = "-Xmx32g", stringsAsFactors = F) 
setwd(dir = "D:/LSG/Data/Microbiome/2019/190814 Figure_1526 samples/")
library(ggplot2); library(reshape2); library(MASS); library(ggpubr); library(vegan); library(ggforce); library(ggsignif); library(ggthemes) ; library("extrafont") ; library(Hmisc)
library(gtable); library(grid) ; library(caret) ; library(ROCR) ; library(pROC) ;library(plotROC) ; library(RJSONIO) ; library(reshape) ;library(NbClust) ;library("RColorBrewer")


## Comparison set
comp_sex <- list(c("Male", "Female")) ; comp_con <- list(c("IBD", "nonIBD")) ; comp_con_3 <- list(c("Random", "nonIBD"), c("nonIBD", "IBD"), c("Random", "IBD"))
comp_type <- list(c("UC", "CD"), c("UC", "nonIBD"), c("CD", "nonIBD")) ; comp_type_2<- list(c("UC", "nonIBD"), c("CD", "nonIBD"))
comp_act <- list(c("Remission", "nonIBD"), c("Border", "nonIBD"), c("Active", "nonIBD"))
comp_prox <- list(c("Distal", "Undetected"), c("Proximal", "Undetected"), c("Detected", "Undetected")) ; comp_prox_one <- list(c("Distal", "Proximal"))
comp_exp <- list(c("Non", "Exp"))
comp_pp <- list(c("Prior", "Posterior"), c("Negative", "Posterior")) ; comp_pp_2 <- list(c("Negative", "Prior"), c("Negative", "Posterior")) ; comp_pp_3 <- list(c("Prior", "Posterior"), c("Negative", "Posterior"), c("Negative", "Prior"))


## Color set
# 1) two color
col2 <- c("#3399FF","#CC3333") # Blue vs. Red  
col2.1 <- c("gray81", "#FC4E07") # Gray vs. Red
col2.2 <- c("#CC3333","#3399FF") # Red vs. Blue
col2.3 <- c("#9966CC", "#009999") # purple vs. blue

# 2) three colors
col3 <- c("gray71", "#00AFBB", "#FC4E07")
col3.1 <- c("gray81", "#009E73", "#FC4E07") # gray for Fn.detection
col3.2 <- c("#009E73", "#FC4E07", "Darkred") # green + red + darkred 
col3.3 <- c("gray","#00AFBB", "#FC4E07") # gray + cyan + orange
col3.4 <- c("gray","Darkred", "#FC4E07") # gray + darkred + orange
col3.5 <- c("#009999", "#990099", "003399") # for cluster_3

# 3) four colors
col4 <- c("#00AFBB", "#009E73", "#E7B800", "#FC4E07") # cyan + green + yellow + red (nonIBD/Remission/Border/Active)
col4.1 <- c("gray81", "#009E73", "#E7B800", "#FC4E07") # for Fn.detection 
col4.2 <- c("#FFCF48", "#FC2847", "#93DFB8", "#158078")  # "Sunglow", "Scarlet", "Sea Green", "Pine Green"
col4.3 <- c("gray71", "gray81", "#00AFBB", "#FC4E07") # Two gray 4 colors

#### 2. Input ####
# Input
bugs <- read.csv("input/bugs_1526.csv", header = T, row.names = 1) # 533 microbes x 1526 samples
  NOD <- rowSums(bugs>0) # Number of detection for every microbes across whole samples
meta <- read.csv("input/meta_1526.csv", header = T, row.names = 1) # Metadata for 1526 samples
mark <- read.csv("input/lefse.marker.csv", header = T, row.names = 1) # LEfSE results; contains 14 non-IBD + 12 IBD marker species 
  ibd.mark <- rownames(mark)[mark$Class=="IBD"] ; non.mark <- rownames(mark)[mark$Class=="nonIBD"] ; mark.list <- c(ibd.mark, non.mark)

ph <- t(read.csv("input/phylum.csv", header = T, row.names = 1)) ; ph <- as.data.frame(ph[order(rownames(ph)),]) # Phylum level composition
core_CRC <- read.csv("input/core.species.CRC.csv", header=T) ; core_CRC <- as.character(core_CRC$core.CRC)  # Core signature of CRC
  
lef <- read.csv("input/lefse.marker.csv", header = T) ; lef_order <- lef$Microbes[order(lef$Class, lef$Log10_LDA, decreasing = T)] ; head(lef) # Import LEfSe results file
re <- read.csv("input/replicated.csv", header = T, row.names = 1) # replicated sample

# Deduplicated metadata
meta_2 <- meta[order(meta$Participant.ID, meta$day.from.first, decreasing = T),]
dedup.meta <- meta_2[!duplicated(meta_2$Participant.ID),]

  # Average age
  mean(dedup.meta$consent_age, na.rm = T) # 27.61905
  mean(dedup.meta$consent_age[dedup.meta$Condition=="nonIBD"]) # 29.65385
  mean(dedup.meta$consent_age[dedup.meta$Condition=="IBD"], na.rm = T) # 26.94937

  # Average days of collection
  mean(dedup.meta$day.from.first) # 296.217
  mean(dedup.meta$day.from.first[dedup.meta$Condition=="nonIBD"]) # 307.5385
  mean(dedup.meta$day.from.first[dedup.meta$Condition=="IBD"]) # 292.5375
  
  # Standard deviation days of collection
  sd(dedup.meta$day.from.first) # 39.61938
  sd(dedup.meta$day.from.first[dedup.meta$Condition=="nonIBD"]) # 38.406061
  sd(dedup.meta$day.from.first[dedup.meta$Condition=="IBD"]) # 39.52597
  

# Ordering factors 
meta$Condition <- factor(meta$Condition, levels=c("nonIBD", "IBD")) 
meta$Type <- factor(meta$Type, levels=c("nonIBD", "CD", "UC"))
mark$Class <- factor(mark$Class, levels = c("nonIBD", "IBD"))
meta$Activity <- factor(meta$Activity, levels=c("nonIBD", "Remission", "Border", "Active", "NA")) # Before NA removal
meta$Fn.experience <- factor(meta$Fn.experience, levels=c("Non", "Exp"))
meta$Fn.proximity <- factor(meta$Fn.proximity, levels=c("Undetected", "Distal", "Proximal", "Detected")) 

## Here, we adopted three different criteria for F. nucleatum prior/posterior classification (1-Original ; 2-Strict ; 3-Simple)
  ## 1) Original 
    ## All F. nucleatum detected samples = positive, 
    ## And, For the samples after the F. nucleatum detection, temporal distance (week) was calculated from the nearest F. nucleatum-detected samples
    meta$Fn.pre.post <- factor(meta$Fn.pre.post, levels=c("Negative", "Prior", "Positive", "Posterior")) # Prior --> Positive --> Posterior 

  ## 2) Strict
    ## Only the very first F. nucleatum-detected samples was considered as F. nucleatum-positive, the latter ones were "posterior"
    ## And, the temporal proximity (week) was calculated from the very first F. nucleatum-detected samples
    ## Usually, this is the common choice for analysis unless undescribed. 
    meta$Fn.pre.post.once <- factor(meta$Fn.pre.post.once, levels=c("Negative", "Prior", "Positive", "Posterior")) 
    
  ## 3) Simple 
    ## Almost same to strict version, but it considered the first F. nucleatum-detected samples as "posterior" --> so, no F. nucleatum-detected samnples (no positive)
    ## This is to increase the number of posterior samples, and simplify the analysis
    meta$Fn.pre.post.simple <- factor(meta$Fn.pre.post.simple, levels=c("Negative", "Prior", "Posterior")) 

# Sorting metadata which contain activity (SCCAI, HBI) information
meta.act <- meta[which(meta$Activity!="NA"),] ; meta.act$Activity <- droplevels(meta.act$Activity) # remove NA in meta; 1487 survived


# Matching sample order in bugs profiling data and metadata
tbugs <- as.data.frame(t(bugs)) # Upon transposition, it is converted into matrix
tbugs.order <- tbugs[order(rownames(tbugs)),]
meta.order <- meta[order(meta$External.ID, decreasing=F),]
bugs.order <- bugs[,order(colnames(bugs))]


# Histogram for QC filtering
## We excluded samples having less than 18 microbes per sample because too low number of species can distort the entire statistics a lot. 

#qc <- read.csv("input/QC_filtering.csv", header = T)

#tiff("figures/qc_filtering.tiff", res=200, width = 300, height = 300, units = "px") ; ggplot(qc, aes(species)) + geom_histogram(stat = "count", color="gray40", fill="gray40") + theme_classic() +
#  xlab("")+ylab("")+  theme(axis.text = element_blank(), axis.ticks = element_blank()) ; dev.off()


#### 3. PCoA & Heatmap ####
# 1) PCoA using Bray-Curtis distance
# Calculate distance between samples & principle coordinates
# Log transformation to reduce the effect by a few dominant microbes
tbugs.log <- log10(tbugs.order+1e-05)+5 # 1e-05 is the half of minimum microbial abundance detected; Add 5 to remove negative value 


# Bray-Curtis distance
d.log <- vegdist(tbugs.log, method = "bray")

# PCoA
pcoa.log <- cmdscale(d.log, 2, eig = T) ; identical(rownames(pcoa.log$points), meta.order$External.ID) # TRUE

# Merging with ordered metadata
pcoa.meta.log <- data.frame(pcoa.log$points, meta.order) ; colnames(pcoa.meta.log)[1:2] <- c("PC1", "PC2")
  pcoa.meta.log$Sex <- factor(pcoa.meta.log$Sex, levels=c("Female", "Male"))
  pcoa.meta.log$Firmicutes <- ph$Firmicutes ; pcoa.meta.log$Bacteroidetes <- ph$Bacteroidetes
  pcoa.meta.act <- pcoa.meta.log[!is.na(pcoa.meta.log$Activity),] ; pcoa.meta.act$Activity <- droplevels(pcoa.meta.act$Activity) # Activity sorting
  
  ## Two eigen values: the first 2 principle coordinates
  pcoa.log$eig[1:2]
  
  ## Calculate variance explained by the two coordinates
  sum((pcoa.log$eig[1:2])^2)/sum((pcoa.log$eig)^2) 
  x.var <- round((sum((pcoa.log$eig[1])^2)/sum((pcoa.log$eig)^2)*100), 1)
  y.var <- round((sum((pcoa.log$eig[2])^2)/sum((pcoa.log$eig)^2)*100), 1)
  
#-----------------PCOA plot-----------------------#
pcoa.plot <- function(data, class) {
  attach(data)
  ggplot(data, aes(PC1, PC2)) + geom_point(aes(color=class), alpha=0.85) + 
    xlab(paste0("PC1 (", x.var, "%)")) + ylab(paste0("PC2 (", y.var, "%)")) + 
    theme_few() + 
    theme(axis.text = element_blank(), 
          axis.title.x = element_text(vjust=-1.5, size=rel(1.2)), 
          axis.title.y = element_text(vjust=1.7, size=rel(1.2)), 
          legend.text = element_text(size=rel(0.8)), 
          legend.title = element_blank(), 
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.background = element_rect(color = "black", size=0.5),
          legend.key = element_rect(fill = "white"),
          legend.key.width = unit(0.8, "line"),
          legend.key.height = unit(0.8,"line"), 
          legend.spacing.y = unit(0, "mm"), # Spacing between legend title and items
          legend.spacing.x = unit(0.6, "mm"))
}

#### 3-1. Draw PCoA ####
# Subject information
tiff(file = "figures/PCOA.site.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, Site) ; dev.off()
tiff(file = "figures/PCOA.sex.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, Sex) +  theme(legend.text = element_text(size=rel(0.9)), legend.key.width=unit(0.8, "line")); dev.off()
tiff(file = "figures/PCOA.patient.tiff", width = 800, height = 800, units = "px", res = 200) ;  pcoa.plot(pcoa.meta.log, Participant.ID) + theme(legend.position = "none") ; dev.off()
tiff(file = "figures/PCOA.condition.tiff", width = 800, height = 800, units = "px", res = 200) ;  pcoa.plot(pcoa.meta.log, Condition) + 
  scale_fill_manual(values = col2) + scale_color_manual(values = col2) + stat_ellipse(type="t", aes(group=Condition, color=Condition), alpha=0.25, size=1.5) +
  theme(legend.text = element_text(size=rel(0.9)), legend.key.width=unit(0.8, "line")) ; dev.off()
tiff(file = "figures/PCOA_IBDtype.tiff", width = 800, height = 800, units = "px", res = 200) ;  pcoa.plot(pcoa.meta.log, Type) + 
  scale_fill_manual(values = col3.2) + scale_color_manual(values = col3.2) +  stat_ellipse(type="t", aes(group=Type, color=Type), alpha=0.25, size=1.5) +
  theme(legend.text = element_text(size=rel(0.9)), legend.key.width=unit(0.8, "line")) ; dev.off()
tiff(file = "figures/PCOA.Activity.tiff", width = 800, height = 800, units = "px", res = 200) ;  pcoa.plot(pcoa.meta.act, Activity) + scale_fill_manual(values = col4.1) + scale_color_manual(values = col4.1); dev.off()

# F. nucleatum-dynamics
tiff(file = "figures/PCOA.Fn.exp.tiff", width = 800, height = 800, units = "px", res = 200);  pcoa.plot(pcoa.meta.log, Fn.experience) + scale_fill_manual(values = col2.1) + scale_color_manual(values = col2.1) + theme(legend.text = element_text(size=rel(0.9)), legend.key.width=unit(0.8, "line")) ; dev.off()
tiff(file = "figures/PCOA.Fn.proximity.tiff", width = 800, height = 800, units = "px", res = 200 );  pcoa.plot(pcoa.meta.log, Fn.proximity) + scale_fill_manual(values = col4.1) + scale_color_manual(values = col4.1); dev.off()
tiff(file = "figures/PCOA.Fn.pre_post.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, Fn.pre.post) + scale_fill_manual(values = col4.2) + scale_color_manual(values = col4.1); dev.off()
tiff(file = "figures/PCOA.Fn.pre_post_strict.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, Fn.pre.post.once) + scale_fill_manual(values = col4.2) + scale_color_manual(values = col4.1); dev.off()
tiff(file = "figures/PCOA.Fn.pre_post_simple.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, Fn.pre.post.simple) + scale_fill_manual(values = col3.1) + scale_color_manual(values = col3.1); dev.off()

# Phlyum level abundance 
tiff(file = "figures/PCOA.Firmicutes.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, Firmicutes) + 
  scale_color_continuous(low = "black", high = "yellow", breaks=c(25, 50, 75), label=c(25, 50, 75)); dev.off()
tiff(file = "figures/PCOA.Bacteroidetes.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, Bacteroidetes) + 
  scale_color_continuous(low = "black", high = "yellow", breaks=c(25, 50, 75), label=c(25, 50, 75)); dev.off()

# To get P-value at PCoA using Principal coordinates (PC1 & PC2)
  # Using logistic regression models,
  # We will test PC1 & PC2 can successfully discriminate two groups. 

eig_test <- data.frame(External.ID=rownames(pcoa.log$points), pcoa.log$points, eig=pcoa.log$eig, # eigen value for each sample
                       Participant.ID=meta.order$Participant.ID, Type=meta.order$Type, Condition=meta.order$Condition, Sex=as.factor(meta.order$Sex),
                       Site=as.factor(meta.order$Site), Fn.experience=meta.order$Fn.experience, 
                       Fn.pre.post.simple=meta.order$Fn.pre.post.simple) ; rownames(eig_test) <- NULL ; colnames(eig_test)[2:3] <- c("PC1", "PC2"); head(eig_test)

# For pairwise comparison between IBD subtypes
eig_UC.CD <- eig_test[eig_test$Type!="nonIBD",] ; eig_UC.CD$Type <- droplevels(eig_UC.CD$Type) 
eig_UC.NON <- eig_test[eig_test$Type!="CD",] ; eig_UC.NON$Type <- droplevels(eig_UC.NON$Type) 
eig_CD.NON <- eig_test[eig_test$Type!="UC",] ; eig_CD.NON$Type <- droplevels(eig_CD.NON$Type)
eig_Fn_pre_post <- eig_test[eig_test$Fn.pre.post.simple!="Negative",] ; eig_Fn_pre_post$Type <- droplevels(eig_Fn_pre_post$Type) 


# Null-hypothesis: IBD and nonIBD are indistinguishable based on Principal coordinates
fit_IBD <- glm(formula = Condition~PC1+PC2, data = eig_test, family = "binomial")
summary(fit_IBD) # PC1: <2e-16 ; PC2: 2.1e-13
anova(fit_IBD, test = "Chisq") # PC1 <2e-16 ; PC2 5.506e-14

# Null-hypothesis: UC and CD are are indistinguishable based on Principal coordinates
  # Regression model to predict 
fit_UC.CD <- glm(formula = Type~PC1+PC2, data = eig_UC.CD, family = "binomial")
summary(fit_UC.CD) # PC1: 0.1495 ; PC2 0.0991
anova(fit_UC.CD, test = "Chisq") # PC1 0.1726 ; PC2 0.0988

# Null-hypothesis: F.n.-innocent subjects and -experienced subjects are indistinguishable based on Principal coordinates
fit_Fn_exp <- glm(formula = Fn.experience~PC1+PC2, data = eig_test, family = "binomial")
summary(fit_Fn_exp) # PC1: 5.45e-14 ; PC2 3.96e-07
anova(fit_Fn_exp, test = "Chisq") # PC1 2.694e-15 ; PC2 2.550e-07

# Null-hypothesis: F.n.-prior samples and -posterior samples are indistinguishable based on Principal coordinates
fit_Fn_pre_post <- glm(formula = Fn.pre.post.simple~PC1+PC2, data = eig_Fn_pre_post, family = "binomial")
summary(fit_Fn_pre_post) # PC1: 5.45e-14 ; PC2 3.96e-07
anova(fit_Fn_pre_post, test = "Chisq") # PC1 2.694e-15 ; PC2 2.550e-07


# Sorting F. nucleatum-experienced subjects
pcoa.meta.log <- pcoa.meta.log[order(pcoa.meta.log$Participant.ID, pcoa.meta.log$visit.num, decreasing = F),]; head(pcoa.meta.log)
pcoa.exp <- pcoa.meta.log[!is.na(pcoa.meta.log$Fn.directionality),] ; dim(pcoa.exp) # 317 x 42 
pcoa.exp_20 <- pcoa.exp[abs(pcoa.exp$Fn.directionality.strict)<=20,]
fn_pat <- unique(pcoa.exp$Participant.ID) ; length(fn_pat) # number of F.nucleatum-experienced subjects; 20

tiff("figures/PCoA.subject.path.tiff", res = 200, width = 800, height = 3200, units = "px"); ggplot(pcoa.exp_20, aes(PC1, PC2)) + 
  geom_path(aes(color=Fn.directionality.strict), arrow = NULL, lineend = "round" , size=1.2)+
  xlab(paste0("PC1 (", x.var, "%)")) + ylab(paste0("PC2 (", y.var, "%)")) + 
  facet_grid(Participant.ID~Condition) + 
  geom_hline(yintercept = 0, linetype="dotted", color="gray40", size=1, alpha=0.5)+
  geom_vline(xintercept = 0, linetype="dotted", color="gray40", size=1, alpha=0.5)+
  scale_color_gradient2(low = "green", mid = "black", midpoint = 0, high = "red") + 
  theme_few() + 
  theme(axis.text = element_text(size=rel(0.7)), 
        axis.title.x = element_text(vjust=-1.5, size=rel(1.2)), 
        axis.title.y = element_text(vjust=1.7, size=rel(1.2)),
        strip.text.x = element_text(size=rel(2.3)),
        strip.text.y = element_text(size=rel(1.1)),
        legend.text = element_text(size=rel(0.8)), 
        legend.title = element_blank(), 
        legend.position = "right",
        legend.justification = c(1, 0),
        legend.background = element_rect(color = "black", size=0.5),
        legend.key = element_rect(fill = "white"),
        legend.key.width = unit(0.8, "line"),
        legend.key.height = unit(0.8,"line"), 
        legend.spacing.y = unit(0, "mm"), # Spacing between legend title and items
        legend.spacing.x = unit(0.6, "mm")) ; dev.off()


#### 3-2. Heatmap ####

# Log transformation 
# If not transformed, most of tile (bacterial abundance per sample) looked dark in heatmap.
# Add pseudo-abundance (0.00001) to avoid log0 ; 0.00001 is a half value of the minimum detected bacterial abundance.
bugs2 <- log10(bugs+0.00001) ; bugs2$microbe <- rownames(bugs2) 

# Ordering heatmap axis (Y-axis=Bugs ; X-axis=Samples)
order1 <- rownames(bugs2)[order(rowMeans(log10(bugs+0.00001)), decreasing = F)] # Y-axis ; Decreasing order of bugs mean abundnace
order2 <- meta$External.ID[order(meta$Type, meta$Participant.ID, meta$Week, decreasing = F)] # X-axis ; samples
  rownames(meta) <- meta$External.ID ; meta_order2 <- meta[order2,] 
  meta_order2$Type # nonIBD --> CD --> UC
  meta_order2$Participant.ID; meta_order2$Participant.ID <- factor(meta_order2$Participant.ID, levels=unique(meta_order2$Participant.ID))
  meta_order2$Week # Increasing order ; 


# Melting
all_melt <- melt(bugs2, id.vars = "microbe") 
  ## Add condition information for molten data
  all_melt$Condition <- "nonIBD"
  all_melt$Condition[all_melt$variable %in% meta$External.ID[meta$Type=="CD"]] <- "CD"
  all_melt$Condition[all_melt$variable %in% meta$External.ID[meta$Type=="UC"]] <- "UC"
  all_melt$Condition <- factor(all_melt$Condition, levels=c("CD", "UC", "nonIBD")) ; head(all_melt)


# Draw heatmap
tiff("figures/heatmap_log.tiff", res = 200, width = 2400, height = 800, units = "px") ; ggplot(all_melt) + 
  geom_tile(aes(variable, microbe, fill=value), color="transparent") +  xlab("") + ylab("") + 
  scale_fill_gradient(low = "black", high = "yellow") + 
  scale_y_discrete(limit=order1) +  scale_x_discrete(limit=order3) +
  theme_few() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text.align = 0.5,
        legend.background = element_rect(color="black", size=0.5),
        legend.text = element_text(size=rel(0.6)),
        legend.justification = c(0, 0), # adjust y-coordinate of legend
        legend.key.height = unit(0.8, units = "line"),
        legend.key.width = unit(0.4, units = "line"),
        legend.spacing.x = unit(1.5, units = "mm")); dev.off() # spcaing x = distance between legend bar and text 

# Add clustering information over heatmap
set.seed(12345)
bugs_mat <- as.matrix(log10(bugs+0.00001)) ; bugs_mat <- bugs_mat[,order2] 
  # kmean clustering of samples using log transformed bugs abundance matrix (3 clusters)
  km <- kmeans(x = t(bugs_mat), centers = 3, iter.max = 1000)
    identical(meta_order2$External.ID, names(km$cluster)) # TRUE
    meta_order2$cluster_3 <- factor(km$cluster, levels = c("1", "2", "3"))

# Draw bar graph for samples  
    ## Ordering for bar samples
    order3 <- meta$External.ID[order(meta$Type, meta$Participant.ID, meta$Week, decreasing = T)] # Reversed order2 

tiff("figures/tile_Type.tiff", res = 300, width = 3200, height = 80, units = "px"); ggplot(meta_order2, aes(x=External.ID, y=1, fill=Type)) + geom_tile() +
  scale_x_discrete(limits=order3) +  scale_fill_manual(values=col3.2) + theme_void() + theme(legend.position = "none") ; dev.off()

tiff("figures/tile_cluster.tiff", res = 300, width = 3200, height = 80, units = "px"); ggplot(meta_order2, aes(x=External.ID, y=1, fill=cluster_3)) + geom_tile() +
  scale_x_discrete(limits=order3) +  scale_fill_manual(values=col3.5) + theme_void() + theme(legend.position = "none") ; dev.off()

tiff("figures/tile_subject.tiff", res = 300, width = 3200, height = 80, units = "px"); ggplot(meta_order2) + 
  geom_tile(aes(x=External.ID, y=1, fill=Participant.ID), color="transparent") + 
  #geom_rect(color="black", aes(xmin=1.5, xmax=1526, ymin=0.5, ymax=1.5), size=0.15, fill="transparent")+
  scale_x_discrete(limits=order3) + scale_fill_manual(values=rep(c("gray32", "gray70"), 53)) + theme_void() + theme(legend.position = "none") ; dev.off()

# Clustering significance
table(meta_order2$cluster_3) # 1=665, 2= 318, 3=543
table(meta_order2$Type) # nonIBD= 407, CD= 702, UC=417
tab <- xtabs(~meta_order2$Type + meta_order2$cluster_3)
print(tab)
 # nonIBD: 1= 102, 2=56, 3=249
 # CD: 1=309, 2=194, 3=199
 # UC: 1= 254, 2=68, 3=95

# nonIBD enrichment in each cluster
  # 1) Cluster 1
  cl <- c("1", "1", "non1", "non1"); con <- c("nonIBD", "other", "nonIBD", "other"); count <- c(102, 563, 305, 556) ; non.dat <- data.frame(cl, con, count)
  non.tab <- xtabs(count~cl+con, data=non.dat) ; fisher.test(non.tab, alternative = "two.sided") # p-value < 2.2e-16 ; odd ratio 0.3305022
  
  # 2) Cluster 2
  cl <- c("2", "2", "non2", "non2") ; con <- c("nonIBD", "other", "nonIBD", "other") ; count <- c(56, 262, 351, 857) ; non.dat <- data.frame(cl, con, count)
  non.tab <- xtabs(count~cl+con, data=non.dat) ; fisher.test(non.tab, alternative = "two.sided") # p-value < 2.2e-16 ; odd ratio 0.5220721
  
  # 3) Cluster 3
  cl <- c("3", "3", "non3", "non3") ; con <- c("nonIBD", "other", "nonIBD", "other") ; count <- c(249, 294, 158, 825) ; non.dat <- data.frame(cl, con, count)
  non.tab <- xtabs(count~cl+con, data=non.dat) ; fisher.test(non.tab, alternative = "two.sided") # p-value < 2.2e-16 ; odd ratio 4.417385

# UC enrichment in each cluster
  # 1) Cluster 1
  cl <- c("1", "1", "non1", "non1"); con <- c("UC", "other", "UC", "other"); count <- c(254, 411, 163, 698) ; UC.dat <- data.frame(cl, con, count)
  UC.tab <- xtabs(count~cl+con, data=UC.dat) ; fisher.test(UC.tab, alternative = "two.sided") # p-value < 2.2e-16 ; odd ratio 0.3781199
  
  # 2) Cluster 2
  cl <- c("2", "2", "non2", "non2") ; con <- c("UC", "other", "UC", "other") ; count <- c(68, 250, 349, 859) ; UC.dat <- data.frame(cl, con, count)
  UC.tab <- xtabs(count~cl+con, data=UC.dat) ; fisher.test(UC.tab, alternative = "two.sided") # p-value = 0.007209 ; odd ratio 1.493327
  
  # 3) Cluster 3
  cl <- c("3", "3", "non3", "non3") ; con <- c("UC", "other", "UC", "other") ; count <- c(95, 448, 322, 661) ; UC.dat <- data.frame(cl, con, count)
  UC.tab <- xtabs(count~cl+con, data=UC.dat) ; fisher.test(UC.tab, alternative = "two.sided") # p-value =7.346e-11 ; odd ratio 2.296072

# CD enrichment in each cluster
  # nonIBD: 1= 102, 2=56, 3=249
  # CD: 1=309, 2=194, 3=199
  # UC: 1= 254, 2=68, 3=95
  
  # 1) Cluster 1
  cl <- c("1", "1", "non1", "non1"); con <- c("CD", "other", "CD", "other"); count <- c(309, 356, 393, 468) ; CD.dat <- data.frame(cl, con, count)
  CD.tab <- xtabs(count~cl+con, data=CD.dat) ; fisher.test(CD.tab, alternative = "two.sided") # p-value = 0.7562 ; odd ratio 1.033617
  
  # 2) Cluster 2
  cl <- c("2", "2", "non2", "non2") ; con <- c("CD", "other", "CD", "other") ; count <- c(194, 124, 508, 700) ; CD.dat <- data.frame(cl, con, count)
  CD.tab <- xtabs(count~cl+con, data=CD.dat) ; fisher.test(CD.tab, alternative = "two.sided") # p-value = 2.243e-09 ; odd ratio 2.154743
  
  # 3) Cluster 3
  cl <- c("3", "3", "non3", "non3") ; con <- c("CD", "other", "CD", "other") ; count <- c(199, 344, 503, 480) ; CD.dat <- data.frame(cl, con, count)
  CD.tab <- xtabs(count~cl+con, data=CD.dat) ; fisher.test(CD.tab, alternative = "two.sided") # p-value = 5.662e-08 ; odd ratio 0.5522441

# Compare cluster enrichment between UC and CD
  # nonIBD are enriched in a particular cluster very significantly (in this trial cluster 3),
  # However, this enrichment can affect p-value of UC and CD in fisher's exact test reducing the count of cluster 3 in other condition
  # So, we tried to compare the ratio of cluster 1 and 2 in UC and CD (2 clusters that are not enriched in nonIBD; somehow reflects dysbiotic microbiome) 

  cl <- c("1", "1", "all", "all") ; con <- c("UC", "CD", "UC", "CD") ; count <- c(254, 309, 417, 702) ; compare.data <- data.frame(cl, con, count)
  compare.mat <- xtabs(count~cl+con, data = compare.data) ; fisher.test(compare.mat, alternative = "two.sided") # p-value = 0.002201, odd ratio 0.7228025
  
  cl <- c("2", "2", "all", "all") ; con <- c("UC", "CD", "UC", "CD") ; count <- c(68, 194, 417, 702) ; compare.data <- data.frame(cl, con, count)
  compare.mat <- xtabs(count~cl+con, data = compare.data) ; fisher.test(compare.mat, alternative = "two.sided") # p-value = 0.0005427, odd ratio 1.694077
  
  
  
#### 4. Alpha Diversity  ####

# Diversity function
div.plot  <- function(data, x, y, color, comparisons, label.y) {
  attach(data)
  ggplot(data, aes(x = x, y=y, color=x, fill=x)) + xlab("") + ylab("")+
    #geom_sina(aes(color=x), alpha=0.2, size=0.8)+
    geom_violin(alpha=0.6, draw_quantiles = 0.5, alpha=0.4) +
    #geom_boxplot(alpha=0.5, outlier.alpha = 0.1) +
    scale_fill_manual(values = color) + scale_color_manual(values = color)+
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label="p.format", tip.length = 0, label.y = label.y) +
    theme_classic() + 
    theme(axis.text.x=element_text(size=rel(1.15), color="black"), 
          axis.text.y = element_text(size=rel(1.15), color="black"), 
          axis.ticks.x = element_blank(),
          legend.position = "none")
}

# 4-1) Sex
# col2.2 = Red vs. Blue
tiff(file = "figures/Diversity_Shannon_sex.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Sex, Shannon, col2.3, comp_sex, 3.7) ; dev.off() 
tiff(file = "figures/Diversity_Evenness_sex.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Sex, Evenness, col2.3, comp_sex, 0.85) ; dev.off() 
tiff(file = "figures/Diversity_Richness_sex.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Sex, Richness, col2.3, comp_sex, 110) ; dev.off() 

# 4-2) Condition
tiff(file = "figures/Diversity_Shannon_Condition.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Condition, Shannon, col2, comp_con, 3.7) ; dev.off() 
tiff(file = "figures/Diversity_Evenness_Condition.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Condition, Evenness, col2, comp_con, 0.85) ; dev.off() 
tiff(file = "figures/Diversity_Richness_Condition.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Condition, Richness, col2, comp_con, 110) ; dev.off() 

# 4-3) Type 
tiff(file = "figures/Diversity_Shannon_Type.tiff", width = 450, height = 700, units = "px", res = 200) ; div.plot(meta, Type, Shannon, col3.2, comp_type, c(3.6, 3.8, 4)) ; dev.off() 
tiff(file = "figures/Diversity_Evenness_Type.tiff", width = 450, height = 700, units = "px", res = 200) ; div.plot(meta, Type, Evenness, col3.2, comp_type, c(0.85, 0.9, 0.95)) ; dev.off() 
tiff(file = "figures/Diversity_Richness_Type.tiff", width = 450, height = 700, units = "px", res = 200) ; div.plot(meta, Type, Richness, col3.2, comp_type, c(106, 112, 118)) ; dev.off() 

# 4-4) Diseases severity (Activity: 4 categories)
tiff(file = "figures/Diversity_Shannon_Activity.tiff", width = 600, height = 700, units = "px", res = 200) ; div.plot(meta.act, Activity, Shannon, col4, comp_act, c(3.6, 3.8, 4)) ; dev.off() 
tiff(file = "figures/Diversity_Evenness_Activity.tiff", width = 650, height = 700, units = "px", res = 200) ; div.plot(meta.act, Activity, Evenness, col4, comp_act, c(0.85, 0.9, 0.95)) ; dev.off() 
tiff(file = "figures/Diversity_Richness_Activity.tiff", width = 650, height = 700, units = "px", res = 200) ; div.plot(meta.act, Activity, Richness, col4, comp_act, c(106, 112, 118)) ; dev.off() 

# 4-5) F.nucleatum-experience
# F. nucleatum "experience" means that the subjects have experienced F. nucleatum at least one time during long-term fecal surveilence 
tiff(file = "figures/Diversity_Shannon_Fn_experience.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.experience, Shannon, col2.1, comp_exp, 3.7) ; dev.off() 
tiff(file = "figures/Diversity_Evenness_Fn_experience.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.experience, Evenness, col2.1, comp_exp, 0.85) ; dev.off() 
tiff(file = "figures/Diversity_Richness_Fn_experience.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.experience, Richness, col2.1, comp_exp, 110) ; dev.off() 

# 4-6) F.nucleatum-proximity
tiff(file = "figures/Diversity_Shannon_Fn_proximity.tiff", width = 600, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.proximity, Shannon, col4, comp_prox, c(3.6, 3.8, 4)) ; dev.off() 
tiff(file = "figures/Diversity_Evenness_Fn_proximity.tiff", width = 650, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.proximity, Evenness, col4, comp_prox, c(0.85, 0.9, 0.95)) ; dev.off() 
tiff(file = "figures/Diversity_Richness_Fn_proximity.tiff", width = 650, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.proximity, Richness, col4, comp_prox, c(106, 112, 118)) ; dev.off() 

# 4-7) F.nucleatum-prior/posterior
tiff(file = "figures/Diversity_Shannon_Fn_pre-post.tiff", width = 450, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.pre.post.simple, Shannon, col3.1, comp_pp_3, c(3.6, 3.8, 4)) ; dev.off() 
tiff(file = "figures/Diversity_Evenness_Fn_pre-post.tiff", width = 450, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.pre.post.simple, Evenness, col3.1, comp_pp_3, c(0.85, 0.9, 0.95)) ; dev.off() 
tiff(file = "figures/Diversity_Richness_Fn_pre-post.tiff", width = 450, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.pre.post.simple, Richness, col3.1, comp_pp_3, c(106, 112, 118)) ; dev.off() 


# All-in-one graph using melting
  # 1) For Condition + Type + Activity
div.melt <- melt(data = meta, id.vars = c("Shannon", "Evenness", "Richness", "Human.Read.Fraction"), 
                 measure.vars = c("Condition", "Type", "Activity"), na.rm = T) 
div.melt <- div.melt[which(div.melt$variable!="Type" | div.melt$value!="nonIBD"),]
div.melt <- div.melt[which(div.melt$variable!="Activity" | div.melt$value!="nonIBD"),]
div.melt$value <- factor(div.melt$value, levels=c("nonIBD", "IBD", "CD", "UC", "Remission", "Border", "Active"))

col7 <- c("#3399FF", "#CC3333", "#FF6633", "#FF6666", "#FFCC99", "#FF9933", "#CC3333")
comp_all_1 <- list(c("nonIBD", "IBD"), c("nonIBD", "CD"), c("nonIBD", "UC"), c("Remission", "Active"), c("Border", "Active"))

melt_draw <- function(data, x, y, color, comparison, label.y) {
  attach(data)
  ggplot(div.melt, aes(x = x, y = y, color=x, fill=x)) + 
    geom_violin(alpha=0.4, draw_quantiles = 0.5, na.rm = T, scale = "area", width=0.7) + 
    scale_fill_manual(values = color) +  scale_color_manual(values = color) +
    stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", tip.length = 0, bracket.size = 0.5, label.y = label.y, size=rel(4))+ 
    xlab("") + ylab("") + theme_pubr() +
    theme(legend.position = "none",
          axis.text.x = element_text(size=rel(1.4)),axis.text.y = element_text(size=rel(1.1)))
}

tiff("figures/diversity_Shannon_all.tiff", res = 200, width = 1600, height = 900, units = "px") ; melt_draw(div.melt, value, Shannon, col7, comp_all_1, 
                                                                                                            label.y = c(3.6, 3.8, 4, 3.79, 3.6)) ; dev.off()
tiff("figures/diversity_Evenness_all.tiff", res = 200, width = 1600, height = 900, units = "px") ; melt_draw(div.melt, value, Evenness, col7, comp_all_1, 
                                                                                                             label.y = c(0.83, 0.88, 0.93, 0.88, 0.83)) ; dev.off()
tiff("figures/diversity_Richness_all.tiff", res = 200, width = 1600, height = 900, units = "px") ; melt_draw(div.melt, value, Richness, col7, comp_all_1, 
                                                                                                             label.y = c(103, 108, 113, 108, 103)) ; dev.off()
tiff("figures/Human.read.fraction_all.tiff", res = 200, width = 1600, height = 900, units = "px") ; melt_draw(div.melt, value, log10(Human.Read.Fraction), col7, comp_all_1, 
                                                                                                              label.y = c(-0.05, 0.2, 0.45, 0.2, -0.05)) ; dev.off()

# 2) For F.nucleatum experience + Proximity + Prior/Posterior
div.melt <- melt(data = meta, id.vars = c("Shannon", "Evenness", "Richness", "Human.Read.Fraction"), 
                 measure.vars = c("Fn.experience", "Fn.proximity", "Fn.pre.post.simple"), na.rm = T) 
div.melt <- div.melt[which(div.melt$value!="Undetected"),] ; div.melt <- div.melt[which(div.melt$value!="Negative"),] ; div.melt <- div.melt[which(div.melt$value!="Detected"),]
div.melt$value <- factor(div.melt$value, levels=c("Non", "Exp", "Distal", "Proximal", "Prior", "Posterior"))


col6 <- c("gray60", "#FC4E07", "#009E73", "darkred", "#00AFBB", "#FC4E07")
comp_all_2 <- list(c("Non", "Exp"), c("Distal", "Proximal"), c("Non", "Proximal"), c("Non", "Posterior"), c("Prior", "Posterior"))

tiff("figures/diversity_Shannon_all_Fn_related.tiff", res = 200, width = 1400, height = 900, units = "px") ; melt_draw(div.melt, value, Shannon, col6, comp_all_2, 
                                                                                                                       label.y = c(3.6, 3.6, 3.8, 4, 3.6)) ; dev.off()
tiff("figures/diversity_Evenness_all_Fn_related.tiff", res = 200, width = 1400, height = 900, units = "px") ; melt_draw(div.melt, value, Evenness, col6, comp_all_2, 
                                                                                                                     label.y = c(0.84, 0.84, 0.88, 0.92, 0.84)) ; dev.off()
tiff("figures/diversity_Richness_all_Fn_related.tiff", res = 200, width = 1400, height = 900, units = "px") ; melt_draw(div.melt, value, Richness, col6, comp_all_2, 
                                                                                                                        label.y = c(103, 103, 108, 113, 103)) ; dev.off()
tiff("figures/Human.read.fraction_all_Fn_related.tiff", res = 200, width = 1400, height = 900, units = "px") ; melt_draw(div.melt, value, log10(Human.Read.Fraction), col6, comp_all_2, 
                                                                                                                        label.y = c(-0.05, -0.05, 0.2, 0.45, -0.05)) ; dev.off()



  # 4-8) F.nucleatum distirubtion over time
# X-axis: Temporal proximity to the observation point of F. nucleatum (week)
# Y-axis: Longitudianl Shannon diversity of F. nucleatum-experienced subjects
# horizontal line indicates median Shannon diversity of F. nucleatum-inexpereinced subjects (Blue = nonIBD ; Red = IBD)

tiff(file = "figures/Diversity_Shannon_Scatter_week.tiff", width = 1000, height = 800, units = "px", res = 200) ; ggplot(meta, aes(x = Fn.directionality.strict, y=Shannon, group=Condition, color=Condition)) + 
  geom_smooth(method = "loess", span=0.8, alpha=0.2, se = T) + xlab("") + ylab("") +  scale_color_manual(values = col2)+
  geom_hline(yintercept = c(median(meta$Shannon[meta$Fn.experience=="Non" & meta$Condition=="nonIBD"]), median(meta$Shannon[meta$Fn.experience=="Non" & meta$Condition=="IBD"])), 
             alpha=0.3, size=1, color=col2, linetype="dotted")+
  geom_vline(xintercept = 0, color="gray", size=1.2, alpha=0.4)+
  theme_classic() + 
  theme(legend.position = c(1,0.85),
        legend.text = element_text(size=rel(1.4)), 
        legend.title = element_blank(),
        legend.justification = c(1,0),
        legend.background = element_rect(color = "black", size=0.5),
        legend.key = element_rect(fill = "white", color="transparent"),
        legend.key.width = unit(0.8, "line"),
        legend.key.height = unit(0.8,"line"), 
        legend.spacing.y = unit(0, "mm"), 
        legend.spacing.x = unit(1, "mm"),
        axis.text.x = element_text(size=rel(1.2), color="black"),
        axis.text.y = element_text(size=rel(1.2), color="black"), 
        axis.title = element_text(size=rel(1.6))) ; dev.off() 

# evenness
tiff(file = "figures/Diversity_Evenness_Scatter_week.tiff", width = 1000, height = 800, units = "px", res = 200) ; ggplot(meta, aes(x = Fn.directionality.strict, y=Evenness, group=Condition, color=Condition)) + 
  geom_smooth(method = "loess", span=0.8, alpha=0.2, se = T) + xlab("") + ylab("") +  scale_color_manual(values = col2)+
  geom_hline(yintercept = c(median(meta$Evenness[meta$Fn.experience=="Non" & meta$Condition=="nonIBD"]), median(meta$Evenness[meta$Fn.experience=="Non" & meta$Condition=="IBD"])), 
             alpha=0.3, size=1, color=col2, linetype="dotted")+
  geom_vline(xintercept = 0, color="gray", size=1.2, alpha=0.4)+
  theme_classic() + 
  theme(legend.position = c(1,0.85),
        legend.text = element_text(size=rel(1.4)), 
        legend.title = element_blank(),
        legend.justification = c(1,0),
        legend.background = element_rect(color = "black", size=0.5),
        legend.key = element_rect(fill = "white", color="transparent"),
        legend.key.width = unit(0.8, "line"),
        legend.key.height = unit(0.8,"line"), 
        legend.spacing.y = unit(0, "mm"), 
        legend.spacing.x = unit(1, "mm"),
        axis.text.x = element_text(size=rel(1.2), color="black"),
        axis.text.y = element_text(size=rel(1.2), color="black"), 
        axis.title = element_text(size=rel(1.6))) ; dev.off() 

# richness
tiff(file = "figures/Diversity_Richness_Scatter_week.tiff", width = 1000, height = 800, units = "px", res = 200) ; ggplot(meta, aes(x = Fn.directionality.strict, y=Richness, group=Condition, color=Condition)) + 
  geom_smooth(method = "loess", span=0.8, alpha=0.2, se = T) + xlab("") + ylab("") +  scale_color_manual(values = col2)+
  geom_hline(yintercept = c(median(meta$Richness[meta$Fn.experience=="Non" & meta$Condition=="nonIBD"]), median(meta$Richness[meta$Fn.experience=="Non" & meta$Condition=="IBD"])), 
             alpha=0.3, size=1, color=col2, linetype="dotted")+
  geom_vline(xintercept = 0, color="gray", size=1.2, alpha=0.4)+
  theme_classic() + 
  theme(legend.position = c(1,0.85),
        legend.text = element_text(size=rel(1.4)), 
        legend.title = element_blank(),
        legend.justification = c(1,0),
        legend.background = element_rect(color = "black", size=0.5),
        legend.key = element_rect(fill = "white", color="transparent"),
        legend.key.width = unit(0.8, "line"),
        legend.key.height = unit(0.8,"line"), 
        legend.spacing.y = unit(0, "mm"), 
        legend.spacing.x = unit(1, "mm"),
        axis.text.x = element_text(size=rel(1.2), color="black"),
        axis.text.y = element_text(size=rel(1.2), color="black"), 
        axis.title = element_text(size=rel(1.6))) ; dev.off() 



#### 4-1. Human read fraction ####

# Human read fraction is the ratio of read count mapped to human genome over total read count
# This value indicates extent of intestinal harsh situation. 
# Indeed, the human read fraction positively correlated with\severity scores of UC and CD

# SCCAI (Ulcerative colitis) ----
tiff(file = "figures/Human_Read_Fraction_SCCAI.tiff", width = 1000, height = 500, units = "px", res = 200) ; ggplot(meta, aes(SCCAI, log10(Human.Read.Fraction))) + 
  xlab("") + ylab("") + geom_point(alpha=0.05) + geom_smooth(method = "loess", span=0.9, color="gray30", alpha=0.4, se = T) + 
  stat_cor(label.x = 0, label.y = 0.6)+ theme_classic() + theme(axis.title = element_text(size=rel(1.2), vjust=5)) ; dev.off()

# HBI (Crohn's disease) ----
tiff(file = "figures/Human_Read_Fraction_HBI.tiff", width = 1000, height = 500, units = "px", res = 200) ; ggplot(meta, aes(HBI, log10(Human.Read.Fraction))) + 
  xlab("") + ylab("") + geom_point(alpha=0.05) + geom_smooth(method = "loess", span=0.9, color="gray30", alpha=0.4, se = T) + 
  stat_cor(label.x = 0, label.y = 0.6)+ theme_classic() + theme(axis.title = element_text(size=rel(1.2), vjust=5)) ; dev.off()


# Draw figures (using diversity plot funciton)----
tiff(file = "figures/Human_Read_Fraction_sex.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Sex, log10(Human.Read.Fraction), col2.3, comp_sex, 0.1) ; dev.off()
tiff(file = "figures/Human_Read_Fraction_Condition.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Condition, log10(Human.Read.Fraction), col2, comp_con, 0.1) ; dev.off()
tiff(file = "figures/Human_Read_Fraction_Type.tiff", width = 450, height = 700, units = "px", res = 200) ; div.plot(meta, Type, log10(Human.Read.Fraction), col3.2, comp_type_2, c(0, 0.3)) ; dev.off()
tiff(file = "figures/Human_Read_Fraction_Activity.tiff", width = 650, height = 700, units = "px", res = 200) ; div.plot(meta.act, Activity, log10(Human.Read.Fraction), col4, comp_act, c(-0.1, 0.2, 0.5)) ; dev.off()

tiff(file = "figures/Human_Read_Fraction_Fn_experience.tiff", width = 300, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.experience, log10(Human.Read.Fraction), col2.1, comp_exp, 0.1) ; dev.off()
tiff(file = "figures/Human_Read_Fraction_Fn_proximity.tiff", width = 650, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.proximity, log10(Human.Read.Fraction), col4.1, comp_prox, c(-0.05, 0.25, 0.55)) ; dev.off()
tiff(file = "figures/Human_Read_Fraction_Fn_pre_post.tiff", width = 450, height = 700, units = "px", res = 200) ; div.plot(meta, Fn.pre.post.simple, log10(Human.Read.Fraction), col3.1, comp_pp_3, c(-0.1, 0.2, 0.5)) ; dev.off()

# Distribution of human read fraction over time (based on F. nucleatum proximity)
# X-axis: Temporal proximity to the observation point of F. nucleatum (week)
# Y-axis: Human read fraction

tiff(file = "figures/Human_Read_Fraction_Scatter_week.tiff", width = 1000, height = 800, units = "px", res = 200) ; ggplot(meta, aes(x = Fn.directionality.strict, y=log10(Human.Read.Fraction), group=Condition, color=Condition)) + 
  geom_smooth(method = "loess", span=0.8, alpha=0.2, se = T) + xlab("") + ylab("") +  scale_color_manual(values = col2)+
  geom_hline(yintercept = c(median(log10(meta$Human.Read.Fraction[meta$Fn.experience=="Non" & meta$Condition=="nonIBD"])), median(log10(meta$Human.Read.Fraction[meta$Fn.experience=="Non" & meta$Condition=="IBD"]))), 
             alpha=0.3, size=1, color=col2, linetype="dotted")+
  geom_vline(xintercept = 0, color="gray", size=1, alpha=0.3)+
  theme_classic() + 
  theme(legend.position = c(1,0.85),
        legend.text = element_text(size=rel(1.4)), 
        legend.title = element_blank(),
        legend.justification = c(1,0),
        legend.background = element_rect(color = "black", size=0.5),
        legend.key = element_rect(fill = "white", color="transparent"),
        legend.key.width = unit(0.8, "line"),
        legend.key.height = unit(0.8,"line"), 
        legend.spacing.y = unit(0, "mm"), 
        legend.spacing.x = unit(1, "mm"),
        axis.text.x = element_text(size=rel(1.2), color="black"),
        axis.text.y = element_text(size=rel(1.2), color="black"), 
        axis.title = element_text(size=rel(1.6))) ; dev.off() 



#### 5. Fusobacterium nucleatum ####

# 5-1) Profiling F. nucleatum-detected samples (and subjects)

# Select samples that actually detected F. nucleatum
Fn_pos <- meta$External.ID[meta$Fn.proximity=="Detected"] 
Fn_meta <- as.data.frame(t(bugs["Fusobacterium_nucleatum",Fn_pos])) ; Fn_meta$External.ID <- rownames(Fn_meta) ; Fn_meta <- Fn_meta[order(Fn_meta$External.ID),] # Getting F.nucleatum abundance data
Fn_meta <- cbind(Fn_meta, subset(meta.order, select = c("Condition", "Type", "Activity", "Sex", "Site", "Participant.ID", "Human.Read.Fraction", "Shannon"))[meta.order$Fn.proximity=="Detected",])
  head(Fn_meta)

# Sorting F. nucleatum-expereinced subject 
Fn_subj <- Fn_meta[!duplicated(Fn_meta$Participant.ID),] ; meta_subj <- meta[!duplicated(meta$Participant.ID),]

# Profiling sample information by F. nucleatum observation
    # Sex
    table(Fn_meta$Sex) # F 18 ; M 23
    table(meta$Sex) # F 772 ; M 754
    table(Fn_subj$Sex) # F 9 ; M 10
    table(meta_subj$Sex) # F 52 ; M 54
    
    # Type
    table(Fn_meta$Type) #nonIBD 7 : UC 6 : CD 28
    table(meta$Type) # nonIBD 407 : UC 417 : CD 702
    table(Fn_count$Type) # nonIBD 4 : UC 4 : CD 10
    table(meta_subj$Type) # nonIBD 26 : UC 30 : CD 50
    
    # Site
    table(Fn_meta$Site) # Cedar 12 : Cincinnati 13 : Emory 1 : MGH 13 : MGH Pediatrics 2
    table(meta$Site) # Cedar 390 : Cincinnati 448 : Emory 49 : MGH 426 : MGH Pediatrics 213
    table(Fn_count$Site) # Cedar 3 : Cincinnati 9 ; Emory 1 ; MGH 4 ; MGH Pediatrics 2
    table(meta_subj$Site) # Cedar 27 : Cincinnati 31 ; Emory 4 ; MGH 28 ; MGH Pediatrics 16
    
    # Condition
    table(Fn_meta$Condition) # nonIBD 7 ; IBD 34
    table(meta$Condition) # nonIBD 407 ; IBD 1119
    table(Fn_count$Condition) # nonIBD 4 ; IBD 15
    table(meta_subj$Condition) # nonIBD 26 ; IBD 80


# 5-2) Distribution of F. nucleatum with or without IBD 
tiff("figures/pie_condition.tiff", res = 400, width = 600, height = 600) ; ggplot(meta, aes(x=1, y = Condition, fill=Condition)) + geom_bar(stat = "identity", alpha=0.6) + xlab("") + ylab("") +
  scale_fill_manual(values = col2) + coord_polar(theta="y") + theme_void() + theme(legend.position = "none") ; dev.off()

tiff("figures/pie_condition_upon_detection.tiff", res = 400, width = 600, height = 600) ; ggplot(Fn_meta, aes(x=1, y = Condition, fill=Condition)) + geom_bar(stat = "identity", alpha=0.9) + xlab("") + ylab("") +
  scale_fill_manual(values = col2) + coord_polar(theta="y") + theme_void() + theme(legend.position = "none") ; dev.off()

# 5-3) Distribution of F.nucleatum upon detection
Fn_meta <- Fn_meta[order(Fn_meta$Fusobacterium_nucleatum, decreasing = T),]  ; Fn.order <- Fn_meta$External.ID # Decreasing F. nucleatum abundance order
tiff(file = "figures/Fn_upon_detection.tiff", width = 800, height = 600, units = "px", res = 200) ; ggplot(Fn_meta, aes(External.ID, log10(Fusobacterium_nucleatum))) + xlab("") + ylab("") + 
  geom_point(alpha=0.6, aes(color=Condition)) +
  scale_x_discrete(limits=Fn.order) + scale_y_continuous(limits = c(-4, 2)) + scale_color_manual(values = col2) + 
  theme_classic() +
  theme(panel.background = element_rect(color="black", linetype = "solid", size = 1),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = "gray95"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.80, 0.80),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.box.background = element_rect(color="transparent", fill = "transparent"),
        legend.key.height = unit(1,"line"),
        legend.text = element_text(size=rel(1.3))) ; dev.off()

# 5-4) F. nucleatum enrichment test in IBD condition
# 1) Differential observation p-value
# Contingency table: F. nucleatum observation by IBD
Fn.2x2_IBD <- data.frame(F.n=c("Observed", "None"), nonIBD=c(7, 400), IBD=c(34, 1085)) ; Fn.2x2mat <- as.matrix(Fn.2x2_IBD[,2:3])
fisher.test(x = Fn.2x2mat, conf.int = 0.95, alternative = "less") # Odd ratio=0.5586331 ; p-value = 0.1062

# 2) Upon detection, Wilcoxon test
wilcox.test(Fn_meta$Fusobacterium_nucleatum~Fn_meta$Condition, alternative="less") # 0.02891

# 3) Upon detection, t-test
t.test(Fn_meta$Fusobacterium_nucleatum~Fn_meta$Condition, alternative="less") # 0.005788


#### 6. IBD/non-IBD marker analysis ####
#### 6-1) Screening by LEfSe ####

tiff("figures/lefse_marker.tiff", res = 400, width = 1000, height = 1000, units = "px") ; ggplot(lef, aes(Microbes, Log10_LDA, color=Class)) + 
  xlab("") + ylab("") + 
  #geom_bar(stat = "identity") + 
  geom_segment(aes(x = Microbes, xend=Microbes, y = 0, yend=Log10_LDA), color="gray12", linetype="dotted", size=0.3, alpha=0.4) +
  geom_point() +
  scale_x_discrete(limits=lef_order) +  
  scale_y_continuous(breaks = seq(1, 4.5, 1)) + scale_color_manual(values=col2.2, name=c("nonIBD", "IBD")) +  
  theme_few() +
  theme(axis.text.x = element_text(angle = 90, color="black", hjust=0.95, vjust=0.2, face="italic", size=rel(0.4)),
        axis.text.y = element_text(color="black", size=rel(0.7)),
        legend.position = c(0.91, 0.88),
        legend.background = element_rect(fill = "transparent"),
        legend.key.width = unit(0.4, "line"),
        legend.key.height = unit(0.25, "line"),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.35))) ; dev.off()


# Merging metadata with IBD/non-IBD marker species profiling data
# 1526 x 26 (12 IBD markers + 14 non-IBD markers)
tbugs.ibd <- tbugs[,ibd.mark]; tbugs.non <- tbugs[,non.mark]
tbugs.mark <- cbind(tbugs.ibd, tbugs.non)
tbugs.mark$External.ID <- rownames(tbugs.mark) # To merge using External.ID 

mark.meta <- merge(tbugs.mark, meta, all = T, by = "External.ID")
mark.meta <- mark.meta[order(mark.meta$Participant.ID, mark.meta$visit.num),]

# Add marker information
# 1) Sum of species
mark.meta$specIBD <- rowSums(mark.meta[,ibd.mark]!=0) ; mark.meta$specNON <- rowSums(mark.meta[,non.mark]!=0) 

# 2) Sum of abundance
mark.meta$sumIBD <- rowSums(mark.meta[,ibd.mark]) ; mark.meta$sumNON <- rowSums(mark.meta[,non.mark])

# 3) Ratio between markers: IBD over nonIBD
mark.meta$specRatio <- (mark.meta$specIBD+1)/(mark.meta$specNON+1) ; mark.meta$sumRatio <- (mark.meta$sumIBD+0.00001)/(mark.meta$sumNON+0.00001)
  #mark.meta.act <- mark.meta[!is.na(mark.meta$Activity),] ; mark.meta.act$Activity <- droplevels(mark.meta.act$Activity) 
  #write.csv(mark.meta, "mark.meta.csv") 
  #mark.meta <- read.csv("input/mark.meta.csv", header = T, row.names = 1) ; rownames(mark.meta) <- NULL ; head(mark.meta) # Import csv files containing IBD/non-IBD markers + metadata 
  

# F. nucleatum proximity ~ Abundance/Number of species
mark.tmp <- mark.meta[, c(mark.list, "External.ID", "Week", "Condition", "Fn.directionality", "Fn.directionality.strict", "Shannon", "Human.Read.Fraction")]
mark.melt <- melt(mark.tmp, id.vars = c("External.ID","Condition", "Week", "Fn.directionality", "Fn.directionality.strict"), measure.vars = mark.list)
mark.melt$Mark.type <- "nonIBD" ; mark.melt$Mark.type[mark.melt$variable %in% ibd.mark] <- "IBD" ; mark.melt$Mark.type <- factor(mark.melt$Mark.type, levels = c("nonIBD", "IBD")) ; colnames(mark.melt)[6:7] <- c("microbes", "Abundance")

mark.plot <- function(data, x, y, comparisons, color, label.y) {
  attach(data)
  ggplot(data, aes(x = x, y = y, fill=x)) + xlab("") + ylab("") + 
    geom_boxplot(alpha=0.6, outlier.alpha = 0.02) + 
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format", tip.length = 0, label.y=label.y) +
    scale_fill_manual(values=color) +
    theme_classic() +
    theme(axis.text.x=element_text(size=rel(1.3), vjust = 0.65, hjust = 0.5, color="black"), 
          axis.text.y=element_text(size=rel(1.1), color="black"), 
          axis.ticks = element_blank(),
          legend.position = "none")
}

# Profiling the distribution of IBD/non-IBD marker species----
scat_box.plot <- function(data, x, y, variable, color) {
  attach(data)
  # 1) Scatter plot at center
  scat1 <- ggplot(data, aes(x, y, fill=variable, color=variable)) + 
    xlab("") + ylab("")+ geom_jitter(alpha=0.3) + 
    # X-axis adjustment
    scale_x_continuous(expand = c(0,0)) +
    expand_limits(x=c(min(x) - 0.1*diff(range(x)),max(x) + 0.1*diff(range(x)))) +
    # Y-axis adjustment
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y=c(min(y) - 0.1*diff(range(y)),max(y) + 0.1*diff(range(y)))) +
    scale_color_manual(values=color) + scale_fill_manual(values=color) +
    theme_classic2() +
    theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"),
          axis.text = element_text(size=rel(0.8), color="black"), axis.line = element_line(size=0.5),
          legend.position = "bottom", legend.margin = margin(20, 0, 0, 0), legend.title = element_blank(), legend.key = element_rect(fill="white", linetype = "blank"), legend.text = element_text(size=rel(1.4)))
  
  # 2) Vertical marginal boxplot (at right)
  box1 <- ggplot(data, aes(variable, y, fill=variable)) + 
    geom_boxplot(outlier.color = NA, alpha=0.7) + 
    scale_fill_manual(values=color) + scale_y_continuous(expand = c(0,0)) +
    stat_compare_means(comparisons = list(c("nonIBD", "IBD")), method = "wilcox.test", label="p.signif", tip.length = 0, label.y = 12.5, size=1.5) +
    expand_limits(y=c(min(y) - 0.1*diff(range(y)),max(y) + 0.1*diff(range(y)))) +
    theme_void() +
    theme(legend.position = "none", axis.text=element_blank(),axis.title = element_blank(), axis.ticks = element_blank(),
          plot.margin = unit(c(0.2, 1, 0.5, -0.5), "lines")) 
  
  
  # 3) Horizontal marginal boxplot (at top)
  box2 <- ggplot(data, aes(variable, x, fill=variable)) + 
    geom_boxplot(outlier.color = NA, alpha=0.7) + 
    scale_fill_manual(values=color) + scale_y_continuous(expand = c(0,0)) +
    stat_compare_means(comparisons = list(c("nonIBD", "IBD")), method = "wilcox.test", label="p.signif", tip.length = 0, label.y = 14.5, size=1.5) +
    expand_limits(y=c(min(x) - 0.1*diff(range(x)), max(x) + 0.1*diff(range(x)))) +
    coord_flip()+ theme_void() +
    theme(legend.position = "none",
          axis.text=element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(c(1,0.2, -0.5, 0.5), "lines"))
  
  gt1 <- ggplot_gtable(ggplot_build(scat1)) ; gt2 <- ggplot_gtable(ggplot_build(box2)) ; gt3 <- ggplot_gtable(ggplot_build(box1))
  
  # Get maximum widths and heights
  gWidth <- unit.pmax(gt1$widths, gt2$widths) ; gHeight <- unit.pmax(gt1$heights, gt3$heights)
  
  # Set the maximums in the gtables for gt1, gt2, and gt3
  gt1$widths <- as.list(gWidth) ; gt2$widths <- as.list(gWidth) ; gt1$heights <- as.list(gHeight) ; gt3$heights <- as.list(gHeight)
  
  # Creat a new gtable 
  gt <- gtable(widths = unit(c(6.5,1), "null"), heights = unit(c(1,8), "null"))
  
  # Insert gt1, gt2, gt3 into the new gtable
  gt <- gtable_add_grob(gt, gt1, 2, 1) ; gt <- gtable_add_grob(gt, gt2, 1, 1) ; gt <- gtable_add_grob(gt, gt3, 2, 2)
  
  # Render the plot
  grid.newpage() ; grid.draw(gt)
}
tiff("figures/IBD_nonIBD_marker_scatter.tiff", width = 800, height = 900, units = "px", res=200) ; scat_box.plot(mark.meta, specNON, specIBD, Condition, col1) ; dev.off()

# by IBD/nonIBD 
# 1) Sum of abundance
tiff(file = "figures/Marker_Condition_sumIBD.tiff", width = 300, height = 700, units = "px", res = 200) ; mark.plot(mark.meta, Condition, log10(sumIBD+0.00001), comp_con, col2, 2.05)  ; dev.off() 
tiff(file = "figures/Marker_Condition_sumNON.tiff", width = 300, height = 700, units = "px", res = 200) ; mark.plot(mark.meta, Condition, log10(sumNON+0.00001), comp_con, col2, 2.05)  ; dev.off() 

# 2) Number of species
tiff(file = "figures/Marker_Condition_specIBD.tiff", width = 300, height = 700, units = "px", res = 200) ; mark.plot(mark.meta, Condition, specIBD, comp_con, col2, 12.5)  ; dev.off() 
tiff(file = "figures/Marker_Condition_specNON.tiff", width = 300, height = 700, units = "px", res = 200) ; mark.plot(mark.meta, Condition, specNON, comp_con, col2, 14.5)  ; dev.off() 

# 3) Ratio of IBD over non-IBD
tiff(file = "figures/Marker_Condition_specRatio.tiff", width = 330, height = 700, units = "px", res = 200) ; mark.plot(mark.meta, Condition, specRatio, comp_con, col2, 9.5) +
  geom_hline(yintercept = 1, size=1, color="gray20", alpha=0.2) ; dev.off()
tiff(file = "figures/Marker_Condition_sumRatio.tiff", width = 330, height = 700, units = "px", res = 200) ; mark.plot(mark.meta, Condition, log10(sumRatio), comp_con, col2, 4.2) +
  geom_hline(yintercept = 0, size=1, color="gray20", alpha=0.2); dev.off() 


# by F.nucleatum exposure
# 1) IBD
tiff(file = "figures/Marker_Fn_experience_sumIBD_by_Condition.tiff", width = 500, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.experience, log10(sumIBD+0.00001), comp_exp, col2.1, 2.05) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")); dev.off() 
tiff(file = "figures/Marker_Fn_experience_specIBD_by_Condition.tiff", width = 500, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.experience, specIBD, comp_exp, col2.1, 12.5) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")); dev.off() 

# 2) nonIBD
tiff(file = "figures/Marker_Fn_experience_sumNON_by_Condition.tiff", width = 500, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.experience, log10(sumNON+0.00001), comp_exp, col2.1, 2.05) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")); dev.off() 
tiff(file = "figures/Marker_Fn_experience_specNON_by_Condition.tiff", width = 500, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.experience, specNON, comp_exp, col2.1, 14.5) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")); dev.off() 

# 3) Ratio
tiff(file = "figures/Marker_Fn_experience_sumRatio_by_Condition.tiff", width = 500, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.experience, log10(sumRatio), comp_exp, col2.1, 7.5) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")) +
  geom_hline(yintercept = 0, size=0.8, color="gray20", alpha=0.2, linetype="dotted"); dev.off() 
tiff(file = "figures/Marker_Fn_experience_specRatio_by_Condition.tiff", width = 500, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.experience, log10(specRatio), comp_exp, col2.1, 1) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")) +
  geom_hline(yintercept = 0, size=0.8, color="gray20", alpha=0.2, linetype="dotted") ; dev.off()


# by Prior/Posterior
# 1) IBD
tiff(file = "figures/Marker_Fn_pre_post_sumIBD_by_Condition.tiff", width = 900, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.pre.post.simple, log10(sumIBD+0.00001), comp_pp_2, col3, c(2.05, 2.5)) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")); dev.off() 
tiff(file = "figures/Marker_Fn_pre_post_specIBD_by_Condition.tiff", width = 900, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.pre.post.simple, specIBD, comp_pp_2, col3, c(12.5, 13.5)) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")); dev.off() 

# 2) nonIBD
tiff(file = "figures/Marker_Fn_pre_post_sumNON_by_Condition.tiff", width = 900, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.pre.post.simple, log10(sumNON+0.00001), comp_pp_2, col3, c(2.05, 2.55)) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")); dev.off() 
tiff(file = "figures/Marker_Fn_pre_post_specNON_by_Condition.tiff", width = 900, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.pre.post.simple, specNON, comp_pp_2, col3, c(14.2, 15.15)) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")); dev.off() 

# 3) Ratio
tiff(file = "figures/Marker_Fn_pre_post_sumRatio_by_Condition.tiff", width = 900, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.pre.post.simple, log10(sumRatio), comp_pp_2, col3, c(7.4, 8.1)) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent")) +
  geom_hline(yintercept = 0, size=0.8, color="gray20", alpha=0.2, linetype="dotted"); dev.off() 
tiff(file = "figures/Marker_Fn_pre_post_specRatio_by_Condition.tiff", width = 900, height = 800, units = "px", res = 200) ; mark.plot(mark.meta, Fn.pre.post.simple, log10(specRatio), comp_pp_2, col3, c(1, 1.15)) + facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.1), margin = margin(1.5, 0, 0.2, 0, "line")), strip.background = element_rect(color="transparent"))+
  geom_hline(yintercept = 0, size=0.8, color="gray20", alpha=0.2, linetype="dotted"); dev.off() 


# Directionality & marker species--------------------
direct.plot <- function(data, x, y) {
  attach(data)
  ggplot(data) + xlab("") + ylab("") + 
    geom_smooth(aes(x = x, y = y, group=variable, color=variable), method="loess", alpha=0.3, fullrange=T, span=0.9) +
    scale_x_continuous(limits = c(-50, 50)) + 
    scale_color_manual(values = col2)+
    facet_wrap(~Condition, ncol=2) +
    geom_vline(xintercept = 0, alpha=0.4, color="gray80") +
    theme_classic() + 
    theme(axis.text.x=element_text(size=rel(1.2), color="black"), 
          axis.text.y=element_text(size=rel(1.2), color="black"), 
          strip.text = element_text(size=rel(1.5)), 
          axis.title.x = element_text(size=rel(1.5), vjust = -1),
          panel.spacing = unit(1.5, "lines"),
          strip.background = element_rect(color="transparent"),
          legend.position = "none")
}

# Distribution of IBD marker species over time 
  # Fn.directionality.strict considered the very first F. nucleatum detection as a 0-week. 

# Melting data
melt_sum <- melt(mark.meta, id.vars = c("External.ID", "Fn.directionality", "Fn.directionality.strict", "Condition"), measure.vars = c("sumIBD", "sumNON"))
melt_spec <- melt(mark.meta, id.vars = c("External.ID", "Fn.directionality", "Fn.directionality.strict", "Condition"), measure.vars = c("specIBD", "specNON"))
melt_sum$variable <- factor(melt_sum$variable, levels=c("sumNON", "sumIBD")) ; melt_spec$variable <- factor(melt_spec$variable, levels=c("specNON", "specIBD")) 

  # Marker species information only in F.nucleatum-"innocent" subject
  mark.meta$sumIBD_woFN <- NA ; mark.meta$sumIBD_woFN[mark.meta$Fn.experience=="Non"] <- rowSums(mark.meta[mark.meta$Fn.experience=="Non", ibd.mark])
  mark.meta$sumNON_woFN <- NA ; mark.meta$sumNON_woFN[mark.meta$Fn.experience=="Non"] <- rowSums(mark.meta[mark.meta$Fn.experience=="Non", non.mark])
  mark.meta$specIBD_woFN <- NA ; mark.meta$specIBD_woFN[mark.meta$Fn.experience=="Non"] <- rowSums(mark.meta[mark.meta$Fn.experience=="Non", ibd.mark]>0)
  mark.meta$specNON_woFN <- NA ; mark.meta$specNON_woFN[mark.meta$Fn.experience=="Non"] <- rowSums(mark.meta[mark.meta$Fn.experience=="Non",non.mark]>0)
  
  # Sum of abundance at F. nucleatum-innocent subjects
  median(mark.meta$sumIBD_woFN[mark.meta$Condition=="IBD"], na.rm = T) # IBD markers at IBD ; 18.41024 
  median(mark.meta$sumIBD_woFN[mark.meta$Condition=="nonIBD"], na.rm = T) # IBD markers at nonIBD ; 8.78698
  
  median(mark.meta$sumNON_woFN[mark.meta$Condition=="IBD"], na.rm = T) # nonIBD markers at IBD ; 9.71617 
  median(mark.meta$sumNON_woFN[mark.meta$Condition=="nonIBD"], na.rm = T) # nonIBD markers at nonIBD ; 17.75386
  
  # Number of detected markers at F. nucleatum-innocent subjects
  median(mark.meta$specIBD_woFN[mark.meta$Condition=="IBD"], na.rm = T) # IBD markers at IBD ; 8 
  median(mark.meta$specIBD_woFN[mark.meta$Condition=="nonIBD"], na.rm = T) # IBD markers at nonIBD ; 6 
  
  median(mark.meta$specNON_woFN[mark.meta$Condition=="IBD"], na.rm = T) # nonIBD markers at IBD ; 7
  median(mark.meta$specNON_woFN[mark.meta$Condition=="nonIBD"], na.rm = T) # nonIBD markers at nonIBD ; 11


# dataframe for horizontal lines   
line_sum <- data.frame(Condition = rep(c("nonIBD", "IBD"), each=2), Z = c(17.75386, 8.78698, 18.41024, 9.71617)); print(line_sum) 
line_spec <- data.frame(Condition = rep(c("nonIBD", "IBD"), each=2), Z = c(11, 6, 8, 7)); print(line_spec)
line_sum$Condition <- factor(line_sum$Condition, levels=c("nonIBD", "IBD")) ; line_spec$Condition <- factor(line_spec$Condition, levels=c("nonIBD", "IBD"))


# Distribution of marker species over time facetted by IBD/nonIBD conditions
tiff(file = "figures/Marker_sum.tiff", width = 900, height = 800, units = "px", res = 200) ; direct.plot(melt_sum, x=Fn.directionality.strict, y=value) +
  geom_hline(data=line_sum, aes(yintercept = Z), color=c(col2, rev(col2)), alpha=0.5, size=1.2, linetype="dotted"); dev.off()

tiff(file = "figures/Marker_spec.tiff", width = 900, height = 800, units = "px", res = 200); direct.plot(melt_spec, x=Fn.directionality.strict, y=value) +
  geom_hline(data=line_spec, aes(yintercept = Z), color=c(col2, rev(col2)), alpha=0.5, size=1.2, linetype="dotted"); dev.off()


#### 6-2) Spearman correlation of all microbes with F. nucleatum ####
t.mat <- as.matrix(tbugs)

# Calculate spearman correlation coefficient between microbes
cor.mat <- rcorr(t.mat, type = "spearman") ; str(cor.mat) # r (correlation) + n (1526; number of samples) + p (p-value)
cor.mat.fn <- cor.mat$r["Fusobacterium_nucleatum",] ; cor.mat.fn.p <- cor.mat$P["Fusobacterium_nucleatum",] ; cor.mat.fn.fdr <- p.adjust(cor.mat.fn.p, method="fdr", n=length(cor.mat.fn.p))
  bugs.spearman <- data.frame(Spearman=cor.mat.fn, p.value=cor.mat.fn.p, fdr=cor.mat.fn.fdr, NOD=NOD) 
  #bugs.spearman <- bugs.spearman[!is.na(bugs.spearman$p.value),] # Remove Fusobacterium nucleautm itself ; 532 now

# Add marker information
bugs.spearman$Class <- "Neither" ; bugs.spearman$Class[rownames(bugs.spearman) %in% ibd.mark] <- "IBD" ; bugs.spearman$Class[rownames(bugs.spearman) %in% non.mark] <- "nonIBD"
bugs.spearman$Class <- factor(bugs.spearman$Class, levels = c("Neither", "nonIBD", "IBD")) ; head(bugs.spearman)
  #write.csv(bugs.spearman, "input/spearman.fn.csv")
  #bugs.spearman <- read.csv("input/spearman.fn.csv", header=T, row.names=1) 

tiff(file = "figures/Spearman_Mark.tiff", width = 500, height = 700, units = "px", res = 200) ; div.plot(bugs.spearman, Class, Spearman, col3, comp_con, 0.18) ; dev.off()

# 2) P-value for differential enrichment for IBD/nonIBD
tb_con <- cbind(tbugs.order, Condition=meta.order$Condition) ; dim(tb_con) # 1526 x 534 

mic_con <- as.data.frame(matrix(nrow = 533, ncol=9)) ; for (i in 1:533) {
  x <- tb_con[tb_con$Condition=="IBD", i] ; y <- tb_con[tb_con$Condition=="nonIBD", i]
  
  # Wilcoxon test
  mic_con[i,1] <- wilcox.test(x, y, paired=F, alternative="greater")$p.value    # IBD-enriched
  mic_con[i,2] <- wilcox.test(x, y, paired=F, alternative="less")$p.value       # nonIBD-enriched
  mic_con[i,3] <- wilcox.test(x, y, paired=F, alternative="two.sided")$p.value  # two.sided
  
  # Log transformation & multiple test correction by FDR
  mic_con[i,4] <- -log10(mic_con[i,1])
  mic_con[i,5] <- -log10(mic_con[i,2])
  mic_con[i,6] <- -log10(mic_con[i,3])
  mic_con[i,7] <- p.adjust(mic_con[i,1], method = "fdr", n=533)
  mic_con[i,8] <- p.adjust(mic_con[i,2], method = "fdr", n=533)
  mic_con[i,9] <- p.adjust(mic_con[i,2], method = "fdr", n=533)
} ; rownames(mic_con) <- colnames(tb_con)[1:533] ; colnames(mic_con) <- c("IBD.Pvalue", "nonIBD.Pvalue", "Diff.Pvalue", "logP.IBD", "logP.nonIBD", "logP.Diff", "fdr.IBD", "fdr.nonIBD", "fdr.Diff")
  #write.csv(mic_con, "input/p.value_IBD.nonIBD.bugs.csv")
  #mic_con <- read.csv("input/p.value_IBD.nonIBD.bugs.csv", header=T, row.names=1)

# 3) P-value for differential enrichment for F.nucleatum-exposed samples
tb_exp <- cbind(tbugs.order, Fn.experience=meta.order$Fn.experience)

mic_exp <- as.data.frame(matrix(nrow = 533, ncol=6)) ; for (i in 1:533) {
  x <- tb_exp[tb_exp$Fn.experience=="Exp", i] ;  y <- tb_exp[tb_exp$Fn.experience=="Non", i]
  
  # Wilcoxon test
  mic_exp[i,1] <- wilcox.test(x, y, paired=F, alternative="greater")$p.value    # F.nucleatum-experienced
  mic_exp[i,2] <- wilcox.test(x, y, paired=F, alternative="less")$p.value       # F.nucleatum-innocent
  mic_exp[i,3] <- wilcox.test(x, y, paired=F, alternative="two.sided")$p.value  # two.sided
  
  # FDR
  mic_exp[i,4] <- p.adjust(mic_exp[i,1], method = "fdr", n=533)
  mic_exp[i,5] <- p.adjust(mic_exp[i,2], method = "fdr", n=533)
  mic_exp[i,6] <- p.adjust(mic_exp[i,3], method = "fdr", n=533)
} ; rownames(mic_exp) <- colnames(tb_exp)[1:533]  ; colnames(mic_exp) <- c("Exp.Pvalue", "Non.Pvalue", "Diff.Pvalue", "Exp.fdr", "Non.fdr", "Diff.fdr")
  #write.csv(mic_exp, "input/Fnucleatum_exposed_Pvalue.csv")
  #mic_exp <- read.csv("input/Fnucleatum_exposed_Pvalue.csv", header=T, row.names=1)


# 4) P-value for differential enrichment for F.nucleatum-posterior 
tb_pp <- cbind(tbugs.order, Fn.pre.post.simple=meta.order$Fn.pre.post.simple) 
  # Sorting F. nucleatum-experienced subject 
  tb_pp <- tb_pp[tb_pp$Fn.pre.post.simple!="Negative", ] ; NOD5_Fn <- colSums(tb_pp[,1:533]>0)>=5
  tb_pp <- tb_pp[,NOD5_Fn] ; tb_pp$Fn.pre.post.simple <- meta.order$Fn.pre.post.simple[meta.order$Fn.pre.post.simple!="Negative"]
  tb_pp$Fn.pre.post.simple <- droplevels(tb_pp$Fn.pre.post.simple)
  

mic_pp <- as.data.frame(matrix(nrow = 258, ncol=6)) ; for (i in 1:258) {
  
  x <- tb_pp[tb_pp$Fn.pre.post.simple=="Posterior", i]
  y <- tb_pp[tb_pp$Fn.pre.post.simple=="Prior", i]
  
  # Wilcoxon test
  mic_pp[i,1] <- wilcox.test(x, y, paired=F, alternative="greater")$p.value   # F.nucleatum-posterior
  mic_pp[i,2] <- wilcox.test(x, y, paired=F, alternative="less")$p.value      # F.nucleatum-prior
  mic_pp[i,3] <- wilcox.test(x, y, paired=F, alternative="two.sided")$p.value # two.sided
  
  # FDR
  mic_pp[i,4] <- p.adjust(mic_pp[i,1], method = "fdr", n=258)
  mic_pp[i,5] <- p.adjust(mic_pp[i,2], method = "fdr", n=258)
  mic_pp[i,6] <- p.adjust(mic_pp[i,3], method = "fdr", n=258)
} ; rownames(mic_pp) <- colnames(tb_pp)[1:258] ; colnames(mic_pp) <- c("Post.Pvalue", "Pre.Pvalue", "Diff.Pvalue", "Post.fdr", "Pre.fdr", "Diff.fdr")
  #write.csv(mic_pp, "input/Fnucleatum_pre_post_Pvalue.csv")
  #mic_pp <- read.csv("input/Fnucleatum_pre_post_Pvalue.csv", header=T, row.names=1) ; head(mic_pp, n=10)

# Merging correlation & enrichment data
  #bugs.spearman <- read.csv("input/spearman.fn.csv", header=T, row.names=1) 
  #mic_con <- read.csv("input/p.value_IBD.nonIBD.bugs.csv", header = T, row.names = 1) ; head(mic_con) # IBD/nonIBD enrichment

#which(rownames(mic_con)=="Fusobacterium_nucleatum") ; mic_con <- mic_con[-251,]  # To remove F. nucleatum
bugs.spearman$microbes <- rownames(bugs.spearman) ; mic_con$microbes <- rownames(mic_con)
bugs.spearman <- merge(bugs.spearman, mic_con, by = "microbes"); bugs.spearman$Class <- factor(bugs.spearman$Class, levels = c("nonIBD", "IBD", "Neither")) ; head(bugs.spearman)

#  Drawing Volcano plot----
# X-axis: Spearman correlation with F. nucleatum
# Y-axis: p-value for IBD/non-IBD differential enrichment (two-sided)
tiff(file = "figures/Volcano_Markers_Enrichment_vs_Correlation.tiff", width = 1200, height = 1000, units = "px", res = 200) ; ggplot(bugs.spearman) + xlab("") + ylab("") + 
  geom_point(aes(Spearman, logP.Diff, color=Class, size=NOD), alpha=0.7) +
  scale_x_continuous(limits = c(-0.25, 0.25)) +  scale_color_manual(values=c("#3399FF", "#CC3333", "gray92")) +
  geom_vline(xintercept = 0, alpha=0.6, size=1.5, color="gray80") +
  stat_smooth(aes(Spearman, logP.Diff), geom = "smooth", method="loess", se=T, span=0.9, alpha=0.3, color="gray60")+
  stat_cor(method = "spearman", aes(x = abs(Spearman), y = logP.Diff), label.x = -0.25, label.y = 32, size=7)+
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=rel(1.1)),
        legend.position = c(0.95,0.98), 
        legend.justification = c(0.95,0.95),
        legend.title = element_text(size=rel(1.1)),
        legend.background = element_rect(color = "black", size=0.5),
        legend.key.height=unit(1,"line"), # default = 3 ; spacing between the key items
        legend.spacing.y = unit(0.2, "line"),
        legend.text = element_text(size=rel(1)),
        legend.key = element_rect(fill = "transparent", color="transparent")) ; dev.off()


#  Distribution of IBD/non-IBD marker species over time
  # Sub-divided into  into 3 categories
  # 1) non-IBD markers
  # 2) IBD marker 1 = IBD markers, positively correlated with F. nucleatum
  # 3) IBD marker 2 = IBD markers, non-positively correlated with F. nucleatum
    # Bacteroides fragilis showed 0.005672 spearman correlation coefficient with F. nucleatum (p-value 0.8247)
    # Thus, the bug was included in IBD marker 2, so called Non-positively associated IBD markers
    ibd.nonPos <- as.character(bugs.spearman[bugs.spearman$Class=="IBD" & bugs.spearman$Spearman < 0.01, ]$microbes)  # 6 bugs
    ibd.Pos <- as.character(bugs.spearman[bugs.spearman$Class=="IBD" & bugs.spearman$Spearman > 0.01, ]$microbes)     # 6 bugs
    
    # Add information 
    mark.meta$nonPos_IBD.sum <- rowSums(mark.meta[,ibd.nonPos]) ; mark.meta$nonPos_spec <- rowSums(mark.meta[,ibd.nonPos]>0)
    mark.meta$Pos_IBD.sum <- rowSums(mark.meta[,ibd.Pos]) ; mark.meta$Pos_spec <- rowSums(mark.meta[,ibd.Pos]>0)
    
    # Melting
    mark.melt <- melt(mark.meta, id.vars = c("Condition", "Fn.directionality.strict"), measure.vars = c("sumNON", "nonPos_IBD.sum", "Pos_IBD.sum")) ; mark.melt <- mark.melt[is.na(mark.melt$Fn.directionality.strict)==F,]
    mark.melt$Strength[mark.melt$variable=="sumNON"] <- "nonIBD"
    mark.melt$Strength[mark.melt$variable=="nonPos_IBD.sum"] <- "Non-P"
    mark.melt$Strength[mark.melt$variable=="Pos_IBD.sum"] <- "Positive"
    mark.melt$Strength <- factor(mark.melt$Strength, levels = c("nonIBD", "Positive", "Non-P"))


tiff("figures/Markers_correlation_over_time.tiff", width = 1000, height = 1000, units = "px", res = 200) ; ggplot(mark.melt, aes(Fn.directionality.strict, log10(value), color=Strength)) +
  geom_jitter(alpha=0.02) +
  geom_vline(xintercept = 0, color="gray10", alpha=0.1, size=1) +
  stat_smooth(method="loess", span=0.9, alpha=0.3, se=T, size=1.4, aes(group=Condition)) +
  scale_x_continuous(limits = c(-40, 40)) + scale_y_continuous(limits = c(-3, 3)) + 
  scale_color_manual(values = col2.3) +
  facet_grid(Condition~Strength) +
  stat_cor(method = "spearman", color="black", size=3.3, alpha=0.8, label.x = -37)+
  theme_classic() +
  theme(strip.text = element_text(size=rel(1.1)), 
        axis.text = element_text(color="black", size=rel(0.8)),
        axis.title = element_blank(),
        legend.position = "none") ; dev.off()


# 1) Distribution of IBD/non-IBD marker species
# Type of markers
all_markers <- c(ibd.mark, non.mark) ; all_markers <- data.frame(microbe=all_markers) ; all_markers$Condition <- ifelse(all_markers$microbe %in% non.mark, yes = "nonIBD", no = "IBD")
  melt_marker <- melt(mark.meta, id.vars = c("Participant.ID", "External.ID","Week", "Condition", "Fn.directionality.strict"), measure.vars = c(ibd.mark, non.mark))
  melt_marker_Exp <- melt_marker[!is.na(melt_marker$Fn.directionality.strict),] # F. nucleatum-experienced subject
  melt_marker_Non <- melt_marker[is.na(melt_marker$Fn.directionality.strict), ] # F. nucleatum-inexperienced subject

draw_microbes <- function(data, date){
  attach(data)
  ggplot(data, aes(x=date, y=log10(value+0.0001), color=Condition)) +
    xlab("") + ylab("")+
    stat_smooth(geom = "smooth", method = "loess", alpha=0.2, se = T, span=0.9, size=1.5) + 
    stat_cor(method = "spearman", size=6, label.y = c(2.1, 1.55)) +
    scale_color_manual(values = col2) + scale_y_continuous(limits = c(-4.05, 2.3))+
    theme_classic() +
    theme(legend.position = "none", 
          axis.text = element_text(color="black", size = rel(1.2)))
}


for(i in 1:length(all_markers$microbe)){
  tmp_bug <- all_markers$microbe[i]
  tmp_condition <- all_markers$Condition[all_markers$microbe==tmp_bug]
  
  # 1) Exposed
  assign(paste0(tmp_bug, "_marker"), melt_marker_Exp[melt_marker_Exp$variable==tmp_bug,])
  tmp_figure <- draw_microbes(get(paste0(tmp_bug, "_marker")), Fn.directionality.strict) + 
    geom_vline(xintercept = 0, color="gray35", alpha=0.1, size=0.7) +
    scale_x_continuous(limits = c(-40, 40))
  ggsave(tmp_figure, path = "figures/05_IBD_nonIBD_marker/4_all markers/Exp/", filename = paste0(tmp_bug, "_", tmp_condition, "_Fn_exp_subject.tiff"), width = 110, height = 110, units = "mm", dpi = 200)
  
  # 2) Non-exposed
  assign(paste0(tmp_bug, "_marker"), melt_marker_Non[melt_marker_Non$variable==tmp_bug,])
  tmp_figure <- draw_microbes(get(paste0(tmp_bug, "_marker")), Week)
  ggsave(tmp_figure, path = "figures/05_IBD_nonIBD_marker/4_all markers/Non/", filename = paste0(tmp_bug, "_", tmp_condition, "_Fn_inexp_subject.tiff"), width = 110, height = 110, units = "mm", dpi = 200)
}


#### 7. Dissimilarity analysis #####
# Beta-diversity
# Bray-Curtis distance using microbial abundance/composition

# Generating distance matrix
d.log <- as.matrix(vegdist(log10(tbugs.order+1e-05)+5, method = "bray", upper = T)) # Distance between samples

#### 7-1) Intra-personal dissimilarity: random selection ####

# Step.1 = Select a random patient
# Step.2 = Select two random samples in the same patient (selected step 1)
# Step 3 = Comparing the samples
  # Step.3-1) Calculate microbial distance  (Bring the pre-calculated distance out of previously-generated distance matrix)
  # Step.3-2) Calculate Temporal distance (week) 
  # Step 3-3) Get other information of subjects (differences in human read fraction ; inflammatory condition)

# Setting
set.seed(150) 
dist.df <- as.data.frame(matrix(ncol = 9)) # Container
pat <- as.character(meta$Participant.ID) # Participant ID 
  pat.fn <- as.character(meta$Participant.ID[meta$Fn.experience=="Exp"]) # F. nucleatum-experienced subjects
  identical(colnames(d.log), as.character(meta.order$External.ID)) # TRUE

for (i in 1:10000){
  pat.tmp <- sample(pat, 1)  # Step.1 : Choose a subject
  tmp <- sample(meta.order$External.ID[meta.order$Participant.ID==pat.tmp],2) # Step.2: Choose two samples
  idx1 <- which(meta.order$External.ID==tmp[1]) ; idx2 <- which(meta.order$External.ID==tmp[2])
  
# dist.df = distance data.frame 
  dist.df[i,] <- data.frame(as.character(meta.order$Condition[idx1]), # Condition
                            as.character(meta.order$Fn.experience[idx1]), # Fn.experience
                            d.log[tmp[1], tmp[2]], # Microbial dissimilarity (beta-diversity) stored at matrix
                            meta.order$Week[idx1]-meta.order$Week[idx2], abs(meta.order$Week[idx1]-meta.order$Week[idx2]), # Week difference 
                            meta.order$Human.Read.Fraction[idx1]-meta.order$Human.Read.Fraction[idx2], abs(meta.order$Human.Read.Fraction[idx1]-meta.order$Human.Read.Fraction[idx2]), # Human.Read.Fractino difference
                            meta.order$Shannon[idx1]-meta.order$Shannon[idx2], abs(meta.order$Shannon[idx1]-meta.order$Shannon[idx2])) # Shannon index difference
} ; rm(tmp, idx1, idx2, pat.tmp)

colnames(dist.df) <- c("Condition", "Fn.experience", "Distance", "diff.Week", "diff.Week.abs",  "diff.HumanFraction", "diff.HumanFraction.abs", "diff.Shannon", "diff.Shannon.abs") 
dist.df$Condition <- factor(dist.df$Condition, levels = c("nonIBD", "IBD")) ; dist.df$Fn.experience <- factor(dist.df$Fn.experience, levels = c("Non", "Exp"))

# Merge two columns into one class
dist.df$Class <- "IBD & F.n.-innocent"
dist.df$Class[dist.df$Condition=="IBD" & dist.df$Fn.experience=="Exp"] <- "IBD & F.n.-exp"
dist.df$Class[dist.df$Condition!="IBD" & dist.df$Fn.experience!="Exp"] <- "nonIBD & F.n.-innocent"
dist.df$Class[dist.df$Condition!="IBD" & dist.df$Fn.experience=="Exp"] <- "nonIBD & F.n.-exp"
dist.df$Class <- factor(dist.df$Class, levels=c("IBD & F.n.-innocent", "IBD & F.n.-exp", "nonIBD & F.n.-innocent", "nonIBD & F.n.-exp"))

# Distance scatter plot ----
dist.plot <- function(data, x, y, group, color) {
  attach(data)
  ggplot(data, aes(x, y, group=group, color=group)) + 
    scale_x_continuous(limits = c(0, 40))+
    scale_fill_manual(values=color) + scale_color_manual(values=color) +
    stat_smooth(geom = "smooth", method = "loess", span=0.8, alpha=0.1, se=T, size=1.2) + 
    xlab("") + ylab("") + 
    theme_classic() +
    theme(axis.text.x = element_text(size=rel(1), color="black"), 
          axis.text.y = element_text(size=rel(1), color="black"),
          legend.position = c(0.26,0.88), 
          legend.title = element_blank(),
          legend.text = element_text(size=rel(1)), 
          legend.key = element_rect(fill = "transparent", color="transparent"),
          legend.key.width = unit(0.7, "line"),
          legend.key.height = unit(0.7,"line"), 
          legend.spacing.y = unit(0, "mm"), 
          legend.spacing.x = unit(0.4, "mm"))
}

# All-in-one figure
tiff(file = "figures/Intra_personal_dist_Condition-Fn_exp_all-in-one.tiff", width = 900, height = 700, units = "px", res = 200) ; dist.plot(dist.df, diff.Week.abs, Distance, Class, col4.2) ; dev.off() 

  # Facetted by Condition
  tiff(file = "figures/Intra_personal_dist_Condition-Fn_exp.tiff", width = 1100, height = 600, units = "px", res = 200) ; dist.plot(dist.df, diff.Week.abs, Distance, Condition, col2) + 
    facet_wrap(~Fn.experience) + theme(strip.text = element_text(size=rel(1.3)), legend.position = "none"); dev.off() 
  
  # Facetted by F. nucleatum experience
  tiff(file = "figures/Intra_personal_dist_Fn_exp-Condition.tiff", width = 1100, height = 600, units = "px", res = 200) ; dist.plot(dist.df, diff.Week.abs, Distance, Fn.experience, col2) + 
    facet_wrap(~Condition) +  theme(strip.text = element_text(size=rel(1.3)), legend.position = "none"); dev.off() 
  
  # Differences in Human read fraction by temporal distance 
  tiff(file = "figures/Human_read_fraction_over_time_by_Fn_exp.tiff", width = 1100, height = 600, units = "px", res = 200) ; dist.plot(dist.df, diff.Week.abs, dist.df$diff.HumanFraction.abs, Fn.experience, col2) + 
    facet_wrap(~Condition) + theme(strip.text = element_text(size=rel(1.3)), legend.position = "none"); dev.off() 



#### 7-2) Intra-persoanl dissimilarity: fixed initial point ####
# Previous trial at 7-1 only considered absolute temporal distance
# Here, we tried to see the directionality toward F. nucleatum obsevation point
# As a control in F. nucleatum-innocent subject,  Random initinal point was selected.

# Participants 
pat.exp <- as.character(unique(meta.order$Participant.ID[meta.order$Fn.experience=="Exp"])) # 20 ; F. nucleatum-experienced subjects
pat.inexp <- as.character(unique(meta.order$Participant.ID[meta.order$Fn.experience=="Non"])) # 86 ; F. nucleatum-innocent subjects

# Setting
set.seed(12345)
Inter.exp <- as.data.frame(matrix(nrow = 10000, ncol = 4)) # Container

# Remove M2083; because the subject have experienced F. nucleatum, but the actual F. nucleatum-detected point was excluded during quality control step due to shortage of total species numbers
# We have used the samples from the M2083 subject for F. nucleatum-expereinced cases 
pat.exp2 <- pat.exp[pat.exp!="M2083"] 


for (i in 1:10000){
  # Select one participants respectively
  tmp.exp <- sample(pat.exp2, 1, replace = F) # F. nucleatum-experienced subjects
  tmp.non <- sample(pat.inexp, 1, replace = F) # F. nucleatum-innocent subjects
  
  # Select two samples from a given F. nucleatum-experiencd subject
  exp1 <- as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.exp[1] & meta.order$Fn.directionality.strict==0]) # Fixed F. nucleatum-detected point 
  exp2 <- sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.exp & meta.order$Fn.directionality.strict!=0]),1) # Randomly selected another samples from the same subject
  
  # Select two random samples from a randomly chosen F. nucleatum-innocent subject
  non.pair <- sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.non[1]]),2, replace = F)
  
  # Container
  Inter.exp[i,] <- data.frame(d.log[exp1,exp2], 
                              d.log[non.pair[1],non.pair[2]],
                              meta.order$Week[meta.order$External.ID==exp1]-
                                meta.order$Week[meta.order$External.ID==exp2],
                              meta.order$Week[meta.order$External.ID==non.pair[1]]-
                                meta.order$Week[meta.order$External.ID==non.pair[2]])
} ; head(Inter.exp, n=30) ; colnames(Inter.exp) <- c("Distance", "Distance", "Week", "Week")



int.exp.tmp <- Inter.exp[,c(1,3)] ; int.exp.tmp$Condition <- "Exp"
int.non.tmp <- Inter.exp[,c(2,4)] ; int.non.tmp$Condition <- "Non"
int.combined <- rbind(int.exp.tmp, int.non.tmp) ; head(int.combined)
int.combined$Condition <- factor(int.combined$Condition, levels=c("Non", "Exp"))

# F. nucleatum-experienced
int_w20 <- int.exp.tmp[abs(int.exp.tmp$Week)<20, ]; dim(int_w20) # 1748 x 3
int_m20 <- int.exp.tmp[abs(int.exp.tmp$Week)>20 , ]; dim(int_m20) # 1115 x 3

int_w20$Class <- ifelse(int_w20$Week>0, yes = "Post", no = "Pre")
int_m20$Class <- ifelse(int_m20$Week>0, yes = "Post", no = "Pre")

wilcox.test(int_w20$Distance~int_w20$Class, alternative="two.sided")
wilcox.test(int_m20$Distance~int_m20$Class, alternative="two.sided")

# F. nucleatum-innocent
int_w20 <- int.non.tmp[abs(int.non.tmp$Week)<20, ]; dim(int_w20) # 1748 x 3
int_m20 <- int.non.tmp[abs(int.non.tmp$Week)>20 , ]; dim(int_m20) # 1115 x 3

int_w20$Class <- ifelse(int_w20$Week>0, yes = "Post", no = "Pre")
int_m20$Class <- ifelse(int_m20$Week>0, yes = "Post", no = "Pre")

wilcox.test(int_w20$Distance~int_w20$Class, alternative="less")
wilcox.test(int_m20$Distance~int_m20$Class, alternative="less")


tiff(file = "figures/Intra_personal_distance_Fixed.tiff", width = 650, height = 600, units = "px", res = 200) ; ggplot(int.combined, aes(Week, Distance, color=Condition)) + 
  xlab("") + ylab("") +
  stat_smooth(geom = "smooth", method="loess", span=0.8, se = T, alpha=0.2) + 
  geom_vline(xintercept = 0, color="gray75", alpha=0.3, size=1, linetype="dotted") +
  scale_fill_manual(values = col2) +  scale_color_manual(values = col2) +
  scale_x_continuous(limits = c(-40, 40)) +
  theme_classic() + 
  theme(legend.position = "none",
        strip.text = element_text(size=rel(1.3))); dev.off()


tiff(file = "figures/Intra_personal_distance_Fixed_box.tiff", width = 220, height = 580, units = "px", res = 200) ; ggplot(int.combined, aes(Condition, Distance)) + 
  geom_sina(aes(color=Condition), alpha=0.01) + 
  geom_boxplot(aes(color=Condition), alpha=0.6, outlier.alpha = 0) + 
  theme_classic() + xlab("") + ylab("") +
  scale_fill_manual(values = col2) +  scale_color_manual(values = col2) +
  stat_compare_means(comparisons = comp_exp, method = "wilcox.test", label = "p.format", tip.length = 0, label.y = 0.95, size=rel(3)) +
  theme(legend.position = "none", 
        axis.text = element_blank(), axis.ticks = element_blank()) ; dev.off()

#### 7-3) Inter-personal distance: Random sampling & distance analysis   ####

# 1) by Condition (IBD vs non-IBD)---------------------

# Dividing "distance matrix" by condition
IBDm <- d.log[meta.order$Condition=="IBD", meta.order$Condition=="IBD"] # 1119 x 1119
nonIBDm <- d.log[meta.order$Condition=="nonIBD", meta.order$Condition=="nonIBD"] # 407 x 407

# Pairwise conversion by melt()
# remove value=0
melt.IBD <- subset(melt(IBDm), value!=0) ; melt.non <- subset(melt(nonIBDm), value!=0) ;rm(melt.IBD, melt.non, IBDm, nonIBDm)

# Drawing function
inter.plot <- function(data, distance, variable, color) {
  attach(data)
  ggplot(data, aes(distance)) + 
    ylab("") + xlab("") + 
    geom_density(alpha=0.4, aes(fill=variable), color=NA) + 
    theme_classic() + 
    scale_fill_manual(values = color)+
    theme(axis.title = element_blank(),
          axis.text.x = element_text(size=rel(1), color="black"), 
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.text = element_text(size=rel(1.1)), 
          legend.title = element_blank(),
          legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_blank(),
          legend.key = element_rect(color="transparent", fill="transparent"),
          legend.key.width = unit(1, "line"),
          legend.key.height = unit(1,"line"), 
          legend.spacing.y = unit(0, "cm"), 
          legend.spacing.x = unit(1, "mm"))
}
box.plot <- function(data, x, y, color) {
  attach(data)
  ggplot(data, aes(x = x, y=y, fill=x)) + 
    xlab("") + ylab("")+
    geom_boxplot(alpha=0.5, outlier.alpha = 0.02) + 
    scale_fill_manual(values = color)+
    #stat_compare_means(comparisons = comparisons, size=1.5, method = "wilcox.test", label="p.format", tip.length = 0, label.y = label.y)+
    # We are going to draw manually 
    theme_void() + 
    coord_flip() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
}

### Procedure ###
# Step.1 Select two participants
  # Step.1-1) Out of same category
  # Step.1-2) Out of different category (control)
# Step.2 Setect one sample for each participants
# Step.3 Calculate distance between the two separately selected samples
# Step.4 Calculate proximity difference/shannon diversity difference/human read fraction difference

# Step.0) Assigning participants into each category 
# Condition
pat.ibd <- as.character(unique(meta.order$Participant.ID[meta.order$Condition=="IBD"])) # 80 patients
pat.non <- as.character(unique(meta.order$Participant.ID[meta.order$Condition=="nonIBD"])) # 26 patients

# 1) by Condition
set.seed(1234)
Inter.condition <- as.data.frame(matrix(nrow = 1000, ncol = 3)) # Condtainer

for (i in 1:1000){
  # Select two participants
  tmp.ibd <- sample(pat.ibd, 2, replace = F) 
  tmp.non <- sample(pat.non, 2, replace = F)
  tmp.mix <- sample(unique(pat), 2, replace = F)
  
  # Bring distance from pre-calculated distance matrix
  # Randomly select one sample from the previously selected subjects
  Inter.condition[i,] <- data.frame(d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.ibd[1]]),1), 
                                          sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.ibd[2]]),1)], # IBD 
                                    d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.non[1]]),1), 
                                          sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.non[2]]),1)], # nonIBD
                                    d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.mix[1]]),1), 
                                          sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.mix[2]]),1)]) # Random
} ; colnames(Inter.condition) <- c("IBD", "nonIBD", "Random")

int.con.melt <- melt(Inter.condition) ; colnames(int.con.melt) <- c("Condition", "Distance")
int.con.melt$Condition <- factor(int.con.melt$Condition, levels=c("Random", "nonIBD", "IBD"))

#tiff(file = "figures/Inter_personal_distance_by_Condition_density only.tiff", width = 800, height = 600, units = "px", res = 200) ; inter.plot(int.con.melt, Distance, Condition, col3) ; dev.off()
#tiff(file = "figures/Inter_personal_distance_by_Condition_box only.tiff", width = 650, height = 800, units = "px", res = 200) ; box.plot(int.con.melt, Condition, Distance, col3) +
# stat_compare_means(comparisons = comp_con_3, method = "wilcox.test", tip.length = 0, label.y=c(0.6, 0.8, 1.0)) ; dev.off()

# Combining into one graph
g1 <- inter.plot(int.con.melt, Distance, Condition, col3) ; g1 <- ggplot_gtable(ggplot_build(g1)) 
g2 <- box.plot(int.con.melt, Condition, Distance, col3) ; g2 <- ggplot_gtable(ggplot_build(g2))

# Get maximum widths and heights
gWidth <- unit.pmax(g1$widths[2:3], g2$widths[2:3])
gHeight <- unit.pmax(g1$heights[4:5], g2$heights[4:5])

# Set the maximums in the gtables for g1 and g2
g1$widths[2:3] <- as.list(gWidth) ; g2$widths[2:3] <- as.list(gWidth)
g1$heights[4:5] <- as.list(gHeight) ; g2$heights[4:5] <- as.list(gHeight)

# Creat a new gtable  & insert gt1, gt2, gt3 into the gtable
gt <- gtable(widths = unit(c(100,1), "null"), heights = unit(c(1,5), "null"))
gt <- gtable_add_grob(gt, g1, 2, 1) ; gt <- gtable_add_grob(gt, g2, 1, 1)

# And render the plot
grid.newpage() ; tiff(file = "figures/Inter_personal_distance_by_Condition.tiff", width = 800, height = 800, units = "px", res = 200) ; grid.draw(gt) ; dev.off()

  # Median dissimilarity
  median(int.con.melt$Distance[int.con.melt$Condition=="IBD"]) # 0.5696
  median(int.con.melt$Distance[int.con.melt$Condition=="nonIBD"]) # 0.5000
  median(int.con.melt$Distance[int.con.melt$Condition=="Random"]) # 0.5586

# 2) F.n-experience
# Fn.exposure
pat.exp <- as.character(unique(meta.order$Participant.ID[meta.order$Fn.experience=="Exp"])) # 20 patients
pat.inexp <- as.character(unique(meta.order$Participant.ID[meta.order$Fn.experience=="Non"])) # 86 patients

set.seed(1234)
Inter.exp <- as.data.frame(matrix(nrow = 1000, ncol = 3)); for (i in 1:1000){
  tmp.exp <- sample(pat.exp, 2, replace = F) # Select two participants
  tmp.non <- sample(pat.inexp, 2, replace = F)
  tmp.mix <- sample(unique(pat), 2, replace = F)
  
  Inter.exp[i,] <- data.frame(d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.exp[1]]),1), 
                                    sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.exp[2]]),1)],
                              d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.non[1]]),1), 
                                    sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.non[2]]),1)],
                              d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.mix[1]]),1), 
                                    sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.mix[2]]),1)])
} ; colnames(Inter.exp) <- c("Exp", "Non", "Random")
int.exp.melt <- melt(Inter.exp) ; colnames(int.exp.melt) <- c("Condition", "Distance")
int.exp.melt$Condition <- factor(int.exp.melt$Condition, levels=c("Random", "Non", "Exp"))

#tiff(file = "figures/inter_distance.fn_exp.tiff", width = 800, height = 600, units = "px", res = 200) ; inter.plot(int.exp.melt, Distance, Condition, col3) ; dev.off()
#tiff(file = "figures/inter_distance.box.fn_exp.tiff", width = 550, height = 800, units = "px", res = 200) ; box.plot(int.exp.melt, Condition, Distance, col3) +
#  stat_compare_means(comparisons = list(c("Non", "Exp"), c("Non", "Random"), c("Exp", "Random")), tip.length = 0, method = "wilcox.test", label.y = c(0.3, 0.5, 0.9)); dev.off()

# Combining into one graph
g1 <- inter.plot(int.exp.melt, Distance, Condition, col3) ; g1 <- ggplot_gtable(ggplot_build(g1))
g2 <- box.plot(int.exp.melt, Condition, Distance, col3)  ; g2 <- ggplot_gtable(ggplot_build(g2))

# Get maximum widths and heights
gWidth <- unit.pmax(g1$widths[2:3], g2$widths[2:3])
gHeight <- unit.pmax(g1$heights[4:5], g2$heights[4:5])

# Set the maximums in the gtables for g1 and g2
g1$widths[2:3] <- as.list(gWidth) ; g2$widths[2:3] <- as.list(gWidth)
g1$heights[4:5] <- as.list(gHeight) ; g2$heights[4:5] <- as.list(gHeight)

# Creat a new gtable & insert gt1, gt2, gt3 into the new gtable
gt <- gtable(widths = unit(c(100,1), "null"), heights = unit(c(1,5), "null"))
gt <- gtable_add_grob(gt, g1, 2, 1) ; gt <- gtable_add_grob(gt, g2, 1, 1)

# And render the plot
grid.newpage() ; tiff(file = "figures/Inter_personal_distance_by_Fn_exp.tiff", width = 800, height = 800, units = "px", res = 200) ; grid.draw(gt) ; dev.off()

  # Median dissimilarity
  median(int.exp.melt$Distance[int.exp.melt$Condition=="Exp"]) # 0.5927
  median(int.exp.melt$Distance[int.exp.melt$Condition=="Non"]) # 0.5401
  median(int.exp.melt$Distance[int.exp.melt$Condition=="Random"]) # 0.5586


# 3) F. nucleatum prior vs. posterior---------------
# Samples were divided into two groups: F. nucleatum-prior, F. nucleatum-posterior (using strict criteria)
meta.non <- meta.order[meta.order$Fn.pre.post.simple=="Negative",]
meta.pre <- meta.order[meta.order$Fn.pre.post.simple=="Prior",]
meta.post <- meta.order[meta.order$Fn.pre.post.simple=="Posterior",]
meta.mix <- rbind(meta.non, meta.pre, meta.post)

set.seed(1234)
inter.prepost <- as.data.frame(matrix(nrow = 1000, ncol = 4)) ; for (i in 1:1000){
  tmp.non <- sample(unique(as.character(meta.non$Participant.ID)), 2, replace = F)
  tmp.pre <- sample(unique(as.character(meta.pre$Participant.ID)), 2, replace = F)
  tmp.post <- sample(unique(as.character(meta.post$Participant.ID)), 2, replace = F)
  tmp.mix <- sample(unique(as.character(meta.mix$Participant.ID)), 2, replace = F)
  
  inter.prepost[i,] <- data.frame(d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.non[1]]),1), 
                                        sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.non[2]]),1)],
                                  d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.pre[1]]),1), 
                                        sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.pre[2]]),1)],
                                  d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.post[1]]),1), 
                                        sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.post[2]]),1)],
                                  d.log[sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.mix[1]]),1), 
                                        sample(as.character(meta.order$External.ID[meta.order$Participant.ID==tmp.mix[2]]),1)])
} ; colnames(inter.prepost) <- c("Negative", "Prior", "Posterior", "Random")
int.pp.melt <- melt(inter.prepost) ; colnames(int.pp.melt) <- c("Condition", "Distance")
int.pp.melt$Condition <- factor(int.pp.melt$Condition, levels=c("Random", "Negative",  "Prior", "Posterior"))

#tiff(file = "figures/inter_distance.pre_post_wNegative.tiff", width = 800, height = 600, units = "px", res = 200) ; #inter.plot(int.pp.melt, Distance, Condition, col4.3) ; dev.off()
#tiff(file = "figures/inter_distance.box.pre_post_wNegative.tiff", width = 850, height = 800, units = "px", res = 200) ; box.plot(int.pp.melt, Condition, Distance, col4.3) +
#  stat_compare_means(comparisons = list(c("Posterior", "Prior"), c("Posterior", "Negative"), c("Posterior", "Random")), method = "wilcox.test", tip.length = 0, label.y = c(0.3, 0.5, 0.9)); dev.off()

# Combining into one graph
g1 <- inter.plot(int.pp.melt, Distance, Condition, col4.3)
g2 <- box.plot(int.pp.melt, Condition, Distance, col4.3)

g1 <- ggplot_gtable(ggplot_build(g1))
g2 <- ggplot_gtable(ggplot_build(g2))

# Get maximum widths and heights
gWidth <- unit.pmax(g1$widths[2:3], g2$widths[2:3])
gHeight <- unit.pmax(g1$heights[4:5], g2$heights[4:5])

# Set the maximums in the gtables for g1 and g2
g1$widths[2:3] <- as.list(gWidth) ; g2$widths[2:3] <- as.list(gWidth)
g1$heights[4:5] <- as.list(gHeight) ; g2$heights[4:5] <- as.list(gHeight)

# Creat a new gtable  & insert gt1, gt2, gt3 into the new gtable
gt <- gtable(widths = unit(c(100,1), "null"), heights = unit(c(1,5), "null"))
gt <- gtable_add_grob(gt, g1, 2, 1) ; gt <- gtable_add_grob(gt, g2, 1, 1)

# And render the plot
grid.newpage() ; tiff(file = "figures/Inter_personal_distance_by_Fn_pre_post.tiff", width = 800, height = 800, units = "px", res = 200) ; grid.draw(gt) ; dev.off()

  # Median dissimilarity
  median(int.pp.melt$Distance[int.pp.melt$Condition=="Posterior"]) # 0.5816
  median(int.pp.melt$Distance[int.pp.melt$Condition=="Prior"]) # 0.5372
  median(int.pp.melt$Distance[int.pp.melt$Condition=="Negative"]) # 0.5454
  median(int.pp.melt$Distance[int.pp.melt$Condition=="Random"]) # 0.5605



#### 7-4) Inter-personal distance: F.nucleatum-detected point to other points from different subjects ####
# F. nucleatum-detected samples as a fixed initial point 

pat.exp <- as.character(unique(meta.order$Participant.ID[meta.order$Fn.experience=="Exp"])) # 20 subjects
pat.exp2 <- pat.exp[pat.exp!="M2083"] # It does not contain F.nucleatum-detected point
pat.non <- as.character(unique(meta.order$Participant.ID[meta.order$Fn.experience=="Non"]))

# Initial point for F.nucleatum-exposed subjects
Fn.point <- as.character(na.omit(meta.order$External.ID[meta.order$Fn.directionality.strict==0]))

Inter.exp.fn <- as.data.frame(matrix(nrow = 0, ncol = 2)) ; for (i in 1:length(Fn.point)){
  tmp.initial <- Fn.point[i]
  tmp.subject <- as.character(unique(meta.order$Participant.ID[meta.order$External.ID==tmp.initial]))
  tmp.wo_subject <- pat.exp2[pat.exp2!=tmp.subject]
  tmp.target <- as.character(meta.order$External.ID[meta.order$Participant.ID %in% tmp.wo_subject])
  
  tmp.distance <- d.log[tmp.initial, tmp.target]
  tmp.diff.week <- meta.order$Fn.directionality.strict[meta.order$External.ID %in% tmp.target]
  
  tmp.df <- data.frame(tmp.distance, 
                       tmp.diff.week)
  rownames(tmp.df) <- NULL
  Inter.exp.fn <- rbind(Inter.exp.fn, tmp.df)
  rm(tmp.df)
}

Inter.exp.fn$Class <- "Exp" ; colnames(Inter.exp.fn) <- c("Distance", "Diff.Week", "Class")


set.seed(12345)
Inter.non <- as.data.frame(matrix(nrow = 0, ncol = 2))
innocent.samples <- as.character(meta.order$External.ID[meta.order$Fn.experience=="Non"]) # 1209 samples

for (i in 1:10000){
  tmp.subject <- sample(pat.non, 2, replace=F)
  tmp.sample1 <- sample(meta.order$External.ID[meta.order$Participant.ID==tmp.subject[1]], 1)
  tmp.sample2 <- sample(meta.order$External.ID[meta.order$Participant.ID==tmp.subject[2]], 1)
  
  tmp.distance <- d.log[tmp.sample1, tmp.sample2]
  tmp.diff.week <- meta.order$Week[meta.order$External.ID==tmp.sample1]-meta.order$Week[meta.order$External.ID==tmp.sample2]
  
  tmp.df <- data.frame(tmp.distance, 
                       tmp.diff.week)
  rownames(tmp.df) <- NULL
  Inter.non <- rbind(Inter.non, tmp.df)
  rm(tmp.df)
} ; Inter.non$Class <- "Non" ; colnames(Inter.non) <- c("Distance", "Diff.Week", "Class")

# Merge
Inter.test <- rbind(Inter.non, Inter.exp.fn)

tiff(file = "figures/Inter_personal_distance_fixed_initial_point.tiff", width = 700, height = 600, units = "px", res = 200) ; ggplot(Inter.test, aes(Diff.Week, Distance, group=Class, color=Class)) +
  xlab("") + ylab("") +
  stat_smooth(geom = "smooth", method = "loess", span = 0.8, size=1, alpha=0.3) +
  geom_vline(xintercept = 0, color="gray75", alpha=0.3, size=1.2) +
  scale_color_manual(values = rev(col2)) +  scale_x_continuous(limits = c(-40, 40)) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black", size=rel(1)),
        axis.text.y = element_text(color="black", size=rel(0.8)),
        legend.position = "none") ; dev.off()
        #legend.text = element_text(size=rel(0.9)), 
        #legend.title = element_blank(),
        #legend.position = c(0.75, 1),
        #legend.justification = c(0,1),
        #legend.background = element_rect(color = "black", size=0.5),
        #legend.key = element_rect(fill = "white"),
        #legend.key.width = unit(0.8, "line"),
        #legend.key.height = unit(0.8,"line"), 
        #legend.spacing.y = unit(0, "cm"), 
        #legend.spacing.x = unit(0, "cm")



#### 8. F. nucleatum-Prior/Posterior analysis ####

# We are going to test the distributional dynamics for all microbes,
# based on F. nucleatum-experience, F. nucleatum-prior/posterior, 

# Pre-screened F.nucleatum-prior/posterior enriched microbes
# NOD >=5 in F. nucleautm-experienced subjects
#mic_pp <- read.csv("input/Fnucleatum_pre_post_Pvalue.csv", header=T, row.names=1)
head(mic_pp)

# Select significant microbes
mic_sig <- subset(mic_pp, Pre.Pvalue<0.01 | Post.Pvalue<0.01)
mic_sig$enriched <- "Posterior" ; mic_sig$enriched[mic_sig$Pre.Pvalue<0.01] <- "Prior"
mic_sig$enriched <- factor(mic_sig$enriched, levels=c("Prior", "Posterior"))
mic_sig$NOD <- apply(t(bugs)[, rownames(mic_sig)]!=0, 2, sum)
  summary(mic_sig$enriched) # Prior 43 ; Posterior 36
  
bugs.prior <- rownames(mic_sig)[mic_sig$enriched=="Prior"] ; bugs.posterior <- rownames(mic_sig)[mic_sig$enriched=="Posterior"] 


#### 8-1) Testing distribution of the prior/posterior markers in whole samples ####
# Boxplot 
dis.plot  <- function(data, x, y, color, comparisons, label.y) {
  attach(data)
  ggplot(data, aes(x = x, y = y, fill=x)) + xlab("") + ylab("")+
    geom_boxplot(alpha=0.5, outlier.alpha = 0.02) + 
    scale_fill_manual(values = color) + 
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format", bracket.size = 0.5, tip.length = 0, label.y=label.y)+
    theme_classic() + 
    theme(axis.text.x = element_text(size=rel(1.3), color="black"), 
          axis.text.y = element_text(size=rel(0.9), color="black"),
          legend.position = "none")
}

# 1) Species ditribution 
bugs.tmp1 <- tbugs.order[, bugs.prior] ; bugs.tmp2 <- tbugs.order[, bugs.posterior]
bugs.pp <- cbind(bugs.tmp1, bugs.tmp2) ; bugs.pp <- data.frame(bugs.pp, meta.order) ; rm(bugs.tmp1, bugs.tmp2)
  identical(rownames(bugs.pp), as.character(bugs.pp$External.ID)) # TRUE

# Sum of abundance
bugs.pp$sum.post <- apply(bugs.pp[,bugs.posterior], 1, sum)
bugs.pp$sum.pre <- apply(bugs.pp[,bugs.prior], 1, sum)

# Number of detected species 
bugs.pp$species.post <- apply(bugs.pp[,bugs.posterior]!=0, 1, sum)
bugs.pp$species.pre <- apply(bugs.pp[,bugs.prior]!=0, 1, sum)

# Ratiometric
bugs.pp$ratio_sumPP <- log10((bugs.pp$sum.post+0.00001)/(bugs.pp$sum.pre+0.00001))
bugs.pp$ratio_specPP <- (bugs.pp$species.post+1)/(bugs.pp$species.pre+1)
bugs.pp.non <- bugs.pp[bugs.pp$Fn.experience=="Non",]


# Disttribution of Prior/Posterior markers in F. nucleatum-innocent subjects

# 8-1) Prior-enriched microbes
# Condition
tiff(file = "figures/Prior_Marker_sum_condition.tiff", width = 300, height = 700, units = "px", res = 200) ; dis.plot(bugs.pp.non, Condition, sum.pre, col2, comp_con, 97) ; dev.off()  
tiff(file = "figures/Prior_Marker_spec_condition.tiff", width = 300, height = 700, units = "px", res = 200) ; dis.plot(bugs.pp.non, Condition, species.pre, col2, comp_con, 32) ; dev.off()

# 8-2) Posterior-enriched microbes
# Condition
tiff(file = "figures/Posterior_Marker_sum_condition.tiff", width = 300, height = 700, units = "px", res = 200) ; dis.plot(bugs.pp.non, Condition, sum.post, col2, comp_con, 43) ; dev.off()  
tiff(file = "figures/Posterior_Marker_spec_condition.tiff", width = 300, height = 700, units = "px", res = 200) ; dis.plot(bugs.pp.non, Condition, species.post, col2, comp_con, 17) ; dev.off()


# 8-3) Ratiometric
# Condition
tiff(file = "figures/Ratio_pre_post_Marker_sum_condition.tiff", width = 300, height = 700, units = "px", res = 200) ; dis.plot(bugs.pp.non, Condition, ratio_sumPP, col2, comp_con, 2.2) ; dev.off()
tiff(file = "figures/Ratio_pre_post_Marker_spec_condition.tiff", width = 300, height = 700, units = "px", res = 200) ; dis.plot(bugs.pp.non, Condition, ratio_specPP, col2, comp_con, 2.5) ; dev.off()  



#### 8-2) Screening Classifiers ####

# Screening "Classifiers" that differentiate the appearing point of F. nucleatum 

# input data to use
dim(tb_pp) # 317 samples with 258 species (NOD>=5) + 1 metadata (Fn.pre.post.simple)
  head(NOD5_Fn) ; num5 <- sum(NOD5_Fn==TRUE) # 258 species: detected at least 5 times in 317 samples

  
# Partitioning (1000)
set.seed(12345)
idx <- createDataPartition(tb_pp$Fn.pre.post.simple, times = 1000, p=0.7, list=F) # 222 samples for 1000 trials 
  # n-th trial : n=1000 --> 1000 partitions
  # i-th microbes : i upto 258

# Container
bugs.auc.stat <- as.data.frame(matrix(nrow = 258, ncol = 4)) 
  bugs.names <- colnames(tb_pp)[-259]
  bugs.auc.stat$microbe <- bugs.names 
  bugs.auc.stat <- bugs.auc.stat[,c(5,1:4)]

auc.tmp <- as.data.frame(matrix(nrow = 1000, ncol = 1)) ; auc.tmp$Index <- paste0("#", 1:1000) ; colnames(auc.tmp) <- c("AUC", "Index")
roc.tmp <- list() ; roc.list <- list()

for (i in 1:258){
  # bugs.name[i] = bring i-th microbes names (temporary)
  bugs.tmp <- bugs.names[i] 
  
  for (n in 1:1000){
    # idx[,n] = n-th partition order,
    # bring the partitions from input data (tb_pp)
    pp.tmp <- tb_pp[idx[,n],]  
    auc.tmp[n,1] <- roc(pp.tmp$Fn.pre.post.simple, pp.tmp[,bugs.tmp])$auc # get AUC for each microbes on n-th partition
    roc.tmp[[n]] <- roc(pp.tmp$Fn.pre.post.simple, pp.tmp[,bugs.tmp])
  }
  
  bugs.auc.stat[i,2] <- mean(auc.tmp$AUC)
  bugs.auc.stat[i,3] <- sd(auc.tmp$AUC)
  
  # Calculating p-value if the AUCs are really higher than 0.5
  bugs.auc.stat[i,4] <- pnorm(mean=bugs.auc.stat[i,2], sd=bugs.auc.stat[i,3], q=0.5, lower.tail=T)
  bugs.auc.stat[i,5] <- p.adjust(bugs.auc.stat[i,4], method="fdr", n=258)
  
  roc.list[[i]] <- roc.tmp
} ; bugs.auc.stat$NOD <- NOD[names(NOD) %in% bugs.auc.stat$microbe] ; colnames(bugs.auc.stat) <- c("microbe", "meanAUC", "stdAUC", "p.value", "fdr", "NOD")
    bugs.auc.sig <- bugs.auc.stat[bugs.auc.stat$fdr<0.001,] # 41 species

#### 8-3) Classifier analysis #### 
# Log transformation of p-values
bugs.auc.stat$logFDR <- -log10(bugs.auc.stat$fdr) ; bugs.auc.stat$logP <- -log10(bugs.auc.stat$p.value)

# Ordering by mean AUC
bugs.auc.stat <- bugs.auc.stat[order(bugs.auc.stat$meanAUC, decreasing = T),] ; rownames(bugs.auc.stat) <- NULL ; bugs.ordered_mAUC <- as.character(bugs.auc.stat$microbe)

# Prior vs. Posterior-enriched microbes
# p-value < 0.01
bugs.auc.stat$pre.post <- "Neither" ; bugs.auc.stat$pre.post[bugs.auc.stat$microbe %in% bugs.posterior] <- "Posterior" ; bugs.auc.stat$pre.post[bugs.auc.stat$microbe %in% bugs.prior] <- "Prior"
  bugs.auc.stat$pre.post <- factor(bugs.auc.stat$pre.post, levels = c("Neither", "Prior", "Posterior"))

# nonIBD vs. IBD markers
bugs.ibd <- rownames(mark)[mark$Class=="IBD"] ; bugs.non <- rownames(mark)[mark$Class=="nonIBD"]
bugs.auc.stat$Condition <- "Neither" ; bugs.auc.stat$Condition[bugs.auc.stat$microbe %in% bugs.ibd] <- "IBD" ; bugs.auc.stat$Condition[bugs.auc.stat$microbe %in% bugs.non] <- "nonIBD"
  bugs.auc.stat$Condition <- factor(bugs.auc.stat$Condition, levels = c("Neither", "nonIBD", "IBD"))

# Add core signature microbes for CRC
bugs.auc.stat$Core <- ifelse(bugs.auc.stat$microbe %in% core_CRC, yes = "Yes", no = "No") ; bugs.auc.stat$Significance <- ifelse(bugs.auc.stat$fdr <0.001 , yes = "Sig", no = "NS")
  #write.csv(bugs.auc.stat, "input/bugs.auc.stat.csv")
  #bugs.auc.stat <- read.csv("input/bugs.auc.stat.csv", header = T, row.names = 1) ; dim(bugs.auc.stat) # 258 x 12

# Distribution of Classifiers in IBD/nonIBD condition
class_post = bugs.auc.stat$microbe[bugs.auc.stat$Significance=="Sig" & bugs.auc.stat$pre.post=="Posterior"]  # 26
class_pre = bugs.auc.stat$microbe[bugs.auc.stat$Significance=="Sig" & bugs.auc.stat$pre.post=="Prior"]  # 15

tb_con$post_spec <- rowSums(tb_con[,class_post]>0) ; tb_con$pre_spec <- rowSums(tb_con[,class_pre]>0)
tb_con$post_sum <- rowSums(tb_con[,class_post]) ; tb_con$pre_sum <- rowSums(tb_con[,class_pre])

melt_con_spec <- melt(tb_con, id.vars = "Condition", measure.vars = c("post_spec", "pre_spec"))
melt_con_sum <- melt(tb_con, id.vars = "Condition", measure.vars = c("post_sum", "pre_sum"))

tiff(file = "figures/classifier_condition_spec.tiff", width = 1000, height = 500, units = "px", res = 200) ; ggplot(melt_con_spec, aes(Condition, value, fill=Condition)) + 
  xlab("")+ylab("")+
  geom_boxplot(alpha=0.7, outlier.alpha = 0.05) +
  scale_fill_manual(values=col2) + 
  coord_flip() +
  stat_compare_means(method = "wilcox.test", comparisons = comp_con, label = "p.signif", label.y= 16) +
  facet_wrap(~variable, ncol=1) + theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size=rel(1.2), color="black"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()); dev.off()



# CRC core signature ~ average AUC of microbes 
# Ordering by p-value (Previosuly, ordered by meanAUC)
mic_pp$microbe <- rownames(mic_pp)
bugs.auc.stat <- merge(bugs.auc.stat, mic_pp, by = "microbe")
bugs.auc.stat <- bugs.auc.stat[order(bugs.auc.stat$logP, decreasing = T),]  ; rownames(bugs.auc.stat) <- NULL
  bugs.auc.stat$Core <- factor(bugs.auc.stat$Core, levels = c("Yes", "No"))
  bugs.list <- bugs.auc.stat$microbe # To get order of microbes (p-value)

  
# Core species ~ Classifying p-values
tiff(file = "figures/Core_CRC_vs_classifying_logP.tiff",  width = 1200, height = 800, units = "px", res = 200) ; ggplot(bugs.auc.stat, aes(color=Core)) + 
  geom_point(aes(x=microbe, y = logP, alpha=Core), size=2.5) +
  geom_hline(yintercept = -log10(0.05), color="gray20", alpha=0.2, size=1, linetype="dotted")+
  xlab("") + ylab("") + 
  scale_alpha_manual(values=c(0.7, 0.3)) +  scale_x_discrete(limits = bugs.list) +  scale_color_manual(values = c("darkred", "gray84"))+
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color="black"), 
        axis.title = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        legend.position = "none") ; dev.off()


# Core species ~ Posterior-enrichment p-values
tiff(file = "figures/Core_CRC_vs_Posterior_logP.tiff",  width = 1200, height = 800, units = "px", res = 200) ; ggplot(bugs.auc.stat, aes(color=Core)) + 
  geom_point(aes(x=microbe, y = -log10(Diff.Pvalue), alpha=Core), size=2.5) +
  geom_hline(yintercept = -log10(0.05), color="gray20", alpha=0.2, size=1, linetype="dotted")+
  xlab("") + ylab("") + 
  scale_alpha_manual(values=c(0.7, 0.3)) +  scale_x_discrete(limits = bugs.list) +  scale_color_manual(values = c("darkred", "gray84"))+
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color="black"), 
        axis.title = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        legend.position = "none") ; dev.off()

# Enrichment figure
tile_draw <- function(data, fill, limits){
  attach(data)
  ggplot(data) +
    geom_tile(aes(x = microbe, y=1, fill=fill), color="white") +
    geom_rect(xmin=0.7, xmax=258.5, ymin=0.5, ymax=1.5, fill="transparent", color="black", size=0.4)+
    scale_x_discrete(limits=limits)+
    theme_void() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")
}

tiff(file = "figures/Enrichment_logP_Core.tiff", res = 300, width = 1000, height = 80, units="px") ; tile_draw(bugs.auc.stat, Core, bugs.ordered_mAUC) + scale_fill_manual(values = c("white", "black")); dev.off()
tiff(file = "figures/Enrichment_logP_IBD-marker.tiff", res = 300, width = 1000, height = 80, units="px") ; tile_draw(bugs.auc.stat, Condition, bugs.ordered_mAUC) + scale_fill_manual(values = c("white", "white", "black")) ; dev.off()
tiff(file = "figures/Enrichment_logP_nonIBD-marker.tiff", res = 300, width = 1000, height = 80, units="px") ; tile_draw(bugs.auc.stat, Condition, bugs.ordered_mAUC) + scale_fill_manual(values = c("white", "black", "white")) ; dev.off()

# Calculating enrichment p-value 
  # Every CRC core microbes, detected at least one time in Fn-exposed samples, had lower p-value than 0.05 for classification (meanAUC>0.5)
  # All CRC core microbes had p-value lower than 0.05  

# Wilcoxon test, p-value ranking ~ CRC core signature
wilcox.test(logP~Core, data = bugs.auc.stat, alternative="two.sided") # 0.0009132

# Wilcoxon test, p-value ranking ~ IBD marker
bugs.test <- as.data.frame(matrix(nrow = 258, ncol = 0))
  bugs.test$IBD <- ifelse(bugs.auc.stat$Condition=="IBD", yes="IBD", no="Non") ;  bugs.test$logP <- bugs.auc.stat$logP
  wilcox.test(logP~IBD, data = bugs.test, alternative="two.sided") # 0.02285

# Wilcoxon test, p-value ranking ~ non-IBD marker
bugs.test <- as.data.frame(matrix(nrow = 258, ncol = 0))
  bugs.test$nonIBD <- ifelse(bugs.auc.stat$Condition=="nonIBD", yes="nonIBD", no="Non");  bugs.test$logP <- bugs.auc.stat$logP
  wilcox.test(logP~nonIBD, data = bugs.test, alternative="two.sided") # 0.07318

# Distribution of classifiers
# Merge classifier data with F.nucleatum-prior/posterior-enrichment data
mic_pp$microbe <- rownames(mic_pp) ; dim(mic_pp) # 258 species x 7
bugs.auc.sig <- bugs.auc.stat[bugs.auc.stat$Significance=="Sig",] 
bugs.auc.sig <- merge(bugs.auc.sig, mic_pp, by = "microbe")
  # Categorize the classifiers based on their enrichment at F.nucleatum Prior vs. Posterior 
  bugs.auc.sig$enriched <- ifelse(test = bugs.auc.sig$Pre.Pvalue < 0.01, yes = "Prior", no = "Posterior")
  bugs.auc.sig$enriched <- factor(bugs.auc.sig$enriched, levels=c("Prior", "Posterior"))
  table(bugs.auc.sig$enriched) # Classifiers= Posterior 26 ; Prior 15
  table(bugs.auc.sig$pre.post) # Classifiers= Posterior 26 ; Prior 15

# Pie chart
pie.draw <- function(data, y, color) {
  attach(data)
  ggplot(data, aes(x=1, y = y, fill=y)) + 
    geom_bar(stat = "identity", alpha=0.8) + 
    xlab("") + ylab("") +
    scale_fill_manual(values = color) + 
    coord_polar(theta="y") + 
    theme_void() + 
    theme(legend.position = "none")
}

# All microbes information
bugs.533 = rownames(bugs) ; bugs.533 <- 
  data.frame(microbe=bugs.533, NOD=rowSums(bugs!=0)) ; rownames(bugs.533) <- NULL

# Add Core signature for CRC
bugs.533$CRC <- ifelse(bugs.533$microbe %in% core_CRC, yes="CRC marker", no="None")
bugs.533$CRC <- factor(bugs.533$CRC, levels=c("CRC marker", "None"))
  table(bugs.533$CRC) # CRC marker 18

# Add IBD/nonIBD markers
bugs.533$Condition <- "Neither"; bugs.533$Condition[bugs.533$microbe %in% ibd.mark] <- "IBD" ; bugs.533$Condition[bugs.533$microbe %in% non.mark] <- "nonIBD"
bugs.533$Condition <- factor(bugs.533$Condition, levels = c("nonIBD", "IBD", "Neither"))
  table(bugs.533$Condition) # nonIBD 14 ; IBD 12 ; Neither 507

# Add Classifiers
bugs.533$Classifier <- "None" ; bugs.533$Classifier[bugs.533$microbe %in%  classifier_sig$microbe] <- "Classifier"
bugs.533$Classifier <- factor(bugs.533$Classifier, levels=c("Classifier", "None"))
  table(bugs.533$Classifier) # Classifier 41

# Add F. nucleatum-prior/posterior enrichment
bugs.533$enriched <- "Neither" ; bugs.533$enriched[bugs.533$microbe %in% bugs.prior] <- "Prior"; bugs.533$enriched[bugs.533$microbe %in% bugs.posterior] <- "Posterior"
bugs.533$enriched <- factor(bugs.533$enriched, levels=c("Prior", "Posterior", "Neither"))
  table(bugs.533$enriched) # Prior 44 ; Posterior 36 ; Neither 453

# Add Classifier + Prior/Posterior
bugs.533$Class <- "None" ; bugs.533$Class[bugs.533$Classifier=="Classifier" & bugs.533$enriched=="Posterior"] <- 
  "Post-Class" ; bugs.533$Class[bugs.533$Classifier=="Classifier" & bugs.533$enriched=="Prior"] <- "Pre-Class"
bugs.533$Class <- factor(bugs.533$Class, levels=c("Pre-Class", "Post-Class", "None"))
  table(bugs.533$Class) # Pre-Class 15 ; Post-Class 26 ;  None 492

tiff("figures/pie_core_CRC_533.tiff", res = 400, width = 600, height = 600) ; pie.draw(bugs.533, CRC, rev(col2.1)) ; dev.off()
tiff("figures/pie_core_CRC_258.tiff", res = 400, width = 600, height = 600) ; pie.draw(bugs.auc.stat, Core, col2.1) ; dev.off()
tiff("figures/pie_core_CRC_in_41_classifiers.tiff", res = 400, width = 600, height = 600) ; pie.draw(bugs.auc.sig, Core, col2.1) ; dev.off()

# Lollipop chart to visualize Classifiers
lol_order <- bugs.auc.sig$microbe[order(bugs.auc.sig$pre.post, bugs.auc.sig$meanAUC, decreasing=T)] ; lol_order_meanAUC  <- bugs.auc.sig$microbe 
rownames(bugs.auc.sig) <- NULL ; dim(bugs.auc.sig) ; head(bugs.auc.sig) # 41 classifiers

# Rename to remove unnecessary "_"
bugs.auc.sig$names <- gsub(bugs.auc.sig$microbe, pattern = "_", replacement = " ")
bugs.auc.sig[bugs.auc.sig$names=="Lachnospiraceae bacterium 3 1 46FAA", "names"] <- "Lachnospiraceae bacterium 3_1_46FAA"
bugs.auc.sig[bugs.auc.sig$names=="Burkholderiales bacterium 1 1 47", "names"] <- "Burkholderiales bacterium 1_1_47"
bugs.auc.sig[bugs.auc.sig$names=="Eubacterium sp 3 1 31", "names"] <- "Eubacterium sp 3_1_31"
bugs.auc.sig[bugs.auc.sig$names=="Erysipelotrichaceae bacterium 5 2 54FAA", "names"] <- "Erysipelotrichaceae bacterium 5_2_54FAA"
bugs.auc.sig[bugs.auc.sig$names=="Bacteroides sp 2 1 22", "names"] <- "Bacteroides sp 2_1_22"

# conversion to Factor
bugs.auc.sig$pre.post<- factor(bugs.auc.sig$pre.post, levels=c("Prior", "Posterior"))
bugs.auc.sig$Condition <- factor(bugs.auc.sig$Condition, levels=c("nonIBD", "IBD", "Neither"))
bugs.auc.sig$Core <- factor(bugs.auc.sig$Core, levels=c("No", "Yes"))
             
# Lollipop graph
tiff("figures/Classifier_lollipop.tiff", res=300, width = 1400, height = 900, units = "px") ; ggplot(bugs.auc.sig, aes(x=microbe, y=meanAUC, group=pre.post)) +
  xlab("") + ylab("") + 
  geom_segment(aes(x=microbe, xend=microbe, y = 0.45, yend=meanAUC), size=0.35, color="gray80", linetype="dotted") +
  geom_point(aes(color=pre.post, shape=Core), size=2) +
  scale_color_manual(values=col2)+
  scale_x_discrete(limits=lol_order, label=bugs.auc.sig$names[order(bugs.auc.sig$pre.post, bugs.auc.sig$meanAUC, decreasing=T)]) + 
  scale_y_continuous(limits=c(0.44, 0.85)) + 
  theme_classic2()+
  theme(axis.text.x.bottom = element_text(size=5.5, color="black", 
                                       angle=90, hjust=1, vjust=0.3, face = "italic"), 
        axis.text.y = element_text(size=6, color="black"), 
        legend.position = "none") ; dev.off()

# P-value for CRC core microbes ~ Classifier
# 5 core CRC species in 26 Posterior-enriched classifiers 
# 11 core in 258 species
dat <- data.frame(Group_marker = c("Core", "All", "Core", "All"), 
                  Group_condition = c("Posterior", "Posterior", "All", "All"),
                  Count_marker = c(5, 26, 11, 258))
tab <- xtabs(Count_marker~Group_marker+Group_condition, data=dat)
fisher.test(tab, alternative = "two.sided") # 0.01635

# P-value for nonIBD markers ~ Classifier
# 4 nonIBD marker species in 15 Prior-enriched classifiers 
# 14 core in 258 species
dat <- data.frame(Group_marker = c("nonIBD", "All", "nonIBD", "All"), 
                  Group_condition = c("Prior", "Prior", "All", "All"),
                  Count_marker = c(4, 15, 14, 258))
tab <- xtabs(Count_marker~Group_marker+Group_condition, data=dat)
fisher.test(tab, alternative = "two.sided") # 0.0222

# P-value for IBD markers ~ Classifier
# 3 IBD marker species in 26 Posterior-enriched classifiers 
# 12 core in 258 species
dat <- data.frame(Group_marker = c("IBD", "All", "IBD", "All"), 
                  Group_condition = c("Posterior", "Posterior", "All", "All"),
                  Count_marker = c(3, 26, 12, 258))
tab <- xtabs(Count_marker~Group_marker+Group_condition, data=dat)
fisher.test(tab, alternative = "two.sided") # 0.1683



# Integrate additional information of microbes
# 1) Spearman correlation with F.nucleatum 
bugs.spearman <- read.csv("input/spearman.fn.csv", header=T, row.names=1)
bugs.spearman$microbe <- rownames(bugs.spearman) ; colnames(bugs.spearman) <- c("Spearman", "Cor.Pvalue", "Cor.fdr", "NOD_whole", "Class", "microbe")
head(bugs.spearman) ; dim(bugs.spearman) 


# Factor
bugs.auc.stat$Condition <- factor(bugs.auc.stat$Condition, levels=c("nonIBD", "IBD", "Neither"))
bugs.auc.stat$Core <- factor(bugs.auc.stat$Core, levels=c("Yes", "No"))
bugs.auc.stat$pre.post <-factor(bugs.auc.stat$pre.post, levels=c("Prior", "Posterior", "Neither"))
bugs.auc.stat$Significance <- factor(bugs.auc.stat$Significance, levels=c("NS", "Sig"))
  str(bugs.auc.stat)

# X-axis: meanAUC
# Y-axis: Posterior enrichment p-value
# Size: Number of detections across 317 samples (F.nucleatum-exposed)

tiff("figures/CRC_core_diff_enrichment_Pvalue.tiff", res = 200, width = 800, height = 800, units = "px") ; ggplot(bugs.auc.stat) + xlab("") + ylab("") +
  geom_point(aes(meanAUC, -log10(Diff.Pvalue), size=NOD, color=Core, alpha=Core)) +
  geom_vline(xintercept = 0.5, size=1, alpha=0.2, color="gray15", linetype="dotted") +  geom_hline(yintercept = -log10(0.05), size=1, alpha=0.2, color="gray15", linetype="dotted") +
  scale_color_manual(values=c("darkred", "gray94")) + scale_alpha_manual(values = c(0.8, 0.6)) +
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.text.x = element_text(color="black", size=rel(1)),
        axis.text.y = element_text(color = "black", size=rel(1)),
        legend.title = element_text(size = rel(0.9)),
        legend.position = c(0.99, 0.2), 
        legend.justification = c(0.95, 0.95),
        legend.spacing.y = unit(0.8, "mm"), 
        legend.spacing.x = unit(0.1, "mm"), 
        legend.background = element_rect(color = "black", size=0.5),
        legend.text = element_text(size=rel(0.8)),
        legend.key.height=unit(0.8,"line"), # default = 3 ; spacing between the key items
        legend.key = element_rect(fill = "white", color="white")) + guides(size=F, color=guide_legend(title="CRC marker"), alpha=F) ; dev.off()

# For IBD/nonIBD marker 
  # X-axis: meanAUC
  # Y-axis: Differential enrichment p-value
  # Size: Number of detections across 317 samples (F.nucleatum-exposed)

bugs.auc.stat$IBD_marker <- ifelse(bugs.auc.stat$Condition=="IBD", yes="IBD", no="Non") ; bugs.auc.stat$IBD_marker <- factor(bugs.auc.stat$IBD_marker, levels=c("IBD", "Non"))
bugs.auc.stat$nonIBD_marker <- ifelse(bugs.auc.stat$Condition=="nonIBD", yes="nonIBD", no="Non"); bugs.auc.stat$nonIBD_marker <- factor(bugs.auc.stat$nonIBD_marker, levels=c("nonIBD", "Non"))

# IBD marker
tiff("figures/IBD_marker_diff_enrichment_Pvalue.tiff", res = 200, width = 800, height = 800, units = "px") ; ggplot(bugs.auc.stat) + xlab("") + ylab("") +
  geom_point(aes(meanAUC, -log10(Diff.Pvalue), size=NOD, color=IBD_marker, alpha=IBD_marker)) +
  geom_vline(xintercept = 0.5, size=1, alpha=0.2, color="gray15", linetype="dotted") +  geom_hline(yintercept = -log10(0.05), size=1, alpha=0.2, color="gray15", linetype="dotted") +
  scale_color_manual(values=c("darkred", "gray94")) + scale_alpha_manual(values = c(0.8, 0.6)) +
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line = element_line(size=0.7),
        axis.text.x = element_text(color="black", size=rel(1)),
        axis.text.y = element_text(color = "black", size=rel(1)),
        legend.title = element_blank(),
        legend.position = c(0.99, 0.2), 
        legend.justification = c(0.95, 0.95),
        legend.spacing.y = unit(0.3, "mm"), 
        legend.spacing.x = unit(0.1, "mm"), 
        legend.background = element_rect(color = "black", size=0.5),
        legend.text = element_text(size=rel(1.5)),
        legend.key.height=unit(0.8,"line"), # default = 3 ; spacing between the key items
        legend.key = element_rect(fill = "white", color="white")) + guides(size=F, alpha=F); dev.off()

# nonIBD marker
tiff("figures/nonIBD_marker_diff_enrichment_Pvalue.tiff", res = 200, width = 800, height = 800, units = "px") ; ggplot(bugs.auc.stat) + xlab("") + ylab("") +
  geom_point(aes(meanAUC, -log10(Diff.Pvalue), size=NOD, color=nonIBD_marker, alpha=nonIBD_marker)) +
  geom_vline(xintercept = 0.5, size=1, alpha=0.2, color="gray15", linetype="dotted") +  geom_hline(yintercept = -log10(0.05), size=1, alpha=0.2, color="gray15", linetype="dotted") +
  scale_color_manual(values=c("#3399FF", "gray94")) + scale_alpha_manual(values = c(0.8, 0.6)) +
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line = element_line(size=0.7),
        axis.text.x = element_text(color="black", size=rel(1)),
        axis.text.y = element_text(color = "black", size=rel(1)),
        legend.title = element_blank(),
        legend.position = c(0.99, 0.2), 
        legend.justification = c(0.95, 0.95),
        legend.spacing.y = unit(0.3, "mm"), 
        legend.spacing.x = unit(0.1, "mm"), 
        legend.background = element_rect(color = "black", size=0.5),
        legend.text = element_text(size=rel(1.5)),
        legend.key.height=unit(0.8,"line"), # default = 3 ; spacing between the key items
        legend.key = element_rect(fill = "white", color="white")) + guides(size=F, alpha=F); dev.off()


# P-value for significant enrichment at AUC>0.5 & p-value <0.05 (differential enrichment in F. nucleautm-prior or posterior samples)
head(bugs.auc.stat)




#### 9. Prediction model for F.nucleatum-associated dysbiosis ####

# The modeling described below is computationally intense, so we performed using our lab server.

#!/data/program/R-3.4.4/bin/Rscript

options(stringsAsFactors=FALSE)
options(java.parameters = "-Xmx64g")
setwd(dir = "/data/jwhuh/metaphlan2/190504-modeling/")

library(MASS) ; library(pROC) ; library(caret) ; library(RJSONIO)

# Input
bugs <- read.csv("input/bugs_1526.csv", header = T, row.names = 1) 
  # Log transformation with pseudo-abundance
  # 0.00001 is the half of the minimum abundance across whole samples
  bugs <- log10(bugs+0.00001) 

meta <- read.csv("input/meta_1526.csv", header = T, row.names = 1) 
auc.sig <- read.csv("input/bugs.auc.stat.csv", header=T, row.names = 1)
  auc.sort <- subset(auc.sig, fdr < 0.001) # 41 species
  auc.sort <- auc.sort[order(auc.sort$meanAUC, decreasing=T),]
  auc.sort <- auc.sort[auc.sort$meanAUC>0.6 & auc.sort$logFDR > 7, ] # Filtering significant features
  auc.sort <- auc.sort[which(auc.sort$microbe!="Fusobacterium_nucleatum"),] # To remove F.nucleatum
  rownames(auc.sort) <- NULL
  numb <- dim(auc.sort)[1]

# F.nucleatum prior vs. posterior criteria

# Matching bugs & metadata: Critical step!! 
# 1) Microbes ordering
tbugs <- t(bugs) ; tb_order <- order(rownames(tbugs))
tbugs.order <- as.data.frame(tbugs[tb_order,]) # Ordered

# 2) Metadata ordering
meta.order <- meta[order(meta$External.ID, decreasing=F),] # Order-Matched
meta.order$Fn.pre.post.simple <- factor(meta.order$Fn.pre.post.simple, levels=c("Negative", "Prior", "Posterior"))
identical(rownames(tbugs.order), as.character(meta.order$External.ID)) # TRUE
rownames(meta.order) <- meta.order$External.ID

# 3) Add metadata to microbial data
tbugs.order$Condition <- meta.order$Condition
tbugs.order$Fn.pre.post.simple <- meta.order$Fn.pre.post.simple
dim(tbugs.order)

#==============================================
#
# Modeling using microbes + sample condition
#
#==============================================

#-----------Sorting dataset-----------#

# 1) Prior/Posterior sorting
# 1-1) Sorting metadata of F.nucleatum-exposed samples
meta.exp <- meta.order[meta.order$Fn.pre.post.simple!="Negative",]
meta.exp$Fn.pre.post.simple <- droplevels(meta.exp$Fn.pre.post.simple)

# 1-2) Sorting microbial abundance data of F.nucleautm-exposed samples 
bugs.exp <- tbugs.order[meta.order$Fn.pre.post.simple!="Negative",]
bugs.exp$Fn.pre.post.simple <- droplevels(bugs.exp$Fn.pre.post.simple)
dim(bugs.exp) # 317 x 535

#### Generate respective datasets having microbial abundance and Condition/Fn.pre.post.once ###

for (i in 1:numb){
  
  # Bugs list
  assign(paste0("auc.top", i), auc.sort[1:i,]) ; tmp <- get(paste0("auc.top", i))
  assign(paste0("list", i), as.character(tmp$microbe)) 
  
  # Dataset for each microbe list
  tmp.tbugs <- bugs.exp[, c(get(paste0("list",i)), "Condition", "Fn.pre.post.simple")]
  assign(paste0("data", i), tmp.tbugs)
  
}

# So, resulting "data"+"number" contained abundance + metadata features with 1 to "numbered" microbe

# 2) F.nucleatum-negative sorting 
# 2-1) Metadata 
meta.inno <- meta.order[meta.order$Fn.pre.post.simple=="Negative", c("Condition", "Fn.pre.post.simple")]
meta.inno$Fn.pre.post.simple <- droplevels(meta.inno$Fn.pre.post.simple)

# 2-2) Bugs sorting 
bugs.inno <- tbugs.order[meta.order$Fn.pre.post.simple=="Negative", ] # 1209 rows (samples)  x 535 (533 + 2 features)
bugs.inno$Fn.pre.post.simple <- droplevels(bugs.inno$Fn.pre.post.simple)
identical(rownames(meta.inno), rownames(bugs.inno)) # TRUE
dim(bugs.inno)

#==============================================
#
# 		Logistic regression 
#
#==============================================

# Control & Data partitioning
set.seed(29382)
trCtrl <- trainControl(method = "repeatedcv",  # Cross-validation
                       number=5, 
                       repeats=10,
                       classProbs = T, 
                       summaryFunction = twoClassSummary) # twoClassSummary gives AUC

idx <- createDataPartition(data1$Fn.pre.post.simple, times = 100, p = 0.7, list = F)
idx[,100]

# Container
mod.list <- list()
mod.list.neg <- list()
mod.list.post <- list()
coef.list <- list()
coef.sig.list <- list()

#------------Generating modeling function--------------#

for (n in 1:numb) {
  
  # Choose dataset 
  data.tmp <- get(paste0("data", n)) # Total feature sets
  
  # Container for model performance results
  mod.tmp <- as.data.frame(matrix(nrow = 100, ncol = 7))
  colnames(mod.tmp) <- c("AUC_cv", "AUC", "Precision", "Sensitivity", "Specificity", "Accuracy", "AIC")
  
  # Container for posterior probability of F.nucleatum-"exposed" samples
  mod.tmp.posterior <- as.data.frame(matrix(nrow = 317, ncol = 100)) 
  rownames(mod.tmp.posterior) <- rownames(data1)  # Rows: sampleID ; Columns: microbes & feature
  
  # Container for posterior probability of F.nucleatum-"innocent" samples 
  mod.tmp.neg <- as.data.frame(matrix(nrow = 1209, ncol = 100)) 
  rownames(mod.tmp.neg) <- rownames(bugs.inno) 
  
  # Container for regression coefficients of each microbe in 100 models
  # Column length of data.tmp (e.g. data1) will be the row number of container
  mod.tmp.coef <- as.data.frame(matrix(nrow = length(data.tmp), ncol = 100)) 
  coef.name <- colnames(data.tmp)[1:length(data.tmp)-1] # Because the last column of data.tmp is "Fn.pre.post.simple"
  rownames(mod.tmp.coef) <- c("Intercept", coef.name)
  
  mod.tmp.sig <- mod.tmp.coef # For p-value
  
  for (i in 1:100){
    # Partitioning data 
    tmp.idx <- idx[,i] # Temporary index for each run
    train <- data.tmp[tmp.idx,]
    test <- data.tmp[-tmp.idx,]
    
    # Logistic regression 
    glm <- train(form=Fn.pre.post.simple~.,
                 data=train,
                 method="glm",
                 trControl=trCtrl,
                 family="binomial",
                 metric="ROC")
    
    # Prediction
    prd.test <- predict(glm, newdata = test, type = "prob")
    prd.roc <- roc(test$Fn.pre.post.simple, prd.test$Posterior)
    prd.negative <- predict(glm, newdata = bugs.inno, type = "prob")
    
    ###  Outcome summary ###
    # Class prediction
    # Criteria = posterior probability
    # p > 0.5 = Posterior
    # p < 0.5 = Prior
    prd.test$pred.class <- ifelse(test = prd.test$Posterior>0.5, yes = "Posterior", no = "Prior")  
    prd.test$real.class <- test$Fn.pre.post.simple
    prd.test$correct <- ifelse(test = prd.test$pred.class==prd.test$real.class, yes="correct", no="incorrect") 
    
    # Fn.pre.post.once = 36 prior + 59 posterior for split test data
    # Actual positive = 59  = True positive + False negative
    # Actual negative = 36  = True negative + False positive
    
    out.positive <- table(prd.test$pred.class)[1] # True positive + False positive
    out.negative <- table(prd.test$pred.class)[2] # True negative + False negative
    
    TT <- table(prd.test$correct)[1] # True positive + True negative
    TP = (59 + out.positive + TT - 95)/2 # {(TP + FN) + (TP + FP) + (TP + TN) - (TP + TN + FP + FN)}/2 = TP
    FP = (out.positive-TP) # False positive = Predicted positive - True positive
    TN = (TT - TP)  # True negative = All correct - True positive
    FN = (59 - TP)  # False negative = Actual positive - True positive
    
    mod.tmp[i,1] <- glm$results[2] # AUC of mod.tmp at the i-th trial  
    mod.tmp[i,2] <- prd.roc$auc
    mod.tmp[i,3] <- TP/out.positive # Precision: Posterior detection
    mod.tmp[i,4] <- TP/59 # Sensitivity
    mod.tmp[i,5] <- TN/36 # Specificity : Prior detection
    mod.tmp[i,6] <- TT/95 # Accuracy
    mod.tmp[i,7] <- glm$finalModel$aic # Model complexity
    
    # Rows (sampleID) x Columns (Posterior Probability) in F.n-exposed samples
    mod.tmp.posterior[rownames(test), i] <- prd.test$Posterior
    
    # Rows (sampleID) x Columns (Posterior Probability) in F.n-innocent samples
    mod.tmp.neg[,i] <- prd.negative$Posterior 
    
    # Store coefficient & p-value of used features
    mod.tmp.coef[,i] <- summary(glm)$coef[,1] # Coefficients
    mod.tmp.sig[,i] <- summary(glm)$coef[,4] # p-value
  }
  
  rownames(mod.tmp.coef) <- names(glm$finalModel$coefficient)
  
  # Calculate mean Posterior probability for F.n-"experienced" participant
  mod.tmp.posterior$meanPost.prob <- rowMeans(mod.tmp.posterior, na.rm=T)  
  mod.tmp.posterior$pred.class <- ifelse(test = mod.tmp.posterior$meanPost.prob>0.5, yes="Posterior", no="Prior")
  mod.tmp.posterior$real.class <-meta.exp$Fn.pre.post.simple
  
  # Calculate mean Posterior probability for F.n-"inexperienced" participants
  mod.tmp.neg$meanPost.prob <- rowMeans(mod.tmp.neg)
  mod.tmp.neg$pred.class <- ifelse(test = mod.tmp.neg$meanPost.prob>0.5, yes="Posterior", no = "Prior")
  
  # Store results in specified spot of each list
  coef.list[[n]] <- mod.tmp.coef
  coef.sig.list[[n]] <- mod.tmp.sig
  mod.list[[n]] <- mod.tmp
  mod.list.neg[[n]] <- mod.tmp.neg
  mod.list.post[[n]] <- mod.tmp.posterior
}

list_names <- 1:numb
names(coef.list) <- list_names ; names(coef.sig.list) <- list_names ; names(mod.list) <- list_names
names(mod.list.neg) <- list_names ; names(mod.list.post) <- list_names

mod.result <- list(mod.list, mod.list.neg, mod.list.post, coef.list, coef.sig.list)

exportJson <- toJSON(mod.result)
write(exportJson, "model.json")



#### 10. Modeling result analysis ####

# Import Modeling results
# Json format
mod.result <- fromJSON("input/model.json")

# Data structure of "mod.result" file

# 1st [[]] : 
  # 1=Model performance 
  # 2=Posterior probability for F.nucleatum-innocent subject
  # 3=Posterior probability for F.nucleatum-experienced samples
  # 4=Feature Coefficient 
  # 5=Feature significance

# 2nd [[]] : model nubmer = 1-40 

# 3rd $
  # for 1st [[1]] : AUC_cv (training set) / AUC (test set) / Precision / Sensitivity / Specificity / Accuracy / AIC
  # for 1st [[2,3]] : Posterior probabilities from 100 different models ; 101th column = meanPostProbability ; 102th column = summarized class (>0.5 or not)
  # for 1st [[4]] : Rows mean features used for modeling ; Columns mean models (100 models)



#### 10-1) Select the best model ####
  # 1) Satisfactory performance
  # 2) Low complexity (less microbial features used)


# Extract model performance results
numb <- 13 # Number of feature combinations

tmp <- list() ; mod.names <- c("mod.AUC_train", "mod.AUC_test", 
                               "mod.precision", "mod.sensitivity", "mod.specificity", "mod.accuracy",
                               "mod.AIC") ; for (i in 1:7){
  tmp.df <- as.data.frame(matrix(nrow = 100, ncol = numb))
  
  for (n in 1:numb){
    tmp.df[,n] <- mod.result[[1]][[n]][[i]]
    colnames(tmp.df)[n] <- n
  }
  
  #colnames(tmp.df)[40] <- "Features"
  tmp[[i]] <- melt(tmp.df)
  assign(mod.names[i], tmp[[i]])
}

# Model performance comparison figure
perf.figure <- function(data, varaible, value) {
  attach(data)
  ggplot(data, aes(variable, value)) + 
    xlab("") + ylab("") + 
    geom_boxplot(aes(fill=variable), alpha=0.4, outlier.alpha = 0.02, size=rel(0.3)) + 
    theme_classic() + 
    theme(legend.position = "none",
          axis.line = element_line(size=rel(0.5)),
          axis.text.x = element_text(color="gray10", size=rel(0.9)),
          axis.text.y = element_text(color="black", size=rel(0.8)),
          axis.title.y = element_blank())
}



# AUC
tiff(file="figures/Compare.AUC_train.tiff", width = 500, height = 450, units = "px", res = 200) ; perf.figure(mod.AUC_train, variable, value) ; dev.off()
tiff(file="figures/Compare.AUC_test.tiff", width = 500, height = 450, units = "px", res = 200) ; perf.figure(mod.AUC_test, variable, value) ; dev.off()

# Precision
# True positive / Predicted positive (True positive + False negative)
tiff(file="figures/Compare.precision.tiff", width = 500, height = 450, units = "px", res = 200) ; perf.figure(mod.precision, variable, value) ; dev.off()

# Sensitivity
tiff(file="figures/Compare.sensitivity.tiff", width = 500, height = 450, units = "px", res = 200) ; perf.figure(mod.sensitivity, variable, value) ; dev.off()

# Specificity
tiff(file="figures/Compare.specificity.tiff", width = 500, height = 450, units = "px", res = 200) ; perf.figure(mod.specificity, variable, value) ; dev.off()

# Accuracy
tiff(file="figures/Compare.accuracy.tiff", width = 500, height = 450, units = "px", res = 200) ; perf.figure(mod.accuracy, variable, value) ; dev.off()

# AIC
tiff(file="figures/Compare.AIC.tiff", width = 500, height = 450, units = "px", res = 200) ; perf.figure(mod.AIC, variable, log10(value)) ; dev.off()


# Choose the best models
  # We are going to use "0.5" as a criteria for posterior prediction.
  # So, we ranked mean AUC, Precision, Sensitivity, Specificity, Accuracy, AIC (the lower is the better) and averaged them all.
  
# Get average stat
model.comp <- as.data.frame(matrix(nrow=numb, ncol=7)) ; for (i in 1:numb){
  model.comp[i,1] <- mean(mod.AUC_train$value[mod.AUC_train$variable==i]) 
  model.comp[i,2] <- mean(mod.AUC_test$value[mod.AUC_test$variable==i]) 
  model.comp[i,3] <- mean(mod.accuracy$value[mod.accuracy$variable==i]) 
  model.comp[i,4] <- mean(mod.sensitivity$value[mod.sensitivity$variable==i])
  model.comp[i,5] <- mean(mod.specificity$value[mod.specificity$variable==i])
  model.comp[i,6] <- mean(mod.precision$value[mod.precision$variable==i]) 
  model.comp[i,7] <- mean(mod.AIC$value[mod.AIC$variable==i]) 
} ; colnames(model.comp) <- c("AUC_train", "AUC_test", "Accuracy", "Sensitivity", "Specificity", "Precision", "AIC")

# Ranking stat
model.comp$Rank_AUC_train <- rank(-model.comp$AUC_train)
model.comp$Rank_AUC_test <- rank(-model.comp$AUC_test)
model.comp$Rank_Accuracy <- rank(-model.comp$Accuracy)
model.comp$Rank_Sensitivity <- rank(-model.comp$Sensitivity)
model.comp$Rank_Specificity <- rank(-model.comp$Specificity)
model.comp$Rank_Precision <- rank(-model.comp$Precision)
model.comp$Rank_AIC <- rank(model.comp$AIC)
  head(model.comp) ; str(model.comp) 
model.comp$Rank_average <- apply(model.comp[,c(8:14)], 1, mean)
  min(model.comp$Rank_average) ; best.model <- rownames(model.comp)[model.comp$Rank_average==min(model.comp$Rank_average)] ; print(best.model)
  model.comp$Color <- "gray" ; model.comp$Color[model.comp$Rank_average==min(model.comp$Rank_average)] <- "hotpink"
  # Model 10 is the best performer
  

# Model summary figure
  model.comp$Rank_average[10]  # 1.43
  model.comp$Rank_average[11]  # 2.43
  
tiff("figures/Average_rank_models.tiff", res = 300, width = 600, height = 500) ; ggplot(model.comp) + 
  xlab("") + ylab("") + 
  geom_segment(aes(x = 1:numb, xend=1:numb, y = numb, yend=Rank_average, color=Color), linetype="dotted", size=0.3)+
  geom_point(size=2.5, aes(x = 1:numb, y = Rank_average, color=Color)) + 
  geom_text(aes(x = 1:numb, y = Rank_average, label=as.numeric(rownames(model.comp))), color="white", size=1.5) + 
  scale_y_continuous(trans="reverse") + scale_color_manual(values=c("gray60", "#F8766D")) +
  theme_pubr() + 
  theme(axis.line = element_line(size=rel(0.5), color="black"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") ; dev.off()


#### 10-2) Estimating "Posterior probability" in F.nucleatum-experiencd subject ####
  
# Bring posterior probability from F. nucleatume-experienced subjects
Prob_Exp <- as.data.frame(matrix(nrow = 317, ncol = numb)) ; for (n in 1:numb){
  Prob_Exp[,n] <- mod.result[[3]][[n]]$meanPost.prob # mod.result[[3]] contains posterior probability of F. nucleatum-experienced subjects
  colnames(Prob_Exp)[n] <- n
}

# Adding metadata
meta.exp <- meta.order[meta.order$Fn.experience=="Exp",]
  meta.exp$Fn.experience <- droplevels(meta.exp$Fn.experience);  meta.exp$Fn.pre.post.simple <- droplevels(meta.exp$Fn.pre.post.simple)
    identical(rownames(tbugs.order[meta.order$Fn.experience=="Exp",]), meta.exp$External.ID)  # TRUE
  meta.exp$specIBD <- rowSums(tbugs.order[meta.order$Fn.experience=="Exp",ibd.mark]>0)
  meta.exp$specNON <- rowSums(tbugs.order[meta.order$Fn.experience=="Exp",non.mark]>0)
  meta.exp$specRatio <- (meta.exp$specIBD+1)/(meta.exp$specNON+1)
  
  # core CRC microbes detected in F. nucleatum-experienced samples
  bugs_in_all <- colnames(tbugs.order)
  bugs_in_Fn_exp <- bugs_in_all[colSums(tbugs.order[meta.order$Fn.experience=="Exp", ])!=0]
  bugs_in_Fn_inno <- bugs_in_all[colSums(tbugs.order[meta.order$Fn.experience=="Non", ])!=0]
    core_CRC_in_all <- core_CRC[core_CRC %in% bugs_in_all] # 18
    core_CRC_in_Fn_exp <- core_CRC[core_CRC %in% bugs_in_Fn_exp] # 17
    core_CRC_in_Fn_inno <- core_CRC[core_CRC %in% bugs_in_Fn_inno] # 17
  
  # Add the number of detected core CRC microbes in F. nucleatum-experiencd subjects
  meta.exp$specCRC <- rowSums(tbugs.order[meta.order$Fn.experience=="Exp", core_CRC_in_Fn_exp]>0)
  
  
  # Select the best model
  Prob_Exp <- cbind(Prob_Exp[, best.model], subset(meta.exp, select = c("External.ID", "Participant.ID", "Week",  "Condition", "Fn.pre.post.simple", 
                                                              "Shannon", "Evenness", "Richness", "Fn.directionality.strict", 
                                                              "specIBD", "specNON", "specRatio", "specCRC")))
  rownames(Prob_Exp) <- NULL ; colnames(Prob_Exp)[c(1, 11:14)] <- c("Post_prob","IBD", "nonIBD", "IBD/nonIBD", "CRC") ; head(Prob_Exp)  
  
  # Assign putative F. nucleatum-prior/posterior based on posterior probability
  Prob_Exp$Predicted <- ifelse(test = Prob_Exp$Post_prob > 0.5, yes = "Posterior", no = "Prior") ; Prob_Exp$Predicted <- factor(Prob_Exp$Predicted, levels=c("Prior", "Posterior")) ; head(Prob_Exp)
  
# Timelapse change of Posterior probability 
prob_function <- function(data, x, y, color, values) {
  attach(data)
  ggplot(data, aes(x, y, color=color, group=color)) + 
    xlab("") + ylab("") + 
    stat_smooth(geom = "smooth", method = "loess", span=0.8, alpha=0.15, se = T, size=1.2)+
    geom_vline(xintercept = 0, size=0.9, color="gray40", alpha=0.2, linetype="dotted")+
    geom_hline(yintercept = 0.5, size=0.9, color="gray40", alpha=0.2, linetype="dotted")+
    scale_color_manual(values = values) + scale_x_continuous(limits = c(-40, 40)) + 
    theme_pubr() +
    theme(axis.title = element_text(size=rel(0.9), vjust=-1),
          axis.text.x = element_text(size=rel(1.1), color="black"), 
          axis.text.y = element_text(size=rel(1), color="black"),
          axis.line = element_line(size=0.6), 
          legend.text = element_text(size=rel(1.5)), 
          legend.title = element_blank(),
          legend.position = c(0.95,0.05),
          legend.justification = c(1,0),
          legend.background = element_rect(color = "black", size=0.5),
          legend.key = element_rect(fill = "white", color="transparent"),
          legend.key.width = unit(1.1, "line"),
          legend.key.height = unit(1.1,"line"), 
          legend.spacing.y = unit(0, "cm"), 
          legend.spacing.x = unit(0, "cm"))
}

tiff(file="figures/Post_Prob_Exp.tiff", width = 1000, height = 800, units = "px", res = 200); prob_function(Prob_Exp, Fn.directionality.strict, Post_prob, Condition, col2) ; dev.off() 

# Melting for each model
melt_exp <- melt(Prob_Exp, id.vars = c("Predicted", "Condition"), measure.vars = c("Shannon", "Evenness", "Richness", "IBD", "nonIBD", "CRC", "IBD/nonIBD"))


# Mean Posterior Probability from 100 models in F.n-experienced samples
Div_facet <- function(data, class, value, color) {
  attach(data)
  ggplot(data, aes(class, value, fill=class)) + 
    xlab("") + ylab("") + 
    #geom_boxplot(alpha=0.7, outlier.alpha = 0.05) + 
    geom_violin(alpha=0.6, draw_quantiles = 0.5) +
    facet_grid(Condition~variable, scales = "free_x") + 
    stat_compare_means(comparisons = list(c("Prior","Posterior")), 
                       method = "wilcox.test", tip.length = 0, size=rel(2.5), label="p.signif", bracket.size = 0.5)+
    scale_x_discrete(label=c("Pre", "Post")) + 
    scale_fill_manual(values=color) + 
    theme_classic() +
    coord_flip() +
    theme(strip.text = element_text(size=rel(1.5)),
          strip.background = element_rect(color="transparent"),
          axis.text.x = element_text(color="black", size=rel(1)),
          axis.text.y = element_text(color="black", size=rel(1.5)),
          axis.title.y = element_blank(),
          legend.position = "none")
}


tiff(file = "figures/Compare_pseudo_pre_post_Fn_exp.tiff",  width = 2800, height = 600, units = "px", res = 200) ; Div_facet(melt_exp, Predicted, value, col2) ; dev.off()


#### 10-3) Estimating posterior probability in F.nucleatum-innocent samples ####
# Bring posterior probability from F. nucleatume-experienced subjects
Prob_Inno <- as.data.frame(matrix(nrow = 1209, ncol = numb)) ; for (n in 1:numb){
  Prob_Inno[,n] <- mod.result[[2]][[n]]$meanPost.prob
  colnames(Prob_Inno)[n] <- n
}

# Adding metadata
meta.inno <- meta.order[meta.order$Fn.experience=="Non",]
  meta.inno$Fn.pre.post.simple <- droplevels(meta.inno$Fn.pre.post.simple) ; meta.inno$Fn.experience <- droplevels(meta.inno$Fn.experience)

  # Detected core CRC microbes in F. nucleatum-innocent subjects
  meta.inno$specCRC <- rowSums(tbugs.order[meta.order$Fn.experience=="Non", core_CRC_in_Fn_inno]>0)
  meta.inno$specIBD <- rowSums(tbugs.order[meta.order$Fn.experience=="Non", ibd.mark]>0)
  meta.inno$specNON <- rowSums(tbugs.order[meta.order$Fn.experience=="Non", non.mark]>0)
  meta.inno$specRatio <- (meta.inno$specIBD+1)/(meta.inno$specNON+1)
    identical(meta.inno$External.ID, names(rowSums(tbugs.order[meta.order$Fn.experience=="Non", core_CRC_in_Fn_inno]>0))) # TRUE
    identical(rownames(tbugs.order), meta.order$External.ID) # TRUE
  
# Select best model
  Prob_Inno <- data.frame(Prob_Inno[, best.model], subset(meta.inno, select = c("External.ID", "Participant.ID", "Week",  "Condition", "Fn.pre.post.simple", 
                                                                                "Shannon", "Evenness", "Richness",
                                                                                "specIBD", "specNON", "specCRC", "specRatio")))
  rownames(Prob_Inno) <- NULL ; colnames(Prob_Inno)[c(1, 10:13)] <- c("Post_prob", "IBD", "nonIBD", "CRC", "IBD/nonIBD") ; head(Prob_Inno)  
  
  # Assign putative F. nucleatum-prior/posterior based on posterior probability
  Prob_Inno$Predicted <- ifelse(test = Prob_Inno$Post_prob > 0.5, yes = "Posterior", no = "Prior") ; Prob_Inno$Predicted <- factor(Prob_Inno$Predicted, levels=c("Prior", "Posterior"))
  melt_inno <- melt(Prob_Inno, id.vars = c("Predicted", "Condition"), measure.vars = c("Shannon", "Evenness", "Richness", "IBD", "nonIBD", "CRC", "IBD/nonIBD"))
  
# Comparing three diversity indexes & marker species in putative posterior groups in F. n.-innocent subjects
tiff(file = "figures/Compare_pseudo_pre_post_Fn_innocent.tiff", width = 2500, height = 500, units = "px", res = 200) ; Div_facet(melt_inno, Predicted, value, col2) ; dev.off()
  # Significance test for putative posterior in nonIBD vs. putative prior in IBD 
  pred.test <- rbind(Prob_Inno[Prob_Inno$Predicted=="Posterior" & Prob_Inno$Condition=="nonIBD",],
                     Prob_Inno[Prob_Inno$Predicted=="Prior" & Prob_Inno$Condition=="IBD",]) ; rownames(pred.test) <- NULL ;   head(pred.test)
  wilcox.test(data=pred.test, IBD~Condition) # 0.2264
  wilcox.test(data=pred.test, `IBD/nonIBD`~Condition) # 0.07793
  wilcox.test(data=pred.test, CRC~Condition)
  
  

Prob_all <- rbind(subset(Prob_Exp, select = c("External.ID", "Participant.ID", "Week", "Condition", "Predicted", "Post_prob", "IBD", "nonIBD", "CRC", "IBD/nonIBD", "Shannon")),
                  subset(Prob_Inno, select = c("External.ID", "Participant.ID", "Week", "Condition", "Predicted", "Post_prob", "IBD", "nonIBD", "CRC", "IBD/nonIBD", "Shannon"))) ; Prob_all <- Prob_all[order(Prob_all$External.ID),]

cor(x = Prob_all$Post_prob, y = Prob_all$`IBD/nonIBD`, method = "spearman") # 0.5322132
cor(x = Prob_all$Post_prob, y = Prob_all$Shannon, method = "spearman") # -0.2924092

tiff("figures/Post_prob_Marker_Shannon_scatter.tiff", res=250, width = 1000, height = 1000, units = "px"); ggplot(Prob_all, aes(Post_prob, log10(`IBD/nonIBD`), color=Shannon)) + 
  xlab("") + ylab("") +
  scale_color_continuous(low = "#660000", high = "#66FFB2") +
  geom_point(size=rel(1), alpha=0.8) + theme_pubr() + 
  stat_smooth(geom="smooth", method="loess", se = T, span=0.8, color="gray80", alpha=0.4)+
  theme(axis.text = element_text(size=rel(0.8)),
        legend.text = element_text(size=rel(0.6)),
        legend.title = element_text(size=rel(0.8), vjust=1),
        legend.position = c(1, 0.02),
        legend.justification = c(1,0),
        legend.direction = "horizontal",
        legend.background = element_rect(color = "black", size=0.5),
        legend.key = element_rect(fill = "white", color="transparent"),
        legend.key.width = unit(0.7, "line"),
        legend.key.height = unit(0.7,"line"), 
        legend.spacing.y = unit(0, "mm"), 
        legend.spacing.x = unit(1, "mm")); dev.off()

# PCoA using posterior probability
# Matching pcoa dataframe with probability data
pcoa.meta.log <- pcoa.meta.log[order(pcoa.meta.log$External.ID, decreasing = T),] ; Prob_all <- Prob_all[order(Prob_all$External.ID, decreasing = T),]
identical(pcoa.meta.log$External.ID, Prob_all$External.ID) # TRUE

# Merging two dataframes
pcoa.meta.log <- cbind(pcoa.meta.log, subset(Prob_all, select = c("Predicted", "Post_prob", "IBD", "nonIBD", "CRC", "IBD/nonIBD"))) ; head(pcoa.meta.log)
pcoa.meta.log$Predicted_vague <- "Vague" ; pcoa.meta.log$Predicted_vague[pcoa.meta.log$Post_prob>0.7] <- "Posterior"; pcoa.meta.log$Predicted_vague[pcoa.meta.log$Post_prob<0.3] <- "Prior"
  pcoa.meta.log$Predicted_vague <- factor(pcoa.meta.log$Predicted_vague, levels=c("Prior", "Vague", "Posterior")) ;   dim(pcoa.meta.log)


tiff(file = "figures/PCOA.prob.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, Post_prob) + 
  scale_color_gradient2(low = "steelblue", midpoint = 0.5, mid = "white", high = "red") ; dev.off()

tiff(file = "figures/PCOA.posterior_predicted.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, Predicted_vague) + 
  scale_color_manual(values=c("#00AFBB", "gray71", "#FC4E07"), label=c("p<0.3", "0.3<p<0.7", "p>0.7")) +
  stat_ellipse(type="t", aes(group=Predicted_vague, color=Predicted_vague), alpha=0.25, size=1.5) +
  theme(legend.text = element_text(size=rel(1))); dev.off()

tiff(file = "figures/PCOA.CRC.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, CRC) + 
  scale_color_gradient(low = "black", high = "yellow" ) ; dev.off()



#### 10-4) Per Participant analysis ####
Prob_all <- Prob_all[order(Prob_all$Condition, Prob_all$Participant.ID, Prob_all$Week, decreasing = F),] # Ordering
Prob_all$Shannon2 <- Prob_all$Shannon/3.5 
Prob.melt <- melt(Prob_all, id.vars = c("Participant.ID", "Condition", "Week"), measure.vars = c("Post_prob", "Shannon2"))
Post_lm <- lm(Prob_all, formula = Shannon~Post_prob) ; print(Post_lm) ; cor(Prob_all$Shannon, Prob_all$Post_prob, method = "pearson")



tiff("figures/Post_prob_change_per_participant.tiff", res=200, width = 1800, height = 1800, units = "px") ; ggplot(Prob.melt, aes(Week, value, group=variable, fill=variable, color=variable)) + 
  xlab("") + ylab("") +
  geom_point(alpha=0.1) + 
  stat_smooth(geom="line", method = "loess", span=0.8, alpha=0.7, size=1.5)+
  scale_color_manual(values=col2) +
  scale_y_continuous(name = "Posterior probability", limits = c(0, 1), sec.axis = sec_axis(~.*3.5, name = "Shannon diversity")) +
  facet_wrap(~Participant.ID, nrow = 11) +
  geom_hline(yintercept = 0.5, color="gray40", alpha=0.4, linetype="dotted") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(), 
        legend.key = element_rect(fill="white"), 
        legend.text = element_text(size=rel(1)), 
        strip.text = element_text(color="black", size=rel(0.8), 
                                  margin=margin(1,0,1,0, "mm")),
        axis.text = element_text(size=rel(0.75), color="black"), 
        axis.title.y = element_text(size=rel(1), vjust = 1)) ; dev.off()


# Calculate average delta-probability in a given participant # 
head(Prob_all)
Part_name <- names(table(Prob_all$Participant.ID)) ; Part_numb <- table(Prob_all$Participant.ID) 
Part_meanDP <- data.frame(Part_name, delta_prob=0, meanProb=0, stdProb=0, Cor=0) ; tmp_all_df <- as.data.frame(matrix()) ; colnames(tmp_all_df)[1] <- "delta_prob" ; for(i in 1:length(Part_name)){
  
  tmp_name <- Part_name[i] # Name of subject
  tmp_n <- as.numeric(Part_numb[tmp_name]) # Number of samples from the subject
  
  tmp_db <- data.frame(Post_prob=Prob_all[Prob_all$Participant.ID==tmp_name, "Post_prob"],  # Probability difference
                       Week=Prob_all[Prob_all$Participant.ID==tmp_name, "Week"],
                       Shannon=Prob_all[Prob_all$Participant.ID==tmp_name, "Shannon"]) # Week difference
    tmp_df <- as.data.frame(matrix(nrow=tmp_n, ncol=0)) ; tmp_df$delta_prob[1] <- 0
  
  for(n in 1:(tmp_n-1)) {
  
  tmp_df$delta_prob[n+1] <- abs(tmp_db[n,1]-tmp_db[n+1,1])/sqrt(abs(tmp_db[n,2]-tmp_db[n+1,2]))
  
  }  
  
  # Stat for probability
  Part_meanDP$delta_prob[i] <- mean(tmp_df$delta_prob) # average delta probability
  Part_meanDP$meanProb[i] <- apply(tmp_db,2,mean) # average probability
  Part_meanDP$stdProb[i] <- apply(tmp_db,2,sd) # average standard deviation
  Part_meanDP$Cor[i] <- cor(tmp_db$Post_prob, tmp_db$Shannon, method = "pearson")
} ; head(Part_meanDP, n=10)


# Sorting "dynamic subjects" 
dynamic_per70 <- Part_meanDP$Part_name[Part_meanDP$stdProb > quantile(Part_meanDP$stdProb, probs = 0.7)] ; Prob.melt_p70 <- Prob.melt[Prob.melt$Participant.ID %in% dynamic_per70, ] 
dynamic_per90 <- Part_meanDP$Part_name[Part_meanDP$stdProb > quantile(Part_meanDP$stdProb, probs = 0.90)] ; Prob.melt_p90 <- Prob.melt[Prob.melt$Participant.ID %in% dynamic_per90, ] 
dynamic_top12 <- Part_meanDP$Part_name[rank(-Part_meanDP$stdProb)<=12] ; Prob.melt_t12 <- Prob.melt[Prob.melt$Participant.ID %in% dynamic_top12, ] 


# All participants
dy106<- unique(Prob.melt$Participant.ID) ; Part_meanDP$Cor[Part_meanDP$Part_name %in% dy106]
mean(Part_meanDP$Cor[Part_meanDP$Part_name %in% dy106]) # -0.1780158

tiff("figures/Post_prob_change_in_all_participants.tiff", res=200, width = 1800, height = 1400, units = "px") ; ggplot(Prob.melt, aes(Week, value, group=variable, fill=variable, color=variable)) + 
  xlab("") + ylab("") +
  stat_smooth(geom="line", method = "loess", span=0.8, alpha=0.7, size=0.6)+
  scale_color_manual(values=rev(col2)) + scale_y_continuous(name = "Posterior probability", limits = c(0, 1), sec.axis = sec_axis(~.*3.5, name = "Shannon diversity")) +
  facet_wrap(~Participant.ID, ncol = 12) +
  geom_hline(yintercept = 0.5, color="gray40", alpha=0.4, linetype="dotted") +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank(), legend.key = element_rect(fill="transparent", color = "transparent"), 
        legend.text = element_text(size=rel(1.4)),legend.background = element_rect(color="transparent", fill="transparent"),
        strip.text = element_text(color="black", size=rel(0.65), margin=margin(1,0,1,0, "mm")),
        strip.background = element_blank(), axis.text = element_text(size=rel(0.6), color="black"), axis.title.y = element_blank()) ; dev.off()

# Percentile 70 criteria 
dy_p70 <- unique(Prob.melt_p70$Participant.ID)
Part_meanDP$Cor[Part_meanDP$Part_name %in% dy_p70]
Part_meanDP$Part_name[Part_meanDP$Part_name %in% dy_p70]
mean(Part_meanDP$Cor[Part_meanDP$Part_name %in% dy_p70]) # -0.287917


tiff("figures/Post_prob_change_in_dynamic_percentile70.tiff", res=250, width = 1600, height = 1000, units = "px") ; ggplot(Prob.melt_p70, aes(Week, value, group=variable, fill=variable, color=variable)) + 
  xlab("") + ylab("") +
  stat_smooth(geom="line", method = "loess", span=0.7, alpha=0.7, size=rel(0.6))+
  scale_color_manual(values=rev(col2)) +   scale_y_continuous(name = "Posterior probability", limits = c(0, 1), sec.axis = sec_axis(~.*3.5, name = "Shannon diversity")) +
  facet_wrap(~Participant.ID, ncol = 8) +
  geom_hline(yintercept = 0.5, color="gray40", alpha=0.4, linetype="dotted") +
  theme_classic() +
  theme(legend.position = "none", legend.title = element_blank(), 
        legend.key = element_rect(fill="transparent", color = "transparent"), legend.text = element_text(size=rel(1.6)),
        legend.background = element_rect(color="transparent", fill="transparent"),
        strip.text = element_text(color="black", size=rel(0.7), margin=margin(1,0,1,0, "mm")),
        strip.background = element_rect(color="transparent"), 
        axis.text = element_text(size=rel(0.4), color="black"), axis.line = element_line(size=0.3),axis.title.y = element_blank()) ; dev.off()

# Top 12 pearson correlation
dy_t12 <- unique(Prob.melt_t12$Participant.ID) 
Part_meanDP$Cor[Part_meanDP$Part_name %in% dy_t12]
mean(Part_meanDP$Cor[Part_meanDP$Part_name %in% dy_t12]) # -0.2855085


tiff("figures/Post_prob_change_in_dynamic_top12.tiff", res=250, width = 900, height = 600, units = "px") ; ggplot(Prob.melt_t12, aes(Week, value, group=variable, fill=variable, color=variable)) + 
  xlab("") + ylab("") +
  stat_smooth(geom="line", method = "loess", span=0.7, alpha=0.7, size=rel(0.5))+
  scale_color_manual(values=rev(col2)) +  
  scale_x_continuous(labels = c(0,20,40), breaks = c(0,20,40)) + scale_y_continuous(name = "Posterior probability", limits = c(0, 1), sec.axis = sec_axis(~.*3.5, name = "Shannon diversity")) +
  facet_wrap(~Participant.ID, nrow = 2) +
  geom_hline(yintercept = 0.5, color="gray40", alpha=0.4, linetype="dotted") +
  theme_classic() +
  theme(legend.position = "none",legend.title = element_blank(), legend.key = element_rect(fill="transparent", color = "transparent"), 
        legend.text = element_text(size=rel(1.6)),legend.background = element_rect(color="transparent", fill="transparent"),
        strip.text = element_text(color="black", size=rel(0.5), margin=margin(1,0,1,0, "mm")),
        strip.background = element_rect(color="transparent"), 
        axis.text = element_text(size=rel(0.4), color="black"), 
        axis.line = element_line(size=0.4),
        axis.title.y = element_blank()) ; dev.off()

# Dynamic subjecst and microbial change in PCoA plot
# 1) dynamic top-12
head(pcoa.meta.log)
pcoa.d12 <- pcoa.meta.log[pcoa.meta.log$Participant.ID %in% dy_t12,] 
pcoa.d12 <- pcoa.d12[order(pcoa.d12$Condition, pcoa.d12$Participant.ID, pcoa.d12$visit.num, decreasing = F),] ; dim(pcoa.d12) # 174 x 49


tiff("figures/PCoA.d12_subject.path.tiff", res = 200, width = 800, height = 2400, units = "px"); ggplot(pcoa.d12, aes(PC1, PC2)) + 
  geom_path(aes(color=Post_prob), arrow = NULL, lineend = "round" , size=1.2)+
  xlab(paste0("PC1 (", x.var, "%)")) + ylab(paste0("PC2 (", y.var, "%)")) + 
  facet_grid(Participant.ID~Condition) + 
  geom_hline(yintercept = 0, linetype="dotted", color="gray40", size=1, alpha=0.5)+
  geom_vline(xintercept = 0, linetype="dotted", color="gray40", size=1, alpha=0.5)+
  scale_color_gradient2(low = "green", mid = "black", midpoint = 0.5, high = "red") + 
  theme_few() + 
  theme(axis.text = element_text(size=rel(0.7)), 
        axis.title.x = element_text(vjust=-1.5, size=rel(1.2)), 
        axis.title.y = element_text(vjust=1.7, size=rel(1.2)),
        strip.text.x = element_text(size=rel(2.3)),
        strip.text.y = element_text(size=rel(1.1)),
        legend.text = element_text(size=rel(0.8)), 
        legend.title = element_blank(), 
        legend.position = "right",
        legend.justification = c(1, 0),
        legend.background = element_rect(color = "black", size=0.5),
        legend.key = element_rect(fill = "white"),
        legend.key.width = unit(0.8, "line"),
        legend.key.height = unit(0.8,"line"), 
        legend.spacing.y = unit(0, "mm"), # Spacing between legend title and items
        legend.spacing.x = unit(0.6, "mm")) ; dev.off()

# 2) dynamic 70th percentile
pcoa.p70 <- pcoa.meta.log[pcoa.meta.log$Participant.ID %in% dy_p70,] 
pcoa.p70 <- pcoa.p70[order(pcoa.p70$Condition, pcoa.p70$Participant.ID, pcoa.p70$visit.num, decreasing = F),] ; dim(pcoa.p70) # 468 x 49

tiff("figures/PCoA.p70_subject.path.tiff", res = 200, width = 1400, height = 2400, units = "px"); ggplot(pcoa.p70, aes(PC1, PC2)) + 
  geom_path(aes(color=Post_prob), arrow = NULL, lineend = "round" , size=1.2)+
  xlab(paste0("PC1 (", x.var, "%)")) + ylab(paste0("PC2 (", y.var, "%)")) + 
  facet_wrap(~Participant.ID, ncol = 4) + 
  geom_hline(yintercept = 0, linetype="dotted", color="gray40", size=1, alpha=0.5)+
  geom_vline(xintercept = 0, linetype="dotted", color="gray40", size=1, alpha=0.5)+
  scale_color_gradient2(low = "green", mid = "black", midpoint = 0.5, high = "red") + 
  theme_few() + 
  theme(axis.text = element_text(size=rel(0.7)), 
        axis.title.x = element_text(vjust=-1.5, size=rel(1.2)), 
        axis.title.y = element_text(vjust=1.7, size=rel(1.2)),
        strip.text = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(0.8)), 
        legend.title = element_blank(), 
        legend.position = "right",
        legend.justification = c(1, 0),
        legend.background = element_rect(color = "black", size=0.5),
        legend.key = element_rect(fill = "white"),
        legend.key.width = unit(0.8, "line"),
        legend.key.height = unit(0.8,"line"), 
        legend.spacing.y = unit(0, "mm"), # Spacing between legend title and items
        legend.spacing.x = unit(0.6, "mm")) ; dev.off()

#### 11. Clustering microbes based on their longitudinal abundance dynamics ####

# K-mean clustering 

# Features for clustering
  # 1) Spearman correlation coefficient with F. nucleatum
  bugs.spearman <- read.csv("input/spearman.fn.csv", header=T, row.names=1) ; bugs.spearman$microbe <- rownames(bugs.spearman) ; str(bugs.spearman)  # 533 bugs

  # 2) Differential enrichment in IBD condition
  mic_con <- read.csv("input/p.value_IBD.nonIBD.bugs.csv", header=T, row.names=1) ; mic_con$microbe <- rownames(mic_con) ; str(mic_con) # 533 bugs
  
  # 3) Differential enrichment in F.nucleatum-experienced subject
  mic_exp <- read.csv("input/Fnucleatum_exposed_Pvalue.csv", header=T, row.names=1) ; mic_exp$microbe <- rownames(mic_exp) ; str(mic_exp) # 533 bugs
  
  # 4) Differential enrichment in F.nucleatum-posterior samples
  mic_pp <- read.csv("input/Fnucleatum_pre_post_Pvalue.csv", header=T, row.names=1) ; mic_pp$microbe <- rownames(mic_pp) ; str(mic_pp) # 258 bugs
  
  # 5) Classifying p-value 
  bugs.auc.stat <- read.csv("input/bugs.auc.stat.csv", header = T, row.names = 1) ; str(bugs.auc.stat) # 258 x 12
  

# Merging in one dataframe
bc_533 <- data.frame(microbe=bugs.spearman$microbe, 
                     NOD=bugs.spearman$NOD,
                     Cor.Pvalue=bugs.spearman$p.value,
                     IBD.Pvalue=mic_con$IBD.Pvalue, 
                     Exp.Pvalue=mic_exp$Exp.Pvalue) ; head(bc_533)

bc_258 <- data.frame(microbe=bugs.auc.stat$microbe,
                     Post.Pvalue=mic_pp$Post.Pvalue,
                     Class.Pvalue=bugs.auc.stat$p.value,
                     Condition=factor(bugs.auc.stat$Condition, levels=c("Neither", "nonIBD", "IBD")),
                     Core=factor(bugs.auc.stat$Core, labels = c("Non", "CRC"), levels = c("No", "Yes")),
                     pre.post=factor(bugs.auc.stat$pre.post, levels=c("Neither", "Prior", "Posterior")),
                     Classifier=factor(bugs.auc.stat$Significance, labels = c("Non", "Classifier"), levels=c("NS", "Sig"))) ; head(bc_258) 

bc_258 <- merge(bc_258, bc_533, by.x = "microbe") ; bc_258 <- bc_258[,c(1,9:11,2:8)] ; head(bc_258)
  # Top-10 classifiers
  t10 <- bugs.auc.stat[c(1:10),"microbe"] ; bc_t10 <- bc_258[bc_258$microbe %in% t10,]
  bc_258$T10 <- ifelse(test = bc_258$microbe %in% t10, yes = "Top10", no = "Non") ; bc_258$T10 <- factor(bc_258$T10, levels=c("Non", "Top10"))

bugs_258 <- bc_258$microbe

#write.csv(bc_258, "input/merge.microbe_simple.csv")
#bc_258 <- read.csv("input/merge.microbe_simple.csv", header = T, row.names=1) ; head(bc_258)

bc_km <- bc_258[,c(2:6)] ; rownames(bc_km) <- bc_258$microbe; head(bc_km) 
bc_km["Fusobacterium_nucleatum","Cor.Pvalue"] <- 0 # NA to 0

#### 11-1) Determining number of clusters (k-mean) ####
set.seed(12345)
tmp.df <- as.data.frame(matrix(nrow = 20, ncol = 6)) ; for (i in 1:20) {
  tmp <- kmeans(bc_km, i)
  tmp.df[i,1] <- max(tmp$withinss) # Maximum withinss
  tmp.df[i,2] <- min(tmp$withinss) # Minimum withinss
  tmp.df[i,3] <- mean(tmp$withinss) # Mean withinss
  tmp.df[i,4] <- tmp$tot.withinss
  tmp.df[i,5] <- tmp$betweenss
  tmp.df[i,6] <- tmp$totss
} ; colnames(tmp.df) <- c("Max.withinss", "Min.withinss", "Avg.withinss", "tot.withinss", "betweenss", "totss") ; tmp.df$kmean <- 1:20 ; tmp.df$totss[1] # 199.0665

tmp_melt <- melt(tmp.df, id.vars = "kmean",  measure.vars = c("Max.withinss", "Min.withinss", "Avg.withinss", "tot.withinss", "betweenss"))

# To visualize different set of sum of square in a one plot
tiff("figures/kmean_number.tiff", units = "px", res = 200, width = 750, height = 450) ; ggplot(tmp_melt, aes(kmean, value, group=variable, color=variable)) + 
  xlab("")+ylab("") + geom_line(alpha=0.6) +   geom_point() + 
  geom_hline(yintercept = 199.0665, color="gray70", size=1.5, alpha=0.5) +
  theme_classic() + theme(legend.title = element_blank(), legend.text = element_text(rel(0.5))) ; dev.off()   

# Chosing the best number of clusters 
set.seed(12345) ; NbClust(as.matrix(bc_km), distance = "manhattan", min.nc = 3, max.nc = 10, method = "kmeans")
# 9 cluster was chosen by the majority rule 
  # 8 proposed 9 as the best number of clusters 
  # Next was cluster number 3 chosen by 4 
  # According to the majority rule, the best number of clutser is 9 

#### 11-2) K-mean clustering & PCoA plot ####
km <- kmeans(bc_km, centers = 9, nstart = 100, iter.max = 100, algorithm = "Lloyd")
  identical(rownames(bc_km), names(km$cluster)) # TRUE
  km$cluster["Fusobacterium_nucleatum"] # 6
  km$cluster["Dorea_longicatena"] # 4

# PCoA plot for bugs 
# 1) By abundance
bugs_258 <- log10(bugs[bc_258$microbe,]+1e-05)+5
db <- vegdist(bugs_258, method = "bray") ; pcoa.bugs <- cmdscale(db, k = 2, eig = T)

  # variance explained
  pcoa.bugs$eig[1:2] ;  sum(pcoa.bugs$eig[1:2])^2/sum((pcoa.bugs$eig))^2
  x.var <- round((sum(pcoa.bugs$eig[1]^2))/sum((pcoa.bugs$eig)^2)*100, 1)
  y.var <- round((sum(pcoa.bugs$eig[2]^2))/sum((pcoa.bugs$eig)^2)*100, 1)

# 2) By dynamics
dd <- vegdist(bc_km, method = "euclidean") ; pcoa.bugs.dynamic <- cmdscale(dd, k = 2, eig = T)
  
  # variance explained
  pcoa.bugs.dynamic$eig[1:2] ;  sum(pcoa.bugs.dynamic$eig[1:2])^2/sum((pcoa.bugs.dynamic$eig))^2
  x.var <- round((sum(pcoa.bugs.dynamic$eig[1]^2))/sum((pcoa.bugs.dynamic$eig)^2)*100, 1)
  y.var <- round((sum(pcoa.bugs.dynamic$eig[2]^2))/sum((pcoa.bugs.dynamic$eig)^2)*100, 1)


pcoa.plot <- function(data, class) {
  attach(data)
  ggplot(data, aes(PC1, PC2)) + geom_point(aes(color=class), alpha=0.85) + 
    xlab(paste0("PC1 (", x.var, "%)")) + ylab(paste0("PC2 (", y.var, "%)")) + 
    theme_few() + 
    theme(axis.text = element_blank(), 
          axis.title.x = element_text(vjust=-1.5, size=rel(1.2)), 
          axis.title.y = element_text(vjust=1.7, size=rel(1.2)), 
          legend.text = element_text(size=rel(1)), 
          legend.title = element_blank(),
          legend.position = c(0.02, 0.02),
          legend.justification = c(0,0),
          legend.background = element_rect(color = "black", size=0.5),
          legend.key = element_rect(fill = "white"),
          legend.key.width = unit(0.8, "line"),
          legend.key.height = unit(0.8,"line"), 
          legend.spacing.y = unit(0, "mm"), # Spacing between legend title and items
          legend.spacing.x = unit(0.6, "mm"))
}

# Abundance PCoA
pcoa.bugs.df <- data.frame(microbe=rownames(pcoa.bugs$points),
                           PC1=pcoa.bugs$points[,1],
                           PC2=pcoa.bugs$points[,2],
                           NOD=bc_258$NOD,
                           Classifier=bc_258$Classifier,
                           Core=bc_258$Core,
                           pre.post=bc_258$pre.post,
                           Condition=bc_258$Condition, 
                           Cluster=factor(km$cluster, levels=1:length(km$size))) ; rownames(pcoa.bugs.df) <- NULL; head(pcoa.bugs.df) 

tiff("figures/PCOA_bugs_Core.tiff", res = 200, units = "px", width = 800, height = 800) ; pcoa.plot(pcoa.bugs.df, Core) + scale_color_manual(values=col2.1) ; dev.off()
tiff("figures/PCOA_bugs_Classifier.tiff", res = 200, units = "px", width = 800, height = 800) ; pcoa.plot(pcoa.bugs.df, Classifier) + scale_color_manual(values=col2.1) ; dev.off()
tiff("figures/PCOA_bugs_Condition.tiff", res = 200, units = "px", width = 800, height = 800) ; pcoa.plot(pcoa.bugs.df, Condition) + scale_color_manual(values=col3.1) ; dev.off()
tiff("figures/PCOA_bugs_pre_post.tiff", res = 200, units = "px", width = 800, height = 800) ; pcoa.plot(pcoa.bugs.df, pre.post) + scale_color_manual(values=col3.1) ; dev.off()
tiff("figures/PCOA_bugs_kmean_clusters.tiff", res = 200, units = "px", width = 900, height = 800) ; pcoa.plot(pcoa.bugs.df, Cluster) + 
  theme(legend.position = "right", legend.background = element_rect(color="transparent"), legend.title = element_blank()) ; dev.off()


# Dynamic PCoA
head(pcoa.bugs.dynamic)
pcoa.dynamic.df <- data.frame(microbe=rownames(pcoa.bugs.dynamic$points),
                              PC1=pcoa.bugs.dynamic$points[,1],
                              PC2=pcoa.bugs.dynamic$points[,2],
                              NOD=bc_258$NOD,
                              Classifier=bc_258$Classifier,
                              Core=bc_258$Core,
                              pre.post=bc_258$pre.post,
                              Condition=bc_258$Condition, 
                              #Top10=bc_258$T10,
                              Cluster=factor(km$cluster, levels=1:length(km$size))) ; rownames(pcoa.dynamic.df) <- NULL; head(pcoa.dynamic.df) 

tiff("figures/PCOA_dynamic_Core.tiff", res = 200, units = "px", width = 800, height = 800) ; pcoa.plot(pcoa.dynamic.df, Core) + scale_color_manual(values=col2.1) ; dev.off()
tiff("figures/PCOA_dynamic_Classifier.tiff", res = 200, units = "px", width = 800, height = 800) ; pcoa.plot(pcoa.dynamic.df, Classifier) + scale_color_manual(values=col2.1) ; dev.off()
tiff("figures/PCOA_dynamic_Condition.tiff", res = 200, units = "px", width = 800, height = 800) ; pcoa.plot(pcoa.dynamic.df, Condition) + scale_color_manual(values=col3.1) ; dev.off()
tiff("figures/PCOA_dynamic_pre_post.tiff", res = 200, units = "px", width = 800, height = 800) ; pcoa.plot(pcoa.dynamic.df, pre.post) + scale_color_manual(values=col3.1) ; dev.off()
#tiff("figures/PCOA_dynamic_Top10.tiff", res = 200, units = "px", width = 800, height = 800) ; pcoa.plot(pcoa.dynamic.df, Top10) + scale_color_manual(values=col2.1) ; dev.off()

color_9 <- rev(brewer.pal(n = 9, name = rev("Spectral")))
tiff("figures/PCOA_dynamic_kmean_clusters.tiff", res = 200, units = "px", width = 950, height = 800) ; pcoa.plot(pcoa.dynamic.df, Cluster) + 
#  stat_ellipse(aes(group=Cluster), alpha=0.3, type = "t")+
  scale_color_manual(values=color_9) +
  theme(legend.position = "right", legend.background = element_rect(color="transparent"), legend.text = element_text(size=rel(1.4)), legend.title = element_blank()) ; dev.off()


# Draw histogram
cluster_histogram <- function(data, clusters, class, color){
  attach(data)
  ggplot(data, aes(clusters, fill=class)) + 
    xlab("") + ylab("") +
    geom_histogram(stat = "count", alpha=0.7) + 
    scale_fill_manual(values = color) + 
    theme_classic() +
    theme(legend.position = "right",
          legend.justification = c(0, 0),
          legend.title = element_blank(),
          legend.text = element_text(size=rel(0.6)),
          legend.key.height = unit(1.5, "mm"),
          legend.key.width = unit(3, "mm"),
          legend.spacing.y = unit(0.1, "mm"), 
          legend.spacing.x = unit(0.3, "mm"),
          axis.text.x = element_text(color="black", size=rel(0.9)),
          axis.text.y = element_text(color="black", size=rel(0.8)))
}


tiff("figures/histogram_cluster_condition.tiff", res = 200, width = 600, height = 600, units = "px") ; cluster_histogram(pcoa.dynamic.df, Cluster, Condition, col3.1); dev.off()
tiff("figures/histogram_cluster_pre_post.tiff", res = 200, width = 600, height = 600, units = "px") ; cluster_histogram(pcoa.dynamic.df, Cluster, pre.post, col3.1); dev.off()
tiff("figures/histogram_cluster_Core.tiff", res = 200, width = 600, height = 600, units = "px") ; cluster_histogram(pcoa.dynamic.df, Cluster, Core, col2.1); dev.off()
tiff("figures/histogram_cluster_classifier.tiff", res = 200, width = 600, height = 600, units = "px") ; cluster_histogram(pcoa.dynamic.df, Cluster, Classifier, col2.1); dev.off()


# Combinding IBD + core 
pcoa.dynamic.df$Pathogenic <- "None" ; pcoa.dynamic.df$Pathogenic[pcoa.dynamic.df$Core=="CRC" & pcoa.dynamic.df$Condition == "IBD"] <- 
  "CRC+IBD" ; pcoa.dynamic.df$Pathogenic[pcoa.dynamic.df$Core=="CRC" & pcoa.dynamic.df$Condition == "nonIBD"] <- 
  "CRC+nonIBD" ; pcoa.dynamic.df$Pathogenic[pcoa.dynamic.df$Core=="CRC" & pcoa.dynamic.df$Condition == "Neither"] <- 
  "CRC" ; pcoa.dynamic.df$Pathogenic[pcoa.dynamic.df$Core!="CRC" & pcoa.dynamic.df$Condition == "IBD"] <- 
  "IBD" ; pcoa.dynamic.df$Pathogenic[pcoa.dynamic.df$Core!="CRC" & pcoa.dynamic.df$Condition == "nonIBD"] <- 
  "nonIBD" ; pcoa.dynamic.df$Pathogenic <- factor(pcoa.dynamic.df$Pathogenic, levels=c("CRC+IBD", "CRC", "IBD", "CRC+nonIBD", "nonIBD", "None")) ; table(pcoa.dynamic.df$Pathogenic)
  write.csv(pcoa.dynamic.df, "input/pcoa.dynamic.df.csv")
tiff("figures/histogram_cluster_pathogenic.tiff", res = 200, width = 600, height = 500, units = "px") ; cluster_histogram(pcoa.dynamic.df, Cluster, Pathogenic, c("darkred", "red", "orange", "steelblue", "skyblue", "gray81"))  ; dev.off()


#### 11-3) Temporal distribution by clusters ####
str(pcoa.dynamic.df)
dat1 <- melt(pcoa.dynamic.df, id.vars = "Cluster", measure.vars = c("Condition", "Core", "Classifier", "pre.post", "Pathogenic"))


# Assign clusters to bugs
for (i in 1:length(km$size)){
  assign(paste0("Cluster", i), names(km$cluster[km$cluster==i]))
}



# bugs + meta data.frame
tbugs.order_ID <- cbind(tbugs.order, meta.order$External.ID) ; colnames(tbugs.order_ID)[length(tbugs.order_ID)] <- "External.ID" ; head(tbugs.order_ID) 
bm <- merge(tbugs.order_ID, meta.order, by = "External.ID") # 1526 x 571


# Get "sum of abundance" or "number of detection" for each cluster across whole samples
km.size <- length(km$size)
cl_df <- as.data.frame(matrix(nrow = 1526, ncol = km.size*3)) ; for (i in 1:km.size){
  tmp_cluster <- get(paste0("Cluster",i))
  cl_df[,i] <- rowSums(bm[ ,tmp_cluster]) # Sum of abundance
  cl_df[,i+km.size] <- rowSums(bm[ , tmp_cluster]>0) # Number of detection
  cl_df[,i+km.size*2] <- round(cl_df[,i+km.size]/length(tmp_cluster)*100, 2) # % of total cluster number
} ; head(cl_df)

cl_df <- cbind(cl_df, bm$External.ID, bm$Fn.experience, bm$Fn.directionality.strict, bm$Fn.pre.post.simple, bm$Condition)
colnames(cl_df)[c((km.size*3+1):(km.size*3+5))] <- c("External.ID", "Fn.experience", "Fn.directionality.strict", "Fn.pre.post.simple", "Condition")
cl_df <- merge(cl_df, subset(Prob_all, select =  colnames(Prob_all)[-which(colnames(Prob_all)=="Condition")]), by = "External.ID") ; head(cl_df)

# Sum of abundance
cl_sum <- melt(cl_df, id.vars = c("External.ID", "Condition", "Fn.directionality.strict", "Fn.pre.post.simple", "Fn.experience", "Post_prob"), measure.vars = paste0("V", 1:km.size)) 
cl_sum$variable <- factor(cl_sum$variable, labels = c(1:km.size)) ; cl_sum$variable

# Number of detected cluster components
cl_spec <- melt(cl_df, id.vars = c("External.ID", "Condition", "Fn.directionality.strict", "Fn.pre.post.simple", "Fn.experience", "Post_prob"), measure.vars = paste0("V", (km.size+1):(km.size*2)))
cl_spec$variable <- factor(cl_spec$variable, labels = c(1:km.size)) ; cl_spec$variable

# Percentage for detected cluster components
cl_spec_perc <- melt(cl_df, id.vars = c("External.ID", "Condition", "Fn.directionality.strict","Fn.pre.post.simple", "Fn.experience", "Post_prob"), measure.vars = paste0("V", (km.size*2+1):(km.size*3)))
cl_spec_perc$variable <- factor(cl_spec_perc$variable, labels = c(1:km.size)) ; cl_spec_perc$variable

draw_cluster <- function(data, x, y, color, label.y){
  attach(data)
  ggplot(data, aes(x, y, color=color)) + 
    xlab("") + ylab("") +
    geom_vline(xintercept = 0, linetype="dotted", size=0.4, color="gray40", alpha=0.5) +
    stat_smooth(method = "loess", span=0.8, se = T, size=0.6, alpha=0.15) +
    stat_cor(method = "spearman", label.y=label.y, size=rel(2)) +
    facet_wrap(~variable, nrow = 3)+ 
    scale_x_continuous(limits = c(-45, 45))+
    scale_color_manual(values=col2) +
    theme_classic() +
    theme(legend.position = "none",
          #legend.justification = c(1,0),
          #legend.title = element_blank(),
          #legend.text = element_text(size=rel(0.5)),
          #legend.key.size = unit(1, "line"),
          strip.text = element_text(size=rel(1), margin = margin(0.6,0,0.6,0, "mm")),
          strip.background = element_rect(color="transparent"),
          axis.text.x = element_text(color="black", size=rel(0.7)),
          axis.text.y = element_text(color="black", size=rel(0.7)))
}

tiff("figures/Cluster_sum_abundance_direct.tiff", res = 200, units = "px", width = 800, height = 800) ; draw_cluster(cl_sum, Fn.directionality.strict, log10(value+0.001), Condition, c(-2.3, -1.6)) ; dev.off()
tiff("figures/Cluster_species_number_direct.tiff", res = 200, units = "px", width = 800, height = 800) ; draw_cluster(cl_spec, Fn.directionality.strict, value, Condition, c(19, 21.5)); dev.off() 
tiff("figures/Cluster_species_percent_direct.tiff", res = 200, units = "px", width = 800, height = 800) ; draw_cluster(cl_spec_perc, Fn.directionality.strict, value, Condition, c(0, 6.5)) ; dev.off() 

#### 11-4) Posterior probability ~ Clutsers ####
head(cl_spec)
tiff("figures/Cluster_Post_prob.tiff", res = 200, units = "px", width = 800, height = 900) ; ggplot(cl_spec, aes(Post_prob, value, group=variable, color=variable)) + 
  geom_jitter(alpha=0.1, size=rel(0.5)) +
  xlab("") + ylab("") + 
  stat_smooth(geom="line", method = "loess", span=0.8, alpha=0.7, size=1.1) + 
  #stat_cor(method = "spearman", label.y=20, color="gray10", size=rel(1.7)) +
  scale_color_manual(values=color_9)+
  scale_x_continuous(labels = c(0, 0.5, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_wrap(~variable, nrow = 3) + theme_classic() +
  theme(legend.position = "none", 
        #legend.justification = c(0,0),
        #legend.direction = "vertical", 
        #legend.title = element_blank(), 
        #legend.text = element_text(size=rel(0.7)),
        #legend.key.height = unit(0.1, "line"),
        #legend.spacing.y = unit(0, "line"),
        strip.background = element_rect(color="transparent"),
        strip.text = element_text(size=rel(1)),
        axis.text.x = element_text(size=rel(0.9)),
        axis.text.y = element_text(size=rel(0.8))) ; dev.off()


#### 11-5) On PCoA plot ####
str(pcoa.meta.log)
Cluster1 # 25
Cluster6 # 46
identical(rownames(tbugs.order), pcoa.meta.log$External.ID) # TRUE

# Sum of abundance per sample
# Cluster 1 (beneficial)
pcoa.meta.log$C1_sum <- rowSums(tbugs.order[,Cluster1])
pcoa.meta.log$C1_spec <- rowSums(tbugs.order[,Cluster1]>0)
pcoa.meta.log$C6_sum <- rowSums(tbugs.order[,Cluster6])
pcoa.meta.log$C6_spec <- rowSums(tbugs.order[,Cluster6]>0)

tiff("figures/pcoa.C1_sum_log.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, class = log10(C1_sum+1)) + scale_color_gradient(low = "black", high = "yellow") +
  xlab("") + ylab(""); dev.off()
tiff("figures/pcoa.C1_spec.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, class = C1_spec) + scale_color_gradient(low = "black", high = "yellow") +
  xlab("") + ylab(""); dev.off()

tiff("figures/pcoa.C6_sum_log.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, class = log10(C6_sum+1)) + scale_color_gradient(low = "black", high = "yellow") +
  xlab("") + ylab(""); dev.off()
tiff("figures/pcoa.C6_spec.tiff", width = 800, height = 800, units = "px", res = 200) ; pcoa.plot(pcoa.meta.log, class = C6_spec) + scale_color_gradient(low = "black", high = "yellow")+
  xlab("") + ylab(""); dev.off()


#### 11-6) Each cluster ####
Cluster4

#### 12. Replicate analysis ####

# Count the number of paired detection for each microbe and compare along with average abundance
colnames(re) ; dim(re) # name of replicated 44 samples (22 pairs) x 218 species

re$mean_Abundance <- round(rowMeans(re), 4)
re$NOD <- rowSums(re[,1:44]>0)
re$detected_Abundance <- round(rowSums(re[,1:44])/re$NOD, 4)

re.list <- list() ; df1 <- as.data.frame(matrix(nrow = 218, ncol = 22)) ; df2 <- as.data.frame(matrix(nrow = 218, ncol = 22)) ; for (i in 1:22){
  
  n1 = 2*i-1
  n2 = 2*i
  
  # mean abundance per pair
  tmp.re <- re[,c(n1, n2)]
  df1[,i] <- round(rowMeans(tmp.re), 4)
  
  # Matched 
  tmp.rule <- rowSums(tmp.re[,1:2]>0)
  tmp.re$Matched[tmp.rule==0] <- "Undetected"
  tmp.re$Matched[tmp.rule==1] <- "Once"
  tmp.re$Matched[tmp.rule==2] <- "Twice"
  
  df2[,i] <- tmp.re$Matched
  colnames(df1)[i] <- colnames(re)[n2]
  colnames(df2)[i] <- colnames(re)[n2]
  
} ; rownames(df1) <- rownames(re) ; rownames(df2) <- rownames(re) ; df1$microbes <- rownames(df1) ; df2$microbes <- rownames(df2) ; head(df1) ; head(df2)

# Melting 
melt1 <- melt(df1, id.vars = "microbes") ; melt1$Matched <- melt(df2, id.vars = "microbes", measure.vars = colnames(df2)[1:22])$value ;
colnames(melt1) <- c("Microbe", "Pair.ID", "mean_Abundance", "Matched")
melt1$Matched <- factor(melt1$Matched, levels = c("Undetected", "Once", "Twice"))

# Name
b.name <- unique(melt1$Microbe) # 218
id.name <- as.character(unique(melt1$Pair.ID)) # 22


pair.df <- as.data.frame(matrix(nrow = 22, ncol=8)) ; for (n in 1:length(id.name)){
  
  # Pair ID
  pair.df[n,1]   <- id.name[n]
  
  # Temporary files
  tmp <- melt1[melt1$Pair.ID==id.name[n],]
  tmp_twice <- tmp[tmp$Matched=="Twice",]
  tmp_once <- tmp[tmp$Matched=="Once",]
  tmp_un <- tmp[tmp$Matched=="Undetected",]
  
  # Number of species per recovery
  pair.df[n,2] <- length(rownames(tmp_twice)) 
  pair.df[n,3] <- length(rownames(tmp_once)) 
  pair.df[n,4] <- pair.df[n,3]/(pair.df[n,2] + pair.df[n,3]) # Proportion of unrecovered microbes
  
  # Mean abundance per microbe group
  pair.df[n,5] <- mean(tmp_twice$mean_Abundance)
  pair.df[n,6] <- mean(tmp_once$mean_Abundance)
  
  # Median abundance per microbe group
  pair.df[n,7] <- median(tmp_twice$mean_Abundance)
  pair.df[n,8] <- median(tmp_once$mean_Abundance)
  
} ; colnames(pair.df) <- c("Pair.ID", "Full_species", "Half_species",  "Ratio", "Mean_twice", "Mean_once", "Med_twice", "Med_once") ; head(pair.df)


# Melting
melt_1_number <- melt(pair.df, id.vars = "Pair.ID", measure.vars = c("Full_species", "Half_species"))
melt_1_number$variable <- ifelse(test = melt_1_number$variable=="Full_species", yes = "Twice", no = "Once")
melt_1_number$variable <- factor(melt_1_number$variable, levels = c("Twice", "Once"))

melt_2_mean <- melt(pair.df, id.vars = "Pair.ID", measure.vars = c("Mean_twice", "Mean_once"))
melt_2_mean$variable <- ifelse(test = melt_2_mean$variable=="Mean_twice", yes = "Twice", no = "Once")
melt_2_mean$variable <- factor(melt_2_mean$variable, levels = c("Twice", "Once"))

melt_3_median <- melt(pair.df, id.vars = "Pair.ID", measure.vars = c("Med_twice", "Med_once"))
melt_3_median$variable <- ifelse(test = melt_3_median$variable=="Med_twice", yes = "Twice", no = "Once")
melt_3_median$variable <- factor(melt_3_median$variable, levels = c("Twice", "Once"))

tiff(filename = "figures/Average_detection_number.tiff", width = 350, height = 600, res = 200, units = "px") ; ggplot(melt_1_number, aes(variable, value)) + 
  xlab("") + ylab("") +
  geom_jitter(aes(color=variable), alpha=0.1) +
  geom_boxplot(aes(fill=variable), alpha=0.4, outlier.alpha = 0) +
  theme_classic2() +
  theme(legend.position = "none",
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(size=rel(1), color="black"),
        axis.text.y = element_text(size=rel(0.8), color="black")) ; dev.off()

tiff(filename = "figures/Average_abundance.tiff", width = 350, height = 600, res = 200, units = "px") ; ggplot(melt_2_mean, aes(variable, log10(value))) + 
  xlab("") + ylab("") +
  geom_jitter(aes(color=variable), alpha=0.1) +
  geom_boxplot(aes(fill=variable), alpha=0.4, outlier.alpha = 0) +
  theme_classic2() +
  theme(legend.position = "none",
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(size=rel(1), color="black"),
        axis.text.y = element_text(size=rel(0.8), color="black")) ; dev.off()

tiff(filename = "figures/Median_abundance.tiff", width = 350, height = 600, res = 200, units = "px") ; ggplot(melt_3_median, aes(variable, log10(value))) + 
  xlab("") + ylab("") +
  geom_jitter(aes(color=variable), alpha=0.1) +
  geom_boxplot(aes(fill=variable), alpha=0.4, outlier.alpha = 0) +
  theme_classic2() +
  theme(legend.position = "none",
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(size=rel(1), color="black"),
        axis.text.y = element_text(size=rel(0.8), color="black")) ; dev.off()

tiff(filename = "figures/Once_among_all_microbe.tiff", width = 200, height = 600, res = 200, units = "px") ; ggplot(pair.df, aes(1, Ratio)) + 
  xlab("") + ylab("")+
  geom_jitter(size=1.2, alpha=0.15, width = 0.3) +
  geom_boxplot(fill="gray20", alpha=0.15, width=0.5) +
  scale_y_continuous(limits = c(0, 1))+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        axis.line = element_line(size=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=rel(0.6), color="black")) ; dev.off()

order1 <- rownames(re)[order(re$mean_Abundance, decreasing = T)]
order2 <- rownames(re)[order(re$NOD, decreasing = T)]
order3 <- rownames(re)[order(re$detected_Abundance, decreasing = T)]

which(order1=="Fusobacterium_nucleatum") # 128
which(order3=="Fusobacterium_nucleatum") # 118

melt1_1 <- melt1 
melt1_1$Matched[melt1_1$Matched=="Undetected"] <- NA
melt1_1$Matched <- factor(melt1_1$Matched, levels = c("Once", "Twice"))
head(melt1_1)

col1 <- c("gray50", "#CC3333", "#3399FF")
tiff("figures/Sup_00_replicated.tiff", width = 1000, height = 600, res = 200, units = "px"); ggplot(melt1_1, aes(Microbe, log10(mean_Abundance))) + 
  geom_point(aes(color=Matched), na.rm = T, alpha=0.2) +
  geom_vline(xintercept = 118, alpha=0.3)+
  scale_color_manual(values = c("#CC3333", "gray50"), name="Detection# \nin a pair", limits=c("Once", "Twice"))+
  xlab("") + ylab("") +
  theme_classic2()+
  scale_x_discrete(limit=order3)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size=0.4),
        legend.key = element_rect(fill="white", color="white"),
        legend.position = c(0.95,0.95), 
        legend.justification = c(0.95,0.95),
        axis.title.x = element_text(size=rel(1.1), vjust=-1.5), 
        legend.spacing.y = unit(0.2, "cm"), # Spacing between legend title and items
        legend.background = element_rect(color = "black", size=0.3),
        legend.title = element_text(size=rel(0.85)), 
        legend.key.height=unit(0.7,"line"), # default = 3 ; spacing between the key items
        legend.text = element_text(size=rel(0.6))) + guides(size=FALSE) ; dev.off()

# Mis-matched pair 
# Matched pair 

df2$Matched <- rowSums(df2=="Twice")
df2$Mis_matched <- rowSums(df2=="Once")
df2$Match_ratio <- round(df2$Matched/(df2$Matched + df2$Mis_matched)*100, 2)
head(df2)
dim(df2)

re$detected_pair <- rowSums(df2[,1:22]!="Undetected")
re$Match_percent <- df2$Match_ratio


which(rownames(re)=="Fusobacterium_nucleatum") # 127

re$Fn <- NA
re$Fn[127] <- "F. nucleatum"

re["Fusobacterium_nucleatum",]
tiff("figures/Sup_00_replicated_Microbes.tiff", width = 1000, height = 600, res = 200, units = "px") ; ggplot(re, aes(log10(detected_Abundance), detected_pair)) + 
  geom_jitter(aes(color=Match_percent), alpha=0.5) +
  #geom_text(aes(label=Fn)) +
  scale_color_gradient(name="Recovery (%)", low="#FF67A4", high = "#00C0AF") + 
  scale_y_continuous(breaks=seq(0, 22,3))+
  #stat_cor(method = "pearson", label.x = -3, label.y = 9, size=3) +
  #stat_smooth(geom = "line", method="loess", alpha=0.5, color="black", se = T) +
  xlab("") + ylab("") +
  theme_classic2()+
  theme(axis.title.x = element_text(size=rel(1.1), vjust=-1.5),
        axis.text = element_text(color="black"),
        axis.line = element_line(size=0.4),
        legend.title = element_text(size=rel(0.7)),
        legend.text.align = 0.5,
        legend.position = c(0.15, 0.7), 
        legend.direction = "vertical",
        legend.title.align = 0.35,
        legend.spacing.y = unit(0.25, "cm"), # Spacing between legend title and items
        legend.background = element_rect(color="black", size=0.3),
        legend.key = element_rect(fill = "white"),
        legend.key.height=unit(1,"line"), # default = 3 ; spacing between the key items
        legend.text = element_text(size=rel(0.6))) ; dev.off()


#### 13) Yogurt or Probiotics with microbes ####

# We tried checking the abundance of Bifidobacterium genus, commonly accepted probiotics, including B. bifidum, B. longum, B. catenulatum, B. dentidum. 
  # Interesting, they had positive correlation with F. nucleatum, and differentially enriched in F. nucleatum-experienced subjects, especially prior to F. nucleatum detection.
  # B. dentium and B. catenulatum were even enrich in IBD condition. 

# But, this result was not included in manuscript. 


# Index
  # 0 = No, I did not consume thse products in the last 7 days
  # 1 = Within the past 4 to 7 days
  # 2 = Within the past 2 to 3 days
  # 4 = Yesterday, 1 to 2 times
  # 5 = Yesterday, 3 or more times

yogurt <- read.csv("input/yogurt.take.csv", header=T) ; head(yogurt)
meta.yrt <- merge(meta, yogurt, all.x = T) ; meta.yrt$Combined <- apply(meta.yrt[,43:44], 1, sum) ; meta.yrt <- merge(meta.yrt, Prob_all, all.x = T) ; meta.yrt <- meta.yrt[order(meta.yrt$External.ID),]; head(meta.yrt)


# Shannon
rcorr(meta.yrt$Shannon, meta.yrt$Combined, type = "spearman") # 0.08 ; P-value = 0.002
rcorr(meta.yrt$Shannon, meta.yrt$Yogurt_index, type = "spearman") # 0.05 ; P-value = 0.0371
rcorr(meta.yrt$Shannon, meta.yrt$Probiotic_index, type = "spearman") # 0.06 ; P-value = 0.0209

# Evenness
rcorr(meta.yrt$Evenness, meta.yrt$Combined, type = "spearman") # 0.06 ; P-value = 0.0242

# Richness
rcorr(meta.yrt$Richness, meta.yrt$Combined, type = "spearman") # 0.06 ; P-value = 0.013

# Human read fraction
rcorr(meta.yrt$Human.Read.Fraction, meta.yrt$Combined, type = "spearman") # 0 ; P-value = 0.9609

# Posterior probability
rcorr(meta.yrt$Post_prob, meta.yrt$Combined, type = "spearman") # -0.02 ; P-value = 0.516


# with Bugs ; espeically Bifidobacterium groups
identical(tbugs.order_ID$External.ID, meta.yrt$External.ID) # TRUE
tbugs.yrt <- merge(tbugs.order_ID, meta.yrt, by = "External.ID") ; dim(tbugs.yrt) # 1526 x 584

tbugs.yrt.cor <- as.matrix(data.frame(tbugs.order, 
                                      Index=meta.yrt$Combined)); head(tbugs.yrt.cor)

Index.cor <- rcorr(tbugs.yrt.cor, type = "spearman")
bugs.index <- data.frame(Spearman=Index.cor$r[,"Index"][1:533],
                         P.value=Index.cor$P[,"Index"][1:533], 
                         FDR=p.adjust(Index.cor$P[,"Index"], method = "fdr", n = 533)[1:533],
                         NOD)
bugs.index <- bugs.index[order(bugs.index$FDR, decreasing = F),] ; bugs.index$microbe <- rownames(bugs.index) # 301 species
bugs.index <- merge(bugs.index, bugs.auc.stat[,c(1, 9:12)], by.y = "microbe") ; head(bugs.index, n = 30)
bugs.index$pre.post <- factor(bugs.index$pre.post, levels = c("Neither", "Prior", "Posterior"))
bugs.index$Condition <- factor(bugs.index$Condition, levels = c("Neither", "nonIBD", "IBD"))
#write.csv(bugs.index, "input/probiotics.ingestion.csv")

# CRC
ggplot(bugs.index, aes(Spearman, -log10(P.value))) + geom_point(aes(color=Core)) +  scale_color_manual(values = col2.1) +
  geom_hline(yintercept = -log10(0.05), color="gray41", linetype="dotted", size=1, alpha=0.5) +  theme_classic() 

# IBD
ggplot(bugs.index, aes(Spearman, -log10(P.value))) + geom_point(aes(color=Condition)) +  scale_color_manual(values = col3.1) +
  geom_hline(yintercept = -log10(0.05), color="gray41", linetype="dotted", size=1, alpha=0.5) +  theme_classic() 

# Classifier
ggplot(bugs.index, aes(Spearman, -log10(P.value))) + geom_point(aes(color=Classifier)) +  scale_color_manual(values = col2.1) +
  geom_hline(yintercept = -log10(0.05), color="gray41", linetype="dotted", size=1, alpha=0.5) +  theme_classic() 

# Pre vs. Post
ggplot(bugs.index, aes(Spearman, -log10(P.value))) + geom_point(aes(color=pre.post)) +  scale_color_manual(values = col3.1) +
  geom_hline(yintercept = -log10(0.05), color="gray41", linetype="dotted", size=1, alpha=0.5) +  theme_classic() 


