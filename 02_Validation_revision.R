# Revision 

# Validation of prediction model into additional healthy individuals



#### 1. Setting ####

## Requirements
options(java.parameters = "-Xmx64g", stringsAsFactors = F) 
setwd(dir = "D:/0_Previous/2018-20_LSG/Data/Microbiome_LSG/2019/190814 Figure_1526 samples/")

library(vegan) ;library(MASS) ; library(pROC) ; library(caret) ; library(RJSONIO)

rev.bugs <- read.csv("D:/0_Previous/2018-20_LSG/Publication/Paper writing/200526/revision_bugs.csv", row.names = 1)
  rev.log <- log10(rev.bugs+1e-5) 
  
mark <- read.csv("input/lefse.marker.csv", row.names = 1) # LEfSE results; contains 14 non-IBD + 12 IBD marker species 
  ibd.mark <- rownames(mark)[mark$Class=="IBD"]
  non.mark <- rownames(mark)[mark$Class=="nonIBD"]
  mark.list <- c(ibd.mark, non.mark)

core_CRC <- read.csv("input/core.species.CRC.csv") ; core_CRC <- as.character(core_CRC$core.CRC)  # Core signature of CRC

lef <- read.csv("input/lefse.marker.csv") ; lef_order <- lef$Microbes[order(lef$Class, lef$Log10_LDA, decreasing = T)] ; head(lef) # Import LEfSe results file





##### Modeling ####

bugs <- read.csv("input/bugs_1526.csv", header = T, row.names = 1) 
  bugs <- log10(bugs+1e-5) 
  # Log transformation with pseudo-abundance
  # 0.00001 is the half of the minimum abundance across whole samples
  
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


#------------- 1) Model construction
# 1-1) Sorting metadata of F.nucleatum-exposed samples
meta.exp <- meta.order[meta.order$Fn.pre.post.simple!="Negative",]
meta.exp$Fn.pre.post.simple <- droplevels(meta.exp$Fn.pre.post.simple)


# 1-2) Sorting microbial abundance data of F.nucleautm-exposed samples 
bugs.exp <- tbugs.order[meta.order$Fn.pre.post.simple!="Negative",]
bugs.exp$Fn.pre.post.simple <- droplevels(bugs.exp$Fn.pre.post.simple)
dim(bugs.exp) # 317 x 535


# Model #10 was the best performer
auc.top10 <- auc.sort[1:10,] # list of classifiers
list10 <- auc.top10$microbe 
data10 <- bugs.exp[,c(list10, "Condition", "Fn.pre.post.simple")] # Dataset ;F.n.-experienced


#------------- 2) Validation healthy control 

rev.log <- as.data.frame(t(rev.log)) # 82 samples x 372 species
rev.log$Condition <- "nonIBD"

# Control & Data partitioning
set.seed(29382)
trCtrl <- trainControl(method = "repeatedcv",  # Cross-validation
                       number=5, 
                       repeats=10,
                       classProbs = T, 
                       summaryFunction = twoClassSummary) # twoClassSummary gives AUC

idx <- createDataPartition(data10$Fn.pre.post.simple, 
                           times = 100, p = 0.7, list = F)


#------------Generating modeling function--------------#


# 1) Container for posterior probability of F.nucleatum-"exposed" samples
mod.posterior <- as.data.frame(matrix(nrow = 317, ncol = 100)) 
  rownames(mod.posterior) <- rownames(data10) 
  # sampleID x (microbes + feature)
  
# 3) Container for "validation set"
mod.rev <- as.data.frame(matrix(nrow = 82, ncol = 100)) 



for (i in 1:100){ # Repeated 100 times with random partitioned dataset
  
# Partitioning 100-times
tmp.idx <- idx[,i] ; train <- data10[tmp.idx,] ; test <- data10[-tmp.idx,]
    
# Logistic regression 
glm <- train(form=Fn.pre.post.simple~.,
             data=train,
             method="glm",
             trControl=trCtrl,
             family="binomial",
             metric="ROC")
    
# Posterior Probability
prd.test <- predict(glm, newdata = test, type = "prob") # F.n-experienced 7:3
mod.posterior[rownames(test), i] <- prd.test$Posterior # F.n-exposed samples

# Validatino set
prd.rev <- predict(glm, newdata = rev.log, type = "prob")
mod.rev[,i] <- prd.rev$Posterior # Healthy validation set
    

} 
  
# Calculate mean Posterior probability for Validation dataset
mod.rev$meanPost.prob <- rowMeans(mod.rev[,c(1:100)])
mod.rev$pred.class <- ifelse(test = mod.rev$meanPost.prob>0.5, 
                             "Posterior", "Prior")
  
#write.csv(mod.rev, "output/Revision_validation_HMP_healthy.csv")


# PCoA and Posterior probability
rev.input <- (rev.log[,-c(369:373)])+5
bugs.dist <- vegdist(rev.input, method = "bray")
bugs.pca <- cmdscale(bugs.dist, eig = T, k = 2)

rev.pca <- data.frame(bugs.pca$points) ; colnames(rev.pca) <- c("PC1", "PC2")
rev.var <- round(bugs.pca$eig^2/sum(bugs.pca$eig^2)*100,1)
  mod.rev$SampleID <- rownames(rev.pca)
  rev.summary <- cbind(mod.rev, rev.pca)

tiff("figures/revision_01_PCoA.tiff", res = 180,
     units = "px", width = 600, height = 600); ggplot(rev.summary, 
                                                      aes(PC1, PC2)) +
  xlab(paste0("PC1 (", rev.var[1], ")%")) + 
  ylab(paste0("PC2 (", rev.var[2], ")%")) +
  geom_point(alpha=0.8, aes(color=meanPost.prob)) +
  scale_color_gradient2( midpoint = 0.5, 
                         low = "navyblue", mid = "white", high = "red")+
  coord_fixed()+
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(vjust=-1.5, size=rel(1)), 
        axis.title.y = element_text(vjust=1.7, size=rel(1)), 
        #legend.position = "bottom",
        legend.justification = c(1,0),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.6)),
        legend.key = element_rect(fill = "white"),
        legend.key.width = unit(0.65, "line"),
        legend.key.height = unit(0.85,"line"), 
        legend.spacing.y = unit(0, "mm"), # Spacing between legend title and items
        legend.spacing.x = unit(0.45, "mm")) ; dev.off()

tiff("figures/revision_01_PCoA_CAT.tiff", res = 180,
     units = "px", width = 600, height = 600); ggplot(rev.summary, 
                                                      aes(PC1, PC2)) +
  xlab(paste0("PC1 (", rev.var[1], ")%")) + 
  ylab(paste0("PC2 (", rev.var[2], ")%")) +
  geom_point(alpha=0.8, aes(color=pred.class)) +
  coord_fixed()+
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(vjust=-1.5, size=rel(1)), 
        axis.title.y = element_text(vjust=1.7, size=rel(1)), 
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.8)),
        legend.key = element_rect(fill = "white"),
        legend.key.width = unit(0.8, "line"),
        legend.key.height = unit(0.8,"line"), 
        legend.spacing.y = unit(0, "mm"), # Spacing between legend title and items
        legend.spacing.x = unit(0.6, "mm")) ; dev.off()

rownames(rev.bugs)

rev.input <- t(rev.bugs[-which(rownames(rev.bugs) %in% c("unclassified", "sum", "Richness", "Bacterial.abundance")),
                      ]) ; tail(rev.input)
rev.div <- data.frame(Shannon=diversity(rev.input),
                      Richness=as.data.frame(t(rev.bugs["Richness",]))$Richness)
rev.div$Evenness <- rev.div$Shannon/(log(rev.div$Richness))

rev.div$IBD.num <- rowSums(rev.input[,ibd.mark]>0)
rev.div$NON.num <- rowSums(rev.input[,non.mark]>0)

core_CRC_detected <- core_CRC[core_CRC %in% colnames(rev.input)]
rev.div$CRC <- rowSums(rev.input[,core_CRC_detected]>0)
rev.div$IBD.NON.ratio <- (rev.div$IBD.num+1)/(rev.div$NON.num+1)
rev.div$PostProb <- rev.summary$meanPost.prob
rev.div$Predicted <- rev.summary$pred.class 
colnames(rev.div)[4:7] <- c("IBD", "nonIBD", "CRC", "IBD/nonIBD")
  write.csv(rev.div, "input/revision_summary.csv")

rev.melt <- reshape::melt(rev.div, id=c("Predicted", "PostProb"))

head(rev.melt)

rev.melt$variable <- factor(rev.melt$variable, 
                            levels=c("Richness", "Evenness", "Shannon",
                                     "IBD", "nonIBD", "IBD/nonIBD",
                                     "CRC"))
rev.melt$Predicted <- factor(rev.melt$Predicted, 
                            levels=c("Prior", "Posterior"))

tiff("figures/revision_Post_prob_vs_diversity2.tiff", 
     res=160, units = "px",
     width = 600, height = 600);ggplot(rev.div, 
                                        aes(PostProb, Shannon)) +
  xlab("Posterior probability") + ylab("Shannon diveristy")+
  geom_point(alpha=0.2)+
  stat_cor(method = "spearman", label.x.npc = 0.01, label.y.npc = 0.05, 
           size=5.5)+
  geom_smooth(method = "lm", span=0.95)+
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=rel(1.1), color="black"),
        axis.title = element_text(size=rel(1.5)));dev.off()

tiff("figures/revision_02_summary2.tiff", 
     res=150, units = "px",
     width = 1600, height = 700);ggplot(rev.melt, 
                                        aes(Predicted, value, fill=Predicted)) +
  xlab("") + ylab("")+
  geom_violin(alpha=0.7, draw_quantiles = 0.5) +
  scale_fill_manual(values = c("#3399FF","#CC3333")) +
  stat_compare_means(comparisons = list(c("Prior","Posterior")),
                     method = "wilcox.test", tip.length = 0, size=3)+
  theme_bw() +
  facet_wrap(~variable, nrow = 1, scales = "free_y")+
  theme(axis.text.x = element_text(color="black", size=rel(1.25)),
        axis.text.y = element_text(color="black", size=rel(0.9),
                                   angle=90, hjust=0.5),
        legend.position = "none",
        strip.text = element_text(size=rel(1.4)));dev.off()


# Cluster 1 & 6

clu <- read.csv("input/cluster.csv")

good.b <- clu$microbe[clu$Cluster==1]
bad.b <-clu$microbe[clu$Cluster==6]

cl.info <- data.frame(C1_sum=colSums(rev.bugs[good.b,]),
                      C6_sum=colSums(rev.bugs[bad.b,], na.rm = T),
                      C1_pre=colSums(rev.bugs[good.b,]>0),
                      C6_pre=colSums(rev.bugs[bad.b,]>0, na.rm = T))

identical(rownames(rev.pca), rownames(cl.info)) # TRUE

rev.pca <- cbind(rev.pca, cl.info)
head(rev.pca)
colnames(rev.pca)[5:6] <- c("Cluster 1", "Cluster 6")
rev.melt <- reshape::melt(rev.pca, id=c("PC1", "PC2"), measure=c("Cluster 1", "Cluster 6"))

head(rev.melt)
tiff("figures/revision_02_Cluster_prevelance.tiff", res = 180,
     units = "px", width = 1200, height = 400); ggplot(rev.melt, 
                                                      aes(PC1, PC2)) +
  xlab(paste0("PC1 (", rev.var[1], ")%")) +
  ylab(paste0("PC2 (", rev.var[2], ")%")) +
  geom_point(alpha=0.8, aes(color=value)) +
  scale_color_gradient(low = "black", high = "yellow")+
  coord_fixed()+
  facet_wrap(~variable)+
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(vjust=-1.5, size=rel(1)), 
        axis.title.y = element_text(vjust=1.7, size=rel(1)), 
        #legend.position = "none",
        legend.justification = c(1,0),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.6)),
        legend.key = element_rect(fill = "white"),
        legend.key.width = unit(0.65, "line"),
        legend.key.height = unit(0.85,"line"), 
        legend.spacing.y = unit(0, "mm"), # Spacing between legend title and items
        legend.spacing.x = unit(0.45, "mm")) ; dev.off()
