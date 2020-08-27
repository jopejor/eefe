#Libraries
library("tidyverse")
library("ggpubr")
library("pheatmap")
library("ComplexHeatmap")
library("FactoMineR")
library("factoextra")
library("UpSetR")
library("pca3d")
library("gplots")
library("cowplot")
library("latex2exp")
library("ggcorrplot")
library("magrittr")

setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/eefe/src/notebooks")
source("../r/eefe_functions.R")

## PCA of biolog data:
#Read and prepare the data
data606 = read.csv("../../data/normalized606Wells86_94RemovedBiologData.csv")
data607 = read.csv("../../data/normalized607Wells86_94RemovedBiologData.csv")
wells = read.csv("../../data/gen_biolog.csv")
wells =  wells %>% filter(ï..Well %in% colnames(data606))
# Read the strategies file
strategy = read.csv("StrategyTable.csv") %>% dplyr::select(Name, Strategy)
data606 = inner_join(data606,strategy, by = "Name") %>% dplyr::select(1:5,Strategy,everything())
data607 = inner_join(data607,strategy, by = "Name") %>% dplyr::select(1:5,Strategy,everything())
colnames(data606)[7:98] = as.character(wells$Assay)
colnames(data607)[7:98] = as.character(wells$Assay)

## 606P
# PCA plot with concentration ellipses
data606$Merged =  paste(data606$Treatment,data606$Condition,sep="") 
data606 = data606 %>% dplyr::select(Merged,everything())
ColFactor =  mutate(data606,Colors = ifelse(Treatment == "Ancestor", Merged, Condition))
ColFactorStrategyColors =  mutate(data606,StrategyColors = Strategy)
# Log transform
dat = data606[,-c(1:7)]
colnames(dat) = wells$Assay
rownames(dat) = make.unique(paste(data606$Name,data606$Condition,sep="-"),sep="-")
logdat = log1p(dat)
# PCA
res.pca = prcomp(logdat, scale = F)
# write scores to potentially use later for biolog restoration analysis
res.pca.pca = PCA(logdat, scale.unit = F, graph = F)
scores =cbind(data606[,c(1:7)],res.pca.pca$ind$coord)
write.csv(scores,"606_scores.csv")

# Talus plot for the informative principal components
eigenvalues = res.pca$eig[,1]
cut = 4
#pdf("Talus_606.pdf")
TalPlot(eigenvalues,cut)
#dev.off()

# PCA plot
#pdf("606_pca.pdf")
p <- fviz_pca_ind(res.pca,
             axes = c(1,2),
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = "black", # color by groups
             palette = c("navy","orange","firebrick","grey53","grey52","grey51"), 
             addEllipses = T, label = "var",
             col.var = "black", repel = TRUE,
             legend.title = "Factor",
             pointshape = 21,
             pointsize = 1.5,
             fill.ind = as.character(ColFactor$Colors),  ellipse.level = 0.8
)
p <- fviz_add(p, res.pca$x[(data606$Strategy=="Generalist" & (data606$Condition=="15" | data606$Condition=="43")),], color = "blue", addlabel = FALSE)
fviz_add(p, res.pca$x[(data606$Strategy=="Specialist" & (data606$Condition=="15" | data606$Condition=="43")),], color = "green", addlabel = FALSE)
#dev.off()

fviz_pca_ind(res.pca)


############################################
## 607P
# PCA plot with concentration ellipses
data607$Merged =  paste(data607$Treatment,data607$Condition,sep="") 
data607 = data607 %>% dplyr::select(Merged,everything())
ColFactor =  mutate(data607,Colors = ifelse(Treatment == "Ancestor", Merged, Condition))
# Log transform
dat = data607[,-c(1:7)]
colnames(dat) = wells$Assay
rownames(dat) = make.unique(paste(data607$Name,data607$Condition,sep="-"),sep="-")
logdat = log1p(dat)
# PCA
res.pca = prcomp(logdat, scale = F)
# write scores to potentially use later for biolog restoration analysis
res.pca.pca = PCA(logdat, scale.unit = F, graph = F)
scores =cbind(data607[,c(1:7)],res.pca.pca$ind$coord)
write.csv(scores,"607_scores.csv")

# Talus plot for the informative principal components
eigenvalues = res.pca$eig[,1]
#cut = 4
#pdf("Talus_607.pdf")
#TalPlot(eigenvalues,cut)
#dev.off()

#pdf("607_pca.pdf")
p <- fviz_pca_ind(res.pca,
                  axes = c(1,2),
                  geom.ind = "point", # show points only (nbut not "text")
                  col.ind = "black", # color by groups
                  palette = c("navy","orange","firebrick","grey53","grey52","grey51"), 
                  addEllipses = T, label = "var",
                  col.var = "black", repel = TRUE,
                  legend.title = "Factor",
                  pointshape = 21,
                  pointsize = 1.5,
                  fill.ind = as.character(ColFactor$Colors),  ellipse.level = 0.7
)
#dev.off()
p <- fviz_add(p, res.pca$x[(data607$Strategy=="Generalist" & (data607$Condition=="15" | data607$Condition=="43")),], color = "blue", addlabel = FALSE)
fviz_add(p, res.pca$x[(data607$Strategy=="Specialist" & (data607$Condition=="15" | data607$Condition=="43")),], color = "green", addlabel = FALSE)



## For PCA of each treaetment separately:
treatmentTemp = "43"
logdat = log1p(dat)[data607$Condition==treatmentTemp,]
# PCA
res.pca = prcomp(logdat, scale = F)
eigenvalues = res.pca$eig[,1]
#pdf("607_pca.pdf")
p <- fviz_pca_ind(res.pca,
                  axes = c(1,2),
                  geom.ind = "point", # show points only (nbut not "text")
                  col.ind = "black", # color by groups
                  palette = c("firebrick","grey53"), 
                  addEllipses = T, label = "var",
                  col.var = "black", repel = TRUE,
                  legend.title = "Factor",
                  pointshape = 21,
                  pointsize = 1.5,
                  fill.ind = as.character(ColFactor$Colors[(ColFactor$Colors == treatmentTemp | ColFactor$Colors == paste("Ancestor",treatmentTemp,sep=""))]),  ellipse.level = 0.7
)
#dev.off()
ind.sub.coord.generalist = data607$Strategy[data607$Condition==treatmentTemp] == "Generalist"
ind.sub.coord.specialist = data607$Strategy[data607$Condition==treatmentTemp] == "Specialist"
p <- fviz_add(p, res.pca$x[ind.sub.coord.generalist,], color = "blue", addlabel = FALSE)
fviz_add(p, res.pca$x[ind.sub.coord.specialist,], color = "green", addlabel = FALSE)


#################################################################################################
#################################################################################################
## Both 606P and 607P strains together:
allData <- rbind(data606,data607)
allData$Merged =  paste(allData$Strain,allData$Treatment,allData$Condition,sep="") 
allData = allData %>% dplyr::select(Merged,everything())
ColFactor =  mutate(allData,Colors = ifelse(Treatment == "Ancestor", Merged, Condition))
# Log transform
dat = allData[,-c(1:7)]
colnames(dat) = wells$Assay
rownames(dat) = make.unique(paste(allData$Name,allData$Condition,sep="-"),sep="-")
logdat = log1p(dat)
# PCA
res.pca = prcomp(logdat, scale = F)
# For writing scores from PCA:
res.pca.pca = PCA(logdat, scale.unit = F, graph = F)
scores =cbind(allData[,c(1:7)],res.pca.pca$ind$coord)
#write.csv(scores,"both606And607_scores.csv")
# Talus plot for the informative principal components
eigenvalues = res.pca.pca$eig[,1]
cut = 4
########## RUN THIS HERE ###########
pdf("Talus_Both606and607.pdf")
TalPlot(eigenvalues,cut)
dev.off()
####################################
#pdf("Both606and607_pca.pdf")
PCsToPlot = c(1,2)
p <- fviz_pca_ind(res.pca,
                  axes = PCsToPlot,
                  geom.ind = "point", # show points only (nbut not "text")
                  col.ind = "black", # color by groups
                  palette = c("navy","orange","firebrick","grey85","grey86","grey87","grey53","grey52","grey51"), 
                  addEllipses = T, label = "var",
                  col.var = "black", repel = TRUE,
                  legend.title = "Factor",
                  pointshape = 21,
                  pointsize = 1.5,
                  fill.ind = as.character(ColFactor$Colors),  ellipse.level = 0.7
)
#p <- fviz_add(p, res.pca$x[(allData$Strain=="606"),], axes = PCsToPlot, color = "blue", addlabel = FALSE)
#fviz_add(p, res.pca$x[(allData$Strain=="607"),], axes = PCsToPlot, color = "green", addlabel = FALSE)
#dev.off()

p <- fviz_add(p, res.pca$x[(allData$Strategy=="Generalist" & (allData$Condition=="43" | allData$Condition=="15")),], axes = PCsToPlot, color = "coral2", addlabel = FALSE)
fviz_add(p, res.pca$x[(allData$Strategy=="Specialist" & (allData$Condition=="43" | allData$Condition=="15")),], axes = PCsToPlot, color = "lightsalmon1", addlabel = FALSE)


p1 = fviz_contrib(res.pca, choice = "var", axes =1, top = 10) + theme(axis.text=element_text(size=6))
p2 = fviz_contrib(res.pca, choice = "var", axes =2, top = 10) + theme(axis.text=element_text(size=6))
p3 = fviz_contrib(res.pca, choice = "var", axes =3, top = 10) + theme(axis.text=element_text(size=6))
p4 = fviz_contrib(res.pca, choice = "var", axes =4, top = 10) + theme(axis.text=element_text(size=6))
#pdf("606_Contributions.pdf")
plot_grid(p1, p2, p3,p4,labels = c('', '',"",""))
#dev.off()
vint1 = p1$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vint2 = p2$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vint3 = p3$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vint4 = p4$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vintmes = unique(c(as.character(vint1$name),as.character(vint2$name),
                   as.character(vint3$name),as.character(vint4$name)))
correlacions = res.pca.pca$var$cor[vintmes,1:4]
# Plot
#pdf("Both606and607_TopCorrelation.pdf")
ggcorrplot(t(correlacions), tl.cex  = 10)
#dev.off()


##############################################################################################################
test2 = dat
annotations = as.data.frame(allData[,1:7])
rownames(annotations) = rownames(dat)
annotations = annotations[,4:7]
annotations$Condition = as.factor(annotations$Condition)
annotations$Strain = as.factor(annotations$Strain)
annotation_colors = list(Treatment = c(Ancestor="gray",Fast = "#00A08A", Random = '#F98400',Slow ='#5BBCD6'),
                         Condition = c(`15`="navy",`37`="orange",`43`="firebrick"),
                         Strategy = c(Generalist = "coral1",Other = "darkgray",Specialist = "bisque1"),
                         Strain = c('606' = "#446455", '607' = "#C7B19C"))

#pdf("Both606and607_heatmap.pdf")
pheatmap(t(test2), scale = "column", annotation_col = annotations, 
               border_color = T,fontsize_row=4,fontsize_col=2, annotation_colors = annotation_colors)
#dev.off()
#save_pheatmap_pdf(hmp, filename = paste("606and607Combined_hmp.pdf",sep=""))# Test heatmap results




## Test heatmap when only pull out one condition:
conditionToUse = "43"
test2 = dat[grepl(conditionToUse, rownames(dat)),]
annotations = as.data.frame(allData[allData$Condition==conditionToUse,1:7])
rownames(annotations) = rownames(test2)
annotations = annotations[,4:7]
annotations$Condition = as.factor(annotations$Condition)
annotations$Strain = as.factor(annotations$Strain)
annotation_colors = list(Treatment = c(Ancestor="gray",Fast = "#00A08A", Random = '#F98400',Slow ='#5BBCD6'),
                         Condition = c(`15`="navy",`37`="orange",`43`="firebrick"),
                         Strategy = c(Generalist = "coral1",Other = "darkgray",Specialist = "bisque1"),
                         Strain = c('606' = "green", '607' = "blue"))
pheatmap(t(test2), scale = "column", annotation_col = annotations, 
         border_color = F,fontsize_row=4,fontsize_col=2, annotation_colors = annotation_colors)




## Test heatmap when only pull out one condition and one lineage:
conditionToUse = "43"
lineageToUse = "607"
dataToUse = allData
test2 = dat[grepl(conditionToUse, rownames(dat)),]
test2 = test2[grepl(lineageToUse, rownames(test2)),]
annotations = as.data.frame(dataToUse[(dataToUse$Condition==conditionToUse & dataToUse$Strain==lineageToUse),1:7])
rownames(annotations) = rownames(test2)
annotations = annotations[,4:7]
annotations$Condition = as.factor(annotations$Condition)
annotations$Strain = as.factor(annotations$Strain)
annotation_colors = list(Treatment = c(Ancestor="gray",Fast = "#00A08A", Random = '#F98400',Slow ='#5BBCD6'),
                         Condition = c(`15`="navy",`37`="orange",`43`="firebrick"),
                         Strategy = c(Generalist = "coral1",Other = "darkgray",Specialist = "bisque1"),
                         Strain = c('606' = "green", '607' = "blue"))
pheatmap(t(test2), scale = "column", annotation_col = annotations, 
         border_color = F,fontsize_row=4,fontsize_col=2, annotation_colors = annotation_colors)




## Test when only use chemical inhibitor wells:
conditionToUse = "37"
#lineageToUse = "607"
dataToUse = allData
chemicalInhibitorWells = c("pH6", "pH5", "1NaCl", "4NaCl", "8NaCl", "1NaLactate","Fusidic Acid","D-Serine","Troleandomycin","Rifamycin","Minocycline","Lincomycin","Guanidine HCl","Niaproof 4","Vancomycin","Tetrazolium Violet","Tetrazolium Blue","Nalidixic Acid","LiCl","Potassium Tellurite","Sodium Butyrate","Sodium Bromate")
chemicalInhibitorData <- dat[,chemicalInhibitorWells]
test2 = chemicalInhibitorData[grepl(conditionToUse, rownames(chemicalInhibitorData)),]
#test2 = test2[grepl(lineageToUse, rownames(test2)),]
#annotations = as.data.frame(dataToUse[(dataToUse$Condition==conditionToUse & dataToUse$Strain==lineageToUse),1:7])
annotations = as.data.frame(dataToUse[(dataToUse$Condition==conditionToUse),1:7])
rownames(annotations) = rownames(test2)
annotations = annotations[,4:7]
annotations$Condition = as.factor(annotations$Condition)
annotations$Strain = as.factor(annotations$Strain)
annotation_colors = list(Treatment = c(Ancestor="gray",Fast = "#00A08A", Random = '#F98400',Slow ='#5BBCD6'),
                         Condition = c(`15`="navy",`37`="orange",`43`="firebrick"),
                         Strategy = c(Generalist = "coral1",Other = "darkgray",Specialist = "bisque1"),
                         Strain = c('606' = "green", '607' = "blue"))
#pdf(paste("606and607Combined_heatmap_chemicalInhibtors",conditionToUse,".pdf",sep=""))
pheatmap(t(test2), scale = "column", annotation_col = annotations, 
         border_color = F,fontsize_row=4,fontsize_col=2, annotation_colors = annotation_colors)
#dev.off()



## Test when only use carbon source wells:
conditionToUse = "37"
#lineageToUse = "606"
dataToUse = allData
chemicalInhibitorWells = c("pH6", "pH5", "1NaCl", "4NaCl", "8NaCl", "1NaLactate","Fusidic Acid","D-Serine","Troleandomycin","Rifamycin","Minocycline","Lincomycin","Guanidine HCl","Niaproof 4","Vancomycin","Tetrazolium Violet","Tetrazolium Blue","Nalidixic Acid","LiCl","Potassium Tellurite","Sodium Butyrate","Sodium Bromate")
carbonSourceData <- dat[, -which(names(dat) %in% chemicalInhibitorWells)]
test2 = carbonSourceData[grepl(conditionToUse, rownames(carbonSourceData)),]
#test2 = test2[grepl(lineageToUse, rownames(test2)),]
#annotations = as.data.frame(dataToUse[(dataToUse$Condition==conditionToUse & dataToUse$Strain==lineageToUse),1:7])
annotations = as.data.frame(dataToUse[(dataToUse$Condition==conditionToUse),1:7])
#annotations = as.data.frame(dataToUse[,1:7])
rownames(annotations) = rownames(test2)
annotations = annotations[,4:7]
annotations$Condition = as.factor(annotations$Condition)
annotations$Strain = as.factor(annotations$Strain)
annotation_colors = list(Treatment = c(Ancestor="gray",Fast = "#00A08A", Random = '#F98400',Slow ='#5BBCD6'),
                         Condition = c(`15`="navy",`37`="orange",`43`="firebrick"),
                         Strategy = c(Generalist = "coral1",Other = "darkgray",Specialist = "bisque1"),
                         Strain = c('606' = "green", '607' = "blue"))
#pdf(paste("606and607Combined_heatmap_caronSources",conditionToUse,".pdf",sep=""))
pheatmap(t(test2), scale = "column", annotation_col = annotations, 
         border_color = F,fontsize_row=4,fontsize_col=2, annotation_colors = annotation_colors)

#dev.off()
