#Libraries
library("tidyverse")
library("ggpubr")
library("ComplexHeatmap")

setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/eefe/src/notebooks")
source("../r/eefe_functions.R")

#Read and prepare the data
data=read.csv("../../data/mutations_populations.csv")
data = data %>% arrange(Population)
data = data %>% mutate(Strain = str_extract(Population,'[0-9][0-9][0-9]'))
head(data)

# Create summary 
resum = as.data.frame(table(data$Population))
resum$Environment  = resum$Var1
resum$Lineage = resum$Var1
#
colnames(resum)[1]="Population"
colnames(resum)[2]="Mutations"
#Recode
resum$Environment = sub(resum$Environment,pattern="F.*",replacement = "Fast")
resum$Environment = sub(resum$Environment,pattern="R.*",replacement = "Random")
resum$Environment = sub(resum$Environment,pattern="S.*",replacement = "Slow")
resum$Lineage = sub(resum$Lineage,pattern="*.606.*",replacement = "606")
resum$Lineage = sub(resum$Lineage,pattern="*.607.*",replacement = "607")
resum$Both = paste(resum$Environment,resum$Lineage,sep=" ")
head(resum)

compare_means(Mutations ~ Both,  data = resum,method = "t.test")

#pdf("NumberMutations_boxplot.pdf")
ggboxplot(resum, x = "Both", y = "Mutations",
          color = "Environment", palette = c("#00A08A",'#F98400','#5BBCD6')
          ,add="jitter")+
  ylab("Number of Mutations")+
  xlab("Linage+Environment")+
  font("x.text", size = 11)+font("y.text", size = 11)+theme_xp()
#dev.off()

# Add Environment
data$Environment = data$Population
data$Environment = sub(data$Environment,pattern="F.*",replacement = "Fast")
data$Environment = sub(data$Environment,pattern="R.*",replacement = "Random")
data$Environment = sub(data$Environment,pattern="S.*",replacement = "Slow")
head(data)

# Plot by mutation class
data$Class = factor(data$Class, levels = c("Non_Synonymous", "Small_indel", "Intergenic_snp",
                                           'Stop','Synonymous','Large_deletion',"Large_amplification"))
data$ClassPlot = data$Class

#pdf("BarplotMutationTypes.pdf")
ggplot(data) + aes(x = Environment, fill = Class)+
  geom_bar(position = "fill") + 
  scale_fill_manual(values= c("#A6CEE3", "#FB9A99",'#FF7F00',
                              '#6A3D9A','#FDBF6F','#E31A1C',"#1F78B4"),
                    labels= c("Nonsynonymous", "Small indel", "Intergenic snv",
                              'Stop','Synonymous','Large deletion',"Large amplification")
  ) + ylab("Relative Abundance") + 
  xlab("Environment")+facet_wrap("Strain")+
  theme_xp()
#dev.off()

# Plot by mutation categories:
data$Category = factor(data$Category, levels = c("Metabolic","Heat_Shock","Membrane_and_Cell_Wall",
                                                 "Other","Regulatory","Cold_Shock","Transporter"))
data$CategoryPlot = data$Category
#pdf("BarplotMutationCategories.pdf")
ggplot(data) + aes(x = Environment, fill = Category)+
  geom_bar(position = "fill") + 
  scale_fill_manual(values= c("cadetblue2", "brown1",'coral',
                              'blueviolet','gold',"dodgerblue",'hotpink'),
                    labels= c("Metabolic","Heat Shock","Membrane and Cell Wall",
                              "Other","Regulatory","Cold Shock","Transporter")
  ) + ylab("Relative Abundance") + 
  xlab("Environment")+facet_wrap("Strain")+
  theme_xp()
#dev.off()

### Oncoprint:
## Prepare the input
prepared = data %>% group_by(.dots=c("Population","Gene")) %>% mutate(collapsed=paste(Class, collapse = ';'))
prepared = prepared %>% select(Gene, Population, collapsed)
prepared = prepared %>% distinct()
prepmat = prepared %>% spread(Population,collapsed,fill="",drop=FALSE)
## Prepare the matrix
mat=as.matrix(prepmat)
rownames(mat) = mat[,1]
mat = mat[, -1]
#Order the matrix
logimat=!(mat == "")
#mat=mat[order(rowSums(logimat),decreasing=T),]
#Save the mutation
prepared$state="1"
mut_mat = prepared %>% spread(Population,state,fill="0",drop=T)
head(mut_mat)

#Get the populations information
pheno = read.csv("../../data/sample_annotations.csv")
Environment=pheno$Environment
Linage = pheno$Linage
ha = HeatmapAnnotation(Environment = Environment, Linage= Linage,
                       col = list(Linage = c("606"='#446455',"607"='#C7B19C'),
                                  Environment = c("Fast" = "#00A08A", "Random" = '#F98400',
                                                  "Slow"='#5BBCD6')
                       ),
                       annotation_height = unit(c(5, 5, 15), "mm"),
                       annotation_legend_param = list(legend_position = "bottom",Environment = list(title = "Environment"),Linage = list(title = "Linage")))

# Orders
StrainOrder = c(grep(colnames(mat), pattern = "R607"),
                grep(colnames(mat), pattern = "F607"),
                grep(colnames(mat), pattern = "S607"),
                grep(colnames(mat), pattern = "R606"),
                grep(colnames(mat), pattern = "F606"),
                grep(colnames(mat), pattern = "S606"))

#OncoPrint
#pdf("oncoPrintFigure2.pdf")
dev.new()
oncoPrint(mat, column_order = StrainOrder,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE,
          row_names_gp = gpar(fontsize = 7),
          heatmap_legend_param = list(title = "Mutation Classes", 
                                      at = c("Intergenic_snp",'Synonymous',"Non_Synonymous", 'Stop', "Small_indel" ,
                                             'Large_deletion',"Large_amplification") , 
                                      labels = c("Intergenic SNP",'Synonymous',
                                                 "Nonsynonymous",'Stop', 
                                                 "Small Indels",'Large Deletion',
                                                 'Large Amplification')),
          bottom_annotation = ha, pct_gp = gpar(fontsize = 7),
          split = sample(c("606","607"), nrow(mat), replace = TRUE)
          
)
#dev.off()



#####################################
## Oncoprint for category of mutation:
# Prepare the input
prepared = data %>% group_by(.dots=c("Population","Gene")) %>% mutate(collapsed=paste(Category, collapse = ';'))
prepared = prepared %>% select(Gene, Population, collapsed)
prepared = prepared %>% distinct()
prepmat = prepared %>% spread(Population,collapsed,fill="",drop=FALSE)
## Prepare the matrix
mat=as.matrix(prepmat)
rownames(mat) = mat[,1]
mat = mat[, -1]
#Order the matrix
logimat=!(mat == "")
#mat=mat[order(rowSums(logimat),decreasing=T),]
#Save the mutation
prepared$state="1"
mut_mat = prepared %>% spread(Population,state,fill="0",drop=T)
head(mut_mat)

#Get the populations information
pheno = read.csv("../../data/sample_annotations.csv")
Environment=pheno$Environment
Linage = pheno$Linage
ha = HeatmapAnnotation(Environment = Environment, Linage= Linage,
                       col = list(Linage = c("606"='#446455',"607"='#C7B19C'),
                                  Environment = c("Fast" = "#00A08A", "Random" = '#F98400',
                                                  "Slow"='#5BBCD6')
                       ),
                       annotation_height = unit(c(5, 5, 15), "mm"),
                       annotation_legend_param = list(legend_position = "bottom",Environment = list(title = "Environment"),Linage = list(title = "Linage")))

# Orders
StrainOrder = c(grep(colnames(mat), pattern = "R607"),
                grep(colnames(mat), pattern = "F607"),
                grep(colnames(mat), pattern = "S607"),
                grep(colnames(mat), pattern = "R606"),
                grep(colnames(mat), pattern = "F606"),
                grep(colnames(mat), pattern = "S606"))

#OncoPrint
#pdf("oncoPrintFigure2.pdf")
oncoPrint(mat, column_order = StrainOrder,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE,
          row_names_gp = gpar(fontsize = 7),
          heatmap_legend_param = list(title = "Mutation Categories", 
                                      at = c("Metabolic","Heat_Shock","Membrane_and_Cell_Wall",
                                             "Other","Regulatory","Cold_Shock","Transporter") , 
                                      labels = c("Metabolic","Heat Shock","Membrane and Cell Wall",
                                                 "Other","Regulatory","Cold Shock","Transporter")),
          bottom_annotation = ha, pct_gp = gpar(fontsize = 7),
          split = sample(c("606","607"), nrow(mat), replace = TRUE)
          
)
#dev.off()

