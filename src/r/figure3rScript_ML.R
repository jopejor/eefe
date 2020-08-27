setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/eefe/src/notebooks")
# libraries
library(tidyverse)
library(multcomp)
library(ggpubr)
library(latex2exp)
#library(broom)
#library(purrr)
#library(car)
#library(data.table)
library(formattable)
library(cowplot)

# Read the data
ancestral =  read.csv("../../data/growth_ancestrals.csv")
head(ancestral,n=2)

# For each experiment
exp43 = ancestral %>% filter(Temperature == 43)
exp15 = ancestral %>% filter(Temperature == 15)
exp37 = ancestral %>% filter(Temperature == 37)

my_comparisons <- list(c("606P","607P") ,c("REL606","REL607"),c("606P","REL606"),c("607P","REL607"))

## Growth at 43C
# Plot
#pdf("FigureS1_A.pdf")
ggboxplot(exp43, x = "Strain", y = "Max_deriv",
          color = "Strain", palette = "jco",add="jitter")+
  stat_compare_means(comparisons = my_comparisons,method ="t.test",
                     method.args = list(alternative = "greater"),
                     p.adjust.method = "fdr",var.equal=F,
                     label="p.signif")+ylab("Maximum Growth Rate (OD/min)")+
  font("x.text", size = 11)+font("y.text", size = 11)
#dev.off()

## Growth at 15C
# Plot
#pdf("FigureS1_C.pdf")
ggboxplot(exp15, x = "Strain", y = "Max_deriv",
          color = "Strain", palette = "jco",add="jitter")+
  stat_compare_means(comparisons = my_comparisons,method ="t.test",p.adjust.method = "fdr",var.equal
                     =F,label="p.signif",method.args = list(alternative = "greater"))+ylab("Maximum Growth Rate (OD/min)")+ 
  scale_y_continuous(labels = function(x) format(x, scientific = F))+
  font("x.text", size = 11)+font("y.text", size = 11)
#dev.off()

## Growth at 37C
# Plot
#pdf("FigureS1_B.pdf")
ggboxplot(exp37, x = "Strain", y = "Max_deriv",
          color = "Strain", palette = "jco",add="jitter")+
  stat_compare_means(comparisons = my_comparisons,method ="t.test",
                     method.args = list(alternative = "greater"),
                     p.adjust.method = "fdr",var.equal
                     =F,label="p.signif")+ylab("Maximum Growth Rate (OD/min)")+ 
  scale_y_continuous(labels = function(x) format(x, scientific = F)) +  
  font("x.text", size = 11)+font("y.text", size = 11)
#dev.off()



#######################################################################
## Figure S2:
# Read the data
dependent =  read.csv("../../data/dependent_15.csv")
dependent$X = factor(paste(dependent$Strain, dependent$Type,sep="_"))
head(dependent,n=2)
table(dependent$X)
levels(dependent$X)
comp_dep = list(c("607P_Acclimation","607P_Condition"), 
                c("606P_Acclimation","606P_Condition"),
                c("F606-2_Acclimation","F606-2_Condition"),
                c("S606-2_Acclimation","S606-2_Condition"),
                c("S607-1_Acclimation","S607-1_Condition"),
                c("S607-2_Acclimation","S607-2_Condition"))
#pdf("FigureS2_dependent_15.pdf")
ggboxplot(dependent, x = "X", y = "Max_deriv",
          color = "Strain", palette = "jco",add="jitter")+
  stat_compare_means(comparisons = comp_dep,method ="t.test",
                     method.args = list(alternative = "greater"),
                     p.adjust.method = "fdr",var.equal
                     =F,label="p.signif")+ylab("Maximum Growth Rate (OD/min)")+ 
  xlab("Strain")+
  scale_y_continuous(labels = function(x) format(x, scientific = F)) +  
  font("x.text", size = 11)+font("y.text", size = 11)+
  theme(axis.text.x=element_text(angle=45, hjust=1))
#dev.off()



##########################################################
## Table 1:
# Find outliers function:
FindOutliers <- function(data) {
  lowerq = quantile(data)[2]
  upperq = quantile(data)[4]
  iqr = upperq - lowerq #Or use IQR(data)
  # mild outliers:
  mild.threshold.upper = (iqr * 1.5) + upperq
  mild.threshold.lower = lowerq - (iqr * 1.5)
  # we identify extreme outliers
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)
  result <- which(data > mild.threshold.upper | data < mild.threshold.lower)
}

#### Testing
clones =  read.csv('../../data/growth_strains_new.csv')
head(clones,n=5)


# Calculate the average
sum_clones = clones %>% group_by(Linage,Strain,Treatment,Experiment,Temperature,Type)  %>%
  summarise(Median = median(Max_deriv), Std = sd(Max_deriv), Cv = sd(Max_deriv)/median(Max_deriv), n= n() ) %>% ungroup()
sum_clones$Experiment = as.character(sum_clones$Experiment)
head(sum_clones)
# Significant
clones = clones %>% distinct() %>% mutate(Group  = factor(paste(Strain,Experiment,Temperature,sep ="_")))
references = clones %>% filter(Type == "Reference")
cloners = clones %>% filter(Type == "Clone") %>% droplevels()
Strain = factor(unique(cloners$Group))
significant = data.frame(vars=double(),Strain=factor(),Temperature=integer())
tempTest = rep(NA, length(factor(unique(cloners$Group))));
tempRef = rep(NA, length(factor(unique(cloners$Group))));
sdTest = rep(NA, length(factor(unique(cloners$Group))));
sdRef = rep(NA, length(factor(unique(cloners$Group))));

for (i in 1:length(factor(unique(cloners$Group)))){
  strain = Strain[i]
  test = cloners %>% filter(Group == strain)
  ref = references %>% filter(Experiment == test$Experiment[1] & Type == "Reference" & Linage ==  test$Linage[1])
  # Find any outliers and remove:
  # use the function to identify outliers
  #tempTest[i] <- FindOutliers(test$Max_deriv)
  #tempRef[i] <- FindOutliers(ref$Max_deriv)
  sdTest[i] <- sd(test$Max_deriv)
  sdRef[i] <- sd(ref$Max_deriv)
}

#### End testing

# Read the data
clones =  read.csv('../../data/growth_strains_new.csv')
head(clones,n=5)


# Calculate the average
sum_clones = clones %>% group_by(Linage,Strain,Treatment,Experiment,Temperature,Type)  %>%
  summarise(Median = median(Max_deriv), Std = sd(Max_deriv), Cv = sd(Max_deriv)/median(Max_deriv), n= n() ) %>% ungroup()
sum_clones$Experiment = as.character(sum_clones$Experiment)
head(sum_clones)
# Significant
clones = clones %>% distinct() %>% mutate(Group  = factor(paste(Strain,Experiment,Temperature,sep ="_")))
references = clones %>% filter(Type == "Reference")
cloners = clones %>% filter(Type == "Clone") %>% droplevels()
Strain = factor(unique(cloners$Group))
significant = data.frame(vars=double(),Strain=factor(),Temperature=integer())
for (i in 1:length(factor(unique(cloners$Group)))){
  strain = Strain[i]
  test = cloners %>% filter(Group == strain)
  ref = references %>% filter(Experiment == test$Experiment[1] & Type == "Reference" & Linage ==  test$Linage[1])
  # Find any outliers and remove:
  # use the function to identify outliers
  tempTest <- FindOutliers(test$Max_deriv)
  tempRef <- FindOutliers(ref$Max_deriv)
  # remove the outliers
  if (length(tempTest) > 0) {
    test <- test[-tempTest,]
    print("i = ", i)
  }
  if (length(tempRef > 0)) {
    ref <- ref[-tempRef]
    print("ref i = ", i)
  }
  #
  d = rbind(ref,test)
  # Sanity plot
  filename = paste(strain,test$Temperature[1],sep="_")
  pdf(filename)
  boxplot(d$Max_deriv~d$Type)
  dev.off()
  # t.test
  var_q = bartlett.test(Max_deriv~Type,d)
  var_T = (var_q$p.value>0.05)
  sig = d %>% summarise_each(funs(t.test(.[Type == "Clone"], .[Type == "Reference"],
                                         var.equal=var_T,alternative="greater")$p.value), 
                             vars = Max_deriv)
  sig$Strain = test$Strain[1]
  sig$Temperature = test$Temperature[1]
  significant =  rbind(significant,sig)
}
# Reshape
Growth = significant %>% distinct() %>% spread(Temperature, vars)
Growth$Strain = factor(Growth$Strain, levels = c("R606-1", "R606-2", "R606-3", "R606-4","F606-1", "F606-2", "F606-3", "F606-4", "S606-1", "S606-2", "S606-3", "S606-4",
                                                 "R607-1", "R607-2", "R607-3", "R607-4","F607-1", "F607-2", "F607-3", "F607-4", "S607-1", "S607-2", "S607-3", "S607-4"))
Growth = Growth %>% arrange(Strain)
# Replicate table 1
sum_clones$RelFit = 0
for (i in 1:nrow(sum_clones)){
  row_df = sum_clones[i,]
  reference =sum_clones %>% 
    filter(Experiment == as.character(row_df["Experiment"]) & 
             Type == "Reference" & 
             Temperature == as.numeric(row_df["Temperature"]) & 
             Linage == as.character(row_df["Linage"]) )
  sum_clones$RelFit[i] = row_df$Median/reference$Median
  
  
}
fitness_table = sum_clones %>% arrange(Linage, Temperature) %>% filter(Type == "Clone")
table1 = fitness_table %>% dplyr::select(Strain,Temperature,RelFit,Treatment,Linage)  %>% distinct() %>% spread(Temperature, RelFit)
table1$Strain = factor(table1$Strain, levels = c("R606-1", "R606-2", "R606-3", "R606-4","F606-1", "F606-2", "F606-3", "F606-4", "S606-1", "S606-2", "S606-3", "S606-4",
                                                 "R607-1", "R607-2", "R607-3", "R607-4","F607-1", "F607-2", "F607-3", "F607-4", "S607-1", "S607-2", "S607-3", "S607-4"))
table1 = table1 %>% arrange(Strain)
colnames(Growth) = c("Strain","Sig15","Sig43")
# Add the pvalues
table1 = full_join(table1,Growth, by = "Strain")
# We will assume that a not significant pvalue implies a ration of 1
table2 = table1
# First deal with the internalized cases
table2$Sig15[table2$Strain == "S607-2"] = (1- table1$Sig15[table1$Strain == "S607-2"])
table2$Sig15[table2$Strain == "S606-2"] = (1- table1$Sig15[table1$Strain == "S606-2"])
table2$`15`[table2$Sig15 >0.05]  = 1
table2$`43`[table2$Sig43 >0.05]  = 1
# Signs
table2$`15` = -1*table2$`15`
## Figure 3 B:
# Classify accordingly
StrategyTable = table2 %>% mutate(., Strategy = if_else(`15`<(-1) & `43`>1,
                                                        "Generalist",if_else(`15`<(-1) & `43`<=1,"Specialist",
                                                                             if_else(`15`>=(-1) & `43`>1,"Specialist",
                                                                                     "Other"))))
# Add ancestral data to strategy table:
ancestor607Df = data.frame("607P","Ancestor",607,0,0,0,0,"Other")
names(ancestor607Df) = colnames(StrategyTable)
ancestor606Df = data.frame("606P","Ancestor",606,0,0,0,0,"Other")
names(ancestor606Df) = colnames(StrategyTable)
StrategyTableToSave = rbind(ancestor606Df,ancestor607Df,StrategyTable)
colnames(StrategyTableToSave)[1] = "Name"
write.csv(StrategyTableToSave,"StrategyTable.csv")
#


SigTable = gather(StrategyTable,Temperature,Fitness,c(`15`,`43`))

# create color vector for generalist, specialist and other labels:
StrategyTable = StrategyTable %>% mutate(., StrategyCol = if_else(`15`<(-1) & `43`>1,
                                                        "coral1",if_else(`15`<(-1) & `43`<=1,"bisque2",
                                                                             if_else(`15`>=(-1) & `43`>1,"bisque2",
                                                                                     "darkgray"))))

## Figure 3 A:
#pdf("Figure3_A.pdf")
ggdotchart(SigTable, x = "Strain", y = "Fitness",
           color = "Treatment",
           fill = "Temperature",
           # Color by groups
           palette = c("#00A08A",'#F98400','#5BBCD6'), # Custom color palette
           #sorting = "asc",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "Strain",                                # Order by groups
           dot.size = 7,                                 # Large dot size
           label = abs(round(SigTable$Fitness,1)),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 10, 
                             vjust = 0.4),               # Adjust label parameters
           ggtheme = theme_pubr(),
           ylim = c(-1.8,1.8),# ggplot2 theme,
           
)+
  geom_hline(yintercept = 1, linetype = 2, color = "firebrick",size=1.5,alpha=0.7)+
  geom_hline(yintercept = 0, linetype = 1, color = "black",size=0.5,alpha=1)+
  geom_hline(yintercept = -1, linetype = 2, color = "navy",size=1.5,alpha=0.7)+
  theme(axis.text.x=element_text(angle=45, hjust=1, color = StrategyTable$StrategyCol))+labs(x="")
#dev.off()


## Figure 3 B:
#pdf("Figure3B.pdf")
ggplot(StrategyTable) + aes(x = Treatment, fill =Strategy)+
  geom_bar(position = "fill")+ 
  scale_fill_manual(values= c("coral1","darkgray","bisque2") )+ ylab("Relative Abundance") + 
  xlab("Environment")+facet_wrap("Linage")
#dev.off()