
#Libraries
library('tidyverse')
library("rlang")

setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/eefe/src/notebooks")

#Functions
compare_to_ancestral=function(data,T_ref, T_test,refname){
  #Empty dataframe
  comparison = data.frame(clone = character(), median=numeric(),pval=numeric(),logistest=character(),
                          referencia=numeric(),PC=character())
  #Set the two datasets to be compared
  datatest=subset(data,data$Condition==T_test) #Test temperature
  dataref=subset(data,data$Condition==T_ref) #Reference temperature    
  pcs = colnames(data)[9:12]
  #Loop
  for (i in 1:length(pcs)) {
    pcs_t = pcs[i]
    #Reference medians
    median_PC = datatest %>% select(clone,Condition,pcs_t) %>% 
      group_by(clone) %>% summarise(median=median(!!sym(pcs_t))) %>% 
      add_column(pval=0,logi_test=0,referencia = T_ref,PC=pcs_t) 
    median_PCref = dataref %>% select(clone,Condition,pcs_t) %>% 
      group_by(clone) %>% summarise(median=median(!!sym(pcs_t))) %>% 
      add_column(pval=0,logi_test=0,referencia = T_ref,PC=pcs_t) 
    #For each component
    for (j in 1:nrow(median_PC)) {
      nom=levels(factor(median_PC$clone[j]))
      temporal = datatest %>%  filter(clone==nom) %>% select(pcs_t)
      temporal_ref = dataref %>%  filter(clone==refname) %>% select(pcs_t)
      ref_mitjana = median_PCref %>% filter(clone == refname) %>% select(median)
      test_mitjana = median_PC %>% filter(clone == nom) %>% select(median)
      if(test_mitjana<ref_mitjana) {
        testtype="less"
      } else if(test_mitjana>ref_mitjana) {
        testtype="greater"
      } else {
        testtype="two.sided"
      }
      median_PC$logi_test[j]=testtype
      test=t.test(temporal, y = temporal_ref,
                  alternative = testtype,
                  mu = 0, paired = FALSE, var.equal=FALSE, exact = TRUE, correct = TRUE,
                  conf.int = FALSE, conf.level = 0.95)
      median_PC$pval[j]=test$p.value
      
    }
    comparison=rbind(comparison,median_PC)
    
  }
  
  return(comparison)
}



#Read the PCA transformmed data
data606=read.csv('606_scores.csv')
data607=read.csv('607_scores.csv')
allData=read.csv("both606And607_scores.csv")

## Format columns
data606$Condition=as.factor(data606$Condition)
data606$Treatment=as.factor(data606$Treatment)
data606$Name=as.factor(data606$Name)
data607$Condition=as.factor(data607$Condition)
data607$Treatment=as.factor(data607$Treatment)
data607$Name=as.factor(data607$Name)
data606$clone = data606$Name
data607$clone = data607$Name
## For both 606 and 607 combined:
allData$Condition=as.factor(allData$Condition)
allData$Treatment=as.factor(allData$Treatment)
allData$Name=as.factor(allData$Name)
allData$clone = allData$Name

#comparisons
refnameVec = c("606P","607P")
for (k in 1:length(refnameVec)) {
  if (refnameVec[k] == "606P"){
    data = data606
    dataToAddTo = data606
  } else {
    data = data607
    dataToAddTo = data607
  }
  T_refVec = c(37, 15, 37, 43)
  T_testVec = c(15, 15, 43, 43)
  # Loop through all the condition comparisons to get data to use to calculate restoration in next section:
  for (i in 1:length(T_refVec)) {
    T_ref = T_refVec[i]
    T_test = T_testVec[i]
    refname = refnameVec[k]
    comp_file=compare_to_ancestral(data, T_ref, T_test ,refname)
    file_name=paste(refname,T_test,"vs",T_ref,".csv",sep= "_")
    write.csv(comp_file,file_name)
  }
  
  ##########################################################3
  #To compare
  T_refVec = c(15, 43)
  T_testVec = c(15, 43)
  # loop through both 15 and 43 conditions:
  for (i in 1:length(T_refVec)) {
    T_ref = T_refVec[i]
    T_test = T_testVec[i]
    refname = refnameVec[k]
    #Load the data and the reference set
    data_name=paste(refname,T_test,"vs",T_ref,".csv",sep= "_")
    data=read.csv(data_name)
    alpha=0.05
    control_name=paste(refname,T_test,"vs","37",".csv",sep= "_")
    control=read.csv(control_name)
    #Rename Columns
    colnames(control)[4]="pvalref"
    colnames(data)[4]="pvaltest"
    data$logi_test=as.character(data$logi_test)
    control$logi_test=as.character(control$logi_test)
    #Multiple comparison correction
    control$pvalref=p.adjust(control$pvalref,method="fdr")
    data$pvaltest=p.adjust(data$pvaltest,method="fdr")
    data$logi_test[data$pvaltest>=alpha] = 'equal'
    control$logi_test[control$pvalref>=alpha] = 'equal'
    
    #Add PCdir to data
    references = control %>% filter((control$clone==refname))
    row.names(references)=references$PC
    data=data %>%  mutate(PCdir=references[as.character(PC),'logi_test'])
    
    #filter controls
    data = data %>% filter(!(data$clone==refname))
    control = control %>% filter(!(control$clone==refname))
    data$comparison=""
    
    # Test Classifyier
    for (i in 1:nrow(data)) {
      
      if (control$logi_test[i]=="less" && data$logi_test[i]=="greater" && data$PCdir[i]=="less") {
        data$comparison[i]="Partially Restored"} 
      if(control$logi_test[i]=="greater" && data$logi_test[i]=="less"  && data$PCdir[i]=="greater") {
        data$comparison[i]="Partially Restored"}
      if(control$logi_test[i]=="less" && data$logi_test[i]=="less" && data$PCdir[i]=="less") {
        data$comparison[i]="Reinforced"}
      if(control$logi_test[i]=="greater" && data$logi_test[i]=="greater" && data$PCdir[i]=="greater") {
        data$comparison[i]="Reinforced"}
      if(control$logi_test[i]=="less" && data$logi_test[i]=="less" && data$PCdir[i]=="greater") {
        data$comparison[i]="Over-restored"}
      if(control$logi_test[i]=="greater" && data$logi_test[i]=="greater" && data$PCdir[i]=="less") {
        data$comparison[i]="Over-restored"}
      if(control$logi_test[i]=="greater" && data$logi_test[i]=="equal" && data$PCdir[i]=="greater") {
        data$comparison[i]="Unrestored"}
      if(control$logi_test[i]=="less" && data$logi_test[i]=="equal"&& data$PCdir[i]=="less") {
        data$comparison[i]="Unrestored"}
      if(control$logi_test[i]=="equal" && data$logi_test[i]=="greater"  && data$PCdir[i]=="less") {
        data$comparison[i]="Restored"}
      if(control$logi_test[i]=="equal" && data$logi_test[i]=="less"  && data$PCdir[i]=="greater") {
        data$comparison[i]="Restored"}
      if(control$logi_test[i]=="less" && data$logi_test[i]=="less"  && data$PCdir[i]=="equal") {
        data$comparison[i]="Novel"}
      if(control$logi_test[i]=="greater" && data$logi_test[i]=="greater"  && data$PCdir[i]=="less") {
        data$comparison[i]="Novel"}
      if(control$logi_test[i]=="equal" && data$logi_test[i]=="less"  && data$PCdir[i]=="equal") {
        data$comparison[i]="Novel"}
      if(control$logi_test[i]=="greater" && data$logi_test[i]=="greater"  && data$PCdir[i]=="equal") {
        data$comparison[i]="Novel"}
      if(control$logi_test[i]=="less" && data$logi_test[i]=="less"  && data$PCdir[i]=="greater") {
        data$comparison[i]="Novel"}
      if(control$logi_test[i]=="less" && data$logi_test[i]=="less"  && data$PCdir[i]=="greater") {
        data$comparison[i]="Novel"}
      if(control$logi_test[i]=="less" && data$logi_test[i]=="equal"  && data$PCdir[i]=="equal") {
        data$comparison[i]="Uninformative"}
      if(control$logi_test[i]=="greater" && data$logi_test[i]=="equal"  && data$PCdir[i]=="equal") {
        data$comparison[i]="Uninformative"}
      if(control$logi_test[i]=="equal" && data$logi_test[i]=="equal"  && data$PCdir[i]=="equal") {
        data$comparison[i]="Uninformative"}
      if(control$logi_test[i]=="equal" && data$logi_test[i]=="equal"  && data$PCdir[i]=="less") {
        data$comparison[i]="Uninformative"}
      if(control$logi_test[i]=="equal" && data$logi_test[i]=="equal"  && data$PCdir[i]=="greater") {
        data$comparison[i]="Uninformative"}
    }
    
    # save
    data$Name = data$clone
    data$comparison = factor(data$comparison, levels = c(
      'Unrestored','Partially Restored','Reinforced','Uninformative','Novel','Restored'))
    restore =  inner_join(data,dataToAddTo, by = "Name")
    #Save
    file_name=paste(refname,T_test,"directionality.csv",sep= "_")
    write.csv(restore,file_name)
  }
}



######################################################
# Read the results:
files = list.files(pattern = "*directionality.csv")
direccio = read.csv(files[1])

for (i in 2:length(files)){
  directe = read.csv(files[i]) 
  direccio = rbind(direccio,directe)
}
head(direccio)

#pdf("BiologRestorationByConditionAndStrain_Fig.pdf")
ggplot(direccio) + aes(x = PC, fill =comparison)+
  geom_bar(position = "fill")+ ylab("Relative Abundance") + 
  xlab("Environment")+facet_grid(Treatment~referencia~Strain~Strategy)
#dev.off()