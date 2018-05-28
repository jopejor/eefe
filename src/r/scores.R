


#### PCA loadings #######
mydata=data[5:92]
fit=prcomp(mydata, center = FALSE, scale = FALSE, retx = TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

###################### Compute the paiwise Wilcox test for each PC and Temperature ###############################

data=read.csv(name)


## Format columns
data$Condition=as.factor(data$Condition)
data$Treatment=as.factor(data$Treatment)
data$Strain=as.factor(data$Strain)
#Set the two datasets to be compared
datatest=subset(data,data$Condition=="43")
dataref=subset(data,data$Condition=="43")

### PC1
median43_PC1=aggregate(x=datatest$PC1, by=list(datatest$Strain), FUN=median)
median43_PC1$pval=1:length(median43_PC1$Group.1);
median43_PC1$logi_test=1:length(median43_PC1$Group.1);
median43_PC1$referencia=rep("43",length(median43_PC1$Group.1))
referencia=subset(dataref$PC1,dataref$Strain==dataref$Strain[1])

for (i in 1:length(median43_PC1$Group.1)) {

	
	nom=levels(factor(median43_PC1$Group.1[i]))
	temporal=subset(datatest$PC1,datatest$Strain==nom)

	if(median(temporal)<median(referencia)) {
		testtype="less"
	} else if(median(temporal)>median(referencia)) {
		testtype="greater"
	} else {
		testtype="two.sided"
            }
 	
    median43_PC1$logi_test[i]=testtype

	test=t.test(temporal, y = referencia,
            alternative = testtype,
            mu = 0, paired = FALSE, var.equal=FALSE, exact = TRUE, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95)


	median43_PC1$pval[i]=test$p.value
}

median43_PC1$PC=rep("PC1",length(median43_PC1$Group.1))

########################################### PC2####################################################

median43_PC2=aggregate(x=datatest$PC2, by=list(datatest$Strain), FUN=median)


median43_PC2$pval=1:length(median43_PC2$Group.1);
median43_PC2$logi_test=1:length(median43_PC2$Group.1);
median43_PC2$referencia=median43_PC1$referencia

referencia=subset(dataref$PC2,dataref$Strain==dataref$Strain[1])
for (i in 1:length(median43_PC2$Group.1)) {

	nom=levels(factor(median43_PC2$Group.1[i]))
	temporal=subset(datatest$PC2,datatest$Strain==nom)

	if(median(temporal)<median(referencia)) {
		testtype="less"
	} else if(median(temporal)>median(referencia)) {
		testtype="greater"
	} else {
		testtype="two.sided"
            }
 	
    median43_PC2$logi_test[i]=testtype

	test=t.test(temporal, y = referencia,
            alternative = testtype,
            mu = 0, paired = FALSE, var.equal=FALSE, exact = TRUE, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95)


	median43_PC2$pval[i]=test$p.value
}

median43_PC2$PC=rep("PC2",length(median43_PC1$Group.1))


########################################### PC3####################################################

median43_PC3=aggregate(x=datatest$PC3, by=list(datatest$Strain), FUN=median)


median43_PC3$pval=1:length(median43_PC3$Group.1);
median43_PC3$logi_test=1:length(median43_PC3$Group.1);
median43_PC3$referencia=median43_PC1$referencia

referencia=subset(dataref$PC3,dataref$Strain==dataref$Strain[1])
for (i in 1:length(median43_PC3$Group.1)) {

	nom=levels(factor(median43_PC3$Group.1[i]))
	temporal=subset(datatest$PC3,datatest$Strain==nom)

	if(median(temporal)<median(referencia)) {
		testtype="less"
	} else if(median(temporal)>median(referencia)) {
		testtype="greater"
	} else {
		testtype="two.sided"
            }
 	
    median43_PC3$logi_test[i]=testtype

	test=t.test(temporal, y = referencia,
            alternative = testtype,
            mu = 0, paired = FALSE, var.equal=FALSE, exact = TRUE, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95)


	median43_PC3$pval[i]=test$p.value
}

median43_PC3$PC=rep("PC3",length(median43_PC1$Group.1))

########################################### PC4####################################################

median43_PC4=aggregate(x=datatest$PC4, by=list(datatest$Strain), FUN=median)


median43_PC4$pval=1:length(median43_PC4$Group.1);
median43_PC4$logi_test=1:length(median43_PC4$Group.1);
median43_PC4$referencia=median43_PC1$referencia

referencia=subset(dataref$PC4,dataref$Strain==dataref$Strain[1])
for (i in 1:length(median43_PC4$Group.1)) {

	nom=levels(factor(median43_PC4$Group.1[i]))
	temporal=subset(datatest$PC4,datatest$Strain==nom)

	if(median(temporal)<median(referencia)) {
		testtype="less"
	} else if(median(temporal)>median(referencia)) {
		testtype="greater"
	} else {
		testtype="two.sided"
            }
 	
    median43_PC4$logi_test[i]=testtype

	test=t.test(temporal, y = referencia,
            alternative = testtype,
            mu = 0, paired = FALSE, var.equal=FALSE, exact = TRUE, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95)


	median43_PC4$pval[i]=test$p.value
}

median43_PC4$PC=rep("PC4",length(median43_PC1$Group.1))

################### Final ####

median43=rbind(median43_PC1,median43_PC2,median43_PC3,median43_PC4)

write.csv(median43,"607_test43vs43.csv")

################################################################ correction and classification ###############################

data=read.csv(name)
alpha=0.05
control=read.csv("607_control43vs37.csv")

for (i in 1:nrow(data)) {
	for (j in 1:nrow(control)) {
	if(data$PC[i]==control$PC[j]) {data$PCdir[i]=as.character(control$logi_test[j])}
}

}

# Multiple comparison correction

data$pvalref=p.adjust(data$pvalref,method="fdr")
data$pvaltest=p.adjust(data$pvalref,method="fdr")

for (i in 1:nrow(data)) {

	if (data$pvalref[i]<=alpha) {
		data$sigref[i]=1
	} else { data$sigref[i]=0}

	if (data$pvaltest[i]<=alpha) {
		data$sigtest[i]=1
	} else { data$sigtest[i]=0}

	}
# Logical reassessment

for (i in 1:nrow(data)) {

	if (data$sigref[i]==0) {
		data$siglogiref[i]="equal"
	} else {data$siglogiref[i]=as.character(data$logi_ref[i])}

	if (data$sigtest[i]==0) {
		data$siglogitest[i]="equal"
	} else {data$siglogitest[i]=as.character(data$logi_test[i])}

	}

# Test Classifyier

for (i in 1:nrow(data)) {

	if (data$siglogiref[i]=="less" && data$siglogitest[i]=="greater" && data$PCdir[i]=="less") {
		data$comparison[i]="Partially Restored"} 
	if(data$siglogiref[i]=="greater" && data$siglogitest[i]=="less"  && data$PCdir[i]=="greater") {
		data$comparison[i]="Partially Restored"}
	if(data$siglogiref[i]=="less" && data$siglogitest[i]=="less" && data$PCdir[i]=="less") {
		data$comparison[i]="Reinforced"}
	if(data$siglogiref[i]=="greater" && data$siglogitest[i]=="greater" && data$PCdir[i]=="greater") {
		data$comparison[i]="Reinforced"}
	if(data$siglogiref[i]=="less" && data$siglogitest[i]=="less" && data$PCdir[i]=="greater") {
		data$comparison[i]="Over-restored"}
	if(data$siglogiref[i]=="greater" && data$siglogitest[i]=="greater" && data$PCdir[i]=="less") {
		data$comparison[i]="Over-restored"}
	if(data$siglogiref[i]=="greater" && data$siglogitest[i]=="equal" && data$PCdir[i]=="greater") {
		data$comparison[i]="Unrestored"}
	if(data$siglogiref[i]=="less" && data$siglogitest[i]=="equal"  && data$PCdir[i]=="less") {
		data$comparison[i]="Unrestored"}
	if(data$siglogiref[i]=="equal" && data$siglogitest[i]=="less"  && data$PCdir[i]=="less") {
		data$comparison[i]="Restored"}
	if(data$siglogiref[i]=="equal" && data$siglogitest[i]=="greater"  && data$PCdir[i]=="greater") {
		data$comparison[i]="Restored"}
	if(data$siglogiref[i]=="less" && data$siglogitest[i]=="less"  && data$PCdir[i]=="equal") {
		data$comparison[i]="Novel"}
	if(data$siglogiref[i]=="greater" && data$siglogitest[i]=="greater"  && data$PCdir[i]=="equal") {
		data$comparison[i]="Novel"}
	if(data$siglogiref[i]=="equal" && data$siglogitest[i]=="equal"  && data$PCdir[i]=="equal") {
		data$comparison[i]="Uninformative"}

	}


write.csv(data,"607_43_comparisons.csv")

compares.freq = table(data$comparison) 
pie(compares.freq)  

################### 

mdata <- melt(mydata, id=c("id","time"))

for (i in 1:nrow(mut)) {
	gene=as.character(mut$Mutation[j])
	strain=colnames(mut[])
	for (j in 1:nrow(case)) {
		if(case$)
		




} 
}
