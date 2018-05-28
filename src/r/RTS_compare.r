 #!/usr/bin/Rscript

#####################################################################
####      Compare the acclimation and acclimated strains           ##
#####################################################################

#Input

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Libraries
library(ggplot2)

# Reading files

group_1=args[1]
group_2=args[2]

data_1=read.csv(group_1,header=TRUE);
data_2=read.csv(group_2,header=TRUE);

# Creating factors 

for (i in 1:length(data_1$Strain)){

	data_1$Temperature[i]=unlist(strsplit(as.character(data_1$Strain[i]),'-'))[1];
	data_1$REL[i]=unlist(strsplit(as.character(data_1$Strain[i]),'-'))[2];
	data_1$Repetion[i]=unlist(strsplit(as.character(data_1$Strain[i]),'-'))[3];

}


for (i in 1:length(data_2$Strain)){

	data_2$Temperature[i]=unlist(strsplit(as.character(data_2$Strain[i]),'-'))[1];
	data_2$REL[i]=unlist(strsplit(as.character(data_2$Strain[i]),'-'))[2];
	data_2$Repetion[i]=unlist(strsplit(as.character(data_2$Strain[i]),'-'))[3];

}

data_1$Pass="First"
data_2$Pass="Second"

# Merge

total=rbind(data_1,data_2)
total$Temperature=as.factor(total$Temperature)

for (i in 1:length(levels(total$Temperature))) {

	temp=subset(total, total$Temperature==levels(total$Temperature)[i])
	box=boxplot(rate~REL*Pass,data=temp)

	prefix="boxplot.pdf";
    name_file=paste(temp$Temperature[1],prefix,sep="_");

    pdf(file=name_file)
 		print(
			boxplot(rate~REL*Pass,data=temp)
  				)
	dev.off()

#General linear model, two factors with interaction

results=lm(rate~REL*Pass,data=temp);

anova_results=anova(results);

write.csv(anova_results,file="two_way_anova.csv");

#General linear model, two factors no interaction


results=lm(rate~REL+Pass,data=temp);

anova_results=anova(results);

write.csv(anova_results,file="two_way_anova_additive.csv");


}
















