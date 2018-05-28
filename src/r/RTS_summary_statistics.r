 #!/usr/bin/Rscript

#####################################################################
####      McFarland Calibration of the RTS instruments             ##
#####################################################################

#Input

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Libraries
library(ggplot2)

# Reading files


for (i in 1:length(args)){

 name=args[i];


 #Reading datasets and filtering

 data=read.csv(name,header=TRUE);


#Plots

 prefix="growth_rate.pdf";

 name_file=paste(name,prefix,sep="_");

pdf(file=name_file)
 print(
    ggplot(data, aes(factor(Strain), rate))+ geom_boxplot()+geom_jitter()+
    xlab("Fluctuation")+ylab("Growth rate")
      )
dev.off()

#Plots

 prefix="carrying_rate.pdf";

 name_file=paste(name,prefix,sep="_");

pdf(file=name_file)
 print(
    ggplot(data, aes(factor(Strain), carrying))+ geom_boxplot()+geom_jitter()+
    xlab("Fluctuation")+ylab("Carrying capacity")
      )
dev.off()

#Plots

 prefix="lag.pdf";

 name_file=paste(name,prefix,sep="_");

pdf(file=name_file)
 print(
    ggplot(data, aes(factor(Strain), lag))+ geom_boxplot()+geom_jitter()+
    xlab("Fluctuation")+ylab("Lag")
      )
dev.off()


}
