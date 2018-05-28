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

 prefix="growth.pdf";

 name_file=paste(name,prefix,sep="_");

pdf(file=name_file)
 print(
  ggplot(data,aes(x = Temp,y = rate, colour=factor(Strain))) + geom_point(shape=1,) +
   stat_smooth(method = 'nls', formula = 'y~(b*(x-L)*(1-exp(c*(x-H))))^2', start = list(b = 1,c=1,L=10,H=50),se=FALSE)+
  ylab("Growth rate") + xlab("Temperature(Celsius)")
  )
dev.off()

 prefix="carrying.pdf";

 name_file=paste(name,prefix,sep="_");

pdf(file=name_file)
 print(
  ggplot(data,aes(x = Temp,y = carrying,colour=factor(Strain))) + geom_point(shape=1) + stat_smooth(method = 'nls', formula = 'y~b*(1-exp(c*(x-H)))', start = list(b = 8,c=1,H=45),se=FALSE)+
  ylab("Carrying capacity") + xlab("Temperature(Celsius)")
  )
dev.off()

prefix="lag.pdf";

 name_file=paste(name,prefix,sep="_");

pdf(file=name_file)
 print(
  ggplot(data,aes(x = Temp,y = lag,colour=factor(Strain))) + geom_point(shape=1) + stat_smooth(method = 'nls', formula = 'y~exp(b/(x-c))', start = list(b = 23,c=2),se=FALSE)+
  ylab("lag") + xlab("Temperature(Celsius)")
  )
dev.off()





}