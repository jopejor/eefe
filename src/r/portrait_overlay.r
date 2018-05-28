 #!/usr/bin/Rscript

#####################################################################
####      Phase portrait comparison						           ##
#####################################################################

#Input

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Libraries
library(ggplot2)
library(gtable)
library(grid)
library("RColorBrewer")
library(reshape)

# Reading files

portrait=read.csv("derivative_37.csv")

for (i in 1:length(args)){

 	name=args[i];

 	data=read.csv(name)

 	prefix="derivative_37_portrait.pdf"

 	name_file=paste(name,prefix,sep="_");


   	pdf(file=name_file)
   		plot(portrait$spl.y,portrait$deriv.y)
   		points(data$spl.y,data$deriv.y,col=2)
	 dev.off()



}