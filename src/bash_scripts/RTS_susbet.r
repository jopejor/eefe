 #!/usr/bin/Rscript

#####################################################################
####      McFarland Calibration of the RTS instruments             ##
#####################################################################

#Input

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)


for (i in 1:length(args)){

 name=args[i];


 #Reading datasets

 data=read.csv(name,skip=1,header=FALSE);

 	OD=as.numeric(as.character(data$V7));

 	Lambda=factor(data$V11);

 	Instrument=factor(data$V2);

 	Time=data$V4

 	Temperature=as.numeric(as.character(data$V6))

 data_raw=data.frame(OD,Lambda,Instrument,Time,Temperature);

 #Writing 

 prefix="filtered";

 name_file=paste(prefix,name,sep="_");

 write.csv(data_raw,name_file, row.names=FALSE, col.names=TRUE);

}
