 #!/usr/bin/Rscript

####################################################################################
####      Data processing for the experimental evolution by Ximo Pechuan          ##
####################################################################################

# Libraries
library(ggplot2)
library(gtable)
library(grid)
library("RColorBrewer")
library(reshape)
library(pspline)
library(grofit)

#Input and data processing 
name="filtered_.csv"
data=read.csv(name,header=TRUE);
data_baseline=data;
data_baseline=subset(data_baseline,data_baseline$OD>-0.2); #To remove negative values
data_baseline=subset(data_baseline,data_baseline$OD<1); #To remove negative values
#### Measure the Baseline #####
bases=1:length(levels(data_baseline$Instrument));
Strain=1:length(levels(data_baseline$Instrument));
for (i in 1:length(levels(data_baseline$Instrument))){
	liner=subset(data_baseline,Instrument==levels(data_baseline$Instrument)[i]);
	bases[i]=mean(liner$OD[1:10])
	Strain[i]=tolower(as.character(levels(data_baseline$Instrument)[i]);
} 
bases_lines=data.frame(Strain,bases);
#Correct the baseline 
prefix="basline.csv"
name_file=paste(prefix,name,sep="_");
write.csv(bases_lines,name_file, row.names=FALSE, col.names=TRUE);

data_baseline=data_baseline[which(is.na(data_baseline$text)), ]

for (i in 1:length(levels(data$Instrument))) {

to_Bottom=subset(data,data$Instrument==levels(data$Instrument)[i])
Bottom=min(to_Bottom$OD);

for (j in 1:nrow(to_Bottom)){

	if(Bottom>0){
		to_Bottom$OD[j]=to_Bottom$OD[j]-Bottom;
	} 

	if(Bottom<=0){
		to_Bottom$OD[j]=to_Bottom$OD[j]+abs(Bottom);
	}
}

data_baseline=rbind(data_baseline,to_Bottom)
}

#Inititializing time variables

hores=1:length(data_baseline$Time);
minuts=1:length(data_baseline$Time);
dates=1:length(data_baseline$Time);
#Time handling

for (i in 1:length(data_baseline$Time)){

temp=unlist(strsplit(as.character(data_baseline$Time[i]),' '));

hores[i]=unlist(strsplit(as.character(temp),":"))[5];

minuts[i]=unlist(strsplit(as.character(data_baseline$Time[i]),' '))[6];

dates[i]=unlist(strsplit(as.character(data_baseline$Time[i]),' '))[1];
}

dies=1:length(dates);

days=sort(unique(as.factor(dates)));

for (i in 1:length(dates)) {

for(j in 1:length(days)){

if (dates[i]==days[j]) {

	dies[i]=j;} 
}
}

hores=as.numeric(hores);
minuts=as.numeric(minuts)+hores*60;


temps=dies*24*60+minuts;

temps=temps-min(temps);

data_baseline$temps=temps;

#Data filtering
data_1=data_baseline

#Plot to detect anomalies
ggplot(data_1,aes(x = temps,y = OD,color=Instrument)) + geom_point() + 
#geom_smooth(method = 'nls', formula = 'y~a*exp(-exp((g*exp(0)/a)*(l-x)+1))', start =list(a = 0.3,g=0.0002,l=100),se=FALSE)+
xlab("Time (minutes)") + ylab("OD 850")+facet_wrap(~Instrument)

#Manual fitering options
data_1=subset(data_1,data_1$OD<-0.25); # By OD
data_1$temps=data_1$temps-min(data_1$temps); #Zero out the time
data_1=subset(data_1,data_1$temps>500); #By time
#Keep only the factors
data_1=subset(data_1,data_1$Instrument=="west"|data_1$Instrument=="zoot"|data_1$Instrument=="sunfyre")
#To remove instruments, repeat for each level of the factor you want to remove 
data_1=subset(data_1,data_1$Instrument!="")
prefix="fitplot.pdf";
name_file=paste(name,prefix,sep="_");
Temp=as.character(Temp)
#Plot of the processed data
pdf(file=name_file)
 print(
  ggplot(data_1,aes(x = temps,y = OD,color=Instrument)) + geom_point() + 
  #geom_smooth(method = 'nls', formula = 'y~a*exp(-exp((g*exp(0)/a)*(l-x)+1))', start =list(a = 0.3,g=0.0002,l=100),se=FALSE)+
  xlab("Time (minutes)") + ylab("OD 850")+facet_wrap(~Instrument)
  )
dev.off()

#Keeps the processed data
data_1$Instrument=factor(data_1$Instrument) #Drops the removed factors
prefix="data.csv"
name_file=paste(prefix,name,sep="_");
write.csv(data_1,name_file, row.names=FALSE, col.names=TRUE);

#Measurements of the growth parameters
#Functions
neg<-function(x) -x 
#Initialize variables
maximum_growth=1:length(levels(data_1$Instrument));
soca=1:length(levels(data_1$Instrument));
temperatura=1:length(levels(data_1$Instrument));
groFit_mu=1:length(levels(data_1$Instrument))
lags=1:length(levels(data_1$Instrument));
carrying=1:length(levels(data_1$Instrument));
#Measurements per factor
for (i in 1:length(levels(data_1$Instrument))) {
	to_Bottom=subset(data_1,data_1$Instrument==levels(data_1$Instrument)[i])
	spl <- smooth.spline(to_Bottom$temps, to_Bottom$OD,spar=0.85)
	deriv=predict(spl, deriv=1)
	#Correct splittting 0.85 for 15 0.65 for 43
	lol=predict(spl,to_Bottom$temps)
	difference=to_Bottom$OD-lol$y
	to_Bottom$difference=difference
	to_Bottom=subset(to_Bottom, to_Bottom$difference>neg(0.02))
	spl <- smooth.spline(to_Bottom$temps, to_Bottom$OD,spar=0.85)
	prefix="phenotype.pdf"
	name_file=paste(levels(data_1$Instrument)[i],prefix,sep="_");
	#Plot
	pdf(file=name_file)
		par(mar=c(5, 4, 4, 8) + 0.1)
 		plot(to_Bottom$temps,to_Bottom$OD,xlab="Time(minutes)", ylab="OD 850 nm")
 		lines(predict(spl),col=2)
 		par(new = T)
 		plot(deriv$x,deriv$y,pch=16, axes=F, xlab=NA, ylab= "" , cex=1.2)
 		axis(side = 4)
 		mtext("Derivative (OD units/minute)", side=4, line=3, cex.lab=1)
 	dev.off()
	#groFit_mu
	test=gcFitSpline(to_Bottom$temps,to_Bottom$OD)

	prefix="groFit_mu_phenotype.pdf"
	name_file=paste(levels(data_1$Instrument)[i],prefix,sep="_");
	#Groofit plot
	pdf(file=name_file)
			plot(test)
	dev.off()

	maximum_growth[i]=max(deriv$y)
	groFit_mu[i]=as.numeric(test$parameters[2]);
	soca[i]=name_file
	temperatura[i]=mean(to_Bottom$Temperature)
	lags[i]=as.numeric(test$parameters[3]);
	carrying[i]=as.numeric(test$parameters[1]);

	}

#Save to a dataframe
parex=data.frame(soca,maximum_growth,temperatura,groFit_mu,lags,carrying)
prefix="parametres.csv";
name_file=paste(name,prefix,sep="_");
write.csv(parex,name_file, row.names=FALSE, col.names=TRUE);