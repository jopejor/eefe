 #!/usr/bin/Rscript

#####################################################################
####      McFarland Calibration of the RTS instruments             ##
#####################################################################

# Libraries
library(ggplot2)
library(gtable)
library(grid)
library("RColorBrewer")
library(reshape)

neg<-function(x) -x 

############################# S-Phenotype ##################################

	maximum_growth=1:length(levels(data_1$Instrument));
    soca=1:length(levels(data_1$Instrument));
    temperatura=1:length(levels(data_1$Instrument));
    slopeses=1:length(levels(data_1$Instrument));
    lags=1:length(levels(data_1$Instrument));
    carrying=1:length(levels(data_1$Instrument));

 for (i in 1:length(levels(data_1$Instrument))) {

 	to_Bottom=subset(data_1,data_1$Instrument==levels(data_1$Instrument)[i])


	spl <- smooth.spline(to_Bottom$temps, to_Bottom$OD)
		deriv=predict(spl, deriv=1)


	temp2=subset(to_Bottom,to_Bottom$temps>550)
	temp2=subset(temp2,temp2$temps<600)
	slopes=lm(temp2$OD ~ temp2$temps)

 	lol=predict(spl,to_Bottom$temps)
	difference=to_Bottom$OD-lol$y
	to_Bottom$difference=difference

	to_Bottom=subset(to_Bottom, to_Bottom$difference>neg(0.017))

	spl <- smooth.spline(to_Bottom$temps, to_Bottom$OD)
 	deriv=predict(spl, deriv=1)

 	 prefix="pilot.pdf"

 	name_file=paste(levels(data_1$Instrument)[i],prefix,sep="_");


   	pdf(file=name_file)

   		par(mar=c(5, 4, 4, 8) + 0.1)
	 	plot(to_Bottom$temps,to_Bottom$OD,xlab="Time(minutes)", ylab="OD 850 nm")
	 	lines(predict(spl),col=2)
	 	abline(lm(temp2$OD ~ temp2$temps))
	 	par(new = T)
	 	plot(deriv$x,deriv$y,pch=16, axes=F, xlab=NA, ylab= "" , cex=1.2)
	 	axis(side = 4)
	 	mtext("Derivative (OD units/minute)", side=4, line=3, cex.lab=1)

	 dev.off()



	maximum_growth[i]=max(deriv$y)
	soca[i]=prefix
	temperatura[i]=mean(to_Bottom$Temperature)
	slopeses[i]=slopes$coefficients[2]
	lags[i]=-slopes$coefficients[1]/slopes$coefficients[2]
	carrying[i]=max(to_Bottom$OD)


 	}
 
 name_file="derivative_lyanna.csv"

 write.csv(data_2,name_file, row.names=FALSE, col.names=TRUE);

 parex=data.frame(soca,maximum_growth,temperatura,slopeses,lags,carrying)

 name_file="parametres_experimentals_pilot.csv";

write.csv(parex,name_file, row.names=FALSE, col.names=TRUE);


####################################################



 	to_Bottom=subset(data_1,data_1$Instrument==levels(data_1$Instrument)[i])


	to_Bottom=subset(to_Bottom,to_Bottom$temps<2200)

	 	spl <- smooth.spline(to_Bottom$temps, to_Bottom$OD)

 	deriv=predict(spl, deriv=1)


	
	 	plot(to_Bottom$temps,to_Bottom$OD,xlab="Time(minutes)", ylab="OD 850 nm")
	 	lines(predict(spl),col=2)
	 	abline(lm(temp2$OD ~ temp2$temps))
	 	par(new = T)
	 	plot(deriv$x,deriv$y,pch=16, axes=F, xlab=NA, ylab= "Derivative (OD units/minute)" , cex=1.2)
	 	points(to_Bottom$temps,predict(new.spl,to_Bottom$temps,nderiv=1),col=6)
	 	axis(side = 4)
	  

	cold=subset(to_Bottom,to_Bottom$Temperature<19)
	hot=subset(to_Bottom,to_Bottom$Temperature>36)

	plot(cold$temps, cold$OD, xlim=c(0,2500),ylim=c(0,0.7))
	points(hot$temps,hot$OD)

	max(cold$OD)
	min(hot$OD)

	correct=max(cold$OD)-min(hot$OD)

	hot$OD=hot$OD+correct

 	spl <- smooth.spline(cold$temps, cold$OD, spar=0.72)

 	deriv=predict(spl, deriv=1)

	 	plot(cold$temps,cold$OD,xlab="Time(minutes)", ylab="OD 850 nm")
	 	lines(predict(spl),col=2)
	 	par(new = T)
	 	plot(deriv$x,deriv$y,pch=16, axes=F, xlab=NA, ylab= "Derivative (OD units/minute)" , cex=1.2)
	 	points(cold$temps,predict(new.spl,cold$temps,nderiv=1),col=6)
	 	axis(side = 4)


 prefix="derivative_cold.csv"

 name_file=paste(levels(data_1$Instrument)[i],prefix,sep="_");

 deriva_df=data.frame(deriv$y,spl$y,spl$x)

 write.csv(deriva_df,name_file, row.names=FALSE, col.names=TRUE);

 spl <- smooth.spline(hot$temps, hot$OD, spar=0.8)

 	deriv=predict(spl, deriv=1)

	 	plot(hot$temps,hot$OD,xlab="Time(minutes)", ylab="OD 850 nm")
	 	lines(predict(spl),col=2)
	 	par(new = T)
	 	plot(deriv$x,deriv$y,pch=16, axes=F, xlab=NA, ylab= "Derivative (OD units/minute)" , cex=1.2)
	 	points(hot$temps,predict(new.spl,hot$temps,nderiv=1),col=6)
	 	axis(side = 4)


 prefix="derivative_hot.csv"

 name_file=paste(levels(data_1$Instrument)[i],prefix,sep="_");

 deriva_df=data.frame(deriv$y,spl$y,spl$x)

 write.csv(deriva_df,name_file, row.names=FALSE, col.names=TRUE);

i=i+1



 	to_Bottom=subset(data_1,data_1$Instrument==levels(data_1$Instrument)[i])
