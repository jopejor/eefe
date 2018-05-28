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
############################### 18 ##############################


data_2=data_baseline[which(is.na(data_baseline$text)), ]

 for (i in 1:length(levels(data_1$Instrument))) {

	to_Bottom=subset(data_1,data_1$Instrument==levels(data_1$Instrument)[i])



	spl <- smooth.spline(to_Bottom$temps, to_Bottom$OD, spar=0.7)
	deriv=predict(spl, deriv=1)

	prefix="derivative.pdf"

	name_file=paste(levels(data_1$Instrument)[i],prefix,sep="_");


	pdf(file=name_file)
	plot(to_Bottom$temps,to_Bottom$OD,xlab="Time(minutes)", ylab="OD 850 nm")
	lines(predict(spl),col=2)
	par(new = T)
	plot(deriv$x,deriv$y,pch=16, axes=F, xlab=NA, ylab= "Derivative (OD units/minute)" , cex=1.2)
	axis(side = 4)
	dev.off()




	instrument=rep(levels(data_1$Instrument)[i],length(deriv$y));

	deriva_df=data.frame(deriv$y,spl$y,spl$x,instrument)

	data_2=rbind(data_2,deriva_df)

 	}

name_file="derivative_lyanna.csv"

write.csv(data_2,name_file, row.names=FALSE, col.names=TRUE);




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

