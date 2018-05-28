#####################################################################
####      McFarland Calibration of the RTS instruments             ##
#####################################################################

# Libraries
install.packages("ggplot2")
library(ggplot2)


# Functions

lm_eqn = function(df){
    m = lm(df$x ~ df$Group.1, df);
    eq <- substitute(y == a + b*x*","~~r^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 3), 
              b = format(coef(m)[2], digits = 3), 
             r2 = format(summary(m)$r.squared, digits = 4)))
    as.character(as.expression(eq));                 
}

# Reading files

data=read.csv("combined");

spec=read.csv("spectro.csv");

# Formatting spectrophotometer reads

Group.2=rep("Spec",length(spec$Volume));

spec=data.frame(spec,Group.2);

dev.x=rep(0.001,length(spec$Volume));

spec=data.frame(spec,dev.x);

names(spec)[2]="Group.1";

names(spec)[1]="x";

# Filtter 17

data=subset(data, Volume!=17);


#  Boxplot

pdf(file="boxplot.pdf")
p <- ggplot(data, aes(factor(Volume), OD))
p + geom_boxplot(aes(fill = factor(Instrument)))+ 
    xlab("Volume (mL)") + ylab("OD 850")+ggtitle("Behaviour of the OD with the volume variation")+
    scale_fill_discrete(name="Instruments", breaks=c("COM3", "COM4", "COM5"), labels=c("RTS 3 ", "RTS 2", "RTS 4"))
dev.off()

#  Scatterplot

means=aggregate(data$OD, by=list(data$Volume, data$Instrument), FUN=mean);

# Standard Error

std <- function(x) sd(x)/sqrt(length(x))

dev=aggregate(data$OD, by=list(data$Volume, data$Instrument), FUN=std);

means=data.frame(means,dev$x);

means=rbind(means, spec);

# Linear regression

factor_levels=levels(means$Group.2);

for (i in 1:length(factor_levels)){

 sub=subset(means, Group.2==factor_levels[i]);

 assign(paste('sub', factor_levels[i], sep=''), sub);

}


# Error bars

limits <- aes(ymax = means$x + means$dev.x, ymin=means$x - means$dev.x)

pdf(file="scatterplot.pdf")
p=ggplot(means, aes(x=means$Group.1, y=means$x, color=means$Group.2)) + geom_point(shape=1)+scale_colour_hue(l=50) 
p + geom_smooth(method=lm,se=FALSE) + geom_errorbar(limits, width=0.2) + 
    xlab("Volume (mL)") + ylab("OD 850")+ggtitle("Behaviour of the OD with the volume variation")+ 
    scale_color_discrete(name="Instruments", breaks=c("COM3", "COM4", "COM5", "Spec"), labels=c("RTS 3 ", "RTS 2", "RTS 4", "spectrophotometer"))+
    annotate("text", x = 14, y = 2.2, label = lm_eqn(subCOM3), parse=TRUE)+
    annotate("text", x = 14, y = 2, label = lm_eqn(subCOM4), parse=TRUE)+
    annotate("text", x = 14, y = 1.8, label = lm_eqn(subCOM5), parse=TRUE)+
    annotate("text", x = 20, y = 0.5, label = lm_eqn(subSpec), parse=TRUE)
dev.off()
