 #!/usr/bin/Rscript

#####################################################################
####      Parametric Curves Fitting Comparison.                    ##
#####################################################################

#Input

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Libraries
library(ggplot2)
library(gtable)
library(systemfit)

for (i in 1:length(args)){

	name=args[i];

	#Reading datasets and filtering

	data=read.csv(name,header=TRUE);

	#Fitting the data to different simple growth models

	fit1=nls(OD~a*exp(-exp((g*exp(0)/a)*(l-temps)+1)),data, start = list(a = 0.3,g=0.002,l=300)) #Gompertz model

	fit2=nls(OD~0.0001+(a-0.0001)/(1+exp(g*(l-temps))), data, start = list(a = 0.3,g=0.002,l=300))  #Isotherm Zweitering

	fit3=nls(OD~b+(a-b)/((c+q*exp(g*(l-temps))^(1/nu))),data, start=list(a=0.5,g=0.002,c=1,nu=0.2,b=0.0000001,q=0.0001,l=300)) #Generalised logistic model

	fit4=nls(OD~0.0001+(a-0.0001)*exp(-exp(g*exp(0)/(a-0.0001)*(lambda-temps)+1)),data,start=list(a=0.3,g=0.002,lambda=300))



	#Fitting the data to different compound growth models

	fit4=nls(OD~b+a-log(exp(b)+(exp(a)-exp(b))*exp(-g*(temps+1/4*log((1+exp(4*(lambda-temps)))/(1+exp(4*lambda)))))),data,start=list(b=0.001,a=0.3,g=0.002,lambda=300)) #Huang Model

	fit5=nls(OD~0.0001+g*(temps+1/g*log(exp(-g*temps)+exp(-h)-exp(-g*temps-h)))-log(1+(exp(g*(temps+1/g*log(exp(-g*temps)+exp(-h)-exp(-g*temps-h))))-1)/(exp(a-0.0001))), data,start=list(a=0.3,g=0.002,h=30))

	A=temps+1/g*log(exp(-g*temps)+exp(-h)-exp(-g*temps-h))


	#plots

	prefix="fitplot.pdf";
	name_file=paste(name,prefix,sep="_");
	Temp=as.character(Temp)


	pdf(file=name_file)
		print(
			ggplot(data,aes(x = temps,y = OD,color=Instrument)) + geom_point() + 
			geom_smooth(method = 'nls', formula = 'y~a*exp(-exp((g*exp(0)/a)*(l-x)+1))', start = list(a = 0.3,g=0.002,l=300),se=FALSE,color='blue')+
			geom_smooth(method = 'nls', formula = 'y~0.0001+(a-0.0001)/(1+exp(g*(l-x)))', start = list(a = 0.3,g=0.002,l=300),se=FALSE,color='green')+
			geom_smooth(method = 'nls', formula = 'y~0.0001+(a-0.0001)*exp(-exp(g*exp(0)/(a-0.0001)*(lambda-x)+1))', start = list(a = 0.3,g=0.002,lambda=300),se=FALSE,color='red')+
			xlab("Time (minutes)") + ylab("OD 850")+ facet_wrap(~Instrument)
		)
	dev.off()
	
}
