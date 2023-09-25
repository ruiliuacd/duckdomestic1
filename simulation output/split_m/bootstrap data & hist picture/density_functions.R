#!/usr/bin/R

# Functions for R and ggplot2 to create log-scale density plots

# Main function: density.log(...)
# Produce a dataframe used to produce a density plot
# Input:
#  data: dataframe, 
#  var: variable to plot on y-axis
#  split: factor to split on
#  n=512: See ?density
#  adjust=1: See ?density
#  title="": x-axis title for plot
#
# Output:
#  A list with three named values
#    points: a dataframe of the points to be on the density plot
#    summary: mean, standard deviation, max, min, median, 25th-percentile, 75th-percentile of var broken out by split
#    plot: a ggplot2 object of the log density plot
# 
# Example code
if (FALSE) {#This is just to make it easier to copy and paste the example code. Do *not* change to TRUE
 source("density_functions.R") #Reads in functions in this file. Make sure this file exists with the same filename

 ex.data<-data.frame(factor=c("A","B"),x=sample(1:10^5,1000,replace=T),y=sample(1:1000,1000,replace=T))
 density.output<-density.log(ex.data,var="y",split="factor",title="Example plot")
 print(density.output$summary)
 print(density.output$plot)

}# End of example code

# Author: Scott A. Hale (http://www.scotthale.net/)
# License: GPLv2
# If you use this in support of an academic publication, please cite:
#
#    Hale, S. A. (2014) Global Connectivity and Multilinguals in the Twitter Network. 
#    In Proceedings of the 2014 ACM Annual Conference on Human Factors in Computing Systems, 
#    ACM (Montreal, Canada).
#
# More details, related code, and the original academic paper using this code
# is available at http://www.scotthale.net/pubs/?chi2014
#

density.dataframe<-function(data,var,split=NA,n=512,adjust=1) {
	dfOutput<-data.frame()
	if (!is.na(split)) {
		if (is.factor(data[, split])) {
			vals<-levels(data[, split])
		} else {
			vals<-unique(data[, split])
		}
		for (val in vals) {
			tmp<-density(log10(data[data[,split]==val,var]),n=n,adjust=adjust,na.rm=TRUE)
			d<-data.frame(x=tmp$x,y=tmp$y,split=val)
			dfOutput<-rbind(dfOutput,d)
		}
	} else {
		tmp<-density(log10(data[,var]),n=n,adjust=adjust,na.rm=TRUE)
		dfOutput<-data.frame(x=tmp$x,y=tmp$y)
	}
	dfOutput$x10<-10^dfOutput$x
	
	#summary(log10(dfEditCountDensity$x10))
	#summary(dfEditCountDensity$x)
	
	return(dfOutput)
}

density.summary<-function(dataf,var,split) {	
	require(plyr)
	dfSummary <- ddply(dataf, c(split), function(df) {
		return(
		c(
			mean=mean(df[, var]),
			sd=sd(df[, var]),
			max=max(df[, var]),
			min=min(df[, var]),
			median=median(df[, var]),
			p25=summary(df[, var])[2],
			p75=summary(df[, var])[5]
		))})
	names(dfSummary)<-c("split","mean","sd","max","min","median","p25","p75")
	return(dfSummary)
}

density.plot<-function(dfDensity,dfSummary,title) {
	require(ggplot2)
	require(scales)

	plot <- ggplot(dfDensity,aes(x=x10,y=y,color=split,linetype=split)) + geom_path()
	plot <- plot + scale_x_log10(title,labels=comma) + scale_y_continuous("Density")
	plot <- plot + geom_vline(data=dfSummary, aes(xintercept=mean,color=split,linetype=split),size=1)#,linetype="dashed"
	plot <- plot +  theme_bw() + 
		theme(legend.title=element_blank(),legend.direction = "horizontal",legend.position = "bottom",
		legend.text=element_text(size=16),
		axis.title.x=element_text(size=18),axis.text.x=element_text(size=16),
		axis.title.y=element_text(size=18),axis.text.y=element_text(size=16))
	return(plot)
}

density.log<-function(data,var,split=NA,n=512,adjust=1,title="") {
	tmp.points <- density.dataframe(data,var,split,n,adjust)
	tmp.summary <- density.summary(data,var,split)
	return(list(
		points=tmp.points,
		summary=tmp.summary,
		plot=density.plot(tmp.points,tmp.summary,title)
	));
}