
#***************************************************************
#		Lab 3.1 Transplant Center 
#			Observational Analysis
#***************************************************************

#***************************************************************
#
#  Read in the data
#
#***************************************************************

#Set working directory

#Read data

r11xplant <- read.table("R11xplant.csv", sep = ",", header = T)
r11donor<-read.table("R11donor.csv", sep = ",", header = T)

uva <- read.table("UVAxplant.csv", sep = ",", header = T)
duke <- read.table("Dukexplant.csv", sep = ",", header = T)
mcv <- read.table("MCVxplant.csv", sep = ",", header = T)
unc <- read.table("UNC.csv", sep = ",", header = T)

#***************************************************************
#
# Soure the scatter plot matrix code and transplant plotting code
#
#***************************************************************

source("SPM_Panel.R")
source("Transplant.plots.R")

#***************************************************************
#
# Scatter plot matrix 
#
#***************************************************************

# Scatter plot matrix for liver transplants
liver <- data.frame(uva$Liver,duke$Liver,mcv$Liver,unc$Liver, r11donor$Liver)
uva.pairs(as.matrix(liver))

# Scatter plot matrix for kidney transplants
kidney <- data.frame(uva$Kidney,duke$Kidney,mcv$Kidney,unc$Kidney, r11donor$Kidney)
uva.pairs(as.matrix(kidney))

# Scatter plot matrix for pancreas transplants
pancreas<-data.frame(uva$Pancreas,duke$Pancreas,mcv$Pancreas,unc$Pancreas, r11donor$Pancreas)
uva.pairs(as.matrix(pancreas))

heart <- data.frame(uva$Heart,duke$Heart,mcv$Heart,unc$Heart, r11donor$Heart)
uva.pairs(as.matrix(heart))


#***************************************************************
#
#  donortype.plot
#
#***************************************************************

# Note in all the plots that follow we remove the 23rd observation
# (2010) since the data were not complete for that year.

# Lung transplants
donortype.plot(cbind(r11xplant$Lung_DD[-23], r11xplant$Lung_LD[-23], r11donor$Lung_DD[-23],r11donor$Lung_LD[-23]), title = "Lung", Year = seq(1988, 2009, 1)) 

# Liver transplants
donortype.plot(cbind(r11xplant$Liver_DD[-23], r11xplant$Liver_LD[-23], r11donor$Liver_DD[-23],r11donor$Liver_LD[-23]), title = "Liver", Year = seq(1988, 2009, 1)) 

# DD means deceased donor; LD means living donor
donortype.plot(cbind(r11xplant$Heart_DD[-23], r11xplant$Heart_LD[-23], r11donor$Heart_DD[-23],r11donor$Heart_LD[-23]), title = "Heart", Year = seq(1988, 2009, 1)) 
donortype.plot(cbind(r11xplant$Liver_DD[-23], r11xplant$Liver_LD[-23], r11donor$Liver_DD[-23],r11donor$Liver_LD[-23]), title = "Liver", Year = seq(1988, 2009, 1)) 


#***************************************************************
#
#  region.plot
#
#***************************************************************


#   Heart
region.plot(cbind(r11xplant$Heart[-23], r11donor$Heart[-23], uva$Heart[-23], unc$Heart[-23], mcv$Heart[-23], duke$Heart[-23]), title = "Heart", Year = seq(1988, 2009, 1))

#   Liver
region.plot(cbind(r11xplant$Liver[-23], r11donor$Liver[-23], uva$Liver[-23], unc$Liver[-23], mcv$Liver[-23], duke$Liver[-23]), title = "Liver", Year = seq(1988, 2009, 1))

#   Lung
region.plot(cbind(r11xplant$Lung[-23], r11donor$Lung[-23], uva$Lung[-23], unc$Lung[-23], mcv$Lung[-23], duke$Lung[-23]), title = "Liver", Year = seq(1988, 2009, 1))


#***************************************************************
#
#  center.plot
#
#***************************************************************

# Pancreas
center.plot(cbind( uva$Pancreas[-23], unc$Pancreas[-23], mcv$Pancreas[-23], duke$Pancreas[-23]), title = "Pancreas", Year = seq(1988, 2009, 1))

# Heart
center.plot(cbind( uva$Heart[-23], unc$Heart[-23], mcv$Heart[-23], duke$Heart[-23]), title = "Heart", Year = seq(1988, 2009, 1))

# Liver
center.plot(cbind( uva$Liver[-23], unc$Liver[-23], mcv$Liver[-23], duke$Liver[-23]), title = "Liver", Year = seq(1988, 2009, 1))

#***************************************************************
#
#		Lab 3.2	 Transplant Center 
#			Bootstrapping
#
#***************************************************************

#***************************************************************
#
#  Read in the data
#
#***************************************************************

# Set working directory

# Read data

r11xplant <- read.table("R11xplant.csv", sep = ",", header = T)
r11donor<-read.table("R11donor.csv", sep = ",", header = T)
uva <- read.table("UVAxplant.csv", sep = ",", header = T)
duke <- read.table("Dukexplant.csv", sep = ",", header = T)
mcv <- read.table("MCVxplant.csv", sep = ",", header = T)
unc <- read.table("UNC.csv", sep = ",", header = T)

# Source the bootstrapping functions

library(boot) #If you don't have this library, install it by: install.packages('boot')
source("TSbootfunctions.R")

#***************************************************************
#
# Bootstrap the differences
#
#***************************************************************

# UVa-MCV Differences in Kidney Transplants 
uva.kidney <- uva$Kidney
mcv.kidney <- mcv$Kidney
kid.diff <- ts(uva.kidney-mcv.kidney,1988,2010)

# recall paired t-test:
t.test(uva.kidney, mcv.kidney,paired=T)

# plot the differences:
ts.plot(kid.diff,ylab='UVa-MCV',main = "Difference in Number of Transplants, UVA-MCV")

#Bootstrap the differences: 
bs.mean <- function(x,i)
{
	return(mean(x[i]))
}

bs.kid.diff<-boot(kid.diff,bs.mean,R=2000)

plot(bs.kid.diff,index=1) 

boot.ci(bs.kid.diff,0.95,type=c('bca','perc'))

#***************************************************************
#
# Bootstrap Regression and Time Series
#
#***************************************************************

#Kidney Model:

#Linear Model:
uva.kid.lm <- lm(uva$Kidney~r11donor$Kidney)

summary(uva.kid.lm)

#Diag
par(mfrow=c(2,2))
plot(uva.kid.lm)
par(mfrow=c(1,1))

# Bootstrapping the LM

# Get the fitted values from the regression model
uva.kfit <- fitted(uva.kid.lm)

# Get the residuals from the regression model
uva.ke <- residuals(uva.kid.lm)

# Get the regression model 
uva.mod <- model.matrix(uva.kid.lm)

	
# Use the RTSB function to obtain the bootstrap distribution for the coefficients. This may take a few minutes.
# Be patient.	

uva.kid.boot <- RTSB(uva$Kidney, r11donor$Kidney, uva.kfit, uva.ke, uva.mod,2000)


# The estimates
uva.kid.boot
summary(uva.kid.lm)

# The 99% CI
boot.ci(uva.kid.boot, .99)

# The 95% CI

boot.ci(uva.kid.boot, .95)
	
# Plot the results for the coeffiecient for region 11 donors
plot(uva.kid.boot, index = 2) 

# A set of configurable plots
par(mfrow = c(1,2))
hist(uva.kid.boot$t[,2], main = "Region 11 Donors",xlab ="Coefficient Values",   col = "steelblue", breaks = 50)
qqnorm(uva.kid.boot$t[,2])
qqline(uva.kid.boot$t[,2])
par(mfrow = c(1,1))

# Bootstrapping TS

# Evaluating residual correlation
par(mfrow = c(1,2))
acf(uva.kid.lm$residuals)
pacf(uva.kid.lm$residuals)
par(mfrow = c(1,1))

# Fit an ar model to the residuals
ar(uva.kid.lm$residuals, method = "yule-walker") #AR(2)

# Add the ar model of the residuals to regression, linear model.
uva.kid2<- lm(uva$Kidney[3:23]~r11donor$Kidney[3:23]+ uva.kid.lm$residuals[2:22] + uva.kid.lm$residuals[1:21])

# The summary still show siginificance for the donor time series
summary(uva.kid2)

# Bootstrap the above time series model:

# Get the fitted values from the regression model
uva.kfit2 <- fitted(uva.kid2)
	
# Get the residuals from the regression model
uva.ke2 <- residuals(uva.kid2)

# Get the regression model 
uva.mod2 <- model.matrix(uva.kid2)
	
# Use the RTSB function to obtain the bootstrap  distribution for the coefficients. This may take a few minutes.
# Be patient.	
uva.boot2 <- RTSB(uva$Kidney[3:23], r11donor$Kidney[3:23], uva.kfit2, uva.ke2, uva.mod2,2000)
	
# The estimates
uva.boot2

# Plot the results for the coeffiecient for region 11 donors
plot(uva.boot2, index = 2)	

# Plot the results for the coeffiecient for time series components
plot(uva.boot2, index = 3)
plot(uva.boot2, index = 4)

# Lab 3.3

#***************************************************************
#
#  Read in the data
#
#***************************************************************

uva <- read.table("UVAxplant.csv", sep = ",", header = T)
mcv <- read.table("MCVxplant.csv", sep = ",", header = T)

# load the boot package

library(boot) 

# Source class code 

source("Transplant.plots.R")
source("TSbootfunctions.R")

#********************************************************************
#
# Difference on the kidney donors at UVA and MCV
#
#********************************************************************

# donor plot
uva.donors <- uva$Kidney_DD + uva$Kidney_LD
mcv.donors <- mcv$Kidney_DD + mcv$Kidney_LD

plot(uva$Year[-23], uva.donors[-23], main = "Kidney Donors", xlab = "Years", type = "l", col ="orange", ylab = "No. of Donors", ylim = c(10, max(mcv.donors)))
lines(uva$Year[-23], mcv.donors[-23])
legend(2004, 45, legend = c("UVA", "MCV"), col = c("orange", "black"), lwd = 2)

# Creating a times series for the diffence in number of donors
# at each medical center from 1988 to 2009.

kid.diff <- ts((uva.donors - mcv.donors), 1988, 2009)

# Statistics

mean(kid.diff)
sd(kid.diff)
var(kid.diff)

# The t test of the difference

t.test(uva.donors, mcv.donors)

# Wilcoxon test of the difference

wilcox.test(uva.donors, mcv.donors)

# Difference series

# A ts plot of the difference.

ts.plot(diff.kid, main = "Difference between Kidney Donors at MCV and UVA", ylab = "MCV - UVA")
abline(h=0)

# ACF and PACF of the difference time series.

par(mfrow = c(1,2))
acf(diff.kid[1:22], main = "ACF for MCV-UVA Kidney")
pacf(diff.kid[1:22], main = "PACF for MCV-UVA Kidney")
par(mfrow = c(1,1))

# Choosing an AR model

diff.ar <- ar(diff.kid[1:22], method = "yule-walker") #uses yule walker
diff.ar <- ar(diff.kid[1:22], method = "mle") # uses MLE

# How good is the model?

# Test for serial correlatoin. 

acf(na.omit(diff.ar$resid), main = "ACF for AR(1) of Center Differences")

# Test for Gaussian distribution.

qqnorm(na.omit(diff.ar$resid)); qqline(na.omit(diff.ar$resid))

# See how well the chosen model performs on AIC vs. others.

barplot(diff.ar$aic, main = "AIC for AR models of center differences", xlab = "Lag", ylab = "AIC", col = "steelblue")

# Get the predictions for 2010 and 2011.

diff.pred <- predict(diff.ar, n.ahead = 2)

# Plot the predictions with confidence intervals.

plot( c(diff.kid)~uva$Year[1:22], xlim = c(1988, 2012), ylim = c(-80, 30), type = "l", xlab = "Years", ylab = "MCV - UVA", main = "Forecasts for Differences in Kidney Transplants")
abline(h = 0)
lines(c(2009, 2010), c(diff.kid[22], diff.pred$pred[1]))
lines(c(2010, 2011), diff.pred$pred, type = "b", pch = 19)
segments(2010, diff.pred$pred[1], 2010, diff.pred$pred[1] + 2*diff.pred$se[1], lty = "dashed" )
segments(2010, diff.pred$pred[1], 2010, diff.pred$pred[1] - 2*diff.pred$se[1], lty = "dashed" )
segments(2011, diff.pred$pred[2], 2011, diff.pred$pred[2] + 2*diff.pred$se[2], lty = "dashed" )
segments(2011, diff.pred$pred[2], 2011, diff.pred$pred[2] - 2*diff.pred$se[2], lty = "dashed" )

#  Bootstrapping the difference

diff.boot <- TSB(diff.kid, diff.kid[22], 1000)

# Plotting the bootstrap results
	
par(mfrow = c(1,2))
	hist(diff.boot$t[,1], main = "UVA - MCV 2010 Forecast",xlab ="Forecast Values",   col = "steelblue", breaks = 50)
	qqnorm(diff.boot$t[,1])
	qqline(diff.boot$t[,1])
	par(mfrow = c(1,1))

# 95% bootstrap confidence interval for the forecast

tsboot.ci(diff.boot)
# Lab 3.4

#***************************************************************
#
#  Read in the data
#
#***************************************************************

uva <- read.table("UVAxplant.csv", sep = ",", header = T)
mcv <- read.table("MCVxplant.csv", sep = ",", header = T)
uvaeth <- read.table("UVAethnic.csv", sep = ",", header =T)
mcveth <- read.table("MCVethnic.csv", sep = ",", header =T)

# Source 
source("DemographicFunctions.R")
source("TSbootfunctions.R")
library(boot)

#*******************************************
#
#		Ethnic Diversity Analysis
#
#*******************************************

# Form data sets for the kidney transplants

uvaketh <- subdata("Kidney", uvaeth)
mcvketh <- subdata("Kidney", mcveth)


# Remove year 2010 and combine all ethnic groups other than white 
# into one category

uvake <- data.frame(uvaketh[-23, "Kidney.W"], Kidney.O = apply(uvaketh[-23,which(colnames(uvaketh) != "Kidney.W")], 1, sum))
mcvke <- data.frame(mcvketh[-23, "Kidney.W"], Kidney.O = apply(mcvketh[-23,which(colnames(mcvketh) != "Kidney.W")], 1, sum))


# Remove year 2010 and combine all ethnic groups other than white 
# into one category

uvake <- data.frame(Kidney.W = uvaketh[-23, "Kidney.W"], Kidney.O = apply(uvaketh[-23,which(colnames(uvaketh) != "Kidney.W")], 1, sum))
mcvke <- data.frame(Kidney.W = mcvketh[-23, "Kidney.W"], Kidney.O = apply(mcvketh[-23,which(colnames(mcvketh) != "Kidney.W")], 1, sum))

# UVA Ethnic plot

par(mfrow = c(1,2))
plot(uva$Year[-23], uvake[,1], type = "l", ylim = c(min(uvake), max(uvake)), xlab = "Year", ylab = "No. of Transplants", main = "UVA Transplants by Ethnic Group")
lines(uva$Year[-23], uvake[,2], col = "green")
legend(1990, 70, legend = c("White", "Other"), lwd = 2, col = c("black", "green"))

plot(mcv$Year[-23], mcvke[,1], type = "l", ylim = c(min(mcvke), max(mcvke)), xlab = "Year", ylab = "No. of Transplants", main = "MCV Transplants by Ethnic Group")
lines(mcv$Year[-23], mcvke[,2], col = "green")
legend(1990, 80, legend = c("White", "Other"), lwd = 2, col = c("black", "green"))
par(mfrow = c(1,1))

# Percent white 

uvapw <- (uvake[,1] * 100)/(apply(uvake[,1:2], 1, sum))
mcvpw <- (mcvke[,1] * 100)/(apply(mcvke[,1:2], 1, sum))

# Time series plots at both medical centers

plot(uva$Year[-23], uvapw, type = "l", ylim = c(10, 100), col = "orange", main = "Percent White", ylab = "Percent", xlab = "Year")
lines(uva$Year[-23], mcvpw, type = "l", col = "black")
legend(2001, 100, legend = c("UVA", "MCV"), lwd = 2, col = c("orange", "black"))

# ACF analysis

par(mfrow = c(1,2))
acf(uvapw, main = "ACF UVA PW")
acf(mcvpw, main = "ACF MCV PW")
par(mfrow = c(1,1))

# Classical tests

t.test(uvapw, mcvpw)
wilcox.test(uvapw, mcvpw)


#*******************************************
#
#		Response Surface Models
#
#*******************************************

# response surface for mcv

mcv.rs <- lm(mcv$Kidney[-23]~mcvpw)
summary(mcv.rs)

par(mfrow = c(2,2))
plot(mcv.rs)
par(mfrow = c(1,1))

par(mfrow = c(1,2))
acf(mcv.rs$resid)
pacf(mcv.rs$resid)
par(mfrow = c(1,1))

ar(mcv.rs$resid, method = "yule-walker")

# clearly a time series is appropriate

mcv.rs2 <- lm(mcv$Kidney[-c(1,23)] ~ mcvpw[-1] + mcv.rs$resid[1:21])
summary(mcv.rs2)

par(mfrow = c(2,2))
plot(mcv.rs2)
par(mfrow = c(1,1))

# check acf

acf(mcv.rs2$resid)

# response surface for uva

uva.rs <- lm(uva$Kidney[-23]~uvapw)
summary(uva.rs)

par(mfrow = c(2,2))
plot(uva.rs)
par(mfrow = c(1,1))

par(mfrow = c(1,2))
acf(uva.rs$resid, main = "ACF for RS")
pacf(uva.rs$resid, main = "PACF for RS")
par(mfrow = c(1,1))

# look for estimated model

ar(uva.rs$resid)

# clearly a time series is appropriate

uva.rs2 <- lm(uva$Kidney[-c(1,23)] ~ uvapw[-1] + uva.rs$resid[1:21])
summary(uva.rs2)

par(mfrow = c(2,2))
plot(uva.rs2)
par(mfrow = c(1,1))

# check acf

par(mfrow = c(1,2))
acf(uva.rs2$resid, main = "ACF Time Series RS")
pacf(uva.rs2$resid, main = "PACF Time Series RS")
par(mfrow = c(1,1))


#*******************************************
#
#  Optimization with Response Surface Models
#
#*******************************************

# Choose two design points, say 70 and 60. 
# what does that do for our predicted no. of xplants?

# Create variables for the model

res2 <- uva$Kidney[-c(1,23)]
v1 <- uvapw[-1]
v2 <- uva.rs$resid[1:21]

# Run the RS model

uva.rs2 <- lm(res2~v1 + v2)

# Get the predictions at the two design points

uva.pred1 <- predict(uva.rs2, newdata = data.frame(v1 = 80, v2 = uva.rs$resid[22]), se.fit = T)
uva.pred2 <- predict(uva.rs2, newdata = data.frame(v1 = 50, v2 = uva.rs$resid[22]), se.fit = T)

# Plot of predictions for each value of the pw and the confidence intervals

plot( uva$Kidney[1:22]~uva$Year[1:22], xlim = c(1988, 2011), type = "l", xlab = "Years", ylab = "No. of Kidney", main = "Predicted Numbers of UVA Kidney Transplants")

points(2010, uva.pred1$fit, pch = 19)
segments(2010, uva.pred1$fit,  2010, uva.pred1$fit + 2*uva.pred1$se.fit[1], lty = "dashed" )
segments(2010, uva.pred1$fit,  2010, uva.pred1$fit - 2*uva.pred1$se.fit[1], lty = "dashed" )

points(2010.5, uva.pred2$fit, pch = 19, col = "orange")
segments(2010.5, uva.pred2$fit,  2010.5, uva.pred2$fit + 2*uva.pred2$se.fit[1], lty = "dashed", col = "orange" )
segments(2010.5, uva.pred2$fit,  2010.5, uva.pred2$fit - 2*uva.pred2$se.fit[1], lty = "dashed" , col = "orange")

legend(1990, 100, legend = c("PW = 80", "PW = 50"), pch = 19, col = c("black", "orange"))

# Getting bootstrapped confidence interval for the predicted values

#**********************************************************************
#
#        Bootstrapped Estimates for the Response Surface Models
#
#**********************************************************************

# Getting the bootstrap estimate for the regression forecast
# using the RFB() function in the class package TSbootfunctions.R

uva.rs2.boot <- RFB(cbind(res2, v1, v2), model = uva.rs2, ndata = data.frame(v1 = 50, v2 = uva.rs$resid[22]), num = 2000 )

# Bootstrap plot 

plot(uva.rs2.boot, index = 1)

# The mean and CI for the bootstrap

mq <- mean(uva.rs2.boot$t[,1])
bq1 <- quantile(uva.rs2.boot$t[,1], .05)
bq2 <- quantile(uva.rs2.boot$t[,1], .95)

# The plot of the bootstrap estimate with the parametric model-based estimate

plot( uva$Kidney[1:22]~uva$Year[1:22], xlim = c(1988, 2011), type = "l", xlab = "Years", ylab = "No. of Kidney", main = "Predicted Numbers of UVA Kidney Transplants")

points(2010, mq, pch = 19)
segments(2010, mq,  2010, bq1, lty = "dashed" )
segments(2010, mq,  2010, bq2, lty = "dashed" )

points(2010.5, uva.pred2$fit, pch = 19, col = "orange")
segments(2010.5, uva.pred2$fit,  2010.5, uva.pred2$fit + 2*uva.pred2$se.fit[1], lty = "dashed", col = "orange" )
segments(2010.5, uva.pred2$fit,  2010.5, uva.pred2$fit - 2*uva.pred2$se.fit[1], lty = "dashed" , col = "orange")

legend(1990, 100, legend = c("Bootstrap", "Gaussian"), pch = 19, col = c("black", "orange"))

Lab 3 Compiled Code
# Lab 3 Compiled Code

#***************************************************************
#
#  Kidney Transplant EISE RSM Analysis
#		Iteration 2
#
#***************************************************************

#***************************************************************
#
#  Read in the data
#
#***************************************************************

uva <- read.table("UVAxplant.csv", sep = ",", header = T)
mcv <- read.table("MCVxplant.csv", sep = ",", header = T)

uvaeth <- read.table("UVAethnic.csv", sep = ",", header =T)
mcveth <- read.table("MCVethnic.csv", sep = ",", header =T)

# Source 
source("Transplant.plots.R")
source("DemographicFunctions.R")
source("TSbootfunctions.R")

# We begin by viewing the data. Explain our purpose in doing this?
# To view the ethnic composition of the kidney transplant patients we need to form new data sets for the kidney transplants patient population.
# We use the subdata functions to do this. (see lab 3.4 notes)

uvaketh <- subdata("Kidney", uvaeth)
mcvketh <- subdata("Kidney", mcveth)

# Remove year 2010 and combine all ethnic groups other than white
# into one category

uvake <- data.frame(uvaketh[-23, "Kidney.W"], Kidney.O = apply(uvaketh[-23,which(colnames(uvaketh) != "Kidney.W")], 1, sum))
mcvke <- data.frame(mcvketh[-23, "Kidney.W"], Kidney.O = apply(mcvketh[-23,which(colnames(mcvketh) != "Kidney.W")], 1, sum))


# Remove year 2010 and combine all ethnic groups other than white into one category

uvake <- data.frame(Kidney.W = uvaketh[-23, "Kidney.W"], Kidney.O = apply(uvaketh[-23,which(colnames(uvaketh) != "Kidney.W")], 1, sum))
mcvke <- data.frame(Kidney.W = mcvketh[-23, "Kidney.W"], Kidney.O = apply(mcvketh[-23,which(colnames(mcvketh) != "Kidney.W")], 1, sum))

# UVA Ethnic plot

par(mfrow = c(1,2))
plot(uva$Year[-23], uvake[,1], type = "l", ylim = c(min(uvake), max(uvake)), xlab = "Year", ylab = "No. of Transplants", main = "UVA Transplants by Ethnic Group")
lines(uva$Year[-23], uvake[,2], col = "green")
legend(1990, 70, legend = c("White", "Other"), lwd = 2, col = c("black", "green"))

plot(mcv$Year[-23], mcvke[,1], type = "l", ylim = c(min(mcvke), max(mcvke)), xlab = "Year", ylab = "No. of Transplants", main = "MCV Transplants by Ethnic Group")
lines(mcv$Year[-23], mcvke[,2], col = "green")
legend(1990, 80, legend = c("White", "Other"), lwd = 2, col = c("black", "green"))
par(mfrow = c(1,1))

# Percent white

uvapw <- (uvake[,1] * 100)/(apply(uvake[,1:2], 1, sum))
mcvpw <- (mcvke[,1] * 100)/(apply(mcvke[,1:2], 1, sum))

# Time series plots at both medical centers

plot(uva$Year[-23], uvapw, type = "l", ylim = c(10, 100), col = "orange", main = "Percent White", ylab = "Percent", xlab = "Year")
lines(uva$Year[-23], mcvpw, type = "l", col = "black")
legend(2001, 100, legend = c("UVA", "MCV"), lwd = 2, col = c("orange", "black"))

# Based on your observations what conjecture do you make to start iteration 2? What is the associated hypothesis?
# How will we test this hypothesis? Use Lab3.4.R to conduct these tests.

# Classical tests

t.test(uvapw, mcvpw)
wilcox.test(uvapw, mcvpw)

# Based on your test results from interation 2
# what conjecture do you make to start
# iteration 3? What is the associated hypothesis?
# How will you test this hypothesis?
#
#######

 
#######
#
# In Lab 3.4 and in lecture 24 we formulated two models for our
# test. In lecture 25 we formulated one model. We will do that here.
# To use one model we need first to remove any serial correlation in
# kidney transplant series. Do this by finding the appropriate models
# for the number of kidney transplants at UVA (one time series model) and
# the number of kidney transplants at MCV (another time series model)
# Test your models to see if they have removed the serial correlation.
# Finally, use var.test() to  see if the variances for centers are the same.
#
#######
uva.ar <- ar(uva$Kidney[-23])

mcv.ar <- ar(mcv$Kidney[-23])
par(mfrow = c(1,2))
acf(uva.ar$resid[-1])
acf(mcv.ar$resid[-1])
par(mfrow = c(1,1))

# Test the variances of the residuals for each center.

var.test(uva.ar$resid[-1], mcv.ar$resid[-1])

#######
#
# Next create a data frame that concatenates the residuals from the
# time series series model (UVA first) for one variable (call it res).
# Add the variable pw for percent white for the observations. Then add
# a variable (UVA) which takes the value 1 for UVA observations and
# 0 for MCV. Call this data frame eth.
#
#######

eth <- data.frame(res = c(uva.ar$resid, mcv.ar$resid), pw = c(uvapw, mcvpw), UVA = c(rep(1,22), rep(0,22)))

#######
#
# Next create a linear model with res as the response and pw and UVA as predictors.
# Evaluate this model.
#
#######

par(mfrow = c(2,2))
plot(eth.rs2)
par(mfrow = c(1,1))

par(mfrow = c(1,2))
acf(eth.rs2$resid)
pacf(eth.rs2$resid)
par(mfrow = c(1,1))

#######
#
# Now use your model to predict the number of kidney transplants
# with pw values of 70 and 60. Then plot these predictions with their
# associated confidence intervals. Use these results to test you conjectured
# hypothesis.
#
#######

eth.p1 <- predict(eth.rs2, newdata = data.frame(pw = 70, UVA = 1), se = T)


eth.p2 <- predict(eth.rs2, newdata = data.frame(pw = 60, UVA = 1), se = T)

eth.pred1 <- predict(uva.ar, newdata = uva$Kidney[22])

eth.pred2 <- predict(uva.ar, newdata = uva$Kidney[22])

ep1 <- eth.pred1$pred + eth.p1$fit
ep2 <- eth.pred2$pred + eth.p2$fit

# Plot of predictions for each value of the pw and the confidence intervals

plot( uva$Kidney[1:22]~uva$Year[-23], xlim = c(1988, 2011), type = "l", xlab = "Years", ylab = "No. of Kidney", main = "Predicted Numbers of UVA Kidney Transplants", ylim = c(20,120))

points(2010, ep1, pch = 19)

segments(2010, ep1,  2010, ep1 + 2*(eth.p1$se.fit+eth.p1$se.fit), lty = "dashed" )
segments(2010, ep1,  2010, ep1 - 2*(eth.p1$se.fit+eth.p1$se.fit), lty = "dashed" )

points(2010.5, ep2, pch = 19, col = "orange")

segments(2010.5, ep2,  2010.5, ep2 + 2*(eth.p2$se.fit+eth.p2$se.fit), lty = "dashed", col = "orange" )
segments(2010.5, ep2,  2010.5, ep2 - 2*(eth.p2$se.fit+eth.p2$se.fit), lty = "dashed" , col = "orange")

legend(1990, 100, legend = c("PW = 70", "PW = 60"), pch = 19, col = c("black", "orange"))

#######
#
# Next use bootstrapping to predict the number of kidney transplants
# with pw values of 70 and 60 and their associated confidence intervals.
# Use these results to test you conjectured
# hypothesis. For this use the RFB() function.
# Finish iteration 3 with your conclusions.
#
#######

eth.boot1 <- RFB(na.omit(eth), model = eth.rs2, ndata = data.frame(pw = 70, UVA =1), num = 2000 )

eth.boot2 <- RFB(na.omit(eth), model = eth.rs2, ndata = data.frame(pw = 60, UVA =1), num = 2000 )

# Bootstrap plot

plot(eth.boot1, index = 1)

# Get the predictions at the two design points

eth.p1 <- mean(eth.boot1$t[,1])

eth.p2 <- mean(eth.boot2$t[,1])

eth.pred1 <- predict(uva.ar, newdata = uva$Kidney[22])

eth.pred2 <- predict(uva.ar, newdata = uva$Kidney[22])

ep1 <- eth.pred1$pred + eth.p1
ep2 <- eth.pred2$pred + eth.p2

blq1 <- eth.p1 - quantile(eth.boot1$t[,1], .05)
buq1 <- quantile(eth.boot1$t[,1], .95) - eth.p1

blq2 <- eth.p2 - quantile(eth.boot2$t[,1], .05)
buq2 <- quantile(eth.boot2$t[,1], .95) - eth.p2

# Plot of predictions for each value of the pw and the confidence intervals

plot( uva$Kidney[1:22]~uva$Year[-23], xlim = c(1988, 2011), type = "l", xlab = "Years", ylab = "No. of Kidney", main = "Bootstrap Predicted Numbers of UVA Kidney Transplants", ylim = c(20,130))

points(2010, ep1, pch = 19)

segments(2010, ep1,  2010, ep1 + 2*eth.pred1$se[1] + buq1, lty = "dashed" )
segments(2010, ep1,  2010, ep1 - 2*eth.pred1$se[1] - blq1, lty = "dashed" )

points(2010.5, ep2, pch = 19, col = "orange")

segments(2010.5, ep2,  2010.5, ep2 + 2*eth.pred2$se[1]+ buq2, lty = "dashed", col = "orange" )
segments(2010.5, ep2,  2010.5, ep2 - 2*eth.pred2$se[1] - blq2, lty = "dashed" , col = "orange")

legend(1990, 100, legend = c("PW = 70", "PW = 60"), pch = 19, col = c("black", "orange"))

#***************************************************************
#
#  Liver Transplant EISE RSM Analysis
#   	 Iteration 1
#
#***************************************************************

#***************************************************************
#
#  Read in the additional data
#
#***************************************************************

r11donor <- read.table("R11donor.csv", sep = ",", header = T)


r11xplant <- read.table("R11xplant.csv", sep = ",", header = T)

duke <- read.table("Dukexplant.csv", sep = ",", header = T)

unc <- read.table("UNC.csv", sep = ",", header = T)


# Start by observing the data using a region plot and
# a center plot.

par(mfrow = c(1,2))
acf(uva$Liver[1:17], main = "UVA 1988-2004")
acf( uva$Liver[18:22], main = "UVA 2005-2009")
acf(mfrow = c(1,1))

t.test(uva$Liver[1:17], uva$Liver[18:22])

wilcox.test(uva$Liver[1:17], uva$Liver[18:22])

#######
#
# To test your hypothesis you will need a linear model with variable for
# the Roanoke center. Do you understand why?
# Create this variable (Roan) which takes on the value 1 for years 2005-2009
# and 0 for years 1988-2004.
# Use this variable as the predictor in a linear model to predict then of liver
# transplants at UVA until 2009. Perform diagnostic tests of your model.
#  
#######

Roan <- c(rep(0, 17), rep(1, 5))

summary(uva.rs1)

par(mfrow = c(2,2))
plot(uva.rs1)
par(mfrow = c(1,1))

par(mfrow = c(1,2))
acf(uva.rs1$resid, main = "ACF for Resid from RS1")
pacf(uva.rs1$resid, main = "PACF for Resid from RS1")
acf(mfrow = c(1,1))

#######
#
# If your model showed serial correlation, then modify it to eliminate this
# correlation. Then test your new model.
#  
#######

uva.rs2 <- lm(uva$Liver[2:22]~Roan[2:22] + uva.rs1$resid[1:21])

summary(uva.rs2)

par(mfrow = c(2,2))
plot(uva.rs2)
par(mfrow = c(1,1))


par(mfrow = c(1,2))
acf(uva.rs2$resid, main = "ACF for Resid from RS2")
pacf(uva.rs2$resid, main = "PACF for Resid from RS2")
par(mfrow = c(1,1))



#######
#
# Now create a time series model to estimate the
# number of transplants UVA would have performed without the
# Roanoke center.
#  
#######

uva.ar <- ar(uva$Liver[1:17])

barplot(uva.ar$aic)

# Getting the 5 period prediction

uva.pred1 <- predict(uva.ar, n.ahead = 5)

# Plot the prediction

plot(uva$Year[1:17], uva$Liver[1:17], xlim = c(1988,2009), ylim = c(0, uva.pred1$pred[5] + 2*uva.pred1$se[5]), main = "Prediction without the Roanoke Center", type = "l", ylab = "Number of Transplants", xlab = "Year")

points(2005:2009, uva.pred1$pred, pch = 19)

segments(2005:2009, uva.pred1$pred,  2005:2009, uva.pred1$pred + 2*uva.pred1$se, lty = "dashed" )
segments(2005:2009, uva.pred1$pred,  2005:2009, uva.pred1$pred - 2*uva.pred1$se, lty = "dashed" )

# Plotting the difference

plot(uva$Year[1:22], uva$Liver[1:22], main = "Estimated Added Value of the New Center", type = "l", xlim = c(1988, 2014), ylim = c(0, uva$Liver[22] + diff[5] + 2*uva.pred1$se[5]),  ylab = "Number of Transplants", xlab = "Year")
points(2010:2014, uva$Liver[22]+ diff, pch = 19, type = "b")

segments(2010:2014, uva$Liver[22]+ diff,  2010:2014, uva$Liver[22]+ diff + 2*uva.pred1$se, lty = "dashed" )
segments(2010:2014, uva$Liver[22]+ diff,  2010:2014, uva$Liver[22]+ diff - 2*uva.pred1$se, lty = "dashed" )


#######
#
# Finally obtain the bootstrap estimates (you will need to use
# the TSF function) for the additional
# transplants that a new center could provide. Plot these
# estimates to show the value of a new center like the Roanoke center.
# What are the results of your hypothesis test for this iteration?
# What do you conclude about your conjecture?
#  
#######


uva.boot.mean <- TSF(uva$Liver[1:17], uva$Liver[17], 1000, n.ahead = 5)
uva.boot.se <- TSF(uva$Liver[1:17], uva$Liver[17], 1000, n.ahead = 5, func = "S")

mq <- apply(uva.boot.mean$t, 2, mean)
bse <- apply(uva.boot.se$t, 2, mean)
diff <- uva$Liver[18:22] - mq

plot(uva$Year[1:22], uva$Liver[1:22], main = "Bootstrap Estimated Added Value of the New Center", type = "l", xlim = c(1988, 2014), ylim = c(0, uva$Liver[22] + diff[5] + 2*uva.pred1$se[5]),  ylab = "Number of Transplants", xlab = "Year")
points(2010:2014, uva$Liver[22]+ diff, pch = 19, type = "b")

segments(2010:2014, uva$Liver[22]+ diff,  2010:2014, uva$Liver[22]+ diff + 2*bse, lty = "dashed" )
segments(2010:2014, uva$Liver[22]+ diff,  2010:2014, uva$Liver[22]+ diff - 2*bse, lty = "dashed" )
