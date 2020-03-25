library(foreign)
library(rgdal)
library(RColorBrewer)
library(tmap)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(pscl)
library(Matrix)
library(lme4)

#Set workspace
setwd("D:/SARS-CoV-2/MasterFile")
#Read shapefile
shp <- readOGR(dsn="D:/SARS-CoV-2/MasterFile/RKI_with_All.shp")
#shpdf <- data.frame(shp)
#Read data
d <- read.dbf("RKI_with_All.dbf")

#independent variable: crude incidence rate
  d$y <- log2(d$cs__100)
#Change negative infinite values to 0
d$y[mapply(is.infinite, d$y)] <- 0
#Select independent variables
d <- d[c(29,48:317)]
#Remove unwanted variables
d <- subset(Xa, select=-c(S_LK))
#Check for factors, due to decimal-comma conversion
names(Filter(is.factor, d))

#Map case counts
tm_shape(shp) + tm_polygons(col="cases", title="Cumulative incidence, 23.03", breaks=c(0,2,5,10,20,50,100,500,750,1000,1500,2000), palette="Greys") + tm_layout(inner.margins=0.03, legend.outside = T, bg.color = "white")

#Poisson Model
m <- lmer(y ~ S_117 + S_051 + S_119 + S_108 + S_125 + S_004 + hair_pp + S_153 + cc_km2 + nigh_pp + EWZ + S_022 + S_130 + S_104 + S_011 + sch_pop + S_142 + S_054 + ff_pop + (1|S_109), 
          data=d, na.action = na.omit)
tab_model(m, show.se=T, show.std=T, p.style = "a", show.r2 = F, col.order = c("est", "se", "ci", "std.est", "std.se"))

#Compute pseudo R^2 measure for poisson models
#pseudoR2 = 1-(m$deviance/m$null.deviance)

#Post-hoc evaluation
#Response residuals - variance should increase with fitted value
pred <- predict.glm(m, type="response")
resp_resid <- (shp$cases-pred)
hist(resp_resid, breaks=50, main="Response Residuals", xlab="Fitted Value", ylab = "P(Response Residual)", freq = FALSE)
plot(pred, resp_resid, main="Response Residuals", xlab="Fitted Value", ylab="Response Residual")
plot(pred, resp_resid, main="Response Residuals, zoomed-in", xlab="Fitted Value", ylab="Response Residual", ylim=c(-100,100), xlim=c(0,100), abline(h=0, col="Red"))
#Pearson Redsiduals
pearsonResid <- resid(m, type="pear")
hist(pearsonResid, breaks=50, main="Pearson Residuals", xlab="Fitted Value", ylab = "P(Pearson Residual)", xlim=c(-20,120), freq = FALSE)
plot(pred, pearsonResid, main="Pearson Residuals", xlab="Fitted Value", ylab="Pearson Residual Residual")
plot(pred, pearsonResid, main="Pearson Residuals, zoomed-in", xlab="Fitted Value", ylab="Pearson Residual", ylim=c(-10,15), xlim=c(0,125), abline(h=0,col="Red"))

#Map residuals, examine results for clustering
shp$pearsonResid <- pearsonResid
tm_shape(shp) + tm_polygons(col="pearsonResid", title="Pearson Residuals", style="sd" , palette="Greys") + tm_layout(inner.margins=0.03, legend.outside = T, bg.color = "white")

#Print pseudo-R^2 statistics in Console window
round(pR2(m),digits=3)

#GLM with mixed effects
#m2 <- glmer(shpdf$cases ~ Einwohner_.6 + Fertilitae + Wahlbeteil+ Abhaengige.1 + (1|BEZ), data = shpdf, family = poisson(link="log"), glmerControl(optimizer = "bobyqa", nAGQ=0, optCtrl = list(maxfun=2e6)))
#tab_model(m2)