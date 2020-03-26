#set to use 30 GB memory
options(java.parameters = "-Xmx30g")

library(foreign)
library(bartMachine)
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

#set to run on 5 cores
set_bart_machine_num_cores(5)

setwd("C:/Users/local-admin/Documents/COVID-19")
#Read shapefile
shp <- readOGR(dsn="RKI_with_All.shp")
#Read DBF
#d <- read.dbf("RKI_with_All.dbf")
#View(shp)

#Seperate dependent and independent variables
y <- log(shp$cs__100)
shp$y = y
hist(y, breaks=50)
#summary(y)
#hist(y,breaks=50)
#Change negative infinite values [log(0)=-Inf] to 0
y[mapply(is.infinite, y)] <- 0

#Select independent variables
Xa <- data.frame(shp)
X <- Xa[c(25,34,40:363)]
rm(Xa)
#Remove unwanted variables
X$S_109 <- as.factor(X$S_109)
rm(Xb)
#Check for factors
names(Filter(is.factor, X))

#Map case counts
#tm_shape(shp) + tm_polygons(col="y", title="Cumulative incidence (log2), 25.03", breaks=c(0,1,2,3,4,5,6,7,8,9,10), palette="Greys") + tm_layout(inner.margins=0.03, legend.outside = T, bg.color = "white")

#bm <- bartMachine(X = X, y = y)
#summary(bm)
#Model will exhibit serious overfitting! Proceed to assess with k-folds CV.

#k-folds CV
#k_fold <- k_fold_cv(X=X, y=y, k_folds=5)
#k_fold$PseudoRsq

#select best-performing variables
#windows()
#varselect <- var_selection_by_permute(bm)
#print(varselect$important_vars_local_names)
#print(varselect$important_vars_global_max_names)
#print(varselect$important_vars_global_se_names) 

X2 <- subset(X, select=c(
  S_117, 
  S_119, 
  S_051, 
  S_104, 
  Rch_den, 
  S_125, 
  nigh_dn, 
  Ostdtsc, 
  S_018, 
  EWZ, 
  S_004, 
  S_092, 
  S_054, 
  BL))

#BARTmachine with selected variables, built using CV parameterisation
bm_cv <- bartMachineCV(X=X2, y=y)
summary(bm_cv)
plot_convergence_diagnostics(bm_cv)
assump <- check_bart_error_assumptions(bm_cv)
plot_y_vs_yhat(bm_cv,credible_intervals = T)
get_sigsqs(bm_cv, plot_hist = T, plot_sigma = T)
vi_splits <- investigate_var_importance(bm_cv, type="splits", num_replicates_for_avg = 50, num_trees_bottleneck = 20)
vi_trees <- investigate_var_importance(bm_cv, type="trees", num_replicates_for_avg = 50, num_trees_bottleneck = 20)

intinv <- interaction_investigator(bm_cv)
#Plot interaction importance as heat map

######## BAR CHART ##############
#Create a nice-looking bar chart of vi_splits and vi_trees! Include Hessen, but then include only one bar for the mean value for the remaining 17 Bundesländer.

########## PD PLOTS ############
pd1 <- pd_plot(bm_cv, "S_117")

############ Residual Map ###################
# join residuals to shapefile, then map residuals
shp$resid <- bm_cv$residuals
tm_shape(shp) + tm_polygons(col="resid", title="BART Machine Residuals (log incidence rate)", breaks=c(-4,-3,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,4), midpoint=NA, palette="RdGy") + tm_layout(inner.margins=0.03, legend.outside = T, bg.color = "white")
#RColorBrewer::brewer.pal.info

save.image("BART_COVID19_Final.RData")
