#set to use 30 GB memory
options(java.parameters = "-Xmx30g")

library(foreign)
library(bartMachine)
library(foreign)
library(rgdal)
library(RColorBrewer)
library(digest)
library(tmap)
library(pscl)

#set to run on 5 cores
set_bart_machine_num_cores(12)

setwd("C:/Users/local-admin/Documents/COVID-19")
#Read shapefile
shp <- readOGR(dsn="RKI_with_All_withDblTime.shp")
#Read DBF
#d <- read.dbf("RKI_with_All.dbf")


#Seperate dependent and independent variables
y <- log(as.numeric(gsub(",", ".", gsub("\\.", "", shp$DBL_TIM))))
y2 <- log(shp$cs__100)
#Change negative infinite values [log(0)=-Inf] to 0
y[mapply(is.infinite, y)] <- 0
y2[mapply(is.infinite, y2)] <- 0
shp$y = y
shp$y2 = y2
hist(y, breaks=50)
summary(y)
hist(y2, breaks=50)
summary(y2)


#Select independent variables
X<- data.frame(shp[c(25,34,40:364)])
#convert to factors
X$S_109 <- as.factor(X$S_109)
X$Ostdtsc <- as.factor(X$Ostdtsc)
#Check for factors
names(Filter(is.factor, X))




############### MAPPING ###################
#Map crude rates
tm_shape(shp) + 
  tm_polygons(col="y2", title="Crude Incidence Rate (log)", 
  #breaks=c(0,1,2,3,4,5,6,7,8,9,10), 
  palette="Greys") +
  tm_layout(inner.margins=0.03, legend.outside = T, bg.color = "white")

# #Calculate spatial lag term
# library(spdep)
# #queen contiguity
# nb <- poly2nb(shp,queen=T)
# lw <- nb2listw(nb, style = "W", zero.policy = T)
# 
# #distance band
coo <- coordinates(shp)
s.dist  <-  dnearneigh(coo, 0, 0.5) #units=decimal degrees
lw <- nb2listw(s.dist, style = "W", zero.policy = T)
# 
# time.lag <- lag.listw(lw, shp$y)
# shp$time.lag <- time.lag
# 



# #Map lagged rates
# tm_shape(shp) + tm_polygons(col="time.lag", title="Lag Doubling Time (log-days)", 
#   #breaks=c(0,1,2,3,4,5,6,7,8,9,10), 
#   palette="Greys") + 
#   tm_layout(inner.margins=0.03, legend.outside = T, bg.color = "white")
# 
# #Compare lagged rates to rates
# lagmodel <- lm(shp$time.lag ~ shp$y)
# summary(lagmodel)
# plot(shp$time.lag.lag, shp$y, pch=19, las=1, xlim=c(0,6), ylim=c(0,6), xlab="Lagged log-Rate", ylab="Observed log-Rate")
# 
# #Moran's I test for spatial autocorrelation, with monte carlo
# mc <- moran.mc(y,lw, 1000)
# mc
# plot(mc, main="", las=1)




################ BART MACHINE ###############

#k-folds CV
#k_fold <- k_fold_cv(X=X, y=y, k_folds=5)
#k_fold$PseudoRsq

#select best-performing variables
#windows()
#varselect <- var_selection_by_permute(bm_cv)
#print(varselect$important_vars_local_names)
#print(varselect$important_vars_global_max_names)
#print(varselect$important_vars_global_se_names) 




#Select variables for final model
X2 <- subset(X, select=c(
  #political units
  S_109, #rural/urban
  #population characteristics
  EWZ, #population
  Pop_Den, #population density
  S_004, #youth unemployment
  S_035, #residents ages 65+
  S_018, #pct of unemployed who are older
  S_051, #voter participation
  #S_054, #apprenticeship positions
  S_092, #persons in need of care
  S_104, #income tax
  S_069, #median income
  #S_113, #population density
  S_119, #persons 65+ receiving special assistance
  #S_124, #proportion of persons in community dependence
  #S_125, #proportion of persons in community dependence who with 3+ children
  #Access to transport infrastructure
  S_129, #airport access, driving time
  S_130, #IC train station access
  S_131, #Oberzentren access
  S_132, #Mittelzentren access
  S_142, #incoming commuters
  #built environment characteristics
  Rch_den, #church density
  hair_dn, #hairdresser density
  kid_den, #kita density
  play_dn, #playground density
  hosp_dn, #hosp density
  doc_den, #doctor office density
  cc_den, #community centre density
  sch_den, #school density
  nigh_dn, #nightclub density
  pub_den #pub density
  ))

##BARTmachine with selected variables, built using CV parameterisation
bm_cv <- bartMachineCV(X=X2, y=y2)
summary(bm_cv)

# Winning bartMachine: k: 2 nu, q: 10, 0.75 m: 200
#bm_cv <- bartMachine(X = X2, y = y, k=2, nu=3, q=0.75, num_trees=200, serialize = T)
#summary(bm_cv)
plot_convergence_diagnostics(bm_cv)
assump <- check_bart_error_assumptions(bm_cv)

#plot model fit
plot_y_vs_yhat(bm_cv,credible_intervals = T)
get_sigsqs(bm_cv, plot_hist = T, plot_sigma = T)

vi_splits <- investigate_var_importance(bm_cv, type="splits", num_replicates_for_avg = 50, num_trees_bottleneck = 20)
vi_trees <- investigate_var_importance(bm_cv, type="trees", num_replicates_for_avg = 50, num_trees_bottleneck = 20)

intinv <- interaction_investigator(bm_cv)
#Plot interaction importance as heat map 

######## BAR CHART ##############
#Create a nice-looking bar chart of vi_splits and vi_trees! Include Hessen, but then include only one bar for the mean value for the remaining 17 Bundesländer.

########## PD PLOTS ############
#SES
#pd_EWZ <- pd_plot(bm_cv, "EWZ")
pd_PopDen <- pd_plot(bm_cv, "Pop_Den") 
pd_S_004 <- pd_plot(bm_cv, "S_004") 
pd_S_035 <- pd_plot(bm_cv, "S_035")
pd_S_051 <- pd_plot(bm_cv, "S_051")
#pd_S_054 <- pd_plot(bm_cv, "S_054")
#pd_S_092 <- pd_plot(bm_cv, "S_092")
#pd_S_104 <- pd_plot(bm_cv, "S_104")
pd_S_069 <- pd_plot(bm_cv, "S_069")
pd_S_113 <- pd_plot(bm_cv, "S_113")
pd_S_119 <- pd_plot(bm_cv, "S_119")
#pd_S_124 <- pd_plot(bm_cv, "S_124")
#pd_S_125 <- pd_plot(bm_cv, "S_125")
#Transportation Infrastructure
pd_S_129 <- pd_plot(bm_cv, "S_129")
#pd_S_130 <- pd_plot(bm_cv, "S_130")
#Ober-/Mittelzentren
pd_S_131 <- pd_plot(bm_cv, "S_131")
pd_S_142 <- pd_plot(bm_cv, "S_142")
#pd_S_132 <- pd_plot(bm_cv, "S_132")
pd_Rch_den <- pd_plot(bm_cv, "Rch_den")
#pd_hair_dn <- pd_plot(bm_cv, "hair_dn") #wrong directionality!
pd_kid_den <- pd_plot(bm_cv, "kid_den") 
pd_play_dn <- pd_plot(bm_cv, "play_dn") #wrong directionality!
pd_hosp_dn <- pd_plot(bm_cv, "hosp_dn")
pd_doc_den <- pd_plot(bm_cv, "doc_den")
pd_cc_den <- pd_plot(bm_cv, "cc_den") #strange pattern...
pd_sch_den <- pd_plot(bm_cv, "sch_den")
pd_nigh_dn <- pd_plot(bm_cv, "nigh_dn")
#pd_pub_den <- pd_plot(bm_cv, "pub_den")
#pd_bar_den <- pd_plot(bm_cv, "bar_km2") #wrong directionality!

############ Residual Map ###################
# join residuals to shapefile, then map residuals
#RColorBrewer::brewer.pal.info
shp$resid <- bm_cv$residuals
tm_shape(shp) + tm_polygons(col="resid", title="BART Machine Residuals (log-days)", style="sd", midpoint=NA, palette="RdGy") + tm_layout(inner.margins=0.03, legend.outside = T, bg.color = "white")

#Moran's I test for residual clustering, monte carlo
mc_resid <- moran.mc(shp$resid,lw, 1000)
mc_resid
plot(mc_resid, main="", las=1)

save.image("BART_COVID19_Final.RData")



