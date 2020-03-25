options(java.parameters = "-Xmx30g")

library(foreign)
library(bartMachine)
set_bart_machine_num_cores(5)

setwd("D:/SARS-CoV-2/MasterFile")

#Read data
d <- read.csv("Export_output_2.txt")

#Seperate dependent and independent variables
y <- d$cases
Xa <- d[c(48:347)]
#Remove unwanted variables
X <- subset(Xa, select=-c(S_011, S_111, S_112, S_115, S_LK, S_Ag, S_068, S_077, S_082, S_085, S_086, S_105, S_107, S_169))
#Check for factors, due to decimal-comma conversion
names(Filter(is.factor, X))

bm <- bartMachine(X = X, y = y)
summary(bm)
#Model will exhibit serious overfitting! Proceed to assess with k-folds CV.

#k-folds CV
k_fold <- k_fold_cv(X=X, y=y, k_folds=5)

windows()
varselect <- var_selection_by_permute(bm, plot = T, num_reps_for_avg = 5, num_permute_samples = 10, num_trees_for_permute = 10)
print(varselect$important_vars_local_names)
print(varselect$important_vars_global_max_names)
print(varselect$important_vars_global_se_names) 

X2 <- subset(Xa, select=c(tot_P, tot_R, S_109, Rskh_k2, Rskh_pp, RmSh_k2, S_078, S_018, S_030, S_176, Rjw_km2, RmSh_pp, Rhn_km2, S_113, S_075, RchE_k2, S_044, Rta_km2, sch_km2, ff_km2, Rtao_pp, cc_km2, S_159, Rch_pop, S_095, S_124))
#BARTmachine with selected variables, built using CV parameterisation

bm_cv <- bartMachineCV(X=X2, y=y)
summary(bm_cv)
assump <- check_bart_error_assumptions(bm_cv)
plot_convergence_diagnostics(bm_cv)
vi <- investigate_var_importance(bm_cv, num_replicates_for_avg = 10)
