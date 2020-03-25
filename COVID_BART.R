options(java.parameters = "-Xmx30g")

library(foreign)
library(bartMachine)
set_bart_machine_num_cores(5)

setwd("D:/SARS-CoV-2/MasterFile")

d <- read.dbf("RKI_with_All.dbf")

y <- d$cases
Xa <- d[c(47:346)]
X <- Xa[-c(319,)]

bm <- bartMachine(X = X, y = y)
#k_fold <- k_fold_cv(X=X, y=y, k_folds=5)
#bm_cv <- bartMachineCV(X=X, y=y)
#assump <- check_bart_error_assumptions(bm)
#plot_convergence_diagnostics(bm_cv)
#vi <- investigate_var_importance(bm, num_replicates_for_avg = 10)
varselect <- var_selection_by_permute(bm, plot = TRUE)
summary(bm)
