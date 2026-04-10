#
# library(LCPA)
#
# 235455
#
# set.seed(56756765)
#
# ################### NNE ###################
# hidden.layers <- c(16, 16)
# activation.function <- "tanh"
# d.model=8
# nhead=2
# dim.feedforward=16
# eps=1e-8
# maxiter <- 1000
# maxiter.early <- 100
# maxcycle <- 20
# initial.temperature <- 1000
# cooling.rate <- 0.5
# maxiter.sa <- 1000
# threshold.sa <- 1e-10
# plot.interval <- 500
# lambda <- 1e-5
# lr = 0.025
# scheduler.patience = 10
# scheduler.factor = 0.80
# par.ini <- "random"
# nrep <- 20
# starts <- 100
# maxiter.wa <- 20
# device <- "CPU"
# vis <- TRUE
#
# ################## Mplus ##################
# files.path <- "inst"
# files.clean <- TRUE
#
# ################################# indicators for LPA ###################################
# library(tidyverse)
# library(openxlsx)
# #### fMRI ####
# data.original <- read.xlsx("W:\\_FilesShare\\LCPA_ Research_and_LCPA_Package\\code\\LCPA research (OSF)\\results\\empirical_studies\\schaefer17_CET.xlsx")
# response <- as.matrix(scale(data.original[, -1]))
#
# ############################################### LPA ###############################################
# set.seed(56756765)
# L.max <- 10
# results <- vector("list", L.max)
# for (i in 1:L.max) {
#   results[[i]] <- list()
# }
# for(l in 1:L.max){
#
#   # Frequently read and save analysis results,
#   # and check whether a certain type of analysis has already been done.
#
#   # This approach ensures that even if the program is interrupted,
#   # it can resume from where it left off when run again,
#   # thereby saving computational resources/costs.
#
#   if(file.exists(paste0("results/empirical_studies/results_LPA.rds"))){
#     results <- readRDS(paste0("results/empirical_studies/results_LPA.rds"))
#     if(!is.null(results[[l]]$res.Mplus) && !is.null(results[[l]]$res.NNE)){
#       next
#     }
#   }
#
#   cat("=========================== Starting: L =", l, " ===========================\n")
#   # Run NNE
#   cat("    Running NNE...\n")
#   time.NNE <- system.time({
#     res.NNE <- LPA(response,
#                    L = l, par.ini = par.ini,
#                    method="NNE", nrep = nrep, starts=starts, maxiter.wa=maxiter.wa, vis = vis,
#                    control.NNE=list(hidden.layers=hidden.layers, activation.function=activation.function,
#                                     d.model=d.model, nhead=nhead, dim.feedforward=dim.feedforward, eps=eps,
#                                     initial.temperature=initial.temperature, cooling.rate=cooling.rate, maxiter.sa=maxiter.sa, threshold.sa=threshold.sa,
#                                     maxiter=maxiter, maxiter.early=maxiter.early, maxcycle=maxcycle,
#                                     lr = lr, scheduler.patience = scheduler.patience, scheduler.factor = scheduler.factor,
#                                     plot.interval=plot.interval,
#                                     device=device))
#   })
#   cat("    NNE done. Elapsed:", round(time.NNE["elapsed"], 2), "sec\n")
#
#   # Run Mplus
#   cat("    Running Mplus...\n")
#   time.Mplus <- system.time({
#     res.Mplus <- tryCatch({
#       temp <- LPA(response,
#                   L=l, par.ini=par.ini, method="Mplus", nrep=nrep, starts=starts, maxiter.wa=maxiter.wa, vis=vis,
#                   control.Mplus=list(maxiter=2000, tol=1e-4, files.path=paste0(files.path, paste0("/L=", l)), files.clean=files.clean))
#
#       if (is.null(temp) ||
#           (!is.list(temp)) ||
#           (is.list(temp) && length(temp) == 0) ||
#           (is.list(temp) && !is.null(temp$error))) {
#         stop("Mplus run failed silently")
#       }
#       temp
#     },
#     error = function(e) {
#       warning("Mplus failed (", conditionMessage(e), "); falling back to NNE")
#       res.Mplus <- res.NNE
#     })
#   })
#
#   cat("    Mplus done. Elapsed:", round(time.Mplus["elapsed"], 2), "sec\n")
#
#   print(compare.model(res.NNE, res.Mplus))
#
#   results[[l]]$res.Mplus <- res.Mplus
#   results[[l]]$res.NNE <- res.NNE
#   results[[l]]$time.NNE <- time.NNE["elapsed"]
#   results[[l]]$time.Mplus <- time.Mplus["elapsed"]
#
#   saveRDS(results, paste0("results/empirical_studies/results_LPA.rds"))
# }
#
# ####################################################### Bootstrap LR test #######################################################
# # Flexible read/write mechanism that allows:
# # - NNE bootstrap likelihood ratio test
# # - Mplus bootstrap likelihood ratio test
# # to run independently of each other in the same file,
# # while saving their results promptly into the same .rds file
# # (results_LPA.rds)
#
# ################## Bootstrap LR test for NNE ##################
# library(LCPA)
# set.seed(56756765)
# L.max <- 10
# n.Bootstrap <- 100
# for(l in 2:L.max){
#   cat("=========================== Starting: L =", l, " ===========================\n")
#   results.NNE <- readRDS(paste0("results/empirical_studies/results_LPA.rds"))
#
#   if(is.null(results.NNE[[l]]$bootstrap.obj.NNE)){
#     results.NNE[[l]]$bootstrap.obj.NNE <- LRT.test.Bootstrap(results.NNE[[l-1]]$res.NNE,
#                                                              results.NNE[[l]]$res.NNE,
#                                                              n.Bootstrap = n.Bootstrap,
#                                                              vis = TRUE)
#
#     temp <- readRDS(paste0("results/empirical_studies/results_LPA.rds"))
#     temp[[l]]$bootstrap.obj.NNE <- results.NNE[[l]]$bootstrap.obj.NNE
#     saveRDS(temp, paste0("results/empirical_studies/results_LPA.rds"))
#   }
# }
#
# ################## Bootstrap LR test for Mplus ##################
# library(LCPA)
# set.seed(56756765)
# L.max <- 10
# n.Bootstrap <- 100
# for(l in 2:L.max){
#   cat("=========================== Starting: L =", l, " ===========================\n")
#   results.Mplus <- readRDS(paste0("results/empirical_studies/results_LPA.rds"))
#
#   if(is.null(results.Mplus[[l]]$bootstrap.obj.Mplus)){
#     results.Mplus[[l]]$bootstrap.obj.Mplus <- LRT.test.Bootstrap(results.Mplus[[l-1]]$res.Mplus,
#                                                                  results.Mplus[[l]]$res.Mplus,
#                                                                  n.Bootstrap = n.Bootstrap,
#                                                                  vis = TRUE)
#
#     temp <- readRDS(paste0("results/empirical_studies/results_LPA.rds"))
#     temp[[l]]$bootstrap.obj.Mplus <- results.Mplus[[l]]$bootstrap.obj.Mplus
#     saveRDS(temp, paste0("results/empirical_studies/results_LPA.rds"))
#   }
# }
#
# ######################################## compare model ############################################
# library(LCPA)
# library(LCPA)
# set.seed(56756765)
# L.max <- 10
# n.Bootstrap <- 100
# for(l in 2:L.max){
#   cat("=========================== Starting: L =", l, " ===========================\n")
#   results.LPA <- readRDS(paste0("results/empirical_studies/results_LPA.rds"))
#
#   results.LPA[[l]]$compare.obj.Mplus <- compare.model(results.LPA[[l-1]]$res.Mplus,
#                                                       results.LPA[[l]]$res.Mplus,
#                                                       n.Bootstrap = 0)
#   results.LPA[[l]]$compare.obj.Mplus$LRT.Bootstrap.obj <- results.LPA[[l]]$bootstrap.obj.Mplus
#
#   results.LPA[[l]]$compare.obj.NNE <- compare.model(results.LPA[[l-1]]$res.NNE,
#                                                     results.LPA[[l]]$res.NNE,
#                                                     n.Bootstrap = 0)
#   results.LPA[[l]]$compare.obj.NNE$LRT.Bootstrap.obj <- results.LPA[[l]]$bootstrap.obj.NNE
#
#   temp <- readRDS(paste0("results/empirical_studies/results_LPA.rds"))
#   temp[[l]]$compare.obj.Mplus <- results.LPA[[l]]$compare.obj.Mplus
#   temp[[l]]$compare.obj.NNE <- results.LPA[[l]]$compare.obj.NNE
#   saveRDS(temp, paste0("results/empirical_studies/results_LPA.rds"))
# }
#
# ######################################## analysis ############################################
# results.LPA <- readRDS(paste0("results/empirical_studies/results_LPA.rds"))
# L.max <- 10
# fit.index.NNE <- fit.index.Mplus <- NULL
# for(l in 1:L.max){
#   ##### NNE
#   fit.index.obj.NNE <- get.fit.index(results.LPA[[l]]$res.NNE)
#   fit.index.NNE.cur <- c(fit.index.obj.NNE$npar, fit.index.obj.NNE$Log.Lik, fit.index.obj.NNE$AIC,
#                          fit.index.obj.NNE$BIC, fit.index.obj.NNE$SIC, fit.index.obj.NNE$CAIC,
#                          fit.index.obj.NNE$AWE, fit.index.obj.NNE$SABIC)
#   fit.index.NNE.cur <- c(get.AvePP(results.LPA[[l]]$res.NNE)[l+1, l+1],
#                          get.entropy(results.LPA[[l]]$res.NNE), fit.index.NNE.cur)
#   if(l > 1){
#     fit.index.NNE.cur <- c(fit.index.NNE.cur,
#                            c(results.LPA[[l]]$compare.obj.NNE$BF,
#                              results.LPA[[l]]$compare.obj.NNE$LRT.obj$statistic,
#                              results.LPA[[l]]$compare.obj.NNE$LRT.obj$parameter,
#                              results.LPA[[l]]$compare.obj.NNE$LRT.obj$p.value,
#                              results.LPA[[l]]$compare.obj.NNE$LRT.VLMR.obj$p.value,
#                              results.LPA[[l]]$compare.obj.NNE$LRT.Bootstrap.obj$p.value))
#   }else{
#     fit.index.NNE.cur <- c(fit.index.NNE.cur, c(NA, NA, NA, NA, NA, NA))
#   }
#   fit.index.NNE <- rbind(fit.index.NNE, fit.index.NNE.cur)
#
#   ##### Mplus
#   fit.index.obj.Mplus <- get.fit.index(results.LPA[[l]]$res.Mplus)
#   fit.index.Mplus.cur <- c(fit.index.obj.Mplus$npar, fit.index.obj.Mplus$Log.Lik, fit.index.obj.Mplus$AIC,
#                            fit.index.obj.Mplus$BIC, fit.index.obj.Mplus$SIC, fit.index.obj.Mplus$CAIC,
#                            fit.index.obj.Mplus$AWE, fit.index.obj.Mplus$SABIC)
#   fit.index.Mplus.cur <- c(get.AvePP(results.LPA[[l]]$res.Mplus)[l+1, l+1],
#                            get.entropy(results.LPA[[l]]$res.Mplus), fit.index.Mplus.cur)
#   if(l > 1){
#     fit.index.Mplus.cur <- c(fit.index.Mplus.cur,
#                              c(results.LPA[[l]]$compare.obj.Mplus$BF,
#                                results.LPA[[l]]$compare.obj.Mplus$LRT.obj$statistic,
#                                results.LPA[[l]]$compare.obj.Mplus$LRT.obj$parameter,
#                                results.LPA[[l]]$compare.obj.Mplus$LRT.obj$p.value,
#                                results.LPA[[l]]$compare.obj.Mplus$LRT.VLMR.obj$p.value,
#                                results.LPA[[l]]$compare.obj.Mplus$LRT.Bootstrap.obj$p.value))
#   }else{
#     fit.index.Mplus.cur <- c(fit.index.Mplus.cur, c(NA, NA, NA, NA, NA, NA))
#   }
#   fit.index.Mplus <- rbind(fit.index.Mplus, fit.index.Mplus.cur)
# }
#
# rownames(fit.index.NNE) <- rownames(fit.index.Mplus) <- paste0("class", 1:L.max)
# colnames(fit.index.NNE) <- colnames(fit.index.Mplus) <- c("AvePP", "entropy", "npar", "Log.Lik", "AIC",
#                                                           "BIC", "SIC", "CAIC", "AWE", "SABIC", "BF",
#                                                           "LR statistic", "df", "p.value", "VLMR.p.value",
#                                                           "bootstrap.p.value")
# ########################################## excel ##########################################
# library(openxlsx)
# wb <- createWorkbook()
# addWorksheet(wb, "NNE")
# writeData(wb, sheet = "NNE", as.table(fit.index.NNE))
# addWorksheet(wb, "Mplus")
# writeData(wb, sheet = "Mplus", as.table(fit.index.Mplus))
# saveWorkbook(wb, file = "results/empirical_studies/LPA_fit.xlsx", overwrite = TRUE)
#
# ####################################### plot #######################################
# l <- 3
# NNE.obj <- results.LPA[[l]]$res.NNE
# Mplus.obj <- adjust.model(object1=NNE.obj, object2=results.LPA[[l]]$res.Mplus)
#
# mean(NNE.obj$Z == Mplus.obj$Z)
#
# plot(NNE.obj, ncol=3)
# plot(Mplus.obj, ncol=3)
#
# # compare.model(NNE.obj, Mplus.obj)
#

