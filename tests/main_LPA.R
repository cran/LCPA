# library(LCPA)
# set.seed(610893)
#
# ################################################# data conditions ###################################################
# # # "E", "V", "E0", "V0", "EE", "VE", "EV", "VV"
# distribution <- c("random")
# constraint <- list(
#                     # VV="VV",
#                     # VE="VE",
#                     # EV="EV",
#                     # EE="EE",
#                     # V0="V0",
#                     # E0="E0",
#                     cus=list(c(1, 2), c(3, 3), c(6, 5), c(4, 4), c(3, 4))
#                   )
# L <- c(3, 6) ## 4, 8
# I <- c(6, 12)
# N <- c(200, 1000)
#
# ################### NNE ###################
# activation.function <- "tanh"
# maxiter <- 1000
# maxiter.early <- 100
# maxcycle <- 20
# get.SE <- FALSE
# vis <- TRUE
# initial.temperature <- 1000
# cooling.rate <- 0.5
# maxiter.sa <- 1000
# threshold.sa <- 1e-10
# plot.interval <- 200
# lambda <- 1e-5
# lr = 0.025
# scheduler.patience = 10
# scheduler.factor = 0.80
# par.ini <- "random"
# nrep <- 20
# starts <- 100
# maxiter.wa <- 20
#
# trials.length <- length(constraint) * length(distribution) * length(L) * length(I) * length(N)
# trial.list <- list(constraint=constraint, distribution=distribution, L=L, I=I, N=N)
# condition.names <- c("constraint", "distribution", "L", "I", "N")
# # methods.names <- c("TEM", "REM", "Mplus", "NNE")
# methods.names <- c("TEM", "REM")
# # methods.names <- c("TEM", "REM", "NNE")
# results.names <- c(c("MSE.m", "MSE.c", "acc"), "time")
#
# times <- 20
# times.interval <- 20
# round.n <- 4
# results.pre <- matrix(0, trials.length, length(methods.names)*length(results.names)+1)
# conditions <- matrix(0, trials.length, length(trial.list))
#
# ######################## load or create result file ########################
# names.file <- "res_LPA"
# res.path <- paste0("tests/results/", names.file, ".rds")
# dimnames.this <- list(paste0("constraint=", names(constraint)), paste0("distribution=", distribution), paste0("L=", L),
#                       paste0("I=", I), paste0("N=", N),
#                       results.names, methods.names, paste0("times=", 1:times))
# if(file.exists(res.path)){
#   results <- readRDS(res.path)
#   if(!identical(dimnames(results), dimnames.this) | dim(results)[length(dim(results))] < times){
#     results <- array(dim = c(length(constraint), length(distribution), length(L), length(I),
#                              length(N), length(results.names), length(methods.names), times),
#                      dimnames = dimnames.this)
#   }
# }else{
#   results <- array(dim = c(length(constraint), length(distribution), length(L), length(I),
#                            length(N), length(results.names), length(methods.names), times),
#                    dimnames = dimnames.this)
# }
#
# ##################################################### pre-load #######################################################
# posi.start <- 1
# t.start <- 1
# while(posi.start <= trials.length){
#
#   runs.cur <- LCPA:::get.runs(posi.start, trial.list)
#   constraint.value.cur <- constraint[runs.cur[1]]
#   constraint.cur <- names(constraint)[runs.cur[1]]
#   distribution.cur <- distribution[runs.cur[2]]
#   L.cur <- L[runs.cur[3]]
#   I.cur <- I[runs.cur[4]]
#   N.cur <- N[runs.cur[5]]
#
#   time.sum <- results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                       paste0("I=", I.cur), paste0("N=", N.cur), 1, "TEM", ]
#   if(any(is.na(time.sum))){
#     t.start <- min(which(is.na(time.sum)))
#     break
#   }else if(any(is.nan(time.sum))){
#     t.start <- min(which(is.nan(time.sum)))
#     break
#   }
#
#   conditions[posi.start, ] <- c(constraint.cur, distribution.cur, L.cur, I.cur, N.cur)
#   results.pre[posi.start, ] <- round(c(apply(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                                                      paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                                      "TEM", ], 1, mean, na.rm = TRUE),
#                                        apply(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                                                      paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                                      "REM", ], 1, mean, na.rm = TRUE),
#                                        # apply(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                                        #               paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                        #               "Mplus", ], 1, mean, na.rm = TRUE),
#                                        # apply(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                                        #               paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                        #               "NNE", ], 1, mean, na.rm = TRUE),
#                                        0), round.n)
#   posi.start <- posi.start + 1
# }
#
# cat("finish load...\n")
# posi.print <- c(max(posi.start-19, 1):min(max(posi.start, 2), trials.length))
# if(all(posi.print == 1)){
#   res_preview <- data.frame(cbind(t(conditions[posi.print, ]), t(results.pre[posi.print, ])))
# }else{
#   res_preview <- data.frame(conditions[posi.print, ], results.pre[posi.print, ])
# }
# rownames(res_preview) <- posi.print
# colnames(res_preview) <- c(condition.names,
#                            as.vector(sapply(methods.names, function(x){
#                              paste0(x, "-", results.names)
#                            })),
#                            "Mean_Time")
# print(res_preview)
#
# ################################################# continue running ###################################################
# posi <- posi.start
# while(posi <= trials.length){
#
#   runs.cur <- LCPA:::get.runs(posi, trial.list)
#   constraint.value.cur <- constraint[[runs.cur[1]]]
#   constraint.cur <- names(constraint)[runs.cur[1]]
#   distribution.cur <- distribution[runs.cur[2]]
#   L.cur <- L[runs.cur[3]]
#   I.cur <- I[runs.cur[4]]
#   N.cur <- N[runs.cur[5]]
#
#   conditions[posi, ] <- c(constraint.cur, distribution.cur, L.cur, I.cur, N.cur)
#   time.posi <- 0
#   if(posi == posi.start){
#     t <- t.start
#   }else{
#     t <- 1
#   }
#
#   while(t <= times){
#     cat("===============================", paste0(posi, "/", trials.length), ":",
#         paste0(t, "/", times), "===============================\n")
#
#     res.t <- data <- sim.LPA(N=N.cur, I=I.cur, L=L.cur, constraint = constraint.value.cur,
#                              distribution=distribution.cur)
#     res.t$covs
#
#     response <- data$response
#
#     time.cur <- system.time({
#
#       time.TEM <- system.time({
#         res.TEM <- LPA(response,
#                        L = L.cur, par.ini = res.t, constraint = constraint.value.cur,
#                        method="EM", nrep = 1, vis = vis,
#                        control.EM=list(maxiter=2000, tol=1e-4))
#       })
#
#       time.REM <- system.time({
#
#         res.REM <- LPA(response,
#                        L = L.cur, par.ini = par.ini, constraint = constraint.value.cur,
#                        method="EM", nrep = nrep, vis = vis,
#                        control.EM=list(maxiter=2000, starts=starts, maxiter.wa=maxiter.wa, tol=1e-4))
#
#       })
#
#       # time.Mplus<- system.time({
#       #   tryCatch({
#       #     res.Mplus <- LPA(response,
#       #                      L = L.cur, par.ini = par.ini, constraint = constraint.value.cur,
#       #                      method="Mplus", nrep = nrep, vis = vis,
#       #                      control.Mplus=list(maxiter=2000, starts=starts, maxiter.wa=maxiter.wa, tol=1e-4, clean.files=TRUE))
#       #   }, error = function(e){
#       #     res.Mplus <- res.REM
#       #   })
#       # })
#       #
#       # time.NNE <- system.time({
#       #   hidden.layers <- c(ceiling(log10(N.cur)) * ceiling(log2(L.cur)),
#       #                      ceiling(log10(N.cur)) * ceiling(log2(L.cur)))
#       #   res.NNE <- LPA(response,
#       #                  L = L.cur, par.ini = par.ini, constraint = constraint.value.cur,
#       #                  method="NNE", nrep = nrep, vis = vis,
#       #                  control.NNE=list(hidden.layers=hidden.layers, activation.function=activation.function,
#       #                                   initial.temperature=initial.temperature, cooling.rate=cooling.rate, maxiter.sa=maxiter.sa, threshold.sa=threshold.sa,
#       #                                   maxiter=maxiter, maxiter.early=maxiter.early, maxcycle=maxcycle,
#       #                                   starts=starts, maxiter.wa=maxiter.wa,
#       #                                   lr = lr, scheduler.patience = scheduler.patience, scheduler.factor = scheduler.factor,
#       #                                   plot.interval=plot.interval))
#       #
#       #   res.NNE$params$covs
#       #   res.NNE$params$means
#       #
#       # })
#     })
#     results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#             paste0("I=", I.cur), paste0("N=", N.cur), ,
#             "TEM", t] <- c(LCPA:::get.index.LPA(res.t, res.TEM), time.TEM[3])
#     results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#             paste0("I=", I.cur), paste0("N=", N.cur), ,
#             "REM", t] <- c(LCPA:::get.index.LPA(res.t, res.REM), time.REM[3])
#     # results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#     #         paste0("I=", I.cur), paste0("N=", N.cur), ,
#     #         "Mplus", t] <- c(LCPA:::get.index.LPA(res.t, res.Mplus), time.Mplus[3])
#     # results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#     #         paste0("I=", I.cur), paste0("N=", N.cur), ,
#     #         "NNE", t] <- c(LCPA:::get.index.LPA(res.t, res.NNE), time.NNE[3])
#
#     ####### preview result ########
#     time.posi <- time.cur[3] + time.posi
#     cat("mean of time cost in", paste0(posi, "/", trials.length), ":", round(time.posi / t, round.n), "\n\n")
#     results.cur.TEM <- matrix(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                                      paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                      "TEM", 1:t], nrow = length(results.names))
#     results.cur.REM <- matrix(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                                      paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                      "REM", 1:t], nrow = length(results.names))
#     # results.cur.Mplus <- matrix(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#     #                                  paste0("I=", I.cur), paste0("N=", N.cur), ,
#     #                                  "Mplus", 1:t], nrow = length(results.names))
#     # results.cur.NNE <- matrix(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#     #                                   paste0("I=", I.cur), paste0("N=", N.cur), ,
#     #                                   "NNE", 1:t], nrow = length(results.names))
#     results.pre[posi, 1:(length(methods.names)*length(results.names)+1)] <- round(c(apply(results.cur.TEM, 1, mean, na.rm = TRUE),
#                                                                                     apply(results.cur.REM, 1, mean, na.rm = TRUE),
#                                                                                     # apply(results.cur.Mplus, 1, mean, na.rm = TRUE),
#                                                                                     # apply(results.cur.NNE, 1, mean, na.rm = TRUE),
#                                                                                     time.posi/t), round.n)
#     posi.print <- c(max(posi-19, 1):min(max(posi, 2), trials.length))
#     if(all(posi.print == 1)){
#       res_preview <- data.frame(cbind(t(conditions[posi.print, ]), t(results.pre[posi.print, ])))
#     }else{
#       res_preview <- data.frame(conditions[posi.print, ], results.pre[posi.print, ])
#     }
#     rownames(res_preview) <- posi.print
#     colnames(res_preview) <- c(condition.names,
#                                as.vector(sapply(methods.names, function(x){
#                                  paste0(x, ".", results.names)
#                                })),
#                                "Mean_Time")
#     cat("------------------- current accurace pre-view -------------------\n")
#     print(res_preview)
#     cat("\n")
#     cat("=============================================================================\n")
#     cat("\n")
#
#     if(t %% times.interval == 0){
#       saveRDS(results, paste0("tests/results/", names.file, ".rds"))
#     }
#     t <- t + 1
#   }
#
#   posi <- posi + 1
#   saveRDS(results, paste0("tests/results/", names.file, ".rds"))
# }
