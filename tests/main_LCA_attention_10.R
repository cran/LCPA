# library(here)
# setwd(here())
#
# library(LCPA)
# set.seed(1237894)
#
# ################################################# data conditions ###################################################
# distribution <- c("random")
# L <- c(10)
# I <- c(100)
# IQ <- c(0.8)
# length.poly <- c(5)
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
# device <- "CPU"
#
# files.path <- "inst"
# files.clean <- TRUE
#
# trials.length <- length(distribution) * length(L) * length(IQ) * length(length.poly) * length(I) * length(N)
# trial.list <- list(distribution=distribution, L=L, IQ=IQ, poly=length.poly, I=I, N=N)
# condition.names <- c("distribution", "L", "IQ", "poly", "I", "N")
# methods.names <- c("NNE.Att", "NNE.noAtt")
# results.names <- c("MSE", "acc", "conv", "time")
#
# times <- 200
# times.interval <- 1
# round.n <- 4
# results.pre <- matrix(0, trials.length, length(methods.names)*length(results.names)+1)
# conditions <- matrix(0, trials.length, length(trial.list))
#
# ######################## load or create result file ########################
# names.file <- paste0("res_LCA_", L, "_", IQ, "_att")
# res.path <- paste0("tests/results/", names.file, ".rds")
# dimnames.this <- list(paste0("distribution=", distribution), paste0("L=", L), paste0("IQ=", IQ),
#                       paste0("poly=", length.poly), paste0("I=", I), paste0("N=", N),
#                       results.names, methods.names, paste0("times=", 1:times))
# if(file.exists(res.path)){
#   results <- readRDS(res.path)
#   if(!identical(dimnames(results), dimnames.this) | dim(results)[length(dim(results))] < times){
#     results <- array(dim = c(length(distribution), length(L), length(IQ), length(length.poly), length(I),
#                              length(N), length(results.names), length(methods.names), times),
#                      dimnames = dimnames.this)
#   }
# }else{
#   results <- array(dim = c(length(distribution), length(L), length(IQ), length(length.poly), length(I),
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
#   distribution.cur <- distribution[runs.cur[1]]
#   L.cur <- L[runs.cur[2]]
#   IQ.cur <- IQ[runs.cur[3]]
#   length.poly.cur <- length.poly[runs.cur[4]]
#   I.cur <- I[runs.cur[5]]
#   N.cur <- N[runs.cur[6]]
#
#   time.sum <- results[paste0("distribution=", distribution.cur), paste0("L=", L.cur), paste0("IQ=", IQ.cur), paste0("poly=", length.poly.cur),
#                       paste0("I=", I.cur), paste0("N=", N.cur), 1, "NNE.Att", ]
#   if(any(is.na(time.sum))){
#     t.start <- min(which(is.na(time.sum)))
#     break
#   }else if(any(is.nan(time.sum))){
#     t.start <- min(which(is.nan(time.sum)))
#     break
#   }
#
#   conditions[posi.start, ] <- c(distribution.cur, L.cur, IQ.cur, length.poly.cur, I.cur, N.cur)
#   results.pre[posi.start, ] <- round(c(apply(results[paste0("distribution=", distribution.cur), paste0("L=", L.cur), paste0("IQ=", IQ.cur),
#                                                      paste0("poly=", length.poly.cur), paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                                      "NNE.Att", ], 1, mean, na.rm = TRUE),
#                                        apply(results[paste0("distribution=", distribution.cur), paste0("L=", L.cur), paste0("IQ=", IQ.cur),
#                                                      paste0("poly=", length.poly.cur), paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                                      "NNE.noAtt", ], 1, mean, na.rm = TRUE),
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
#   distribution.cur <- distribution[runs.cur[1]]
#   L.cur <- L[runs.cur[2]]
#   IQ.cur <- IQ[runs.cur[3]]
#   length.poly.cur <- length.poly[runs.cur[4]]
#   I.cur <- I[runs.cur[5]]
#   N.cur <- N[runs.cur[6]]
#
#   conditions[posi, ] <- c(distribution.cur, L.cur, IQ.cur, length.poly.cur, I.cur, N.cur)
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
#     res.t <- NULL
#     data <- sim.LCA(N=N.cur, I=I.cur, L=L.cur, poly.value=length.poly.cur, IQ=IQ.cur, distribution=distribution.cur)
#     response <- data$response
#     res.t$P.Z.Xn <- data$P.Z.Xn
#     res.t$P.Z <- data$P.Z
#     res.t$par <- data$par
#     res.t$Z <- data$Z
#     res.t$poly.value <- data$poly.value
#
#     time.cur <- system.time({
#
#       time.NNE.Att <- system.time({
#         hidden.layers <- c(16,  16)
#         res.NNE.Att <- LCA(response,
#                            L = L.cur, par.ini = par.ini,
#                            method="NNE", nrep = nrep, starts=starts, maxiter.wa=maxiter.wa, vis = vis,
#                            control.NNE=list(hidden.layers=hidden.layers, activation.function=activation.function,
#                                             use.attention = TRUE,
#                                             d.model=8, nhead=2, dim.feedforward=16, eps=1e-8,
#                                             initial.temperature=initial.temperature, cooling.rate=cooling.rate, maxiter.sa=maxiter.sa, threshold.sa=threshold.sa,
#                                             maxiter=maxiter, maxiter.early=maxiter.early, maxcycle=maxcycle,
#                                             lr = lr, scheduler.patience = scheduler.patience, scheduler.factor = scheduler.factor,
#                                             plot.interval=plot.interval,
#                                             device=device))
#       })
#       res.NNE.Att$par <- res.NNE.Att$params$par
#       res.NNE.Att$P.Z <- res.NNE.Att$params$P.Z
#
#       time.NNE.noAtt <- system.time({
#         hidden.layers <- c(16,  16)
#         res.NNE.noAtt <- LCA(response,
#                              L = L.cur, par.ini = par.ini,
#                              method="NNE", nrep = nrep, starts=starts, maxiter.wa=maxiter.wa, vis = vis,
#                              control.NNE=list(hidden.layers=hidden.layers, activation.function=activation.function,
#                                               use.attention = FALSE,
#                                               d.model=8, nhead=2, dim.feedforward=16, eps=1e-8,
#                                               initial.temperature=initial.temperature, cooling.rate=cooling.rate, maxiter.sa=maxiter.sa, threshold.sa=threshold.sa,
#                                               maxiter=maxiter, maxiter.early=maxiter.early, maxcycle=maxcycle,
#                                               lr = lr, scheduler.patience = scheduler.patience, scheduler.factor = scheduler.factor,
#                                               plot.interval=plot.interval,
#                                               device=device))
#       })
#       res.NNE.noAtt$par <- res.NNE.noAtt$params$par
#       res.NNE.noAtt$P.Z <- res.NNE.noAtt$params$P.Z
#
#     })
#     results[paste0("distribution=", distribution.cur), paste0("L=", L.cur), paste0("IQ=", IQ.cur),
#             paste0("poly=", length.poly.cur), paste0("I=", I.cur), paste0("N=", N.cur), ,
#             "NNE.Att", t] <- c(LCPA:::get.index.LCA(res.t, res.NNE.Att),
#                                length(res.NNE.Att$Log.Lik.history) < maxiter * maxcycle, time.NNE.Att[3])
#     results[paste0("distribution=", distribution.cur), paste0("L=", L.cur), paste0("IQ=", IQ.cur),
#             paste0("poly=", length.poly.cur), paste0("I=", I.cur), paste0("N=", N.cur), ,
#             "NNE.noAtt", t] <- c(LCPA:::get.index.LCA(res.t, res.NNE.noAtt),
#                                  length(res.NNE.noAtt$Log.Lik.history) < maxiter * maxcycle, time.NNE.noAtt[3])
#
#     ####### preview result ########
#     time.posi <- time.cur[3] + time.posi
#     cat("mean of time cost in", paste0(posi, "/", trials.length), ":", round(time.posi / t, round.n), "\n\n")
#     results.cur.NNE.Att <- matrix(results[paste0("distribution=", distribution.cur), paste0("L=", L.cur), paste0("IQ=", IQ.cur),
#                                           paste0("poly=", length.poly.cur), paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                           "NNE.Att", 1:t], nrow = length(results.names))
#     results.cur.NNE.noAtt <- matrix(results[paste0("distribution=", distribution.cur), paste0("L=", L.cur), paste0("IQ=", IQ.cur),
#                                             paste0("poly=", length.poly.cur), paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                             "NNE.noAtt", 1:t], nrow = length(results.names))
#     results.pre[posi, 1:(length(methods.names)*length(results.names)+1)] <- round(c(apply(results.cur.NNE.Att, 1, mean, na.rm = TRUE),
#                                                                                     apply(results.cur.NNE.noAtt, 1, mean, na.rm = TRUE),
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
#       saveRDS(results, res.path)
#     }
#     t <- t + 1
#   }
#
#   posi <- posi + 1
#   saveRDS(results, res.path)
# }
