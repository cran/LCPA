# library(LCPA)
# set.seed(610893)
#
# ################################################# data conditions ###################################################
# # # "E", "V", "E0", "V0", "EE", "VE", "EV", "VV"
# distribution <- c("random")
# constraint <- list(VV="VV")
# L <- c(10)
# I <- c(18)
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
# device = "CPU"
#
# ################## Mplus ##################
# files.path <- "inst"
# files.clean <- TRUE
#
# trials.length <- length(constraint) * length(distribution) * length(L) * length(I) * length(N)
# trial.list <- list(constraint=constraint, distribution=distribution, L=L, I=I, N=N)
# condition.names <- c("constraint", "distribution", "L", "I", "N")
# methods.names <- c("NNE.Att", "NNE.noAtt")
# results.names <- c("MSE.m", "MSE.c", "acc", "conv", "time")
#
# times <- 200
# times.interval <- 1
# round.n <- 4
# results.pre <- matrix(0, trials.length, length(methods.names)*length(results.names)+1)
# conditions <- matrix(0, trials.length, length(trial.list))
#
# ######################## load or create result file ########################
# names.file <- paste0("res_LPA_", names(constraint)[1], "_", L, "_att")
# res.path <- paste0("results/", names.file, ".rds")
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
#                       paste0("I=", I.cur), paste0("N=", N.cur), 1, "NNE.Att", ]
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
#                                                      "NNE.Att", ], 1, mean, na.rm = TRUE),
#                                        apply(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                                                      paste0("I=", I.cur), paste0("N=", N.cur), ,
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
#   constraint.value.cur <- constraint[[runs.cur[1]]]
#   constraint.cur <- names(constraint)[runs.cur[1]]
#   distribution.cur <- distribution[runs.cur[2]]
#   L.cur <- L[runs.cur[3]]
#   I.cur <- I[runs.cur[4]]
#   N.cur <- N[runs.cur[5]]
#
#   if(constraint.cur == "cus"){
#     if(I.cur == 6){
#       constraint.value.cur <- list(c(1, 1), c(4, 4),
#                                    c(1, 2), c(5, 6))
#     }else if(I.cur == 12){
#       constraint.value.cur <- list(c(1, 1), c(4, 4), c(7, 7), c(10, 10),
#                                    c(1, 2), c(5, 6), c(7, 8), c(11, 12))
#     }else if(I.cur == 18){
#       constraint.value.cur <- list(c(1, 1), c(4, 4), c(7, 7), c(10, 10), c(13, 13), c(16, 16),
#                                    c(1, 2), c(5, 6), c(7, 8), c(11, 12), c(13, 14), c(17, 18))
#     }
#   }
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
#       time.NNE.Att <- system.time({
#         hidden.layers <- c(16, 16)
#         res.NNE.Att <- LPA(response,
#                            L = L.cur, par.ini = par.ini, constraint = constraint.value.cur,
#                            method="NNE", nrep = nrep, starts=starts, maxiter.wa=maxiter.wa, vis = vis,
#                            control.NNE=list(hidden.layers=hidden.layers, activation.function=activation.function, use.attention=TRUE,
#                                             initial.temperature=initial.temperature, cooling.rate=cooling.rate, maxiter.sa=maxiter.sa, threshold.sa=threshold.sa,
#                                             maxiter=maxiter, maxiter.early=maxiter.early, maxcycle=maxcycle,
#                                             lr = lr, scheduler.patience = scheduler.patience, scheduler.factor = scheduler.factor,
#                                             plot.interval=plot.interval,
#                                             device=device))
#
#         res.NNE.Att$params$covs
#         res.NNE.Att$params$means
#
#       })
#
#       time.NNE.noAtt <- system.time({
#         hidden.layers <- c(16, 16)
#         res.NNE.noAtt <- LPA(response,
#                              L = L.cur, par.ini = par.ini, constraint = constraint.value.cur,
#                              method="NNE", nrep = nrep, starts=starts, maxiter.wa=maxiter.wa, vis = vis,
#                              control.NNE=list(hidden.layers=hidden.layers, activation.function=activation.function, use.attention=FALSE,
#                                               initial.temperature=initial.temperature, cooling.rate=cooling.rate, maxiter.sa=maxiter.sa, threshold.sa=threshold.sa,
#                                               maxiter=maxiter, maxiter.early=maxiter.early, maxcycle=maxcycle,
#                                               lr = lr, scheduler.patience = scheduler.patience, scheduler.factor = scheduler.factor,
#                                               plot.interval=plot.interval,
#                                               device=device))
#
#         res.NNE.noAtt$params$covs
#         res.NNE.noAtt$params$means
#
#       })
#
#     })
#     results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#             paste0("I=", I.cur), paste0("N=", N.cur), ,
#             "NNE.Att", t] <- c(LCPA:::get.index.LPA(res.t, res.NNE.Att),
#                                length(res.NNE.Att$Log.Lik.history) < maxiter * maxcycle,
#                                time.NNE.Att[3])
#     results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#             paste0("I=", I.cur), paste0("N=", N.cur), ,
#             "NNE.noAtt", t] <- c(LCPA:::get.index.LPA(res.t, res.NNE.noAtt),
#                                  length(res.NNE.noAtt$Log.Lik.history) < maxiter * maxcycle,
#                                  time.NNE.noAtt[3])
#
#     ####### preview result ########
#     time.posi <- time.cur[3] + time.posi
#     cat("mean of time cost in", paste0(posi, "/", trials.length), ":", round(time.posi / t, round.n), "\n\n")
#     results.cur.NNE.Att <- matrix(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                                           paste0("I=", I.cur), paste0("N=", N.cur), ,
#                                           "NNE.Att", 1:t], nrow = length(results.names))
#     results.cur.NNE.noAtt <- matrix(results[paste0("constraint=", constraint.cur), paste0("distribution=", distribution.cur), paste0("L=", L.cur),
#                                             paste0("I=", I.cur), paste0("N=", N.cur), ,
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
#       saveRDS(results, paste0("results/", names.file, ".rds"))
#     }
#     t <- t + 1
#   }
#
#   posi <- posi + 1
#   saveRDS(results, paste0("results/", names.file, ".rds"))
# }
