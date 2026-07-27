# 15/10/2018: This file contains the modified rERR functions, to be able to copmute correctly the lags
library(stats4)
library(dplyr)

# If the ID variable contains any characters, update line 53 to use lapply(as.numeric) excluding the ID column.

dyn.load("Y:/llhood_linear_v3.dll") # location of "llhood_linear_v3.dll"

#///////////////////////////////////////////////////////////////////////////////////////////////
  
f_fit_linERR <- function (formula, data, rsets, n_lin_vars, n_loglin_vars, id_name, 
                             time_name) 
{
  # canvis per Jordi F:
  #	- a linia 76 he afegit la condicio "& data$n_pe != -1" per no contar la primera linia ficticia de cada subjecte com a tac
  
  formula_sv <- formula[[2]]
  v_id <- eval(parse(text = paste0("data$", id_name)))
  v_n_pe <- eval(parse(text = paste0("data$", "n_pe")))
  v_entry <- eval(parse(text = paste0("data$", formula_sv[[2]])))
  v_exit <- eval(parse(text = paste0("data$", formula_sv[[3]])))
  v_outcome <- eval(parse(text = paste0("data$", formula_sv[[4]])))
  v_time <- eval(parse(text = paste0("data$", time_name)))
  nrow_cases <- which(v_outcome == 1)
  failtimes <- v_exit[nrow_cases]
  id_cases <- v_id[nrow_cases]
  nrow_cases <- names(rsets)
  beta <- rep(0.1, n_lin_vars + n_loglin_vars)
  beta <- as.list(beta)
  names(beta) <- paste0("x", 1:length(beta))
  add.arguments <- function(f, n) {
    t = paste("arg <- alist(", paste(sapply(1:n, function(i) paste("x", 
                                                                   i, "=", sep = "")), collapse = ","), ")", sep = "")
    formals(f) <- eval(parse(text = t))
    f
  }
  p.est <- function() {
    beta3 <- vector()
    for (i in 1:length(beta)) {
      beta3[i] <- eval(parse(text = paste0("x", i)))
    }
    beta_2 <- as.list(beta3)
    beta_2 <- unlist(lapply(beta_2, as.numeric))
    constr_ind <- as.integer(0)
    res <- .Call("llhood_linear_v3", beta_2, length(beta_2), 
                 length(rsets), rsets, nrseti, data, nrow_cases, n_lin_vars, 
                 n_loglin_vars, constr_ind)
    res <- res[-length(res)]
    res1 <- log(res)
    return(-sum(log(res)))
  }
  p.est <- add.arguments(p.est, length(beta))
  data <- lapply(data, as.numeric)
  rsets <- lapply(rsets, as.numeric)
  nrseti <- lapply(rsets, length)
  nrseti <- lapply(nrseti, as.numeric)
  n_lin_vars <- as.integer(n_lin_vars)
  n_loglin_vars <- as.integer(n_loglin_vars)
  nrow_cases <- lapply(nrow_cases, as.numeric)
  constr <- vector()
  if (n_lin_vars > 0) 
    for (i in 0:(n_lin_vars - 1)) {
      beta_2 <- unlist(lapply(beta, as.numeric))
      res <- .Call("llhood_linear_v3", beta_2, length(beta_2), 
                   length(rsets), rsets, nrseti, data, nrow_cases, 
                   n_lin_vars, n_loglin_vars, as.integer(i))
      constr[i + 1] <- res[length(res)]
    }
  llim <- c(constr, rep(-Inf, length(beta) - n_lin_vars))
  llim2 <- c(constr, rep(0.1, length(beta) - n_lin_vars))
  reduction <- 0.999
  p.est_num <- function(x) {
    x <- as.list(x)
    names(x) <- paste0("x", 1:length(x))
    return(do.call(p.est, x))
  }
  suppressWarnings(while (any(is.na(diag(numDeriv::hessian(p.est_num, 
                                                           llim2 * reduction))))) reduction <- reduction^2)
  res <- mle(minuslogl = p.est, start = beta, method = "L-BFGS-B", 
             lower = llim * reduction)
  names <- names(data)[8:(8 + n_lin_vars + n_loglin_vars - 
                            1)]
  names(attr(res, "coef")) <- names
  names(attr(res, "fullcoef")) <- names
  attr(res, "desc") <- list(desc = list(n_events = length(rsets), 
                                        n_observations = length(data$n_pe[which(data$n_pe != 
                                                                                  0 & data$n_pe != -1)]), n_subjects = length(data$n_pe[which(data$n_pe == 
                                                                                                                                                1)])))
  if (n_lin_vars > 0) {
    lrt_ci_mat <- matrix(nrow = n_lin_vars, ncol = 2)
    rownames(lrt_ci_mat) <- names[1:n_lin_vars]
    colnames(lrt_ci_mat) <- c("lower .95", "upper .95")
    for (i in 1:n_lin_vars) {
      attr(res, "details")
      beta_good <- attr(res, "coef")
      llik <- attr(res, "details")$value
      prob <- 0.95
      llim <- constr * 0.999
      g <- function(x) {
        aux <- beta_good
        aux <- as.list(aux)
        aux[[i]] <- x
        names(aux) <- paste0("x", 1:length(aux))
        y <- do.call(p.est, aux)
        return(-y + llik + qchisq(prob, 1)/2)
      }
      l1_t <- try(uniroot(g, lower = llim[i], upper = beta_good[i], 
                          extendInt = "no")$root)
      l2_t <- try(uniroot(g, lower = beta_good[i], upper = 1e+05, 
                          extendInt = "yes")$root)
      if (is.numeric(l1_t)) 
        l1 <- l1_t
      else l1 <- -Inf
      if (is.numeric(l2_t)) 
        l2 <- l2_t
      else l2 <- Inf
      if (beta_good[i] > l2) 
        l2 <- Inf
      lrt_ci_mat[i, 1] <- l1
      lrt_ci_mat[i, 2] <- l2
    }
    attr(res, "lrt_ci") <- lrt_ci_mat
    sum <- summary(res)
    coef <- attr(sum, "coef")[1:n_lin_vars, 1]
    se <- attr(sum, "coef")[1:n_lin_vars, 2]
    coef_mat <- matrix(nrow = n_lin_vars, ncol = 4)
    rownames(coef_mat) <- names[1:n_lin_vars]
    colnames(coef_mat) <- c("coef", "se(coef)", "z", "Pr(>|z|)")
    coef_mat[, 1] <- coef
    coef_mat[, 2] <- se
    coef_mat[, 3] <- coef_mat[, 1]/coef_mat[, 2]
    coef_mat[, 4] <- 2 * (1 - pnorm(abs(coef_mat[, 3])))
    attr(res, "lin_coef") <- coef_mat
  }
  if (n_loglin_vars > 0) {
    sum <- summary(res)
    coef <- attr(sum, "coef")[(n_lin_vars + 1):(n_lin_vars + 
                                                  n_loglin_vars), 1]
    se <- attr(sum, "coef")[(n_lin_vars + 1):(n_lin_vars + 
                                                n_loglin_vars), 2]
    coef_mat <- matrix(nrow = n_loglin_vars, ncol = 5)
    rownames(coef_mat) <- names[(n_lin_vars + 1):(n_lin_vars + 
                                                    n_loglin_vars)]
    colnames(coef_mat) <- c("coef", "exp(coef)", "se(coef)", 
                            "z", "Pr(>|z|)")
    coef_mat[, 1] <- coef
    coef_mat[, 2] <- exp(coef)
    coef_mat[, 3] <- se
    coef_mat[, 4] <- coef_mat[, 1]/coef_mat[, 3]
    coef_mat[, 5] <- 2 * (1 - pnorm(abs(coef_mat[, 4])))
    conf_mat <- matrix(nrow = n_loglin_vars, ncol = 4)
    colnames(conf_mat) <- c("exp(coef)", "exp(-coef)", "lower .95", 
                            "upper .95")
    rownames(conf_mat) <- names[(n_lin_vars + 1):(n_lin_vars + 
                                                    n_loglin_vars)]
    conf_mat[, 1] <- exp(coef)
    conf_mat[, 2] <- exp(-coef)
    conf_mat[, 3] <- exp(coef_mat[, 1] - qnorm(0.975, mean = 0, 
                                               sd = 1) * se)
    conf_mat[, 4] <- exp(coef_mat[, 1] + qnorm(0.975, mean = 0, 
                                               sd = 1) * se)
    attr(res, "wald_ci") <- conf_mat
    attr(res, "loglin_coef") <- coef_mat
  }
  lrt_stat <- 2 * (-attr(res, "details")$value - sum(log(1/unlist(lapply(rsets, 
                                                                         length)))))
  lrt_pval <- 1 - pchisq(lrt_stat, 1, lower.tail = T)
  attr(res, "formula") <- formula
  attr(res, "lrt") <- list(lrt_stat = lrt_stat, lrt_pval = lrt_pval)
  attr(res, "se") <- attr(summary(res), "coef")[, 2]
  attr(res, "AIC") <- AIC(res)
  class(res) <- "rERR"
  return(res)
}

#///////////////////////////////////////////////////////////////////////////////////////////////

f_risksets <- function (formula, data, lag, id_name, time_name) {
  # canvis per Jordi F:
  #	- comentada la linia 80 que elimina tacs
  #	- comentada linia 82 i substituida per linies 83-88. depenent si failtime - lag és <0 o >=0
  #			sel.leccionarem tacs reals o la linia ficticia n_pe = -1
  
  id <- NULL
  n_row <- NULL
  formula_sv <- formula[[2]]
  formula_terms <- sum(gregexpr("+", paste0(as.character(formula[[3]]), 
                                            collapse = ""), fixed = TRUE)[[1]] > 0) + 1
  if (formula_terms == 1) {
    if (any(grepl("lin", formula[[3]]))) 
      formula_lin <- formula[[3]]
    if (any(grepl("logl", formula[[3]]))) 
      formula_loglin <- formula[[3]]
    if (any(grepl("strata", formula[[3]]))) 
      formula_strat <- formula[[3]]
    if (any(grepl("logl", formula_lin))) 
      rm(formula_lin)
  }
  if (formula_terms == 2) {
    for (i in 2:3) {
      if (any(grepl("lin", formula[[3]][[i]]))) 
        formula_lin <- formula[[3]][[i]]
      if (any(grepl("logl", formula[[3]][[i]]))) 
        formula_loglin <- formula[[3]][[i]]
      if (any(grepl("strata", formula[[3]][[i]]))) 
        formula_strat <- formula[[3]][[i]]
    }
    
    if (any(grepl("logl", formula_lin))) 
      rm(formula_lin)
  }
  if (formula_terms == 3) {
    if (any(grepl("lin", formula[[3]][[3]]))) 
      formula_lin <- formula[[3]][[3]]
    if (any(grepl("logl", formula[[3]][[3]]))) 
      formula_loglin <- formula[[3]][[3]]
    if (any(grepl("strata", formula[[3]][[3]]))) 
      formula_strat <- formula[[3]][[3]]
    for (i in 2:3) {
      if (any(grepl("lin", formula[[3]][[2]][[i]]))) 
        formula_lin <- formula[[3]][[2]][[i]]
      if (any(grepl("logl", formula[[3]][[2]][[i]]))) 
        formula_loglin <- formula[[3]][[2]][[i]]
      if (any(grepl("strata", formula[[3]][[2]][[i]]))) 
        formula_strat <- formula[[3]][[2]][[i]]
    }
    if (any(grepl("logl", formula_lin))) 
      rm(formula_lin)
  }
  if (exists("formula_lin")) {
    lin_vars <- unlist(strsplit(as.character(formula_lin)[2:length(formula_lin)], 
                                split = "+", fixed = T))
    lin_vars <- gsub(" ", "", lin_vars)
  }
  if (exists("formula_loglin")) {
    loglin_vars <- unlist(strsplit(as.character(formula_loglin)[2:length(formula_loglin)], 
                                   split = "+", fixed = T))
    loglin_vars <- gsub(" ", "", loglin_vars)
  }
  if (exists("formula_strat")) {
    strata_vars <- unlist(strsplit(as.character(formula_strat)[2:length(formula_strat)], 
                                   split = "+", fixed = T))
    strata_vars <- gsub(" ", "", strata_vars)
  }
  v_id <- eval(parse(text = paste0("data$", id_name)))
  v_n_pe <- eval(parse(text = paste0("data$", "n_pe")))
  v_entry <- eval(parse(text = paste0("data$", formula_sv[[2]])))
  v_exit <- eval(parse(text = paste0("data$", formula_sv[[3]])))
  v_outcome <- eval(parse(text = paste0("data$", formula_sv[[4]])))
  v_time <- eval(parse(text = paste0("data$", time_name)))
  nrow_cases <- which(v_outcome == 1)
  failtimes <- v_exit[nrow_cases]
  id_cases <- v_id[nrow_cases]
  rsets <- list()
  nrows_cases_ <- vector()
  for (i in 1:length(failtimes)) {
    dt <- data[which(v_entry < failtimes[i] & v_exit >= failtimes[i]), 1:7]
    #dt <- dt[which(dt[, 7] <= failtimes[i] - lag), ]
    names(dt)[c(2, 4:6, 7)] <- c("id", "entry", "exit", "outcome", "time")
    #dt <- dt %>% group_by(id) %>% dplyr::summarize(row = max(n_row))
    if ((failtimes[i] - lag) < 0){
      dt <- dt[which(dt[, 7] <=0 & dt[, 3] == -1), ] %>% group_by(id) %>% dplyr::summarize(row = max(n_row))
    } else {
      dt <- dt[which(dt[, 7] <= failtimes[i] - lag), ] %>% group_by(id) %>% 
        dplyr::summarize(row = max(n_row))
    }
    nrows_cases_[i] <- dt$row[which(dt$id == id_cases[i])]
    dt$IN <- T
    if (exists("formula_strat")) 
      for (j in 1:length(strata_vars)) {
        x <- eval(parse(text = paste0("data$", strata_vars[j])))
        x_case <- x[nrow_cases[i]]
        dt$x <- x[dt$row]
        dt$IN <- dt$x == x_case
        dt <- dt[dt$IN, ]
      }
    rsets[[i]] <- dt$row
  }
  names(rsets) <- nrows_cases_
  rsets_2 <- list()
  for (i in 1:length(rsets)) {
    rsets_2[[i]] <- rsets[[order(failtimes)[i]]]
    names(rsets_2)[i] <- names(rsets)[order(failtimes)[i]]
  }
  return(rsets_2)
}

#///////////////////////////////////////////////////////////////////////////////////////////////

f_to_event_table_ef_v2 <- function (id, start, stop, outcome, data, times, doses, covars) 
{
  # canvis per Jordi F:
  # - noves linies: de 25 a 29 - afegim una 1a linia per tots els subjectes and age i dosi a 0
  # - modificació linia 33: incoporem n_pe per ordenar data, deixant nova linia n_pe = -1 la primera
  # - nova linia 38: posem les dosis dels tacs que cauen en lagged time a 0
  
  id_ <- NULL
  time_aux <- NULL
  call <- match.call()
  id_name <- eval(id)
  dose_name <- eval(doses)
  stop_name <- eval(stop)
  time_name <- eval(times)
  data$id_ <- eval(parse(text = paste0(call$data, "$", id_name)))
  #data$id_ <- eval(parse(text = paste0("data$", id_name)))
  data <- data %>% group_by(id_) %>% dplyr::mutate(n_pe = 1:length(id_))
  dt <- data[which(data$n_pe == 1), ]
  dt$n_pe <- 0
  dt[, which(names(dt) == dose_name)] <- 0
  dt[, which(names(dt) == time_name)] <- eval(parse(text = paste0("dt$",stop_name)))
  data[, which(names(data) == outcome)] <- 0
  
  data <- rbind(data, dt)
  ## jordi: add new row, unexposed
  dt <- data[which(data$n_pe == 1), ]
  dt$n_pe <- -1
  dt[, which(names(dt) == dose_name)] <- 0
  dt[, which(names(dt) == time_name)] <- 0
  data <- rbind(data, dt)
  ##
  data$time_aux <- eval(parse(text = paste0("data$", time_name)))
  #data <- plyr::arrange(data, id_, time_aux)
  data <- plyr::arrange(data, id_, time_aux, n_pe)
  data <- data[, -dim(data)[2]]
  dose_num <- eval(parse(text = paste0("data$", dose_name)))
  data$dose_num <- dose_num
  # jordi: doses to 0 in lagged CTs
  data[eval(parse(text = paste0("data$", stop_name))) - eval(parse(text = paste0("data$", time_name))) < lag, c("dose_num")] <- 0
  data <- data %>% group_by(id_) %>% dplyr::mutate(dose_cum = cumsum(dose_num))
  a <- plyr::count(data$id_)
  data$id1 <- unlist(lapply(seq_along(1:dim(dt)[1]), function(x) rep(x, 
                                                                     each = a$freq[x])))
  return(data)
}

#///////////////////////////////////////////////////////////////////////////////////////////////

f_fit_linERR_all <- function (formula, data, id_name, time_name, lag) 
{
  dt2 <- f_to_model_data(formula, data, id_name, time_name)
  n_lin_vars <- attr(dt2, "n_lin_vars")
  n_loglin_vars <- attr(dt2, "n_loglin_vars")
  rsets <- f_risksets(formula, data = dt2, lag, id_name, time_name)
  fit <- f_fit_linERR(formula, data = dt2, rsets, n_lin_vars, 
                      n_loglin_vars, id_name, time_name)
  return(fit)
}

#///////////////////////////////////////////////////////////////////////////////////////////////

f_to_model_data <- function (formula, data, id_name, time_name) 
{
  id_aux <- NULL
  time_aux <- NULL
  formula_sv <- formula[[2]]
  formula_terms <- sum(gregexpr("+", paste0(as.character(formula[[3]]), 
                                            collapse = ""), fixed = TRUE)[[1]] > 0) + 1
  if (formula_terms == 1) {
    if (any(grepl("lin", formula[[3]]))) 
      formula_lin <- formula[[3]]
    if (any(grepl("logl", formula[[3]]))) 
      formula_loglin <- formula[[3]]
    if (any(grepl("strata", formula[[3]]))) 
      formula_strat <- formula[[3]]
    if (any(grepl("logl", formula_lin))) 
      rm(formula_lin)
  }
  if (formula_terms == 2) {
    for (i in 2:3) {
      if (any(grepl("lin", formula[[3]][[i]]))) 
        formula_lin <- formula[[3]][[i]]
      if (any(grepl("logl", formula[[3]][[i]]))) 
        formula_loglin <- formula[[3]][[i]]
      if (any(grepl("strata", formula[[3]][[i]]))) 
        formula_strat <- formula[[3]][[i]]
    }
    if (any(grepl("logl", formula_lin))) 
      rm(formula_lin)
  }
  if (formula_terms == 3) {
    if (any(grepl("lin", formula[[3]][[3]]))) 
      formula_lin <- formula[[3]][[3]]
    if (any(grepl("logl", formula[[3]][[3]]))) 
      formula_loglin <- formula[[3]][[3]]
    if (any(grepl("strata", formula[[3]][[3]]))) 
      formula_strat <- formula[[3]][[3]]
    for (i in 2:3) {
      if (any(grepl("lin", formula[[3]][[2]][[i]]))) 
        formula_lin <- formula[[3]][[2]][[i]]
      if (any(grepl("logl", formula[[3]][[2]][[i]]))) 
        formula_loglin <- formula[[3]][[2]][[i]]
      if (any(grepl("strata", formula[[3]][[2]][[i]]))) 
        formula_strat <- formula[[3]][[2]][[i]]
    }
    if (any(grepl("logl", formula_lin))) 
      rm(formula_lin)
  }
  if (exists("formula_lin")) {
    lin_vars <- unlist(strsplit(as.character(formula_lin)[2:length(formula_lin)], 
                                split = "+", fixed = T))
    lin_vars <- gsub(" ", "", lin_vars)
  }
  if (exists("formula_loglin")) {
    loglin_vars <- unlist(strsplit(as.character(formula_loglin)[2:length(formula_loglin)], 
                                   split = "+", fixed = T))
    loglin_vars <- gsub(" ", "", loglin_vars)
  }
  if (exists("formula_strat")) {
    strata_vars <- unlist(strsplit(as.character(formula_strat)[2:length(formula_strat)], 
                                   split = "+", fixed = T))
    strata_vars <- gsub(" ", "", strata_vars)
  }
  v_id <- eval(parse(text = paste0("data$", id_name)))
  v_n_pe <- eval(parse(text = paste0("data$", "n_pe")))
  v_entry <- eval(parse(text = paste0("data$", formula_sv[[2]])))
  v_exit <- eval(parse(text = paste0("data$", formula_sv[[3]])))
  v_outcome <- eval(parse(text = paste0("data$", formula_sv[[4]])))
  v_time <- eval(parse(text = paste0("data$", time_name)))
  dt <- data.frame(v_id, v_n_pe, v_entry, v_exit, v_outcome, 
                   v_time)
  names(dt) <- c(id_name, "n_pe", as.character(formula_sv[[2]]), 
                 as.character(formula_sv[[3]]), as.character(formula_sv[[4]]), 
                 time_name)
  n_lin_vars <- 0
  if (exists("formula_lin")) {
    for (i in 1:length(lin_vars)) {
      is_factor <- F
      if (grepl("factor", lin_vars[i])) {
        lin_vars[i] <- substr(lin_vars[i], 8, nchar(lin_vars[i]) - 
                                1)
        is_factor <- T
      }
      x <- eval(parse(text = paste0("data$", lin_vars[i])))
      if (is.factor(x) | is_factor) {
        x <- as.factor(x)
        levels <- levels(x)
        for (j in 2:length(levels)) {
          n_lin_vars <- n_lin_vars + 1
          x_lev <- as.numeric(x == levels[j])
          dt <- cbind(dt, x_lev)
          names(dt)[dim(dt)[2]] <- paste0(lin_vars[i], 
                                          "_", levels[j])
        }
      }
      else {
        n_lin_vars <- n_lin_vars + 1
        dt <- cbind(dt, x)
        names(dt)[dim(dt)[2]] <- lin_vars[i]
      }
    }
  }
  n_loglin_vars <- 0
  if (exists("formula_loglin")) {
    for (i in 1:length(loglin_vars)) {
      is_factor <- F
      if (grepl("factor", loglin_vars[i])) {
        loglin_vars[i] <- substr(loglin_vars[i], 8, nchar(loglin_vars[i]) - 
                                   1)
        is_factor <- T
      }
      x <- eval(parse(text = paste0("data$", loglin_vars[i])))
      if (is.factor(x) | is_factor) {
        x <- as.factor(x)
        levels <- levels(x)
        for (j in 2:length(levels)) {
          n_loglin_vars <- n_loglin_vars + 1
          x_lev <- as.numeric(x == levels[j])
          dt <- cbind(dt, x_lev)
          names(dt)[dim(dt)[2]] <- paste0(loglin_vars[i], 
                                          "_", levels[j])
        }
      }
      else {
        n_loglin_vars <- n_loglin_vars + 1
        dt <- cbind(dt, x)
        names(dt)[dim(dt)[2]] <- loglin_vars[i]
      }
    }
  }
  if (exists("formula_strat")) {
    for (i in 1:length(strata_vars)) {
      x <- eval(parse(text = paste0("data$", strata_vars[i])))
      dt <- cbind(dt, x)
      names(dt)[dim(dt)[2]] <- strata_vars[i]
    }
  }
  dt <- cbind(dt, dt[, c(1, 6)])
  names(dt)[(dim(dt)[2] - 1):dim(dt)[2]] <- c("id_aux", "time_aux")
  dt <- plyr::arrange(dt, id_aux, time_aux)
  n_row <- data.frame(n_row = 1:(dim(dt)[1]))
  dt <- cbind(n_row, dt)
  dt <- dt[, c(1:(dim(dt)[2] - 2))]
  attr(dt, "n_lin_vars") <- n_lin_vars
  attr(dt, "n_loglin_vars") <- n_loglin_vars
  return(dt)
}

#///////////////////////////////////////////////////////////////////////////////////////////////

f_parse_formula <- function (formula) {
  formula_sv <- formula[[2]]
  if (exists("formula_lin")) 
    rm(formula_lin)
  if (exists("formula_loglin")) 
    rm(formula_loglin)
  if (exists("formula_strat")) 
    rm(formula_strat)
  if (exists("lin_vars")) 
    rm(lin_vars)
  if (exists("loglin_vars")) 
    rm(loglin_vars)
  if (exists("strata_vars")) 
    rm(strata_vars)
  formula_terms <- sum(gregexpr("+", paste0(as.character(formula[[3]]), 
                                            collapse = ""), fixed = TRUE)[[1]] > 0) + 1
  if (formula_terms == 1) {
    if (any(grepl("lin", formula[[3]]))) 
      formula_lin <- formula[[3]]
    if (any(grepl("logl", formula[[3]]))) 
      formula_loglin <- formula[[3]]
    if (any(grepl("strata", formula[[3]]))) 
      formula_strat <- formula[[3]]
    if (exists("formula_lin") & exists("formula_loglin")) 
      if (formula_lin == formula_loglin) 
        rm(formula_lin)
  }
  if (formula_terms == 2) {
    for (i in 2:3) {
      if (any(grepl("lin", formula[[3]][[i]]))) 
        formula_lin <- formula[[3]][[i]]
      if (any(grepl("logl", formula[[3]][[i]]))) 
        formula_loglin <- formula[[3]][[i]]
      if (any(grepl("strata", formula[[3]][[i]]))) 
        formula_strat <- formula[[3]][[i]]
    }
    if (exists("formula_lin") & exists("formula_loglin")) 
      if (formula_lin == formula_loglin) 
        rm(formula_lin)
  }
  if (formula_terms == 3) {
    if (any(grepl("lin", formula[[3]][[3]]))) 
      formula_lin <- formula[[3]][[3]]
    if (any(grepl("logl", formula[[3]][[3]]))) 
      formula_loglin <- formula[[3]][[3]]
    if (any(grepl("strata", formula[[3]][[3]]))) 
      formula_strat <- formula[[3]][[3]]
    for (i in 2:3) {
      if (any(grepl("lin", formula[[3]][[2]][[i]]))) 
        formula_lin <- formula[[3]][[2]][[i]]
      if (any(grepl("logl", formula[[3]][[2]][[i]]))) 
        formula_loglin <- formula[[3]][[2]][[i]]
      if (any(grepl("strata", formula[[3]][[2]][[i]]))) 
        formula_strat <- formula[[3]][[2]][[i]]
    }
  }
  if (exists("formula_lin")) {
    lin_vars <- unlist(strsplit(as.character(formula_lin)[2:length(formula_lin)], 
                                split = "+", fixed = T))
    lin_vars <- gsub(" ", "", lin_vars)
  }
  if (exists("formula_loglin")) {
    loglin_vars <- unlist(strsplit(as.character(formula_loglin)[2:length(formula_loglin)], 
                                   split = "+", fixed = T))
    loglin_vars <- gsub(" ", "", loglin_vars)
  }
  if (exists("formula_strat")) {
    strata_vars <- unlist(strsplit(as.character(formula_strat)[2:length(formula_strat)], 
                                   split = "+", fixed = T))
    strata_vars <- gsub(" ", "", strata_vars)
  }
  res <- list()
  res[[1]] <- list(entry = formula_sv[[2]], exit = formula_sv[[3]], 
                   outcome = formula_sv[[4]])
  if (exists("formula_lin")) {
    res[[2]] <- lin_vars
  }
  else res[[2]] <- character(0)
  if (exists("formula_loglin")) {
    res[[3]] <- loglin_vars
  }
  else res[[3]] <- character(0)
  if (exists("formula_strat")) {
    res[[4]] <- strata_vars
  }
  else res[[4]] <- character(0)
  res[[5]] <- character(0)
  res[[6]] <- character(0)
  if (exists("formula_lin")) {
    res[[5]] <- character(0)
    for (i in 1:length(lin_vars)) {
      if (grepl("factor", lin_vars[i])) 
        res[[5]] <- append(res[[5]], substr(lin_vars[i], 
                                            8, nchar(lin_vars[i]) - 1))
    }
  }
  if (exists("formula_loglin")) {
    res[[6]] <- character(0)
    for (i in 1:length(loglin_vars)) {
      if (grepl("factor", loglin_vars[i])) 
        res[[6]] <- append(res[[6]], substr(loglin_vars[i], 
                                            8, nchar(loglin_vars[i]) - 1))
    }
  }
  for (i in 1:length(res)) if (is.null(res[[i]]) | length(res[[i]]) == 
                               0) 
    res[[i]] <- NA
  names(res) <- c("Surv", "lin_vars", "loglin_vars", "strata_vars", 
                  "lin_factor", "loglin_factor")
  if (exists("formula_lin")) 
    rm(formula_lin)
  if (exists("formula_loglin")) 
    rm(formula_loglin)
  if (exists("formula_strat")) 
    rm(formula_strat)
  if (exists("lin_vars")) 
    rm(lin_vars)
  if (exists("loglin_vars")) 
    rm(loglin_vars)
  if (exists("strata_vars")) 
    rm(strata_vars)
  return(res)
}

#///////////////////////////////////////////////////////////////////////////////////////////////

f_to_event_table_ef_all <- function (formula, data, id_name, dose_name, time_name, covars_names) 
{
  form <- f_parse_formula(formula)
  entry_name <- as.character(form$Surv$entry)
  exit_name <- as.character(form$Surv$exit)
  #exit_name <- as.character("exit_age")
  outcome_name <- as.character(form$Surv$outcome)
  dt1 <- f_to_event_table_ef_v2(id = id_name, start = entry_name, 
                                stop = exit_name, outcome = outcome_name, data = data, 
                                doses = dose_name, times = time_name, covars = covars_names)
  return(dt1)
}
