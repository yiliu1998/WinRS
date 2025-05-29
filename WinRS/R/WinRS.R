#' Estimation and inference for the win ratio for two hierarchical endpoints using the S score estimator (nonparametric maximum likelihood estimator)
#' @param df_a data frame for data of treatment group a (treated group)
#' @param df_b data frame for data of treatment group a (treated group), the structure and names of variables should match df_a
#' @param cov.adjust whether to adjust covariate for missingness of the second endpoint Y2 (the quality-of-life outcome), default is FALSE
#' @param covar.name a vector for names of columns of covariates in df_a and df_b; if cov.adjust=TRUE, one must specify a set of covariates
#' @param t.horizon time horizon, a positive finite constant specified as the administrative end time of the study
#' @param out1.name column name of the primary survival outcome in both df_a and df_b
#' @param delta1.name column name of the primary event indicator in both df_a and df_b (1 indicates event observed, while 0 indicates right censored)
#' @param out2.name column name of the second outcome (non-survival) in both df_a and df_b
#' @param boot whether to conduct bootstrap for variance estimation, default is TRUE
#' @param n_boot number of bootstrap replicates if boot=TRUE, default is 500
#' @param seed seed for generating random numbers used in set.seed() function when splitting data; default is 4399
WinRS <- function(df_a, df_b,
                  cov.adjust=FALSE, covar.name=NULL,
                  t.horizon, out1.name, delta1.name,
                  out2.name,
                  boot=TRUE, n_boot=500, seed=4399) {

  set.seed(seed=seed)
  h <- t.horizon
  data_a <- get_S_df(data=df_a, t.horizon=h, covar.name, out1.name, delta1.name, out2.name)
  data_b <- get_S_df(data=df_b, t.horizon=h, covar.name, out1.name, delta1.name, out2.name)

  if(!cov.adjust) {
    km_a <- fit_kaplan_meier_S(data_a)
    km_b <- fit_kaplan_meier_S(data_b)

    common_times <- sort(unique(c(km_a$time, km_b$time)))
    K <- length(common_times)

    S_a <- approx(c(0,km_a$time), c(1,km_a$surv), xout=common_times, rule=2, method="constant")$y
    S_b <- approx(c(0,km_b$time), c(1,km_b$surv), xout=common_times, rule=2, method="constant")$y

    diff_S_b <- diff(c(0, 1-S_b))
    diff_S_a <- diff(c(0, 1-S_a))

    N <- sum(S_a * diff_S_b)
    D <- sum(S_b * diff_S_a)
    theta <- N / D

    K_a <- nrow(data_a)
    K_b <- nrow(data_b)
    n <- K_a+K_b
    F_i <- numeric(n)

    denom_a <- sapply(common_times[-K], function(t_k) mean(data_a$S >= t_k))*(K_a/n)
    denom_b <- sapply(common_times[-K], function(t_k) mean(data_b$S >= t_k))*(K_b/n)
    denom_a[denom_a==0] <- 0.0001
    denom_b[denom_b==0] <- 0.0001

    sum_mask <- upper.tri(matrix(1, K-1, K-1), diag=TRUE)

    compute_F_vectorized <- function(data, km, denom, group) {
      num_obs <- nrow(data)
      F_temp <- numeric(num_obs)
      time_vec <- data$S
      delta_vec <- data$Delta_S

      if (group=="a") {
        q_k <- S_a / c(1, S_a[-K])
        weight_sums <- sum_mask %*% (S_a[1:(K-1)] / D^2 *
                                       (  N * (S_b[1:(K-1)]-S_b[2:K]) +
                                            D * (c(1, S_b[1:(K-2)])-S_b[1:(K-1)])))
      } else {
        q_k <- S_b / c(1, S_b[-K])
        weight_sums <- sum_mask %*% (S_b[1:(K-1)] / D^2 *
                                       ( -D * (S_a[1:(K-1)]-S_a[2:K]) -
                                           N * (c(1, S_a[1:(K-2)])-S_a[1:(K-1)])))
      }

      Y_mat <- outer(time_vec, common_times[-K], FUN=">=")
      A_mat <- outer(time_vec, common_times[-K], FUN=function(t, tk) (t > tk | (t==tk & delta_vec==0)))
      for (i in 1:num_obs) {
        E_i_k <- (Y_mat[i, ] * (A_mat[i, ]-q_k[-K])) / denom
        F_temp[i] <- sum((1 / q_k[-K]) * E_i_k * weight_sums, na.rm=T)
      }
      return(F_temp)
    }
    F_i[1:K_a] <- compute_F_vectorized(data_a, km_a, denom_a, group="a")
    F_i[(K_a+1):n] <- compute_F_vectorized(data_b, km_b, denom_b, group="b")

    theta.sd <- sqrt(mean(F_i^2, na.rm=TRUE) / n)

    boot.sd <- boot.qt.lwr <- boot.qt.upr <- NA
    if(boot) {
      if(is.na(n_boot)) { stop("Please specify the number of bootstrap replicates!") }
      bootstrap_wr <- numeric(n_boot)
      for (b in 1:n_boot) {
        resampled_a <- data_a[sample(1:nrow(data_a), replace=TRUE), ]
        resampled_b <- data_b[sample(1:nrow(data_b), replace=TRUE), ]
        bootstrap_wr[b] <- WinRS_est(resampled_a, resampled_b, cov.adjust=FALSE, covar.name=NULL, h=h)
      }
      boot.sd=sd(bootstrap_wr)
      boot.qt.lwr=quantile(bootstrap_wr, 0.025)
      boot.qt.upr=quantile(bootstrap_wr, 0.975)
    }
  }

  if(cov.adjust) {
    km_a1 <- fit_kaplan_meier_Y1(data_a)
    km_b1 <- fit_kaplan_meier_Y1(data_b)

    cdf_a2 <- fit_cdf_Y2(data_a, h=h, covar=covar)
    cdf_b2 <- fit_cdf_Y2(data_b, h=h, covar=covar)

    F1a_h <- 1-km_a1$surv[km_a1$time==max(km_a1$time[km_a1$time!=(h+1)])]
    F1b_h <- 1-km_b1$surv[km_b1$time==max(km_b1$time[km_b1$time!=(h+1)])]

    cdf_Sa <- c(1-km_a1$surv[km_a1$time<=h],F1a_h + (1-F1a_h)*cdf_a2$F2)
    cdf_Sb <- c(1-km_b1$surv[km_b1$time<=h],F1b_h + (1-F1b_h)*cdf_b2$F2)

    common_times <- sort(unique(c(km_a1$time,cdf_a2$time2+h+1,km_b1$time,cdf_b2$time2+h+1)))
    K <- length(common_times)

    S_a <- approx(c(km_a1$time,cdf_a2$time2+h+1),c(1,1-cdf_Sa),xout=common_times,rule=2,method="constant")$y
    S_b <- approx(c(km_b1$time,cdf_b2$time2+h+1),c(1,1-cdf_Sb),xout=common_times,rule=2,method="constant")$y

    N <- sum(S_a*diff(c(0,1-S_b)))
    D <- sum(S_b*diff(c(0,1-S_a)))
    theta <- N / D

    K_a <- nrow(data_a)
    K_b <- nrow(data_b)
    n <- K_a+K_b
    F_i <- numeric(n)

    ind_a <- which(data_a$Y1>h)
    ind_b <- which(data_b$Y1>h)
    ipw_a <- rep(1,K_a)
    ipw_b <- rep(1,K_b)
    ipw_a[ind_a] <- cdf_a2$ipw
    ipw_b[ind_b] <- cdf_b2$ipw

    denom_a <- sapply(common_times[-K],function(t_k) mean((data_a$S>=t_k)*ipw_a))*(K_a/n)
    denom_b <- sapply(common_times[-K],function(t_k) mean((data_b$S>=t_k)*ipw_b))*(K_b/n)
    denom_a[denom_a==0] <- 0.0001
    denom_b[denom_b==0] <- 0.0001

    sum_mask <- upper.tri(matrix(1,K-1,K-1),diag=TRUE)

    compute_F_vectorized <- function(data, denom, group) {
      num_obs <- nrow(data)
      F_temp <- numeric(num_obs)
      time_vec <- data$S
      delta_vec <- data$Delta_S

      if (group=="a") {
        q_k <- S_a / c(1,S_a[-K])
        weight_sums <- sum_mask %*% (S_a[1:(K-1)] / D^2 *
                                       (  N * (S_b[1:(K-1)]-S_b[2:K]) +
                                            D * (c(1,S_b[1:(K-2)])-S_b[1:(K-1)])))

        Y_mat <- outer(time_vec,common_times[-K],FUN=">=")
        A_mat <- outer(time_vec,common_times[-K],FUN=function(t,tk) (t > tk | (t==tk & delta_vec==0)))
        for (i in 1:num_obs) {
          E_i_k <- (Y_mat[i,] * (A_mat[i,]-q_k[-K]) * ipw_a[i]) / denom
          F_temp[i] <- sum((1 / q_k[-K]) * E_i_k * weight_sums, na.rm=T)
        }
      } else {
        q_k <- S_b / c(1,S_b[-K])
        weight_sums <- sum_mask %*% (S_b[1:(K-1)] / D^2 *
                                       ( -D * (S_a[1:(K-1)]-S_a[2:K]) -
                                           N * (c(1,S_a[1:(K-2)])-S_a[1:(K-1)])))

        Y_mat <- outer(time_vec,common_times[-K],FUN=">=")
        A_mat <- outer(time_vec,common_times[-K],FUN=function(t,tk) (t > tk | (t==tk & delta_vec==0)))
        for (i in 1:num_obs) {
          E_i_k <- (Y_mat[i,] * (A_mat[i,]-q_k[-K]) * ipw_b[i]) / denom
          F_temp[i] <- sum((1 / q_k[-K]) * E_i_k * weight_sums, na.rm=T)
        }
      }
      return(F_temp)
    }
    F_i[1:K_a] <- compute_F_vectorized(data_a, denom_a, group="a")
    F_i[(K_a+1):n] <- compute_F_vectorized(data_b, denom_b, group="b")

    theta.sd <- sqrt(mean(F_i^2,na.rm=TRUE) / n)

    boot.sd <- boot.qt.lwr <- boot.qt.upr <- NA
    if(boot) {
      if(is.na(n_boot)) { stop("Please specify the number of bootstrap replicates!") }
      bootstrap_wr <- numeric(n_boot)
      for (b in 1:n_boot) {
        resampled_a <- data_a[sample(1:nrow(data_a), replace=TRUE), ]
        resampled_b <- data_b[sample(1:nrow(data_b), replace=TRUE), ]
        bootstrap_wr[b] <- WinRS_est(resampled_a, resampled_b, cov.adjust=TRUE, covar.name=covar.name, h=h)
      }
      boot.sd=sd(bootstrap_wr)
      boot.qt.lwr=quantile(bootstrap_wr, 0.025)
      boot.qt.upr=quantile(bootstrap_wr, 0.975)
    }
  }

  return(list(WR.est=theta,
              IF.sd=theta.sd,
              IF.lwr=theta-theta.sd*1.96,
              IF.upr=theta+theta.sd*1.96,
              boot.sd=boot.sd,
              boot.Wald.lwr=theta-boot.sd*1.96,
              boot.Wald.upr=theta+boot.sd*1.96,
              boot.qt.lwr=boot.qt.lwr,
              boot.qt.upr=boot.qt.upr))
}

WinRS_est <- function(data_a, data_b,
                      cov.adjust=FALSE, covar.name=NULL, h) {

  if(!cov.adjust) {
    km_a <- fit_kaplan_meier_S(data_a)
    km_b <- fit_kaplan_meier_S(data_b)

    common_times <- sort(unique(c(km_a$time, km_b$time)))
    S_a <- approx(c(0,km_a$time), c(1,km_a$surv), xout=common_times, rule=2, method="constant")$y
    S_b <- approx(c(0,km_b$time), c(1,km_b$surv), xout=common_times, rule=2, method="constant")$y
  }

  if(cov.adjust) {
    km_a1 <- fit_kaplan_meier_Y1(data_a)
    km_b1 <- fit_kaplan_meier_Y1(data_b)

    cdf_a2 <- fit_cdf_Y2(data_a, h=h, covar=covar.name)
    cdf_b2 <- fit_cdf_Y2(data_b, h=h, covar=covar.name)

    F1a_h <- 1-km_a1$surv[km_a1$time==max(km_a1$time[km_a1$time!=(h+1)])]
    F1b_h <- 1-km_b1$surv[km_b1$time==max(km_b1$time[km_b1$time!=(h+1)])]

    cdf_Sa <- c(1-km_a1$surv[km_a1$time<=h],F1a_h + (1-F1a_h)*cdf_a2$F2)
    cdf_Sb <- c(1-km_b1$surv[km_b1$time<=h],F1b_h + (1-F1b_h)*cdf_b2$F2)

    common_times <- sort(unique(c(km_a1$time, cdf_a2$time2+h+1, km_b1$time, cdf_b2$time2+h+1)))

    S_a <- approx(c(km_a1$time,cdf_a2$time2+h+1), c(1,1-cdf_Sa), xout=common_times, rule=2, method="constant")$y
    S_b <- approx(c(km_b1$time,cdf_b2$time2+h+1), c(1,1-cdf_Sb), xout=common_times, rule=2, method="constant")$y
  }

  N <- sum(S_a*diff(c(0, 1-S_b)))
  D <- sum(S_b*diff(c(0, 1-S_a)))
  theta <- N / D

  return(theta)
}

get_S_df <- function(data,
                     t.horizon,
                     covar.name,
                     out1.name,
                     delta1.name,
                     out2.name) {

  h <- t.horizon
  T1 <- data[,out1.name]
  Y1 <- pmin(T1, h) + I(T1>h)
  Delta1 <- data[,delta1.name]
  Delta1.tilde <- I(Delta1==1 & Y1<=h) + I(Y1>h)
  Y2 <- data[,out2.name]
  Y2[Y1<=h] <- 0
  R2 <- as.numeric(!is.na(Y2))
  S <- Y1 + I(Y1>h)*ifelse(R2==0, 0, Y2)
  Delta_S <- ifelse((Delta1==1 & Y1<=h) | (R2==1 & Y2>h), 1, 0)
  X <- data[,covar.name]
  df <- as.data.frame(cbind(X,
                            data.frame(Y1=Y1, Delta1=Delta1,
                                       Y2=Y2, R2=R2,
                                       S=S, Delta_S=Delta_S)))
  return(df)
}

fit_kaplan_meier_S <- function(data) {
  survfit(Surv(S, Delta_S) ~ 1, data=data)
}

fit_kaplan_meier_Y1 <- function(data) {
  survfit(Surv(Y1, Delta1) ~ 1, data=data)
}

fit_cdf_Y2 <- function(data, h, covar) {
  dat.Y2 <- data[data$Y1>h,]
  X <- as.matrix(dat.Y2[,covar])
  R2 <- dat.Y2$R2
  Y2 <- dat.Y2$Y2
  Y2[R2==0] <- -999
  ps <- glm(R2~X,family=binomial(link="logit"))
  ipw <- 1/ps$fitted.values
  time2 <- unique(Y2[Y2!=-999]) %>% sort()
  F2 <- c()
  for(i in 1:length(time2)) {
    F2[i] <- sum(R2*ipw*I(Y2<=time2[i])) / sum(R2*ipw)
  }
  return(list(time2=time2,F2=F2,ipw=ipw))
}
