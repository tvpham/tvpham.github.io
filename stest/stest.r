#######################################################################
#
# Author: Thang V. Pham, t.pham@vumc.nl
#
# All rights reserved.
#
# Citation for the s-test:
#
# Pham TV, Jimenez CR.
# Simulated linear test applied to quantitative proteomics.
# Bioinformatics. 2016 Sep 1; 32(17):i702-i709.
#
# Software version: 2.0
#
#######################################################################

stest <- list()

stest$ms_intensity_data <- function(dat, d.log2 = NULL, pdffile = "_stest.pdf", moderate_factor = 1.0, min_factor = 1.0) {

    # dat is the nonlog data for technical estimate. Provide d.log2 for the real data, if data for estimation of technical variation (a, b, cc) is different.

    pdf(pdffile, width = 7, height = 7)

    ret <- list()

    # technical variation
    cc <- stest$estimate_cv(dat, draw = TRUE)
    ab <- stest$estimate_ab(dat, draw = TRUE)

    tvar <- log1p(cc$cc^2)

    ret$technical_variation <- list(s2 = tvar, ab = ab)


    if (is.null(d.log2)) {
        ret$d.log2 <- log2(dat)
        ret$d.log2[ret$d.log2 == -Inf] <- NA
    } else {
        ret$d.log2 <- d.log2
    }

    # sd and degree of freedom
    ret$n <- rowSums(!is.na(ret$d.log2))
    sigma <- apply(ret$d.log2, 1, sd, na.rm = TRUE)
    ret$sn1 <- sigma * (ret$n-1)


    # estimation via the gamma distribution over precision
    precision = 1.0 / sigma^2
    precision <- precision[!is.na(precision)]
    precision <- precision[ precision < quantile(precision, probs = 0.95)]
    a_rate  <- stest$fit_gamma(precision)
    ret$d0 <-  moderate_factor * 2 * a_rate$a
    ret$s0_2 <-  a_rate$rate / a_rate$a

    # plot histogram
    x <- seq(1e-6, 10, length.out = 1000)

    pxa <- exp(ret$d0/2*log(ret$s0_2*ret$d0/2) - lgamma(ret$d0/2) + (-ret$s0_2*ret$d0/2 / x) - (1+ret$d0/2)*log(x))
    #pxa2 <- exp(0.5*d0/2*log(s0_2*0.5*d0/2) - lgamma(0.5*d0/2) + (-s0_2*0.5*d0/2 / x) - (1+0.5*d0/2)*log(x))

    h <- hist(sigma^2, 40, freq = FALSE, ylim = c(0, 2))
    #lines(x, px, col = "red")
    lines(x, pxa, col = "green")
    #lines(x, pxa2, col = "orange")

    # plot estimated variances
    m <- apply(ret$d.log2, 1, mean, na.rm = TRUE)
    plot(m, sigma^2, type = "p", col = "lightblue", pch = 16)
    var_post <- rep(NA, length(sigma))
    for (i in 1:length(var_post)) {
        var_post[i] <- stest$update_chisquare(ret$d.log2[i,], matrix(1, nrow = ncol(ret$d.log2), ncol = 1), ret$d0, ret$s0_2)$s2
    }
    points(m, var_post, col = "red")
    legend("topright", legend = c("Empirical variance", "Squeezed variance"),
           pch = c(16, 1),
           col = c("lightblue", "red"))


    # imputation
    ret$imputation <- list(m = min_factor * min(ret$d.log2, na.rm = TRUE),
                           s2 = mean(apply(ret$d.log2, 1, sd, na.rm = TRUE), na.rm = TRUE)^2)


    dev.off()

    return(ret)
}

stest$ms_intensity_matrix <- function(d.log2, group, sdata, s2 = NULL, cv = NULL, id = NULL, n.threads = 1, index = 1:nrow(d.log2)) {

    # s2 and cv will overwrite estimation in sdata

    if (length(group) != ncol(d.log2)) {
        stop("length of 'group' must be equal to ncol of matrix\n")
    }

    if (length(levels(factor(group))) < 2) {
        stop("minimum 2 groups\n")
    }

    # testing
    d <- d.log2

    if (is.null(id)) {
        X <- model.matrix(~ factor(group))
    } else {
        X <- model.matrix(~ factor(group) + factor(id))
        X_restricted <- model.matrix(~ factor(id))
    }


    tech_var <- sdata$technical_variation$s2

    if (!is.null(s2)) {
        tech_var <- s2
    } else {
        if (!is.null(cv)) {
            tech_var <- log1p(cv^2)
        }
    }

    pval <- rep(1, length(index))
    log_fc <- rep(1, length(index))

    #-- single thread or multiple threads
    if (n.threads == 1) {
        for (i in index) {
            y  <- as.numeric(d[i,])
            if (is.null(id)) {
                out <- stest$ms_intensity_vector(y, X, tech_var, sdata$technical_variation$ab, sdata$imputation$m, sdata$imputation$s2,
                                                 sdata$d0, sdata$s0_2)
            } else {
                out <- stest$ms_intensity_vector(y, X, tech_var, sdata$technical_variation$ab, sdata$imputation$m, sdata$imputation$s2,
                                                 sdata$d0, sdata$s0_2, X_restricted = X_restricted)
            }
            pval[i] <- out$pval
            log_fc[i] <- out$beta[2]
        }
    } else {
        require(parallel)

        nt <- ifelse(n.threads > 1, n.threads, detectCores(logical = FALSE) + n.threads)

        cat("Using", nt, "core(s)\n")

        cl <- makeCluster(nt)
        clusterExport(cl=cl, varlist=c("stest"))

        outer <- parLapply(cl, as.list(index), function(i) {
            y  <- as.numeric(d[i,])

            if (is.null(id)) {
                out <- stest$ms_intensity_vector(y, X, tech_var, sdata$technical_variation$ab, sdata$imputation$m, sdata$imputation$s2,
                                                 sdata$d0, sdata$s0_2)
            } else {
                out <- stest$ms_intensity_vector(y, X, tech_var, sdata$technical_variation$ab, sdata$imputation$m, sdata$imputation$s2,
                                                 sdata$d0, sdata$s0_2, X_restricted = X_restricted)
            }
            out
        })
        stopCluster(cl)

        for (i in 1:length(index)) {
            pval[i] <- outer[[i]]$pval
            log_fc[i] <- outer[[i]]$beta[2]
        }
    }

    pval.BH <- pval

    pval.BH[pval < 1] <- p.adjust(pval[pval < 1], method = "BH")

    if (length(levels(factor(group))) == 2) {
        return (list(dd = d,
                     log_fc = log_fc,
                     pval = pval,
                     pval.BH = pval.BH,
                     n_detected = rowSums(!is.na(d))))
    } else {
        return (list(dd = d,
                     pval = pval,
                     pval.BH = pval.BH,
                     n_detected = rowSums(!is.na(d))))
    }
}

stest$ms_intensity_vector <- function(y, X, techvar2, ab, null_mu, null_sigma2, d0, s0_2,
                                      X_restricted = matrix(1, nrow = length(y), ncol = 1)) {


    if (all(is.na(y))) {
        return(list(pval = 1, beta = rep(0, ncol(X))))
    } else if (is.na(null_mu)) {
        y_ok <- !is.na(y)
        y_new <- y[y_ok]
        X_new <- X[y_ok,, drop = FALSE]
        X_restricted_new <- X_restricted[y_ok,, drop = FALSE]
    } else {
        y_new <- y
        X_new <- X
        X_restricted_new <- X_restricted
    }

    # TODO: for multiple row, this is common, thus can be further optimized, computationally.
    A <- stest$check_M(t(matrix(X_new, ncol = ncol(X_new))))
    A_restricted <- stest$check_M(t(matrix(X_restricted_new, ncol = ncol(X_restricted_new))))
    
    if (is.na(A) || is.na(A_restricted)) {
        return(list(pval = 1, beta = rep(0, ncol(X_new))))
    }
    
    para1 <- stest$update_chisquare(y_new, X_new, d0, s0_2, na_v = null_mu)
    #para1 <- stest$update_chisquare(y, X, d0, s0_2)
    #para2 <- stest$update_chisquare(y, X_restricted, d0, s0_2)

    #s2_hat <- (para1$d*para1$s2 + para2$d*para2$s2) / (para1$d + para2$d)
    s2_hat <- para1$s2

    # ----- simulate
    Lw <- NULL
    U <- NULL

    for (k in 1:length(y_new)) {
        s <-  stest$weighted_sampling_pyu(y_new[k], ab, techvar2, null_mu, null_sigma2)
        #s <-  stest$weighted_sampling_pyu(y[k], ab, techvar2, min(y_ok)-sqrt(para$s2), null_sigma2)
        Lw  <- cbind(Lw, s$log_lambda)
        U  <- cbind(U, s$u)
    }


    #% ----- unrestricted MLE
    uLogL <-  stest$max_LL(Lw, U, t(matrix(X_new, ncol = ncol(X_new))), y_new, s2_hat)

    #% ----- restrictied MLE
    rLogL <-  stest$max_LL(Lw, U, t(matrix(X_restricted_new, ncol = ncol(X_restricted_new))), y_new, s2_hat)


    stat  <-  2 * (uLogL$L - rLogL$L)

    if (stat < 0) {
        stat  <-  0
    }

    #stat

    pval <- 1 - pchisq(stat, ncol(X_new)-ncol(X_restricted_new))

    return(list(pval = pval, beta = uLogL$beta))
}

stest$check_M <- function(M) {
    A <- NA
    tryCatch({
        A  <- solve(M %*% t(M), M)
    }, error = function(e) {
        A <- NA
    })
    A
}

stest$max_LL <- function(Lw, U, X, y, s2) {
    
    
    maxT <- 1e4
    f_tol <- 1e-8
    beta_tol <- 1e-8
    
    
    M <- nrow(U)
    K <- ncol(U)
    
    A  <- solve(X %*% t(X), X)
    
    yy  <-  y
    yy[is.na(y)]  <-  mean(y[!is.na(y)]) # for starting values
    
    beta  <-  A %*% as.matrix(yy)
    
    xbeta  <-  t(t(X) %*% beta)
    
    L <- -Inf
    
    
    for (t in 1:maxT) {
        
        #matplot(t(U), pch=16)
        #cat(t, " -> ", L, "\n")
        #if (readline() == 'q') stop()
        
        
        old_L  <-  L
        old_beta  <-  beta
        
        R  <- U - matrix(xbeta, nrow = M, ncol = length(xbeta), byrow = TRUE)
        R2  <- R * R
        
        H_mk  <-  Lw - 0.5 * log(2 * pi * s2) - R2 / (2 * s2)
        
        Lk  <-  stest$lse(H_mk)
        
        L  <-  sum(Lk)
        
        W_mk  <- exp(H_mk - matrix(Lk, nrow = M, ncol = length(Lk), byrow = TRUE))
        
        
        #matplot(t(U), pch='.')
        #print(colSums(W_mk))
        #symbols(rep(1:K, M), as.vector(t(U)), circles = as.vector(t(W_mk)), inches=0.05, bg="seagreen")
        #cat(t, " -> ", L, "\n")
        #if (readline() == 'q') stop()
        
        v  <- colSums(U * W_mk)
        
        beta  <- A %*% as.matrix(v)
        xbeta  <- t(t(X) %*% beta)
        
        if (max(abs(old_L - L)) < f_tol && max(abs(old_beta - beta)) < beta_tol) {
            break
        }
        
    }
    
    R  <- U - matrix(xbeta, nrow = M, ncol = length(xbeta), byrow = TRUE)
    R2  <- R * R
    
    H_mk  <-  Lw - 0.5 * log(2 * pi * s2) - R2 / (2 * s2)
    
    Lk  <-  stest$lse(H_mk)
    
    L  <-  sum(Lk)
    
    list(L = L, beta = beta)
}


stest$plot_mle <- function(y, techvar2, ab, null_mu, null_sigma2, d0, s0_2) {

    y_ok <- y[!is.na(y)]

    if (all(is.na(y))) {
        return(1.0)
    }

    # ----- simulate
    Lw <- NULL
    U <- NULL

    for (k in 1:length(y)) {
        s <-  stest$weighted_sampling_pyu(y[k], ab, techvar2, null_mu, null_sigma2)
        Lw  <- cbind(Lw, s$log_lambda)
        U  <- cbind(U, s$u)
    }

    #% ----- restricted MLE

    X <- matrix(1, nrow = length(y), ncol = 1)

    para <- stest$update_chisquare(y, X, d0, s0_2)
    logl <-  stest$max_LL(Lw, U, t(X), y, para$s2)


    x <- seq(0, 40, length.out = 1000)
    #px <-function(x, xbeta, d0, s0_2) {
    #    exp(lgamma((d0+1.0)/2.0) - lgamma(d0/2.0) - 0.5*log(pi*d0*s0_2) - 0.5*(d0+1.0)*log1p((x-xbeta)^2/(d0*s0_2)))
    #}
    #cat(matrix(1, nrow = length(y), ncol = 1) %*% rLogL$beta)
    cat("beta=", logl$beta, "\n")
    fx <- x
    for (i in 1:length(fx)) fx[i] <- dnorm(x[i], logl$beta, sqrt(para$s2))
    plot(x, fx, col = "red", type='l', ylim=c(0,1))
    points(x, dnorm(x, mean = mean(y_ok, na.rm = TRUE), sd = sd(y_ok, na.rm = TRUE)) , col = "blue", type = "l")
    points(y, rep(0, length(y)), col = "green", type = "p")
    legend("topleft", legend = c("Empirical", "Moderated", "Data points"),
           lty = c(1,1,NA),
           pch = c(NA,NA, 1),
           col = c("blue", "red", "green"))

}

stest$update_chisquare <- function(y, X, d0, s0_2, na_v = NA) {

    fit <- stest$estimate_chisquare(y, X, na_v)

    if (fit$d0 > 1) {
        return(list(d = d0 + fit$d0, s2 = (d0 * s0_2 + fit$s2 * fit$d0)/(d0 + fit$d0)))
    } else {
        return(list(d = d0, s2 = s0_2))
    }
}


stest$estimate_chisquare <- function(y, X, na_v = NA) {

    s2 <- NA
    d0 <- 0

    if (is.na(na_v)) {
        y_ok <- !is.na(y)

        if (sum(y_ok) > 0) {
            XX <- X[y_ok, , drop = FALSE]
            yy <- y[y_ok]

            out <- lm.fit(XX, yy)

            d0 <- out$df.residual
            
            if (d0 > 0) {
                s2 <- mean(out$effects[-(1:out$rank)]^2)
            }
        }
    } else {

        yy <- y
        yy[is.na(y)] <- na_v
        
        XX <- X
            
        out <- lm.fit(XX, yy)
            
        d0 <- out$df.residual
            
        if (d0 > 0) {
            s2 <- mean(out$effects[-(1:out$rank)]^2)
        }
    
    }
    return(list(d0 = d0, s2 = s2))
}


stest$fit_gamma <- function(x) { # x = precision

    tmp <- mean(log(x)) - log(mean(x))

    a <- 0.5 / (-tmp) # eq. 13
    for (i in 1:1000) {
        old_a <- a
        a <- 1.0 / (1.0/a + (tmp + log(a) - digamma(a))/(a*a*(1/a-trigamma(a))))
        if (abs(a-old_a) < 1e-10 && i > 100) {
            break
        }
    }

    return(list(a = a, rate = a/mean(x)))

}


stest$log_pyu <- function(y, u, ab, sigma2_hat) {

    phi <- exp(u * ab$a + ab$b)

    if (!is.na(y)) {
        ll <- -log(1.0 + phi) - log(sqrt(2.0 * pi * sigma2_hat)) - (y - u)*(y - u) / (2.0 * sigma2_hat)
    } else {
        ll <- log(phi / (1.0 + phi))
    }

    return(ll)
}

stest$weighted_sampling_pyu <- function(y, ab, techvar2, mu, sigma2) {

    z  <- c(-6.995680124, -6.275078705, -5.673961445, -5.133595577, -4.631559506, -4.156271756,
            -3.700743403, -3.260320732, -2.831680453, -2.412317705, -2.000258549, -1.59388586,
            -1.191826998, -0.792876977, -0.395942736, 0, 0.395942736, 0.792876977, 1.191826998,
            1.59388586, 2.000258549, 2.412317705, 2.831680453, 3.260320732, 3.700743403,
            4.156271756, 4.631559506, 5.133595577, 5.673961445, 6.275078705, 6.995680124)

    wz <- c(0.829310817, 0.644938481, 0.565491089, 0.518694458, 0.487223526, 0.46448379,
            0.447333229, 0.434058004, 0.423635472, 0.415416223, 0.408969796, 0.404003106,
            0.400314539, 0.397766974, 0.396271629, 0.395778556, 0.396271629, 0.397766974,
            0.400314539, 0.404003106, 0.408969796, 0.415416223, 0.423635472, 0.434058005,
            0.447333229, 0.464483791, 0.487223526, 0.518694459, 0.565491089, 0.644938482, 0.829310817)

    #mu, sigma2 for -Inf case

    if (!is.na(y)) {
        mu_hat <- y
        sigma2_hat <- techvar2
    } else {
        mu_hat <- mu
        sigma2_hat <- sigma2
    }


    u  <-  z * sqrt(2) * sigma2_hat + mu_hat

    log_lambda <- log(sqrt(2) * sigma2_hat * wz) + stest$log_pyu(y, u, ab, sigma2_hat)

    #log_lambda <- log(sqrt(2) * sigma2_hat * wz) + log_pyu_v2(mu_hat, u, ab, sigma2_hat); # effectively imputation

    list(log_lambda = log_lambda , u = u)
}


stest$estimate_cv <- function(dat, draw = FALSE) {

    # dat = normalized, nonlog data
    thres  <- ncol(dat) / 2

    m <- NULL
    s <- NULL

    for (i in 1:nrow(dat)) {
        x  <- dat[i, ]
        x <- x[!is.na(x)]
        if (length(x) > 0 && sum(x > 1e-8) > thres) {
            x <- x[x > 1e-8]
            m <- c(m, mean(x))
            n <- length(x)
            s <- c(s, sqrt((n - 1) * var(x)/n))
        }
    }
    cc <- median(s / m)

    # robust estimation
    cc_robust <- as.numeric(quantile(s / m, 0.25))

    if (draw) {
        plot(log2(m), log2(s), pch = 20, col = "gray", main = "mean variation plot", xlab = "log2 mean", ylab = "log2 standard deviation",
             xlim = c(min(min(log2(m)), min(log2(s)))-1, max(max(log2(m)), max(log2(s))) + 1), ylim = c(min(min(log2(m)), min(log2(s)))-1, max(max(log2(m)), max(log2(s)) + 1)))

        abline(a = 0, b = 1, col = "red", lwd = 2)

        #cv <- sqrt(exp(cc)-1)
        #c <- -log2(cv)

        xl <- seq(min(min(log2(m)), min(log2(s))) - 2, max(max(log2(m)), max(log2(s))) + 2, length.out = 1000)
        points(xl, xl + log2(cc), type = 'l', lwd = 3, col = "green")
        points(xl, xl + log2(cc_robust), type = 'l', lwd = 3, col = "firebrick")

        points(xl, xl + log2(0.1), type = 'l', lwd = 3, col = "blue", lty=3)
        points(xl, xl + log2(0.15), type = 'l', lwd = 3, col = "seagreen", lty=3)
        points(xl, xl + log2(0.3), type = 'l', lwd = 3, col = "gold", lty=3)

        legend("topleft", legend = c("median CV", paste0("robust CV = ", sprintf("%0.2f", cc_robust)), "CV=10%", "CV=15%", "CV=30%"),
               lty = c(1, 1, 3, 3, 3),
               lwd = c(3, 3, 3, 3, 3),
               col = c("green", "firebrick", "blue", "seagreen", "gold"))
    }

    list(cc = cc_robust, m = m, s = s)
}


stest$estimate_ab <- function(dat, draw = FALSE) {

    # dat = normalized, nonlog data
    dd <- dat
    dd[dd < 1e-8] <- NA

    N  <- nrow(dd)

    # calculate row means in log space
    m  <-  rep(1, N)

    for (i in 1:N) {
        x  <- dd[i,]
        x <- x[!is.na(x)]
        if (length(x) > 0) {
            m[i] <- mean(log2(x))
        } else {
            m[i] <- NA
        }
    }

    x <- seq(from = min(m, na.rm  = TRUE), to = max(m, na.rm  = TRUE), length.out = 21)

    n <- hist(m, breaks = x, plot = FALSE)



    aa <- NULL
    bb <- NULL

    to_plots <- list()
    max_y <- -Inf
    min_y <- Inf


    cc <- 0

    for (i in 1:ncol(dd)) {

        n1 <- hist(m[is.na(dd[, i])], breaks = x, plot = FALSE)

        ratios  <-  n1$counts / n$counts

        ok  <- (ratios < 1) & (ratios > 0) & !is.na(ratios)

        if (sum(ok) == 0) {
            cat(i, ": no missing data, perhaps this is not the best test\n")
        } else {

            p <- lsfit(x = n$mids[ok], y =  log(ratios[ok] / (1 - ratios[ok])))

            aa <- c(aa, p$coefficients[2])
            bb <- c(bb, p$coefficients[1])

            cc <- cc + 1
            to_plots[[cc]] <- list(x = n$mids[ok], y =  log(ratios[ok] / (1 - ratios[ok])))
            max_y <- max(max_y, max(to_plots[[i]]$y))
            min_y <- min(min_y, min(to_plots[[i]]$y))
        }

    }

    if (cc == 0) {
        return(list(a = 0, b = 0))
    }

    # for plot
    xx <- seq(from = (min(m, na.rm  = TRUE)-1), to = (max(m, na.rm  = TRUE)+1), length.out = 1000)
    if (draw) {
        plot(NULL, xlim=c(min(xx), max(xx)), ylim=c(min_y - 0.1, max_y + 0.1), ylab="logistics transformed of misdetection rate", xlab="log2 mean")
        for (i in 1:cc) {
            points(to_plots[[i]]$x, to_plots[[i]]$y, type = 'l', col = palette()[i %% 8 + 1])
        }

        plot(NULL, xlim=c(min(xx), max(xx)), ylim=c(0, 1), ylab="phi_u", xlab="log2 mean")
        for (i in 1:cc) {
            tmp <- exp(xx * aa[i] + bb[i])
            points(xx, tmp/(1.0+tmp), type = 'l', col = palette()[i %% 8 + 1], lty = 3, lwd = 2)
        }

        tmp <- exp(xx * median(aa) + median(bb))
        points(xx, tmp/(1.0+tmp), type = 'l', col = "steelblue", lwd = 3)

        legend("topright", legend = c("Median detection rate"),
               lty = c(1),
               lwd = c(3),
               col = c("steelblue"))

    }


    list( a = median(aa), b = median(bb))
}


stest$lse <- function(a) {

    #y = column max
    y <- apply(a, 2, max)

    # substract the max
    a  <-  a - matrix(y, nrow = nrow(a), ncol = ncol(a), byrow = TRUE)

    return(y + log(colSums(exp(a))))

}


imputation_log2 <- function(log2_data){

    ret <- list()
    ret$original <- log2_data
    ret$seed <- 1203

    set.seed(ret$seed)

    ret$normal_imputed_data <- log2_data

    ret$halfmin_imputation <- ret$normal_imputed_data

    ret$global_mean <- global_mean <- min(ret$normal_imputed_data, na.rm = TRUE)
    ret$global_sd <- global_sd <- mean( apply(ret$normal_imputed_data, 1, sd, na.rm = TRUE), na.rm = TRUE)

    ret$normal_imputed_data[is.na(ret$normal_imputed_data)] <- rnorm(length(ret$normal_imputed_data[is.na(ret$normal_imputed_data)]), mean = global_mean, sd = global_sd)

    ret$halfmin_imputation[is.na(ret$halfmin_imputation)] <- 0.5 * global_mean

    cat("global mean =", global_mean,"; sd = ", global_sd, "\n")
    return(ret)
}


tp_limma_2g <- function(dat, group1, group2) {

    require(limma)
    require(Biobase)

    N <- nrow(dat)


    myeset <- ExpressionSet(assayData = as.matrix(dat)) # use all data

    groups <- colnames(dat)
    for (i in 1:length(groups)) {
        if (groups[i] %in% group1) {
            groups[i] <- "g1"
        } else {
            if (groups[i] %in% group2) {
                groups[i] <- "g2"
            } else {
                groups[i] <- "other"
            }
        }
    }

    if (sum(groups == "other") == 0) {
        groups <- factor(groups, levels = c("g1", "g2"))
        design <- model.matrix(~ 0 + groups)
        colnames(design) <- c("g1", "g2")
    } else {
        groups <- factor(groups, levels = c("g1", "g2", "other"))
        design <- model.matrix(~ 0 + groups)
        colnames(design) <- c("g1", "g2", "other")
    }

    contrast.matrix <- makeContrasts("g2-g1", levels = design)

    fit <- lmFit(myeset, design)

    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)

    a <- topTable(fit2, sort="none",n=Inf)

    return (list(dd = dat[, c(group1,  group2)],
                 logFC = a[, "logFC"],
                 pval = a[, "P.Value"],
                 pval.BH = a[, "adj.P.Val"],
                 design = design))
}
