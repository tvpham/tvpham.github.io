
estimate_cv <- function(dat) {

    # dat = normalized, nonlog data
    thres  <- ncol(dat) / 2

    m <- NULL
    s <- NULL

    for (i in 1:nrow(dat)) {
        x  <- dat[i, ]
        if (sum(x > 0) > thres) {
            x  <- x[x > 0];
            m  <- c(m, mean(x));
            n <- length(x)
            s  <- c(s, sqrt((n - 1) * var(x)/n));
        }
    }
    cc  <-  median(s / m);

    list(cc = cc, m = m, s = s)
}


estimate_ab <- function(dat) {

    # dat = normalized, nonlog data

    N  <- nrow(dat)

    # calculate row means in log space
    m  <-  rep(1, N)

    for (i in 1:N) {
        x  <- dat[i,];
        x  <- x[x > 0];
        m[i] <- mean(log(x))
    }

    x <- seq(from = min(m, na.rm  = TRUE), to = max(m, na.rm  = TRUE), length.out = 21)

    n <- hist(m, breaks = x, plot = FALSE)


    aa <- NULL
    bb <- NULL

    for (i in 1:ncol(dat)) {

        n1 <- hist(m[dat[, i] == 0], breaks = x, plot = FALSE)

        ratios  <-  n1$counts / n$counts;
        ok  <- (ratios < 1) & (ratios > 0)

        p <- lsfit(x = n$mids[ok], y =  log(ratios[ok] / (1 - ratios[ok])))

        aa <- c(aa, p$coefficients[2])
        bb <- c(bb, p$coefficients[1])
    }

    list( a = median(aa), b = median(bb))
}


estimate_mean_variance_in_log_space <- function(ldat) {

    # dat = normalized, nonlog data

    thres <- ncol(ldat) / 2

    N <- nrow(ldat);

    # calculate row means in log space
    l_m  <- NULL
    l_s2 <- NULL

    for (i in 1:N) {

        x  <- as.numeric(ldat[i, ])

        if (any(is.finite(x))) {

            x  <- x[is.finite(x)]; # detected entries

            n <- length(x)

            if (n > thres) {
                l_m <- c(l_m, mean(x))
                l_s2 <- c(l_s2, (n - 1) * var(x) / n)
            }
        }
    }

    list(l_m = l_m, l_s2 = l_s2)
}


weighted_sampling_pyu <- function(y, ab, techvar2, mu, sigma2) {

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

    if (is.finite(y)) {
        mu_hat <- y
        sigma2_hat <- techvar2
    } else {
        mu_hat <- mu
        sigma2_hat <- sigma2
    }

    u  <-  z * sqrt(2) * sigma2_hat + mu_hat;

    log_lambda <- log(sqrt(2) * sigma2_hat * wz) + log_pyu(y, u, ab, techvar2);

    list(log_lambda = log_lambda , u = u)
}


log_pyu <- function(y, u, ab, c2) {

    phi <- exp(u * ab$a + ab$b);

    if (is.finite(y)) {
        ll <- -log(1 + phi) - log(sqrt(2.0 * pi * c2)) - (y - u)*(y - u) / (2.0 * c2)
    } else {
        ll <- log(phi / (1.0 + phi))
    }

    return(ll)
}


estimate_var_per_group <- function(y, techvar2, df0) {


    y_ok  <- y[is.finite(y)]

    if (all(is.infinite(y_ok)) || length(y_ok) < 2) {
        var2 <- techvar2
        df <- df0;
    } else {
        n <- length(y_ok)
        var2  <-  (n-1) * var(y_ok) / n

        if (n > 1 && var2 < techvar2) {
            var2 <- techvar2
            df <- df0
        } else {
            df <- length(y_ok)
        }
    }

    list(var2 = var2, df = df)
}


lse <- function(a) {

    #y = column max
    y <- apply(a, 2, max)

    # substract the max
    a  <-  a - matrix(y, nrow = nrow(a), ncol = ncol(a), byrow = TRUE)

    return(y + log(colSums(exp(a))))

}


max_LL <- function(Lw, U, X, y, s2) {

    maxT <- 1e4;
    f_tol <- 1e-8;
    beta_tol <- 1e-8;

    M <- nrow(U)
    K <- ncol(U)

    A  <- solve(X %*% t(X), X)

    yy  <-  y;
    yy[is.infinite(yy)]  <-  mean(yy[is.finite(y)]); # for starting values

    beta  <-  A %*% as.matrix(yy)

    xbeta  <-  t(t(X) %*% beta)

    L  <-  -Inf;

    for (t in 1:maxT) {

        old_L  <-  L
        old_beta  <-  beta

        R  <- U - matrix(xbeta, nrow = M, ncol = length(xbeta), byrow = TRUE)
        R2  <- R * R

        H_mk  <-  Lw - 0.5 * log(2 * pi * s2) - R2 / (2 * s2);

        Lk  <-  lse(H_mk)

        L  <-  sum(Lk)

        W_mk  <- exp(H_mk - matrix(Lk, nrow = M, ncol = length(Lk), byrow = TRUE))

        v  <- colSums(U * W_mk)

        beta  <- A %*% as.matrix(v)
        xbeta  <- t(t(X) %*% beta)

        if (max(abs(old_L - L)) < f_tol && max(abs(old_beta - beta)) < beta_tol) {
            break;
        }
    }

    list(L = L, beta = beta)
}


s_test <- function(y1, y2, techvar2, df0, ab) {

    y <- c(y1, y2)
    X <- rbind(rep(1, length(y)),
               c(rep(1, length(y1)), rep(2, length(y2))))

    y_ok <- y[is.finite(y)]

    if (all(is.infinite(y))) {
        return(1.0)
    }

    null_mu <- min(y_ok)
    null_sigma2 <- techvar2

    # ----- simulate
    Lw <- NULL
    U <- NULL

    for (k in 1:length(y)) {
        s <-  weighted_sampling_pyu(y[k], ab, techvar2, null_mu, null_sigma2)
        Lw  <- cbind(Lw, s$log_lambda)
        U  <- cbind(U, s$u)
    }


    #% ----- unrestricted MLE

    vd1 <- estimate_var_per_group(y1, techvar2, df0)
    vd2 <- estimate_var_per_group(y2, techvar2, df0)

    sigma2_hat = (vd1$df * vd1$var2 + vd2$df * vd2$var2) / (vd1$df + vd2$df)

    uLogL <-  max_LL(Lw, U, X, y, sigma2_hat)

    #% ----- restrictied MLE
    rLogL <-  max_LL(Lw, U, matrix(1, nrow =1, ncol = length(y)) , y, sigma2_hat);

    #%[h, pValue, stat] = lratiotest(uLogL, rLogL, 1);

    stat  <-  2 * (uLogL$L - rLogL$L)

    if (stat < 0)
        stat  <-  0;
    end

    #stat

    pval <- 1 - pchisq(stat, 1)

    return(pval)

}


# iterface
s.test <- function(d, g1, g2) {

    cc <- estimate_cv(d)
    ab <- estimate_ab(d)

    tvar2 <- log1p(cc$cc * cc$cc)

    ldat <- log(d[, c(g1, g2)])

    v  <- estimate_mean_variance_in_log_space(ldat)

    mean(v$l_s2)
    precision = 1.0 / v$l_s2

    require(MASS)

    shape_rate  <- fitdistr(precision, "gamma")

    alpha <- as.numeric(shape_rate$estimate[1])

    df0  <-  2 * alpha

    N <- nrow(d)
    pval  <- rep(1, N)

    ind1 <- 1:length(g1)
    ind2 <- (length(g1) + 1):(length(g1) + length(g2))

    for (i in 1:N) {
        y  <- as.numeric(ldat[i,])
        pval[i] <- s_test(as.numeric(y[ind1]), as.numeric(y[ind2]), tvar2, df0, ab)
    }

    return(pval)
}
