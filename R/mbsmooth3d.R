
mb.setup.eigvals <- function(s, n)
{
    if (length(n) == 2) {
        Lambda <- kronecker( X = (2 - 2 * cos((0:(n[1] - 1)) * pi / n[1]))*s[1],
                            Y = (2 - 2 * t(cos((0:(n[2] - 1)) * pi / n[2])))*s[2], 
                            FUN = "+")
    }
    else {
        if (length(n) == 3) {
            Lambda <- array(kronecker( X = (2 - 2 * cos((0:(n[1] - 1)) * pi / n[1]))*s[1],
                                      Y = (2 - 2 * t(cos((0:(n[2] - 1)) * pi / n[2])))*s[2], 
                                      FUN = "+"), dim = n) +
                array(s[3]*rep(2 - 2 * cos((0:(n[3] - 1)) * pi / n[3]), each = n[1]*n[2]), dim = n)
        }
    }
    Lambda
}

neg.llhd.general <- function(s, DCT3y) {
    if ((min(s) <= 1e-5) | (max(s[-1]) > 2))
        Inf     else {
                    n <- dim(DCT3y)
                    Lambda <- mb.setup.eigvals(s[-1], n)
                    Gamma <- s[1] + ifelse(Lambda == 0, 0, 1/Lambda)
                    sum(DCT3y^2/Gamma + log(Gamma))
                }
}

mb.smooth3d.general <- function(y, initval)
{
    if ((length(dim(y))+1) != length(initval))
        cat("Error: use of this function only makes sense for higher dimensions.\n")
    else {
        n <- dim(y)
        dct3y <- DCT3(y, inverse = FALSE)
        
        par.val <- optim(par = initval, fn = neg.llhd.general, DCT3y = dct3y)
        
        shat <- par.val$par

        lambda <- mb.setup.eigvals(shat[-1], n)
        gamma <- 1/(shat[1]*lambda + 1)
        c(list(im.smooth = DCT3(gamma * dct3y, inverse = TRUE), par.val = par.val))
    }
}

smooth3d.general.mb <- function(y, shat) {
    n <- dim(y)
    dct3y <- DCT3(y, inverse = FALSE)
    lambda <- mb.setup.eigvals(shat[-1], n)
    gamma <- 1/(1 + shat[1]*lambda)
    c(list(im.smooth = DCT3(gamma * dct3y, inverse = TRUE), par.val = shat))
}


