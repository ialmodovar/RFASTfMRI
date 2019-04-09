biweight.scale.est <- function(x, center = median(x), w = 6){
    y <- x-center
    n <- length(y)
    stilde <- mad(x)
    u <- y/w/stilde
    num <- sum(y^2*(1-u^2)^4*(u^2 < 1))
    den <- sum((1-u^2)*(1-5*u^2)*(u^2 < 1))
    sqrt(n * num)/abs(den)
}
