.name <- function(len, pfx=NULL, max=NULL, sfx=NULL, alpha=NULL)
{
    if(is.null(alpha))
        alpha <- c(0:9, LETTERS)
    b <- length(alpha) # base
    if(is.null(max))
        max <- len

    n <- seq(len)
    r <- list()
    q <- 1
    while(q <= max)
    {
        r <- c(list(alpha[n %% b + 1]), r)
        n <- n %/% b
        q <- q * b
    }
    paste0(pfx, do.call(paste0, r), sfx)
}
