#' Line Count
#'
#' Count the occurance of '\\n', much faster than \code{readLine} and
#' \code{read.table}. Although slower than the unix command "wc -l",
#' it upholds platform independency.
#'
#' @param f the file name, or a connection.
#' @noRd
lc <- function(f)
{
    f <- file(f, open="rb")                        
    n <- 0L
    p <- as.raw(10L)                    # '\n'
    while (length(chunk <- readBin(f, "raw", 65536)) > 0)
    {
        n <- n + sum(chunk == p)
    }
    if(isOpen(f))
       close(f)
    n
}


i2b <- function(i, n=1)
{
    n <- max(1, min(n, 4)) * 8
    b <- intToBits(12)[1:n]
    paste(sapply(strsplit(paste(b),""), `[[`, 2), collapse="")
}

r2b <- function(x)
{
    b <- sapply(x, function(.)
    {
        . <- rawToBits(.)
        paste(sapply(strsplit(paste(.),""), `[[`, 2), collapse="")
    })
    attributes(b) <- attributes(x)
    b
}
