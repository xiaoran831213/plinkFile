#' number of rows
#'
#' Count the occurance of '\\n', much faster than \code{readLine} and
#' \code{read.table}. Although slower than the unix command "wc -l",
#' it upholds platform independency.
#'
#' @param f the file name, or a connection.
#' @noRd
nr <- function(f)
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

#' number of fields
#'
#' Count the columns according to the first line.
#'
#' @param f the file name, or a connection
#' @noRd
nf <- function(f)
{
    length(scan(f, "\t", nlines = 1, quiet = TRUE))
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
