#' Line Count
#'
#' Count the occurance of '\\n', much faster than \code{readLine} and
#' \code{read.table}. Although slower than the unix command "wc -l",
#' it upholds platform independency.
#'
#' @param f the file name, or a connection.
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
