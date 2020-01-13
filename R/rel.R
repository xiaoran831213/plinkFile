#' read binary symmetric matrix
#'
#' read binary symmetric  matrix from PLINK files {p}.bin and  {p}.id, where {p}
#' is a shared prefix.
#'
#' {p}.bin stores the matrix, which can be:
#'
#'   * the N x N symmetric matrix itself
#'   * the lower triangle w/t diagonal
#'   * the lower triangle w/o diagonal
#' 
#' saved as either single or double precision.
#'
#' {p}.id contains  family ID (FID)  and individual ID  (IID) in two  columns of
#' text, by default only the 2nd column (IID) is used.
#'
#' @param pfx shared prefix for two PLINK files
#' @param dgn diagonal values for lower triangle without diagonal (def=1.0)
#' @param fid read FID as well (def=FALSE)
#' @param sep separator put between fid and iid to form the final ID
#'
#' @return N by N symmetric matrix loaded from file.
#' @export
readBSM <- function(pfx, dgn=1, fid=FALSE, sep=".")
{
    ## ID in {p}.id
    ids <- matrix(scan(paste0(pfx, ".id"), quiet=TRUE, comment.char = "#"), 2)
    if(isTRUE(fid))
        ids <- paste(ids[1, ], ids[2, ], sep=sep)
    else
        ids <- ids[1, ]
    N <- length(ids)
    
    ## data matrix
    bin <- paste0(pfx, ".bin")
    S <- file.size(bin) # file size
    R <- matrix(.0, N, N, dimnames=list(ids, ids))
    H <- "UNK" # shape, start with unknown

    ## try: lower triangle w/t diagonal
    if(H == "UNK")
    {
        L <- N * (N + 1.0) / 2.0
        U <- S / L
        if(U == 4 || U == 8)
        {
            R[upper.tri(R, 1)] <- readBin(bin, .0, L, U)
            R[lower.tri(R, 0)] <- t(R)[lower.tri(R, 0)]
            H="LWD"
        }
    }

    ## try: lower triangle, no diagonal
    if(H == "UNK")
    {
        L <- N * (N - 1.0) / 2.0 # number of entries
        U <- S / L               # unit size
        if(U == 4 || U == 8)
        {
            R[upper.tri(R, 0)] <- readBin(bin, .0, L, U)
            R[lower.tri(R, 0)] <- t(R)[lower.tri(R, 0)]
            diag(R) <- dgn # diagnal
            H="LND"
        }
    }

    ## try: a squre
    if(H == "UNK")
    {
        L <- 1.0 * N * N
        U <- S / L
        if(U == 4 || U == 8)
        {
            R[ , ] <- readBin(bin, .0, L, U)
            H="SQR"
        }
    }

    ## fail or return
    print(data.frame(pfx=pfx, size=S, entries=L, unit=U, shape=H, N=N))
    if(H == "UNK")
        stop("unknown storage type.")
    R
}

#' Infer Sample ID from a symmetric matrix
#'
#' Exam the row name family and individual id.
#'
#' For matrices without rowname, id are automatically generated.
#' 
#' By common practice, the row names or a relatedness matrix are
#' in the form of [FID.]IID. Samples without family ID are given
#' one identical to their individual ID.
#'
#' @param x  matrix
#' 
#' @return data.frame of inferred family ID and individual ID.
.inf.id <- function(x, sep=".")
{
    N <- nrow(x)
    i <- row.names(x)
    if(is.null(i))
        i <- sprintf('I%d05', seq(N))

    sep <- paste0('[', sep, ']')
    i <- lapply(strsplit(i, sep), rep, l=2)
    i <- data.frame(do.call(rbind, i), stringsAsFactors=FALSE)
    names(i) <- c('FID', 'IID')
    i
}

#' Save Symmetric Matrix to Binary
#'
#' save the matrix to binary {p}.bin, and  the row names text file {p}.id, where
#' {p} is a given file prefix.
#' 
#' @param pfx prefix of output files
#' @param x matrix to be saved
#' @param l store the lower triangle? def=TRUE
#' @param d store the diagnal? def=TRUE
#' @param u numerical unit, def=4 (single precision)
#' @export
saveBSM <- function(pfx, x, l=TRUE, d=TRUE, u=4L)
{
    ## save subject IDs
    id <- .inf.id(x)
    write(t(id), paste0(pfx, '.id'), 2L, sep='\t')
    
    ## R is colume majored, to save the lower TRI, use upper TRI instead.
    if(isTRUE(l))
        x <- x[upper.tri(x, d)]
        
    ## save matrix data
    writeBin(x, paste0(pfx, '.bin'), u)
}


#' Read PLINK Binary IBS matrix
#'
#' A PLINK IBS (Identity by State) matrix is represented by
#'
#'   * {p}.mibs.bin   : IBS matrix in binary
#'   * {p}.mibs.id    : samples FID and IID in text
#'
#' IBS matrix is a result of PLINK --distance ibs
#'
#' @param pfx shared prefix of the three files
#' @return matrix of relatedness
#' @export
readIBS <- function(pfx) readBSM(paste0(pfx, ".mibs"))

#' Read PLINK Binary REL matrix
#'
#' A PLINK REL (Relatedness) matrix is represented by
#'
#'   * {p}.rel.bin   : IBS matrix in binary
#'   * {p}.rel.id    : samples FID and IID in text
#'
#' REL matrix is a result of PLINK --distance ibs
#'
#' @param pfx shared prefix of the three files
#' @return a relatedness matrix, with row and column names set to subject IDs.
#' @export
readREL <- function(pfx) readBSM(paste0(pfx, ".rel"))

#' Read GRM binary of GCTA
#'
#' GRM (genetic relatedness matrix) is the main formt of GCTA, basically a PLINK
#' REL with differnent surfix:
#'
#'   * .grm.bin : binary matrix,
#'   * .grm.id  : samples,
#'
#' and it always uses single precision.
#' 
#' GRM  also comes  with another  binary matrix  ".grm.N.bin" counting  variants
#' contributed to each relatedness.
#'
#' @param pfx shared prefix of data files
#' @return N by N matrix of relatedness
#' @export
readGRM <- function(pfx) readBSM(paste0(pfx, ".grm"))

saveGRM <- function(pfx, grm)
{
    ## get file names
    fn.rmx <- paste0(pfx, ".grm.bin")
    fn.N <- paste0(pfx, ".grm.N.bin")
    fn.id <- paste0(pfx, ".grm.id")

    ## complete id and N
    if(is.matrix(grm))
    {
        grm <- list(rmx=grm, id=.inf.id(grm), N=1.0)
    }

    with(grm,
    {
        ## upper.tri of col major = lower.tri of row major
        idx <- upper.tri(diag(nrow(id)), T)
        
        ## genomic relatedness matrix
        rmx <- rmx[idx]
        writeBin(rmx, fn.rmx, 4L)

        ## genomic variant count matrix
        N <- N[idx]
        writeBin(N, fn.N, 4L)

        ## subject IDs
        write(t(id), fn.id, 2, sep='\t')
    })
}

#' Test Genetic Relatedness Matrix Reader
#'
#' Compare the read from genetic relatedness matrix created
#' from the same genome segment but stored in different shapes
#' and types.
#' @export
relTest <- function()
{
    one <- function(pfx)
    {
        fdr <- file.path(system.file("extdata", package="plinkBED"), pfx)
        
        ## standardized relatedness, lower trangle with diagonal, single
        st4 <- readREL(paste0(fdr, ".st4"))
        ## standardized relatedness, lower trangle with diagonal, double
        st8 <- readREL(paste0(fdr, ".st8"))
        ## standardized relatedness, square matrix, single
        ss4 <- readREL(paste0(fdr, ".ss4"))
        ## standardized relatedness, square matrix, double
        ss8 <- readREL(paste0(fdr, ".ss8"))
        ## GCTA genetic relatedness matrix
        grm <- readGRM(fdr)
        if(!all(st4 == ss4, st8 == ss8))
            stop("Failed REL reading test:", pfx, ", full-square vs low-triangle.")
        if(!all.equal(st4, st8, tolerance=1e-7))
            stop("Failed REL reading test:", pfx, ", single vs double precision.")
        if(!all.equal(ss4, ss8, tolerance=1e-7))
            stop("Failed REL reading test:", pfx, ", single vs double precision.")
        if(!all(st4 == grm, ss4 == grm))
            stop("Failed GRM reading test:", pfx, ", single vs double precision.")

        cat("Pass PLINK relatedness matrix reading test: ", pfx, "\n", sep="")
    }
    one('m20')
}
