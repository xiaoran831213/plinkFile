#' Read a list of ID
#'
#' The typical file subject ID has no header, and two columns, family
#' ID (FID) and within family individual ID (IID).
#'
#' For data exclusively composed of unrelated individuals, FID equals
#' IID in most cases. 
#'
#' Because of the usual practice of assigning uniqe IID across the entire
#' cohort even if the individuals form families, the default behavior is
#' fething only the IID. If both FID and IID are read, the final ID will
#' have the form "FID.IID".
#' 
#' @param fn file name of hte ID list
#' @param only.IID Individuals ID only, ignore Family ID?
#'
#' @return an vector of sample ID
.get.id <- function(fn, only.IID=TRUE)
{
    . <- matrix(scan(fn, "", quiet=TRUE), ncol=2, byrow=TRUE)
    . <- as.data.frame(., stringsAsFactors=FALSE)
    names(.) <- c('FID', 'IID')
    if(only.IID)
        .$IID
    else
        with(., paste(FID, IID, Read='.'))
}

#' sep PLINK Binary Square Matrix
#'
#' With the binary data file, and a list of N subject IDs,
#' this function figure out the storage type of the matrix
#' and retrieve it as a standard R matrix.
#'
#' PLINK users can save the genetic relatedness matrices in
#' various forms:
#'
#' * an N by N square matrix as it is
#' * the lower triangle with the diagonal
#' * the lower triangle without the diagonal
#'
#' Aside from the differences in shape, they can choose
#' between single or double precision numbers, occupying 4
#' or 8 bytes per entry, respectively.
#' 
#' @param fn filename of the binary square matrix
#' @param id vector of N subject ID
#' @param dg diagonal value for lower triangle matrices
#' without diagonal, the default is 1.
#'
#' @return an N by N square matrix (GRM) load from file.
.get.bm <- function(f, id, dg=1)
{
    M <- length(id)
    S <- file.size(f)                  # file size
    R <- matrix(.0, M, M, dimnames=list(id, id))

    bin <- sub("[.].*$", "", basename(f))

    ## try: lower triangle with diagonal
    L <- M * (M + 1.0) / 2.0
    U <- S / L
    if(U == 4 || U == 8)
    {
        R[upper.tri(R, 1)] <- readBin(f, .0, L, U)
        R[lower.tri(R, 0)] <- t(R)[lower.tri(R, 0)]
        print(data.frame(bin=bin, size=S, num.entries=L, unit.bytes=U, shape='LWD'))
        return(R)
    }

    ## try: lower triangle without diagonal
    L <- M * (M - 1.0) / 2.0            # number of entries
    U <- S / L                          # unit size
    if(U == 4 || U == 8)
    {
        R[upper.tri(R, 0)] <- readBin(f, .0, L, U)
        R[lower.tri(R, 0)] <- t(R)[lower.tri(R, 0)]
        diag(R) <- dg                   # assigned diagnal
        print(data.frame(bin=bin, size=S, num.entries=L, unit.bytes=U, shape='LND'))
        return(R)
    }

    ## try: a squre
    L <- 1.0 * M * M
    U <- S / L
    if(U == 4 || U == 8)
    {
        R[ , ] <- readBin(f, .0, L, U)
        print(data.frame(bin=bin, size=S, num.entries=L, unit.bytes=U, shape='SQR'))
        return(R)
    }

    ## fail: return
    stop("can not figure out the shape of recoded entries in the relatedness matrix.")
}

#' Infer Sample ID from a Relatedness Matrix
#'
#' Exam the row name a relatedness matrix for family IDs and
#' individual IDs, which form the sample IDs.
#'
#' The IDs are automatically generate for symmetric matrices
#' without a row name.
#' 
#' By common practice, the row names or a relatedness matrix are
#' in the form of [FID.]IID. Samples without family ID are given
#' one identical to their individual ID.
#'
#' @param rmx relatedness matrix
#' 
#' @return data.frame of inferred family ID and individual ID.
.inf.id <- function(rmx)
{
    N <- nrow(rmx)
    i <- row.names(rmx)
    if(is.null(i))
        i <- sprintf('I%d05', seq(N))

    i <- lapply(strsplit(i, '[.]'), rep, l=2)
    i <- data.frame(do.call(rbind, i), stringsAsFactors=FALSE)
    names(i) <- c('FID', 'IID')
    i
}

#' Write Squre Matrix to PLINK Binary
#'
#' Both the matrix data and subject ID will be written to
#' files with a given prefix.
#' 
#' @param f character prefix of target files
#' @param x matrix to be saved
#' @param l logical only store use lower triangle? def=TRUE
#' @param d logical store diagnal? def=TRUE
#' @param u numerical unit size, def=4 (single precision)
.put.bm <- function(f, x, l=TRUE, d=TRUE, u=4L)
{
    ## save subject IDs
    id <- .inf.id(x)
    write(t(id), paste0(f, '.id'), 2L, sep='\t')
    
    ## R is colume major, to save the lower TRI, one extract
    ## upper TRI instead.
    if(isTRUE(l))
        x <- x[upper.tri(x, d)]
        
    ## save matrix data
    writeBin(x, paste0(f, '.bin'), u)
}


#' Read PLINK Binary IBS matrix
#'
#' A PLINK IBS (Identity by State) matrix is represented by
#'
#' * .mibs.bin   : the binary file for N by N relatedness matrix;
#' * .mibs.id    : list of samples IDs
#'
#' @param pfx shared prefix of the three files
#' @return matrix of relatedness
#' @export
readIBS <- function(pfx)
{
    .get.bm(paste0(pfx, ".mibs.bin"), .get.id(paste0(pfx, ".mibs.id")))
}

#' Read PLINK Binary REL matrix
#'
#' PLINK REL (relatedness) matrix comes with two files:
#'
#' * .grm.bin : an N by N matrix of relatedness;
#' * .grm.id  : a list of sample ID.
#'
#' @param pfx shared prefix of the three files
#' @return a relatedness matrix, with row and column names set to subject IDs.
#' @export
readREL <- function(pfx)
{
    .get.bm(paste0(pfx, ".rel.bin"), .get.id(paste0(pfx, ".rel.id")))
}

saveREL <- function(pfx, grm)
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

#' Read GRM binary of GCTA
#'
#' GRM (genetic relatedness matrix) is the main formt for GCTA
#' developed by Jian Yang et. al. A GRM has the surfix:
#'
#' * .grm.bin : an N by N matrix of relatedness;
#' * .grm.id  : a list of sample ID.
#'
#' which is essentially the same with PLINK REL.
#'
#' A GRM comes with an addtional binary matrix ".grm.N.bin" that counts
#' the genomic variants contributed to each relatedness, and it always
#' adapts single precision and lower triangle with diagnonal. These
#' specifics however, does not affect the behavior of the reader.
#' 
#'
#' @param pfx shared prefix of data files
#' @return N by N matrix of relatedness
#' @export
readGRM <- function(pfx)
{
    .get.bm(paste0(pfx, ".grm.bin"), .get.id(paste0(pfx, ".grm.id")))
}
