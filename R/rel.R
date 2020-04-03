#' Read Binary Symmetric Matrix (BSM)
#'
#' Read BSM represented by a pair of files suffixed by ".bin" and ".id",
#' usually produced by PLINK and GCTA.
#'
#' The ".bin" is a binary file storing the matrix entries, which can be
#'
#' \itemize{
#'   \item{the N x N symmetric matrix in full}
#'   \item{the lower triangle with diagonal}
#'   \item{the lower triangle w/o diagonal}
#' },
#' saved as either single or double precision.
#'
#' The ".id"  a text  file of  family ID (FID)  and individual  ID (IID)  in two
#' columns. by default, IID is used as matix row and column names.
#'
#' PLINK option  \code{--make-red bin},  \code{--distance bin}, and  GCTA option
#' \code{--make-grm} all creats binary symmetric matrices, widely used in linear
#' mixed model or kernel based models for genetics.
#'
#' @param pfx prefix of data files {pfx}.id and {pfx}.bin
#' @param dgv diagonal value for matrix without a diagonal (def=1.0)
#' @param fid separator between FID and IID (def=NULL, use IID only)
#' @param bin use bin file instead of the default \code{{pfx}.bin}
#' @param id use id file instead of the default \code{{pfx}.id}
#'
#' @return  symmetric matrix  loaded from file,  with sample ID  in the  row and
#'     column names.
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20.rel")
#' (readBSM(pfx, fid=":"))
#' 
#' @export
readBSM <- function(pfx, dgv=1, fid=NULL, id=NULL, bin=NULL)
{
    ## ID in {p}.id
    if(is.null(id))
        id <- paste0(pfx, ".id")
    ids <- matrix(scan(id, quiet=TRUE, comment.char = "#"), 2)
    if(is.null(fid))
        ids <- ids[1, ]
    else
        ids <- paste(ids[1, ], ids[2, ], sep=fid)
    N <- length(ids)

    ## data matrix
    if(is.null(bin))
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
            R[upper.tri(R, 1)] <- readBin(bin, 0., L, U)
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
            R[upper.tri(R, 0)] <- readBin(bin, 0., L, U)
            R[lower.tri(R, 0)] <- t(R)[lower.tri(R, 0)]
            diag(R) <- dgv # diagnal
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
            R[ , ] <- readBin(bin, 0., L, U)
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
#' Exam the row name for family and individual id.
#'
#' For matrices without rowname, id are automatically generated.
#' 
#' By common practice, the  row names or a matrix are in  the form of [FID.]IID.
#' Samples without family ID are given one identical to their individual ID.
#'
#' @param x  matrix
#' @param sep separator between FID and IID forming the sample ID
#' @return data.frame of inferred family ID and individual ID.
gid <- function(x, sep=".")
{
    N <- nrow(x)
    i <- row.names(x)
    if(is.null(i))
    {
        i <- sprintf(paste0('I%d0', nchar(N)), seq(N))
    }
    sep <- paste0('[', sep, ']')
    i <- lapply(strsplit(i, sep), rep, l=2)
    i <- data.frame(do.call(rbind, i), stringsAsFactors=FALSE)
    names(i) <- c('FID', 'IID')
    i
}

#' Save Symmetric Matrix to Binary
#'
#' Save symmetric matrix  to a binary core  file (.bin), and a text  file of IDs
#' (.id), recognizable by PLINK.
#' 
#' @param pfx prefix of output files
#' @param x symmetric matrix to save
#' @param ltr store the lower triangle only? (def=TRUE)
#' @param diag save diagnal? (def=TRUE) ignored if \code{ltr} is FALSE.
#' @param unit numerical unit, (def=4, single precision)
#' @param fid separator between FID and IID (def=".").
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20.rel")
#' rel <- readBSM(pfx)  # relatedness kernel matrix
#' re2 <- rel^2         # 2nd order polynomial kernel
#'
#' tmp <- tempdir()
#' dir.create(tmp, FALSE)
#' out <- file.path(tmp, 'm20.re2')
#' saveBSM(out, re2)    # save the polynomial kernel
#' dir(tmp)             # show new files, then clean up
#' unlink(tmp, recursive=TRUE)
#' 
#' @export
saveBSM <- function(pfx, x, ltr=TRUE, diag=TRUE, unit=4L, fid=".")
{
    ## save subject IDs
    id <- gid(x, sep=fid)
    write(t(id), paste0(pfx, '.id'), 2L, sep='\t')
    
    ## R is colume majored, to save the lower TRI, use upper TRI instead.
    if(length(ltr) == 0 || ltr)
        x <- x[upper.tri(x, diag)]
        
    ## save matrix data
    writeBin(x, paste0(pfx, '.bin'), unit)
}


#' Read PLINK Binary IBS matrix
#'
#' A PLINK IBS (Identity by State) matrix is represented by
#' \itemize{
#'   \item {.mibs.bin:}{IBS matrix in binary}
#'   \item {.mibs.id :}{FID and IID in text}
#' }
#' A binary IBS matrix is the result of PLINK \code{--distance ibs bin}
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' (readIBS(pfx))
#' 
#' @param pfx prefix of the IBS file set.
#' @param fid seperate after family ID (def=NULL, use IID only)
#' @return IBS matrix with row and column names set to sample ID.
#' @export
readIBS <- function(pfx, fid=".") readBSM(paste0(pfx, ".mibs"), fid=fid)


#' Read PLINK Binary REL matrix
#'
#' A PLINK REL (Relatedness) matrix is represented by
#' \itemize{
#'   \item {.rel.bin:}{REL matrix in binary}
#'   \item {.rel.id :}{FID and IID in text}
#' }
#' A binary REL matrix is the result of PLINK \code{--make-rel bin}
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' (readREL(pfx))
#' 
#' @param pfx prefix of the REL file set
#' @param fid separate after family ID. (def=NULL, use IID only)
#' @return relatedness matrix with row and column names set to sample ID.
#' @export
readREL <- function(pfx, fid=".") readBSM(paste0(pfx, ".rel"), fid=fid)

#' Read Genetic Related Matrix (GRM) of GCTA
#'
#' GRM is the  core formt of GCTA,  which is an binary symmetric  matrix with an
#' extra variant count  matrix (VCM), this function reads  the binary sysmmetric
#' matrix.
#'
#' GCTA GRM is represented by a set of three files:
#' 
#' \itemize{
#'   \item {.grm.bin   :}{GRM matrix in binary}
#'   \item {.grm.id    :}{sample FID and IID in text}
#'   \item {.grm.N.bin :}{number of valid variants for each GRM entry}}
#' 
#' and it always  uses single precision (4  bytes per entry).
#'
#' To read the extra the extra VCM (grm.N.bin), use \code{\link{readVCM}}.
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' (readGRM(pfx))
#'
#' @param pfx prefix of GRM file set
#' @param fid separator after family ID (def=NULL, use IID only)
#' @return matrix of relatedness with sample ID in row and column names.
#' @export
readGRM <- function(pfx, fid=".") readBSM(paste0(pfx, ".grm"), fid=fid)

#' Read Variant Count Matrix (VCM) accompanying a GCTA GRM
#'
#' GRM (Genetic Relatedness Matrix) is the core formt of GCTA, which is a PLINK
#' binary  symmetric matrix  with  an  extra variant  count  matrix (VCM),  this
#' function reads the VCM.
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' (readVCM(pfx))
#'
#' @param pfx prefix of GRM file set
#' @param fid seperate after family ID (def=NULL, use IID only)
#' @return matrix of variant count with sample ID in row and column names.
#' @export
readVCM <- function(pfx, fid=NULL) readBSM(paste0(pfx, ".grm.N"), fid=fid, id=paste0(pfx, ".grm.id"))

#' Save symmetic matrix to GCTA GRM format.
#'
#' GRM (Genetic  Relatedness Matrix) is  the core  formt of GCTA,  this function
#' saves a R symmetric matrix to a file set recgnizable by GCTA.
#'
#' Three files will be saved:
#' 
#' \itemize{
#'   \item {.grm.bin   :}{genetic relatedness matrix in binary}
#'   \item {.grm.id    :}{FID and IID for N individuals in text}
#'   \item {.grm.N.bin :}{variant count matrix (VCM) in binary}
#' }
#'
#' FID and IID will be generated if the \code{grm} to be saved has no row names.
#' 
#' When save the  \code{vcm}, if a single  number is given, this  number is used
#' as the variant count for all entries in the GRM.
#'
#' \code{saveGRM} is useful in exporting  customized kinship matrices (such as a
#' Gaussian or a  Laplacian kernel) to a  GRM acceptable by GCTA,  which are not
#' supported by GCTA's own GRM builder.
#'
#' @param pfx prefix of data files
#' @param grm genome relatedness matrix to save
#' @param vcm variant counts matrix to save (def=1).
#' @param fid separator after family ID. (def=".")
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' gmx <- readBED(pfx)  # read genotype matrix from PLINK BED.
#' gmx <- scale(gmx)    # standardize
#' tmp <- tempdir()     # for example outputs
#' dir.create(tmp, FALSE)
#' 
#' # kinship matrix as Gaussian kernel, built from the first 10 variants
#' gmx.gau <- gmx[, +(1:10)]                 # the first 10 variants
#' not.na.gau <- tcrossprod(!is.na(gmx.gau)) # variant count matrix
#' kin.gau <- exp(as.matrix(-dist(gmx.gau, "euc")) / not.na.gau)
#' print(kin.gau)                            # the Gaussian kernel
#' out.gau <- file.path(tmp, "m20.gau")
#' saveGRM(out.gau, kin.gau, not.na.gau)     # gau.grm.* should appear
#'
#' # kinship matrix as Laplacian kernel, built without the first 10 variants
#' gmx.lap <- gmx[, -(1:10)]                 # drop the first 10 variants
#' not.na.lap <- tcrossprod(!is.na(gmx.lap)) # variant count matrix
#' kin.lap <- exp(as.matrix(-dist(gmx.lap, "man")) / not.na.lap)
#' out.lap <- file.path(tmp, "m20.lap")
#' print(kin.lap)                            # the Laplacian kernel
#' saveGRM(out.lap, kin.lap, not.na.lap)     # lap.grm.* should appear
#'
#' # merge kinship in R language for a radius based function kernel matrix
#' not.na.rbf <- not.na.gau + not.na.lap
#' kin.rbf <- (kin.gau * not.na.gau + kin.lap * not.na.lap) / not.na.rbf
#' print(kin.rbf)
#' out.rbf <- file.path(tmp, "m20.rbf")
#' saveGRM(out.rbf, kin.rbf, not.na.rbf)     # rbf.grm.* should appear
#'
#' # show saved matrices, then clean up
#' dir(tmp, "(gau|lap|rbf)")
#' unlink(tmp, recursive=TRUE)
#' 
#' @export
saveGRM <- function(pfx, grm, vcm=NULL, fid=".")
{
    ## get file names
    fn.rmx <- paste0(pfx, ".grm.bin")
    fn.N <- paste0(pfx, ".grm.N.bin")
    fn.id <- paste0(pfx, ".grm.id")

    ## complete id and N
    id <- gid(grm, sep=fid)

    if(is.null(vcm))
        vcm <- 1
    if(length(vcm) != length(grm) && length(vcm) != 1)
        stop("wrong dimension for variant count matrix.")

    ## upper.tri of col major = lower.tri of row major
    u <- upper.tri(diag(nrow(id)), TRUE)
    grm <- grm[u]
    writeBin(grm, fn.rmx, 4L)

    ## variant count matrix
    if(length(vcm) > 1)
        vcm <- vcm[u]
    writeBin(vcm, fn.N, 4L)

    ## individual ID table
    write(t(id), fn.id, 2, sep='\t')

    invisible(NULL)
}

#' Test Genetic Relatedness Matrix Reader
#'
#' Compare the read from genetic relatedness matrix created from the same genome
#' segment but stored in different shapes and types.
testReadBSM <- function()
{
    one <- function(pfx)
    {
        fdr <- file.path(system.file("extdata", package="plinkFile"), pfx)
        
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
