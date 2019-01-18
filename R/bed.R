#' Decompress Byte Data
#'
#' For each SNP (i.e., a logical row in the BED file), a byte in the BED
#' file can encodes the genotype of up to 4 consecutive samples (2 bits
#' each).
#'
#' The function decondens the bytes read from a BED into one byte per sample
#' bases, then decode the bytes into alternative allele dosage counts or NA.
#'
#' @param bd byte data in R "raw" mode
#' "integer", "numeric", "double", and "character".
#' @return a vector of dosage data in desired type, untruncated.
db1 <- function(bd)
{
    ## fill the genomic data vector
    .m <- as.raw(0x03)                       # 0x03 = 11000000

    cat("converting ...\n")
    cat("sample 1, 5,  9, ...\n")
    g1 <- as.integer(bd & .m)                # sample 1, 5, ...
    cat("sample 2, 6, 10, ...\n")
    g2 <- as.integer(rawShift(bd, -2L) & .m) # sample 2, 6, ...
    cat("sample 3, 7, 11, ...\n")
    g3 <- as.integer(rawShift(bd, -4L) & .m) # sample 3, 7, ...
    cat("sample 4, 8, 12, ...\n")
    g4 <- as.integer(rawShift(bd, -6L) & .m) # sample 4, 8, ...
    
    ## decode to dosage
    cat("decoding ...\n")
    ## 00->2: homozygous A1
    ## 01->3: missing
    ## 10->1: heterozygous
    ## 11->0: homozygous A2
    dc <- c(2L, NA, 1L, 0L)
    cat("sample 1, 5,  9, ...\n")
    g1 <- dc[g1 + 1L]
    cat("sample 2, 6, 10, ...\n")
    g2 <- dc[g2 + 1L]
    cat("sample 3, 7, 11, ...\n")
    g3 <- dc[g3 + 1L]
    cat("sample 4, 8, 12, ...\n")
    g4 <- dc[g4 + 1L]

    cat("matrix packing ...\n")
    rt <- as.vector(rbind(g1, g2, g3, g4), "integer")
    rt
}

#' Read the BED file of a set of PLINK Triplets
#'
#' Read genotype matrix from a PLINK BED file set. This function is meant
#' for in-memory examination of small PLINK BEDs. For out-of-memory computation
#' tasks involving giagantic BEDs, see \code{\link{scanBED}} instead.
#' 
#' PLINK BED (\emph{binary biallelic genotype table}) is comprised of three files
#' (usually) sharing identical prefix:
#' \itemize{
#'   \item {pfx}.fam: table of N typed individuals
#'   \item {pfx}.bim: table of P typed genomic variants (i.e., SNPs);
#'   \item {pfx}.bed: genotype matrix of N rows and P columns stored in condensed
#'   binary format.
#' }
#'
#' The three files are commonly referred by their common prefix, e.g.:
#'
#' chrX.bed, chrX.fam, and chrX.bim, are jointly specified by "chrX".
#'  
#' @param pfx  character prefix of a PLINK BED file set.
#' @param row  rows and names, NULL by default; set 1 to use IID (2nd column in a FAM)
#' as row names, or 2 to use FID.IID (the 1st and 2nd column in a FAM) as row names;
#' alternatively, one could supply a vector as row names, but must make sure its size
#' being identical to the actual number of individuals recorded by the PLINK files.
#' 
#' @param col  columns and names, NULL by default; set to 1 to use SNP rs ID (2nd column
#' in a BIM) as column names, or set to 2 to use automatically generated "CHR.POS" as
#' column names; alternatively, one could supply a vector as column names, but must be
#' sure its size equals the actual number of variants in the PLINK files.
#'
#' @return matrix of N rows of samples rows, and P columns of biallelic
#' genomic features (i.e., SNPs).
#'
#' @seealso{
#' scanBED}
#' 
#' @examples{
#'
#' ## read PLINK BED 1, 10 individuals, 20  SNPs
#' f1 <- system.file("extdata", paste0('m20', '.bed'), package="plinkBED")
#' print(f1)
#'
#' ## get prefix and read
#' g1 <- readBED(sub("[.]bed$", "", f1))
#' print(g1)
#'
#'
#' ## read PLINK BED 2, 10 individuals, 14K SNPs
#' (f2 <- system.file("extdata", paste0('i10', '.bed'), package="plinkBED"))
#'
#' ## get prefix and read
#' g2 <- readBED(sub("[.]bed$", "", f2))
#' str(g2)
#'
#' }
#' 
#' @export
readBED <- function(pfx, row=NULL, col=NULL)
{
    pfx <- sub("[.]bed", "", pfx)
    ## the triplets
    bedFile <- paste0(pfx, '.bed')
    famFile <- paste0(pfx, '.fam')
    bimFile <- paste0(pfx, '.bim')

    ## read FAM and BIM if necessary, count samples and variants.
    if(is.numeric(row) && length(row) == 1 && row > 0)
    {
        fam <- readFAM(pfx)
        N <- nrow(fam)
        if(row == 1)                    # IID as row name
            row <- fam[, 2]             
        else                            # FID.IID as row name
            row <- paste(fam[, 1], fam[, 2], sep='.')
        rm(fam)
    }
    else if(length(row) > 1)
    {
        row <- as.character(row)
        N <- length(row)
    }
    else
    {
        row <- NULL
        N <- lc(famFile)
    }

    if(is.numeric(col) && col > 0)
    {
        bim <- readBIM(pfx)
        P <- nrow(bim)
        if(col == 1)                    # rs ID as column name
            col <- bim[, 2]
        else                            # CHR.POS as column name
        {
            col <- sprintf("%02d.%09d", bim[, 1], bim[, 3])
        }
        rm(bim)
    }
    else if(length(col) > 1)
    {
        col <- as.character(col)
        P <- length(col)
    }
    else
    {
        col <- NULL
        P <- lc(bimFile)
    }

    ## open the the BED file
    fp <- file(bedFile, open="rb")

    ## first 3 bytes must be 6c 1b 01
    if(any(readBin(fp, "raw", 3L) != as.raw(c(0x6c, 0x1b, 0x01))))
    {
        if(isOpen(fp, "rb"))
        {
            ## print(paste("Close", bedFile, fp))
            close(fp)
        }
        stop("wrong magic numbers, BED triplets may be outdated.")
    }

    ## scan the entire data
    dl <- file.size(bedFile) - 3L        # data length
    bd <- readBin(fp, "raw", dl)         # byte data
    if(isOpen(fp, "rb"))
    {
        ## print(paste("Close", bedFile, fp))
        close(fp)
    }

    ## decompress the byte data
    rt <- db1(bd)

    ## form matrix, and truncate
    dim(rt) <- c(length(rt) %/% P, P)
    rt <- rt[seq.int(N), ]

    ## assign row and column names
    if(!is.null(row))
        rownames(rt) <- row
    if(!is.null(col))
        colnames(rt) <- col
    rt
}

#' Read FAM file
#'
#' @param pfx prefix of a set of PLINK triplets, or the fullname of
#' the bed, fam or bim in that set.
#' @return a data frame describing the individuals, loaded from the
#' *fam* file of the triplets.
#' @export
readFAM <- function(pfx)
{
    pfx <- sub("[.]bed", "", pfx)
    fn <- paste0(pfx, '.fam')
    hdr <- c("fid", "iid", "mom", "dad", "sex", "phe")
    clz <- c(rep("character", 4), rep("integer", 2))
    utils::read.table(fn, FALSE, col.names=hdr, colClasses=clz)
}

#' Read BIM file
#'
#' @param pfx prefix of a set of PLINK triplets, or the fullname of
#' the bed, fam or bim in that set.
#' @return a data frame of describing genomic variants, loaded from
#' the *bim* file of the triplets.
#' @export
readBIM <- function(pfx)
{
    fn <- paste0(pfx, '.bim')
    hdr <- c("chr", "id", "cm", "pos", "a1", "a2")
    clz <- c("integer", "character", "integer", "integer", "character", "character")
    utils::read.table(fn, FALSE, col.names=hdr, colClasses=clz)
}

#' Read all 3 files in a set of PLINK triplets
#' 
#' Read the genotype matrix, individual table, and variant table from a PLINK BED
#' file set. This function is meant for in memory operation of small BED. For out
#' of memory computation tasks involving huge BED, see \code{\link{scanBED}}
#' instead.
#' 
#' A PLINK BED (\emph{binary biallelic genotype table}) is comprised of three files
#' (preferably) using identical prefix:
#' \itemize{
#'   \item {pfx}.fam: table of N typed individuals
#'   \item {pfx}.bim: table of P typed genomic variants (i.e., SNPs);
#'   \item {pfx}.bed: genotype matrix of N rows and P columns stored in condensed
#'   binary format.
#' }
#'
#' The surfix is not mandatory, but is prefered, because in that way all three can
#' be commonly referred by their shared prefix, e.g.:
#'
#' chrX.bed, chrX.fam, and chrX.bim, are jointly specified by "chrX".
#' 
#' @param pfx prefix of a PLINK BED
#' @return a list containing the genomic matrix (bed), data frame of individuals
#' (fam), and data frame of variants (bim).
#' @export
readPLK <- function(pfx)
{
    bim <- readBIM(pfx)
    fam <- readFAM(pfx)
    bed <- readBED(pfx, fam[, 2], bim[, 2])
    list(bed=bed, bim=bim, fam=fam)
}


#' Scan genotypes in PLINK BED(s)
#'
#' Visit one variant at a time, through one or more PLINK BED file sets. This is
#' meant for out-of-memory screening of huge PLINK BED, such as a GWAS study.
#'
#' To read an entire BED into a R matrix, see \code{\link{readBED}} instead.
#' 
#' A BED (\emph{binary biallelic genotype table}) is comprised of three files
#' (usually) sharing identical prefix:
#' \itemize{
#'   \item {pfx}.fam: table of N typed individuals
#'   \item {pfx}.bim: table of P typed genomic variants (i.e., SNPs);
#'   \item {pfx}.bed: genotype matrix of N rows and P columns stored in condensed
#'   binary format.
#' }
#'
#' The three files are commonly referred by their common prefix, e.g.:
#'
#' chrX.bed, chrX.fam, and chrX.bim, are jointly specified by "chrX".
#'  
#' @param pfx  character one or more prefix of PLINK BED to go through.
#' @param FUN  function to used on each genomic variant, it must take the first
#' argument as the vector of variant's values.
#' @param ...  additional argument to pass to \emph{\code{FUN}}.
#'
#' @return list of all computation results with respect to each SNP.
#' genomic features (i.e., SNPs).
#'
#' @seealso {readBED}
#' @export
scanBED <- function(pfx, FUN=summary, ...)
{
    ## the triplets
    bedFile <- paste0(pfx, '.bed')
    famFile <- paste0(pfx, '.fam')
    bimFile <- paste0(pfx, '.bim')

    ## number of samples and features
    N <- lc(famFile)
    P <- lc(bimFile)

    ## open the the BED file
    fp <- file(bedFile, open="rb")

    ## first 3 bytes must be 6c 1b 01
    if(any(readBin(fp, "raw", 3L) != as.raw(c(0x6c, 0x1b, 0x01))))
    {
        if(isOpen(fp, "rb"))
        {
            ## print(paste("Close", bedFile, fp))
            close(fp)
        }
        stop("wrong magic numbers, PLINK BED may be outdated or damaged.")
    }

    ## scan the data
    dl <- file.size(bedFile) - 3L       # bytes remaining
    bpv <- dl %/% P                     # bytes per variant
    csz <- 2L^20L %/% bpv               # chunk size (1MB)
    vdx <- 0L                           # variant index
    ret <- list()
    while (length(.g <- readBin(fp, "raw", csz)) > 0)
    {
        .g <- db1(.g)                   # genotype chunk
        dim(.g) <- c(length(.g) %/% P, P)
        .g <- .g[seq.int(N), ]
        for(j in seq.int(NCOL(.g)))
        {
            vdx <- vdx + 1L
            ret[[vdx]] <- FUN(.g[, j], ...)  # one variant
        }
    }

    ## close the BED file
    if(isOpen(fp, "rb"))
    {
        ## print(paste("Close", bedFile, fp))
        close(fp)
    }
    simplify2array(ret)
}

#' Test BED Reader
#'
#' Let the bedRead take "i10.{bed,bim,fam}" and "m20.{bed,bim,fam}"
#' under "extdata", and compare the integer output with texture files
#' ("i10.txt" and "m20.txt") that were manully converted from the
#' triplets using PLINK 1.9.
#' 
#' @export
bedTest <- function()
{
    fun <- function(pfx)
    {
        txt <- system.file("extdata", paste0(pfx, '.txt'), package="plinkBED")
        ctrl <- scan(txt, 1L, quiet=TRUE)
        bed <- system.file("extdata", paste0(pfx, '.bed'), package="plinkBED")
        read <- readBED(sub("[.]bed", "", bed))
        if(any(t(read) != ctrl, na.rm=TRUE) || any(is.na(t(read)) != is.na(ctrl)))
            stop("Failed BED reading test:", pfx)
        cat("Passed BED reading test: ", pfx, "\n", sep="")
        invisible(TRUE)
    }

    fun('m20')
    fun('i10')
}
