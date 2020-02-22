#' Decompress Byte Data
#'
#' For each SNP  (i.e., a row in the  BED), a byte encodes the up  to 4 genotype
#' samples (2 bits each).
#'
#' The function decodes bytes read from a BED to allele dosage or NA.
#'
#' @param bd byte data in R "raw" mode
#' @return a vector of dosage data in integer, untruncated.
db1 <- function(bd, N, P, quiet=TRUE)
{
    M <- as.integer(ceiling(N / 4) * 4L)
    r <- matrix(0L, M, P) # data
    ## 00->2: homozygous A1
    ## 01->3: missing
    ## 10->1: heterozygous
    ## 11->0: homozygous A2
    C <- c(2L, NA, 1L, 0L)

    ## fill the genomic data vector
    m <- as.raw(0x03) # bit mask: 0x03 = 11000000
    if(!quiet)
        cat("converting ...\n")

    if(!quiet)
        cat("sample 1, 5,  9, ...\n")
    r[seq.int(1L, M, 4L), ] <- C[as.integer(bd                & m) + 1L]

    if(!quiet)
        cat("sample 2, 6, 10, ...\n")
    r[seq.int(2L, M, 4L), ] <- C[as.integer(rawShift(bd, -2L) & m) + 1L]

    if(!quiet)
        cat("sample 3, 7, 11, ...\n")
    r[seq.int(3L, M, 4L), ] <- C[as.integer(rawShift(bd, -4L) & m) + 1L]

    if(!quiet)
        cat("sample 4, 8, 12, ...\n")
    r[seq.int(4L, M, 4L), ] <- C[as.integer(rawShift(bd, -6L) & m) + 1L]
    
    r <- r[seq.int(1L, N), ]
    r
}


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
    rt <- db1(bd, N, P)

    ## assign row and column names
    ## if(!is.null(row))
    ##     rownames(rt) <- row
    ## if(!is.null(col))
    ##     colnames(rt) <- col
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
scanBED <- function(pfx, FUN, ...)
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
    dl <- file.size(bedFile) - 3L         # bytes remaining
    bpv <- dl %/% P                       # bytes per variant
    vpc <- min(P, max(2L^20L %/% bpv, 1)) # variant per chunk (chunk size = 1MB)
    vdx <- 0L                             # variant index
    ret <- list()
    while (length(.g <- readBin(fp, "raw", vpc)) > 0)
    {
        .p <- length(.g) / N
        .g <- db1(.g, N, P) # genotype chunk
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

    r <- simplify2array(ret)
    if(is.matrix(r))
        r <- t(r)
    r
}

#' Test BED Reader
#'
#' Read m20 (bed, bim, and fam) under  "extdata" and compare with the content in
#' text file "i10.txt" converted from m20 by PLINK.
test.readBed <- function()
{
    pfx <- sub("[.]bed$", "", system.file("extdata", 'm20.bed', package="plinkBED"))
    txt <- system.file("extdata", paste0(pfx, '.txt'), package="plinkBED")

    bed <- readBED(pfx)
    txt <- scan(txt, 1L, quiet=TRUE)

    if(any(t(bed) != txt, na.rm=TRUE) || any(is.na(t(bed)) != is.na(txt)))
        stop("Failed BED reading test:", pfx)
    cat("Passed BED reading test: ", pfx, "\n", sep="")
    invisible(TRUE)
}


#' Test BED Scanner
#'
#' Go through  i20.* (bed, bim,  and fam) under  "extdata", summerize every SNP.
test.scanBed <- function()
{
    pfx <- sub("[.]bed$", "", system.file("extdata", 'i10.bed', package="plinkBED"))
    bed <- readBED(pfx)
    
    ret <- scanBED(pfx, function(g)
    {
        af <- mean(g, na.rm=TRUE) / 2
        maf <- min(af, 1 - af)
        c(mu=mean(g, na.rm=TRUE), sd=sd(g, na.rm=TRUE), maf=maf, nas=sum(is.na(g)))
    })
    ret
}
