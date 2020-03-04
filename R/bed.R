#' Decompress Byte Data
#'
#' For each SNP  (i.e., a row in the  BED), a byte encodes the up  to 4 genotype
#' samples (2 bits each).
#'
#' The function decodes bytes read from a BED to allele dosage or NA.
#'
#' @param B byte data in R "raw" mode
#' @param N number of individuals in the byte data.
#' @param quiet do not report (def=TRUE)
#' 
#' @return a N x P matrix of genotype, where P is the number of variants.
dbd <- function(B, N, quiet=TRUE)
{
    M <- as.integer(ceiling(N / 4) * 4L) # round up to 4 persons
    P <- length(B) / M * 4              # infer the number of variants

    ## genotype matrix to be returned
    r <- matrix(0L, M, P)

    ## Dictionary of 2 bits codings:
    ## 00->2: homozygous A1
    ## 01->3: missing
    ## 10->1: heterozygous
    ## 11->0: homozygous A2
    C <- c(2L, NA, 1L, 0L)

    ## fill the genotype matrix
    m <- as.raw(0x03) # bit mask: 0x03 = 11000000
    if(!quiet)
        cat("converting ...\n")

    if(!quiet)
        cat("sample 1, 5,  9, ...\n")
    r[seq.int(1L, M, 4L), ] <- C[as.integer(B                & m) + 1L]

    if(!quiet)
        cat("sample 2, 6, 10, ...\n")
    r[seq.int(2L, M, 4L), ] <- C[as.integer(rawShift(B, -2L) & m) + 1L]

    if(!quiet)
        cat("sample 3, 7, 11, ...\n")
    r[seq.int(3L, M, 4L), ] <- C[as.integer(rawShift(B, -4L) & m) + 1L]

    if(!quiet)
        cat("sample 4, 8, 12, ...\n")
    r[seq.int(4L, M, 4L), ] <- C[as.integer(rawShift(B, -6L) & m) + 1L]

    ## truncate padded, non-existing individuals
    if (N < M)
        r <- r[seq.int(1L, N), ]
    r
}


#' Read BED file
#'
#' Read a BED file into a R matrix. This is meant for in-of-memory process of moderate
#' to small sized genotype.
#'
#' To scan a huge BED one varant at time without reading it into the memoty, see
#' \code{\link{scanBED}} instead.
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
#' @param pfx prefix of PLINK file set, or the fullname of a BED file.
#' @param row  the row names: 1 =  use individual ID, 2 =  family and individual
#'     ID, def = NULL.
#' @param col the column names: 1 =  use variant ID (i.e., rsID), 2 = CHR:POS, 3
#'     = CHR:POS_A1_A2
#' @param quiet suppress screen printing? (def=TRUE)
#' @return genotype matrix with row individuals and column variants.
#' 
#' @examples
#' bed <- system.file("extdata", 'm20.bed', package="plinkFile")
#' pfx <- sub("[.]bed$", "", bed)
#' bed <- readBED(pfx, quiet=FALSE)
#'
#' @seealso {readBED}
#' @export
readBED <- function(pfx, row=NULL, col=NULL, quiet=TRUE)
{
    pfx <- sub("[.]bed", "", pfx)
    ## the triplets
    bedFile <- paste0(pfx, '.bed')
    famFile <- paste0(pfx, '.fam')
    bimFile <- paste0(pfx, '.bim')
    
    ## read FAM and BIM if necessary, count samples and variants.
    if(is.numeric(row) && row > 0)
    {
        fam <- readFAM(pfx)
        N <- nrow(fam)
        if(row == 1) # IID as row name
            row <- fam[, 2]             
        else # FID.IID as row name
            row <- paste(fam[, 1], fam[, 2], sep='.')
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
        if(col == 1) # rs ID as column name
            col <- bim[, 2]
        else if(col == 2) # CHR:POS as column name
            col <- sprintf("%02d:%09d", bim[, 1], bim[, 4])
        else # CHR:POS_A1_A2 as column name
            col <- sprintf("%02d:%09d_%s_%s", bim[, 1], bim[, 4], bim[, 5], bim[, 6])
        rm(bim)
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
    rt <- dbd(bd, N, quiet=quiet)

    ## assign row and column names
    rownames(rt) <- row
    colnames(rt) <- col
    rt
}

#' Read FAM file
#'
#' @param pfx prefix of a PLINK file set.
#' @return data frame describing individuals,  loaded from the FAM file.
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' bed <- readBED(pfx, row=1, col=1, quiet=FALSE)
#' bed
#' 
#' @export
readFAM <- function(pfx)
{
    fn <- paste0(pfx, '.fam')
    hdr <- c("fid", "iid", "mom", "dad", "sex", "phe")
    clz <- c(rep("character", 4), rep("integer", 2))
    utils::read.table(fn, FALSE, col.names=hdr, colClasses=clz)
}

#' Read BIM file
#'
#' @param pfx prefix of a PLINK file set.
#' @return data frame describing genome variants, loaded from the BIM file.
#' @export
readBIM <- function(pfx)
{
    fn <- paste0(pfx, '.bim')
    hdr <- c("chr", "id", "cm", "pos", "a1", "a2")
    clz <- c("integer", "character", "integer", "integer", "character", "character")
    utils::read.table(fn, FALSE, col.names=hdr, colClasses=clz)
}


#' Scan genotypes in PLINK BED(s)
#'
#' Go through a BED file set and visit  one variant at a time. This is meant for
#' out-of-memory screening of huge genotype, such as a GWAS study.
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
#' @param pfx prefix of PLINK BED.
#' @param FUN the function to process each variant.
#' @param ...  additional argument to pass to \emph{\code{FUN}}.
#' @param simplify TRUE to simplify the result as an array.
#'
#' @return an array with each row  corresponding to a variant if simplify is set
#'     to TRUE; otherwise,  a list with each element corresponding  to a variant
#'     is returned.
#'
#' A  context  vaiable  ".i"  is  assigned to  the  environment  of  \code{FUN},
#' therefore, one can  access the index of variant current  being processed from
#' within the body of \code{FUN}.
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "000")
#' ret <- scanBED(pfx, function(g)
#' {
#'     af <- mean(g, na.rm=TRUE) / 2
#'     maf <- min(af, 1 - af)
#'     c(idx=.i, mu=mean(g, na.rm=TRUE), maf=maf, nas=sum(is.na(g)))
#' })
#' print(ret[1:5, ])
#'
#' @seealso {readBED}
#' @export
scanBED <- function(pfx, FUN, ..., simplify=TRUE)
{
    ## PLINK file set
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
    csz <- 2^19                         # chunk size: 512KB
    dl <- file.size(bedFile) - 3L       # bytes remaining
    bpv <- dl %/% P                     # bytes per variant
    vpc <- min(P, max(csz %/% bpv, 1))  # variants per chunk
    bpc <- bpv * vpc                    # bytes per chunk
    idx <- 0L                           # variant index
    cdx <- 1L                           # chunk index
    ret <- list()
    while (length(chk <- readBin(fp, "raw", bpc)) > 0)
    {
        ## cat("Chunk #", cdx, "\n") # process a chunk
        gmx <- dbd(chk, N)        # genotype matrix
        for(j in seq(ncol(gmx)))  # go through variants
        {
            idx <- idx + 1L
            environment(FUN)[[".i"]] <- idx  # contex: index
            ret[[idx]] <- FUN(gmx[, j], ...) # one variant
        }
        cdx <- cdx + 1
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
testReadBED <- function()
{
    pfx <- sub("[.]bed$", "", system.file("extdata", 'm20.bed', package="plinkFile"))
    txt <- system.file("extdata", "m20.txt", package="plinkFile")

    bed <- readBED(pfx, quiet=FALSE)
    txt <- unname(as.matrix(utils::read.table(txt)))
    
    if(any(bed != txt, na.rm=TRUE) || any(is.na(bed) != is.na(txt)))
        stop("Failed BED reading test:", pfx)
    cat("Passed BED reading test: ", pfx, "\n", sep="")
    invisible(TRUE)
}


#' Test BED Scanner
#'
#' Go through file set "000" under "extdata", summerize every SNP.
testScanBED <- function()
{
    pfx <- sub("[.]bed$", "", system.file("extdata", '000.bed', package="plinkFile"))
    ret <- scanBED(pfx, function(g)
    {
        af <- mean(g, na.rm=TRUE) / 2
        maf <- min(af, 1 - af)
        c(mu=mean(g, na.rm=TRUE), maf=maf, nas=sum(is.na(g)))
    })
    ret
}
