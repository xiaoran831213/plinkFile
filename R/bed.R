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

#' Read Genotype in PLINK BED Triplets
#'
#' The function take in the shared prefix of the BED triplets, return the
#' P genomic features of N individuals organized in P by N matrix.
#' 
#' PLINK BED comes with three files: a *.fam as the table of N individuals,
#' a *.bim as the table of P genomic features, and a *.bed archiving the
#' actual data matrix in condensed binary form. The triplets usually share
#' the prefix to ease referencing, for example:
#'
#' the triplets chrX.bed, chrX.fam, and chrX.bim, can be collectively
#' refered by their shared prefix "chrX".
#'  
#' @param pfx  the prefix shared by the PLINK BED triplets.
#' @param type change the storage type of read values.
#'
#' @return a matrix of N individual samples in the rows, and P genomic
#' features in the columns.
#'
#' @export
readBED <- function(pfx, type="integer")
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
            print(paste("Close", bedFile, fp))
            close(fp)
        }
        stop("wrong magic numbers, BED triplets may be outdated.")
    }
    ## scan the data
    L <- file.size(bedFile) - 3L        # bytes remaining
    bd <- readBin(fp, "raw", L)
    if(isOpen(fp, "rb"))
    {
        print(paste("Close", bedFile, fp))
        close(fp)
    }

    ## fill the genomic data vector
    gd <- raw(L * 4L)
    .i <- seq.int(1L, by=4L, length.out=L)
    .m <- as.raw(0x03)                    # 0x03 = 11000000
    gd[.i     ] <- bd & .m                # sample 1, 5, ...
    gd[.i + 1L] <- rawShift(bd, -2L) & .m # sample 2, 6, ...
    gd[.i + 2L] <- rawShift(bd, -4L) & .m # sample 3, 7, ...
    gd[.i + 3L] <- rawShift(bd, -6L) & .m # sample 4, 8, ...
    
    ## recode to dosage
    rt <- raw(L * 4L)
    rt[gd==as.raw(0L)] <- as.raw(2L)    # 00->2: homozygous A1
    rt[gd==as.raw(1L)] <- as.raw(3L)    # 01->3: missing
    rt[gd==as.raw(2L)] <- as.raw(1L)    # 10->1: heterozygous
    rt[gd==as.raw(3L)] <- as.raw(0L)    # 11->0: homozygous A2

    ## storage type
    fc <- switch(
        type,
        integer=as.integer, numeric=as.numeric,
        double=as.double,
        character=as.character,
        raw=as.raw,
        as.integer)
    rt <- fc(rt)

    dim(rt) <- c(length(rt) %/% P, P)
    rt <- rt[seq.int(N), ]
    rt
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
        if(!all(t(read) == ctrl))
            stop("Failed BED reading test:", pfx)
        cat("Passed BED reading test: ", pfx, "\n", sep="")
        invisible(TRUE)
    }

    fun('m20')
    fun('i10')
}
