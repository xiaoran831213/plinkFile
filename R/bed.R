#' Decode Byte Data
#'
#' Decodes bytes read from a BED into allele dosage or NA.
#'
#' For each SNP (a logical line in the  BED), each byte decode into to 4 dosage
#' samples. (2 bits each).
#'
#' @param B byte data in R "raw" mode
#' @param N number of individuals in the byte data.
#' @param quiet do not report (def=TRUE)
#' 
#' @return a N x P matrix of genotype, where P is the number of variants.
#' @noRd
dbd <- function(B, N, quiet=TRUE)
{
    M <- as.integer(ceiling(N / 4) * 4L) # round to 4 persons
    P <- length(B) / M * 4               # number of variants

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
        r <- r[seq.int(1L, N), , drop=FALSE]
    r
}

#' Encode Byte Data
#'
#' #' Decodes bytes read from a BED into allele dosage or NA.
#'
#' For each SNP (a logical line in the  BED), each byte decode into to 4 dosage
#' samples. (2 bits each).
#'
#' @param D N x P dosage genotype data
#' @param N number of individuals
#' @param quiet do not report (def=TRUE)
#' 
#' @return a length `ceiling(N/4)` vector of bytes in R "raw" mode.
#' @noRd
ebd <- function(D, N, quiet=TRUE)
{
    M <- as.integer(ceiling(N / 4)) # round to 4 persons per byte
    P <- length(D) / N              # number of variants
    D[is.na(D)] <- 3L               # 3 <- NA
    
    ## bytes required
    r <- matrix(raw(M * P), M, P)
    
    ## Dictionary of 2 bits codings:
    ## 0->11: homozygous A2
    ## 1->10: heterozygous
    ## 2->00: homozygous A1
    ## 3->01: missing
    C <- as.raw(c(3, 2, 0, 1))

    ## fill the bytes
    m <- as.raw(0x03) # bit mask: 0x03 = 11000000
    if(!quiet)
        cat("converting ...\n")

    if(!quiet)
        cat("sample 1, 5,  9, ...\n")
    i <- seq.int(1L, N , 4L); j <- seq_along(i)
    r[j, ] <- r[j, ] | rawShift(C[D[i, ] + 1], 0)

    if(!quiet)
        cat("sample 2, 6, 10, ...\n")
    i <- seq.int(2L, N , 4L); j <- seq_along(i)
    r[j, ] <- r[j, ] | rawShift(C[D[i, ] + 1], 2)

    if(!quiet)
        cat("sample 3, 7, 11, ...\n")
    i <- seq.int(3L, N , 4L); j <- seq_along(i)
    r[j, ] <- r[j, ] | rawShift(C[D[i, ] + 1], 4)

    if(!quiet)
        cat("sample 4, 8, 12, ...\n")
    i <- seq.int(4L, N , 4L); j <- seq_along(i)
    r[j, ] <- r[j, ] | rawShift(C[D[i, ] + 1], 6)
    
    dim(r) <- NULL
    r
}

#' Navigate Byte Data
#'
#' Jump to the start of a variant in a BED data.
#'
#' Each variant is a logical line in  the BED file, and each byte encode allele
#' dosage of 4 individuals  (2 bits each).  The beginning of  the K th. variant
#' is the K * up(N / 4) th. byte of the BED data stream.
#'
#' @param C opened file pointing to a BED
#' @param N number of individuals in the byte data.
#' @param j integer where to jump?
#' 
#' @return a N x P matrix of genotype, where P is the number of variants.
#' @examples
#' bed <- system.file("extdata", "000.bed", package="plinkFile")
#' fam <- system.file("extdata", "000.fam", package="plinkFile")
#' num <- length(readLines(fam))
#' 
#' con <- file(bed, open="rb")
#' close(con)
#' @noRd
nbd <- function(C, N, j=0)
{
    ## bytes per variant, rounded to 4 persion per byte
    M <- as.integer(ceiling(N / 4))

    ## move the file pointer
    seek(C, origin="current", where=M * j)
}

#' Read BED file
#'
#' Read a PLINK BED file into a R matrix.
#'
#' This is meant  for genotype that can  fit into system memory; the  size of R
#' matrix is 16 times the size of the  BED file. To scan a large BED one varant
#' at time without loading it entirely into the memoty, see [scanBED()].
#'
#' A PLINK1 binary fileset has three files,
#' \itemize{
#'   \item{`pfx.fam`}: {text table of `N` individuals.}
#'   \item{`pfx.bim`}: {text table of `P` genomic variants (i.e., SNPs).}
#'   \item{`pfx.bed`}: {N x P genotype matrix in condensed binary format.}}
#' 
#' The three files  comprising a genotype data are typically  referred by their
#' common  prefix,  for  example,  the X  chromosome  genotype  represented  by
#' `chrX.bed`, `chrX.fam`, and `chrX.bim` are jointly refered by `chrX`.
#' 
#' @param pfx prefix of PLINK file set, or the fullname of a BED file.
#' @param row row names: 1=IID, 2=FID.IID, def=NULL.
#' @param col col names: 1=variant ID, 2= CHR:POS, 3=CHR:POS_A1_A2, def=NULL
#' @param vfr from which variant? (def=1st)
#' @param vto till which variant? (def=EOF)
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
readBED <- function(pfx, row=NULL, col=NULL, vfr=NULL, vto=NULL, quiet=TRUE)
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
        N <- nr(famFile)
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
        P <- nr(bimFile)
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


#' Save BED file
#'
#' Save a R matrix into a PLINK BED file.
#'
#' This is meant for genotype small enough  to fit into system memory. The size
#' of R matrix is 16  times the size of the BED file.
#' 
#' @param pfx prefix of the output file set, in PLINK1 BED format.
#' @param bed N x P genotype matrix
#' @param quiet do not report (def=TRUE)
saveBED <- function(pfx, bed, quiet=TRUE)
{
    pfx <- sub("[.]bed", "", pfx)
    bedFile <- paste0(pfx, '.bed')

    ## open the BED file
    fp <- file(bedFile, open="wb")

    ## first 3 bytes must be 6c 1b 01
    mb <- as.raw(c(0x6c, 0x1b, 0x01))
    writeBin(mb, fp, 1L)
    
    ## rest of the data, binary encoded
    bd <- ebd(bed, nrow(bed), quiet=quiet) # byte data
    writeBin(bd, fp, 1L)

    ## close the BED file
    close(fp)

    invisible(NULL)
}


#' Read FAM file
#'
#' Read sample meta-data form the *fam* file of a PLINK1 BED fileset.
#'
#' There are six columns in a *bim* file
#' * fid: family ID;
#' * iid: individual ID, default row name used by `[readBED]`;
#' * mom: maternal ID;
#' * dad: paternal ID;
#' * sex: individual sex.
#' * phe: phenotype, not often used;
#'
#' The PLINK1 *bim* file has no header line, this is changed in PLINK2.
#'
#' The columns "sex" and "phe" are  mostly the legency of early GWAS, nowerdays
#' it is common to provide sex, among other covariates, and multiple phenotypes
#' in a separate file.
#' 
#' @param pfx prefix of a PLINK file set.
#' @return meta data frame describing individuals, loaded from "{_pfx_.fam}".
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' fam <- readFAM(pfx)
#' fam
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
#' Get variant meta-data form the *bim* file of a PLINK1 BED fileset.
#'
#' @details
#' There are six columns in a *bim* file
#' * chr: chromosme of the variant
#' * vid: variant id, such as an RS number;
#' * cmg: position by centimorgan;
#' * bps: position by basepairs;
#' * al1: allele 1, the one counted as dosage.
#' * al2: allele 2.
#' 
#' @param pfx prefix of a PLINK file set.
#' @return meta data frame describing the variants, loaded from "{_pfx_.bim}".
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' bim <- readBIM(pfx)
#' bim
#' @export
readBIM <- function(pfx)
{
    fn <- paste0(pfx, '.bim')
    hdr <- c("chr", "vid", "cmg", "bps", "al1", "al2")
    clz <- c("character", "character", "integer", "integer", "character", "character")
    utils::read.table(fn, FALSE, col.names=hdr, colClasses=clz)
}

#' read variant ID
#'
#' Generate variant ID automatically, or based on a *bim* file.
#'
#' @details
#' The option (`opt`) can be:
#' *  1 = the 2nd column in _pfx_.bim (default),
#' *  2 = formated as %CHR(02d):%BPS(09d),
#' *  3 = formated as %CHR(02d):%BPS(09d)_AL1(s)_AL2(s)
#' *  0 = nothing
#' * -1 = numbering of variants, decimal
#' * -2 = numbering of variants, zero-padded, fixed length decimal
#' * -3 = numbering of variants, zero-padded, fixed length hexedemical
#' * or, a vector of IDs to use.
#' 
#' @param pfx prefix of a PLINK file set.
#' @param opt option (def=1: the 2nd column in _pfx_.bim).
#' @param bim use existing table loaded via `readBIM(pfx)` (def=NULL).
#' @return a vector of variant ID
#'
#' @examples
#' # read variant ID
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#'
#' # opt=1: 2nd column in the BED file (default)
#' vid <- readVID(pfx, 1); head(vid); tail(vid)
#'
#' # opt=2: format by position
#' vid <- readVID(pfx, 2); head(vid); tail(vid)
#'
#' # opt=3: format by position and alleles
#' vid <- readVID(pfx, 3); head(vid); tail(vid)
#' @export
readVID <- function(pfx, opt=NULL, bim=NULL)
{
    if(is.null(opt))
        opt <- 1
    if(is.numeric(opt) && length(opt) == 1)
    {
        if(opt > 0)
        {
            if(is.null(bim))
                bim <- readBIM(pfx)
            vid <- with(bim,
            {
                if(opt == 1)
                    vid
                else if(opt == 2)
                    sprintf("%02d:%09d", as.integer(chr), bps)
                else
                    sprintf("%02d:%09d_%s_%s", as.integer(chr), bps, al1, al2)
            })
            P <- length(vid)
        }
        else
        {
            P <- if(is.null(bim)) nr(paste0(pfx, ".bim")) else nrow(bim)
            if(length(opt) == 1)
            {
                if(opt == 0)
                    vid <- NULL
                else if(opt == -1)
                    vid <- sprintf("%d", seq(P))
                else if(opt == -2)
                    vid <- .name(P, alpha=c(0:9))
                else
                    vid <- .name(P, alpha=c(0:9, letters[1:6]), pfx="0x")
            }
        }
    }
    else if(length(opt) > 1)
    {
        P <- if(is.null(bim)) nr(paste0(pfx, ".bim")) else nrow(bim)
        if(length(opt) != P)
            stop(gettextf("unequal num of IDs and variants (%d and %d).", length(opt), P))
        vid <- opt
    }
    else
        stop(gettextf("fail to get variant IDs from \"%s\" with opt=%s.", pfx, opt))
    vid
}

#' @describeIn bed apply a function to variants in a PLINK1 BED fileset
#'
#' The scripts in `FUN` has no side effects on environment of `scanBED`.
#'
#' @param FUN a function to process each window of variants;
#' @param ... additional argument for *`FUN`* when `scanBED` is used.
#' @param win window size - the number of variants passed to `FUN` per iteration.
#' @param vid variant ID option (def=1, use the 2nd column in _pfx_.bim).
#' @examples
#' ## traverse genotype, apply R function without side effects
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "000")
#' ret <- scanBED(pfx, function(g)
#' {
#'     .af <- colMeans(g, na.rm=TRUE) / 2
#'     maf <- pmin(.af, 1 - .af)
#'     mis <- colSums(is.na(g)) / .N
#'     wnd <- .w
#'     pct <- round(.w / .W * 100, 2)
#'     cbind(buf=.b, wnd=.w, idx=.i, MAF=maf, MIS=mis, PCT=pct)
#' },
#' win=13, simplify=rbind, buf=2^18)
#' head(ret)
#' tail(ret)
#'
#' @export
scanBED <- function(pfx, FUN, ..., win=1, vid=NULL, buf=2^24, simplify=TRUE)
{
    ## PLINK file set
    bedFile <- paste0(pfx, '.bed')
    famFile <- paste0(pfx, '.fam')
    bimFile <- paste0(pfx, '.bim')

    ## number of samples and features
    N <- nr(famFile)
    vid <- readVID(pfx, vid)
    P <- if(is.null(vid)) nr(bimFile) else length(vid)

    ## open the the BED file
    fp <- fp <- .try.open.bed(bedFile)

    ## praperation
    byt <- file.size(bedFile) - 3L  # bytes to go through
    bpv <- byt %/% P                # bytes per variant
    vpb <- max(buf %/% bpv, win)    # variants per chunk
    vpb <- ceiling(vpb / win) * win # round up to full window
    bpc <- bpv * vpb                # bytes per chunk
    wdx <- 0L                       # window index
    bdx <- 0L                       # buffer chunk index
    env <- environment(FUN)         # processing environment
    env[[".P"]] <- P
    env[[".N"]] <- N
    env[[".W"]] <- as.integer(ceiling(P / win))
    env[[".B"]] <- as.integer(ceiling(byt / bpc))
    ## go through the binary data
    ret <- list()
    while (length(chk <- readBin(fp, "raw", bpc)) > 0)
    {
        bdx <- bdx + 1
        env[[".b"]] <- bdx
        ## cat("Chunk #", bdx, "\n") # process a chunk
        gmx <- dbd(chk, N)               # genotype matrix
        if(!is.null(vid))                # vid -> column name
            colnames(gmx) <- vid[seq(wdx * win + 1, len=ncol(gmx))]
        for(j in seq(1, ncol(gmx), win)) # go through windows
        {
            jdx <- seq(j, min(j + win - 1, ncol(gmx)))
            idx <- wdx * win + seq_along(jdx)
            wdx <- wdx + 1L
            ## cat("Window #", wdx, "\n") # process a window
            env[[".i"]] <- idx
            env[[".w"]] <- wdx
            env[[".p"]] <- length(jdx)
            ## process a window of variants
            ret[[wdx]] <- FUN(gmx[, jdx, drop=FALSE], ...)
        }
        ## msg <- sprintf("%4d %3d %4d\n", bdx, ncol(gmx), length(chk))
        ## cat(msg)
    }

    ## close the BED file
    if(isOpen(fp, "rb"))
    {
        ## print(paste("Close", bedFile, fp))
        close(fp)
    }

    ret <- .try.simplify(ret, simplify)
    ret
}

#' @describeIn bed evaluate an expression on variants in a PLINK1 BED
#'
#' The scripts in `EXP` have side effects on the environment of `loopBED`.
#' @param EXP a R expression to evaluate with each window of variants;
#' @param GVR a R variable name to assign the window to (def="g").
#' @param vid option for variant IDs (def=1, the 2nd column in _pfx_.bim).
#' @examples
#' ## traversing genotype, evaluate R expression with side effects
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "000")
#' ret <- list() # use side effect to keep the result of each window.
#' loopBED(pfx,
#' {
#'     af <- colMeans(gt, na.rm=TRUE) / 2
#'     sg <- af * (1 - af)
#'     ret[[.w]] <- cbind(wnd=.w, alf=af, var=sg)
#' },
#' win=13, GVR="gt", vid=3, buf=2^18)
#' head(ret)
#' tail(ret)
#'
#' @export
loopBED <- function(pfx, EXP, GVR="g", win=1, vid=1, buf=2^24, simplify=TRUE)
{
    ## PLINK file set
    bedFile <- paste0(pfx, '.bed')
    famFile <- paste0(pfx, '.fam')
    bimFile <- paste0(pfx, '.bim')

    ## number of samples and features
    N <- nr(famFile)
    vid <- readVID(pfx, vid)
    P <- if(is.null(vid)) nr(bimFile) else length(vid)

    ## open the the BED file
    fp <- .try.open.bed(bedFile)

    ## scan the data
    byt <- file.size(bedFile) - 3L  # bytes to go through
    bpv <- byt %/% P                # bytes per variant
    vpb <- max(buf %/% bpv, win)    # variants per chunk
    vpb <- ceiling(vpb / win) * win # round up to full window
    bpc <- bpv * vpb                # bytes per chunk
    wdx <- 0L                       # window index
    bdx <- 0L                       # buffer index
    EXP <- substitute(EXP)          # expression
    env <- parent.frame()           # processing environment
    env[[".P"]] <- P
    env[[".N"]] <- N
    env[[".W"]] <- as.integer(ceiling(P / win))
    env[[".B"]] <- as.integer(ceiling(byt / bpc))
    ret <- list()
    while (length(chk <- readBin(fp, "raw", bpc)) > 0)
    {
        bdx <- bdx + 1
        env[[".b"]] <- bdx
        ## cat("Chunk #", bdx, "\n") # process a chunk
        gmx <- dbd(chk, N)               # genotype matrix
        if(!is.null(vid))                # vid -> column name
            colnames(gmx) <- vid[seq(wdx * win + 1, len=ncol(gmx))]
        for(j in seq(1, ncol(gmx), win)) # go through windows
        {
            jdx <- seq(j, min(j + win - 1, ncol(gmx)))
            idx <- wdx * win + seq_along(jdx)
            wdx <- wdx + 1L
            ## cat("Window #", wdx, "\n") # process a window
            env[[".i"]] <- idx
            env[[".w"]] <- wdx
            env[[".p"]] <- length(jdx)
            env[[GVR]] <- gmx[, jdx, drop=FALSE]
            ## process a window of variants
            ret[[wdx]] <- eval(EXP, env)
        }
        ## cat(sprintf("%4d %3d %4d\n", bdx, ncol(gmx), length(chk)))
    }
    ## close the BED file
    if(isOpen(fp, "rb"))
    {
        close(fp)
    }

    ret <- .try.simplify(ret, simplify)
    invisible(ret)
}

.try.open.bed <- function(fn)
{
    ## open the the BED file
    fp <- file(fn, open="rb")

    ## first 3 bytes must be 6c 1b 01
    if(any(readBin(fp, "raw", 3L) != as.raw(c(0x6c, 0x1b, 0x01))))
    {
        if(isOpen(fp, "rb"))
        {
            ## print(paste("Close", bedFile, fp))
            close(fp)
        }
        stop(gettextf("wrong magic numbers, \"%s\" compromised.", fn))
    }
    fp
}

.try.simplify <- function(ret, simplify)
{
    .r <- NULL
    ## try provided simplifier
    if(is.character(simplify) || is.function(simplify))
        .r <- try(do.call(simplify, ret), silent = TRUE)
    ## try default simplifier
    if(inherits(.r, 'try-error') || isTRUE(simplify))
        .r <- try(simplify2array(ret), silent = TRUE)
    if(!is.null(.r) && !inherits(.r, 'try-error'))
        ret <- .r
    ret
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
