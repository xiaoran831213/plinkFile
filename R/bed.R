## DBT <- dbd(as.raw(seq(0x00, 0xFF)), 4L)

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

DBD <- function(B, N, quiet=TRUE)
{
    M <- as.integer(ceiling(N / 4) * 4L) # round to 4 persons
    P <- length(B) / M * 4               # number of variants

    ## genotype matrix to be returned
    r <- DBT[, as.integer(B) + 1]
    dim(r) <- c(M, P)

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
    ## bytes per variant, rounded to 4 person per byte
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
#' @param iid IID as row names (def=1, see [readIID()]). 
#' @param vid VID as col names (def=1, see [readVID()]).
#' @param vfr from which variant? (def=1st)
#' @param vto till which variant? (def=EOF)
#' @param quiet suppress screen printing? (def=TRUE)
#' @return genotype matrix with row individuals and column variants.
#' 
#' @examples
#' ## read an entire small data
#' bed <- system.file("extdata", 'm20.bed', package="plinkFile")
#' gmx <- readBED(bed, quiet=FALSE)
#'
#' ## read part of a large data
#' bed <- system.file("extdata", '000.bed', package="plinkFile")
#' U <- readBED(bed, vfr=01, vto=10, quiet=FALSE)
#' V <- readBED(bed, vfr=11, vto=20, quiet=FALSE)
#' W <- cbind(U, V)
#' X <- readBED(bed, vfr=01, vto=20, quiet=FALSE)
#' all.equal(W, X)
#'
#' @seealso {readBED}
#' @export
readBED <- function(pfx, iid=1, vid=1, vfr=NULL, vto=NULL, quiet=TRUE)
{
    pfx <- sub("[.](bim|fam|bed)$", "", bim)
    ## the triplets
    bedFile <- paste0(pfx, '.bed')
    famFile <- paste0(pfx, '.fam')
    bimFile <- paste0(pfx, '.bim')
    
    ## number of samples and features
    iid <- readIID(pfx, iid)
    vid <- readVID(pfx, vid)
    N <- if(is.null(iid)) nr(famFile) else length(iid)
    P <- if(is.null(vid)) nr(bimFile) else length(vid)

    ## get range
    vfr <- if(is.null(vfr)) 1 else min(max(vfr, 0), P) # variant-from
    vto <- if(is.null(vto)) P else min(max(vto, 0), P) # variant-to
    if(vfr < 1)     # in case of fraction
    {
        vfr <- as.integer(vfr * P) + 1
        if(vto == 1)
            vto <- P
    }
    if(vto < 1)
        vto <- as.integer(vto * P)
    if(vto < vfr)
        stop(gettextf("to-variant (%d) preceeds from-variant (%d)", vto, vfr))

    ## praperation
    bpv <- (file.size(bedFile) - 3L)  %/% P # bytes per variant
    vid <- vid[vfr:vto]                     # update VID
    nvr <- vto - vfr + 1                    # number of variants
    byt <- nvr * bpv                        # number of bytes

    fp <- .try.open.bed(bedFile)
    seek(fp, origin="current", where=(vfr - 1L) * bpv)  # move file pointer
    bd <- readBin(fp, "raw", byt)                       # byte data
    
    if(isOpen(fp, "rb"))
    {
        ## print(paste("Close", bedFile, fp))
        close(fp)
    }
    
    ## decompress the byte data
    rt <- dbd(bd, N, quiet=quiet)

    ## assign row and column names
    rownames(rt) <- iid
    colnames(rt) <- vid
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
#' @param fam prefix or name of a PLINK file.
#' @return data frame of individuals, loaded from FAM.
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' fam <- readFAM(pfx)
#' fam
#' @export
readFAM <- function(fam)
{
    pfx <- sub("[.](bim|fam|bed)$", "", fam)
    hdr <- c("fid", "iid", "mom", "dad", "sex", "phe")
    clz <- c(rep("character", 4), rep("integer", 2))
    utils::read.table(paste0(pfx, '.fam'), FALSE, col.names=hdr, colClasses=clz)
}

#' read individual ID
#'
#' Generate individual ID automatically, or based on a *fam* file.
#'
#' @details
#' The option (`opt`) can be:
#' *  1 = the *iid* column in FAM (default),
#' *  2 = formated as *fid*.*iid*,
#' *  0 = nothing
#' * -1 = numbering of individuals, decimal
#' * -2 = numbering of individuals, zero-padded fix-length decimal
#' * -3 = numbering of individuals, zero-padded fix-length hexedemical
#' or, a vector of IDs to use.
#' 
#' @param fam prefix or name of a PLINK file, or data fram from a FAM file.
#' @param opt option (def=1: the 2nd column in FAM file).
#' @return a vector of individual ID
#'
#' @examples
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "m20")
#' readIID(pfx,  1) # opt= 1: IID
#' readIID(pfx,  2) # opt= 2: FID.IID
#' readIID(pfx, -1) # opt=-1: number sequence
#' readIID(pfx, -2) # opt=-2: number sequence, fixed length, decimal
#' readIID(pfx, -3) # opt=-3: number sequence, fixed length, hexidemical
#'
#' @export
readIID <- function(fam, opt=NULL)
{
    opt <- if(is.null(opt)) 1 else opt            # default option = 1
    if(identical(opt, 0L) || identical(opt, 0))   # do nothing
        iid <- NULL
    else if(is.numeric(opt) && length(opt) == 1 && opt > 0)
    {
        if(is.character(fam) && length(fam) == 1) # read FAM from file
            fam <- readFAM(fam)
        if(is.data.frame(fam) && all(c("fid", "iid") %in% names(fam)))
        {
            if(opt == 1) # IID as row name
                iid <- fam$iid
            else         # FID.IID as row name
                iid <- paste(fam$fid, fam$iid, sep='.')
        }
        else
            stop("Bad FAM:", fam)
    }
    else if(is.numeric(opt) && length(opt) == 1 && opt < 0)
    {
        if(is.character(fam) && length(fam) == 1)
            N <- nr(paste0(sub("[.](bed|fam|fam)$", "", fam), ".fam"))
        else
            N <- NROW(fam)
        if(opt == -1)
            iid <- sprintf("%d", seq(N))
        else if(opt == -2)
            iid <- .name(N, alpha=c(0:9))
        else
            iid <- .name(N, alpha=c(0:9, letters[1:6]))
    }
    else if(length(opt) > 1)
    {
        if(is.character(fam) && length(fam) == 1)
            N <- nr(paste0(sub("[.](bed|fam|fam)$", "", fam), ".fam"))
        else
            N <- NROW(fam)
        if(length(opt) != N)
            stop(gettextf("unequal num of IDs and individuals (%d and %d).", length(opt), N))
        iid <- opt
    }
    else
        stop(gettextf("fail to get variant IDs from \"%s\" with opt=%s.", fam, opt))
    iid
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
#' @param bim prefix or name of a PLINK file.
#' @return data frame of variants, loaded from BIM.
#' @examples
#' bed <- file.path(system.file("extdata", package="plinkFile"), "000.bed")
#' bim <- readBIM(bed, 20, 30)
#' bim
#' @export
readBIM <- function(bim, vfr=NULL, vto=NULL)
{
    ## fix filename
    bimFile <- paste0(sub("[.](bed|bim|fam)$", "", bim), ".bim")
    ## read BIM
    hdr <- c("chr", "vid", "cmg", "bps", "al1", "al2")
    ccs <- c("character", "character", "integer", "integer", "character", "character")
    skp <- if(is.null(vfr))  0  else vfr - 1
    nln <- if(is.null(vto)) -1  else vto - skp
    bim <- utils::read.table(bimFile, FALSE, col.names=hdr, colClasses=ccs, nrows=nln, skip=skp)
    ## convert chromosomes to integer
    within(bim, chr <- CHR[chr])
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
#' @param bim prefix or name of a PLINK file, or data frame from a BIM file.
#' @param opt option (def=1: the 2nd column in BIM file).
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
#'
#' # opt=-1: number sequence
#' vid <- readVID(pfx, -1); head(vid); tail(vid)
#'
#' # opt=-2: number sequence, fixed length, decimal
#' vid <- readVID(pfx, -2); head(vid); tail(vid)
#'
#' # opt=-3: number sequence, fixed length, hexidemical
#' vid <- readVID(pfx, -3); head(vid); tail(vid)
#' @export
readVID <- function(bim, opt=NULL, vfr=NULL, vto=NULL)
{
    opt <- if(is.null(opt)) 1 else opt            # default option = 1
    if(identical(opt, 0L) || identical(opt, 0))             # do nothing
        vid <- NULL
    else if(is.numeric(opt) && length(opt) == 1 && opt > 0) # based on BIM table
    {
        if(is.character(bim) && length(bim) == 1) # read BIM from file
            bim <- readBIM(bim, vfr, vto)
        if(is.data.frame(bim) && all(c("chr", "bps", "al1", "al2") %in% names(bim)))
        {
            if(opt == 1)
                vid <- bim$vid
            else if(opt == 2)
                vid <- with(bim, sprintf("%02d:%09d",       chr, bps))
            else
                vid <- with(bim, sprintf("%02d:%09d_%s_%s", chr, bps, al1, al2))
        }
        else
            stop("Bad BIM:", bim)
    }
    else if(is.numeric(opt) && length(opt) == 1 && opt < 0) # by variant counts
    {
        if(is.character(bim) && length(bim) == 1)
            P <- nr(paste0(sub("[.](bed|bim|fam)$", "", bim), ".bim"))
        else
            P <- NROW(bim)
        if(opt == -1)
            vid <- sprintf("%d", seq(P))
        else if(opt == -2)
            vid <- .name(P, alpha=c(0:9))
        else
            vid <- .name(P, alpha=c(0:9, letters[1:6]))
    }
    else if(length(opt) > 1)
    {
        if(is.character(bim) && length(bim) == 1)
            P <- nr(paste0(sub("[.](bed|bim|fam)$", "", bim), ".bim"))
        else
            P <- NROW(bim)
        if(length(opt) != P)
            stop(gettextf("unequal num of IDs and variants (%d and %d).", length(opt), P))
        vid <- opt
    }
    else
        stop(gettextf("fail to get variant IDs from \"%s\" with opt=%s.", bim, opt))
    vid
}

#' @describeIn bed apply a function to variants in a PLINK1 BED fileset
#'
#' The scripts in `FUN` has no side effects on environment of `scanBED`.
#'
#' @param FUN a function to process each window of variants;
#' @param ... additional argument for *`FUN`* when `scanBED` is used.
#' @param win window size - the number of variants passed to `FUN` per iteration.
#' @param iid IID as row names (def=1, see [readIID()]). 
#' @param vid VID as col names (def=1, see [readVID()]).
#' @examples
#' ## traverse genotype, apply R function without side effects
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "000")
#' ret <- scanBED(pfx, function(g)
#' {
#'     .af <- colMeans(g, na.rm=TRUE) / 2
#'     maf <- pmin(.af, 1 - .af)
#'     mis <- colSums(is.na(g)) / .N
#'     pct <- round(.w / .W * 100, 2)
#'     cbind(buf=.b, wnd=.w, idx=.i, MAF=maf, MIS=mis, PCT=pct)
#' },
#' vfr=NULL, vto=NULL, win=13, simplify=rbind, buf=2^18)
#' head(ret)
#' tail(ret)
#'
#' @export
scanBED <- function(pfx, FUN, ..., win=1, iid=1, vid=1, vfr=NULL, vto=NULL, buf=2^24, simplify=TRUE)
{
    ## PLINK file set
    bedFile <- paste0(sub("[.](bed|bim|fam)$", "", pfx), ".bed")
    famFile <- paste0(sub("[.](bed|bim|fam)$", "", pfx), ".fam")
    bimFile <- paste0(sub("[.](bed|bim|fam)$", "", pfx), ".bim")

    ## number of samples and features
    N <- nr(famFile)
    P <- nr(bimFile)

    ## get range
    vfr <- if(is.null(vfr)) 1 else min(max(vfr, 1), P) # variant-from
    vto <- if(is.null(vto)) P else min(max(vto, 1), P) # variant-to
    if(vto < vfr)
        stop(gettextf("to-variant (%d) preceeds from-variant (%d)", vto, vfr))

    ## IDs
    vid <- readVID(pfx, opt=vid, vfr=vfr, vto=vto)
    iid <- readIID(pfx, iid)
    
    ## praperation
    bpv <- (file.size(bedFile) - 3L)  %/% P # bytes per variant
    nvr <- vto - vfr + 1 # number of variants
    byt <- nvr * bpv     # number of bytes

    vpb <- max(buf %/% bpv, win)    # variants per buffer,
    vpb <- ceiling(vpb / win) * win # round up to full window
    bpc <- bpv * vpb                # bytes per buffer chunk
    wdx <- 0L                       # window index
    bdx <- 0L                       # buffer chunk index
    env <- environment(FUN)         # processing environment
    env[[".P"]] <- P
    env[[".N"]] <- N
    env[[".p"]] <- nvr
    env[[".W"]] <- as.integer(ceiling(nvr / win))
    env[[".B"]] <- as.integer(ceiling(byt / bpc))

    ## go through the binary data
    ret <- list()
    fp <- .try.open.bed(bedFile)
    seek(fp, origin="current", where=(vfr - 1L) * bpv)  # move file pointer
    byp <- 0 # bytes passed
    while (length(chk <- readBin(fp, "raw", min(bpc, byt - byp))) > 0)
    {
        bdx <- bdx + 1
        env[[".b"]] <- bdx
        ## cat("Chunk #", bdx, "\n") # process a chunk
        gmx <- dbd(chk, N)           # genotype matrix
        colnames(gmx) <- vid[seq(wdx * win + 1, len=ncol(gmx))]
        rownames(gmx) <- iid
        for(j in seq(1, ncol(gmx), win)) # go through windows
        {
            jdx <- seq(j, min(j + win - 1, ncol(gmx)))
            idx <- wdx * win + seq_along(jdx)
            wdx <- wdx + 1L
            ## cat("Window #", wdx, "\n") # process a window
            env[[".i"]] <- idx
            env[[".I"]] <- vfr + idx - 1
            env[[".w"]] <- wdx
            env[[".J"]] <- length(jdx) # actual window size
            ## process a window of variants
            ret[[wdx]] <- FUN(gmx[, jdx, drop=FALSE], ...)
        }
        byp <- byp + length(chk)
        ## cat(sprintf("%4d %3d %4d\n", bdx, ncol(gmx), length(chk)))
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
#' @param iid IID as row names (def=1, see [readIID()]). 
#' @param vid VID as col names (def=1, see [readVID()]).
#' @examples
#' ## traversing genotype, evaluate R expression with side effects
#' pfx <- file.path(system.file("extdata", package="plinkFile"), "000.bed")
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
loopBED <- function(pfx, EXP, GVR="g", win=1, iid=1, vid=1, vfr=NULL, vto=NULL, buf=2^24, simplify=TRUE)
{
    ## PLINK file set
    bedFile <- paste0(sub("[.](bed|bim|fam)$", "", pfx), ".bed")
    famFile <- paste0(sub("[.](bed|bim|fam)$", "", pfx), ".fam")
    bimFile <- paste0(sub("[.](bed|bim|fam)$", "", pfx), ".bim")
    
    ## number of samples and features
    N <- nr(famFile)
    P <- nr(bimFile)

    ## get range
    vfr <- if(is.null(vfr)) 1 else min(max(vfr, 1), P) # variant-from
    vto <- if(is.null(vto)) P else min(max(vto, 1), P) # variant-to
    if(vto < vfr)
        stop(gettextf("to-variant (%d) preceeds from-variant (%d)", vto, vfr))

    ## IDs
    vid <- readVID(pfx, opt=vid, vfr=vfr, vto=vto)
    iid <- readIID(pfx, iid)

    ## praperation
    bpv <- (file.size(bedFile) - 3L)  %/% P # bytes per variant
    nvr <- vto - vfr + 1 # number of variants
    byt <- nvr * bpv     # number of bytes

    vpb <- max(buf %/% bpv, win)    # variants per buffer,
    vpb <- ceiling(vpb / win) * win # round up to full window
    bpc <- bpv * vpb                # bytes per buffer chunk
    wdx <- 0L                       # window index
    bdx <- 0L                       # buffer chunk index
    EXP <- substitute(EXP)          # expression
    env <- parent.frame()           # processing environment
    env[[".P"]] <- P
    env[[".N"]] <- N
    env[[".p"]] <- nvr
    env[[".W"]] <- as.integer(ceiling(nvr / win))
    env[[".B"]] <- as.integer(ceiling(byt / bpc))

    ## go through the binary data
    ret <- list()
    fp <- .try.open.bed(bedFile)
    seek(fp, origin="current", where=(vfr - 1L) * bpv)  # move file pointer
    byp <- 0 # bytes passed
    while (length(chk <- readBin(fp, "raw", min(bpc, byt - byp))) > 0)
    {
        bdx <- bdx + 1
        env[[".b"]] <- bdx
        ## cat("Chunk #", bdx, "\n") # process a chunk
        gmx <- dbd(chk, N)               # genotype matrix
        colnames(gmx) <- vid[seq(wdx * win + 1, len=ncol(gmx))]
        rownames(gmx) <- iid
        for(j in seq(1, ncol(gmx), win)) # go through windows
        {
            jdx <- seq(j, min(j + win - 1, ncol(gmx)))
            idx <- wdx * win + seq_along(jdx)
            wdx <- wdx + 1L
            ## cat("Window #", wdx, "\n") # process a window
            env[[".i"]] <- idx
            env[[".I"]] <- vfr + idx - 1
            env[[".w"]] <- wdx
            env[[".J"]] <- length(jdx) # actual window size
            env[[GVR]] <- gmx[, jdx, drop=FALSE]
            ## process a window of variants
            ret[[wdx]] <- eval(EXP, env)
        }
        byp <- byp + length(chk)
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
