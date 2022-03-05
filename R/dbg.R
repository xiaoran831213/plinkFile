tst1 <- function(vfr=NULL, vto=NULL)
{
    pfx <- file.path(system.file("extdata", package="plinkFile"), "000")
    ret <- scanBED(pfx, function(g)
    {
        .af <- colMeans(g, na.rm=TRUE) / 2
        maf <- pmin(.af, 1 - .af)
        mis <- colSums(is.na(g)) / .N
        pct <- round(.i / .p * 100, 2)
        cbind(buf=.b, wnd=.w, idx=.I, MAF=maf, MIS=mis, PCT=pct)
    },
    vfr=vfr, vto=vto, win=13, simplify=rbind, buf=2^18)
    ret
}


tst2 <- function(vfr=0, vto=1)
{
    pfx <- file.path(system.file("extdata", package="plinkFile"), "000")
    ret <- loopBED(pfx,
    {
        .af <- colMeans(g, na.rm=TRUE) / 2
        maf <- pmin(.af, 1 - .af)
        mis <- colSums(is.na(g)) / .N
        pct <- round(.i / .p * 100, 2)
        cbind(buf=.b, wnd=.w, idx=.I, MAF=maf, MIS=mis, PCT=pct)
    },
    GVR="g", vfr=vfr, vto=vto, win=13, simplify=rbind, buf=2^18)
    ret
}
