tst1 <- function()
{
    pfx <- file.path(system.file("extdata", package="plinkFile"), "000")
    ret <- scanBED(pfx, function(g)
    {
        .af <- colMeans(g, na.rm=TRUE) / 2
        maf <- pmin(.af, 1 - .af)
        mis <- colSums(is.na(g)) / .N
        wnd <- .w
        pct <- round(.w / .W * 100, 2)
        cbind(buf=.b, wnd=.w, idx=.i, MAF=maf, MIS=mis, PCT=pct)
    },
    win=13, simplify=rbind, buf=2^18)
    ret
}

tst2 <- function()
{
    pfx <- file.path(system.file("extdata", package="plinkFile"), "000")
    ret <- list() # use side effect to keep the result of each window.
    loopBED(pfx,
    {
        af <- colMeans(gt, na.rm=TRUE) / 2
        sg <- af * (1 - af)
        ret[[.w]] <- cbind(wnd=.w, alf=af, var=sg)
    },
    win=13, GVR="gt", vid=3, buf=2^18)
    ret <- do.call(rbind, ret)
    ret
}
