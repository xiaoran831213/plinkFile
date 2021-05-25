#' travers variants in a PLINK1 BED fileset
#'
#' @description
#' Sequentially visits variants in a PLINK1  BED fileset with a stepping window
#' matrix, and process each window matrix  with user scripts either in function
#' or expression form, meant for data to big to fit in the memory.
#'
#' To read the entire BED into a R matrix, use `[readBED]()` instead.
#'
#' @section BED PLINK1 Binary Pedigree fileset:
#' A popular format to store biallelic dosage genotype, with three files,
#' * _pfx_.fam: text table for `N` individuals, detailed in [readFAM];
#' * _pfx_.bim: text table for `P` variants, detailed in [readBIM];
#' * _pfx_.bed: transposed genotype matrix (`P` x `N` ) in binary format.
#' 
#' The triplets are commonly referred by the shared prefix (`pfx`), e.g., the X
#' chromosome represented by "chrX.bed", "chrX.fam", and "chrX.bim" are refered
#' by `"chrX"`.
#'
#' The binary file "_pfx_.bed" represent each dosage value with two bits - just
#' enough to encode all four possiblities: 0, 1, or 2 alleles, or missing.
#'
#' The number of variants (`P`) and samples (`N`) equals to the number of lines
#' in text file "_pfx_.bim" and "_pfx_.fam", respectively.
#'
#' For the detailed specification of PLINK1 BED genotype format, see the lagecy
#' PLINK v1.07 page at: \\
#' <https://zzz.bwh.harvard.edu/plink/binary.shtml>. \\
#' For the modern use and management of PLINK1 BED, see the PLINK v1.9 page: \\
#' <https://www.cog-genomics.org/plink/1.9/input#bed>.
#'
#' @section detailed arguments:
#' * `win`: visiting window size.
#' 
#'   the number of variants per window, that  is, the number of columns in each
#'   window matrix passed to the user script.
#'
#'   For example, a size one window means  the user script will be dealing with
#'   only one variant at  a time, received from in a matrix  of a single column
#'   -- a manner similar to  genome wide association analysis (GWAS).  However,
#'   a larger, multi-variant window coupled with R language's vector and matrix
#'   syntax can significantly boost efficiency.
#'
#'   The default size is 1000 variants / columns per window.
#'
#' * `buf`: buffer size in bytes
#'
#'   a large buffer reduces the frequency of hard disk visits when traversing a
#'   PLINK1 BED file, which in turn reduces non-computation overhead.
#'
#'   The default size is `2^24` bytes, or 16 MB.
#'
#' * `simplify`:
#'
#'   when FALSE: resuts  of user script processing each window  of variants are
#'   returned in a list;
#'
#'   when TRUE,  use `simplify2array` to put  the results into an  array, if it
#'   fails, fallback and return a list.
#'
#'   when a function is specified, it is then used to simplify the results, if
#'   an execption is thrown, fallback and return a list.
#'
#'   e.g., the window script returns a  data frame of estimate, standard error,
#'   t-statistic, and p-value  for each variant, `simplify =  rbind` to combine
#'   results of all windows into one data frame of `P` rows and four columns of
#'   statistics.
#' 
#' @section genotype context:
#' context infomation  such the number of  variants and samples are  updated in
#' the window processing environment to ease user scripting, which includes:
#' 
#' - `.i`: indies of variants in the current visiting window;
#' - `.p`: number of variants in the current visiting window.
#' - `.P`: total number of variants;
#' - `.w`: index of the current window;
#' - `.W`: total number of windows to go through;
#' - `.N`: number of individuals.
#' - `.b`: index of the current buffer.
#' - `.B`: number of buffers to be swapped.
#'
#' e.g. (1) print percentage progress with `print(.w / .W * 100)`; \\
#' e.g. (2) use `inf <- readBIM(pfx)` to  read the table of variants before the
#' window visits,  later use `inf[.i,  ]` to  access meta-data for  variants in
#' each window.
#' 
#' @name bed
#' @param pfx prefix of PLINK BED.
#' @param win reading window size (def=100 variants per window)
#' @param buf buffer size in byptes (def=2^24, or 16 MB).
#' @param simplify try simplifying the results  into an array, or leave them in
#'     a list, or specify a function to simplify the said list.
#' @return results of all windows processed by the user script.
#' @seealso {`[readBED]`}
NULL
