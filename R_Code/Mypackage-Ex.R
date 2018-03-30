pkgname <- "Mypackage"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Mypackage')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Convert_VCF_PLINK")
### * Convert_VCF_PLINK

flush(stderr()); flush(stdout())

### Name: Convert_VCF_PLINK
### Title: COnvert function
### Aliases: Convert_VCF_PLINK
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file, quality_th = 0.3, output) 
{
    Beagle = read.table(gzfile(file), stringsAsFactors = F)
    temp = cbind(unlist(lapply(Beagle[, 4], nchar)), unlist(lapply(Beagle[, 
        5], nchar)))
    nrow(temp[temp[, 1] > 1 | temp[, 2] > 1, ])
    Beagle = Beagle[-which(temp[, 1] > 1 | temp[, 2] > 1), ]
    r2 = unlist(lapply(1:nrow(Beagle), function(x) {
        strsplit(Beagle[x, 8], split = ";")[[1]][2]
    }))
    r2_value = unlist(lapply(1:nrow(Beagle), function(x) {
        strsplit(r2[x], split = "=")[[1]][2]
    }))
    Beagle = Beagle[which(as.numeric(r2_value) > quality_th), 
        ]
    phasebook = data.frame(matrix(0, nrow = 146, ncol = (nrow(Beagle))))
    count = 1
    for (i in 10:ncol(Beagle)) {
        print(count)
        x = 1
        phasebook[count, ] = ifelse(substring(Beagle[, i], 0, 
            3) == "0|0", "1 1", ifelse(substring(Beagle[, i], 
            0, 3) == "1|1", "2 2", ifelse(substring(Beagle[, 
            i], 0, 3) %in% c("0|1", "1|0"), "1 2", "error")))
        count = count + 1
    }
    Extreme = read.table("/root/Cow/PLINK_070415_0630/extremes_pheno.txt", 
        stringsAsFactors = F)
    sample_names = Extreme[, 2]
    Family_names = Extreme[, 1]
    phasebook = cbind(Family_names, sample_names, rep(0, 146), 
        rep(0, 146), rep(-9, 146), rep(-9, 146), phasebook)
    write.table(phasebook, file = paste(output, ".ped", sep = ""), 
        quote = F, row.names = F, col.names = F)
    Map = cbind(rep(11, nrow(Beagle)), paste("X", 1:nrow(Beagle), 
        sep = ""), rep(0, nrow(Beagle)), Beagle[, 2])
    write.table(Map, file = paste(output, ".map", sep = ""), 
        quote = F, row.names = F, col.names = F)
  }



cleanEx()
nameEx("Mypackage-package")
### * Mypackage-package

flush(stderr()); flush(stdout())

### Name: Mypackage-package
### Title: What the package does (short line) ~~ package title ~~
### Aliases: Mypackage-package Mypackage
### Keywords: package

### ** Examples

~~ simple examples of the most important functions ~~



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
