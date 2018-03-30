pkgname <- "CowAnalysis"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('CowAnalysis')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Convert_BEAGLE_PhasedFile_PHASEBOOK")
### * Convert_BEAGLE_PhasedFile_PHASEBOOK

flush(stderr()); flush(stdout())

### Name: Convert_BEAGLE_PhasedFile_PHASEBOOK
### Title: Convert_BEAGLE_PhasedFile_PHASEBOOK
### Aliases: Convert_BEAGLE_PhasedFile_PHASEBOOK
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ped_file, bim_file, chromosome_num = 30, output_folder) 
{
    beagle = read.table(ped_file, header = F, stringsAsFactors = F)
    temp = beagle[, 7:ncol(beagle)]
    Allele1 = temp[, seq(from = 1, to = ncol(temp), by = 2)]
    Allele2 = temp[, seq(from = 2, to = ncol(temp), by = 2)]
    All_Indv = list(NULL)
    Allele1 = t(Allele1)
    Allele2 = t(Allele2)
    All_Indv = lapply(1:ncol(Allele1), function(x) {
        cbind(Allele1[, x], Allele2[, x])
    })
    All_Indv = do.call("cbind", All_Indv)
    beagle_map = read.table(bim_file, header = F)
    Beagle_Geno = cbind(rep("M", nrow(beagle_map)), as.character(beagle_map[, 
        2]), All_Indv)
    sample_names = 1:nrow(temp)
    colnames(Beagle_Geno) = c("I", "id", unlist(lapply(1:nrow(temp), 
        function(i) {
            rep(sample_names[i], 2)
        })))
    for (j in 1:chromosome_num) {
        print(paste("BTA ", j, " started", sep = ""))
        start = sum(summary(as.factor(beagle_map[, 1]))[1:j]) + 
            1
        end = sum(summary(as.factor(beagle_map[, 1]))[1:(j + 
            1)])
        write.table(beagle_map[start:end, c(2, 4:6)], file = paste(output_folder, 
            "Beagle_CHR_", j, ".map", sep = ""), quote = F, row.names = F, 
            col.names = F)
        write.table(Beagle_Geno[start:end, ], file = paste(output_folder, 
            "Beagle_CHR_", j, ".txt", sep = ""), quote = F, row.names = F)
    }
  }



cleanEx()
nameEx("Convert_PLINK_BEAGLE")
### * Convert_PLINK_BEAGLE

flush(stderr()); flush(stdout())

### Name: Convert_PLINK_BEAGLE
### Title: Convert_PLINK_BEAGLE
### Aliases: Convert_PLINK_BEAGLE
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ped_file, bim_file, chromosome_num = 30, output_folder) 
{
    beagle = read.table(ped_file, header = F, stringsAsFactors = F)
    temp = beagle[, 7:ncol(beagle)]
    Allele1 = temp[, seq(from = 1, to = ncol(temp), by = 2)]
    Allele2 = temp[, seq(from = 2, to = ncol(temp), by = 2)]
    All_Indv = list(NULL)
    Allele1 = t(Allele1)
    Allele2 = t(Allele2)
    All_Indv = lapply(1:ncol(Allele1), function(x) {
        cbind(Allele1[, x], Allele2[, x])
    })
    All_Indv = do.call("cbind", All_Indv)
    beagle_map = read.table(bim_file, header = F)
    Beagle_Geno = cbind(rep("M", nrow(beagle_map)), as.character(beagle_map[, 
        2]), All_Indv)
    sample_names = 1:nrow(temp)
    colnames(Beagle_Geno) = c("I", "id", unlist(lapply(1:nrow(temp), 
        function(i) {
            rep(sample_names[i], 2)
        })))
    for (j in 1:chromosome_num) {
        print(paste("BTA ", j, " started", sep = ""))
        start = sum(summary(as.factor(beagle_map[, 1]))[1:j]) + 
            1
        end = sum(summary(as.factor(beagle_map[, 1]))[1:(j + 
            1)])
        write.table(beagle_map[start:end, c(2, 4:6)], file = paste(output_folder, 
            "Beagle_CHR_", j, ".map", sep = ""), quote = F, row.names = F, 
            col.names = F)
        write.table(Beagle_Geno[start:end, ], file = paste(output_folder, 
            "Beagle_CHR_", j, ".txt", sep = ""), quote = F, row.names = F)
    }
  }



cleanEx()
nameEx("Convert_PLINK_PHASEBOOK_unphased")
### * Convert_PLINK_PHASEBOOK_unphased

flush(stderr()); flush(stdout())

### Name: Convert_PLINK_PHASEBOOK_unphased
### Title: Convert_PLINK_PHASEBOOK_unphased
### Aliases: Convert_PLINK_PHASEBOOK_unphased
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file, chromosome_num = 30, output_folder) 
{
    map = read.table(paste(file, ".map", sep = ""), header = F)
    ped = read.table(paste(file, ".ped", sep = ""), header = F)
    CHRs = list(NULL)
    CHR_0 = NULL
    for (i in 0:chromosome_num) {
        if (i == 0) {
            CHR_0 = nrow(map[map[, 1] == i, ])
        }
        else {
            CHRs[[i]] = nrow(map[map[, 1] == i, ])
        }
    }
    CHRs = unlist(CHRs)
    skip = 6 + (2 * CHR_0) + 1
    spaces = c(rep("    ", 9), rep("   ", 90), rep("  ", (nrow(ped) - 
        99)))
    geno = NULL
    for (i in 1:chromosome_num) {
        if (i == 1) {
            from = skip
            to = skip + (CHRs[i] * 2) - 1
            geno = ped[, c(2, from:to)]
        }
        else {
            from = skip + (sum(CHRs[1:(i - 1)]) * 2)
            to = skip + (sum(CHRs[1:i]) * 2) - 1
            geno = ped[, c(2, from:to)]
        }
        geno[, 1] = 1:nrow(geno)
        geno = cbind(spaces, geno)
        cat("Writing Chromosome number ", i, " to folder\n")
        write.table(geno, file = paste(output_folder, "geno_CHR_", 
            i, ".txt", sep = ""), quote = F, row.names = F, col.names = F)
    }
    for (i in 1:30) {
        chr_map = cbind(1:CHRs[i], map[map[, 1] == i, c(2, 4)])
        chr_map[, 3] = chr_map[, 3]/10^6
        cat("Writing Map file of Chromosome number ", i, " to folder\n")
        write.table(chr_map, file = paste(output_folder, "map_CHR_", 
            i, ".txt", sep = ""), quote = F, row.names = F, col.names = F)
    }
  }



cleanEx()
nameEx("Convert_VCF_PHASEBOOK")
### * Convert_VCF_PHASEBOOK

flush(stderr()); flush(stdout())

### Name: Convert_VCF_PHASEBOOK
### Title: Convert_VCF_PHASEBOOK
### Aliases: Convert_VCF_PHASEBOOK
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
    Beagle = Beagle[which(as.numeric(r2_value) > quality), ]
    phasebook = data.frame(matrix(0, nrow = 146 * 2, ncol = nrow(Beagle) + 
        2))
    count = 1
    exclude = NULL
    ex_index = 1
    for (i in 10:ncol(Beagle)) {
        print(count)
        for (j in 1:nrow(Beagle)) {
            print(j)
            if (strsplit(Beagle[j, i], split = ":")[[1]][1] == 
                "0|0") {
                phasebook[count, (j + 2)] = "1"
                phasebook[count + 1, (j + 2)] = "1"
            }
            else if (strsplit(Beagle[j, i], split = ":")[[1]][1] == 
                "0|1") {
                phasebook[count, (j + 2)] = "1"
                phasebook[count + 1, (j + 2)] = "2"
            }
            else if (strsplit(Beagle[j, i], split = ":")[[1]][1] == 
                "1|0") {
                phasebook[count, (j + 2)] = "2"
                phasebook[count + 1, (j + 2)] = "1"
            }
            else if (strsplit(Beagle[j, i], split = ":")[[1]][1] == 
                "1|1") {
                phasebook[count, (j + 2)] = "2"
                phasebook[count + 1, (j + 2)] = "2"
            }
            else {
                print("Error")
            }
        }
        count = count + 2
    }
    sample_names = 1:146
    phasebook[, 1] = unlist(lapply(1:146, function(i) {
        rep(sample_names[i], 2)
    }))
    spaces = c(rep("    ", 9 * 2), rep("   ", 90 * 2), rep("  ", 
        (146 - 99) * 2))
    x = c("1", "2")
    phasebook = cbind(spaces, phasebook[, 1], rep(x, nrow(phasebook)/2), 
        phasebook[, 3:ncol(phasebook)])
    write.table(phasebook, file = paste(output, sep = ""), quote = F, 
        row.names = F, col.names = F)
  }



cleanEx()
nameEx("Convert_VCF_PLINK")
### * Convert_VCF_PLINK

flush(stderr()); flush(stdout())

### Name: Convert_VCF_PLINK
### Title: Convert_VCF_PLINK
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
nameEx("mht_plot_plink")
### * mht_plot_plink

flush(stderr()); flush(stdout())

### Name: mht_plot_plink
### Title: mht_plot_plink
### Aliases: mht_plot_plink
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file, autosome_num = 29, output, Disease_name = "Unknown") 
{
    res = read.table(file, header = T)
    res = res[res$TEST == "ADD", ]
    res = na.omit(res)
    res = res[res$CHR != 0, ]
    res = res[res$CHR != (autosome_num + 1), ]
    chr_no = cbind(1:autosome_num, as.numeric(summary(as.factor(res$CHR))))
    colors = NULL
    ticks = NULL
    for (i in 1:nrow(chr_no)) {
        if (i %in% seq(1, autosome_num, 2)) {
            colors = c(colors, rep("grey", chr_no[i, 2]))
            if (i == 1) {
                ticks[i] = chr_no[i, 2]/2
            }
            else {
                ticks[i] = sum(chr_no[1:(i - 1), 2]) + chr_no[i, 
                  2]/2
            }
        }
        else {
            colors = c(colors, rep("black", chr_no[i, 2]))
            ticks[i] = sum(chr_no[1:(i - 1), 2]) + chr_no[i, 
                2]/2
        }
    }
    png(paste("mht_", output, ".png", sep = ""), width = 290, 
        height = 210, units = "mm", res = 150)
    plot(1:nrow(res), -log10(res$P), col = colors, xaxt = "n", 
        pch = 19, xlab = "Chromosomes", ylab = "-log10(pvalue)", 
        main = paste("MHT plot for ", Disease_name, sep = ""))
    axis(side = 1, at = ticks, labels = 1:29, cex = 3)
    dev.off()
    png(paste("qq_", output, ".png", sep = ""), width = 290, 
        height = 210, units = "mm", res = 150)
    qqunif(res$P, main = paste("qq plot for ", Disease_name, 
        sep = ""))
    dev.off()
    png(paste("hist_", output, ".png", sep = ""), width = 290, 
        height = 210, units = "mm", res = 150)
    hist(res$P, main = paste("Histogram for ", Disease_name, 
        sep = ""))
    dev.off()
  }



cleanEx()
nameEx("mht_plot_plink_chr")
### * mht_plot_plink_chr

flush(stderr()); flush(stdout())

### Name: mht_plot_plink_chr
### Title: mht_plot_plink_chr
### Aliases: mht_plot_plink_chr
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file, chr_num, output, Disease_name = "Unknown") 
{
    res = read.table(file, header = T)
    res = res[res$TEST == "ADD", ]
    res = na.omit(res)
    png(paste("mht_", output, "_chr_", chr_num, ".png", sep = ""), 
        width = 290, height = 210, units = "mm", res = 150)
    plot(res[res$CHR == chr_num, ]$BP/10^6, -log10(res[res$CHR == 
        chr_num, ]$P), pch = 19, xlab = "position in MB", ylab = "-log10(pvalue)", 
        main = paste("MHT plot for CHR ", chr_num, " for", Disease_name, 
            sep = ""))
    dev.off()
  }



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
