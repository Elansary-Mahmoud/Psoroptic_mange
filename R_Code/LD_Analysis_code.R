library('SNPRelate')
setwd('I:\\Projects\\porc\\Plots\\BE_CR\\')
snpgdsBED2GDS('Data.bed', 'Data.fam', 'Data.bim', "BE_CR.gds")
genofile <- openfn.gds("BE_CR.gds")
snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))[(54492 - 20):(54492 + 20)]
rem = c("CASI0005728","ASGA0034469","H3GA0051449","ASGA0080725")
snpset = setdiff(snpset,rem)

L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="r", slide=-1)
colnames(L1$LD) = snpset
rownames(L1$LD) = snpset
write.table((L1$LD),file='BE_CR_r.txt')

L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="dprime", slide=-1)
colnames(L1$LD) = snpset
rownames(L1$LD) = snpset
write.table(L1$LD,file='BE_CR_dprime.txt')

L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="r", slide=-1)
colnames(L1$LD) = snpset
rownames(L1$LD) = snpset
write.table((L1$LD^2),file='BE_CR_r2.txt')


L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="r", slide=-1)
rownames(L1$LD) = snpset
pdf('BE_CR_r2.pdf')
levelplot(t(L1$LD^2), col.regions = terrain.colors)
dev.off()

L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="dprime", slide=-1)
rownames(L1$LD) = snpset
pdf('BE_CR_dprime.pdf')
levelplot(t(L1$LD), col.regions = terrain.colors)
dev.off()


L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="r", slide=-1)
rownames(L1$LD) = snpset
pdf('BE_CR_r.pdf')
levelplot(t(L1$LD), col.regions = terrain.colors)
dev.off()

########################################################

L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="r", slide=-1)
rownames(L1$LD) = snpset
temp = L1$LD
# positions = which(is.na(temp),arr.in=T)
# for(i in 1:nrow(positions))
# {
	# temp[positions[i,1],positions[i,2]] = 0
# }
pdf('BE_CR_r2_nonNA.pdf')
levelplot(na.omit(t(temp^2)), col.regions = terrain.colors)
dev.off()
temp=NULL

L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="dprime", slide=-1)
rownames(L1$LD) = snpset
temp = L1$LD
# positions = which(is.na(temp),arr.in=T)
# for(i in 1:nrow(positions))
# {
	# temp[positions[i,1],positions[i,2]] = 0
# }
pdf('BE_CR_dprime_nonNA.pdf')
levelplot(t(temp), col.regions = terrain.colors)
dev.off()
temp=NULL

L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="r", slide=-1)
rownames(L1$LD) = snplons))
# {
	# temp[positions[i,1],positions[i,2]] = 0
# }
pdf('BE_CR_r_nonNA.pdf')
levelplot(t(temp), col.regions = terrain.colors)
dev.off()


library('SNPRelate')
setwd('I:/Projects/porc/Plots/CR BE + NO/NO/')
snpgdsBED2GDS('Data.bed', 'Data.fam', 'Data.bim', "NO_CR.gds")
sex = read.table('Data.fam')
genofile <- openfn.gds("NO_CR.gds")
snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))
g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(-1,-1))
L1 = snpgdsLDpair(sex[,5],g[,1],method="r")
L2 = snpgdsLDpair(sex[,5],g[,2],method="r")
L3 = snpgdsLDpair(sex[,5],g[,3],method="r")
L4 = snpgdsLDpair(sex[,5],g[,4],method="r")
L1$ld^2
L2$ld^2
L3$ld^2
L4$ld^2
closefn.gds(genofile)
rm(list=ls())

Results = data.frame(NULL)
 for(i in 1:length(Marker1))
 {
    if(Marker1[i] == 0){
	 Results[i,1] = 1
	 Results[i,2] = 1
	} else if (Marker1[i] == 1){
		Results[i,1] = 1
		Results[i,2] = 2
	} else if (Marker1[i] == 2){
		Results[i,1] = 2
		Results[i,2] = 2
	} else {
		Results[i,1] = 0
		Results[i,2] = 0
	}
}

 for(i in 1:nrow(sex))
 {
    if(sex[i,5] == 1){
	 Results[i,3] = 1
	 Results[i,4] = 2
	} else if (sex[i,5] == 2){
		Results[i,3] = 1
		Results[i,4] = 1
	} else {
		print("Error")
	}
}
