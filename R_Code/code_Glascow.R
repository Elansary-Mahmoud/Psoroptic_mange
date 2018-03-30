setwd('/media/HD-E1/Projects/Cow/Annalys/plink-1.07-dos')
snpgdsBED2GDS('Annelies_4900.bed', 'Annelies_4900.fam','Annelies_4900.bim', "Data.gds")
genofile <- openfn.gds("Data.gds")
snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))
RV <- snpgdsPCA(genofile, num.thread=2)
plot(RV$eigenvect[,2], RV$eigenvect[,1], xlab="PC 2", ylab="PC 1", pch=19)



closefn.gds(genofile)




fam = read.table("Annelies_4900.fam")
mat <- snpgdsGetGeno(genofile)
rownames(mat) = fam[,1]
genotypes_labels = c("AA","AB","BB")


system('/home/ulg/genan/mahmoud/plink-1.07-x86_64/plink --noweb --file /home/ulg/genan/mahmoud/Annalys/COW_Analysis/plink --recode12 --out /home/ulg/genan/mahmoud/Annalys/COW_Analysis/Analysis_12')
x = system("awk {'print NF'} Analysis_12.ped | head -1",intern=TRUE)
system(paste("cut -d' ' -f 1,7-",x," Analysis_12.ped > test_12.ped",sep=""))

PFILE = read.table('pheno.txt',header=T,stringsAsFactor=F)
PFILE = data.frame(PFILE[,1],RV$eigenvect[,1:5],PFILE[,3])

MapFILE = read.table('20140630_490s2.map')




beagle = read.table('Beagle_genofile.ped',stringsAsFactors=F)
temp = beagle[,7:ncol(beagle)]

Indv = NULL
All_Indv = NULL
for(i in 1:nrow(beagle))
{
	print(i)
	Indv = cbind(t(Allele1[i,]),t(Allele2[i,]))
	All_Indv = cbind(All_Indv,Indv)
}

beagle_map = read.table('Beagle_genofile.bim')

Beagle_Geno = cbind(rep("M",nrow(beagle_map)),as.character(beagle_map[,2]),All_Indv)

colnames(Beagle_Geno) = c("I","id",unlist(lapply(1:nrow(beagle), function(i) { rep(beagle[i,2],2) } )))

for(j in 1:29)
{
	print(paste("BTA ",j," started",sep=""))
	start = sum(summary(as.factor(beagle_map[,1]))[1:j]) + 1
	end = sum(summary(as.factor(beagle_map[,1]))[1:(j+1)])
	write.table(beagle_map[start:end,c(2,4:6)],file=paste("Beagle_CHR_",j,".map",sep=""),quote=F,row.names=F)	
	write.table(Beagle_Geno[start:end,],file=paste("Beagle_CHR_",j,".txt",sep=""),quote=F,row.names=F)	
}


  Beagle_CHR_3
for(i in 1:1)
{
	system(paste("java -Xmx2000m -jar /media/HD-E1/Projects/porc/GIGA\\ lab/Norway\\ day/Beagle/RG\\ NO/beagle.jar unphased=Beagle_CHR_",i,".txt missing=0 out=CHR",i,sep=''))
}


Indv = list(NULL)
for(chr in 10:29)
{
	print(paste("CHR ",chr,sep=""))
	phased = NULL
	phased = read.table(gzfile(paste("CHR",chr,".Beagle_CHR_",chr,".phased.gz",sep="")), header=T,check.names=F,stringsAsFactor=F)
	for(i in 1:nrow(phased))
	{
		print(i)
		SNP1 = NULL
		SNP2 = NULL
		for(j in 3:ncol(phased))
		{
			if(!is.null(SNP1) & !is.null(SNP2)) {
				break;
			}
			else {
				if(j == 3) {
					SNP1 = phased[i,j]
				} else {
					if(phased[i,j] != SNP1){
						SNP2 = phased[i,j]
					}
				}
			}
		}
		phased[i,which(phased[i,] == SNP1)] = rep(1,length(phased[i,which(phased[i,] == SNP1)]))
		phased[i,which(phased[i,] == SNP2)] = rep(2,length(phased[i,which(phased[i,] == SNP2)]))
		
	}
	phased = t(phased[3:ncol(phased)])
	Indv[[chr]] = rownames(phased)[seq(from=1,to=nrow(phased),by=2)]
	FID = 1:(nrow(phased)/2)
	FID = unlist(lapply(1:(nrow(phased)/2),function(k) { rep(FID[k],2) } ) )
	phased = cbind(FID,phased)
	write.table(phased,file=paste("GLASCOW_phased_CHR_",chr,".txt",sep=""),quote=F,row.names=F,col.names=F)
}

count = 1
for(i in 1:29)
 {
 	print(i)
 for(j in i:29)
 {
 	print(j)
	x[count] = all.equal(Indv[[1]],Indv[[2]])
	count = count + 1
 }
 }

for(chr in 2:29)
{
	print(chr)
	phase = read.table(paste("GLASCOW_phased_CHR_",chr,".txt",sep=""), header=F,check.names=F,stringsAsFactor=F)	
	x = c(1,2)
	x = rep(x,nrow(phase)/2)
	temp = data.frame(phase[,1],x,phase[,2:ncol(phase)])
	write.table(temp,file=paste("GLASCOW_phased_CHR_",chr,".txt",sep=""),quote=F,col.names=F,row.names=F)
}






