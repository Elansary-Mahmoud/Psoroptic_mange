map = read.table('GWAS_genotypes12.map',header=F)
ped = read.table('GWAS_genotypes12.ped',header=F)

CHRs = list(NULL)
for(i in 0:30)
{
	if(i == 0)
	{
		CHR_0 = nrow(map[map[,1] == i,])
	}
	else 
	{
		CHRs[[i]] = nrow(map[map[,1] == i,])
	}
}
CHRs = unlist(CHRs)
skip = 6 + (2*CHR_0) + 1
spaces = c(rep('    ',9),rep('   ',90),rep('  ',384))
for(i in 1:30)
{
	if(i == 1)
	{
		from = skip
		to = skip + (CHRs[i]*2) -1 
		geno = ped[,c(2,from:to)]
	}
	else
	{	
		from = skip + (sum(CHRs[1:(i-1)])*2)
		to = skip + (sum(CHRs[1:i])*2) - 1
		geno = ped[,c(2,from:to)]
	}
	geno = geno[-c(156,200),]
	geno[,1] = 1:nrow(geno)	
	geno=cbind(spaces,geno)
	write.table(geno,file=paste('/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs/QC_genotyped_files/geno_CHR_',i,'.txt',sep=''),quote=F,row.names=F,col.names=F)
}

for(i in 1:30)
{
	chr_map = cbind(1:CHRs[i],map[map[,1] == i,c(2,4)])
	chr_map[,3] = chr_map[,3] / 1000000
	write.table(chr_map,file=paste('/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs/QC_genotyped_files/map_CHR_',i,'.txt',sep=''),quote=F,row.names=F,col.names=F)
}

#####################################	Convert to Beagle Format 	####################################



beagle = read.table('GWAS_genotypes12.ped',header=F,stringsAsFactors=F)
temp = beagle[,7:ncol(beagle)]
temp = temp[-c(156,200),]
Allele1 = temp[,seq(from=1,to=ncol(temp),by=2)]
Allele2 = temp[,seq(from=2,to=ncol(temp),by=2)]
Indv = NULL
All_Indv = NULL
for(i in 1:nrow(temp))
{
	print(i)
	Indv = cbind(t(Allele1[i,]),t(Allele2[i,]))
	All_Indv = cbind(All_Indv,Indv)
}

beagle_map = read.table('GWAS_genotypes12.map',header=F)

Beagle_Geno = cbind(rep("M",nrow(beagle_map)),as.character(beagle_map[,2]),All_Indv)

sample_names = 1:nrow(temp)
colnames(Beagle_Geno) = c("I","id",unlist(lapply(1:nrow(temp), function(i) { rep(sample_names[i],2) } )))

for(j in 30:30)
{
	print(paste("BTA ",j," started",sep=""))
	start = sum(summary(as.factor(beagle_map[,1]))[1:j]) + 1
	end = sum(summary(as.factor(beagle_map[,1]))[1:(j+1)])
	#write.table(beagle_map[start:end,c(2,4:6)],file=paste("Beagle_CHR_",j,".map",sep=""),quote=F,row.names=F)	
	write.table(Beagle_Geno[start:end,],file=paste("Beagle_CHR_",j,".txt",sep=""),quote=F,row.names=F)	
}


#####################################	Convert FROM Beagle Format TO PHASEBOOK	####################################
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
spaces = c(rep('    ',9*2),rep('   ',90*2),rep('  ',384*2))
for(i in 1:30)
{
	print(i)
	#system(paste("gunzip /root/Cow/Beagle_files/CHR",i,".Beagle_CHR_",i,".txt.phased.gz",sep=""))
	Beagle = read.table(paste('CHR',i,'.Beagle_CHR_',i,'.txt.phased',sep=''),stringsAsFactor=F)
	Beagle = t(Beagle)
	Beagle = Beagle[-c(1,2),]
	Beagle = trim(Beagle)
	x=c("1","2")
	Beagle = cbind(spaces,Beagle[,1],rep(x,nrow(Beagle)/2),Beagle[,2:ncol(Beagle)])
	write.table(Beagle,file=paste('phased_CHR_',i,'.txt',sep=''),quote=F,row.names=F,col.names=F)
}




