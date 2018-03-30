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
for(i in 10:10)
{
	system(paste("java -Xmx2000m -jar /media/HD-E1/Projects/porc/GIGA\\ lab/Norway\\ day/Beagle/RG\\ NO/beagle.jar unphased=Beagle_CHR_",i," missing=0 out=CHR",i,sep=''))
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





###########################		LD pruning	#################################

snpgdsBED2GDS('clean-inds-GWA-data.bed', 'clean-inds-GWA-data.fam','clean-inds-GWA-data.bim', "Data_v2.gds")
genofile <- openfn.gds("Data_v2.gds")
put.attr.gdsn(index.gdsn(genofile,"snp.chromosome"),"autosome.end",29)
put.attr.gdsn(index.gdsn(genofile,"snp.chromosome"),"X",30)
put.attr.gdsn(index.gdsn(genofile,"snp.chromosome"),"Y",31)
put.attr.gdsn(index.gdsn(genofile,"snp.chromosome"),"XY",32)
put.attr.gdsn(index.gdsn(genofile,"snp.chromosome"),"MT",33)
put.attr.gdsn(index.gdsn(genofile,"snp.chromosome"),"M",33)
set.seed(10000)
snpset=snpgdsLDpruning(genofile,ld.threshold=0.2)
snpset.id = unlist(snpset)
#ibd=snpgdsIBDMoM(genofile,snp.id=snpset.id)
ibd=snpgdsIBDMLE(genofile,snp.id=snpset.id, verbose=T)
ibd.coeff=snpgdsIBDSelection(ibd)
head(ibd.coeff)
write.table(ibd.coeff,file='IBD_Matrix_SNPRelate_MLE.txt',quote=F)

##################################################################################################

for(i in 1:28)
{
	print(i)
	if(i == 1) {
		new = as.numeric(system(paste("wc -l Beagle_CHR_",i,".map | cut -f1 -d' '",sep=""),intern=TRUE)) - 1
		old = 10

		system(paste("sed -i 's/HAPLO 10 ",as.character(old)," 3/HAPLO 10 ",as.character(new)," 3/g' new_Parameters_GLASCOW.txt",sep=""))
		system(paste("./GLASCOW new_Parameters_GLASCOW.txt -dbg -allowtwins > GLASCOW_CHR_",i,".log",sep=""))
	} else {
		new = as.numeric(system(paste("wc -l Beagle_CHR_",i,".map | cut -f1 -d' '",sep=""),intern=TRUE)) - 1
		old = as.numeric(system(paste("wc -l Beagle_CHR_",i-1,".map | cut -f1 -d' '",sep=""),intern=TRUE)) - 1
		system(paste("sed -i 's/ANAME Analysis_GLASCOW_CHR_",i-1,"/ANAME Analysis_GLASCOW_CHR_",i,"/g' new_Parameters_GLASCOW.txt",sep=""))
		system(paste("sed -i 's/HAPLO 10 ",as.character(old)," 3/HAPLO 10 ",as.character(new)," 3/g' new_Parameters_GLASCOW.txt",sep=""))
		system(paste("sed -i 's/GFILE GLASCOW_phased_CHR_",i-1,".txt GLASCOW_phased_CHR_",i,".txt/GFILE GLASCOW_phased_CHR_",i,".txt GLASCOW_phased_CHR_",i+1,".txt/g' new_Parameters_GLASCOW.txt",sep=""))
		system(paste("./GLASCOW new_Parameters_GLASCOW.txt -dbg -allowtwins > GLASCOW_CHR_",i,".log",sep=""))

	}
}

x = paste("Analysis_GLASCOW_CHR_",1:28,".out",sep="")
file = file.path("/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs",x)
res = lapply(file[file.exists(file)],read.table)
library('qvalue')
library('gap')
CHR = which(file.exists(file))

x = paste("Beagle_CHR_",1:28,".map",sep="")
file = file.path("/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs",x)
map = lapply(file[file.exists(file)],read.table,header=T)


for(i in 1:length(res))
{
	print(CHR[i])
	res[[i]][,8] = pgamma(res[[i]][,4],shape=res[[i]][,2],scale=res[[i]][,3],lower.tail=F)
	
	png(paste("/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs/plots/hist_CHR_",CHR[i],".png",sep=""),width=290,height=210,units="mm",res=100)
	hist(res[[i]][,8],main=paste("Histogram CHR ",CHR[i],sep=""))
	dev.off()
	
	png(paste("/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs/plots/qq_plot_CHR_",CHR[i],".png",sep=""),width=290,height=210,units="mm",res=100)
	qqunif(res[[i]][,8],main=paste("QQ Plot CHR ",CHR[i],sep=""))
	dev.off()
	
	png(paste("/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs/plots/mht_plot_CHR_",CHR[i],".png",sep=""),width=290,height=210,units="mm",res=100)
	plot(map[[CHR[i]]][,2]/1000000,-log10(res[[i]][,8]),main=paste("mht Plot CHR ",CHR[i],sep=""),xlab="position in MB",ylab="-log10(pvalue)")
	dev.off()
}


###############################################################################################
##########################		Remove Extremes			#######################
###############################################################################################

extremes = read.table("Extremes.txt",header=T)
PFILE = read.table("PFILE_GLASCOW_v2.txt")
lookup = read.table("../gcta_1.24.4/GWAS_v2.fam")
names = paste("states_",1:29,sep="")
#names = names[-c(2,4,6)]
files = file.path("/home/ulg/genan/mahmoud/Annalys/GLASCOW_analysis/HMM15",names)
states = lapply(files,read.table)
names(states) = names
extremes_FID = lookup[lookup[,1] %in% intersect(lookup[,1],extremes[,1]),2]
PFILE = PFILE[PFILE[,1] %in% extremes_FID,]
PFILE[,1] = 1:nrow(PFILE)
#write.table(PFILE,"PFILE_GLASCOW_Extremes.txt",quote=F,row.names=F,col.names=F)



FID = 1:length(extremes_FID)
FID = unlist(lapply(1:length(extremes_FID),function(k) { rep(FID[k],2) } ) )
spaces = c(rep('    ',9),rep('   ',90),rep('  ',12))
for(i in 1:length(states))
{
	print(i)
	temp=NULL
	for(j in 1:length(extremes_FID))
	{
		temp=rbind(temp, states[[i]][((extremes_FID[j]*2)-1):(extremes_FID[j]*2),])
	}
	print(all.equal(unique(temp[,1]),extremes_FID))
	temp[,1] = FID
	temp = cbind(spaces,temp)
	write.table(temp,file=paste("Extremes_analysis_HMM15/",names[i],sep=""),quote=F,col.names=F,row.names=F)
}



 
###############################################################################################
##########################		MNHT			#######################
###############################################################################################
res = read.table("Extremes_Hered_scores.out")
chr_no = read.table("../matrix.in")
colors = NULL
ticks = NULL
counter = 1
chr_names = unlist(lapply(1:29,function(i) { rep(i,chr_no[i,2]) }))
for(i in 1:nrow(chr_no))
{
	if(i %in% c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29))
	{
		colors = c(colors,rep("grey",chr_no[i,2]),rep("white",500))
		if(i == 1)
		{
			ticks[i] = chr_no[i,2]/2 
		}
		else
		{
			ticks[i] = sum(chr_no[1:(i-1),2]) + chr_no[i,2]/2 + (500*(i-1))
		}
	}
	else
	{
		colors = c(colors,rep("black",chr_no[i,2]),rep("white",500))
		ticks[i] = sum(chr_no[1:(i-1),2]) + chr_no[i,2]/2 + (500*(i-1))
	}
}
temp = cbind(res,chr_names)
insert = matrix(1,ncol=8,nrow=500)
colnames(insert) = colnames(temp)
new = NULL
for(i in 1:nrow(chr_no))
{
	chr = temp[temp$chr_names == i,]
	new =  rbind(new,chr,insert)
	
}

new[,1] = 1:nrow(new)
png("C_MHT_EXTREMES_CASE_CONTROL.png",width=290,height=210,units="mm",res=150)
plot(new[,1],-log10(new[,5]),col=colors,xaxt = "n",pch=19,xlab="Chromosomes",ylab="-log10(p-value)")
axis(side=1 ,at=ticks,labels = 1:29,cex=3)
dev.off()

###############################################################################################
########################	Save DATA			##############################
###############################################################################################
library("gap")
png("mht_plot_Extreme_Hered_scores.png",width=290,height=210,units="mm",res=150)
plot(new[,1],-log10(new[,5]),col=colors,xaxt = "n",pch=19,xlab="chromosome",ylab="-log10(pvalue)" )#,ylab=expression(-log[10](pvalue)))
axis(side=1 ,at=ticks,labels = 1:29,cex=3)
dev.off()	
png("qq_plot_Extreme_Hered_scores.png",width=290,height=210,units="mm",res=150)
qqunif(new[,5])	
dev.off()
png("hist_plot_Extreme_Hered_scores.png",width=290,height=210,units="mm",res=100)	
hist(new[,5])
dev.off()
rm(list=ls())
###############################################################################################
###############################################################################################

for(i in 1:1)
{
	print(i)
	png(paste("mht_CHR_",i,".png",sep=""),width=290,height=210,units="mm",res=150)
	if(i == 1){
		from = 1
	} else {
		from = sum(chr_no[1:(i-1),2])
	}
	to = sum(chr_no[1:i,2])
	plot(1:length(from:to),-log10(res[from:to,5]))
	dev.off()
}





for(i in 1:nrow(chr_no))
{
	if(i %in% c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29))
	{
		colors = c(colors,rep("grey",chr_no[i,2]),rep("white",300))
		if(i == 1)
		{
			ticks[i] = chr_no[i,2]/2
		}
		else
		{
			ticks[i] = sum(chr_no[1:(i-1),2]) + chr_no[i,2]/2 + 300
		}
	}
	else
	{
		colors = c(colors,rep("black",chr_no[i,2]),rep("white",300))
		ticks[i] = sum(chr_no[1:(i-1),2]) + chr_no[i,2]/2 + 300
	}
}
