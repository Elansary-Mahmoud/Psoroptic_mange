dhfxhdesetwd('/media/HD-E1/Projects/Cow/Annalys/plink-1.07-dos')
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


Allele1 = temp[,seq(from=1,to=ncol(temp),by=2)] 
Allele2 = temp[,seq(from=2,to=ncol(temp),by=2)]        
Indv = NULL
All_Indv = NULL
for(i in 1:nrow(beagle))
{
	print(i)
	Indv = cbind(t(Allele1[i,]),t(Allele2[i,]))
	All_Indv = cbind(All_Indv,Indv)
}

beagle_map = read.table("../PSOROVIS_clean_final_12.bim")

Beagle_Geno = cbind(rep("M",nrow(beagle_map)),as.character(beagle_map[,2]),All_Indv)

change_Names = cbind(1:nrow(beagle),beagle[,2])
write.table(change_Names,"lookup_IId_phasebook_ids.txt",quote=F,row.names=F,col.names=F)
colnames(Beagle_Geno) = c("I","id",unlist(lapply(1:nrow(beagle), function(i) { rep(change_Names[i,1],2) } )))

#beagle_map = beagle_map[beagle_map[,1] != 0,]
#Beagle_Geno = Beagle_Geno[-c(1:414),]

setwd("Beagle/")
for(j in 1:29)
{
	print(paste("BTA ",j," started",sep=""))
	start = sum(summary(as.factor(beagle_map[,1]))[1:j]) + 1
	end = sum(summary(as.factor(beagle_map[,1]))[1:(j+1)])
        output = cbind(1:length(start:end), as.character(beagle_map[start:end,2]) , beagle_map[start:end,4]/(10^6))
	write.table( output,file=paste("Beagle_CHR_",j,".map",sep=""),quote=F,row.names=F,col.names=F)	
	#write.table(Beagle_Geno[start:end,],file=paste("Beagle_CHR_",j,".txt",sep=""),quote=F,row.names=F)	
}


for(i in 2:29)
{
	system(paste("java -Xmx2000m -jar /media/Mahmoud_HD/Projects/porc/GIGA\\ lab/Norway\\ day/Beagle/RG\\ NO/beagle.jar unphased=Beagle_CHR_",i,".txt missing=0 out=CHR",i,sep=''))
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





###########################		LD proning	#################################

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
names = names[-c(2,4,6)]
files = file.path("/home/ulg/genan/mahmoud/Annalys/GLASCOW_analysis",names)
states = lapply(files,read.table)
names(states) = names
extremes_FID = lookup[lookup[,1] %in% intersect(lookup[,1],extremes[,1]),2]
PFILE = PFILE[PFILE[,1] %in% extremes_FID,]
PFILE[,1] = 1:nrow(PFILE)
write.table(PFILE,"PFILE_GLASCOW_Extremes.txt",quote=F,row.names=F,col.names=F)



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
	temp[,1] = FID
	temp = cbind(spaces,temp)
	write.table(temp,file=paste("Extremes_analysis/",names[i],sep=""),quote=F,col.names=F,row.names=F)
}



 
###############################################################################################
##########################		MNHT			#######################
###############################################################################################
res = read.table("Testing_GLASCOW_Extremes_v2.out")

res = read.table("seq_Scores.assoc.linear",header=T)
res = res[res$TEST == "ADD",]
res = na.omit(res)
res = res[res$CHR !=0,]
res = res[res$CHR !=30,]
chr_no = cbind(1:29,as.numeric(summary(as.factor(res$CHR))))
#chr_no = read.table("../matrix.in")
chr_no = cbind(1:29,as.numeric(system("for i in `seq 1 29`; do wc -l ~/Cow_10Apr2015/Map/map_CHR_${i}.txt; done | awk '{print $1}'",intern=TRUE)))
colors = NULL
ticks = NULL
for(i in 1:nrow(chr_no))
{
	if(i %in% c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29))
	{
		colors = c(colors,rep("grey",chr_no[i,2]))
		if(i == 1)
		{
			ticks[i] = chr_no[i,2]/2
		}
		else
		{
			ticks[i] = sum(chr_no[1:(i-1),2]) + chr_no[i,2]/2
		}
	}
	else
	{
		colors = c(colors,rep("black",chr_no[i,2]))
		ticks[i] = sum(chr_no[1:(i-1),2]) + chr_no[i,2]/2
	}
}

png("mht_plot_plink_scores_extremes_chr11_seq_imputed.png",width=290,height=210,units="mm",res=150)
plot(1:nrow(res),-log10(res$P),col=colors,xaxt = "n",pch=19,xlab="Chromosomes",ylab="-log10(pvalue)",main="MHT for Extremes (146) using Binary values")
axis(side=1 ,at=ticks,labels = 1:29,cex=3)
dev.off()	
png("qq_plot_plink_scores_extremes_chr11_seq_imputed.png",width=290,height=210,units="mm",res=150)
qqunif(res$P)	
dev.off()
png("hist_plot_v2.png",width=290,height=210,units="mm",res=100)	
hist(res[,5])
dev.off()

res = read.table("seq_Visit4.assoc.linear",header=T)
png("mht_plot_plink_visit4_extremes_chr11_seq_imputed.png",width=290,height=210,units="mm",res=150)
plot(res$BP/10^6,-log10(res$P),pch=19,xlab="Position in MB",ylab="-log10(pvalue)",main="CHR 11 imputed region from seq using Visit 4 values")
dev.off()
png("qq_plot_plink_visit4_extremes_chr11_seq_imputed.png",width=290,height=210,units="mm",res=150)
qqunif(res$P)
dev.off()


Map = list()
Map = lapply(paste("~/Cow_10Apr2015/Map/map_CHR_",1:29,".txt",sep=""),read.table)
for(i in )
{
	print(i)
	png(paste("mht_CHR_",i,".png",sep=""),width=290,height=210,units="mm",res=150)
	if(i == 1){
		from = 1
	} else {
		from = sum(chr_no[1:(i-1),2]) + 1
	}
	to = sum(chr_no[1:i,2])
	plot(Map[[i]][,3],-log10(res[from:to,5]),main=paste("Chr number ",i,sep=""),xlab="position in BP",ylab="-log10(pvalue)",pch=19,cex=1.1)
	dev.off()
}



pheno_extra = read.xls("~/Cow/PLINK_070415_0630/PhenoExtra.xls",sheet=1,header=T,na.strings=c(""),colClass=rep("character",5))
rownames(pheno_extra)=as.character(pheno_extra[,1])
pheno_extra = pheno_extra[as.character(Alternative[,2]),]
head(pheno_extra)
head(Alternative)
all.equal(as.character(Alternative[,2]),as.character(pheno_extra[,1]))
Alternative = cbind(Alternative,pheno_extra[,2:5])





###################################	imuptation of sequencing data 	#########################
#################################################################################################

system(paste("java -Xmx2000m -jar /media/Mahmoud_HD/Projects/porc/GIGA\\ lab/Norway\\ day/Beagle/RG\\ NO/beagle.jar gt=chr11_peak.vcf ref=belgian_blue.vcf missing=0 out=chr11_peak_700k",sep=''))



index = (unlist(lapply(paste("^",as.character(inverse[,1]),"$",sep=""),grep,x=as.character(map[,2])))*2) + 6 - 1
map_index = unlist(lapply(paste("^",as.character(inverse[,1]),"$",sep=""),grep,x=as.character(map[,2])))
count = 1

for(i in index)
{
	print(count)
	
	temp = ped[,i:(i+1)]
	temp[temp[,1] ==  bim[map_index[count],6],1] = "x"
	temp[temp[,1] ==  bim[map_index[count],5],1] =  bim[map_index[count],6]
	temp[temp[,1] == "x",1] =  bim[map_index[count],5]

	temp[temp[,2] ==  bim[map_index[count],6],2] = "x"
	temp[temp[,2] ==  bim[map_index[count],5],2] =  bim[map_index[count],6]
	temp[temp[,2] == "x",2] =  bim[map_index[count],5]

	ped[,i:(i+1)] = temp
	count = count + 1
}


###################################	COnvert from Beagle imputed formate to phasebook format	######################
######################################################################################################################

Beagle = read.table(gzfile("seq.vcf.gz"),stringsAsFactor=F)

temp = cbind(unlist(lapply(Beagle[,4],nchar)),unlist(lapply(Beagle[,5],nchar)))
nrow(temp[temp[,1] > 1 | temp[,2] > 1,])
Beagle = Beagle[-which(temp[,1] > 1 | temp[,2] > 1),]

r2 = unlist(lapply(1:nrow(Beagle),function(x) { strsplit(Beagle[x,8],split=";")[[1]][2] } ))
r2_value = unlist(lapply(1:nrow(Beagle),function(x) { strsplit(r2[x],split="=")[[1]][2] } ))

Beagle = Beagle[which(as.numeric(r2_value) > 0.3),]



phasebook = data.frame(matrix(0,nrow=146*2,ncol=nrow(Beagle)+2))
count = 1
exclude = NULL
ex_index = 1
for(i in 10:ncol(Beagle))
{
	print(count)
	
	for(j in 1:nrow(Beagle))
	{
		print(j)
		if(strsplit(Beagle[j,i],split=":")[[1]][1] == "0|0") { 
			phasebook[count,(j+2)] = "1"
			phasebook[count+1,(j+2)] = "1" 
		} else if(strsplit(Beagle[j,i],split=":")[[1]][1] == "0|1") {
			phasebook[count,(j+2)] = "1"
			phasebook[count+1,(j+2)] = "2" 
		} else if(strsplit(Beagle[j,i],split=":")[[1]][1] == "1|0") {
			phasebook[count,(j+2)] = "2"
			phasebook[count+1,(j+2)] = "1" 
		} else if(strsplit(Beagle[j,i],split=":")[[1]][1] == "1|1") {
			phasebook[count,(j+2)] = "2"
			phasebook[count+1,(j+2)] = "2" 
		} else {
			print("Error")
			#print(strsplit(Beagle[j,i],split=":"))
			#exclude[ex_index] = Beagle[j,2]
			#ex_index = ex_index + 1
		}
		
	}
	count = count + 2
}

sample_names = 1:146
phasebook[,1] = unlist(lapply(1:146, function(i) { rep(sample_names[i],2) } ))
spaces = c(rep('    ',9*2),rep('   ',90*2),rep('  ',(146 - 99)*2))
x=c("1","2")
phasebook = cbind(spaces,phasebook[,1],rep(x,nrow(phasebook)/2),phasebook[,3:ncol(phasebook)])
write.table(phasebook,file=paste('phased_CHR_11_peak.txt',sep=''),quote=F,row.names=F,col.names=F)


r2 = unlist(lapply(1:nrow(Beagle),function(x) { strsplit(Beagle[x,8],split=";")[[1]][2] } ))
r2_value = unlist(lapply(1:nrow(Beagle),function(x) { strsplit(r2[x],split="=")[[1]][2] } ))
pdf("imputation_seq_quality.pdf")
hist(as.numeric(r2_value),main="imputation Quality")
dev.off()



########################################################################

Beagle = read.table(gzfile("seq.vcf.gz"),stringsAsFactor=F)

temp = cbind(unlist(lapply(Beagle[,4],nchar)),unlist(lapply(Beagle[,5],nchar)))
nrow(temp[temp[,1] > 1 | temp[,2] > 1,])
Beagle = Beagle[-which(temp[,1] > 1 | temp[,2] > 1),]

r2 = unlist(lapply(1:nrow(Beagle),function(x) { strsplit(Beagle[x,8],split=";")[[1]][2] } ))
r2_value = unlist(lapply(1:nrow(Beagle),function(x) { strsplit(r2[x],split="=")[[1]][2] } ))

Beagle = Beagle[which(as.numeric(r2_value) > 0.3),]

phasebook = data.frame(matrix(0,nrow=146,ncol=(nrow(Beagle))))
count = 1
for(i in 10:ncol(Beagle))
{
	print(count)
	x = 1
	phasebook[count,] = ifelse(substring(Beagle[,i],0,3) == "0|0","1 1",ifelse(substring(Beagle[,i],0,3) == "1|1","2 2",ifelse(substring(Beagle[,i],0,3) %in% c("0|1","1|0"),"1 2","error")))
	count = count + 1
}
Extreme = read.table("/root/Cow/PLINK_070415_0630/extremes_pheno.txt",stringsAsFactors=F)
sample_names = Extreme[,2]
Family_names = Extreme[,1]        
phasebook = cbind(Family_names,sample_names,rep(0,146),rep(0,146),rep(-9,146),rep(-9,146),phasebook)
write.table(phasebook,file=paste('seq_data.ped',sep=''),quote=F,row.names=F,col.names=F)



Map = cbind(rep(11,nrow(Beagle)),paste("X",1:nrow(Beagle),sep=""),rep(0,nrow(Beagle)),Beagle[,2])
write.table(Map,file=paste('seq_data.map',sep=''),quote=F,row.names=F,col.names=F)



