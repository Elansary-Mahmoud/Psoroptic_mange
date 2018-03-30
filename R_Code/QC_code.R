####################################################################################################
#################################### Script For Cow Analysis	#####################################
####################################################################################################
###	This scrit will be used to finish the Quality Control steps and Analysis step           ###
###										       ###
####################################################################################################

####################################################################################################
####################################################################################################
##########			The Quality Control Steps				##########
####################################################################################################
####################################################################################################
			######## Prepare the directories for plink ###########
####################################################################################################
plink     = "/home/ulg/genan/mahmoud/plink-1.07-x86_64/plink";
Geno_dir  = "/home/ulg/genan/mahmoud/Annalys/Final_Analysis/";
pheno_dir = "/home/ulg/genan/mahmoud/Annalys/Final_Analysis/";
out_dir   = "/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs/";

Geno_file  = "raw-GWA-data";
pheno_file = "";
out_file   = "raw-GWA-data";

####################################################################################################

	      ######## Identify all markers with an excessive missing data rate ###########
system(paste(plink," --noweb --cow --bfile ",Geno_dir,Geno_file," --geno 0.05 --recode --make-bed --out ",out_dir,"QC-GWA-data-SNP",sep=""))

	  	       ######## Hardy Weinberg equilibrium Test ######## 
	  
system(paste(plink," --noweb --cow --bfile ",out_dir,"QC-GWA-data-SNP"," --hardy --adjust --out ",out_dir,"QC-GWA-data-SNP",sep=""))


HWE = read.table("/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs/QC-GWA-data-SNP.hwe",header=T,stringsAsFactors=F)
HWE = HWE[HWE$TEST == "ALL",]#in case/control: only controls are used because you expect a deviation in hwe in the diseased
HWE$qvalue = p.adjust(HWE$P,method="fdr")#correct for multiple testing , method =  False discovery rate, calculate a p/q-value)
thr = HWE[HWE$qvalue < 0.05 ,]#select all markers with a qvalue < 0.05)

thr = thr[which.max(thr$qvalue),]$P # find the largest of these selected q values, to use it as your thershold value


system(paste(plink," --noweb --cow --bfile ",out_dir,"QC-GWA-data-SNP"," --maf 0.05 --hwe ",thr," --recode --make-bed --out ",out_dir,"QC-GWA-data-SNP-MAF-HWE",sep=""))

 #####################################################################################################################
 #############################################		QC Animals		######################################
 #####################################################################################################################
 

##########################################missingness vs heterozygosity  ###################################
system(paste(plink," --noweb --cow --bfile ",out_dir,"QC-GWA-data-SNP-MAF-HWE"," --missing --out ",out_dir,"QC-GWA-data-SNP-MAF-HWE",sep=""))
system(paste(plink," --noweb --cow --bfile ",out_dir,"QC-GWA-data-SNP-MAF-HWE"," --het --out ",out_dir,"QC-GWA-data-SNP-MAF-HWE",sep=""))

imiss=read.table(paste(out_dir,"QC-GWA-data-SNP-MAF-HWE.imiss",sep=""),h=T)
imiss$logF_MISS = log10(imiss[,6])#F_miss = proportion of missing SNPs, log = for visualisation
het=read.table(paste(out_dir,"QC-GWA-data-SNP-MAF-HWE.het",sep=""),h=T)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.#calculate heterozygosity BY SNP from observed homozygosity and total non-missing values.

library("geneplotter")
colors  <- densCols(imiss$logF_MISS,het$meanHet)
pdf(paste(out_dir,"QC-GWA-data-SNP-MAF-HWE.imiss-vs-het.pdf",sep=""))
plot(imiss$logF_MISS,het$meanHet, col=colors, xlim=c(min(imiss$logF_MISS),0),ylim=c(0,0.5), xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
axis(1,at=c(min(imiss$logF_MISS),-4,-3,-2,-1,0),labels=c(10^min(imiss$logF_MISS),0.0001,0.001,0.01,0.1,1))
abline(h=mean(het$meanHet)-(3*sd(het$meanHet)),col="RED",lty=2)
abline(h=mean(het$meanHet)+(3*sd(het$meanHet)),col="RED",lty=2)
abline(v= -1.30103, col="RED", lty=2)#threshold missingnes: log10 0.05 = -1.3
dev.off()# shuts down the current graphic device#

fail_imisshet_qc = imiss[imiss$F_MISS > 0.05,]  # the names of animals with low call rate less than 95%
write.table(fail_imisshet_qc[,1:2],file=paste(out_dir,"fail-imisshet-qc.txt",sep=""),row.name=F,col.name=F,quote=F)

system(paste(plink," --noweb --cow --bfile ",out_dir,"QC-GWA-data-SNP-MAF-HWE --remove ", out_dir,"fail-imisshet-qc.txt --recode --make-bed --out ",out_dir,"clean-inds-GWA-data",sep=""))

######################################## cluster the data and draw MDS plot#########################

system(paste(plink," --noweb --cow --bfile ",out_dir,"clean-inds-GWA-data --genome --out ",out_dir,"clean-inds-GWA-data",sep=""))
 
system(paste(plink," --noweb --cow --bfile ",out_dir,"clean-inds-GWA-data --read-genome ",out_dir,"clean-inds-GWA-data.genome --cluster --mds-plot 4 --out ",out_dir,"clean-inds-GWA-data",sep=""))

MDS = read.table(paste(out_dir,"clean-inds-GWA-data.mds",sep=""),header=T,stringsAsFactor=F)
pdf(paste(out_dir,"MDS-plot.pdf",sep=""))
plot(MDS$C1,MDS$C2,xlab="PC 1",ylab="PC 2",main="MDS Plot",pch=19)
dev.off()

########################################## Sex Check  ###############################################

map = read.table(paste(out_dir,'clean-inds-GWA-data.bim',sep=""))
chrX = map[map[,1] == 30,] # select the snp's on the X-chromosome
chrX = chrX[chrX[,4] < 137000000,]  # 137 million base pair   or 143451431 for GPPR143 = leave out the pseudoautosomal region

library('SNPRelate')
setwd('/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs/')
snpgdsBED2GDS('clean-inds-GWA-data.bed', 'clean-inds-GWA-data.fam','clean-inds-GWA-data.bim', "Data.gds")
genofile <- openfn.gds("Data.gds")
snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))
Sex = read.table("clean-inds-GWA-data.fam")#fam = 6 columns of PED = FID,IID, SireID, DamID, Sex, Pheno)
Sex_new = read.table("/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Pheno.txt",header=T)
rownames(Sex_new) = Sex_new[,2] # name both equal: necessary for comparison
Sex_new = Sex_new[as.character(Sex[,2]),]#what does this do exactly?

mat <- snpgdsGetGeno(genofile)#function not described in SNPrelate manual? wel on internet: To get a genotype matrix from a specified GDS file
chrX_geno = mat[,as.numeric(rownames(chrX))]# takes the columns from the matrix that contain snp's on the X chrom?
rownames(mat) = Sex[,2]# rownames = IID

discordant_sex = NULL
for(i in 1:nrow(mat))#nrow mat = # individuals
{
	discordant_sex[i] = sum(chrX_geno[i,] == 1) / (ncol(chrX_geno) - sum(chrX_geno[i,] == 3))#??????? when ==3 == missing, 1=heterozygous, 0, 2= homozygous????
}

Sex$discordant_sex = discordant_sex

Sex_new = data.frame(Sex_new,discordant_sex)


pdf(paste(out_dir,"Sex_prediction.pdf",sep=""))
P1 = hist(Sex_new[Sex_new$Sex == 2,]$discordant_sex)
P2 = hist(Sex_new[Sex_new$Sex == 1,]$discordant_sex)
plot( P2, col=rgb(1,0,0,1/4), main="Sex prediction Histogram",ylim=c(0, max(P1$counts)),xlim=c(0,0.5)) 
plot( P1, col=rgb(0,0,1,1/4),add=T,main="Sex prediction Histogram")
legend("topright",col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),legend=c("Females","Males"), bty="n",pch = 19)
dev.off()


temp = rbind(Sex_new[Sex_new$Sex == 1 & Sex_new$discordant_sex > 0.1,1:2],Sex_new[Sex_new$Sex == 2 & Sex_new$discordant_sex < 0.1,1:2])#?
write.table(temp,file="",quote=F,col.name=F,row.name=F)


system(paste(plink," --noweb --cow --bfile ",out_dir,"clean-inds-GWA-data --remove ", out_dir,"Sex_Problem.txt --recode --make-bed --out ",out_dir,"clean-inds-GWA-data",sep=""))

############################################################################################
############################		PCA		####################################
############################################################################################


RV <- snpgdsPCA(genofile, num.thread=2)
plot(RV$eigenvect[,1], RV$eigenvect[,2], xlab="PC 1", ylab="PC 2",pch=19)
temp = data.frame(Sex[,1:2],RV$eigenvect[,1:5])
colnames(temp) = c("FID","IID",paste("PC",1:5,sep=""))
write.table(temp,file='/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Outputs/PCA_file.txt',quote=F,row.names=F)



Pheno = read.table("/home/ulg/genan/mahmoud/Annalys/Final_Analysis/Pheno.txt",header=T)
res = NULL
for(i in 1:nrow(Sex))
{
	if(length(grep(pattern=paste("^",Sex[i,2],"$",sep=""),x=as.character(Pheno[,2]))) == 1) {
		Sex[i,7:13] = Pheno[grep(pattern=paste("^",Sex[i,2],"$",sep=""),x=as.character(Pheno[,2])),]
		
	} else {
		print("Error")
	}
}


#Add PCA-factors as covariates, add the phenotypes to the plink files
system(paste(plink," --noweb --cow --bfile ",out_dir,"clean-inds-GWA-data --covar ", out_dir,"PCA_file.txt --linear --pheno ",out_dir,"Alternative_pheno.txt --all-pheno --out ",out_dir,"Analysis_Files",sep=""))

names = c("PhenoA","PhenoB","PhenoC","PhenoD")
for(i in 4:4)
{
	print(names[i])
	Results = read.table(paste(out_dir,"Analysis_Files.",names[i],".assoc.linear",sep=""),header=T,stringsAsFactor=F)
	Results = Results[Results$TEST == "ADD",]
	Results = Results[Results$CHR != 30,]
	Results = Results[Results$CHR != 31,]
	#mhtplot(Results[,c(1,3,9)],control = mht.control(usepos = FALSE,col=colors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
	png(paste(out_dir,names[i],"/qqplot.png",sep=""))
	qqunif(Results$P,main=names[i])
	dev.off()

	for(j in 1:29)
	{	
		print(paste("CHR",j,sep=""))
		png(paste(out_dir,names[i],"/CHR",j,".png",sep=""))
		plot(Results[Results$CHR == j,]$BP,-log10(Results[Results$CHR == j,]$P),main=paste(names[i]," CHR ",j,sep=""))	
		dev.off()
	}	
	
}

#############################################		Glascow		##############################################
# this is an old version, were we ran glascow on genotypes and not on haplotypes + were PCA was used as random factor and not the correlation matrix

system(paste(plink," --noweb --cow --bfile ",out_dir,"clean-inds-GWA-data --recode12 --out ",out_dir,"GFILE_GLASCOW",sep=""))
x = system("awk {'print NF'} GFILE_GLASCOW.ped | head -1",intern=TRUE)
system(paste("cut -d' ' -f 1,7-",x," GFILE_GLASCOW.ped > GFILE_GLASCOW.txt",sep=""))
PCA = read.table(paste(out_dir,"PCA_file.txt",sep=""),header=T)
Pheno = read.table(paste(out_dir,"Alternative_pheno.txt",sep=""),header=T)

PFILE = data.frame(ID=Pheno[,1],PCA[,3:7],PhenoA=Pheno[,3])
write.table(PFILE,file=paste(out_dir,"PFILE_GLASCOW.txt",sep=""),quote=F,row.names=F)



