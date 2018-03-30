#################################################################################
# Plot Heterozygosity rate vs Missing rate
#################################################################################
het = read.table('/home/mahmoud/plink.het',header=T,sep="\t")
missing = read.table('/home/mahmoud/plink.imiss',header=T,sep="\t")
jpeg('/home/mahmoud/hetvsmissing.jpeg',quality=100)
plot((het$N.NM. - het$O.HOM.)/ het$N.NM. , missing$F_MISS ,xlab="heterozygosity rate",ylab="missing rate",main="missing vs het")
abline(h=0.01)
abline(lm(missing$F_MISS ~ (het$N.NM. - het$O.HOM.)/ het$N.NM. ))
dev.off()


########################################################################################
# Plot Heterozygosity rate vs Missing rate after removing the outlying animals(?) OR by renaming
#################################################################################
temp = cbind( het=(het$N.NM. - het$O.HOM.)/ het$N.NM. , missing=missing$F_MISS)
temp = temp[temp[,2] <= 0.01,]
temp = as.data.frame(temp)
jpeg('/home/mahmoud/hetvsmissing2.jpeg',quality=100)
plot(temp$het , temp$missing ,xlab="heterozygosity rate",ylab="missing rate",main="missing vs het")
abline(lm( missing ~ het,data=temp))
dev.off()


##########################################################################################
#recode sex in ped file as given in the sex_file
##########################################################################################
ped_file = read.table('/home/mahmoud/Analysis.ped',sep="")
map_file = read.table('/home/mahmoud/Analysis.map',sep="")
sex_file = read.table('/home/mahmoud/sex_status.txt',sep="\t",header=T)
index = NULL
for(i in 1:nrow(ped_file)){
	index = NULL
	#grep searches for matches with the argument behind "pattern': search the animals in the sex_file en match ze
	index = grep(pattern = paste("^",as.character(ped_file[i,2]),"$",sep="") , x= as.character(sex_file$IID))
	if(length(index) > 0){
		#recode the sexes from M/F/missing to 1/2/-9
		ped_file[i,5] =   ifelse(as.character(sex_file$sex[index]) == "m",1,ifelse(as.character(sex_file$sex[index]) == "f",2,-9))
	}
	else{
		print("Error")
	}
}

##########################################################################################
#calculate heterozygosity in PAB/ on X chrom???
##########################################################################################
b = (52602 * 2) + 6
e = (53489 * 2) + 6
count = 0
temp = NULL
for(j in 1:nrow(ped_file))
{
	count = 0
	for(i in seq(start=b,end=e,by=2))
	{
		if(as.character(ped_file[j,i]) != as.character(ped_file[j,i+1]))
		{
			count = count + 1
		}
	}
	temp[j] = count 
}

##############################################################################################
## 	Plot AIREML for one marker haplotype
##############################################################################################
state = read.table("solutions",header=TRUE,sep="",na.strings="",stringsAsFactors=FALSE,check.names=FALSE, row.names=NULL)
state = state[,-1]
state = state[-1,]
pheno = read.table("pheno.txt")
x = barplot(summary(as.factor(c(pheno[,3],pheno[,4]))))
png("Test_Top_Marker_haplotype_cluster.png",width=290,height=210,units="mm",res=150)
barplot(summary(as.factor(c(pheno[,3],pheno[,4]))),xlim=c(0.5,25),ylim=c(0,50),col="skyblue",xlab="Haplotype Cluster",ylab="Number of chromosomes")
par(new = TRUE)
y=state[,3]
sd=state[,4]
plot(x,y,xlab="",ylab="",ylim=c(-1,1),xlim=c(0.5,25),axes=F,col="red",pch=19)
arrows(x, y-sd, x, y+sd, length=0.05, angle=90, code=3)
mtext("effect on disease score",side=4,col="red",line=-0.5)
axis(4, ylim=c(-1,1), col="red",col.axis="red",line=-2.5)
dev.off()

Top_Marker = states_11[,c(1,1742)
control = Top_Marker[Top_Marker[,1] %in% pheno[pheno[,3] == 0,1],]
cases = Top_Marker[Top_Marker[,1] %in% pheno[pheno[,3] == 1,1],]
height = rbind(  unlist(lapply(paste("^",1:20,"$",sep=""),function(x) { length(grep(pattern=x,x=control[,2])) } )) ,   unlist(lapply(paste("^",1:20,"$",sep=""),function(x) { length(grep(pattern=x,x=cases[,2])) } )) )
colnames(height) = 1:20
png("Top_Marker_haplotype_cluster_case_control_frq.png",width=290,height=210,units="mm",res=150)                
mp = barplot(height,ylim=c(0,40),xlab="Haplotype Cluster",ylab="Number of chromosomes",main="Top Marker 'ARS-BFGL-NGS-118288' frequency in cases and controls",col=c("blue","red"),beside=TRUE)
text(mp, height, labels = format(height, 4),pos = 3, cex = .75)
legend("topright", fill=c("blue","red"), legend=c("Controls", "Cases"));
dev.off()

for(hap in 1:20)
{
        print(paste("cluster",hap))
        temp = NULL                
        for(i in 1:146)
        {
                temp[i] = sum(Top_Marker[Top_Marker[,1] == i,2] == hap)
        }
        temp = cbind(1:146,temp,scores[,3],pheno[,3])
        #col = rep("black",nrow(temp))
        #col[which(temp[,3] == 0)] = "red"
        temp = cbind(temp,rep("red",nrow(temp)))
        temp[which(temp[,4] == 0),5] = "black" 
        png(paste("frq_haplotype_cluster_",hap,".png",sep=""),width=290,height=210,units="mm",res=150)
        plot(1:146,temp[order(as.numeric(temp[,3])),2],col=temp[order(as.numeric(temp[,3])),5],pch=19,main=paste("Haplotype Cluster ",hap,sep=""),xlab="Individuals sorted by their score",ylab="Number of times the haplotype appeared",ylim=c(0,2))
        dev.off()
        print(sum(as.numeric(temp[temp[,2] > 0 & temp[,5] == "red",2]))) # Cases
        print(sum(as.numeric(temp[temp[,2] > 0 & temp[,5] == "black",2]))) # Control
}       

for(hap in 1:20)
{
       print(paste("cluster",hap))
       temp = NULL                
       for(i in 1:146)
       {
               temp[i] = sum(Top_Marker[Top_Marker[,1] == i,2] == hap)
       }
       temp = cbind(1:146,temp,scores[,3],pheno[,3])
       #col = rep("black",nrow(temp))
       #col[which(temp[,3] == 0)] = "red"
       temp = cbind(temp,rep("red",nrow(temp)))
       temp[which(temp[,4] == 0),5] = "black" 
       png(paste("frq_real_scores_haplotype_cluster_",hap,".png",sep=""),width=290,height=210,units="mm",res=150)
       plot(temp[order(as.numeric(temp[,3])),3],temp[order(as.numeric(temp[,3])),2],col=temp[order(as.numeric(temp[,3])),5],pch=19,main=paste("Haplotype Cluster ",hap,sep=""),xlab="Individuals sorted by their score",ylab="Number of times the haplotype appeared",ylim=c(0,2))
       dev.off()
       print(sum(as.numeric(temp[temp[,2] > 0 & temp[,5] == "red",2]))) # Cases
       print(sum(as.numeric(temp[temp[,2] > 0 & temp[,5] == "black",2]))) # Control
}      
                       
setwd("/media/Mahmoud_HD/Projects/Cow/ubuntu folder/Glascow_summary/peak/")
pca=read.table("PCA_file_farms.txt",header=TRUE)     
rownames(pca)=as.character(pca$IID)
pca = pca[rownames(mat),]
library("SNPRelate")
snpgdsBED2GDS("seq_data.bed","seq_data.fam","seq_data.bim","Data.gds")                   
genofile <- openfn.gds("Data.gds")
snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))
mat <- snpgdsGetGeno(genofile)
                       
####### LD with the highest peak (SNP) ############

LD = unlist(lapply(1:ncol(mat),function(x) { snpgdsLDpair(mat[,grep(pattern="X11834",x=snpset)],mat[,x],method="r")$ld^2 } ))                       
col = rep("black",nrow(res))
col[which(LD == 1)] = "black"
col[which(LD >= 0.8 & LD < 1)] = "blue"
col[which(LD >= 0.5 & LD < 0.8)] = "yellow"
col[which(LD >= 0.20 & LD < 0.5)] = "orange"
pch = rep(23,nrow(res))
pch[which(LD == 1)] = 18
pch[which(LD >= 0.2)] = 18
cex = rep(1.25,nrow(res))
cex[which(LD == 1)] = 4.5
                      
                       
names = read.table("seq_data.fam")                       
rownames(mat) = names[,2]
genotypes_labels = c("AA","AB","BB")
scores = read.table("Alternative_pheno_v2.txt",header=TRUE)
rownames(scores) = as.character(scores[,2])
SNP = data.frame(mat[,grep(pattern="X6319",x=snpset)],as.character(scores[rownames(mat),3]),pca[rownames(mat),3:7],as.character(scores[rownames(mat),8]),check.names=F)
SNP = SNP[SNP[,1] < 3,]
SNP = na.omit(SNP)
png(paste("Top_SNP_1st_peak_boxplot.png",sep=""),width=290,height=210,units="mm",res=150)
plot(as.factor(SNP[,1]),as.numeric(as.character(SNP[,2])),xlab=paste("SNP ","rs385007695",sep=""),ylab="scores phenotype",xaxt='n',main="Top SNP in the first peak")
axis(side=1, labels=c("AA","AB","BB"),at=c(1,2,3),cex.axis=1.3,srt=45)
points(as.factor(SNP[,1]),as.numeric(as.character(SNP[,2])),col=SNP[,8],pch=19)
legend("topright", col=c("black","red"), legend=c("Controls", "Cases"),pch=19);
dev.off()
model = lm(as.numeric(as.character(SNP[,2]))~SNP[,1] + SNP[,3] + SNP[,4] + SNP[,5] + SNP[,6] + SNP[,7])
abline(a=coefficients(model)[1],b=coefficients(model)[2])                       

png(paste("Top_SNP_1st_peak_boxplot_mht.png",sep=""),width=390,height=170,units="mm",res=150)                       
par(mfrow=c(1,2),oma = c(0, 0, 3, 0))
plot(as.factor(SNP[,1]),as.numeric(as.character(SNP[,2])),xlab=paste("SNP ","rs385007695",sep=""),ylab="scores phenotype",xaxt='n',main="Top SNP in the first peak")
axis(side=1, labels=c("AA","AB","BB"),at=c(1,2,3),cex.axis=1.3,srt=45)
points(as.factor(SNP[,1]),as.numeric(as.character(SNP[,2])),col=SNP[,8],pch=19)
legend("topright", col=c("black","red"), legend=c("Controls", "Cases"),pch=19);                       
res = read.table("seq_Scores.assoc.linear",header=TRUE)
res = res[res$TEST == "ADD",]
#res = na.omit(res)
col = rep("black",nrow(res))
col[grep(pattern="X6319",x=as.character(res$SNP))] = "red"
plot(res$BP/10^6,-log10(res$P),pch=19,xlab="position in MB",ylab="-log10(pvalue)",main="manhattan plot of chr11 imputed region using scores as phenpotype",col=col)                     
dev.off()
                       
                       
                       
                       
########################### MHT PLOT ##############################
setwd("/media/elansary/Mahmoud_HD/Projects/Cow/ubuntu folder/Glascow_summary/peak/")
library("SNPRelate")
snpgdsBED2GDS("seq_data.bed","seq_data.fam","seq_data.bim","Data.gds")                   
genofile <- openfn.gds("Data.gds")
snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))
mat <- snpgdsGetGeno(genofile)
pca=read.table("PCA_file_farms.txt",header=TRUE)     
rownames(pca)=as.character(pca$IID)
pca = pca[rownames(mat),]

LD = unlist(lapply(1:ncol(mat),function(x) { snpgdsLDpair(mat[,grep(pattern="X11834",x=snpset)],mat[,x],method="r")[1]^2 } ))  
LD = as.double(LD)
LD[which(is.na(LD))] = 0
col = rep("black",nrow(res))

col[grep(pattern="X11834",x=snpset)] = "black"
col[which(LD >= 0.8 & LD < 0.99)] = "blue"
col[which(LD >= 0.5 & LD < 0.8)] = "yellow"
col[which(LD >= 0.20 & LD < 0.5)] = "orange"
pch = rep(23,nrow(res))
pch[grep(pattern="X11834",x=snpset)] = 18
pch[which(LD >= 0.2)] = 18
cex = rep(1.25,nrow(res))
cex[grep(pattern="X11834",x=snpset)] = 4.5

res = read.table("seq_Scores.assoc.linear",header=TRUE)
res = res[res$TEST == "ADD",]
map = read.table("/media/elansary/Mahmoud_HD/Projects/Cow/ubuntu folder/PLINK_070415_0630/PSOROVIS_clean_final.map")
recomb = map[map[,1] == 11 & map[,4] >= min(res$BP) & map[,4] <= max(res$BP),]


#par(mfrow=c(2,1),oma = c(0, 0.5, 0, 0))
png("mht_plot_peak_gene_anotation.png",width=330,height=210,units="mm",res=150)                       
par(mar=c(5,4,4,5)+.1)
plot(res$BP/10^6,-log10(res$P),yaxt="n",cex.lab=1.2,pch=23,xlab="position in MB",ylab=expression('-log'[10]*'(pvalue)'),main="manhattan plot of chr11 imputed region using scores as phenpotype",col="black",cex=1.28,ylim=c(-6,8))                     
axis(2, at=seq(0,8,2), labels=seq(0,8,2),tck = -0.02)
points(res$BP/10^6,-log10(res$P),pch=pch,col=col,cex=cex)
legend("topleft", fill=c("red","blue","yellow","orange","white"), legend=c("r2 = 1", "0.8 ≤ r2 <1","0.5 ≤ r2 < 0.8","0.2 ≤ r2 < 0.5","r2 < 0.2"));
text((res[which.min(res$P),]$BP/10^6) + 0.25 , 7.5, labels = "rs110723767",pos = 3, cex = 1.3)
text((res[which.min(res$P),]$BP/10^6) + 0.25 , 7, labels = bquote("P = 9.5 x 10"^"-8"*"" ),pos = 3, cex = 1.3)

par(new=TRUE)
plot(recomb[,4], recomb[,3]-107.350,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-2,4))
points(recomb[,4], recomb[,3]-107.350,col="red",pch="x")
axis(4)
mtext("recommbination rate (cM/MB)",side=4,line=3)

geneNames_2 = read.table("/media/Mahmoud_HD/Projects/Cow/Gene Names/mart_export.txt",header=TRUE,stringsAsFactors=FALSE,sep=",")
geneNames_2 = geneNames_2[!duplicated(geneNames_2$Ensembl.Gene.ID),]                 
geneNames_2 = geneNames_2[(geneNames_2$Gene.Start..bp./10^6 >= 106.5 & geneNames_2$Gene.Start..bp./10^6 <= 107)  | (geneNames_2$Gene.Start..bp./10^6 >= 106 & geneNames_2$Gene.Start..bp./10^6 <= 106.3),]
                       
                       
geneNames = read.table("/media/Mahmoud_HD/Projects/Cow/Gene Names/geneNames",header=TRUE,stringsAsFactors=FALSE)
geneNames = geneNames[!(geneNames$txStart/10^6) < 104.5,]
                       
                       
#plot(res$BP/10^6,rep(60,nrow(res)),col="white",axes=FALSE, frame.plot=FALSE,ann=FALSE)
count = seq(-0.6,-6,1*((-6+0.6))/nrow(geneNames_2))
for(i in 1:nrow(geneNames_2))
{
       if(geneNames_2$Strand[i] == 1)
       {
               arrows(x0=geneNames_2$Gene.Start..bp.[i]/10^6,x1=geneNames_2$Gene.End..bp.[i]/10^6,y0=count[i],y1=count[i], angle=15, length=.12, lwd=5,col="red")                       
               text( (geneNames_2$Gene.Start..bp.[i]/10^6) + (geneNames_2$Gene.End..bp.[i] - geneNames_2$Gene.Start..bp.[i])/10^6, count[i], labels = geneNames_2$Associated.Gene.Name[i],pos = 3, cex = .75)
       } else {
               arrows(x0=geneNames_2$Gene.End..bp.[i]/10^6,x1=geneNames_2$Gene.Start..bp.[i]/10^6,y0=count[i],y1=count[i], angle=15, length=.12, lwd=5,col="blue")                       
               text( (geneNames_2$Gene.Start..bp.[i]/10^6) + (geneNames_2$Gene.End..bp.[i] - geneNames_2$Gene.Start..bp.[i])/10^6, count[i], labels = geneNames_2$Associated.Gene.Name[i],pos = 3, cex = .75)
       }      
}
dev.off()
                       
count = seq(-0.3,-4,-0.078)
for(i in 1:nrow(geneNames))
{
        if(geneNames$strand[i] == "+")
        {
                arrows(x0=geneNames$txStart[i]/10^6,x1=geneNames$txEnd[i]/10^6,y0=count[i],y1=count[i], angle=15, length=.12, lwd=5,col="red")                       
        } else {
                arrows(x0=geneNames$txEnd[i]/10^6,x1=geneNames$txStart[i]/10^6,y0=count[i],y1=count[i], angle=15, length=.12, lwd=5,col="red")                       
                
        }      
}
                       
setwd("/media/Mahmoud_HD/Projects/Cow/ubuntu folder/Glascow_summary/extras/Bootstrapping/")
res = read.table("Number_of_times_in_top_10.txt",header=TRUE)
dbSNP = read.table("/media/Mahmoud_HD/Projects/Cow/Gene Names/rs_names_peak.txt",header=TRUE,sep=",")
as.character(unlist(lapply(paste("^",res$BP,"$",sep=""),function(x){ dbSNP[grep(pattern=x,x=as.character(dbSNP[,4])),1]})))
new_res = cbind(res[,1:2],as.character(unlist(lapply(paste("^",res$BP,"$",sep=""),function(x){ dbSNP[grep(pattern=x,x=as.character(dbSNP[,4])),1]})))
,res[3:ncol(res)])
head(new_res)
colnames(new_res)[3] = "dbSNP"
write.table(new_res,"Number_of_times_in_top_10.txt",quote=F,col.names=TRUE,row.names=F)
                       
                       
                       
################################################
                    
emission = read.table("/media/Mahmoud_HD/Projects/Cow/ubuntu folder/Glascow_summary/extras/clusters/emission_prob.txt")
markers = read.table("PSOROVIS_clean_final_12.bim")
alleles = markers[markers[,1] == 11 & markers[,4]/10^6 > 106.43 & markers[,4]/10^6 <= 106.75,5:6]
map = read.table("map_CHR_11.txt")      
states = lapply(1:20,function(x) { emission[emission[,1] == x,c(1,map[which(as.character(map[,2]) %in% as.character(markers[markers[,1] == 11 & markers[,4]/10^6 > 106.43 & markers[,4]/10^6 <= 106.75,2])),1]+1)] }  )


for(j in 1:20)
{
        print(j)
        m=matrix(0,nrow=4,ncol=9)
        rownames(m) = c("A","C","G","T")
        colnames(m) = 1:9
        count = 1
        for(i in 2:10)
        {
                m[alleles[count,1],count]  = states[[j]][1,i]        
                m[alleles[count,2],count] =  states[[j]][2,i]        
                count = count + 1
                
        }
        png(paste("motif_hidden_state_",j,".png",sep=""),width=290,height=210,units="mm",res=150)
        seqLogo(m,ic.scale=FALSE)
        dev.off()
                               
}
setwd("/media/Mahmoud_HD/Projects/Cow/ubuntu folder/Glascow_summary/extras/clusters");
states = read.table("states_11");
phases = read.table("phases_11");
temp = states[,c(1,2,1744)];
Mean = NULL;
freq = NULL;
for(i in 1:20)
{
temp_phases = phases[which(temp[,3] == i),]
Mean[[i]] = unlist(lapply(temp_phases[,3:ncol(temp_phases)],mean))
freq[[i]] =  unlist(lapply(temp_phases[,3:ncol(temp_phases)],function(x) { summary(as.factor(x))["1"]/nrow(temp_phases) } ))
}
#plot(1:101,Mean[[14]][1653:1753],pch=1,ylim=c(1,2),type="n")
#plot(1:21,freq[[14]][1732:1752],pch=1,ylim=c(1,2),type="n")            
for(i in 1:20)
{
        dir.create(paste("State_",i,sep=""))
        print(i)
        for(j in 1:20)
        {
                
                png(paste("State_",i,"/plot_",i,"_vs_",j,".png",sep=""),width=330,height=210,units="mm",res=150)
                plot(1:21,freq[[i]][1732:1752],pch=1,ylim=c(1,2),type="n",ylab="Mean of the Marker",xlab="Subset of Markers",main=paste("Hidden state ",i," vs Hidden state ",j,sep=""))            
                
                freq[[j]][which(is.na(freq[[j]]))] = 0 
                freq[[i]][which(is.na(freq[[i]]))] = 0 
                #text(1:101,Mean[[i]][1653:1753],label=i)
                #text(1:101,Mean[[i]][1732:1752],label=i)
                text(1:101,(1 - freq[[i]][1732:1752]) + 1,label=i)
                text(1:101,(1 - freq[[j]][1732:1752]) + 1,label=j)
        
                abline(v=11,col="red")                       
                dev.off()
        }
}                      
                       
