res = read.table("seq_Scores.assoc.linear",header=TRUE)
res = res[res$TEST == "ADD",]
res = na.omit(res)

res2 = read.table("More_Analysis/Second_peak_added.assoc.linear",header=TRUE)
res2 = res2[res2$TEST == "ADD",]
res2 = na.omit(res2)



png("More_Analysis/mht_plot_plink_effect_of_second_peak.png",width=290,height=210,units="mm",res=150)
plot(res$BP/10^6,-log10(res$P),xlab="Position in MB",ylab="-log10(pvalue)",pch=19,main="CHR 11 peak before and after accouting for second peak")
points(res2$BP/10^6,-log10(res2$P),xlab="Position in MB",ylab="-log10(pvalue)",pch=20,col="red")
dev.off()


#################### Boot Straping #############################
setwd("/media/Mahmoud_HD/Projects/Cow/ubuntu folder/Glascow_summary/peak/")
pheno = read.table("Alternative_pheno_v2.txt",header=TRUE)
ped = read.table("More_Analysis/Bootstrapping/seq_data.ped") #,colClasses=rep("numeric",27744))
Bootstrap_ped_list=NULL
Bootstrap_pheno_list=NULL
for(boot in 1:100)
{
        print(boot)        
        indv = sort(sample(x=1:146, size=146, replace = TRUE))
        Bootstrap_ped = ped[indv,]
        Bootstrap_ped[,2] = as.character(Bootstrap_ped[,2])
        
        duplicates = summary(as.factor(ped[indv,2]))
        Bootstrap_pheno = NULL
        Bootstrap_pheno = pheno[pheno$IID %in% names(duplicates),]
        rownames(Bootstrap_pheno) = 1:nrow(Bootstrap_pheno)
        duplicates = duplicates[duplicates > 1]
        duplicates = duplicates - 1
        
        
        for(i in 1:length(duplicates))
        {
                for(j in 1:duplicates[i])
                {
                        temp = Bootstrap_pheno[Bootstrap_pheno$IID == names(duplicates)[i],]
                        temp[,2] = paste(temp[,2],"_",j,sep="")
                        Bootstrap_pheno[(nrow(Bootstrap_pheno)+1),] = temp                                
                }
        }
        
        for(i in 1:length(duplicates))
        {
                found = grep(pattern=names(duplicates)[i],x=Bootstrap_ped[,2])
                count = 1
                for(j in 2:length(found))
                {
                        Bootstrap_ped[found[j],2] = paste(as.character(names(duplicates)[i]),"_",count,sep="")
                        count = count + 1
                }
        }
        Bootstrap_ped_list[[boot]] = Bootstrap_ped
        Bootstrap_pheno_list[[boot]] = Bootstrap_pheno
        
}
setwd("/media/Mahmoud_HD/Projects/Cow/ubuntu folder/Glascow_summary/peak/More_Analysis/Bootstrapping/")

for(i in 1:length(Bootstrap_ped_list))
{
                print(i)
        write.table(Bootstrap_ped_list[[i]],gzfile(paste("seq_data_boot_",i,".ped",sep="")),row.names=F,col.names=F,quote=F)
        write.table(Bootstrap_pheno_list[[i]],paste("Alternative_pheno_boot_",i,".txt",sep=""),quote=F,row.names=F,col.names=F)     
}


pca=read.table("PCA_file_farms.txt",header=TRUE)
for(i in 1:100)
{
        print(i)
        temp = NULL
        common = NULL
        diff = NULL
        common = intersect(Bootstrap_pheno_list[[i]][,2],pca[,2])
        temp = pca[pca[,2] %in% common,]
        rownames(temp) = 1:nrow(temp)
        diff = setdiff(Bootstrap_pheno_list[[i]][,2],pca[,2])
        for(j in 1:length(diff))
        {
                temp_diff = NULL
                temp_diff = pca[pca[,2] ==  strsplit(diff[j],split="_")[[1]][1], ]
                temp_diff[,2] = diff[j]
                rownames(temp_diff) = (nrow(temp)+1)
                temp[(nrow(temp)+1),] = temp_diff
        }
        write.table(temp,paste("PCA_file_farm_boot_",i,".txt",sep=""),row.names=F,col.names=F,quote=F)
}

