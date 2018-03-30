justin=read.table("/root/Dropbox/IBD Morocan/PLOTS/justin_regions.txt")

DE = c("CD","UC","IBD")
for(j in 1:3)

res = read.table(paste(DE[j],".assoc.logistic",sep=""),header=T)
res = res[res$TEST == "ADD",]
chr_no = data.frame(chr=1:23,No=summary(as.factor(res[,1])))
colors = NULL
ticks = NULL
counter = 1
chr_names = unlist(lapply(1:23,function(i) { rep(i,chr_no[i,2]) }))
space = 2000
for(i in 1:nrow(chr_no))
{
	if(i %in% c(1,3,5,7,9,11,13,15,17,19,21,23))
	{
		colors = c(colors,rep("grey",chr_no[i,2]),rep("white",space))
		if(i == 1)
		{
			ticks[i] = chr_no[i,2]/2 
		}
		else
		{
			ticks[i] = sum(chr_no[1:(i-1),2]) + chr_no[i,2]/2 + (space*(i-1))
		}
	}
	else
	{
		colors = c(colors,rep("black",chr_no[i,2]),rep("white",space))
		ticks[i] = sum(chr_no[1:(i-1),2]) + chr_no[i,2]/2 + (space*(i-1))
	}
}
temp = cbind(res,chr_names)
insert = matrix(10^-7,ncol=10,nrow=space)
colnames(insert) = colnames(temp)
new = NULL
for(i in 1:nrow(chr_no))
{
	chr = temp[temp$chr_names == i,]
	new =  rbind(new,chr,insert)
	
}

new= cbind(new,colors)
new$colors = as.character(new$colors)
for(i in 1:nrow(justin))
{
	found =	new[new$CHR == justin$CHR[i] & new$BP >= (justin$start[i]*1000000) & new$BP <= (justin$end[i]*1000000),]
	if(nrow(found) > 0)
	{
		new[new$CHR == justin$CHR[i] & new$BP >= (justin$start[i]*1000000) & new$BP <= (justin$end[i]*1000000),]$colors = "red"
	}
}
new$pch = 19
new[new$colors == "red",]$pch = 17

#new[,1] = 1:nrow(new)
png(paste("/root/Dropbox/IBD Morocan/PLOTS/",DE[j],".assoc.logistic.with_risk_regions.png",sep=""),width=390,height=210,units="mm",res=450)
plot(1:nrow(new),-log10(new[,9]),col=new$colors,xaxt = "n",yaxt="n",pch=new$pch,xlab="",ylab="",ylim=c(0,8))
axis(side=2, at=1:7,labels=1:7,cex.axis=1.3)
axis(side=1, at=ticks,labels=FALSE)
text(ticks, par("usr")[3] - 0.2, labels = 1:23, pos = 1, xpd = TRUE,cex=c(rep(1.3,17),0.9,1.3,0.9,1.3,0.9,1.3))
title(ylab=expression(paste('-log' [10],'(p)',sep="")),cex.lab=1.7,line=2)
title(xlab="Chromosome",cex.lab=1.7)
abline(h=-log10(0.05/nrow(na.omit(res))),col=2,lty=2,lwd =2,cex.lab=2)		# Red line genomewide Threshold
legend("topright",c("Justin Risk Regions","Non-Risk Regions"),col=c("red","black"),pch=c(17,19))
dev.off()

}


