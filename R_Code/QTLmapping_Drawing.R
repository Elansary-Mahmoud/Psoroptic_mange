library('gap')
setwd("C:\\Users\\Administrator\\Documents\\Visual Studio 2010\\Projects\\porc_Test\\Backups\\");

results <- read.table("TDT_2.txt", header=F)
results = results[,c(2,1,3,7,4,5,8,9,10)]
colnames(results) = c("Marker_ID","CHR","BP","P","T","U","Nonsense","Conflict","Missing")
mycolors <- c("red","blue","green","cyan","brown4","black","magenta","red","blue","green",
            "cyan","brown4","black","magenta","red","blue","green","cyan","black")
results$BP = results$BP/1000000

results = results[-log10(results$P) < 4,]
rownames(results) = seq(1:1,nrow(results))
T = T[-log10(T$P) > 5,]
rownames(T) = seq(1:1,nrow(T))


for(i in 0:18)
{
	win.graph();
	mhtplot(results[ results$CHR == i,c(2,3,4)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3))
	axis(2,at=0:10)
	abline(h=0)
	title(paste("TDT For Umbilical Hernia Disease Belgium CHR ",i,collapse = ""))
}

qqunif(results[,4],type="unif")
abline(h=3.5)


results <- read.table("TDT_6_X_v2.txt", header=F)
results = results[,c(2,1,3,7,4,5,8,9,10)]
chr = "X"

min(results[results$CHR == "2",4])
results = results[results$BP >= 47.000000,2:4] 
results = results[results$BP <= 101.000000,] 


mhtplot(results[,2:4],control = mht.control(usepos = FALSE,cutoffs = 3.5,colors=mycolors,yline=3,xline=3))
textvec <- as.character(c("X","Uwgs"))
xat <- locator()$x
xat <- round(xat)
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=0:10)
abline(h=0)
title("TDT For Cryptorchidism Disease CHR X")


setwd("C:\\Users\\Administrator\\Documents\\Visual Studio 2010\\Projects\\porc_Test\\Significant RG\\Chromosomes RG_V2\\");
significant_RG_v2 = results[-log10(results$P) >= 3.5 ,]
rownames(significant_RG_v2) = seq(1:1,nrow(significant_RG_v2))
write.table(significant_RG_v2, file = "significant_RG_v2.txt", sep  = "\t", quote = FALSE, row.names = FALSE);

0.25 / (8/100000)


##################### Manhattan Plots #########################
####################################################################################
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\");
results <- read.table("TDT_All.txt", header=F)
results = results[,c(1,2,3,7,4,5,8,9,10)]
colnames(results) = c("Marker_ID","CHR","BP","P","T","U","Nonsense","Conflict","Missing")
results$BP = results$BP / 1000000
results$logp = -log10(results$P)
win.graph()
qqunif(results$P)
title("QQ plot Cryptorchidism Disease Controls")
abline(h=5.5)
results2 = results[-log10(results$P) >= 5, ]
rownames(results2) = seq(1:1,nrow(results2))

results = results[results$CHR != 19 , ]
results2 = results2[-log10(results2$P) >= 3.5, ]
rownames(results2) = seq(1:1,nrow(results2))


r = results[,1:4]
write.table(results2, file = "CHR 14 Norway.txt", sep  = "\t", quote = FALSE, row.names = FALSE);

formate(6.94E-06)
format(6.94E-06, sci=F)
-log10(5.60E-05)

qvalue.gui()

(72-65)^2 / (65+72)

results2 = results2[results2$BP >= 7000000, ]
results2 = results2[results2$BP <= 95952589 + 100000, ]
cases = results2
controls = results2
All = results2

rownames(results2) = seq(1:1,nrow(results2))


mycolors <- c("red","blue","green","cyan","brown4","black","magenta","red","blue","green",
            "cyan","brown4","black","magenta","red","blue","green","cyan","black","green")
win.graph();
textvec <- c(as.character(0:18),"X")
textvec <- c(as.character(0:18))
mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=16,by=1))
abline(h=0)
title("TDT For Cryptorchidism Disease controls")

-log10(0.001090835)


xat <- locator()$x
xat <- round(xat)
abline(h=4.5,lty = 3)

setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Cryptorchidism\\plink\\Chromosomes\\");
for(i in 0:19)
{
	win.graph();
	jpeg(paste("CHR ",i,".jpg",collapse = ""), quality = 100)
	mhtplot(results[ results$CHR == i,c(1,3,10)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3))
	axis(2,at=0:14)
	abline(h=0)
	title(paste("TDT For Cryptorchidism Norway CHR ",i,collapse = ""))
	dev.off()
}


win.graph()
r = results[results$CHR != "19",]
qqunif(results$P)
title("QQ Plot TDT For Scrotal Cryptorchidism Disease DC")
abline(h=4.8)

####################################################################################

setwd("C:\\Users\\Administrator\\Documents\\Visual Studio 2010\\Projects\\porc_Test\\Backups")
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Cryptorchidism\\");
results <- read.table("TDT_cases.txt", header=F)
#results = results[,c(1,2,3,9,10,11)] # the last one is 6 or 9
results = results[,c(1,2,3,7,4,5,8,9,10)] # the last one is 6 or 9
colnames(results) = c("Marker_ID","CHR","BP","P","T","U","nonsense","Conflict","Missing")
results$BP = results$BP / 1000000
hist(results[results$Conflict > 0 , 8])
results = results[results$Conflict < 5 , ]

win.graph()
qqunif(results$P)
results2 = results[results$CHR == 23 , ]
rownames(results2) = seq(1:1,nrow(results2))
qqunif(results$P)
title("QQ Cryptorchidism Norway less than 5 conflict")



nrow(results)
nrow(results2)



mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=20,by=1))
abline(h=0)
title("TDT For Cryptorchidism Disease Norway")

results2 = results[results$CHR == 23, ]
results2 = results2[results2$BP >= 49.28828, ]
results2 = results2[results2$BP <= 51.11037, ]

-log10(results$P)



results = results[results$CHR == 14, ]
results = results[-log10(results$P) <= 10, ]
rownames(results) = seq(1:1,nrow(results))

####################################################################################
results <- read.table("f_pseudo_TDT_6.txt", header=F)
results = results[,c(2,1,3,15,12,13)] # the last one is 6 or 9
colnames(results) = c("Marker_ID","CHR","BP","P","T","U")


############################### Noregian Data ##########################
######### Cryptorchidism ##############
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Norway Data Set\\Cryptorchidism\\");
results <- read.table("allmarkers_PB_no.txt", header=T)
xat <- locator()$x
xat <- round(xat)
textvec <- c(as.character(0:18),"X")

mhtplot(results[,c(2:3,16)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=20,by=1))
abline(h=0)
title("TDT For Scrotal Hernia Disease Norway")
max(results$P1df)
win.graph();
qqunif(results$P1df)
title("QQ For Scrotal Hernia Disease Norway")


setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Umbilical\\plink\\");
results <- read.table("porc.tdt", header=T)
results = results[results$CHR == 2, ]
results = results[-log10(na.omit(results$P)) >= 8, ]
rownames(results) = seq(1:1,nrow(results))

win.graph()
qqunif(na.omit(results$P))
title("QQ For Umbilical Hernia Disease Norway data Plink")
xat <- locator()$x
xat <- round(xat)
textvec <- c(as.character(0:18),"X")


mhtplot(results[,c(1,3,10)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=20,by=1))
abline(h=0)
title("TDT For Umbilical Hernia Disease Norway cases Plink")



for(i in 0:19)
{
	win.graph();
	jpeg(paste("CHR ",i,".jpg",collapse = ""), quality = 100)
	mhtplot(results[ results$CHR == i,c(1,3,10)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3))
	axis(2,at=0:10)
	abline(h=0)
	#abline(v=116) #116
	title(paste("TDT For Cryptorchidism Disease Norway CHR ",i,collapse = ""))
	dev.off()
}


######################## Belgian Data ########################
##################### Manhattan Plots #########################
####################################################################################
setwd("C:\\Users\\Administrator\\Documents\\Visual Studio 2010\\Projects\\porc_Test\\Backups\\");
results <- read.table("TDT_2.txt", header=F)
results = results[,c(2,1,3,15,12,13,16,17,18)]
colnames(results) = c("Marker_ID","CHR","BP","P","T","U","Nonsense","Conflict","Missing")
results$BP = results$BP / 1000000
results$logp = -log10(results$P)
results = results[-log10(na.omit(results$P)) <= 22, ]
rownames(results) = seq(1:1,nrow(results))

xat <- locator()$x
xat <- round(xat)
textvec <- c(as.character(0:18),"X")



win.graph()
mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3),srt = 10)
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=23,by=1))
abline(h=0)
title("TDT For Cryptorchidism Disease Belgium")


win.graph()
qqunif(na.omit(results$P))
title("QQ For Cryptorchidism Disease Belgium")


results2 = results[results$CHR == 14 , ]
results2 = results2[results2$BP >= 49.28828, ]
results2 = results2[results2$BP <= 51.11037, ]
rownames(results2) = seq(1:1,nrow(results2))
write.table(results2, file = "CHR 14 Belgium.txt", sep  = "\t", quote = FALSE, row.names = FALSE);


##################### Pesudoautosomal Region Males and Females ##############
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\");
results <- read.table("TDT_cases_Male.txt", header=F)
#results = results[,c(1,2,3,9,10,11)] # the last one is 6 or 9
results = results[,c(1,2,3,7,4,5)] # the last one is 6 or 9
colnames(results) = c("Marker_ID","CHR","BP","P","T","U")
results$BP = results$BP / 1000000
win.graph()
qqunif(results$P)
results2 = results[results$CHR != 19 , ]
rownames(results2) = seq(1:1,nrow(results2))
qqunif(results2$P)



nrow(results)
nrow(results2)


win.graph()
mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=20,by=1))
abline(h=0)
title("TDT For Scrotal Hernia Disease Males")


write.table(results, file = "P_values.txt", sep  = "\t", quote = FALSE, row.names = FALSE);


###################### Norway Scrotal Hernia ########################

setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\");
results <- read.table("TDT_cases.txt", header=F)
results = results[,c(1,2,3,7,4,5,8,9,10)] # the last one is 6 or 9
colnames(results) = c("Marker_ID","CHR","BP","P","T","U","nonsense","Conflict","Missing")
results$BP = results$BP / 1000000
results = results[-log10(results$P) <= 20, ]
results = results[results$Conflict < 5 , ]
rownames(results) = seq(1:1,nrow(results))


win.graph()
hist(results[results$Conflict > 0 , 8])

win.graph()
qqunif(results$P)
results2 = results[results$CHR == 19 , ]
results2 =  results2[-log10(na.omit(results2$P)) >= 4, ]
rownames(results2) = seq(1:1,nrow(results2))

qqunif(results$P)
title("QQ Scrotal Hernia Norway less than 5 conflict")


xat <- locator()$x
xat <- round(xat)
textvec <- c(as.character(0:18),"X")
mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=20,by=1))
abline(h=0)
title("TDT For Scrotal Hernia Disease Norway")

setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\Chromosomes\\");
for(i in 0:19)
{
	win.graph();
	jpeg(paste("CHR ",i,".jpg",collapse = ""), quality = 100)
	mhtplot(results[ results$CHR == i,c(2,3,4)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3))
	axis(2,at=0:14)
	abline(h=0)
	title(paste("TDT For Scrotal Norway CHR ",i,collapse = ""))
	dev.off()
}

################## PLink for Scrotal ############################
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\plink\\");
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Belgium Data Set\\Cryptorchidism\\");
results <- read.table("Family.tdt", header=T)
results$BP = results$BP / 1000000
results = na.omit(results)
rownames(results) = seq(1:1,nrow(results))
nrow(results)
results = results[-log10(results$P) >= 10, ]
write.table(results, file = "Family.tdt", sep  = "\t", quote = FALSE, row.names = FALSE);


xat[12] = 33988 - 300
win.graph()
mhtplot(results[,c(1,3,10)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=30,by=1))
abline(h=0)
title("TDT For Scrotal Disease Belgium data plink")

qqunif(results$P)
title("QQ  Cryptorchidism Belgium data plink")


setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\plink\\Chromosomes\\");
for(i in 0:19)
{
	win.graph();
	jpeg(paste("CHR ",i,".jpg",collapse = ""), quality = 100)
	mhtplot(results[ results$CHR == i,c(1,3,10)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3))
	axis(2,at=0:14)
	abline(h=0)
	title(paste("TDT For Scrotal Norway plink CHR ",i,collapse = ""))
	dev.off()
}


################### Norway Cryptorchidism ########################
###################### Norway Scrotal Hernia ########################

rownamessetwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\");
results <- read.table("TDT_cases.txt", header=F)
results = results[,c(1,2,3,7,4,5,8,9,10)] # the last one is 6 or 9
colnames(results) = c("Marker_ID","CHR","BP","P","T","U","nonsense","Conflict","Missing")
results$BP = results$BP / 1000000
results = results[-log10(results$P) <= 20, ]
results = results[results$Conflict < 5 , ]
rownames(results) = seq(1:1,nrow(results))
write.table(results, file = "TDT_cases_M.txt", sep  = "\t", quote = FALSE, row.names = FALSE);



win.graph()
hist(results[results$Conflict > 0 , 8])

win.graph()
qqunif(results$P)
results = results[results$CHR != 19 , ]
results =  results2[-log10(na.omit(results2$P)) >= 4, ]
(results2) = seq(1:1,nrow(results2))
win.graph()
qqunif(results$P)
title("QQ Scrotal Hernia Norway with X")


xat <- locator()$x
xat <- round(xat)
textvec <- c(as.character(0:18))
mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=20,by=1))
abline(h=0)
title("TDT For Scrotal Hernia Disease Norway without X chr")

setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\Chromosomes\\");
for(i in 0:19)
{
	win.graph();
	jpeg(paste("CHR ",i,".jpg",collapse = ""), quality = 100)
	mhtplot(results[ results$CHR == i,c(2,3,4)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3))
	axis(2,at=0:14)
	abline(h=0)
	title(paste("TDT For Scrotal Norway CHR ",i,collapse = ""))
	dev.off()
}

setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\");
Female <- read.table("Female_modified.txt", header=F)
Male <- read.table("Male_modified.txt", header=F)

win.graph();
hist(Male$V4,xlab="Male Heterozygosity",main= "Male histogram")
hist(Female$V4,xlab="Female Heterozygosity",main= "Female histogram")

win.graph()
plot(x=seq(1:1:nrow(Male)), y =Male$V4, ylab="Male Heterozygosity", xlab = "Markers")
abline(v=141,col="BLACK",lty=2)
title("Male X chr")

win.graph()
plot(x=seq(1:1:nrow(Female)), y =Female$V4, ylab="Female Heterozygosity", xlab = "Markers")
abline(v=141,col="BLACK",lty=2)
title("Female X chr")


win.graph()
plot(x=Male$V3/1000000, y =Male$V4, ylab="Male Heterozygosity", xlab = "BP * million")
abline(v=7,col="BLACK",lty=2)
title("Male X chr")



win.graph()
plot(x=Female$V3/1000000, y =Female$V4, ylab="Male Heterozygosity", xlab = "BP * million")
abline(v=7,col="BLACK",lty=2)
title("Female X chr")



###################### Norway Umbilical Hernia ########################

setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Umbilical\\");
results <- read.table("TDT_cases_new.txt", header=F)
results = results[,c(1,2,3,7,4,5,8,9,10)] # the last one is 6 or 9
colnames(results) = c("Marker_ID","CHR","BP","P","T","U","nonsense","Conflict","Missing")
results$BP = results$BP / 1000000
results = results[-log10(results$P) >= 3.5, ]
results = results[results$Conflict < 5 , ]
rownames(results) = seq(1:1,nrow(results))
write.table(results, file = "TDT_cases_M.txt", sep  = "\t", quote = FALSE, row.names = FALSE);



win.graph()
hist(results[results$Conflict > 0 , 8])
max(results$Conflict)

win.graph()
qqunif(results$P)
results2 = results[results$CHR == 19 , ]
results2 =  results2[-log10(na.omit(results2$P)) >= 4, ]
rownames(results2) = seq(1:1,nrow(results2))

win.graph()
qqunif(results$P)
title("QQ Umbilical Hernia Norway less than 5 conflict")


xat <- locator()$x
xat <- round(xat)
textvec <- c(as.character(0:18),"X")
mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=20,by=1))
abline(h=0)
title("TDT For Umbilical Hernia Disease Norway")

setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Scrotal\\Chromosomes\\");
for(i in 0:19)
{
	win.graph();
	jpeg(paste("CHR ",i,".jpg",collapse = ""), quality = 100)
	mhtplot(results[ results$CHR == i,c(2,3,4)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3))
	axis(2,at=0:14)
	abline(h=0)
	title(paste("TDT For Scrotal Norway CHR ",i,collapse = ""))
	dev.off()
}


####################### Belgium and Norway Scrotal Hernia ##################
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Cryptorchidism\\");
results <- read.table("TDT_DC.txt", header=F)
results = results[,c(1,2,3,9,4,5,6,7)] # the last one is 6 or 9
colnames(results) = c("Marker_ID","CHR","BP","P","Tn","Un","Tb","Ub")
results$BP = results$BP / 1000000

win.graph()
qqunif(results$P)
title("QQ Scrotal Hernia joint analysis")

xat <- locator()$x
xat <- round(xat)
textvec <- c(as.character(0:18),"X")
mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=44,by=1))
abline(h=0)
title("TDT For Scrotal Hernia joint analysis")


results =  results[-log10(na.omit(results$P)) >= 10, ]
rownames(results2) = seq(1:1,nrow(results2))

(127 - 117)^2 / (117 + 127)


setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Cryptorchidism\\Chromosomes joint\\");
for(i in 19:19)
{
	win.graph();
	jpeg(paste("CHR ",i,".jpg",collapse = ""), quality = 100)
	mhtplot(results[ results$CHR == i,c(2,3,4)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3))
	axis(2,at=0:44)
	abline(h=0)
	title(paste("TDT For Cryptorchidism joint CHR ",i,collapse = ""))
	dev.off()
}

####################### Belgium and Cryptorchidism ##################
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Cryptorchidism\\Test\\");
results <- read.table("TDT_DC.txt", header=F)
results = results[,c(1,2,3,9,4,5,6,7)] # the last one is 6 or 9
colnames(results) = c("Marker_ID","CHR","BP","P","Tn","Un","Tb","Ub")
results$BP = results$BP / 1000000

win.graph()
qqunif(results$P)
title("QQ Cryptorchidism joint analysis")

xat <- locator()$x
xat <- round(xat)


textvec <- c(as.character(0:18),"X")
mhtplot(results[,2:4],control = mht.control(colors=mycolors,xline=2,labels=paste("chr",0:19,sep=""),srt=270))
mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=44,by=1))
abline(h=0)
title("TDT For Cryptorchidism joint analysis")

mhtplot



################# PAR in Belgium #######################
setwd("C:\\Users\\Administrator\\Documents\\Visual Studio 2010\\Projects\\porc_Test\\");
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\")
PAR = read.table("PAR.txt", header=T)
Female <- PAR[,c(1,2,4)]
Male <- PAR[,c(1,2,3)]
Male = na.omit(Male)
colnames(Male) = c("Marker_ID","BP","Perc")
rownames(Male) = seq(1:1,nrow(Male))

Markers_Norway =  Male[Male$BP >= 7000000,]
rownames(Markers_Norway) = seq(1:1,nrow(Markers_Norway))


Markers_Belgium =  Male[Male$BP >= 7000000,]
Markers_Belgium = Markers_Belgium[Markers_Belgium$Perc >= 0.1,]
rownames(Markers_Belgium) = seq(1:1,nrow(Markers_Belgium))

Markers_Norway
Markers_Belgium


Female = na.omit(Female)
colnames(Female) = c("Marker_ID","BP","Perc")
rownames(Female) = seq(1:1,nrow(Female))


win.graph();
hist(Male$V4,xlab="Male Heterozygosity",main= "Male histogram")
hist(Female$V4,xlab="Female Heterozygosity",main= "Female histogram")

win.graph()
plot(x=seq(1:1:nrow(Male)), y =Male$Perc, ylab="Male Heterozygosity", xlab = "Markers")
abline(v=170,col="BLACK",lty=2)
title("Male X chr")


win.graph()
plot(x=seq(1:1:nrow(Female)), y =Female$Perc, ylab="Female Heterozygosity", xlab = "Markers")
abline(v=170,col="BLACK",lty=2)
title("Female X chr")


win.graph()
plot(x=Male$BP/1000000, y =Male$Perc, ylab="Male Heterozygosity", xlab = "BP * million",xaxt='n')
abline(v=7,col="BLACK",lty=2)
axis(1,at=seq(from=0,to=150,by=20))
title("Male X chr Belgium")



win.graph()
plot(x=Female$BP/1000000, y =Female$Perc, ylab="Male Heterozygosity", xlab = "BP * million",xaxt='n')
abline(v=7,col="BLACK",lty=2)
axis(1,at=seq(from=0,to=150,by=20))
title("Female X chr Belgium")


############## Belgium data Umbilical Hernia ##################
setwd("C:\\Users\\Administrator\\Documents\\Visual Studio 2010\\Projects\\porc_Test\\");
results <- read.table("TDT_2_v2.txt", header=F)
results = results[,c(2,1,3,7,4,5,8,9,10)]
colnames(results) = c("Marker_ID","CHR","BP","P","T","U","Nonsense","Conflict","Missing")
results$BP = results$BP / 1000000
results = results[-log10(results$P) <= 20,]

win.graph()
qqunif(results$P)
title("QQ Scrotal Hernia belgium data")

mycolors <- c("red","blue","green","cyan","brown4","black","magenta","red","blue","green",
            "cyan","brown4","black","magenta","red","blue","green","cyan","black","red")

xat <- locator()$x
xat <- round(xat)
textvec <- c(as.character(0:18),"X")
mhtplot(results[,2:4],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=44,by=1))
abline(h=0)
title("TDT For Cryptorchidism Belgium Data")

install.packages("ggbio")
library(ggbio)
data(hg19IdeogramCyto, package = "biovizBase")
data(hg19Ideogram, package = "biovizBase")
install.packages(GenomicRanges)
library(GenomicRanges)

###########################################
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Combine P values\\Umbilical\\TDT_total")
results <- read.table("TDT_DC.txt", header=F)

win.graph();
qqunif(results$V11);
title("QQ Umbilical Hernia belgium and Norway TDT")
mhtplot(results[,c(2,3,11)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=44,by=1))
abline(h=0)
title("TDT For Umbilical Hernia belgium and Norway TDT")





setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\All Data\\NO plink BE mine\\")
setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Combine P values\\Scrotal\\TDT_total\\");
results <- read.table("SR_joint.txt", header=T)
results$

results = results[,c(1,2,3,7,4,5,8,9,10)]
colnames(results) = c("Marker_ID","CHR","BP","P","T","U","Nonsense","Conflict","Missing")
results$BP = results$BP / 1000000

results = na.omit(results)
results = results[results$CHR != "19",]

win.graph()
qqunif(results$P)
title("QQ RG + SR Data without X CHR")

xat <- locator()$x
xat <- round(xat)
textvec <- c(as.character(0:18),"X")
mhtplot(results[,c(1,3,9)],control = mht.control(usepos = FALSE,colors=mycolors,yline=3,xline=3,cex = 3),srt = 10,pch = ".")
mtext(textvec, side=1, line=0, at=xat)
axis(2,at=seq(from=0,to=65,by=1))
abline(h=0)
title("TDT For RG + SR ")

for(i in 14:15)
{
	win.graph();
	jpeg(paste("CHR ","X",".jpg",collapse = ""), quality = 100)
	mhtplot(results[results$CHR == 19 ,c(2,3,11)],control = mht.control(usepos = TRUE,colors=mycolors,yline=3,xline=3,cex = 4),srt = 10,pch = ".")
	axis(2,at=0:65)
	abline(h=0)
	title(paste("TDT For Scrotal hernia combined data CHR ","X",collapse = ""))
	dev.off()
}

(135-45)^2 / (135+45)





setwd("D:\\Master of Bioinformatics\\GIGA Lab\\Norway day\\Til Mahmoud\\Cryptorchidism\\Sex SNP\\");
Males <- read.table("Males_RG_Eli.txt", header=F)
Females <- read.table("Females_RG_Eli.txt", header=F)

Males = Males[!is.nan(Males$V2),]
Females = Females[!is.nan(Females$V2),]

win.graph();
hist(Males$V2,xlab="Male Heterozygosity",main= "Male histogram RG")
win.graph();
hist(Females$V2,xlab="Female Heterozygosity",main= "Female histogram RG")


(1.4+0.8-0.1-0.4)/(0.7^2 + 0.4^2 + 0.1^2 + 0.4^2 + 0.5^2)

(1.588785*0.4) + 1