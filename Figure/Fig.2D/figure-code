###figure1  numbaer

barplot(c(11635,14427),names.arg = c("Ctrl-0h","kd-0h"),col=c("red","blue"),las=1,cex.axis=1.5,cex.names = 2)

barplot(c(15308,13818),names.arg = c("Ctrl-6h","kd-6h"),col=c("red","blue"),las=1,cex.axis=1.5,cex.names = 2)

barplot(c(14918,15105),names.arg = c("Ctrl-3d","kd-3d"),col=c("red","blue"),las=1,cex.axis=1.5,cex.names = 2)

barplot(c(16289 ,12038),names.arg = c("Ctrl-7d","kd-7d"),col=c("red","blue"),las=1,cex.axis=1.5,cex.names = 2)

barplot(c(16551,12192),names.arg = c("Ctrl-11d","kd-11d"),col=c("red","blue"),las=1,cex.axis=1.5,cex.names = 2)

###figure2 level

a<-read.table("E:/total-Ctrl-11d.txt")
b<-read.table("E:/total-kd-11d.txt")


par(mfrow=c(1,2))
boxplot(log2(a[,7]),log2(b[,7]),col=2:3,names = c("ctrl","kd"),cex.axis=2,las=1)

plot(ecdf(log2(a[,7])),col=2,lwd=2,main="",xlab="Level",ylab="Density",cex.axis=1.5,las=1,cex.lab=1.5)
lines(ecdf(log2(b[,7])),col=3,lwd=2)
legend("right",c("ctrl","kd"),col=c(2,3),lty=c(1,1))

ks.test(log2(a[,7]),log2(b[,7]))
mean(log2(a[,7]))
mean(log2(b[,7]))


###figure3 »ðÉ½Í¼
###0h
a<-read.table("E:/total-0h")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))

plot(a$fold, -log10(a$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p)", las = 1, main = "Volcano Plot",xlim=c(-5,5),cex.axis=1.5,cex.lab=1.5)
points(a[a[,15] >= 0.26 & a[,16]<0.05,15], -log10(a[a[,15] >= 0.26 & a[,16]<0.05,16]), pch=16, col="red")
points(a[a[,15] <= -0.26 & a[,16]<0.05,15], -log10(a[a[,15] <= -0.26 & a[,16]<0.05,16]), pch=16, col="green")

length(a[a[,15] >= 0.26 & a[,16]<0.05,15])
length(a[a[,15] <= -0.26 & a[,16]<0.05,15])

abline(v=0.26,lty=8)
abline(v=-0.26,lty=8)
abline(h=-log10(0.05),lty=8)

text(4,4, "n=752",cex=2)
text(-4,4, "n=714",cex=2)

write.table(a, "E:/0h-new.txt",sep="\t",quot=F)

###6h
a<-read.table("E:/total-6h")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))

plot(a$fold, -log10(a$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p)", las = 1, main = "Volcano Plot",xlim=c(-5,5),cex.axis=1.5,cex.lab=1.5)
points(a[a[,15] >= 0.26 & a[,16]<0.05,15], -log10(a[a[,15] >= 0.26 & a[,16]<0.05,16]), pch=16, col="red")
points(a[a[,15] <= -0.26 & a[,16]<0.05,15], -log10(a[a[,15] <= -0.26 & a[,16]<0.05,16]), pch=16, col="green")

length(a[a[,15] >= 0.26 & a[,16]<0.05,15])
length(a[a[,15] <= -0.26 & a[,16]<0.05,15])

abline(v=0.26,lty=8)
abline(v=-0.26,lty=8)
abline(h=-log10(0.05),lty=8)

text(4,4, "n=549",cex=2)
text(-4,4, "n=632",cex=2)

write.table(a, "E:/6h-new.txt",sep="\t",quot=F)

###3d
a<-read.table("E:/total-3d")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))

plot(a$fold, -log10(a$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p)", las = 1, main = "Volcano Plot",xlim=c(-5,5),cex.axis=1.5,cex.lab=1.5)
points(a[a[,15] >= 0.26 & a[,16]<0.05,15], -log10(a[a[,15] >= 0.26 & a[,16]<0.05,16]), pch=16, col="red")
points(a[a[,15] <= -0.26 & a[,16]<0.05,15], -log10(a[a[,15] <= -0.26 & a[,16]<0.05,16]), pch=16, col="green")

length(a[a[,15] >= 0.26 & a[,16]<0.05,15])
length(a[a[,15] <= -0.26 & a[,16]<0.05,15])

abline(v=0.26,lty=8)
abline(v=-0.26,lty=8)
abline(h=-log10(0.05),lty=8)

text(4,6, "n=839",cex=2)
text(-4,6, "n=2331",cex=2)

write.table(a, "E:/3d-new.txt",sep="\t",quot=F)

###7d
a<-read.table("E:/total-7d")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))

plot(a$fold, -log10(a$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p)", las = 1, main = "Volcano Plot",xlim=c(-5,5),cex.axis=1.5,cex.lab=1.5)
points(a[a[,15] >= 0.26 & a[,16]<0.05,15], -log10(a[a[,15] >= 0.26 & a[,16]<0.05,16]), pch=16, col="red")
points(a[a[,15] <= -0.26 & a[,16]<0.05,15], -log10(a[a[,15] <= -0.26 & a[,16]<0.05,16]), pch=16, col="green")

length(a[a[,15] >= 0.26 & a[,16]<0.05,15])
length(a[a[,15] <= -0.26 & a[,16]<0.05,15])

abline(v=0.26,lty=8)
abline(v=-0.26,lty=8)
abline(h=-log10(0.05),lty=8)

text(4,4, "n=764",cex=2)
text(-4,4, "n=3363",cex=2)

write.table(a, "E:/7d-new.txt",sep="\t",quot=F)

###11d
a<-read.table("E:/total-11d")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))

plot(a$fold, -log10(a$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p)", las = 1, main = "Volcano Plot",xlim=c(-5,5),cex.axis=1.5,cex.lab=1.5)
points(a[a[,15] >= 0.26 & a[,16]<0.05,15], -log10(a[a[,15] >= 0.26 & a[,16]<0.05,16]), pch=16, col="red")
points(a[a[,15] <= -0.26 & a[,16]<0.05,15], -log10(a[a[,15] <= -0.26 & a[,16]<0.05,16]), pch=16, col="green")

length(a[a[,15] >= 0.26 & a[,16]<0.05,15])
length(a[a[,15] <= -0.26 & a[,16]<0.05,15])

abline(v=0.26,lty=8)
abline(v=-0.26,lty=8)
abline(h=-log10(0.05),lty=8)

text(4,4, "n=1202",cex=2)
text(-4,4, "n=1973",cex=2)

write.table(a, "E:/11d-new.txt",sep="\t",quot=F)



###figure3 deg-level

a<-read.table("E:/m6A-fc")
b<-read.table("E:/non-m6A-fc")


par(mfrow=c(1,2))
boxplot(a[,1],b[,1],col=2:3,names = c("m6A","non-m6A"),cex.axis=2,las=1)

plot(ecdf(a[,1]),col=2,lwd=2,main="",xlab="Level",ylab="Density",cex.axis=1.5,las=1,cex.lab=1.5)
lines(ecdf(b[,1]),col=3,lwd=2)
legend("right",c("m6A","non-m6A"),col=c(2,3),lty=c(1,1))

ks.test(a[,1],b[,1])
mean(a[,1])
mean(b[,1])
