#’€œﬂÕº
par(mfrow=c(2,2),mar=c(3,4,1,1))
all<-read.table("E:/1-p.txt",sep="\t",stringsAsFactors=F,row.names=1,header=F);
all<-log2(all)
matplot(apply(all,1,scale),type="l",lty=1,col="deepskyblue",ylab="Normalized expression",ylim=c(-2,2),las=1,xaxt="n",cex.lab=1.5,cex.axis=1.5);
#axis(1,at=1:5,labels=c("0h","6h","3d","7d","11d"),cex.axis=1.5);
lines(c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),col="black",lwd=3);
text(1,1.5,"gene:737",pos=4,cex=3);
all<-read.table("E:/2-p.txt",sep="\t",stringsAsFactors=F,row.names=1,header=T);
all<-log2(all)
matplot(apply(all,1,scale),type="l",lty=1,col="orange",ylab="",ylim=c(-2,2),las=1,xaxt="n",cex.lab=1.5,cex.axis=1.5);
#axis(1,at=1:5,labels=c("0h","6h","3d","7d","11d"),cex.axis=1.5);
lines(c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),col="black",lwd=3);
text(1,1.5,"gene:893",pos=4,cex=3);
all<-read.table("E:/3-p.txt",sep="\t",stringsAsFactors=F,row.names=1,header=T);
all<-log2(all)
matplot(apply(all,1,scale),type="l",lty=1,col="purple",ylab="Normalized expression",ylim=c(-2,2),las=1,xaxt="n",cex.lab=1.5,cex.axis=1.5);
axis(1,at=1:5,labels=c("0h","6h","3d","7d","11d"),cex.axis=1.5);
lines(c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),col="black",lwd=3);
text(1,1.5,"gene:2779",pos=4,cex=3);
all<-read.table("E:/4-p.txt",sep="\t",stringsAsFactors=F,row.names=1,header=F);
all<-log2(all)
matplot(apply(all,1,scale),type="l",lty=1,col="pink",ylab="",ylim=c(-2,2),las=1,xaxt="n",cex.lab=1.5,cex.axis=1.5);
axis(1,at=1:5,labels=c("0h","6h","3d","7d","11d"),cex.axis=1.5);
lines(c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),col="black",lwd=3);
text(1,1.5,"gene:1976",pos=4,cex=3);
