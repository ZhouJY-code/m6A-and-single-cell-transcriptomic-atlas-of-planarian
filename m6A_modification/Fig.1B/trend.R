#’€œﬂÕº
par(mfrow=c(3,2),mar=c(3,4,1,1))
all<-read.table("E:/1f-7677.txt",sep="\t",stringsAsFactors=F,row.names=1,header=F);
sp=spline(c(1:5),c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),n=100)
plot(sp,col="deepskyblue",ylab="Normalized expression",ylim=c(-2,2),las=1,xlab="",cex.lab=1.5,cex.axis=1.5,type="l",lwd=3);

all<-read.table("E:/2f-3363.txt",sep="\t",stringsAsFactors=F,row.names=1,header=F);
sp=spline(c(1:5),c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),n=100)
plot(sp,col="deepskyblue",ylab="Normalized expression",ylim=c(-2,2),las=1,xlab="",cex.lab=1.5,cex.axis=1.5,type="l",lwd=3);

all<-read.table("E:/3f-6853.txt",sep="\t",stringsAsFactors=F,row.names=1,header=F);
sp=spline(c(1:5),c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),n=100)
plot(sp,col="deepskyblue",ylab="Normalized expression",ylim=c(-2,2),las=1,xlab="",cex.lab=1.5,cex.axis=1.5,type="l",lwd=3);

all<-read.table("E:/4f-6313.txt",sep="\t",stringsAsFactors=F,row.names=1,header=F);
sp=spline(c(1:5),c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),n=100)
plot(sp,col="deepskyblue",ylab="Normalized expression",ylim=c(-2,2),las=1,xlab="",cex.lab=1.5,cex.axis=1.5,type="l",lwd=3);

all<-read.table("E:/5f-7334.txt",sep="\t",stringsAsFactors=F,row.names=1,header=F);
sp=spline(c(1:5),c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),n=100)
plot(sp,col="deepskyblue",ylab="Normalized expression",ylim=c(-2,2),las=1,xlab="",cex.lab=1.5,cex.axis=1.5,type="l",lwd=3);

all<-read.table("E:/6f-3672.txt",sep="\t",stringsAsFactors=F,row.names=1,header=F);
sp=spline(c(1:5),c(median(t(apply(all,1,scale))[,1]),median(t(apply(all,1,scale))[,2]),median(t(apply(all,1,scale))[,3]),median(t(apply(all,1,scale))[,4]),median(t(apply(all,1,scale))[,5])),n=100)
plot(sp,col="deepskyblue",ylab="Normalized expression",ylim=c(-2,2),las=1,xlab="",cex.lab=1.5,cex.axis=1.5,type="l",lwd=3);
