##############################################
control_utr5<-read.table("E:/6h-m6a-utr5",stringsAsFactors=F);
control_cds<-read.table("E:/6h-m6a-cds",stringsAsFactors=F);
control_utr3<-read.table("E:/6h-m6a-utr3",stringsAsFactors=F);

control_utr5[,13]<-(control_utr5[,2]+ control_utr5[,3])/2;
control_utr5[,14]<-control_utr5[,10]+(control_utr5[,13]-control_utr5[,6]+1)/control_utr5[,9];


control_cds[,13]<-(control_cds[,2]+ control_cds[,3])/2;
control_cds[,14]<-control_cds[,10]+(control_cds[,13]-control_cds[,6]+1)/control_cds[,9];

control_utr3[,13]<-(control_utr3[,2]+ control_utr3[,3])/2;
control_utr3[,14]<-control_utr3[,10]+(control_utr3[,13]-control_utr3[,6]+1)/control_utr3[,9];


control_utr5[,15]<-ceiling(control_utr5[,14]* 235);
control_cds[,15]<-ceiling(control_cds[,14]* 1119)+ 235;
control_utr3[,15]<-ceiling(control_utr3[,14]* 260)+235+ 1119;

x.control<-as.numeric(names(table(c(control_utr5[,15], control_cds[,15], control_utr3[,15]))));
y.control<-as.numeric(table(c(control_utr5[,15], control_cds[,15], control_utr3[,15])));

plot(-1,type="l",lwd=5,col="deepskyblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0, 1614),ylim=c(0,40));
lines(smooth.spline(x.control, y.control,df=50)$x,smooth.spline(x.control, y.control,df=50)$y,lwd=5,col="firebrick2",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
abline(v= 235.5,lwd=5,lty=2);
abline(v= 235 + 1080 + 0.5,lwd=5,lty=2); ##长度不对，应为1119


control_utr5<-read.table("E:/6h-kd-m6a-utr5",stringsAsFactors=F);
control_cds<-read.table("E:/6h-kd-m6a-cds",stringsAsFactors=F);
control_utr3<-read.table("E:/6h-kd-m6a-utr3",stringsAsFactors=F);

control_utr5[,13]<-(control_utr5[,2]+ control_utr5[,3])/2;
control_utr5[,14]<-control_utr5[,10]+(control_utr5[,13]-control_utr5[,6]+1)/control_utr5[,9];


control_cds[,13]<-(control_cds[,2]+ control_cds[,3])/2;
control_cds[,14]<-control_cds[,10]+(control_cds[,13]-control_cds[,6]+1)/control_cds[,9];

control_utr3[,13]<-(control_utr3[,2]+ control_utr3[,3])/2;
control_utr3[,14]<-control_utr3[,10]+(control_utr3[,13]-control_utr3[,6]+1)/control_utr3[,9];


control_utr5[,15]<-ceiling(control_utr5[,14]* 235);
control_cds[,15]<-ceiling(control_cds[,14]* 1119)+ 235;
control_utr3[,15]<-ceiling(control_utr3[,14]* 260)+235+ 1119;

x.control<-as.numeric(names(table(c(control_utr5[,15], control_cds[,15], control_utr3[,15]))));
y.control<-as.numeric(table(c(control_utr5[,15], control_cds[,15], control_utr3[,15])));
lines(smooth.spline(x.control, y.control,df=50)$x,smooth.spline(x.control, y.control,df=50)$y,lwd=5,col="deepskyblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
abline(v= 235.5,lwd=5,lty=2);
abline(v= 235 + 1080 + 0.5,lwd=5,lty=2); ##长度不对，应为1119


axis(1,at=c(0, 235.5, 235 + 1080 + 0.5, 1575),label=rep("",4),lwd=5,tcl=-1.2);##长度不对，应为1614
axis(2,at=c(0,20,40),label=rep("",3),lwd=5,tcl=-1.2);
