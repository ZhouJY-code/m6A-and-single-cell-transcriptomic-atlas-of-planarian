number<- c(80.8868292,15.7296474,2.6904241,100-80.8868292-15.7296474-2.6904241)
number1<- c(80.9,15.7,2.7,0.7)
names<-c("n=1","n=2","n=3","n>=4")
barplot(number1,names.arg=names,main="",xlab="",ylab="",las=1,col="#BC3C28",cex.axis=1.2,lwd=5,ylim=c(0,100))
