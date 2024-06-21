number<- c(8298, 9443, 6164,10374,7289)
number1<- c(9057, 12779, 11224, 13362,12362)
names<-c("0d","6h","3d","7d","11d")
barplot(number1,names.arg=names,main="",xlab="stage",ylab="number",las=1,col="#BC3C28",cex.axis=1.2,lwd=5,ylim=c(0,16000))
