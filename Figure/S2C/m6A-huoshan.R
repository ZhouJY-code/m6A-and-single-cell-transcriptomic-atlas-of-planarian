library(ggplot2)

###0h
data<-read.table("E:/0h-new.txt",header=T,row.names=1)


data$change<- as.factor(ifelse(data$pvalue <0.05 & abs(data$fold)>0.26,ifelse(data$fold >0.26, 'UP','DOWN'),'NOT'))

p<-ggplot(data = data, aes(x = fold, y = -log10(pvalue),color = change)) +
  geom_point(alpha=0.8, size = 5) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +xlim(-5,5)+
  geom_hline(yintercept=1.301 ,linetype=2 ,size=1) +
geom_vline(xintercept=c(-0.26,0.26) ,linetype=2 ,size=1,) +
  scale_color_manual(name = "", values = c("#0072B5","#BC3C28","grey"), limits = c("UP", "DOWN", "NOT")) +theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p +theme(axis.ticks.y=element_line(color="black",size=2))+theme(axis.ticks.x=element_line(color="black",size=2))

###6h
data<-read.table("E:/6h-new.txt",header=T,row.names=1)


data$change<- as.factor(ifelse(data$pvalue <0.05 & abs(data$fold)>0.26,ifelse(data$fold >0.26, 'UP','DOWN'),'NOT'))

p<-ggplot(data = data, aes(x = fold, y = -log10(pvalue),color = change)) +
  geom_point(alpha=0.8, size = 5) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +xlim(-5,5)+
  geom_hline(yintercept=1.301 ,linetype=2 ,size=1) +
geom_vline(xintercept=c(-0.26,0.26) ,linetype=2 ,size=1,) +
  scale_color_manual(name = "", values = c("#0072B5","#BC3C28","grey"), limits = c("UP", "DOWN", "NOT")) +theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p +theme(axis.ticks.y=element_line(color="black",size=2))+theme(axis.ticks.x=element_line(color="black",size=2))


###3d
data<-read.table("E:/3d-new.txt",header=T,row.names=1)


data$change<- as.factor(ifelse(data$pvalue <0.05 & abs(data$fold)>0.26,ifelse(data$fold >0.26, 'UP','DOWN'),'NOT'))

p<-ggplot(data = data, aes(x = fold, y = -log10(pvalue),color = change)) +
  geom_point(alpha=0.8, size = 5) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +xlim(-5,5)+
  geom_hline(yintercept=1.301 ,linetype=2 ,size=1) +
geom_vline(xintercept=c(-0.26,0.26) ,linetype=2 ,size=1,) +
  scale_color_manual(name = "", values = c("#0072B5","#BC3C28","grey"), limits = c("UP", "DOWN", "NOT")) +theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p +theme(axis.ticks.y=element_line(color="black",size=2))+theme(axis.ticks.x=element_line(color="black",size=2))

###7d
data<-read.table("E:/7d-new.txt",header=T,row.names=1)


data$change<- as.factor(ifelse(data$pvalue <0.05 & abs(data$fold)>0.26,ifelse(data$fold >0.26, 'UP','DOWN'),'NOT'))

p<-ggplot(data = data, aes(x = fold, y = -log10(pvalue),color = change)) +
  geom_point(alpha=0.8, size = 5) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +xlim(-5,5)+
  geom_hline(yintercept=1.301 ,linetype=2 ,size=1) +
geom_vline(xintercept=c(-0.26,0.26) ,linetype=2 ,size=1,) +
  scale_color_manual(name = "", values = c("#0072B5","#BC3C28","grey"), limits = c("UP", "DOWN", "NOT")) +theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p +theme(axis.ticks.y=element_line(color="black",size=2))+theme(axis.ticks.x=element_line(color="black",size=2))


###11d
data<-read.table("E:/11d-new.txt",header=T,row.names=1)


data$change<- as.factor(ifelse(data$pvalue <0.05 & abs(data$fold)>0.26,ifelse(data$fold >0.26, 'UP','DOWN'),'NOT'))

p<-ggplot(data = data, aes(x = fold, y = -log10(pvalue),color = change)) +
  geom_point(alpha=0.8, size = 5) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +xlim(-5,5)+
  geom_hline(yintercept=1.301 ,linetype=2 ,size=1) +
geom_vline(xintercept=c(-0.26,0.26) ,linetype=2 ,size=1,) +
  scale_color_manual(name = "", values = c("#0072B5","#BC3C28","grey"), limits = c("UP", "DOWN", "NOT")) +theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p +theme(axis.ticks.y=element_line(color="black",size=2))+theme(axis.ticks.x=element_line(color="black",size=2))
