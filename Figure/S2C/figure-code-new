###进入R中
###0h
a<-read.table("./total-0h")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))
write.table(a, "./0h-new.txt",sep="\t",quot=F)

###6h
a<-read.table("./total-6h")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))
write.table(a, "./6h-new.txt",sep="\t",quot=F)

###3d
a<-read.table("./total-3d")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))
write.table(a, "./3d-new.txt",sep="\t",quot=F)

###7d
a<-read.table("./total-7d")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))
write.table(a, "./7d-new.txt",sep="\t",quot=F)

###11d
a<-read.table("./total-11d")
a$fold<-log2(a[,14]/a[,7])
a$pvalue <- as.numeric(apply(a[,1:14],1,function(x){t.test(as.numeric(x[4:6]),as.numeric(x[11:13]),paired=F)$p.value}))
write.table(a, "./11d-new.txt",sep="\t",quot=F)


#出R
cat 0h-new.txt|awk '$16>=0.26 && $17<0.05'|wc
cat 0h-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
cat 6h-new.txt|awk '$16>=0.26 && $17<0.05'|wc
cat 6h-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
cat 3d-new.txt|awk '$16>=0.26 && $17<0.05'|wc
cat 3d-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
cat 7d-new.txt|awk '$16>=0.26 && $17<0.05'|wc
cat 7d-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
cat 11d-new.txt|awk '$16>=0.26 && $17<0.05'|wc
cat 11d-new.txt|awk '$16<=-0.26 && $17<0.05'|wc



[sunbf@logina1 deg]$ cat 11d-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
    556    9451   84003
[sunbf@logina1 deg]$ cat 0h-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
    520    8840   78467
[sunbf@logina1 deg]$ cat 6h-new.txt|awk '$16>=0.26 && $17<0.05'|wc
    425    7224   64374
[sunbf@logina1 deg]$ cat 6h-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
    405    6885   61180
[sunbf@logina1 deg]$ cat 3d-new.txt|awk '$16>=0.26 && $17<0.05'|wc
    488    8295   73340
[sunbf@logina1 deg]$ cat 3d-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
   1263   21471  191721
[sunbf@logina1 deg]$ cat 7d-new.txt|awk '$16>=0.26 && $17<0.05'|wc
    560    9519   84230
[sunbf@logina1 deg]$ cat 7d-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
   2162   36754  330001
[sunbf@logina1 deg]$ cat 11d-new.txt|awk '$16>=0.26 && $17<0.05'|wc
    538    9145   81048
[sunbf@logina1 deg]$ cat 11d-new.txt|awk '$16<=-0.26 && $17<0.05'|wc
    778   13226  118351
