
table=htmltab("1543.html",which=3)
genes1543=table[11]

Reduce(intersect, list(genesstrongFirst,genes37,genes15,genes43,genes4315,genes1543))

vlist=list(genes607,genesstrongFirst,genes37,genes15,genes43,genes4315,genes1543)
countertable=addmargins(table(gene=unlist(vlist), vec=rep(c("genes607","genesstrongFirst","genes37","genes15","genes43","genes4315","genes1543"),times=sapply(vlist,length))),2,list(Count=function(x) sum(x[x>0])))


vlist=list(genes37,genes15)
addmargins(table(gene=unlist(vlist), vec=rep(paste0("v",1:2),times=sapply(vlist,length))),2,list(Count=function(x) sum(x[x>0])))

colnames(genes1543) <- NULL
colnames(genes607) <- NULL
colnames(genes4315) <- NULL
colnames(genesstrongFirst) <- NULL
colnames(genes43) <- NULL
colnames(genes37) <- NULL
colnames(genes15) <- NULL

genes1543=unique(genes1543) 
genes607=unique(genes607) 
genes4315=unique(genes4315) 
genesstrongFirst=unique(genesstrongFirst) 
genes43=unique(genes43) 
genes37=unique(genes37)
genes15=unique(genes15) 


genes1543=unlist(genes1543) 
genes607=unlist(genes607) 
genes4315=unlist(genes4315) 
genesstrongFirst=unlist(genesstrongFirst) 
genes43=unlist(genes43) 
genes37=unlist(genes37)
genes15=unlist(genes15) 


c("genes607","genesstrongFirst","genes37","genes15","genes43","genes4315","genes1543")


