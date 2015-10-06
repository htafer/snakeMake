#!/usr/bin/env Rscript                                                                                                                                                                                                                                                                                                                                                                        
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-c","--countFile"),
                help="Count file as generated for example by featureCounts."
                ),
    make_option(c("-p","--pValue"),
                help="Pvalue threshold for a gene to be reported as differentially expressed",
                default=0.1
                ),
    make_option(c("-l","--logFC"),
                help="Absolute Value of the log Fold change",
                default=0
                ),
    make_option(c("-n","--numberBack"),
                help="Number of differentially expressed genes to return",
                default=1000)
)
opt<- parse_args(OptionParser(option_list=option_list))

#define contrast matrix                                                                                                                                                                                                                                                                                                                                                                       
getContrastDesign<-function(classes){
    list<-combn(classes,2,simplify=F);
    contrast<-vector();
    for(i in 1:length(list)){
        str=paste(list[[i]][1],list[[i]][2],sep="-");
        contrast<-c(contrast,str);
        str=paste(list[[i]][2],list[[i]][1],sep="-");
        contrast<-c(contrast,str);
    }
    return(contrast);
}


library(edgeR)
library(limma)
library(stringr)
library(statmod)
#read data                                                                                                                                                                                                                                                                                                                                                                                    
datafile=read.table(opt$countFile,header=T,row.names=1)
#get name of categories                                                                                                                                                                                                                                                                                                                                                                       
groupName=str_replace(colnames(datafile),"_.+","" )
classes<-unique(groupName)
#Do what do below but by comparing i,j                                                                                                                                                                                                                                                                                                                                                        
lC=length(classes)
lC_1<-lC-1

for(a in 1:lC_1){
    c=a+1;
    for(b in c:lC){
        class<-c(classes[a],classes[b]);
        columnsa<-grep(paste(classes[a],"_",sep=""),colnames(datafile))
        columnsb<-grep(paste(classes[b],"_",sep=""),colnames(datafile))
        columns<-c(columnsa,columnsb)
        gName<-groupName[columns]
        all<-DGEList(data.matrix(datafile[,columns]),group=gName)
        isexpr<-rowSums(cpm(all)>1)>=length(gName)
        all<-all[isexpr,,keep.lib.sizes=FALSE]
        all<-calcNormFactors(all)
        design<-model.matrix(~ 0+factor(match(gName,class)))
        colnames(design)<-unique(gName)
        v<-voom(all,design,plot=F)
        fit <- lmFit(v, design)
        #Contrast Design                                                                                                                                                                                                                                                                                                                                                                      
        contrasts<-getContrastDesign(class)
        contrast.matrix<-makeContrasts(contrasts=contrasts,levels=design)
        #Fit + predictions                                                                                                                                                                                                                                                                                                                                                                    
        fit2<-contrasts.fit(fit,contrast.matrix)
        fit2<-eBayes(fit2)
        #Now extract results                                                                                                                                                                                                                                                                                                                                                                  
        for( i in contrasts){
            fra<-topTable(fit2,coef=i,number=Inf,adjust.method="fdr",lfc=opt$logFC,p.value=opt$pValue,sort.by="P")
            write.table(head(fra[fra$logFC <=  opt$logFC   & fra$adj.P.Val < opt$pValue,],n=opt$numberBack),file=paste(i,"DOWN",sep="."),qmethod="double")
            write.table(head(fra[fra$logFC >=  opt$logFC   & fra$adj.P.Val < opt$pValue,],n=opt$numberBack),file=paste(i,"UP",sep=".")  ,qmethod="double")
        }
    }
}
