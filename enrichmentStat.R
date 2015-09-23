#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
#library("RColorBrewer")
#library("cummeRbund")
#library("VennDiagram")
#library("rPlotter")
#library("EBImage")
#library("devtools")
#library("Category")

option_list <- list(
    make_option(c("-b","--background"),
                help="Background data containing the gene->{GO,PFAM,KEGG,relationship}
                     It must contain a header with name IDs and ANN"),
    make_option(c("-d","--geneFile"),
                help="File contains the genes for gene set enrichment has to be computed"),
    make_option(c("-t","--type"),
                help="type of data: GO or anything else"),
    make_option(c("-p","--pValue"),
                default="0.1",
                help="pValue threshold"),
    make_option(c("-D","--dictionary"),
                default="NA",
                help="dictionary for ID description conversion")
)
opt<- parse_args(OptionParser(option_list=option_list))

bla<-list()
################################################################################
#Function to get the upregulated genes  corresponding to the enriched function
################################################################################



getCorrespondingGene<-function(el){
    IPR<-read.table(opt$background,header=T)
    geneUP<-read.table(opt$geneFile,header=T,row.names=1)
    return(toString(row.names(geneUP)[which(row.names(geneUP)  %in% IPR$IDs[IPR$ANN %in% el])]))
}


enrichmentGO<-function(){
    exoDerGo<-read.table(opt$background,header=T);
    goFrame<-GOFrame(exoDerGo,organism="Exophiala dermatitidis")
    goAllFrame=GOAllFrame(goFrame)
    gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
                                        #generate universe id
    universe<-getGOFrameData(goAllFrame)
    universe<-unique(universe$gene_id)
    diffData<-read.table(opt$geneFile,header=T,row.names=1)
    geneList<-row.names(diffData)
    geneList<-intersect(universe,geneList)
    params<-GSEAGOHyperGParams(name="My Custom GSEA based annot Params", geneSetCollection=gsc, geneIds=geneList, universeGeneIds=universe, ontology=c("BP","CC","MF"),pvalueCutoff=1,conditional=F,testDirection="over")
    over<-hyperGTest(params)
    over<-summary(over)
    over$fdr<-p.adjust(over$Pvalue,method="fdr")
    bla<-over[over$fdr<as.numeric(opt$pValue),]
    write.table(bla,file=paste(opt$geneFile,opt$type,"csv",sep="."))
}

enrichmentElse<-function(){
    exoPFAM <- read.table(opt$background,head=T,stringsAsFactors=FALSE, row.names=NULL)
    sets<-split(exoPFAM$IDs,exoPFAM$ANN)
    gsc <- GeneSetCollection(Map(function(pid, gids) {
        GeneSet(gids, setName=pid, collectionType=PfamCollection(pid))
    }, names(sets), sets))
    universe<-unique(exoPFAM$IDs)
    diffData<-read.table(opt$geneFile,header=T,row.names=1)
    geneList<-row.names(diffData)
    geneList<-intersect(universe,geneList)
    params<-GSEAKEGGHyperGParams(name="Test",geneSetCollection=gsc,geneIds=geneList,universeGeneIds=universe,testDirection="over",pvalueCutoff=1)
    over<-hyperGTest(params)
    over<-summary(over)
    #write.table(ps);
    over$fdr<-p.adjust(over$Pvalue,method="fdr")
    bla<-over[over$fdr<as.numeric(opt$pValue),]
    if(opt$dictionary!="NA"){
        desc<-fread(opt$dictionary,header=F)
        bla$desc<-desc$V2[match(bla$KEGGID,desc$V1)]
    }
    bla$genes<-as.character(lapply(bla$KEGGID,getCorrespondingGene))
    write.table(bla,file=paste(opt$geneFile,opt$type,"csv",sep="."),sep=",")
}

if(opt$type=="GO"){
    library("AnnotationForge")
    library("GOstats")
    library("GSEABase")
    library("xtable")
    enrichmentGO()
}else{
    library("data.table")
    library("AnnotationForge")
    library("GOstats")
    library("GSEABase")
    enrichmentElse()
}

