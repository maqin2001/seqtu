#!/usr/bin/Rscript
args <- commandArgs(TRUE)

#set path
#pathForwardPredictResult <- "/home/wenchichou/project/TU/latestResults/final.c1f"
pathForwardPredictResult <- args[1] 
#pathReversePredictResult <- "/home/wenchichou/project/TU/latestResults/final.c1r"
pathReversePredictResult <- args[2] 
#outFileName <- "dataset1"
outFileName <- args[3]
#read ppt 
ptt=read.table("/home/wenchichou/project/TU/NC_009012.ptt", sep = "\t", skip = 2, header=T, check.names =F,  quote = "")
#dim(ptt)

#read gff
x=read.table("/home/wenchichou/project/TU/NC_009012.gff", sep = "\t", skip = 5)
x1<-x[x$V3=="gene",]
tag=unlist(strsplit(as.character(x1[1,]$V9),";"))[grep("locus_tag=",unlist(strsplit(as.character(x1[1,]$V9),";")))]
tag=unlist(strsplit(tag,"="))[2]
x2<-within(x1[1,], y<-data.frame(do.call('rbind',strsplit(as.character(x1[1,]$V9), ';', fixed=TRUE))))
x3<-within(x2, y<-data.frame(do.call('rbind',strsplit(as.character(x2$y$X2), '=', fixed=TRUE))))
Ngene = data.frame(tag, x3$V4, x3$V5, x3$V7)
for (i in 2:nrow(x1)){
    tag=unlist(strsplit(as.character(x1[i,]$V9),";"))[grep("locus_tag=",unlist(strsplit(as.character(x1[i,]$V9),";")))]
    tag=unlist(strsplit(tag,"="))[2]
    x2<-within(x1[i,], y<-data.frame(do.call('rbind',strsplit(as.character(x1[i,]$V9), ';', fixed=TRUE))))
    x3<-within(x2, y<-data.frame(do.call('rbind',strsplit(as.character(x2$y$X2), '=', fixed=TRUE))))
    x4<-data.frame(tag, x3$V4, x3$V5, x3$V7)
    Ngene = rbind(Ngene,x4)
}
colnames(Ngene)=c("ID", "Start", "End", "Strand") # Set Column Names
Ngene=Ngene[order(Ngene[,2]),] #Order the list by start positions of each gene
#Ngene[92,]
AllForwardGenes=Ngene[Ngene[,4]=="+",] #Select all genes on forward strand
AllReverseGenes=Ngene[Ngene[,4]=="-",] #Select all genes on reverse strand



#Read TU result tables
ForwardPredictResult<-read.table(pathForwardPredictResult, header=F)
ReversePredictResult<-read.table(pathReversePredictResult, header=F)

dim(ForwardPredictResult)
dim(ReversePredictResult)

dim(AllForwardGenes)
dim(AllReverseGenes)
#head(ForwardPredictResult)
#head(ReversePredictResult)

#Generate Target TU SVM format
GenerateOneTargetTU<-function(InputGeneList, i){
    OneTUList=data.frame()
    OneTUList=c("5'Gene", InputGeneList[i:(i+1),][1,2], InputGeneList[i:(i+1),][1,3], as.character(InputGeneList[i:(i+1),][1,4]));
    OneTUList=rbind(OneTUList, c("IntergenicRegion", (InputGeneList[i:(i+1),][1,3]+1), (InputGeneList[i:(i+1),][2,2]-1),as.character( InputGeneList[i:(i+1),][1,4])));
    OneTUList=rbind(OneTUList, c("3'Gene", as.numeric(as.character(InputGeneList[i:(i+1),][2,2])), InputGeneList[i:(i+1),][2,3], as.character(InputGeneList[i:(i+1),][2,4])));
    OneTUList=rbind(OneTUList, c("FullTU", InputGeneList[i:(i+1),][1,2], InputGeneList[i:(i+1),][2,3], as.character(InputGeneList[i:(i+1),][1,4])));
    row.names(OneTUList)=NULL
    OneTUList=as.data.frame(OneTUList)
    OneTUList[,2]=as.numeric(as.character(OneTUList[,2]))
    OneTUList[,3]=as.numeric(as.character(OneTUList[,3]))
    colnames(OneTUList)=c("ID", "Start", "End", "Strand") # Set Column Names
    return(OneTUList)
}


GenerateFinalTUTable<-function(TargetTUsPrediction, SelectedGenes, PTT){
#TargetTUsPrediction=ForwardPredictResult
#SelectedGenes=AllForwardGenes
#PTT=ptt
    colnames(TargetTUsPrediction)[1] = "TUresult"
    TargetTUName1Name2StartEnd=data.frame()
    j=1
    for(i in 1:(nrow(SelectedGenes)-1)){
        TargetTUName1Name2StartEnd=rbind(TargetTUName1Name2StartEnd, cbind(SelectedGenes[i,1], SelectedGenes[(i+1),1],GenerateOneTargetTU(SelectedGenes,i)[4,2:3]))
        j=j+1
    }
    TargetTUsPrediction = (cbind(TargetTUsPrediction, TargetTUName1Name2StartEnd))
#	head(TargetTUsPrediction)

	colnames(TargetTUsPrediction)[2:3]=c("GeneName1", "GeneName2")
	GeneProduct1=NULL
	GeneProduct1=as.character(PTT$Product[match(TargetTUsPrediction$GeneName1, PTT$Synonym)])
	GeneProduct1[which(is.na(GeneProduct1))]=as.character(TargetTUsPrediction$GeneName1[which(is.na(GeneProduct1))])
#	GeneProduct1[1669:1675]
	GeneProduct2=NULL
	GeneProduct2=as.character(PTT$Product[match(TargetTUsPrediction$GeneName2, PTT$Synonym)])
	GeneProduct2[which(is.na(GeneProduct2))]=as.character(TargetTUsPrediction$GeneName2[which(is.na(GeneProduct2))])
	TargetTUsPrediction=data.frame(TargetTUsPrediction, GeneProduct1=GeneProduct1, GeneProduct2=GeneProduct2)

#	head(TargetTUsPrediction)

    FinalTUAll=data.frame()
    FinalTUOne=data.frame()
    FinalTUOneNameStart=as.character(TargetTUsPrediction[1,"GeneName1"])
    FinalTUOneProductStart=as.character(TargetTUsPrediction[1,"GeneProduct1"])
    FinalTUOneNameEnd=as.character("")
    FinalTUOneProductEnd=as.character("")
    FinalTUOnePosStart=TargetTUsPrediction[1,4]
    FinalTUOnePosEnd=TargetTUsPrediction[1,5]
    FinalTUOne=data.frame(FinalTUOneNameStart,FinalTUOneNameEnd,FinalTUOnePosStart,FinalTUOnePosEnd,FinalTUOneProductStart,FinalTUOneProductEnd)
    flag=TargetTUsPrediction[1,1]
    for(i in c(1:nrow(TargetTUsPrediction))){
        if(TargetTUsPrediction[i,1]==1){  #TU result = 1
            FinalTUOne[,2]=paste(FinalTUOne[,2],TargetTUsPrediction[i,"GeneName2"], sep=";")
            FinalTUOne[,6]=paste(FinalTUOne[,6],TargetTUsPrediction[i,"GeneProduct2"], sep=";")
            FinalTUOne[,4]=TargetTUsPrediction[i,5]
            flag=1
        }else{                                   #TU result = -1
            FinalTUAll=rbind(FinalTUAll,FinalTUOne)
            FinalTUOne=NULL
            flag=-1
            FinalTUOneNameStart=as.character(TargetTUsPrediction[i,"GeneName2"])
            FinalTUOneProductStart=as.character(TargetTUsPrediction[i,"GeneProduct2"])
            FinalTUOneNameEnd=as.character("")
    		FinalTUOneProductEnd=as.character("")
            FinalTUOnePosStart=TargetTUsPrediction[i+1,4]
            FinalTUOnePosEnd=TargetTUsPrediction[i,5]
    		FinalTUOne=data.frame(FinalTUOneNameStart,FinalTUOneNameEnd,FinalTUOnePosStart,FinalTUOnePosEnd,FinalTUOneProductStart,FinalTUOneProductEnd)
        }
    }#for
    FinalTUAll=rbind(FinalTUAll,FinalTUOne)
	FinalTUAll$allNames=paste(FinalTUAll[,"FinalTUOneNameStart"],FinalTUAll[,"FinalTUOneNameEnd"],sep="")
	FinalTUAll$allProducts=paste(FinalTUAll[,"FinalTUOneProductStart"],FinalTUAll[,"FinalTUOneProductEnd"],sep="")
	FinalTUAll<-FinalTUAll[,c("FinalTUOnePosStart","FinalTUOnePosEnd","allNames","allProducts")]
    return(FinalTUAll)
} #end of function

ForwardTUtable<-cbind(GenerateFinalTUTable(ForwardPredictResult,AllForwardGenes,ptt),Strand="+")
ReverseTUtable<-cbind(GenerateFinalTUTable(ReversePredictResult,AllReverseGenes,ptt),Strand="-")
#write.csv(ForwardTUtable, file = "ForwardTUtable.csv", quote = FALSE)
#write.csv(ReverseTUtable, file = "ReverseTUtable.csv", quote = FALSE)

AllTUtable<-rbind(ForwardTUtable,ReverseTUtable)
AllTUtable <- AllTUtable[,c("FinalTUOnePosStart","FinalTUOnePosEnd","Strand","allNames","allProducts")]
AllTUtable <- AllTUtable[order(AllTUtable$FinalTUOnePosStart),]
row.names(AllTUtable) <- paste("TU",c(1:nrow(AllTUtable)), sep="")

write.csv(AllTUtable, file = paste(outFileName,"-AllTUtable.csv",sep=""),quote = TRUE)


