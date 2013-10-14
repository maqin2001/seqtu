#!/usr/bin/Rscript
#Identify Transcription Units using strand-specific RNA-seq data
#Authors: 
#	1. Wen-Chi Chou wenchichou@hsl.harvard.edu 
#	2. Qin Ma maqin@uga.edu

#prerequisite:
#libsvm
#	The current release (Version 3.17, April 2013) of LIBSVM can be obtained by downloading the zip file or tar.gz file. You can also check this github directory. Please e-mail us if you have problems to download the file. 
#	http://www.csie.ntu.edu.tw/~cjlin/cgi-bin/libsvm.cgi?+http://www.csie.ntu.edu.tw/~cjlin/libsvm+tar.gz
#	http://www.csie.ntu.edu.tw/~cjlin/libsvm/


#load necessary libraries
library(grid)
library(gridBase)
library(ggplot2)
library(seqinr)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

#set paths
file_mappedStrandSpecificRNAseqSignals <- "/home/wenchichou/project/TU/data/s_1_1_qfilter_20q_95p.fastq.bwa.sam.multi.NA"
file_gff <- "/home/wenchichou/project/TU/NC_009012.gff"
file_genomeSequence <- "/home/wenchichou/project/TU/NC_009012.fna"
FileNamePrefix <- "testChou"
dir_libsvm <- "/home/wenchichou/project/TU/bin/libsvm-3.17/"

#args <- commandArgs(TRUE)
#Clean SVM files
system("ls -ltr")        #Check ls
system("rm SVM.*")       #Remove SVM files
system("rm TargetSVM.*") #Remove target SVM files
system("ls -ltr")        #Check ls

# Remove all variables before running R script 
rm(list=ls(all=TRUE))

# open a statistic file
#sink("stat.txt")

# set some global variables
MeanQuantile = 0.1
MedianValue = 1
FilterMatrixQuantile=0.3

#Read EBCDs of each RNAseq data set (first column:Forward strand; second column:Reverse strand)
SelectedRNAseqData<-read.table(file_mappedStrandSpecificRNAseqSignals, head=F)


#Read Gene Annotation 
x=read.table(file_gff, sep = "\t", skip = 5)
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

#Select all genes on forward strand
AllForwardGenes=Ngene[Ngene[,4]=="+",] 
#Select all genes on reverse strand
AllReverseGenes=Ngene[Ngene[,4]=="-",] 
#dim(AllForwardGenes)
#dim(AllReverseGenes)

#Read bacterial Genome Sequence
seq=read.fasta(file_genomeSequence,seqonly=TRUE,as.string=TRUE)
seq=unlist(strsplit(seq[[1]], split="")) #Extract string
#length(seq) #Check genome length

#Function: Generate intergenic regions according to a given gene list
GenerateIGList<-function(InputGeneList){
	IGList=data.frame()
	for(i in 1:(nrow(InputGeneList)-1)){
	    Name=paste(InputGeneList[i:(i+1),][1,1], InputGeneList[i:(i+1),][2,1])
	    Start=InputGeneList[i:(i+1),][1,3]+1
	    End=InputGeneList[i:(i+1),][2,2]-1
	    Strand=InputGeneList[i:(i+1),][1,4]
	    IGList[i,1]=Name
	    IGList[i,2]=Start
	    IGList[i,3]=End
	    IGList[i,4]=Strand
	}
	colnames(IGList)=c("ID", "Start", "End", "Strand") # Set Column Names
	return(IGList)
} #END Function: Generate intergenic regions according to a given gene list
AllForwardIntergenicRegions=data.frame()
AllForwardIntergenicRegions=GenerateIGList(AllForwardGenes) #All intergenic regions on forward strand
AllReverseIntergenicRegions=data.frame()
AllReverseIntergenicRegions=GenerateIGList(AllReverseGenes) #All intergenic regions on reverse strand 

#Function: Get all expression values of a region
GetRegionExp<-function(SelectedRNAseqData, InputGeneList, i){
        if(InputGeneList[i,4] == "+"){
                SubExp=SelectedRNAseqData[InputGeneList[i,2]:InputGeneList[i,3],1]
        }else{
                SubExp=SelectedRNAseqData[InputGeneList[i,2]:InputGeneList[i,3],2]
        }
    return(SubExp)
}#END Function: Get all expression values of a region

#Function: Get all gap values of a region
GetRegionGap<-function(SelectedRNAseqData, InputGeneList, i){
        if(InputGeneList[i,4] == "+"){
                SubGap=(SelectedRNAseqData[InputGeneList[i,2]:InputGeneList[i,3],1]==0)*1
        }else{
                SubGap=(SelectedRNAseqData[InputGeneList[i,2]:InputGeneList[i,3],2]==0)*1
        }
    return(SubGap)
}#END Function: Get all gap values of a region



#Get mean expression level, mean gaps and lenght for each gene on forward strand;
AllForwardGenesExpressionMean=NULL;
AllForwardGenesExpressionMedian=NULL;
AllForwardGenesGapsMean=NULL;
AllForwardGenesLength=NULL;
for(i in 1:nrow(AllForwardGenes)){;
	AllForwardGenesExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllForwardGenes, i));
	AllForwardGenesExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllForwardGenes, i));
    AllForwardGenesGapsMean[i]=mean(GetRegionGap(SelectedRNAseqData,AllForwardGenes, i));
    AllForwardGenesLength[i]=AllForwardGenes[i,3]-AllForwardGenes[i,2]+1;
};
#Get mean expression level, mean gaps and lenght for each gene on reverse strand;
AllReverseGenesExpressionMean=NULL;
AllReverseGenesExpressionMedian=NULL;
AllReverseGenesGapsMean=NULL;
AllReverseGenesLength=NULL;
for(i in 1:nrow(AllReverseGenes)){;
	AllReverseGenesExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllReverseGenes, i));
	AllReverseGenesExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllReverseGenes, i));
    AllReverseGenesGapsMean[i]=mean(GetRegionGap(SelectedRNAseqData,AllReverseGenes, i));
    AllReverseGenesLength[i]=AllReverseGenes[i,3]-AllReverseGenes[i,2]+1;
};
;
#Get mean expression level, mean gaps and lenght for each intergenic region on forward strand;
AllForwardIntergenicRegionsExpressionMean=NULL;
AllForwardIntergenicRegionsGapsMean=NULL;
AllForwardIntergenicRegionsLength=NULL;
for(i in 1:(nrow(AllForwardIntergenicRegions))){;
	AllForwardIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
    AllForwardIntergenicRegionsGapsMean[i]=mean(GetRegionGap(SelectedRNAseqData,AllForwardIntergenicRegions, i));
    AllForwardIntergenicRegionsLength[i]=AllForwardIntergenicRegions[i,3]-AllForwardIntergenicRegions[i,2]+1;
};
#Get mean expression level, mean gaps and lenght for each intergenic region on reverse strand;
AllReverseIntergenicRegionsExpressionMean=NULL;
AllReverseIntergenicRegionsGapsMean=NULL;
AllReverseIntergenicRegionsLength=NULL;
for(i in 1:(nrow(AllReverseIntergenicRegions))){;
	AllReverseIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
    AllReverseIntergenicRegionsGapsMean[i]=mean(GetRegionGap(SelectedRNAseqData,AllReverseIntergenicRegions, i));
    AllReverseIntergenicRegionsLength[i]=AllReverseIntergenicRegions[i,3]-AllReverseIntergenicRegions[i,2]+1;
};

#plot statistics of ssRNA-seq signals
EG1=EG2=EG12=p12=NULL;
EG1 <- data.frame(E=log2(AllForwardGenesExpressionMean+1), G=AllForwardGenesGapsMean);
EG2 <- data.frame(E=log2(AllReverseGenesExpressionMean+1), G=AllReverseGenesGapsMean);
EG12<-rbind(EG1, EG2);
#cat(summary(EG12),"\n");
p12<-ggplot(EG12,aes(x=EG12$E,y=EG12$G)) +
geom_point(alpha = 0.3) +
#geom_smooth() +;
xlab("log2(Expression)") +
ylab("Gap %") 

EG3=EG4=EG34=p34=NULL;
EG3 <- data.frame(E=log2(AllForwardIntergenicRegionsExpressionMean+1), G=AllForwardIntergenicRegionsGapsMean);
EG4 <- data.frame(E=log2(AllReverseIntergenicRegionsExpressionMean+1), G=AllReverseIntergenicRegionsGapsMean);
EG34<-rbind(EG3, EG4);
p34<-ggplot(EG34,aes(x=EG34$E,y=EG34$G)) +
geom_point(alpha = 0.3) +
#geom_smooth() +
scale_colour_manual(values="blue") +
xlab("log2(Expression)") +
ylab("Gap %")

EG1234<-rbind(cbind(EG12,type="Coding region"), cbind(EG34,type="Intergenic region"))
EG1234$G<-EG1234$G*100
p1234<-ggplot(EG1234,aes(x=E,y=G, group=type,col=type, shape=type)) +
geom_point(alpha = 0.5) +
theme(legend.position=c(1.2,1.2),legend.title = element_blank())+
scale_colour_manual(values=c("blue", "red")) +
xlab("log2(Expression)") +
ylab("Gap %")
d1<-ggplot(EG1234, aes(x=E,y=..count.., fill=type))+geom_density(alpha = .5,size=0.01)+theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(), plot.margin = unit(c(0,0,-2,2.5),"lines")
)+scale_fill_manual(values=c("blue", "red"))

d2<-ggplot(EG1234, aes(x=G,y=..count.., fill=type))+geom_density(alpha=.5, size=0.01)+coord_flip()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),plot.margin = unit(c(0,0,1.7,-2),"lines")
)+scale_fill_manual(values=c("blue", "red"))

pdf(file=paste("BasicStatistics-",FileNamePrefix,"-",format(Sys.time(),"%m%d%y-%k%M"),".pdf",sep=""), onefile=TRUE);
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 5))) # a 5 by 5 grid
print(p1234, vp=vplayout(2:5,1:4)) # the main x/y plot will instead spread across most of the grid
print(d1, vp=vplayout(1,1:4)) # the first density plot will occupy the top of the grid
print(d2, vp=vplayout(2:5,5)) # with the second density plot occupying a narrow vertical strip at the right
dev.off();












getIntergenicInfo<-function(SelectedRNAseqData){
	AllForwardIntergenicRegionsExpressionMean=NULL;
	AllForwardIntergenicRegionsExpressionMedian=NULL;
	AllForwardIntergenicRegionsGapsMean=NULL;
	AllForwardIntergenicRegionsLength=NULL;
	for(i in 1:(nrow(AllForwardIntergenicRegions))){
		AllForwardIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
		AllForwardIntergenicRegionsExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
	    AllForwardIntergenicRegionsGapsMean[i]=mean(GetRegionGap(SelectedRNAseqData,AllForwardIntergenicRegions, i));
	    AllForwardIntergenicRegionsLength[i]=AllForwardIntergenicRegions[i,3]-AllForwardIntergenicRegions[i,2]+1;
	}
	#Get mean expression level, mean gaps and lenght for each intergenic region on reverse strand
	AllReverseIntergenicRegionsExpressionMean=NULL;
	AllReverseIntergenicRegionsExpressionMedian=NULL;
	AllReverseIntergenicRegionsGapsMean=NULL;
	AllReverseIntergenicRegionsLength=NULL;
	for(i in 1:(nrow(AllReverseIntergenicRegions))){
		AllReverseIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
		AllReverseIntergenicRegionsExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
	    AllReverseIntergenicRegionsGapsMean[i]=mean(GetRegionGap(SelectedRNAseqData,AllReverseIntergenicRegions, i));
	    AllReverseIntergenicRegionsLength[i]=AllReverseIntergenicRegions[i,3]-AllReverseIntergenicRegions[i,2]+1;
	}
	MyList<- list("a"=AllForwardIntergenicRegionsExpressionMean, "b"=AllReverseIntergenicRegionsExpressionMean, "c"=AllForwardIntergenicRegionsExpressionMedian, "d"=AllReverseIntergenicRegionsExpressionMedian);
	return(MyList); 
}#end function getIntergenicInfo<-function(SelectedRNAseqData){


#Get length proportions of all forward intergenic regions
ForwardIntergenicRegionProportion=NULL
ForwardIntergenicRegionDeviate=NULL
for(i in 1:(nrow(AllForwardGenes)-1)){
    FullLength=AllForwardGenes[i:(i+1),][2,3]-AllForwardGenes[i:(i+1),][1,2]+1
    GeneLengthDifferent = abs(AllForwardGenes[i:(i+1),][2,3]+AllForwardGenes[i:(i+1),][1,2]-AllForwardGenes[i:(i+1),][2,2]-AllForwardGenes[i:(i+1),][1,3])
    ForwardIntergenicRegionProportion[i]=(AllForwardGenes[i:(i+1),][2,2]-AllForwardGenes[i:(i+1),][1,3]-1)/FullLength
    ForwardIntergenicRegionDeviate[i]=GeneLengthDifferent/FullLength
}
cat ("11.Proportion of intergenic region's length and sum length of flanking two genes on forward strand:",summary(ForwardIntergenicRegionProportion), "\n", sep = "\t")
#summary(ForwardIntergenicRegionProportion)
#summary(ForwardIntergenicRegionDeviate)

#Get length proportions of all reverse intergenic regions
ReverseIntergenicRegionProportion=NULL
ReverseIntergenicRegionDeviate=NULL
for(i in 1:(nrow(AllReverseGenes)-1)){
    FullLength=AllReverseGenes[i:(i+1),][2,3]-AllReverseGenes[i:(i+1),][1,2]+1
    GeneLengthDifferent = abs(AllReverseGenes[i:(i+1),][2,3]+AllReverseGenes[i:(i+1),][1,2]-AllReverseGenes[i:(i+1),][2,2]-AllReverseGenes[i:(i+1),][1,3])
    ReverseIntergenicRegionProportion[i]=(AllReverseGenes[i:(i+1),][2,2]-AllReverseGenes[i:(i+1),][1,3]-1)/FullLength
    ReverseIntergenicRegionDeviate[i]=GeneLengthDifferent/FullLength
}
cat ("12.Proportion of intergenic region's length and sum length of flanking two genes on reverse strand:",summary(ReverseIntergenicRegionProportion), "\n", sep= "\t")
#summary(ReverseIntergenicRegionProportion)
#summary(ReverseIntergenicRegionDeviate)

#generate the ForwardIntergenicRegionProportionNew
ForwardIntergenicRegionProportionNew=NULL
num=1
for (i in 1:(nrow(AllForwardGenes)-1)){ #Why "-1" ?
    if(AllForwardIntergenicRegionsExpressionMean[i]>= quantile(AllForwardIntergenicRegionsExpressionMean, 0.60)){ #capture expressed intergenic regions
		if (ForwardIntergenicRegionProportion[i]>0){
		    ForwardIntergenicRegionProportionNew[num] = ForwardIntergenicRegionProportion[i]
	    	num=num+1    
		}
    }
}
#summary(ForwardIntergenicRegionProportionNew)
cat("13.",length(ForwardIntergenicRegionProportionNew), "forward intergenic regions are kept to simulate TUs, whose expression value is larget than 60 quantile","\n",sep=" ")

#generate the ReverseIntergenicRegionProportionNew
ReverseIntergenicRegionProportionNew=NULL
num=1
for (i in 1:(nrow(AllReverseGenes)-1)){ #Why "-1" ?
    if(AllReverseIntergenicRegionsExpressionMean[i]>= quantile(AllReverseIntergenicRegionsExpressionMean, 0.60)){ #capture expressed intergenic regions
		if (ReverseIntergenicRegionProportion[i]>0){
		    ReverseIntergenicRegionProportionNew[num] = ReverseIntergenicRegionProportion[i]
		    num=num+1    
		}
    }
}
#summary(ReverseIntergenicRegionProportionNew)
cat("14.",length(ReverseIntergenicRegionProportionNew), "reverse intergenic regions are kept to simulate TUs, whose expression value is larget than 60 quantile","\n",sep=" ")


#Function: Generate one simulated TU
GenerateOneSimulatedTU<-function(InputGeneList, IntergenicRegionProportion, IntergenicRegionDeviate, i){
    SimIGLength=as.integer((InputGeneList[i,3]-InputGeneList[i,2]+1)*sample(IntergenicRegionProportion,1)) #Need a real probability
    if(SimIGLength %% 2 == 1){SimIGLength=SimIGLength+1}
    MiddlePoint=as.integer(median(InputGeneList[i,2]:InputGeneList[i,3]))
    #get the maximal AT-rich region in the perturbation
    GCvalue = 1
    for (j in 1:as.integer((InputGeneList[i,3]-InputGeneList[i,2]+1)*median(IntergenicRegionDeviate)/2)){
        SubSeq=seq[(MiddlePoint-(0.5*SimIGLength)+1-j):(SimIGLength+MiddlePoint-(0.5*SimIGLength)-j)]
        if( GC(SubSeq)<=GCvalue & length(SubSeq)>1){
            GCvalue = GC(SubSeq)
            MiddlePointNew=MiddlePoint-j
        }
    }
    for (j in 1:as.integer((InputGeneList[i,3]-InputGeneList[i,2]+1)*median(IntergenicRegionDeviate)/2)){
        SubSeq=seq[(MiddlePoint-(0.5*SimIGLength)+1+j):(SimIGLength+MiddlePoint-(0.5*SimIGLength)+j)]
        if( GC(SubSeq)<=GCvalue & length(SubSeq)>1){
            GCvalue = GC(SubSeq)
            MiddlePointNew=MiddlePoint+j
        }
    }
    MiddlePoint = MiddlePointNew
    LeftSubGeneStart=InputGeneList[i,2]
    LeftSubGeneEnd=MiddlePoint-(0.5*SimIGLength)
    IGSubGeneStart=MiddlePoint-(0.5*SimIGLength)+1
    IGSubGeneEnd=SimIGLength+MiddlePoint-(0.5*SimIGLength)
    RightSubGeneStart=SimIGLength+MiddlePoint-(0.5*SimIGLength)+1
    RightSubGeneEnd=InputGeneList[i,3]
    SubGeneList=NULL
    SubGeneList=c("LeftSubGene", LeftSubGeneStart, LeftSubGeneEnd, as.character(InputGeneList[i,4]));
    SubGeneList=rbind(SubGeneList, c("IGSubGene", IGSubGeneStart, IGSubGeneEnd, as.character(InputGeneList[i,4])));
    SubGeneList=rbind(SubGeneList, c("RightSubGene", RightSubGeneStart, RightSubGeneEnd, as.character(InputGeneList[i,4])));
    SubGeneList=rbind(SubGeneList, c("FullSubGene", InputGeneList[i,2], InputGeneList[i,3], as.character(InputGeneList[i,4])));
    row.names(SubGeneList)=NULL
    SubGeneList=as.data.frame(SubGeneList)
    SubGeneList[,2]=as.numeric(as.character(SubGeneList[,2]))
    SubGeneList[,3]=as.numeric(as.character(SubGeneList[,3]))
    return(SubGeneList)
}#END Function: Generate one simulated TU

#Function: Get gap percentage of a Intergenic region over a full TU
GetGapPercentage<-function(NAExp, InputGeneList, i){
	if(InputGeneList[i,4] == "+"){
		if(sum((NAExp[InputGeneList[4,2]:InputGeneList[4,3],1]<1)*1)==0){
			SubGaps=0
		}else{
			SubGaps=sum((NAExp[InputGeneList[i,2]:InputGeneList[i,3],1]<1)*1)/sum((NAExp[InputGeneList[4,2]:InputGeneList[4,3],1]<1)*1)
		}
	}else{
		if(sum((NAExp[InputGeneList[4,2]:InputGeneList[4,3],2]<1)*1)==0){
			SubGaps=0
		}else{
			SubGaps=sum((NAExp[InputGeneList[i,2]:InputGeneList[i,3],2]<1)*1)/sum((NAExp[InputGeneList[4,2]:InputGeneList[4,3],2]<1)*1)
		}
	}
	return(SubGaps)
}#END Function: Get gap percentage of a Intergenic region over a full TU


#Function: Calculate foldchanges
FoldChange<-function(num1, num2){
    MAX=NULL
    MAX=max(num1,num2)
    MIN=NULL
    if(min(num1,num2) <= 0){
        MIN=0.00001
    }else{
        MIN=min(num1, num2)
    }
    Change=MAX/MIN
    return(Change)
}#END Function: Calculate foldchanges

#Function: Get New Start and End positions for removing the position bias
NewStartEnd<-function(Start, End){
    Range=c(Start:End)
    NewStartInFactor=round(length(Range)/100*5,0)
    NewEndInFactor=length(Range)-NewStartInFactor
    NewStart=Range[NewStartInFactor]
    NewEnd=Range[NewEndInFactor]
    Positions=c(NewStart, NewEnd)
    return(Positions)
}#END Function: Get New Start and End positions for removing the position bias

#Function: Get Expression of a region with its position bias removed
GetRegionExpRemovePositionBias<-function(NAExp, InputGeneList, i){
        if(InputGeneList[i,4] == "+"){
                NewRange=NewStartEnd(InputGeneList[i,2], InputGeneList[i,3])
                NewRange=c(NewRange[1]:NewRange[2])
                SubExp=NAExp[NewRange,1]
        }else{
                SubExp=NAExp[NewRange,2]
        }
    return(SubExp)
}#END Function: Get Expression of a region with its position bias removed

#101413##Get mean expression level, mean gaps and lenght for each gene on forward strand;
#101413#AllForwardGenesExpressionMean=NULL;
#101413#AllForwardGenesGapsMean=NULL;
#101413#AllForwardGenesLength=NULL;
#101413#for(i in 1:nrow(AllForwardGenes)){;
#101413#	AllForwardGenesExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllForwardGenes, i));
#101413#    AllForwardGenesGapsMean[i]=mean(GetRegionGap(SelectedRNAseqData,AllForwardGenes, i));
#101413#    AllForwardGenesLength[i]=AllForwardGenes[i,3]-AllForwardGenes[i,2]+1;
#101413#};
#101413##Get mean expression level, mean gaps and lenght for each gene on reverse strand;
#101413#AllReverseGenesExpressionMean=NULL;
#101413#AllReverseGenesGapsMean=NULL;
#101413#AllReverseGenesLength=NULL;
#101413#for(i in 1:nrow(AllReverseGenes)){;
#101413#	AllReverseGenesExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllReverseGenes, i));
#101413#    AllReverseGenesGapsMean[i]=mean(GetRegionGap(SelectedRNAseqData,AllReverseGenes, i));
#101413#    AllReverseGenesLength[i]=AllReverseGenes[i,3]-AllReverseGenes[i,2]+1;
#101413#};

#Generate forward simulated TU matrix
SimulatedForwardTUMatrix=data.frame()
SimulatedForwardTUIGExpressionMean=NULL
SimulatedForwardTU5endExpressionMean=NULL
SimulatedForwardTUIGGapsMean=NULL
SimulatedForwardTU5endGapsMean=NULL
SimulatedForwardTUIGLength=NULL
num=1
for(i in 1:(nrow(AllForwardGenes))){
		if ((AllForwardGenesExpressionMean[i] > quantile(AllForwardGenesExpressionMean, MeanQuantile)) & (AllForwardGenesExpressionMedian[i] >= MedianValue)){
        OneForwardSimulatedTU=GenerateOneSimulatedTU(AllForwardGenes, ForwardIntergenicRegionProportionNew, ForwardIntergenicRegionDeviate, i)
		SimulatedForwardTUMatrix[num,1]=GetGapPercentage(SelectedRNAseqData, OneForwardSimulatedTU,2) # 2 is the intergenic reiong of a simulated TU
        SimulatedForwardTUMatrix[num,2]=FoldChange(median(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,1)),median(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
		SimulatedForwardTUIGExpressionMean[num]=mean(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,2))
		SimulatedForwardTU5endExpressionMean[num]=mean(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,1))
		SimulatedForwardTUIGGapsMean[num]=mean(GetRegionGap(SelectedRNAseqData, OneForwardSimulatedTU,2))
		SimulatedForwardTU5endGapsMean[num]=mean(GetRegionGap(SelectedRNAseqData, OneForwardSimulatedTU,1))
		SimulatedForwardTUIGLength[num]=OneForwardSimulatedTU[2,3]-OneForwardSimulatedTU[2,2]+1
		num = num +1
    }
}
#Generate reverse simulated TU matrix
SimulatedReverseTUMatrix=data.frame()
SimulatedReverseTUIGExpressionMean=NULL
SimulatedReverseTU5endExpressionMean=NULL
SimulatedReverseTUIGGapsMean=NULL
SimulatedReverseTU5endGapsMean=NULL
SimulatedReverseTUIGLength=NULL
num=1
for(i in 1:(nrow(AllReverseGenes))){
    if ((AllReverseGenesExpressionMean[i] > quantile(AllReverseGenesExpressionMean, MeanQuantile)) & (AllReverseGenesExpressionMedian[i] >= MedianValue)){
        OneReverseSimulatedTU=GenerateOneSimulatedTU(AllReverseGenes, ReverseIntergenicRegionProportionNew, ReverseIntergenicRegionDeviate, i)
        SimulatedReverseTUMatrix[num,1]=GetGapPercentage(SelectedRNAseqData, OneReverseSimulatedTU,2) # 2 is the intergenic reiong of a simulated TU
        SimulatedReverseTUMatrix[num,2]=FoldChange(median(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,1)),median(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
		SimulatedReverseTUIGExpressionMean[num]=mean(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,2))
		SimulatedReverseTU5endExpressionMean[num]=mean(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,1))
		SimulatedReverseTUIGGapsMean[num]=mean(GetRegionGap(SelectedRNAseqData, OneReverseSimulatedTU,2))
		SimulatedReverseTU5endGapsMean[num]=mean(GetRegionGap(SelectedRNAseqData, OneReverseSimulatedTU,1))
		SimulatedReverseTUIGLength[num]=OneReverseSimulatedTU[2,3]-OneReverseSimulatedTU[2,2]+1
	num = num + 1
    }
}

EG9=EG0=EG90=EG1=EG2=EG12=NULL;
EG9 <- data.frame(E=log2(SimulatedForwardTUIGExpressionMean+1), G=SimulatedForwardTUIGGapsMean);
EG0 <- data.frame(E=log2(SimulatedReverseTUIGExpressionMean+1), G=SimulatedReverseTUIGGapsMean);
EG90<-rbind(EG9, EG0);
#EG1 <- data.frame(E=log2(AllForwardGenesExpressionMean+1), G=AllForwardGenesGapsMean);
#EG2 <- data.frame(E=log2(AllReverseGenesExpressionMean+1), G=AllReverseGenesGapsMean);
EG1 <- data.frame(E=log2(SimulatedForwardTU5endExpressionMean+1), G=SimulatedForwardTU5endGapsMean);
EG2 <- data.frame(E=log2(SimulatedReverseTU5endExpressionMean+1), G=SimulatedReverseTU5endGapsMean);
EG12<-rbind(EG1, EG2);

EG1290<-rbind(cbind(EG12,type="Coding region"), cbind(EG90,type="Intergenic region"))
EG1290$G<-EG1290$G*100

p1290<-ggplot(EG1290,aes(x=E,y=G, group=type,col=type, shape=type)) +
geom_point(alpha = 0.5) +
theme(legend.position=c(1.2,1.2),legend.title = element_blank())+
scale_colour_manual(values=c("orange", "green")) +
scale_shape_manual(values=c(16,15)) +
xlab("log2(Expression)") +
ylab("Gap %")

d1<-ggplot(EG1290, aes(x=E,y=..density.., fill=type))+geom_density(alpha = .5,size=0.01)+theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(), plot.margin = unit(c(0,0,-2,2.5),"lines")
)+scale_fill_manual(values=c("orange", "green"))

d2<-ggplot(EG1290, aes(x=G,y=..density.., fill=type))+geom_density(alpha=.5, size=0.01)+coord_flip()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),plot.margin = unit(c(0,0,1.7,-2),"lines")
)+scale_fill_manual(values=c("orange", "green"))



pdf(file=paste("BasicStatistics-","Simulated-",FileNamePrefix,"-",format(Sys.time(),"%m%d%y-%k%M"),".pdf",sep=""), onefile=TRUE);
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 5))) # a 5 by 5 grid
print(p1290, vp=vplayout(2:5,1:4)) # the main x/y plot will instead spread across most of the grid
print(d1, vp=vplayout(1,1:4)) # the first density plot will occupy the top of the grid
print(d2, vp=vplayout(2:5,5)) # with the second density plot occupying a narrow vertical strip at the right
dev.off();


which(is.na(SimulatedForwardTUMatrix[,1])) #check NaN elements
which(is.na(SimulatedForwardTUMatrix[,2])) #check NaN elements
which(is.na(SimulatedReverseTUMatrix[,1])) #check NaN elements
which(is.na(SimulatedReverseTUMatrix[,2])) #check NaN elements
summary(SimulatedForwardTUMatrix)
summary(SimulatedReverseTUMatrix)
OriginalSimulatedForwardTUMatrix=SimulatedForwardTUMatrix
OriginalSimulatedReverseTUMatrix=SimulatedReverseTUMatrix
#SimulatedForwardTUMatrix=(SimulatedForwardTUMatrix[which(SimulatedForwardTUMatrix[,2]<2e+04),])
#SimulatedReverseTUMatrix=(SimulatedReverseTUMatrix[which(SimulatedReverseTUMatrix[,2]<2e+04),])

SimulatedForwardTUMatrix=OriginalSimulatedForwardTUMatrix[which(OriginalSimulatedForwardTUMatrix[,1] <= (quantile(OriginalSimulatedForwardTUMatrix[,1], FilterMatrixQuantile))),]
SimulatedReverseTUMatrix=OriginalSimulatedReverseTUMatrix[which(OriginalSimulatedReverseTUMatrix[,1] <= (quantile(OriginalSimulatedReverseTUMatrix[,1], FilterMatrixQuantile))),]
SimulatedForwardTUMatrix=OriginalSimulatedForwardTUMatrix[which(OriginalSimulatedForwardTUMatrix[,2] <= 10),]
SimulatedReverseTUMatrix=OriginalSimulatedReverseTUMatrix[which(OriginalSimulatedReverseTUMatrix[,2] <= 10),]

summary(SimulatedForwardTUMatrix)
summary(SimulatedReverseTUMatrix)

#Output SVM format
sink("SVM.Forward.Training")
for(i in 1:nrow(SimulatedForwardTUMatrix)){
    cat("+1\t1:", SimulatedForwardTUMatrix[i,1], "\t2:", SimulatedForwardTUMatrix[i,2], "\n", sep="")
}
sink()

sink("SVM.Reverse.Training")
for(i in 1:nrow(SimulatedReverseTUMatrix)){
    cat("+1\t1:", SimulatedReverseTUMatrix[i,1], "\t2:", SimulatedReverseTUMatrix[i,2], "\n", sep="")
}
sink()

#Run system() command in R script
#Check ls
#system("ls -ltr", wait = FALSE)
#system("ls -ltr SVM.*")
#Rename SVM format files
#system("for i in SVM.*; do j=\"_081112\"; mv $i $i$j; done")


#Scale SVM format files
system(paste("for i in SVM.*.Training; do ",dir_libsvm,"svm-scale -s $i.scaled.Ruler $i > $i.scaled; done", sep=""))
ForwardScaledFileName<-system("ls SVM.Forward.Training.scaled", intern=TRUE)
ForwardRulerFileName=paste(ForwardScaledFileName,".Ruler", sep="")
ForwardModelFileName=paste(ForwardScaledFileName,".Model", sep="")
ForwardPredictFileName=paste(ForwardScaledFileName,".Predict", sep="")

ReverseScaledFileName<-system("ls SVM.Reverse.Training.scaled", intern=TRUE)
ReverseRulerFileName=paste(ReverseScaledFileName,".Ruler", sep="")
ReverseModelFileName=paste(ReverseScaledFileName,".Model", sep="")
ReversePredictFileName=paste(ReverseScaledFileName,".Predict", sep="")

#Check scaled SVM format files
system(paste("for i in SVM.*.scaled; do ",dir_libsvm,"tools/checkdata.py $i; done", sep=""))
#Find the best parameters and results
GridResults<-system(paste("for i in SVM.*.scaled; do python ",dir_libsvm,"tools/grid_WCC081812.py -s 2 -t 2 -v 5 -m 300 $i; done",sep=""), intern=TRUE)
#Check Training results
#system("for i in *scaled.out; do echo $i; sort -n -k 3 $i|tail -1; done")

BestParameters=grep("local", GridResults, invert=TRUE, value=TRUE)
ForwardN=as.numeric(strsplit(BestParameters, " ")[[1]][1])
ForwardG=as.numeric(strsplit(BestParameters, " ")[[1]][2])
ReverseN=as.numeric(strsplit(BestParameters, " ")[[2]][1])
ReverseG=as.numeric(strsplit(BestParameters, " ")[[2]][2])


system(paste(dir_libsvm,"svm-train -s 2 -t 2 -n ",ForwardN," -g ",ForwardG," -n ",ForwardN," ",ForwardScaledFileName," ",ForwardModelFileName, sep=""))
system(paste(dir_libsvm,"svm-train -s 2 -t 2 -n ",ReverseN," -g ",ReverseG," ",ReverseScaledFileName," ",ReverseModelFileName, sep=""))
SimulatedForwardTUResult<-system(paste(dir_libsvm,"svm-predict ",ForwardScaledFileName," ",ForwardModelFileName," ",ForwardPredictFileName, sep=""), intern=TRUE)
SimulatedReverseTUResult<-system(paste(dir_libsvm,"svm-predict ",ReverseScaledFileName," ",ReverseModelFileName," ",ReversePredictFileName, sep=""), intern=TRUE)

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

TargetForwardTUMatrix=data.frame()
num=1
for(i in 1:(nrow(AllForwardGenes)-1)){
	OneForwardTargetTU=GenerateOneTargetTU(AllForwardGenes,i)
	TargetForwardTUMatrix[num,1]=GetGapPercentage(SelectedRNAseqData, OneForwardTargetTU,2) # 2 is the intergenic reiong of a target TU
    TargetForwardTUMatrix[num,2]=FoldChange(median(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,1)),median(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
    num = num +1
}
which(is.na(TargetForwardTUMatrix[,1])) #check NaN elements
which(is.na(TargetForwardTUMatrix[,2])) #check NaN elements


TargetReverseTUMatrix=data.frame()
num=1
for(i in 1:(nrow(AllReverseGenes)-1)){
	OneReverseTargetTU=GenerateOneTargetTU(AllReverseGenes,i)
	TargetReverseTUMatrix[num,1]=GetGapPercentage(SelectedRNAseqData, OneReverseTargetTU,2) # 2 is the intergenic reiong of a target TU
    TargetReverseTUMatrix[num,2]=FoldChange(median(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,1)),median(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
    num = num +1
}
which(is.na(TargetReverseTUMatrix[,1])) #check NaN elements
which(is.na(TargetReverseTUMatrix[,2])) #check NaN elements

#Output TargetSVM format
sink("TargetSVM.Forward")
for(i in 1:nrow(TargetForwardTUMatrix)){
    cat("+1\t1:", TargetForwardTUMatrix[i,1], "\t2:", TargetForwardTUMatrix[i,2], "\n", sep="")
}
sink()

sink("TargetSVM.Reverse")
for(i in 1:nrow(TargetReverseTUMatrix)){
    cat("+1\t1:", TargetReverseTUMatrix[i,1], "\t2:", TargetReverseTUMatrix[i,2], "\n", sep="")
}
sink()

#Rename SVM format files
#system("for i in TargetSVM.*; do j=\"_081112\"; mv $i $i$j; done")
#Scale SVM format files
system(paste(dir_libsvm,"svm-scale -r ",ForwardRulerFileName," TargetSVM.Forward > TargetSVM.Forward.scaled", sep=""))
TargetForwardScaledFileName<-system("ls TargetSVM.Forward.scaled", intern=TRUE)
TargetForwardPredictFileName=paste(TargetForwardScaledFileName,".Predict", sep="")

system(paste(dir_libsvm,"svm-scale -r ",ReverseRulerFileName," TargetSVM.Reverse > TargetSVM.Reverse.scaled", sep=""))
TargetReverseScaledFileName<-system("ls TargetSVM.Reverse.scaled", intern=TRUE)
TargetReversePredictFileName=paste(TargetReverseScaledFileName,".Predict", sep="")



TargetForwardTUResult<-system(paste(dir_libsvm,"svm-predict ",TargetForwardScaledFileName," ",ForwardModelFileName," ",TargetForwardPredictFileName, sep=""), intern=TRUE)
TargetReverseTUResult<-system(paste(dir_libsvm,"svm-predict ",TargetReverseScaledFileName," ",ReverseModelFileName," ",TargetReversePredictFileName, sep=""), intern=TRUE)

cat ("SVM Accuracy of simulated forward TUs:",SimulatedForwardTUResult,"\n",sep= "\t")
cat ("SVM Accuracy of simulated reverse TUs:",SimulatedReverseTUResult,"\n",sep= "\t")





