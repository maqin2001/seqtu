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

#041514 version updates:
#1. merge SVM training from two models (positive and negative strands) to one models (merging two strands to one strand before training)
#2. add negative data which are intergenic regions with over 10-fold change of their own two flanking genes.
#3. add 454 data as positive dataset

#020914 version updates: 
#adding
# 1. pass arguments
# 2. checking 454 results and output accuracy
# 3. add lamda as 0 to include all low expressed genes

rm(list=ls(all=TRUE))
source("/home/wc102/project/tu/bin/seqTU_functions_041114.r")
library(seqinr)
#load necessary libraries
args <- commandArgs(TRUE)
#library(grid)
#library(gridBase)
#library(ggplot2)
#vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
# Remove all variables before running R script 

#set paths
# /hms/scratch1/ifarmusc/wenchichou/data/s_7_1_qfilter_20q_95p.fastq.bwa.sam.multi.NA sampling.treatment3. /home/wc102/project/tu/data/header.treatment1.2.3.454.coveredreads.length.changeRatio

#file_RNAseqSignals <- "/groups/ifaromics/wenchichou/seqTU/data/s_1_1_qfilter_20q_95p.fastq.bwa.sam.multi.NA"
file_RNAseqSignals <- args[1]

#fileNamePrefix <- "treatment3"
#fileNamePrefix <- "control"
fileNamePrefix <- args[2]

file_gff <- "/home/wc102/project/tu/data/NC_009012.gff"
#file_gff <- args[2]
file_genomeSequence <- "/home/wc102/project/tu/data/NC_009012.fna"
#file_genomeSequence <- args[3]
dir_libsvm <- "/home/wc102/project/tu/bin/libsvm-3.17/"
#dir_libsvm <- args[5]

#file_454TU <- "/home/wc102/project/tu/data/454.Control.CoveredReads.Length.ChangeRatio"
file_454TU <- args[3]

#Clean SVM files
#052714WCCsystem("ls -ltr")        #Check ls
#052714WCCsystem("rm *SVM.*")       #Remove SVM files
#052714WCCsystem("ls -ltr")        #Check ls
#052714WCC

# open a statistic file
#sink("stat.txt")

# set some global variables
MeanQuantile = 0.1
FilterMatrixQuantile = 0.3
lambda = 0
MinimumMeanExpressionValue = 1
MinimumMedianExpressionValue = 1
MaximumGapProportionValue = 0.1
#Read EBCDs of each RNAseq data set (first column:Forward strand; second column:Reverse strand)
#SelectedRNAseqData<-read.table(file_mappedStrandSpecificRNAseqSignals, head=F)

cat("Reading the plot file ... \n")
SelectedRNAseqData<-read.table(file_RNAseqSignals, head=F)


#Read Gene Annotation 
cat("Reading gene annotations from the gff file ... \n")
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
cat("Done with reading gene annotations from the gff file ... \n")
cat("In this genome, there are\n 1)",nrow(AllForwardGenes),"forward genes and\n 2)",nrow(AllReverseGenes),"reverse genes.\n")

#Read bacterial Genome Sequence
cat("Reading genome sequence for the fna file ... \n")
seq=read.fasta(file_genomeSequence,seqonly=TRUE,as.string=TRUE)
seq=unlist(strsplit(seq[[1]], split="")) #Extract string
#length(seq) #Check genome length
cat("Done with reading genome sequence from the fna file ... \n")
cat("The input genome is",length(seq),"bp long.\n")

cat("Reading intergenic regions ... \n")
AllForwardIntergenicRegions=data.frame()
AllForwardIntergenicRegions=GenerateIGList(AllForwardGenes) #All intergenic regions on forward strand
AllReverseIntergenicRegions=data.frame()
AllReverseIntergenicRegions=GenerateIGList(AllReverseGenes) #All intergenic regions on reverse strand 
cat("Done with reading intergenic regions ... \n")
cat("In this genome, there are\n 1)",nrow(AllForwardIntergenicRegions),"forward intergenic regions and\n 2)",nrow(AllReverseIntergenicRegions),"reverse intergenic regions.\n")


#Get mean expression level, mean gaps and lenght for each gene on forward strand;
AllForwardGenesExpressionMean=NULL;
AllForwardGenesExpressionMedian=NULL;
AllForwardGenesExpressionSD=NULL;
AllForwardGenesGapsProportion=NULL;
AllForwardGenesGapsLongest=NULL;
AllForwardGenesLength=NULL;
for(i in 1:nrow(AllForwardGenes)){
	AllForwardGenesExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllForwardGenes, i));
	AllForwardGenesExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllForwardGenes, i));
	AllForwardGenesExpressionSD[i]=sd(GetRegionExp(SelectedRNAseqData,AllForwardGenes, i));
    AllForwardGenesGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllForwardGenes, i));
    AllForwardGenesGapsLongest[i]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData,AllForwardGenes, i));
    AllForwardGenesLength[i]=AllForwardGenes[i,3]-AllForwardGenes[i,2]+1;
};
#Get mean expression level, mean gaps and lenght for each gene on reverse strand;
AllReverseGenesExpressionMean=NULL;
AllReverseGenesExpressionMedian=NULL;
AllReverseGenesExpressionSD=NULL;
AllReverseGenesGapsProportion=NULL;
AllReverseGenesGapsLongest=NULL;
AllReverseGenesLength=NULL;
for(i in 1:nrow(AllReverseGenes)){;
	AllReverseGenesExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllReverseGenes, i));
	AllReverseGenesExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllReverseGenes, i));
	AllReverseGenesExpressionSD[i]=sd(GetRegionExp(SelectedRNAseqData,AllReverseGenes, i));
    AllReverseGenesGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllReverseGenes, i));
    AllReverseGenesGapsLongest[i]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData,AllReverseGenes, i));
    AllReverseGenesLength[i]=AllReverseGenes[i,3]-AllReverseGenes[i,2]+1;
}

#Get mean expression level, mean gaps and lenght for each intergenic region on forward strand;
AllForwardIntergenicRegionsExpressionMean=NULL;
AllForwardIntergenicRegionsExpressionMedian=NULL;
AllForwardIntergenicRegionsExpressionSD=NULL;
AllForwardIntergenicRegionsGapsProportion=NULL;
AllForwardIntergenicRegionsLength=NULL;
for(i in 1:(nrow(AllForwardIntergenicRegions))){;
	AllForwardIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
	AllForwardIntergenicRegionsExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
	AllForwardIntergenicRegionsExpressionSD[i]=sd(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
    AllForwardIntergenicRegionsGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllForwardIntergenicRegions, i));
    AllForwardIntergenicRegionsLength[i]=AllForwardIntergenicRegions[i,3]-AllForwardIntergenicRegions[i,2]+1;
};
#Get mean expression level, mean gaps and lenght for each intergenic region on reverse strand;
AllReverseIntergenicRegionsExpressionMean=NULL;
AllReverseIntergenicRegionsExpressionMedian=NULL;
AllReverseIntergenicRegionsExpressionSD=NULL;
AllReverseIntergenicRegionsGapsProportion=NULL;
AllReverseIntergenicRegionsLength=NULL;
for(i in 1:(nrow(AllReverseIntergenicRegions))){;
	AllReverseIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
	AllReverseIntergenicRegionsExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
	AllReverseIntergenicRegionsExpressionSD[i]=sd(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
    AllReverseIntergenicRegionsGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllReverseIntergenicRegions, i));
    AllReverseIntergenicRegionsLength[i]=AllReverseIntergenicRegions[i,3]-AllReverseIntergenicRegions[i,2]+1;
};


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

ForwardIntergenicRegionProportionNew=ForwardIntergenicRegionProportion
ReverseIntergenicRegionProportionNew=ReverseIntergenicRegionProportion

#Generate Forward simulated TU matrix
SimulatedForwardTUMatrixExpressionMean=data.frame()
SimulatedForwardTUMatrixExpressionSD=data.frame()
SimulatedForwardTUMatrixGapProportion=data.frame()
SimulatedForwardTUMatrixGapLongest=data.frame()
SimulatedForwardTUGeneFoldChange=NULL
num=1
for(i in 1:(nrow(AllForwardGenes))){
    	if ((AllForwardGenesExpressionMean[i] >= MinimumMeanExpressionValue) & (AllForwardGenesExpressionMedian[i] >= MinimumMedianExpressionValue) & (AllForwardGenesGapsProportion < MaximumGapProportionValue) & (AllForwardGenesGapsLongest[i] < 50)){
	        OneForwardSimulatedTU=GenerateOneSimulatedTU(AllForwardGenes, ForwardIntergenicRegionProportionNew, ForwardIntergenicRegionDeviate, i)
			for(j in 1:4){
				SimulatedForwardTUMatrixExpressionMean[num,j]=mean(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,j))
				SimulatedForwardTUMatrixExpressionSD[num,j]=sd(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,j))
				SimulatedForwardTUMatrixGapProportion[num,j]=mean(GetRegionGap(SelectedRNAseqData, OneForwardSimulatedTU,j))
				SimulatedForwardTUMatrixGapLongest[num,j]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData, OneForwardSimulatedTU,j))

			}
			SimulatedForwardTUGeneFoldChange[num] = FoldChange(mean(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,1)),mean(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
		num = num +1
    }
}

#Generate Reverse simulated TU matrix
SimulatedReverseTUMatrixExpressionMean=data.frame()
SimulatedReverseTUMatrixExpressionSD=data.frame()
SimulatedReverseTUMatrixGapProportion=data.frame()
SimulatedReverseTUMatrixGapLongest=data.frame()
SimulatedReverseTUGeneFoldChange=NULL
num=1
for(i in 1:(nrow(AllReverseGenes))){
    	if ((AllReverseGenesExpressionMean[i] >= MinimumMeanExpressionValue) & (AllReverseGenesExpressionMedian[i] >= MinimumMedianExpressionValue) & (AllReverseGenesGapsProportion < MaximumGapProportionValue) & (AllReverseGenesGapsLongest[i] < 50)){
	        OneReverseSimulatedTU=GenerateOneSimulatedTU(AllReverseGenes, ReverseIntergenicRegionProportionNew, ReverseIntergenicRegionDeviate, i)
			for(j in 1:4){
				SimulatedReverseTUMatrixExpressionMean[num,j]=mean(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,j))
				SimulatedReverseTUMatrixExpressionSD[num,j]=sd(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,j))
				SimulatedReverseTUMatrixGapProportion[num,j]=mean(GetRegionGap(SelectedRNAseqData, OneReverseSimulatedTU,j))
				SimulatedReverseTUMatrixGapLongest[num,j]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData, OneReverseSimulatedTU,j))
			}
			SimulatedReverseTUGeneFoldChange[num] = FoldChange(mean(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,1)),mean(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
		num = num +1
    }
}

#102014WCCwhich(is.na(SimulatedForwardTUMatrix[,1])) #check NaN elements
#102014WCCwhich(is.na(SimulatedForwardTUMatrix[,2])) #check NaN elements
#102014WCCwhich(is.na(SimulatedReverseTUMatrix[,1])) #check NaN elements
#102014WCCwhich(is.na(SimulatedReverseTUMatrix[,2])) #check NaN elements
#102014WCC
#102014WCCsummary(SimulatedForwardTUMatrix)
#102014WCCsummary(SimulatedReverseTUMatrix)
#102014WCC
#102014WCCnrow(SimulatedForwardTUMatrix)
#102014WCCnrow(SimulatedReverseTUMatrix)

#041514WCCOriginalSimulatedForwardTUMatrix=SimulatedForwardTUMatrix
#041514WCCOriginalSimulatedReverseTUMatrix=SimulatedReverseTUMatrix
#OriginalSimulatedTwoStrandsTUMatrix=rbind(SimulatedForwardTUMatrix,SimulatedReverseTUMatrix)
OriginalSimulatedTwoStrandsTUMatrix=cbind(
rbind(SimulatedForwardTUMatrixExpressionMean,SimulatedReverseTUMatrixExpressionMean),
rbind(SimulatedForwardTUMatrixExpressionSD,SimulatedReverseTUMatrixExpressionSD),
rbind(SimulatedForwardTUMatrixGapProportion,SimulatedReverseTUMatrixGapProportion),
rbind(SimulatedForwardTUMatrixGapLongest,SimulatedReverseTUMatrixGapLongest),
c(SimulatedForwardTUGeneFoldChange,SimulatedReverseTUGeneFoldChange)
)
colnames(OriginalSimulatedTwoStrandsTUMatrix)=c(
"ExpressionMean1","ExpressionMean2","ExpressionMean3","ExpressionMean4",
"ExpressionSD1","ExpressionSD2","ExpressionSD3","ExpressionSD4",
"GapProportion1","GapProportion2","GapProportion3","GapProportion4",
"GapLongest1","GapLongest2","GapLongest3","GapLongest4",
"GeneFoldChange"
)

#colnames(OriginalSimulatedTwoStrandsTUMatrix)=c("GapPercentage", "FoldChange", "MeanExpressionIntergenicRegion", "MeanExpressionFiveEnd", "")


#geting real cases from the genome
TargetForwardTUMatrixExpressionMean=data.frame()
TargetForwardTUMatrixExpressionSD=data.frame()
TargetForwardTUMatrixGapProportion=data.frame()
TargetForwardTUMatrixGapLongest=data.frame()
TargetForwardTUGeneFoldChange=NULL
num=1
for(i in 1:(nrow(AllForwardGenes)-1)){
	OneForwardTargetTU=GenerateOneTargetTU(AllForwardGenes,i)
	for(j in 1:4){
		TargetForwardTUMatrixExpressionMean[num,j]=mean(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,j))
		TargetForwardTUMatrixExpressionSD[num,j]=sd(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,j))
		TargetForwardTUMatrixGapProportion[num,j]=mean(GetRegionGap(SelectedRNAseqData, OneForwardTargetTU,j))
		TargetForwardTUMatrixGapLongest[num,j]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData, OneForwardTargetTU,j))
	}
	TargetForwardTUGeneFoldChange[num] = FoldChange(mean(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,1)),mean(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
    num = num +1
}
#which(is.na(TargetForwardTUMatrix[,1])) #check NaN elements
#which(is.na(TargetForwardTUMatrix[,2])) #check NaN elements
TargetReverseTUMatrixExpressionMean=data.frame()
TargetReverseTUMatrixExpressionSD=data.frame()
TargetReverseTUMatrixGapProportion=data.frame()
TargetReverseTUMatrixGapLongest=data.frame()
TargetReverseTUGeneFoldChange=NULL
num=1
for(i in 1:(nrow(AllReverseGenes)-1)){
	OneReverseTargetTU=GenerateOneTargetTU(AllReverseGenes,i)
	for(j in 1:4){
		TargetReverseTUMatrixExpressionMean[num,j]=mean(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,j))
		TargetReverseTUMatrixExpressionSD[num,j]=sd(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,j))
		TargetReverseTUMatrixGapProportion[num,j]=mean(GetRegionGap(SelectedRNAseqData, OneReverseTargetTU,j))
		TargetReverseTUMatrixGapLongest[num,j]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData, OneReverseTargetTU,j))
	}
	TargetReverseTUGeneFoldChange[num] = FoldChange(mean(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,1)),mean(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
    num = num +1
}

#Output TargetSVM format
TargetPositiveTUMatrix=cbind(
TargetForwardTUMatrixExpressionMean,
TargetForwardTUMatrixExpressionSD,
TargetForwardTUMatrixGapProportion,
TargetForwardTUMatrixGapLongest,
TargetForwardTUGeneFoldChange
)
TargetNegativeTUMatrix=cbind(
TargetReverseTUMatrixExpressionMean,
TargetReverseTUMatrixExpressionSD,
TargetReverseTUMatrixGapProportion,
TargetReverseTUMatrixGapLongest,
TargetReverseTUGeneFoldChange
)

TargetTwoStrandsTUMatrix=cbind(
rbind(TargetForwardTUMatrixExpressionMean,TargetReverseTUMatrixExpressionMean),
rbind(TargetForwardTUMatrixExpressionSD,TargetReverseTUMatrixExpressionSD),
rbind(TargetForwardTUMatrixGapProportion,TargetReverseTUMatrixGapProportion),
rbind(TargetForwardTUMatrixGapLongest,TargetReverseTUMatrixGapLongest),
c(TargetForwardTUGeneFoldChange,TargetReverseTUGeneFoldChange)
)
colnames(TargetPositiveTUMatrix)=colnames(TargetNegativeTUMatrix)=colnames(TargetTwoStrandsTUMatrix)=c(
"ExpressionMean1","ExpressionMean2","ExpressionMean3","ExpressionMean4",
"ExpressionSD1","ExpressionSD2","ExpressionSD3","ExpressionSD4",
"GapProportion1","GapProportion2","GapProportion3","GapProportion4",
"GapLongest1","GapLongest2","GapLongest3","GapLongest4",
"GeneFoldChange"
)

#454 data 
TUverifiedBy454<-read.table(file_454TU, header=T)
verifiedTU_Reads_Length_Ratio=subset(TUverifiedBy454,select=c(coveredreads,length,changeRatio))
Temp454PositiveIndex<-which(verifiedTU_Reads_Length_Ratio$coveredreads>=1 & verifiedTU_Reads_Length_Ratio$changeRatio <=2)
Temp454NegativeIndex<-which(verifiedTU_Reads_Length_Ratio$coveredreads<1 | verifiedTU_Reads_Length_Ratio$changeRatio >2)

#get a half of Temp454PositiveIndex into training 
randomizedTemp454PositiveIndex<-sample(Temp454PositiveIndex)
#get the first half of randomized index for training 
Temp454PositiveIndexInTraining <- randomizedTemp454PositiveIndex[1:ceiling(length(Temp454PositiveIndex)/2)]
#get the second half of randomized index for verification
Temp454PositiveIndexInTest <- randomizedTemp454PositiveIndex[(ceiling(length(Temp454PositiveIndex)/2)+1):length(Temp454PositiveIndex)]
sink(paste("Temp454PositiveIndexInTest",fileNamePrefix,sep=""))
for(i in 1:length(Temp454PositiveIndex)){
	cat(i,"\t",Temp454PositiveIndex[i],"\n",sep="")
}
sink()






#only keep fold chang < 10 in pseudo TU
#SimulatedPositiveTUMatrix=SimulatedPositiveTUMatrix[which(SimulatedPositiveTUMatrix[,2] < 10),]
SimulatedPositiveTUMatrix = OriginalSimulatedTwoStrandsTUMatrix
#SimulatedPositiveTUMatrix=rbind(SimulatedPositiveTUMatrix,TargetTwoStrandsTUMatrix[Temp454PositiveIndexInTraining,])
dim(SimulatedPositiveTUMatrix)
head(SimulatedPositiveTUMatrix)

#SimulatedNegativeTUMatrix = TargetTwoStrandsTUMatrix[Temp454NegativeIndex,]
SimulatedNegativeTUMatrix = TargetTwoStrandsTUMatrix
SimulatedNegativeTUMatrix = SimulatedNegativeTUMatrix[which(SimulatedNegativeTUMatrix[,"GeneFoldChange"] >= 7),]   # fold changes > 7
#SimulatedNegativeTUMatrix = SimulatedNegativeTUMatrix[which(SimulatedNegativeTUMatrix[,"ExpressionSD4" ]> 100),] 
#SimulatedNegativeTUMatrix = SimulatedNegativeTUMatrix[which(SimulatedNegativeTUMatrix[,"GapLongest2"] > 100),] 
dim(SimulatedNegativeTUMatrix)
head(SimulatedNegativeTUMatrix)


#select negative data set
#SimulatedNegativeTUMatrix = TargetTwoStrandsTUMatrix[which(TargetTwoStrandsTUMatrix[,2] >=5),]







#Final filtering to prepare training sets
#041514WCCSimulatedForwardTUMatrix=OriginalSimulatedForwardTUMatrix[which(OriginalSimulatedForwardTUMatrix[,1] <= (quantile(OriginalSimulatedForwardTUMatrix[,1], FilterMatrixQuantile))),]
#041514WCCSimulatedReverseTUMatrix=OriginalSimulatedReverseTUMatrix[which(OriginalSimulatedReverseTUMatrix[,1] <= (quantile(OriginalSimulatedReverseTUMatrix[,1], FilterMatrixQuantile))),]
#041514WCCSimulatedForwardTUMatrix=OriginalSimulatedForwardTUMatrix[which(OriginalSimulatedForwardTUMatrix[,2] <= 10),]
#041514WCCSimulatedReverseTUMatrix=OriginalSimulatedReverseTUMatrix[which(OriginalSimulatedReverseTUMatrix[,2] <= 10),]

#041514WCCcat("Summary of SimulatedForwardTUMatrix: ",summary(SimulatedForwardTUMatrix), "\n")
#041514WCCcat("Summary of SimulatedReverseTUMatrix: ",summary(SimulatedReverseTUMatrix), "\n")
#041514WCCcat("Length of SimulatedForwardTUMatrix: ", nrow(SimulatedForwardTUMatrix), "\n")
#041514WCCcat("Length of SimulatedReverseTUMatrix: ", nrow(SimulatedReverseTUMatrix), "\n")




#Output SVM format for Training dataset
iter = c(1:ncol(TargetTwoStrandsTUMatrix))

sink(paste("SVM.TwoStrands.Training.",fileNamePrefix,sep=""))
#SimulatedTwoStrandsTUMatrix=rbind(SimulatedForwardTUMatrix, SimulatedReverseTUMatrix)
for(i in 1:nrow(SimulatedPositiveTUMatrix)){
        cat("+1")
#   for (j in 1:ncol(SimulatedPositiveTUMatrix)){
    for (j in iter){
        cat("\t",j,":",SimulatedPositiveTUMatrix[i,j],sep="")
    }
    cat("\n")
}
for(i in 1:nrow(SimulatedNegativeTUMatrix)){
        cat("-1")
#   for (j in 1:ncol(SimulatedNegativeTUMatrix)){
    for (j in iter){
        cat("\t",j,":",SimulatedNegativeTUMatrix[i,j],sep="")
    }
    cat("\n")
}
sink()


#Output SVM format for Target dataset
iter = c(1:ncol(TargetTwoStrandsTUMatrix))
sink(paste("SVM.TwoStrands.Target.",fileNamePrefix,sep=""))
for(i in 1:nrow(TargetPositiveTUMatrix)){
        cat("+1")
    for (j in iter){
        cat("\t",j,":",TargetPositiveTUMatrix[i,j],sep="")
    }
    cat("\n")
}
for(i in 1:nrow(TargetNegativeTUMatrix)){
        cat("-1")
    for (j in iter){
        cat("\t",j,":",TargetNegativeTUMatrix[i,j],sep="")
    }
    cat("\n")
}
sink()








#042214#Output SVM format
#042214sink("SVM.TwoStrands.Training")
#042214#SimulatedTwoStrandsTUMatrix=rbind(SimulatedForwardTUMatrix, SimulatedReverseTUMatrix)
#042214for(i in 1:nrow(SimulatedPositiveTUMatrix)){
#042214		cat("+1")
#042214#	for (j in 1:ncol(SimulatedPositiveTUMatrix)){
#042214	for (j in 2){
#042214	    cat("\t",j,":",SimulatedPositiveTUMatrix[i,j],sep="")
#042214	}
#042214	cat("\n")
#042214}
#042214for(i in 1:nrow(SimulatedNegativeTUMatrix)){
#042214		cat("-1")
#042214#	for (j in 1:ncol(SimulatedNegativeTUMatrix)){
#042214	for (j in 2){
#042214	    cat("\t",j,":",SimulatedNegativeTUMatrix[i,j],sep="")
#042214	}
#042214	cat("\n")
#042214}
#042214sink()



#041514WCCsink("SVM.Forward.Training")
#041514WCCfor(i in 1:nrow(SimulatedForwardTUMatrix)){
#041514WCC    cat("+1\t1:", SimulatedForwardTUMatrix[i,1], "\t2:", SimulatedForwardTUMatrix[i,2], "\n", sep="")
#041514WCC}
#041514WCCsink()
#041514WCC
#041514WCCsink("SVM.Reverse.Training")
#041514WCCfor(i in 1:nrow(SimulatedReverseTUMatrix)){
#041514WCC    cat("+1\t1:", SimulatedReverseTUMatrix[i,1], "\t2:", SimulatedReverseTUMatrix[i,2], "\n", sep="")
#041514WCC}
#041514WCCsink()

#Run system() command in R script
#Check ls
#system("ls -ltr", wait = FALSE)
#system("ls -ltr SVM.*")
#Rename SVM format files
#system("for i in SVM.*; do j=\"_081112\"; mv $i $i$j; done")

#052714WCC
#052714WCC#Scale SVM format files
#052714WCCsystem(paste("for i in SVM.*.Training; do ",dir_libsvm,"svm-scale -s $i.scaled.Ruler $i > $i.scaled; done", sep=""))
#052714WCCTwoStrandsScaledFileName<-system("ls SVM.TwoStrands.Training.scaled", intern=TRUE)
#052714WCCTwoStrandsRulerFileName=paste(TwoStrandsScaledFileName,".Ruler", sep="")
#052714WCCTwoStrandsModelFileName=paste(TwoStrandsScaledFileName,".Model", sep="")
#052714WCCTwoStrandsPredictFileName=paste(TwoStrandsScaledFileName,".Predict", sep="")
#052714WCCweight = nrow(SimulatedPositiveTUMatrix)/nrow(SimulatedNegativeTUMatrix)
#052714WCC
#052714WCC#Check scaled SVM format files
#052714WCCsystem(paste("for i in SVM.*.scaled; do ",dir_libsvm,"tools/checkdata.py $i; done", sep=""))
#052714WCC
#052714WCC#Find the best parameters and results
#052714WCC#"for i in SVM.*.scaled; do python /home/wc102/project/tu/bin/libsvm-3.17/tools/grid_WCC081812.py -s 2 -t 2 -v 5 -m 300 -svmtrain /home/wc102/project/tu/bin/libsvm-3.17/svm-train $i; done"
#052714WCC#/home/wc102/project/tu/bin/libsvm-3.17/svm-train -s 0 -t 2 -w1 9 -w-1 1 -v 5 -m 300 -c 0.03125 -g 256 SVM.TwoStrands.Training.scaled
#052714WCCGridResults<-system(paste("for i in SVM.*.scaled; do python ",dir_libsvm,"tools/grid_WCC081812.py -s 0 -t 2 -w1 ",weight," -w-1 1 -v 5 -m 300 -svmtrain ",dir_libsvm,"svm-train $i; done",sep=""), intern=TRUE)
#052714WCC
#052714WCC#Check Training results
#052714WCC#system("for i in *scaled.out; do echo $i; sort -n -k 3 $i|tail -1; done")
#052714WCCBestParameters=grep("local", GridResults, invert=TRUE, value=TRUE)
#052714WCCTwoStrandsC=2^(as.numeric(strsplit(BestParameters, " ")[[1]][1]))
#052714WCCTwoStrandsG=2^(as.numeric(strsplit(BestParameters, " ")[[1]][2]))
#052714WCC
#052714WCC#system(paste(dir_libsvm,"svm-train -s 2 -t 2 -n ",TwoStrandsN," -g ",TwoStrandsG," -n ",TwoStrandsN," ",TwoStrandsScaledFileName," ",TwoStrandsModelFileName, sep=""))
#052714WCCsystem(paste(dir_libsvm,"svm-train -s 0 -t 2 -c ",TwoStrandsC," -g ",TwoStrandsG," -w1 ",weight," -w-1 1 ",TwoStrandsScaledFileName," ",TwoStrandsModelFileName, sep=""))
#052714WCCSimulatedTwoStrandsTUResult<-system(paste(dir_libsvm,"svm-predict ",TwoStrandsScaledFileName," ",TwoStrandsModelFileName," ",TwoStrandsPredictFileName, sep=""), intern=TRUE)
#052714WCC
#052714WCCcat ("Best parameters:",BestParameters,"\n",sep= "\t")
#052714WCCcat ("SVM Accuracy of simulated forward TUs:",SimulatedTwoStrandsTUResult,"\n",sep= "\t")
#052714WCC
#052714WCC
#052714WCC#041514WCCsystem(paste(dir_libsvm,"svm-train -s 2 -t 2 -n ",ForwardN," -g ",ForwardG," -n ",ForwardN," ",ForwardScaledFileName," ",ForwardModelFileName, sep=""))
#052714WCC#041514WCCsystem(paste(dir_libsvm,"svm-train -s 2 -t 2 -n ",ReverseN," -g ",ReverseG," ",ReverseScaledFileName," ",ReverseModelFileName, sep=""))
#052714WCC#041514WCC
#052714WCC#041514WCCSimulatedForwardTUResult<-system(paste(dir_libsvm,"svm-predict ",ForwardScaledFileName," ",ForwardModelFileName," ",ForwardPredictFileName, sep=""), intern=TRUE)
#052714WCC#041514WCCSimulatedReverseTUResult<-system(paste(dir_libsvm,"svm-predict ",ReverseScaledFileName," ",ReverseModelFileName," ",ReversePredictFileName, sep=""), intern=TRUE)
#052714WCC
#052714WCC
#052714WCCTargetForwardTUMatrix=data.frame()
#052714WCCnum=1
#052714WCCfor(i in 1:(nrow(AllForwardGenes)-1)){
#052714WCC	OneForwardTargetTU=GenerateOneTargetTU(AllForwardGenes,i)
#052714WCC	TargetForwardTUMatrix[num,1]=GetGapPercentage(SelectedRNAseqData, OneForwardTargetTU,2) # 2 is the intergenic reiong of a target TU
#052714WCC    TargetForwardTUMatrix[num,2]=FoldChange(median(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,1)),median(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
#052714WCC    num = num +1
#052714WCC}
#052714WCCwhich(is.na(TargetForwardTUMatrix[,1])) #check NaN elements
#052714WCCwhich(is.na(TargetForwardTUMatrix[,2])) #check NaN elements
#052714WCC
#052714WCC
#052714WCCTargetReverseTUMatrix=data.frame()
#052714WCCnum=1
#052714WCCfor(i in 1:(nrow(AllReverseGenes)-1)){
#052714WCC	OneReverseTargetTU=GenerateOneTargetTU(AllReverseGenes,i)
#052714WCC	TargetReverseTUMatrix[num,1]=GetGapPercentage(SelectedRNAseqData, OneReverseTargetTU,2) # 2 is the intergenic reiong of a target TU
#052714WCC    TargetReverseTUMatrix[num,2]=FoldChange(median(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,1)),median(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
#052714WCC    num = num +1
#052714WCC}
#052714WCCwhich(is.na(TargetReverseTUMatrix[,1])) #check NaN elements
#052714WCCwhich(is.na(TargetReverseTUMatrix[,2])) #check NaN elements
#052714WCC
#052714WCC#Output TargetSVM format
#052714WCCTargetTwoStrandsTUMatrix=rbind(TargetForwardTUMatrix, TargetReverseTUMatrix)
#052714WCCsink("TargetSVM.TwoStrands")
#052714WCCfor(i in 1:nrow(TargetTwoStrandsTUMatrix)){
#052714WCC    cat("+1\t1:", TargetTwoStrandsTUMatrix[i,1], "\t2:", TargetTwoStrandsTUMatrix[i,2], "\n", sep="")
#052714WCC}
#052714WCCsink()
#052714WCC
#052714WCC
#052714WCC#Rename SVM format files
#052714WCC#system("for i in TargetSVM.*; do j=\"_081112\"; mv $i $i$j; done")
#052714WCC#Scale SVM format files
#052714WCC
#052714WCCsystem(paste(dir_libsvm,"svm-scale -r ",TwoStrandsRulerFileName," TargetSVM.TwoStrands > TargetSVM.TwoStrands.scaled", sep=""))
#052714WCCTargetTwoStrandsScaledFileName<-system("ls TargetSVM.TwoStrands.scaled", intern=TRUE)
#052714WCCTargetTwoStrandsPredictFileName=paste(TargetTwoStrandsScaledFileName,".Predict", sep="")
#052714WCC
#052714WCC
#052714WCC
#052714WCC
#052714WCC
#052714WCC
#052714WCCTargetTwoStrandsTUResult<-system(paste(dir_libsvm,"svm-predict ",TargetTwoStrandsScaledFileName," ",TwoStrandsModelFileName," ",TargetTwoStrandsPredictFileName, sep=""), intern=TRUE)
#052714WCC
#052714WCC
#052714WCC
#052714WCCcat ("SVM Accuracy of simulated forward TUs:",SimulatedTwoStrandsTUResult,"\n",sep= "\t")
#052714WCCcat ("% TUs:", TargetTwoStrandsTUResult, "\n",sep= "\t")
#052714WCC
