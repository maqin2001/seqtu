Instruction for installing and using the scripts of SeqTU 
Installation 
1.	Download and install libSVM:
 The current release (Version 3.17, April 2013) of LIBSVM can be obtained by downloading the zip file or tar.gz file from http://www.csie.ntu.edu.tw/~cjlin/libsvm/ 
2.	Copy "grid_WCC081812.py" into libsvm-3.17/tools  directory. 
3.	You also need to have R installed and four R packages installed. The four R packages are
a.	library(grid)
b.	library(gridBase)
c.	library(ggplot2)
d.	library(seqinr)
Usage
1.	Your strand-specific RNA-seq data need to be mapped to the corresponding genome. In our study, we use Clostridium Thermocellum genome as a reference genome for mapping.  
2.	The mapped results need to be re-formatted to two-column single-base signals. Please see a sample file named "ssRNAseq.forward.reversed.signals". The first line of the file presents the first position of the genome, and the two numbers separated by TAB are RNA-seq signals of forward strand and reversed strand, respectively.
3.	With the "ssRNAseq.forward.reversed.signals" ready, you also need prepare a GFF file and a FASTA file of the corresponding genome for running an R script named "seqTU_101413.r". In our study, we used NC_009012.gff  and NC_009012.fna as the GFF file and the FASTA file.
4.	The TU identification will be performed by running "seqTU_101413.r". The results then can be post-processed by another R script named "print.final.table_101413.r".
Contact
If you have any questions, please contact us 
1.	Wen-Chi Chou  
wcc957@gmail.com 
2.	Qin Ma 
maqin2001@gmail.com 


