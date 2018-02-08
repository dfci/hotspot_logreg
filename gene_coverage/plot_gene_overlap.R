###################################
#'	IntelCCC gene overlap and plotting function: This code accepts a list of mafs. The mafs are expected to have the folloing columns: 
#'		-Hugo_Symbol
#'		-Center
#'		-Tumour_Sample_Barcode
#'
#'	Input: 
#'	list_file: a text file containing a list of maf files to read
#'	output_dir: directory to save output plot and text file
#'	sort: (TRUE) sort the genes according to total sample count, (FALSE) sort genes alphatically 
#'	plot: (TRUE) generates a stacked boxplot, (FALSE) no plot is generated 
#'
#'	Output: a text file that summarizes coverage per gene along with a plot to visualize the gene overlap across centers.
#'	
#'	Date: June 9th 2017
#' 	Fouad Yousif, fouad.yousif@oicr.on.ca


#load required libraries
library(reshape)
library(ggplot2)

generateGeneOverlap<-function(
	list_file,
	output_dir,
	sort,
	plot
	){

#open the file that contains the list of mafs to read
	mafs <- read.delim(list_file,header=T)

#count the number of maf files
	l <- nrow(mafs)

#loop through the files and read them
	for (i in 1:l) {
  		x <- read.delim(as.character(mafs$filename[i]),header=T)
  		if (i==1){
  			genes<-data.frame(Hugo_Symbol= x$Hugo_Symbol,Center=x$Center,Sample=x$Tumor_Sample_Barcode)
  		}else{
  			genes<-rbind(genes,data.frame(Hugo_Symbol= x$Hugo_Symbol,Center=x$Center,Sample=x$Tumor_Sample_Barcode))
  		}
	}


#generate a summary
	genes_collapsed <- data.frame(as.matrix(t(table(genes$Center,genes$Hugo_Symbol))))
	names(genes_collapsed) <- c('Hugo_Symbol','Center','Sample')
	genes_collapsed<-genes_collapsed[order(genes_collapsed$Sample,decreasing=TRUE),]
	if(sort == TRUE){
		genes_collapsed$Hugo_Symbol <- reorder(genes_collapsed$Hugo_Symbol,sort(genes_collapsed$Sample,decreasing=TRUE))
	}else{}


#write the table
	write.table(genes_collapsed,paste(output_dir,"gene_overlap_summary.txt",sep=""),quote=FALSE,sep="\t",row.names=F)

#generate the plot
	if(plot==TRUE){
		plot <- ggplot(genes_collapsed, aes(x = Hugo_Symbol, y = Sample, fill = Center,order=Hugo_Symbol)) + geom_bar(stat = "identity")
		plot + theme(axis.text.x=element_text(angle=90, hjust=1,size=4))
		ggsave(paste(output_dir,"gene_overlap_summary.png",sep=""), width = 12, height = 5)
		}
}


#call the function
generateGeneOverlap("./maf_list.txt","./",TRUE,TRUE)
