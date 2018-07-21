# ROSeq

ROSeq evaluates differential expression between two subpopulations of single cells by modeling the rank distribution of read counts for each gene. 

## Installation

To get started with ROSeq, *git clone* or download the complete ROSeq [repository](https://github.com/malaalam/ROSeq) on to your PC. Next, double click the executable file *ROSeq/06_installation/OS/for_redistribution/MyAppInstaller_mcr.exe* for installing the software on your respective operating system. 

## Tutorial 

This tutorial below, provides instructions to upload the (un) normalized read count file corresponding to single cells coming from the individuals *NA19098* and *NA19101*, in the [tung2017](https://www.ncbi.nlm.nih.gov/pubmed/28045081) data set; and subsequently analyse the two subpopulations for differential expression. 

The number of single cells equals 288 for each individual. After filtering with a threshold equal to six, low quality genes are eliminated and the number of genes considered for differential expression analysis reduces from 19027 to 16087. This is also prompted as a message in the Status Panel. 

Following would constitute as a complete set of steps from uploading read count data to saving results as a '*.csv' file:

1. Click on the radio button - *Upload Raw Read Count file*. This would produce an *Open File* button, for the user to upload the respective file. 
    
2. Click on the *Open File* button, to select the raw read count file *single_counts_no_filter_no_norm.csv*. Prior to successful selection of the file, an import wizard window would appear. Click on the buttons *Next* and *Finish*.
    
3. Next, the number of genes in the uploaded file will be suggested (19027). The user could alter the number of genes, in case you wish to look only at the first few (10) genes.

4. Decide on the Filter Threshold. A filter threshold equal to 6, would mean that all genes with the number of non-zero read counts less than 6, would not be considered during the subsequent differential analaysis, as they would be considered *low-quality genes*.
    
5. Select the starting (1) and ending (288) column index corresponding to the first group (*NA19098*). Select the starting (289) and ending (576) column of the second group similarly (*NA19101*). 

6. If all steps were performed so far, then the *Solve DE* Button becomes enabled. Click on it to run the analysis for differential expression. This would open up a green progress bar for the Subpopulation One, followed by another for Subpopulation Two.
    
7. Once the analysis is complete, a message (*DE analysis completed successfully.*) would appear in the *Status* Panel. At this stage, the ROSeq software would look as in the Figure below.    

8. Lastly, click on *Save Results* button to save the variables [pValue, adjusted pValue and log2(fold count)] as three columns in a *.csv file format.

![Screenshot of ROSeq]({{ "/05_images/ROSeq.PNG" | absolute_url }})

## Downloads

The Paper and Supplementary Material is available for download here.

The R code used for producing differential analysis results with methods like SCDE, ROTS, MAST and Wilcoxon-Rank Sum is available for download here.

## Datasets

The [tung2017](https://www.ncbi.nlm.nih.gov/pubmed/28045081) dataset was downloaded from [jdblischak](https://github.com/jdblischak/singleCellSeq)'s Repository. 

The [trapnell2014](http://www.nature.com/articles/nbt.2859) dataset was downloaded from 

