### Plot ChIP-seq peaks

#####Rscript to plot ChIP-Seq peaks. 

Dependencies: 
[rtracklayer](http://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)
ggplot2

Example usage: 

```
source("plot_chipPeaks.R")
plotPeaks(bigWig = bw,peaks_bed = 'peaks.bed',chr = 'chr10',start = 74020000,end = 74050000,isBigWigFile = T,showGeneStruct = T)
```

`bigWig` = Input file in [bigwig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) format.

`peaks_bed` = Peak regions in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.

`chr` = Chromosome name to visualize.

`start` = start position from which to plot.

`end` = end position till which to plot.

`isBigWigFile` = If input file is is bigwig. LOGICAL - defaut FALSE.

`showGeneStruct` = Whether to plot genes in located in the plot region. Default FALSE.


#####Sample output from above command.


![sample peak](https://github.com/PoisonAlien/Plot_chipPeaks/blob/master/Rplot01.png)

