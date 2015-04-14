plotPeaks = function(bigWig,peaks_bed,chr,start,end,isBigWigFile = FALSE,showGeneStruct = FALSE){
  
  library(rtracklayer)
  library(ggplot2)
  
  ################################# multiplot function ################################################################
  ## from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ ##
  
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  #################################################################################################################################################
  
  if(isBigWigFile){
    #read bigwig file
    print(paste("Reading ",bigWig,sep = ":"))
    bw = import.bw(bigWig)  
  }
  
  #read peaks file (in bed format)
  peaks = read.delim(peaks_bed,header = F,stringsAsFactors = F)
  #convert it into GRanges object.
  peaks = GRanges(seqnames = peaks[,1],ranges = IRanges(start = peaks[,2],end = peaks[,3]))
  
  #form a GRanges area for plotting area.
  plot.region = GRanges(seqnames = chr,ranges = IRanges(start = c(start),end = c(end)))
  
  #Get scores for regions to plot.
  reg = subsetByOverlaps(query = bw,subject = plot.region)
  reg.df = data.frame(midpoint = start(reg)+(end(reg) - start(reg))/2,score = reg$score)
  
  #Get peaks to highlight.
  peaksToShow = as.data.frame(subsetByOverlaps(query = peaks,subject = plot.region))
  
  #plot main peaks
  p1 = ggplot(reg.df, aes(x=midpoint, y=score)) +geom_area(colour="black", fill="blue", alpha=.2)+ylim(0,max(reg.df$score)+2)+xlab("position")+ylab("rpm")+theme_bw()+theme_minimal()+theme(axis.line = element_line(colour = "black"),panel.grid.major=element_blank())+xlim(start,end)
  
  #plot bed regions
  p2 = ggplot(data = peaksToShow,aes(ymin = -0.5,ymax = 0.5,xmin = start,xmax = end,aplha = 0.2))+geom_rect()+theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())
  p2 = p2 + xlim(start,end)
  
  if(showGeneStruct){
    
    #start an ucsc browser session
    mySession = browserSession("UCSC")
    genome(mySession) = "hg19"
    
    #Download all genes which lie in the plot region.
    geneStruct = getTable(ucscTableQuery(mySession, track="RefSeq Genes",range=plot.region, table="refGene"))
    geneStruct = geneStruct[1,] #For now this will only plot the fisrt entry :(
    
    #Get start and end locus of each exon and make a data frame.
    exon.start = as.numeric(unlist(strsplit(x = as.character(geneStruct$exonStarts),split = ",",fixed = T)))
    exon.end = as.numeric(unlist(strsplit(x = as.character(geneStruct$exonEnds),split = ",",fixed = T)))
    exon.df = data.frame(start = as.numeric(exon.start),end = as.numeric(exon.end),gene = as.character(geneStruct$name2),exon = paste("exon",1:length(exon.start),sep = "_"),stringsAsFactors = F)
    
    #plot the gene structure.
    p3 = ggplot(data = exon.df,aes(ymin = -0.5,ymax = 0.5,xmin = start,xmax = end,aplha = 0.2))+geom_rect()+ylim(-2,1)+xlim(start,end)+geom_rect(data = geneStruct,aes(ymin = -0.1,ymax = 0.1,xmin = txStart,xmax = txEnd,aplha = 0.2))+theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())
    p3 = p3+annotate("text",label = as.character(geneStruct$name2),x = start+500,y = -1) #Label it.
    
    #form a layout for plotting.
    layout = matrix(c(1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3), nrow = 8, byrow = TRUE)
    
    multiplot(plotlist = list(p2,p1,p3), layout = layout)
  } else{
    layout = matrix(c(1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2), nrow = 7, byrow = TRUE)
    
    multiplot(plotlist = list(p2,p1), layout = layout)
  }
}
