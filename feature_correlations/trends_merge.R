#!/usr/bin/env Rscript

#Taking in arguments
library("optparse")

option_list = list(
  make_option(c("-i", "--inputfolder"), type="character", default=NULL,
              help="input folder",metavar="character"),
  make_option(c("-f", "--files"), type="character", default=NULL,
              help="filepattern",metavar="character"),
  make_option(c("-y", "--ymax"), type="character", default=4,
              help="y axis limit",metavar="character"),
  make_option(c("-o", "--output_prefix"), type="character", default="out",
              help="output file name [default %default]", metavar="character"),
  make_option(c("-b", "--tick_interval"), type="numeric", default=1,
              help="interval between major Y axis ticks [default %default]", metavar="numeric"),
  make_option(c("-q", "--nolab"), action="store_true", default=FALSE,
              help="suppress axis text labels; append _nolab to filename [default %default]")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
inp=opt$inputfolder; print(inp)
file=opt$files; print(file)
ymax=opt$ymax; if (ymax > 1) {ymax= as.integer(ymax)} else {ymax= round(as.numeric(ymax), digits = 2)};  print(ymax)
out=opt$output_prefix; print(out)
tick_interval=opt$tick_interval; print(tick_interval)

#Checking for required input
writeLines("\n...Checking input...\n")
if (is.null(opt$files)){
  print_help(opt_parser)
  stop("Please specify file basename.n", call.=FALSE)
}
if (is.null(opt$inputfolder)){
  print_help(opt_parser)
  stop("Please specify the input folder to be used to make matrix.n", call.=FALSE)
}


#####Main Script
#1.Calling required Packages
writeLines("\n...Calling required Package(s)...\n")

library(GenomicRanges, quiet = T)
library(plyranges, quiet = T)
library(ggplot2, quiet = T)
library(tidyr, quiet = T)


#2.Defining functions if any
writeLines("\n...Defining required functions...\n")

#3. Preprocessing input files
writeLines("\n...Processing Input files...\n")

#filenames = list.files(inp, pattern=file)
filenames=Filter(function(x) grepl("stats.tsv",x),list.files(inp,pattern=file))

#4.Main code
writeLines("\n...Executing Main code...\n")
colors <- c(
  "rA"="#D55E00",
  "rC"="#0072B2",
  "rG"="#F0C742",
  "rU"="#009E73"
)
shapes <- c(
  "rA"=16,  # circle
  "rC"=17,  # triangle
  "rG"=15,  # square
  "rU"=23   # diamond (shape 23 renders at consistent visual weight with other shapes)
)

#Merge matrix
for (p in filenames) {
  writeLines(paste("Working on", p))
  seq.df=read.table(paste(inp,p,sep="/"),sep='\t', header=T)
  # Extract nucleotide prefix (rA/rC/rG/rU) from the name column
  ribo=strsplit(unique(seq.df$name),split='_')[[1]][1]
  color=unname(colors[ribo])
  colnames(seq.df)[4]=ribo #using ymean
  if (p == filenames[1]) {
    mat=seq.df[c('groups','Featurex')]
  }
  mat=merge(mat,seq.df[c('groups',ribo)])
}

data=gather(mat, key='ribo',value='ymean',-c(Featurex,groups))
suffix <- if (isTRUE(opt$nolab)) "_nolab" else ""
svg(paste(out,'/',file,'_all_linear_trend', suffix, '.svg', sep=""), width = 1.75, height = 1.75)
#png(paste(out,'/',file,'_all_linear_trend.png', sep=""), width = 1.75, height = 1.75, units = "in", res=600, type="cairo", bg="transparent")
#data$groups<-as.factor(data$groups)
dodge=position_dodge(width = 0.8)
plot<-ggplot(data, aes(x=groups,y=ymean, fill=ribo, colour=ribo, shape=ribo, size=ribo, stroke=ribo))+
  geom_point()+
  geom_smooth(method=lm, size = 0.3, se=FALSE)+
  #stat_poly_line() +
  #geom_ribbon(aes(ymin=low,ymax=high),alpha=0.2, colour=NA)+
  scale_fill_manual(values=as.character(colors[unique(data$ribo)]))+
  scale_color_manual(values=as.character(colors[unique(data$ribo)]))+
  scale_shape_manual(values=shapes[unique(data$ribo)])+
  scale_size_manual(values=c("rA"=0.3,"rC"=0.3,"rG"=0.3,"rU"=0.55), guide="none")+
  scale_discrete_manual("stroke", values=c("rA"=0.5,"rC"=0.5,"rG"=0.5,"rU"=0), guide="none")+
  scale_x_continuous(expand = c(0, 0), breaks=mat$groups, labels=as.character(mat$groups),limits=c(min(mat$groups), max(mat$groups)+0.5))+
  #scale_x_continuous(expand = c(0, 0), breaks=unique(mat$groups), labels=mat$label) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,ymax+ymax*0.0035),breaks=seq(0,ymax,by=tick_interval), name = "rNMP Enrichment")+
  #annotate("text",x=6,y=3.3,label=bquote(italic(R) == .(format(TSS_stats[[4]][[1]], digits = 3))~";"~italic(p) == .(format(TSS_stats[[5]], digits = 3))), size = 4)+
  #annotate("text",x=6,y=0.5,label=bquote(italic(R) == .(format(gene_stats[[4]][[1]], digits = 3))~";"~italic(p) == .(format(gene_stats[[5]], digits = 3))), size = 4)+
  #annotate("text",x=6,y=0.9,label=paste(bquote(italic(R)^2) round(gene_stats[[3]], digits=4),"\n p-val=", round(gene_stats[[5]], digit=4),sep=""), size = 3)+
  #stat_cor(label.x = 3, label.y = 34)+
    #stat_summary(fun = "median", geom="text", size = 2, vjust=-0.5, label=round(median$x, 2))+
    #scale_fill_manual(values=get_brewer_pal("Spectral", n= bins, contrast = c(0.3, 0.8), stretch = F, plot = F))+
    #scale_fill_manual(values=colorRampPalette(c("#F8C6CB","#CB4A42"),  bias=1)(10))+
    #scale_fill_manual(values=get_brewer_pal("BuPu", n= bins, contrast = c(0.3, 0.6), stretch = F, plot = F))+
    #guides(colour = "colorbar", size = "legend", shape = "none")+
    #guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10))+
    theme_classic(base_size=20)+
    theme(legend.position = "none",panel.border = element_blank(),
        plot.title = element_text(color="black",size=11,hjust=0.5),
        axis.title= element_blank(),
        axis.ticks = element_line(linewidth=0.28,color = "black"),
        axis.ticks.length = unit(0.075, "cm"),
        axis.line = element_line(linewidth=0.28,color = "black"),
        axis.text = element_text(color="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_blank(),
        #axis.text.x=element_text(size=8,  angle = 90, hjust=0.75),
        axis.text.x=if (isTRUE(opt$nolab)) element_blank() else element_text(size=9),
        #axis.text.x=element_text(size=9),
        axis.text.y=if (isTRUE(opt$nolab)) element_blank() else element_text(size=11),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
print(plot)
dev.off()

# Generate legend SVGs once (on the non-nolab call only to avoid duplicate writes)
if (!isTRUE(opt$nolab)) {
  get_legend_grob <- function(p) {
    g <- ggplotGrob(p)
    idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
    if (length(idx) == 0L) return(NULL)
    g$grobs[[idx]]
  }

  legend_base <- plot +
    guides(
      colour = guide_legend(override.aes = list(size=c(3,3,3,3), stroke=c(0.5,0.5,0.5,0)))
    ) +
    theme(
      legend.position   = "right",
      legend.text       = element_text(size=10, color="black"),
      legend.title      = element_blank(),
      legend.key        = element_rect(fill="transparent"),
      legend.background = element_rect(fill="transparent", color=NA)
    )

  leg <- get_legend_grob(legend_base)
  if (!is.null(leg)) {
    svg(paste0(out, "/legend.svg"), width=1.2, height=1.0)
    grid::grid.draw(leg)
    dev.off()

    leg_nolab <- get_legend_grob(
      legend_base + theme(legend.text = element_blank())
    )
    svg(paste0(out, "/legend_nolab.svg"), width=0.5, height=1.0)
    grid::grid.draw(leg_nolab)
    dev.off()
  }
}
