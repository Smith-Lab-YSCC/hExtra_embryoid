# Panel F left
for i in ../data/human/bedgraph/*.bed; do
    n=`basename ${i}`;
    echo ${n};
    bedtools intersect -wa -wb -a CGIs_hg19.bed -b ${i} | bedtools groupby -g 1,2,3 -c 8 -o mean >CGI_${n};
done

```R
require(vioplot)

# load data
track_paths <- list.files(".", "CGI_.*\\.bed", full.names = TRUE)
tracks <- lapply(track_paths, function(x){read.table(x, header=F, col.names=c('chr','start','end','mean'))})
names(tracks) <- sapply(track_paths, function(x) strsplit(basename(x), ".bed")[[1]][1])

tracks <- lapply(names(tracks), function(x){
    data <- tracks[[x]]
    colnames(data) <- c('chr','start','end', x)
    return(data)
})
names(tracks) <- sapply(track_paths, function(x) strsplit(basename(x), ".bed")[[1]][1])
data <- Reduce(function(x,y) merge(x = x, y = y, all = TRUE), tracks)

dim(data)
[1] 26638     8
dim(na.omit(data))
[1] 25871     8
data <- na.omit(data)

# annotate data
CGI <- subset(data, CGI_placenta <= CGI_hESC+0.1)[,-c(1:3)]
hyperCGI <- subset(data, CGI_placenta > CGI_hESC+0.1)[,-c(1:3)]

dim(CGI)
[1] 22629     5
dim(hyperCGI)
[1] 3242    5

pdf('CGI_placenta_hyperCGI_violin_n22629_n3242.pdf', width=21)
vioplot(CGI, side = "left", col='grey', plotCentre = "line")
vioplot(hyperCGI, side = "right", add = T, col='cornflowerblue', plotCentre = "line")
dev.off()
```


# Panel right
for i in ../data/human/bedgraph/*.bed; do
        n=`basename $i`;
        echo $n;
        bedtools intersect -wa -wb -a hg19_tiles_1kb.bed -b $i | bedtools groupby -g 1,2,3,4 -c 8 -o mean >1kb_tiles_${n};
done

```R
require(vioplot)

# load data
track_paths <- list.files(".", "1kb_tiles_.*\\.bed", full.names = TRUE)
tracks <- lapply(track_paths, function(x){read.table(x, header=F, col.names=c('chr','start','end','tile','mean'))})
names(tracks) <- sapply(track_paths, function(x) strsplit(basename(x), ".bed")[[1]][1])

tracks <- lapply(names(tracks), function(x){
    data <- tracks[[x]]
    colnames(data) <- c('chr','start','end', 'tile', x)
    return(data)
})
names(tracks) <- sapply(track_paths, function(x) strsplit(basename(x), ".bed")[[1]][1])

data <- Reduce(function(x,y) merge(x = x, y = y, all = TRUE), tracks)
dim(data)
[1] 2662818       9
dim(na.omit(data))
[1] 2494542       9dim(data)
data <- na.omit(data)

# vioplots
pdf('1kb_tile_vioplots_n2494542.pdf', width=14)
vioplot(data[,-c(1:4)])
dev.off()
```


# Panel G
```R
library(ComplexHeatmap)
library(circlize)

# load data
track_paths <- list.files("../data/human/bedgraph", "\\.bed", full.names = TRUE)
tracks <- lapply(track_paths, function(x){read.table(x, header=F, col.names=c('chr','start','end','mean'))})
names(tracks) <- sapply(track_paths, function(x) strsplit(basename(x), ".bed")[[1]][1])

tracks <- lapply(names(tracks), function(x){
    data <- tracks[[x]]
    colnames(data) <- c('chr','start','end', x)
    return(data)
})
names(tracks) <- sapply(track_paths, function(x) strsplit(basename(x), ".bed")[[1]][1])

data <- Reduce(function(x,y) merge(x = x, y = y, all = TRUE), tracks)

# correlation between samples
cor.mat <- cor(data[,-c(1:3)], use="pairwise.complete.obs")

pdf('WGBS_samples_pairwise_correlation_CpGs_n26674525.pdf')
Heatmap(cor.mat, cluster_rows=hclust(as.dist(1-cor.mat)), cluster_columns=hclust(as.dist(1-cor.mat)))
dev.off()
```
