library('edgeR')

#Set working directory where output will go
working_dir = "C:/Users/fangyu/Desktop/final/"
setwd(working_dir)

#Read in gene mapping
mapping=read.table("ENSG_ID2Name.txt", header=FALSE, stringsAsFactors=FALSE, row.names=1)
head(mapping)

#Read in count matrix
rawdata=read.table("gene_read_counts_table2_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
head(rawdata)

#Check dimensions
dim(rawdata)

#Require at least 25% of samples to have count > 5
quant <- apply(rawdata,1,quantile,0.75)
keep <- which((quant >= 5) == 1)
rawdata <- rawdata[keep,]
dim(rawdata)

#Make class labels
class <- factor( c( rep("KO",3), rep("rescue",3) ))

#Get common gene names
genes=rownames(rawdata)
gene_names=mapping[genes,1]


# Make DGEList object
y <- DGEList(counts=rawdata, genes=genes, group=class)
nrow(y)

y <- calcNormFactors(y)


barplot(y$samples$lib.size*1e-6, names=1:6, ylab="Library size (millions)")

plotMDS(y)

# Estimate dispersion
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

# Differential expression test
et <- exactTest(y)

# Print top genes
topTags(et, n = 10, adjust.method = "fdr", sort.by = "PValue", p.value = 1)

#Transform topTags object in data frame
x=as.data.frame(topTags(et, n = nrow(et), adjust.method = "fdr", sort.by = "PValue", p.value = 1))

# Print number of up/down significant genes at FDR = 0.05  significance level
summary(de <- decideTestsDGE(et, p=.05))
detags <- rownames(y)[as.logical(de)]

#Plot log-fold change against log-counts per million, with DE genes highlighted:
plotMD(et)
abline(h=c(-1, 1), col="blue") #The blue lines indicate 2-fold changes

#Plot heatmap
heatmap(heatmap(y$pseudo.counts, margins = c(10,5)))

#Plot and save a heatmap
png(filename='heatmap.png', width=800, height=900)
heatmap(y$pseudo.counts)
graphics.off()

#Make a basic volcano plot
with(x,plot(logFC,-log10(PValue),pch=20,main="Volcanoplot"))

#Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(x,FDR<.05),points(logFC,-log10(PValue),pch=20,col="red"))
with(subset(x,abs(logFC)>1),points(logFC,-log10(PValue),pch=20,col="orange"))
with(subset(x,FDR<.05 & abs(logFC)>1),points(logFC,-log10(PValue), pch=20,col="green"))

# Matrix of significantly DE genes
mat <- cbind(
  genes,gene_names,
  sprintf('%0.3f',log10(et$table$PValue)),
  sprintf('%0.3f',et$table$logFC)
)[as.logical(de),]

colnames(mat) <- c("Gene", "Gene_Name", "Log10_Pvalue", "Log_fold_change")

# Order by log fold change
o <- order(et$table$logFC[as.logical(de)],decreasing=TRUE)
mat <- mat[o,]

# Save table
write.table(mat, file="DE_genes.txt", quote=FALSE, row.names=FALSE, sep="\t")
