
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO
pid = "GSE3483"
gset <- getGEO("GSE3483", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL339", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "000011112222"
sml <- strsplit(gsms, split="")[[1]]

ex <- exprs(gset)
ex[which(ex <= 0)] <- NaN

write.csv(ex, paste0("data/", pid, "_rawCounts.csv" ), quote = F)
exprs(gset) <- log2(ex) # log2 transform

exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data
aquantile_normalizedCounts = exprs(gset)
write.csv(aquantile_normalizedCounts, paste0("data/", pid, "_aquantile_NormalizedCounts.csv" ), quote = F)

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("QSC","ASC","NMC"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

nall <- nrow(gset)
gset <- gset[complete.cases(exprs(gset)), ]

# calculate precision weights and show plot of mean-variance trend
v <- vooma(gset, design, plot=T)
# OR weights by group
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # attach gene annotations

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""), paste(groups[1],"-",groups[3],sep=""), paste(groups[2],"-",groups[3],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title","Gene.ID","Chromosome.location","Chromosome.annotation","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession","Platform_CLONEID","Platform_ORF","Platform_SPOTID","GO.Function","GO.Process","GO.Component","GO.Function.ID","GO.Process.ID","GO.Component.ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf, coef = NULL)


hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE3483", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE3483", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 5, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=5", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)


write.table(tT2, paste0("data/", pid, "_DE_Results.tsv" ), quote = F, sep = "\t", col.names = NA)

tT3 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf, coef = 1)
tT4 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf, coef = 2)
tT5 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf, coef = 3)

rmC = c("UniGene.title",
  "UniGene.symbol",
  "UniGene.ID",
  "Nucleotide.Title",
  "GI",
  "GenBank.Accession",
  "Platform_CLONEID",
  "Platform_ORF",
  "Platform_SPOTID",
  "GO.Function.ID",
  "GO.Process.ID",
  "GO.Component.ID")



write.table(tT3, paste0("data/", pid, "_", colnames(fit2$coefficients)[1], "_DE_Results.tsv" ), quote = F, sep = "\t", col.names = NA)
write.table(tT4, paste0("data/", pid, "_", colnames(fit2$coefficients)[2], "_DE_Results.tsv" ), quote = F, sep = "\t", col.names = NA)
write.table(tT5, paste0("data/", pid, "_", colnames(fit2$coefficients)[3], "_DE_Results.tsv" ), quote = F, sep = "\t", col.names = NA)





