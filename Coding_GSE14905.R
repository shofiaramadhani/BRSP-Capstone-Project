# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE14905", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("00000000000000000000011111111111111111111111111111",
               "11111111111111111111111111111111")
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Healthy","Disease"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC",
                          "GenBank.Accession","Platform_SPOTID",
                          "Gene.symbol","Gene.title"))

write.table(tT, file=stdout(), row.names=FALSE, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
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

title <- "Boxplot Expression Value Distribution per Sample"

boxplot(ex[,ord], boxwex=0.6, notch=TRUE, main=title,
        outline=FALSE, las=2, col=gs[ord])

legend("topright", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE14905", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE14905")

# DENSITY PLOT
# Visualizes the global distribution of gene expression
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Gene Expression Value Distribution",
    x = "Log2 Expression Value",
    y = "Density"
  )

# UMAP (DIMENSIONALITY REDUCTION)
# See how samples cluster globally based on thousands of genes
umap_input <- t(ex) # Transpose: Samples must be rows
umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot: Sample Clustering",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  )

# VOLCANO PLOT
# Prepare the data (Using a separate object to keep tT clean)
volcano_data <- data.frame(
  logFC = tT$logFC,
  adj.P.Val = tT$adj.P.Val,
  Gene = tT$Gene.symbol
)

# Visualization
volcano_data$status <- "Not Significant"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "Upregulated"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "Downregulated"

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.5) +
  # Custom colors
  scale_color_manual(values = c("Downregulated" = "blue", 
                                "Not Significant" = "grey", 
                                "Upregulated" = "red")) +
  # Vertical helper lines (logFC cut-offs)
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.4) +
  # Horizontal helper line (Significance cut-off)
  # Note: If your data starts at 9, this line might be "below the floor" 
  # and won't be seen, which is what makes it look like Image 2.
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.4) +
  # THE TRICK: Remove the 'expand' command to let ggplot fit the data naturally
  theme_minimal() +
  labs(
    title = "Volcano Plot: GSE14905",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-Value)"
  )

# HEATMAP (TOP 50 DEGs)
# Sorting tT by significance
tT_sorted <- tT[order(tT$adj.P.Val), ]
top50 <- head(tT_sorted, 50)

# Extract expression matrix for these 50 genes using the Probe ID (ID column)
mat_heatmap <- ex[top50$ID, ]

# Create labels: Use Gene Symbol, fallback to ID if empty
gene_labels <- ifelse(
  is.na(top50$Gene.symbol) | top50$Gene.symbol == "",
  top50$ID,
  top50$Gene.symbol
)
rownames(mat_heatmap) <- gene_labels

# Data Cleaning for Heatmap
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ] # Remove NAs
gene_var <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_var > 0, ] # Remove zero variance

# Setup column annotation
annotation_col <- data.frame(Group = gset$group)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row", # Z-score normalization per gene
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  main = "Heatmap: Top 50 Differentially Expressed Genes"
)

# SAVING RESULTS
write.csv(tT, "GSE14905_DEG_Full_Results.csv", row.names = FALSE)
message("Analysis and Visualizations Complete. Results saved to CSV.")

# ENRICHMENT ANALYSIS
# DATA PREPARATION
# Filter significant genes based on Adjusted P-Value < 0.05
# Using the 'tT' variable from your previous DEG analysis
sig_genes <- subset(tT, adj.P.Val < 0.05)

# Convert Gene Symbols to Entrez IDs (Required for KEGG)
# Note: Using "Gene.symbol" column from your specific dataset
gene_conversion <- bitr(sig_genes$Gene.symbol, 
                        fromType = "SYMBOL", 
                        toType   = "ENTREZID", 
                        OrgDb    = org.Hs.eg.db)

# Store the unique Entrez IDs in a variable
entrez_ids <- gene_conversion$ENTREZID

# GENE ONTOLOGY (GO) ENRICHMENT ANALYSIS
# Analyzing all three categories: Biological Process (BP), 
# Cellular Component (CC), and Molecular Function (MF)
ego_results <- enrichGO(gene          = entrez_ids,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        readable      = TRUE)

# Visualization: Faceted Dotplot for GO
dotplot(ego_results, showCategory = 10, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(strip.text = element_text(face = "bold", size = 9)) +
  labs(title = "Gene Ontology Enrichment Analysis",
       subtitle = "Dataset: GSE14905 (Healthy vs Disease)")

# KEGG PATHWAY ENRICHMENT ANALYSIS
# Ensure you have an active internet connection for this step
kegg_results <- enrichKEGG(gene         = entrez_ids,
                           organism     = 'hsa', 
                           pvalueCutoff = 0.05)

# Visualization: Barplot for KEGG Pathways
barplot(kegg_results, showCategory = 15) + 
  labs(title = "KEGG Pathway Enrichment Analysis",
       subtitle = "Dataset: GSE14905")

# EXPORTING RESULTS TO CSV
write.csv(as.data.frame(ego_results), "Enrichment_Results_GO.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_results), "Enrichment_Results_KEGG.csv", row.names = FALSE)