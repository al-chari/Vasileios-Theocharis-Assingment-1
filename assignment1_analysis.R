#Import expression data
expMat <- read.table("E-MEXP-3715_data.txt", header = TRUE, row.names = 1,
                     as.is = TRUE, sep = "\t")

#Define sample groups
sdrf <- read.table("E-MEXP-3715.sdrf.txt", header = TRUE, as.is = TRUE,
                   row.names = 1, sep = "\t")
group <- sdrf$Characteristics.clinical.information.
sex <- sdrf$Characteristics.sex.
pairs <- sdrf$Characteristics.individual.

#Load package
library(limma)

#Define design matrix
studyDesign = data.frame(Group = group, Sex = sex, Pairs = pairs)
designMatrix <- model.matrix(~ 0 + Group + Sex + Pairs, data = studyDesign)
colnames(designMatrix)[1:3] <- c("Tumor", "Healthy", "Sex")

#Fit linear model
fit <- lmFit(expMat, designMatrix)

#Define contrast matrix
contrastMatrix <- makeContrasts(TumorVsHealthy = Tumor - Healthy,
                                levels = designMatrix)

#Specify contrasts
fit <- contrasts.fit(fit, contrastMatrix)

#Perform statistical testing
fit <- eBayes(fit, robust = TRUE)

#Filter results
results <- topTable(fit, coef = "TumorVsHealthy", number = nrow(expMat),
                    p.value = 0.05, lfc = 1)
resultsAll <- topTable(fit, coef = "TumorVsHealthy", number = nrow(expMat))

#Add gene symbols
library(hugene10sttranscriptcluster.db)
probe2symbol <- mapIds(
  hugene10sttranscriptcluster.db,
  keys = rownames(results),
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first")
results$Symbol <- probe2symbol

#Volcano plot
library(openxlsx)
library(ggplot2)
write.xlsx(results, "Tumor vs healthy DEGs.xlsx", rowNames = TRUE)
resDF <- data.frame(
  logFC = resultsAll$logFC,
  FDR = resultsAll$adj.P.Val
)
volcanoPlot <- ggplot(resDF,
                      aes(x = logFC, y = -log10(FDR),
                          color = ifelse(FDR < 0.05 & abs(logFC) > 1, "darkred", "grey"))) +
  geom_point() +
  xlab(expression("Fold Change, Log"[2]*"")) +
  ylab(expression("Adjusted P value, Log"[10]*"")) +
  geom_vline(
    xintercept = c(-1, 1), linetype = "dotted", size = 1) +
  geom_hline(
    yintercept = -log10(0.05), linetype = "dotted", size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("darkred", "grey", "steelblue"))
print(volcanoPlot)
ggsave(plot = volcanoPlot, filename = "Volcano plot.pdf",
       width = 8, height = 5)
