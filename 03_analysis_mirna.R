#######################################################
# ANALISI miRNA QIAseq (UMI-deduplicated)
# Obiettivo: biomarcatori di prurito
# Gruppi: SANO, PB, PSO, DA
#######################################################

########################
# LIBRERIE
########################
library(edgeR)
library(limma)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(UpSetR)

cat("==== INIZIO ANALISI miRNA QIAseq ====\n")

########################
# 1️⃣ LETTURA COUNTS
########################

cat("Lettura file counts...\n")

files <- list.files("counts",
                    pattern = "_counts.tsv$",
                    full.names = TRUE)

if(length(files) == 0){
  stop("ERRORE: nessun file trovato in counts/")
}

counts_list <- lapply(files, function(f){
  df <- read.delim(f, stringsAsFactors = FALSE)
  colnames(df) <- c("miRNA", basename(f))
  df
})

counts <- Reduce(function(x, y) merge(x, y, by = "miRNA", all = TRUE),
                 counts_list)

rownames(counts) <- counts$miRNA
counts$miRNA <- NULL
counts[is.na(counts)] <- 0

########################
# 2️⃣ METADATA
########################

cat("Lettura metadata...\n")

meta <- read.delim("metadata.tsv", stringsAsFactors = FALSE)

if(!all(meta$sample %in% colnames(counts))){
  stop("ERRORE: metadata non corrisponde ai file counts")
}

counts <- counts[, meta$sample]

group <- factor(meta$condition,
                levels = c("SANO", "PB", "PSO", "DA"))

########################
# (OPZIONALE) BATCH
########################
# Se hai batch di sequenziamento / preparazione
# scommenta queste righe e aggiungi la colonna "batch" a metadata.tsv

# batch <- factor(meta$batch)

########################
# 3️⃣ OGGETTO DGE
########################

dge <- DGEList(counts = counts, group = group)

# filtro espressione bassa
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]

# normalizzazione TMM
dge <- calcNormFactors(dge)

########################
# 4️⃣ CARTELLE OUTPUT
########################

dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)

########################
# 5️⃣ PCA & MDS
########################

cat("PCA e MDS...\n")

logCPM <- cpm(dge, log = TRUE, prior.count = 1)

# PCA
pca <- prcomp(t(logCPM), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Condition = group
)

ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA – miRNA expression",
       subtitle = "QIAseq UMI-deduplicated",
       x = "PC1", y = "PC2")

ggsave("results/plots/PCA_all_samples.png", width = 7, height = 6)

# MDS
png("results/plots/MDS_all_samples.png", width = 800, height = 600)
plotMDS(dge, col = as.numeric(group), pch = 16,
        main = "MDS – miRNA expression")
legend("topright", legend = levels(group),
       col = seq_along(levels(group)), pch = 16)
dev.off()

########################
# 6️⃣ MODELLO LINEARE (voom)
########################

cat("Modello lineare...\n")

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# se usi batch:
# design <- model.matrix(~0 + group + batch)

v <- voomWithQualityWeights(dge, design)
fit <- lmFit(v, design)

########################
# 7️⃣ CONTRASTI
########################

contrasts <- makeContrasts(
  PB_vs_SANO  = PB  - SANO,
  PSO_vs_SANO = PSO - SANO,
  DA_vs_SANO  = DA  - SANO,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

########################
# 8️⃣ RISULTATI DE + PLOT
########################

contrast_list <- list(
  PB_vs_SANO  = "PB_vs_SANO",
  PSO_vs_SANO = "PSO_vs_SANO",
  DA_vs_SANO  = "DA_vs_SANO"
)

for(name in names(contrast_list)){
  
  res <- topTable(fit2,
                  coef = contrast_list[[name]],
                  number = Inf,
                  sort.by = "P")
  
  write.csv(res, paste0("results/", name, ".csv"))
  
  # Volcano
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = "logFC",
                  y = "adj.P.Val",
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  title = name,
                  subtitle = "miRNA biomarkers",
                  caption = "FDR < 0.05")
  
  ggsave(paste0("results/plots/Volcano_", name, ".png"),
         width = 8, height = 6)
  
  # Heatmap top 30
  sig <- rownames(res[res$adj.P.Val < 0.05, ])
  sig <- head(sig, 30)
  
  if(length(sig) >= 2){
    expr <- v$E[sig, ]
    annotation_col <- data.frame(Condition = group)
    rownames(annotation_col) <- colnames(expr)
    
    pheatmap(expr,
             scale = "row",
             annotation_col = annotation_col,
             show_colnames = FALSE,
             main = paste("Top miRNA –", name))
    
    ggsave(paste0("results/plots/Heatmap_", name, ".png"),
           width = 8, height = 6)
  }
}

########################
# 9️⃣ FIRMA GLOBALE PRURITO
########################

cat("Analisi globale PRURITUS vs SANO...\n")

meta$phenotype <- ifelse(meta$condition == "SANO",
                         "SANO", "PRURITUS")
group2 <- factor(meta$phenotype,
                 levels = c("SANO", "PRURITUS"))

dge2 <- DGEList(counts = counts, group = group2)
dge2 <- dge2[filterByExpr(dge2), , keep.lib.sizes = FALSE]
dge2 <- calcNormFactors(dge2)

design2 <- model.matrix(~0 + group2)
colnames(design2) <- levels(group2)

v2 <- voomWithQualityWeights(dge2, design2)
fit_pru <- lmFit(v2, design2)

contrast2 <- makeContrasts(
  PRURITUS_vs_SANO = PRURITUS - SANO,
  levels = design2
)

fit_pru <- contrasts.fit(fit_pru, contrast2)
fit_pru <- eBayes(fit_pru)

res_pruritus <- topTable(fit_pru, number = Inf)
write.csv(res_pruritus, "results/PRURITUS_vs_SANO.csv")

EnhancedVolcano(res_pruritus,
                lab = rownames(res_pruritus),
                x = "logFC",
                y = "adj.P.Val",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Pruritus miRNA signature")

ggsave("results/plots/Volcano_PRURITUS_vs_SANO.png",
       width = 8, height = 6)

########################
# 🔟 UPSET PLOT
########################

sig_PB  <- rownames(subset(
  topTable(fit2, coef="PB_vs_SANO", number=Inf),
  adj.P.Val < 0.05))

sig_PSO <- rownames(subset(
  topTable(fit2, coef="PSO_vs_SANO", number=Inf),
  adj.P.Val < 0.05))

sig_DA  <- rownames(subset(
  topTable(fit2, coef="DA_vs_SANO", number=Inf),
  adj.P.Val < 0.05))

png("results/plots/UpSet_miRNA_pruritus.png",
    width = 1000, height = 800)

upset(fromList(list(
  PB  = sig_PB,
  PSO = sig_PSO,
  DA  = sig_DA)),
  order.by = "freq")

dev.off()

cat("==== FINE ANALISI miRNA QIAseq ====\n")
