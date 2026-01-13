#######################################################
# ANALISI miRNA SEQUENCING - DE, GRAFICI & BIOMARCATORI
# INPUT: counts/*.tsv, metadata.tsv
# OUTPUT: CSV, Volcano plot, Heatmap, MDS, UpSet
# Cartella output: results/
#######################################################

# -------------------------
# Librerie
# -------------------------
library(edgeR)
library(limma)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(UpSetR)

cat("==== INIZIO ANALISI ====\n")

# -------------------------
# 1️⃣ Leggi counts e metadata
# -------------------------
cat("Leggo counts e metadata...\n")
files <- list.files("counts", pattern="tsv$", full.names=TRUE)
if(length(files) == 0){
  stop("ERROR: nessun file counts trovato in counts/")
}

counts <- Reduce(function(x,y) merge(x,y,by="miRNA",all=TRUE),
                 lapply(files, read.delim))
rownames(counts) <- counts$miRNA
counts$miRNA <- NULL
counts[is.na(counts)] <- 0

meta <- read.delim("metadata.tsv")
if(!all(meta$sample %in% colnames(counts))){
  stop("ERROR: nomi dei campioni in metadata.tsv non corrispondono ai file counts")
}
counts <- counts[, meta$sample]
group <- factor(meta$condition)

# -------------------------
# 2️⃣ Creazione cartella results
# -------------------------
dir.create("results", showWarnings = FALSE)
cat("Cartella results/ creata.\n")

# -------------------------
# 3️⃣ MDS plot globale
# -------------------------
cat("Creazione MDS plot...\n")
dge <- DGEList(counts=counts, group=group)
dge <- dge[filterByExpr(dge), , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

png("results/MDS_all_samples.png", width=800, height=600)
plotMDS(dge,
        col = as.numeric(group),
        pch = 16,
        main = "MDS plot – miRNA expression")
legend("topright",
       legend = levels(group),
       col = seq_along(levels(group)),
       pch = 16)
dev.off()

# -------------------------
# 4️⃣ Linear model + voom
# -------------------------
cat("Calcolo DE singole patologie...\n")
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

v <- voomWithQualityWeights(dge, design)
fit <- lmFit(v, design)

contrasts <- makeContrasts(
  DA_vs_CTRL = atopic_dermatitis - control,
  PB_vs_CTRL = pemphigoid - control,
  PSO_vs_CTRL = psoriasis - control,
  levels=design
)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# -------------------------
# 5️⃣ Salva risultati singoli + grafici
# -------------------------
contrast_list <- list(
  "DA_vs_CTRL" = topTable(fit2, coef=1, number=Inf),
  "PB_vs_CTRL" = topTable(fit2, coef=2, number=Inf),
  "PSO_vs_CTRL" = topTable(fit2, coef=3, number=Inf)
)

cat("Salvo CSV e generazione grafici per confronti singoli...\n")
for(name in names(contrast_list)){
  res <- contrast_list[[name]]
  write.csv(res, paste0("results/", name, ".csv"))
  
  # Volcano plot
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = "logFC",
                  y = "adj.P.Val",
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  title = paste0(name, " Volcano"),
                  subtitle = "Differential miRNA expression",
                  caption = "FDR < 0.05")
  ggsave(paste0("results/Volcano_", name, ".png"), width=8, height=6)
  
  # Heatmap top 30 DE miRNA
  top_mirna <- rownames(res[res$adj.P.Val < 0.05, ])
  top_mirna <- head(top_mirna, 30)
  if(length(top_mirna) >= 2){
    expr <- v$E[top_mirna, ]
    annotation_col <- data.frame(Condition = group)
    rownames(annotation_col) <- colnames(expr)
    
    pheatmap(expr,
             scale="row",
             annotation_col=annotation_col,
             show_colnames=FALSE,
             main=paste0("Top DE miRNA – ", name))
    ggsave(paste0("results/Heatmap_top30_", name, ".png"), width=8, height=6)
  }
}

# -------------------------
# 6️⃣ Confronto globale: tutti i malati vs sani
# -------------------------
cat("Calcolo DE globale: tutti i malati vs CTRL...\n")
meta$phenotype <- ifelse(meta$condition=="control", "control", "disease")
meta$phenotype <- factor(meta$phenotype)
group2 <- meta$phenotype

dge2 <- DGEList(counts=counts, group=group2)
dge2 <- dge2[filterByExpr(dge2), , keep.lib.sizes=FALSE]
dge2 <- calcNormFactors(dge2)

design2 <- model.matrix(~0 + group2)
colnames(design2) <- levels(group2)

v2 <- voomWithQualityWeights(dge2, design2)
fit2b <- lmFit(v2, design2)

contrast2 <- makeContrasts(DISEASE_vs_CTRL = disease - control,
                           levels=design2)
fit2b <- contrasts.fit(fit2b, contrast2)
fit2b <- eBayes(fit2b)

res_disease <- topTable(fit2b, number=Inf)
write.csv(res_disease, "results/DISEASE_vs_CTRL.csv")

EnhancedVolcano(res_disease,
                lab = rownames(res_disease),
                x = "logFC",
                y = "adj.P.Val",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "All diseases vs Control",
                subtitle = "Shared pruritus-associated miRNAs")
ggsave("results/Volcano_DISEASE_vs_CTRL.png", width=8, height=6)

top_shared <- rownames(res_disease[res_disease$adj.P.Val < 0.05, ])
top_shared <- head(top_shared, 30)
if(length(top_shared) >= 2){
  expr2 <- v2$E[top_shared, ]
  annotation_col <- data.frame(Phenotype = group2)
  rownames(annotation_col) <- colnames(expr2)
  
  pheatmap(expr2,
           scale="row",
           annotation_col=annotation_col,
           show_colnames=FALSE,
           main="Top shared miRNAs – Disease vs Control")
  ggsave("results/Heatmap_shared_miRNA.png", width=8, height=6)
}

# -------------------------
# 7️⃣ UpSet plot
# -------------------------
cat("Creazione UpSet plot...\n")
sig_DA <- rownames(topTable(fit2, coef=1, number=Inf)[topTable(fit2, coef=1, number=Inf)$adj.P.Val<0.05,])
sig_PB <- rownames(topTable(fit2, coef=2, number=Inf)[topTable(fit2, coef=2, number=Inf)$adj.P.Val<0.05,])
sig_PSO <- rownames(topTable(fit2, coef=3, number=Inf)[topTable(fit2, coef=3, number=Inf)$adj.P.Val<0.05,])

mir_list <- list("DA"=sig_DA, "PB"=sig_PB, "PSO"=sig_PSO)

png("results/UpSet_miRNA.png", width=1000, height=800)
upset(fromList(mir_list),
      order.by="freq",
      main.bar.color="darkblue",
      sets.bar.color="steelblue",
      text.scale=c(2,2,1.5,1.5,1.5,1.5))
dev.off()

# -------------------------
# 8️⃣ Identificazione miRNA candidati biomarcatori di prurito
# -------------------------
cat("Selezione automatica dei biomarcatori di prurito...\n")
sig_global <- rownames(res_disease[res_disease$adj.P.Val < 0.05, ])
sig_shared <- Reduce(intersect, list(sig_DA, sig_PB, sig_PSO))
pruritus_biomarkers <- intersect(sig_global, sig_shared)

biomarker_table <- data.frame(miRNA = pruritus_biomarkers)
write.csv(biomarker_table, "results/miRNA_pruritus_biomarkers.csv", row.names=FALSE)

if(length(pruritus_biomarkers) >= 2){
  expr_biom <- v2$E[pruritus_biomarkers, ]
  annotation_col <- data.frame(Phenotype = group2)
  rownames(annotation_col) <- colnames(expr_biom)
  
  pheatmap(expr_biom,
           scale="row",
           annotation_col=annotation_col,
           show_colnames=FALSE,
           main="Candidate pruritus miRNAs")
  ggsave("results/Heatmap_pruritus_biomarkers.png", width=8, height=6)
}

cat("==== ANALISI COMPLETA! Tutti i file in results/ ====\n")
