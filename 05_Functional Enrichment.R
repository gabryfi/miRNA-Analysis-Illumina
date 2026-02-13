###############################################
# Functional Enrichment dei miRNA DE
###############################################

# Librerie richieste
library(multiMiR)         # per target validati/predetti
library(clusterProfiler)  # per GO/KEGG
library(org.Hs.eg.db)     # annotazione umana
library(ReactomePA)       # per Reactome pathways
library(enrichR)          # alternative database

cat("==== INIZIO ENRICHMENT FUNZIONALE ====\n")

##########################################
# 1️⃣ Estrai miRNA differenzialmente espressi
##########################################

# Lista di miRNA significativi (ad esempio PRURITUS vs SANO)
sig_mir <- rownames(
  subset(res_pruritus,
         adj.P.Val < 0.05 & abs(logFC) > 1)
)

cat("Numero miRNA significativi:", length(sig_mir), "\n")

##########################################
# 2️⃣ Target VALIDATI sperimentalmente
##########################################

# Ottieni i target GENICI per i miRNA
targets <- get_multimir(
  mirna = sig_mir,
  table = "validated",  # richiede evidenza sperimentale
  summary = TRUE
)

# Estrai simboli dei geni
target_genes <- unique(targets@data$target_symbol)

cat("Numero target validati:", length(target_genes), "\n")

##########################################
# 3️⃣ ENRICHMENT GO (Biological Process)
##########################################

ego <- enrichGO(
  gene          = target_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Salva e disegna
write.csv(as.data.frame(ego),
          "results/GO_enrichment_pruritus.csv")

dotplot(ego, showCategory = 25) +
  ggtitle("GO enrichment – BP (target genes)")

ggsave("results/plots/GO_enrichment_pruritus.png",
       width = 10, height = 8)

##########################################
# 4️⃣ ENRICHMENT KEGG
##########################################

# Converti in ENTREZID per KEGG
gene_df <- bitr(target_genes,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

ekegg <- enrichKEGG(
  gene       = gene_df$ENTREZID,
  organism   = "hsa",
  pAdjustMethod = "BH"
)

write.csv(as.data.frame(ekegg),
          "results/KEGG_enrichment_pruritus.csv")

dotplot(ekegg, showCategory = 25) +
  ggtitle("KEGG pathway enrichment")

ggsave("results/plots/KEGG_enrichment_pruritus.png",
       width = 10, height = 8)

##########################################
# 5️⃣ ENRICHMENT REACTOME
##########################################

ereact <- enrichPathway(
  gene          = gene_df$ENTREZID,
  organism      = "human",
  pAdjustMethod = "BH"
)

write.csv(as.data.frame(ereact),
          "results/Reactome_enrichment_pruritus.csv")

dotplot(ereact, showCategory = 25) +
  ggtitle("Reactome pathway enrichment")

ggsave("results/plots/Reactome_enrichment_pruritus.png",
       width = 10, height = 8)

##########################################
# 6️⃣ VALIDAZIONE ALTERNATIVA con enrichR
##########################################

# Scegli database multipli
dbs <- c("GO_Biological_Process_2021",
         "KEGG_2021_Human",
         "Reactome_2022")

enrich_results <- enrichr(target_genes, dbs)

# Salva GO da Enrichr
write.csv(enrich_results$GO_Biological_Process_2021,
          "results/enrichR_GO.csv")

cat("OK: enrichR completato\n")

##########################################
# 7️⃣ Filtraggio pathway infiammatori
##########################################

# Esempio di filtro testo
inflam_go <- subset(as.data.frame(ego),
                    grepl("inflamm|cytokine|immune|NF|interleukin|JAK|STAT",
                          Description,
                          ignore.case = TRUE))

write.csv(inflam_go,
          "results/GO_inflammatory_terms.csv")

cat("==== FINE ENRICHMENT FUNZIONALE ====\n")

