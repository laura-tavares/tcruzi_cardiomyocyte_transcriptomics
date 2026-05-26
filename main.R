############################################################
# ANÁLISE DE EXPRESSÃO DIFERENCIAL
# CARDIOMIÓCITOS HUMANOS INFECTADOS POR T. cruzi
# GSE223600
# 24hpi vs 0hpi
############################################################

############################################################
# INSTALAÇÃO DE PACOTES
############################################################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("DESeq2")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("clusterProfiler")
# BiocManager::install("enrichplot")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("GEOquery")
#
# install.packages("pheatmap")
# install.packages("ggplot2")
# install.packages("gplots")
# install.packages("ggrepel")
# install.packages("dplyr")

############################################################
# CARREGAR BIBLIOTECAS
############################################################

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(gplots)
library(ggrepel)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(GEOquery)
library(dplyr)

############################################################
# BAIXAR DADOS DO GEO
############################################################

getGEOSuppFiles("GSE223600")

############################################################
# EXTRAIR ARQUIVOS
############################################################

untar(
  "GSE223600/GSE223600_RAW.tar",
  exdir = "GSE223600/raw"
)

############################################################
# VER ARQUIVOS
############################################################

list.files("GSE223600/raw")

############################################################
# DEFINIR ARQUIVOS UTILIZADOS
############################################################

files <- c(
  "GSE223600/raw/GSM6965202_0hpi_1_.counts.txt.gz",
  "GSE223600/raw/GSM6965203_0hpi_2_counts.txt.gz",
  "GSE223600/raw/GSM6965204_0hpi_3_counts.txt.gz",
  "GSE223600/raw/GSM6965205_24hpi_1_counts.txt.gz",
  "GSE223600/raw/GSM6965206_24hpi_2_counts.txt.gz",
  "GSE223600/raw/GSM6965207_24hpi_3_counts.txt.gz"
)

############################################################
# NOMES DAS AMOSTRAS
############################################################

sample_names <- c(
  "Controle_1",
  "Controle_2",
  "Controle_3",
  "Infectado_1",
  "Infectado_2",
  "Infectado_3"
)

############################################################
# LER PRIMEIRO ARQUIVO
############################################################

df <- read.delim(
  gzfile(files[1])
)

############################################################
# VER ESTRUTURA
############################################################

head(df)

############################################################
# CRIAR MATRIZ INICIAL
############################################################

count_matrix <- data.frame(
  ENSEMBL = df$X
)

############################################################
# ADICIONAR CONTAGENS À MATRIZ
############################################################

for(i in 1:length(files)){
  
  temp <- read.delim(
    gzfile(files[i])
  )
  
  count_matrix[[sample_names[i]]] <-
    temp[, grep("_Raw$", colnames(temp))]
}

############################################################
# REMOVER VERSÃO DOS ENSG
############################################################

count_matrix$ENSEMBL <- sub(
  "\\..*",
  "",
  count_matrix$ENSEMBL
)

############################################################
# REMOVER GENES DUPLICADOS
############################################################

count_matrix <- count_matrix[
  !duplicated(count_matrix$ENSEMBL),
]

############################################################
# DEFINIR IDs COMO ROW NAMES
############################################################

rownames(count_matrix) <- count_matrix$ENSEMBL

############################################################
# REMOVER COLUNA ENSEMBL
############################################################

count_matrix <- count_matrix[, -1]

############################################################
# VER MATRIZ FINAL
############################################################

head(count_matrix)

dim(count_matrix)

############################################################
# CRIAR METADADOS
############################################################

conditions <- factor(c(
  rep("Controle", 3),
  rep("Infectado", 3)
))

sampleTable <- data.frame(
  sampleName = colnames(count_matrix),
  condition = conditions
)

rownames(sampleTable) <- colnames(count_matrix)

sampleTable

############################################################
# CRIAR OBJETO DESeq2
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sampleTable,
  design = ~ condition
)

############################################################
# FILTRAR BAIXAS CONTAGENS
############################################################

dds <- dds[
  rowSums(counts(dds)) > 10,
]

############################################################
# ANÁLISE DE EXPRESSÃO DIFERENCIAL
############################################################

dds <- DESeq(dds)

############################################################
# PCA
############################################################

rld <- rlogTransformation(
  dds,
  blind = FALSE
)

plotPCA(
  rld,
  intgroup = "condition",
  ntop = nrow(counts(dds))
) +
  ggtitle(
    expression(
      paste(
        "Cardiomiócitos infectados por ",
        italic("T. cruzi")
      )
    )
  ) +
  labs(
    color = "Grupo"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

############################################################
# HEATMAP DE CORRELAÇÃO
############################################################

cU <- cor(as.matrix(assay(rld)))

cols <- c(
  "dodgerblue3",
  "firebrick3"
)[conditions]

heatmap.2(
  cU,
  symm = TRUE,
  col = colorRampPalette(
    c("white", "darkblue")
  )(100),
  labCol = colnames(cU),
  labRow = colnames(cU),
  distfun = function(c) as.dist(1 - c),
  trace = "none",
  Colv = TRUE,
  cexRow = 1,
  cexCol = 1,
  margins = c(6,6),
  key = FALSE,
  font = 2,
  main = "Correlação entre amostras",
  RowSideColors = cols,
  ColSideColors = cols
)

############################################################
# RESULTADOS DE EXPRESSÃO DIFERENCIAL
############################################################

res <- results(
  dds,
  contrast = c(
    "condition",
    "Infectado",
    "Controle"
  )
)

############################################################
# NORMALIZAÇÃO
############################################################

norm.counts <- counts(
  dds,
  normalized = TRUE
)

############################################################
# CONVERTER PARA DATAFRAME
############################################################

all <- data.frame(
  res,
  norm.counts
)

############################################################
# ADICIONAR ENSEMBL
############################################################

all$ENSEMBL <- rownames(all)

############################################################
# CONVERTER ENSG PARA SYMBOL
############################################################

gene_conversion <- bitr(
  all$ENSEMBL,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

############################################################
# JUNTAR ANOTAÇÕES
############################################################

all <- merge(
  all,
  gene_conversion,
  by = "ENSEMBL",
  all.x = TRUE
)

############################################################
# MANTER ENSG QUANDO NÃO EXISTIR SYMBOL
############################################################

all$SYMBOL[
  is.na(all$SYMBOL)
] <- all$ENSEMBL[
  is.na(all$SYMBOL)
]

############################################################
# EXPORTAR RESULTADOS
############################################################

write.table(
  all,
  file = "Tabelas/DESeq2_all_GSE223600.txt",
  sep = "\t"
)

# write.csv(
#   all,
#   file = "DESeq2_all_GSE223600.csv"
# )

############################################################
# RESUMO
############################################################

summary(
  res,
  alpha = 0.05
)

sum(
  res$padj < 0.05 &
    res$log2FoldChange > 1,
  na.rm = TRUE
)

sum(
  res$padj < 0.05 &
    res$log2FoldChange < -1,
  na.rm = TRUE
)

############################################################
# MA PLOT
############################################################

# plotMA(
#   res,
#   ylim = c(-8,8),
#   alpha = 0.05,
#   main = expression(
#     paste(
#       "Infecção por ",
#       italic("T. cruzi")
#     )
#   )
# )

par(
  cex.main = 1.3,
  cex.lab = 1.2,
  cex.axis = 1.1
)

plotMA(
  res,
  ylim = c(-8,8),
  alpha = 0.05,
  colSig = "red",
  
  main = expression(
    paste(
      "Infecção por ",
      italic("T. cruzi")
    )
  ),
  
  xlab = "Média normalizada de contagens",
  ylab = "Log2 Fold Change"
)

abline(
  h = 0,
  col = "blue",
  lty = 2
)

############################################################
# VOLCANO PLOT
############################################################

all$grupo <- "Não significativo"

all$grupo[
  all$padj < 0.05 &
    all$log2FoldChange > 1
] <- "Superexpresso"

all$grupo[
  all$padj < 0.05 &
    all$log2FoldChange < -1
] <- "Subexpresso"

volcano <- na.omit(all)

top_up <- volcano %>%
  filter(
    padj < 0.05 &
      log2FoldChange > 2
  ) %>%
  arrange(padj) %>%
  head(5)

top_down <- volcano %>%
  filter(
    padj < 0.05 &
      log2FoldChange < -2
  ) %>%
  arrange(padj) %>%
  head(5)

top_genes <- rbind(
  top_up,
  top_down
)

ggplot(
  volcano,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = grupo
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 2.5
  ) +
  
  scale_color_manual(
    values = c(
      "Superexpresso" = "firebrick3",
      "Subexpresso" = "dodgerblue3",
      "Não significativo" = "gray75"
    )
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = SYMBOL),
    size = 4,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  # coord_cartesian(
  #   xlim = c(-5, 10),
  #   ylim = c(0, 60)
  # ) +
  
  labs(
    title = expression(
      paste(
        "Cardiomiócitos infectados por ",
        italic("T. cruzi")
      )
    ),
    x = "Log2 Fold Change",
    y = expression(-log[10](padj)),
    color = NULL
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "top",
    legend.text = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    )
  )


############################################################
# TOP 5 GENES
############################################################

sig <- subset(
  all,
  padj < 0.05
)

top5 <- head(
  sig[
    order(
      abs(sig$log2FoldChange),
      decreasing = TRUE
    ),
    c("ENSEMBL", "log2FoldChange", "SYMBOL")
  ],
  5
)

top5

############################################################
# TOP 20 GENES
############################################################

top20gene <- rownames(res)[
  res$padj <= sort(res$padj)[20] &
    !is.na(res$padj)
]

gene_names <- gene_conversion$SYMBOL[
  match(top20gene, gene_conversion$ENSEMBL)
]

gene_names[
  is.na(gene_names)
] <- top20gene[
  is.na(gene_names)
]

mat20top <- assay(rld)[top20gene, ]

rownames(mat20top) <- gene_names

mat20top <- mat20top[
  !duplicated(rownames(mat20top)),
]

mat20top <- mat20top - rowMeans(mat20top)

anno20top <- as.data.frame(
  colData(rld)[c("condition")]
)

colnames(anno20top) <- "Grupo"

ann_colors <- list(
  Grupo = c(
    Controle = "dodgerblue3",
    Infectado = "firebrick3"
  )
)

pheatmap(
  mat20top,
  annotation_col = anno20top,
  annotation_colors = ann_colors,
  fontsize_row = 8,
  main = "Top 20 genes diferencialmente expressos"
)

############################################################
# ENRIQUECIMENTO GO - BP
############################################################

sig_genes_df <- subset(
  all,
  padj < 0.05
)

genes <- sig_genes_df$log2FoldChange

names(genes) <- sig_genes_df$SYMBOL

genes <- na.omit(genes)

genes <- names(genes)[
  abs(genes) > 1
]

original_gene_list <- all$log2FoldChange

names(original_gene_list) <- all$SYMBOL

gene_list <- na.omit(original_gene_list)

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

go_enrich_bp <- enrichGO(
  gene = genes,
  universe = names(gene_list),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  readable = TRUE,
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

barplot(
  go_enrich_bp,
  drop = TRUE,
  showCategory = 10,
  title = "Processos Biológicos enriquecidos",
  font.size = 8
)

dotplot(
  go_enrich_bp,
  showCategory = 10,
  title = "Processos biológicos enriquecidos"
)

cnetplot(
  go_enrich_bp,
  showCategory = 10,
  foldChange = gene_list
)

############################################################
# ENRIQUECIMENTO GO - CC
############################################################

go_enrich_cc <- enrichGO(
  gene = genes,
  universe = names(gene_list),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  readable = TRUE,
  ont = "CC",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

barplot(
  go_enrich_cc,
  drop = TRUE,
  showCategory = 10,
  title = "GO - Componentes Celulares",
  font.size = 8
)

dotplot(
  go_enrich_cc,
  showCategory = 10,
  title = "Componentes celulares enriquecidos"
)

cnetplot(
  go_enrich_cc,
  showCategory = 10,
  foldChange = gene_list
)

############################################################
# ENRIQUECIMENTO GO - MF
############################################################

go_enrich_mf <- enrichGO(
  gene = genes,
  universe = names(gene_list),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  readable = TRUE,
  ont = "MF",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

barplot(
  go_enrich_mf,
  drop = TRUE,
  showCategory = 10,
  title = "GO - Funções Moleculares",
  font.size = 8
)

dotplot(
  go_enrich_mf,
  showCategory = 10,
  title = "Funções moleculares enriquecidas"
)

cnetplot(
  go_enrich_mf,
  showCategory = 10,
  foldChange = gene_list
)


############################################################
# heatmap genes específicos
############################################################
############################################################
# HEATMAP DE GENES ESPECÍFICOS
############################################################

# Lista de genes de interesse
genes_interesse <- c(
  "PFKFB3",
  "VEGFB",
  "PFKP",
  "ENO1",
  "NOCT",
  "BHLHE40",
  "GPI",
  "SLC243",
  "FOS",
  "ENO2",
  "GAPDH",
  "HK2",
  "PGM1",
  "SLC2A1",
  "GLRX",
  "POF",
  "TPI1",
  "HK1",
  "LOHA",
  "MAFF"
)

############################################################
# CONVERTER SYMBOL -> ENSEMBL
############################################################

genes_heatmap <- gene_conversion$ENSEMBL[
  gene_conversion$SYMBOL %in% genes_interesse
]

############################################################
# EXTRAIR MATRIZ RLOG
############################################################

mat_genes <- SummarizedExperiment::assay(rld)[
  genes_heatmap,
]

############################################################
# TROCAR NOMES ENSG POR SYMBOL
############################################################

rownames(mat_genes) <- gene_map$SYMBOL[
  match(
    rownames(mat_genes),
    gene_map$ENSEMBL
  )
]

############################################################
# REMOVER DUPLICADOS
############################################################

mat_genes <- mat_genes[
  !duplicated(rownames(mat_genes)),
]

############################################################
# ESCALAR POR GENE (Z-SCORE)
############################################################

mat_genes <- t(scale(t(mat_genes)))

############################################################
# ANOTAÇÃO DAS COLUNAS
############################################################

annotation_col <- data.frame(
  Grupo = conditions
)

rownames(annotation_col) <- colnames(mat_genes)

############################################################
# CORES DOS GRUPOS
############################################################

ann_colors <- list(
  Grupo = c(
    Controle = "dodgerblue3",
    Infectado = "firebrick3"
  )
)

############################################################
# PLOTAR HEATMAP
############################################################

pheatmap(
  mat_genes,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  fontsize_row = 10,
  fontsize_col = 10,
  
  main = "Hipóxia"
)

























