# Instalar pacotes se necessário
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")

# Carregar bibliotecas
library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(circlize)
library(RColorBrewer)

# 1. Carregar matriz de metilação
beta.m <- read.csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\Normalização dos dados ChAMP\\pré menopausa\\Dados Camila\\dadosnormalizados_CAMILA.csv", 
                   sep = ",", 
                   stringsAsFactors = FALSE, 
                   check.names = FALSE)
rownames(beta.m) <- beta.m[, 1]
beta.m <- beta.m[, -1]
beta.m <- as.matrix(beta.m)
rownames(beta.m) <- trimws(rownames(beta.m))

# 2. Carregar resultado do limma
res_limma <- fread("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\vo2.csv")
colnames(res_limma)[7] <- "CpG"
res_limma$CpG <- trimws(res_limma$CpG)

# 3. Filtrar CpGs com FDR < 0.05
cpgs_signif <- res_limma %>% filter(adj.P.Val < 0.05) %>% pull(CpG)
cat("Total de CpGs significativos:", length(cpgs_signif), "\n")

# 4. Filtrar matriz de beta com os CpGs que realmente estão na matriz
cpgs_comuns <- intersect(cpgs_signif, rownames(beta.m))
cat("CpGs presentes na matriz beta.m:", length(cpgs_comuns), "\n")
beta.f <- beta.m[cpgs_comuns, ]

# 5. Ordenar amostras para pré e pós
amostras_pre <- grep("_pre$", colnames(beta.f), value = TRUE)
amostras_post <- grep("_post$", colnames(beta.f), value = TRUE)
amostras_ordenadas <- c(amostras_pre, amostras_post)
beta.ord <- beta.f[, amostras_ordenadas]

# ⚠️ Converter valores para numérico (evita erro no Heatmap)
beta.ord <- as.data.frame(beta.f[, amostras_ordenadas])
beta.ord <- beta.ord %>% mutate_all(as.numeric)
beta.ord <- as.matrix(beta.ord)
rownames(beta.ord) <- rownames(beta.f)
colnames(beta.ord) <- amostras_ordenadas


# 6. Anotação por tempo
grupo_tempo <- ifelse(grepl("_pre$", amostras_ordenadas), "Pré", "Pós")
ha <- HeatmapAnnotation(
  Timepoint = grupo_tempo, 
  col = list(Timepoint = c("Pré" = "#7a92ac", "Pós" = "#254e70")),
  annotation_name_side = "left",
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9)
  )
)

# Paleta de cores quente (mantida como no seu código)
col_fun <- colorRamp2(c(0, 0.5, 1), c("#f1a340", "white", "#998ec3"))

# Criar o objeto do Heatmap
heatmap_plot <- Heatmap(
  beta.ord,
  name = "Beta",
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  col = col_fun,
  column_title = "Heatmap Time PreM – CpGs FDR < VO2",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "Beta",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(3.5, "cm")
  )
)

# Exportar com proporção apropriada
jpeg("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\HeatmapPreM_Metilacao_VO2_FDR005_FINAL.jpg",
     width = 4000, height = 5000, res = 1200)

draw(
  heatmap_plot,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

width = 4000, height = 5000, res = 1200)

width = 1600, height = 1800, res = 300
