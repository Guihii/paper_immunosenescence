if (!requireNamespace("missMethyl", quietly = TRUE)) BiocManager::install("missMethyl")
if (!requireNamespace("limma", quietly = TRUE)) install.packages("limma")

# Carregar os pacotes
library(missMethyl)
library(readr)
library(ggplot2)
library(dplyr)

# Carregar seus dados
dados <- read_csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\sg.csv")

# Separar CpGs significativos (FDR < 0.05)
dados_sig <- dados %>% filter(adj.P.Val < 0.05)

# Separar hiper e hipo com base no logFC
cpg.hyper <- dados_sig %>% filter(logFC > 0) %>% pull(MarkerName)
cpg.hypo  <- dados_sig %>% filter(logFC < 0) %>% pull(MarkerName)
cpg.all <- dados$MarkerName

# Rodar enriquecimento para cada conjunto
go.hyper <- gometh(sig.cpg = cpg.hyper, all.cpg = cpg.all, collection = "GO", plot.bias = FALSE)
go.hypo  <- gometh(sig.cpg = cpg.hypo,  all.cpg = cpg.all, collection = "GO", plot.bias = FALSE)

# Selecionar top 7 de cada
top_go_hyper <- go.hyper %>%
  slice_min(P.DE, n = 7) %>%
  mutate(Direction = "Hypermethylated")

top_go_hypo <- go.hypo %>%
  slice_min(P.DE, n = 7) %>%
  mutate(Direction = "Hypomethylated")

# Combinar e ajustar TERM
top_go <- bind_rows(top_go_hyper, top_go_hypo) %>%
  mutate(TERM = factor(TERM, levels = rev(unique(TERM))))

# Gráfico combinado
PrevsPostGO <- ggplot(top_go, aes(x = -log10(P.DE), y = TERM, size = N, color = Direction)) +
  geom_point() +
  scale_color_manual(values = c("Hypermethylated" = "#1D3557", "Hypomethylated" = "#E63946")) +
  theme_minimal() +
  labs(
    title = "GO Enrichment - PostM Pre vs Post",
    x = "-log10(P-value)",
    y = "GO Term",
    size = "Num CpGs",
    color = "Direction"
  )

print(PrevsPostGO)

ggsave("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\GO_Enrichment_PostM_PrevsPostGO.jpg", plot = PrevsPostGO, width = 10, height = 6, dpi = 1200)
