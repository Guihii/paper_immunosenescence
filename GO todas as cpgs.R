library(missMethyl)
library(readr)
library(ggplot2)
library(dplyr)

# Carregar seus dados
dados <- read_csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\vo2.csv")

# Selecionar CpGs significativos (FDR < 0.05)
dados_sig <- dados %>% filter(adj.P.Val < 0.05)
cpg.sig <- dados_sig$MarkerName
cpg.all <- dados$MarkerName

# Rodar análise de enriquecimento GO considerando todos os CpGs significativos juntos
go.all <- gometh(sig.cpg = cpg.sig, all.cpg = cpg.all, collection = "GO", plot.bias = FALSE)

# Selecionar top 15 vias com menor p-valor
top_go_all <- go.all %>%
  arrange(P.DE) %>%
  slice(1:15) %>%
  mutate(TERM = factor(TERM, levels = rev(TERM)))  # para ordem no gráfico

# Criar gráfico de barras horizontais
GOplot_all1 <- ggplot(top_go_all, aes(x = -log10(P.DE), y = TERM)) +
  geom_bar(stat = "identity", fill = "#457B9D") +
  theme_minimal() +
  labs(
    title = "Top 15 GO Enrichment – All Significant CpGs VO2",
    x = "-log10(P-value)",
    y = "GO Term"
  )

# Mostrar gráfico
print(GOplot_all1)

# Salvar imagem
ggsave("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\GO_vo2_ALL.jpg",
       plot = GOplot_all1, width = 10, height = 6, dpi = 1200)
