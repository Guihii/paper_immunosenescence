library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# 1. Ler o arquivo
dados <- read_csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\sg.csv")  # Substitua pelo caminho do seu arquivo

# 2. Carregar anotação
annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotation_df <- as.data.frame(annotation)

# 3. Juntar anotação com os dados
dados_annot <- dados %>%
  left_join(annotation_df %>%
              select(MarkerName = Name,
                     IslandContext = Relation_to_Island,
                     GeneContext = UCSC_RefGene_Group),
            by = "MarkerName")

# 4. Filtrar e preparar dados significativos
dados_sig <- dados_annot %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(Direction = ifelse(logFC > 0, "Hypermethylated", "Hypomethylated"))

# 5. Expandir múltiplas categorias
# Parte A: gene context
gene_data <- dados_sig %>%
  separate_rows(GeneContext, sep = ";") %>%
  filter(GeneContext != "") %>%
  group_by(Context = GeneContext, Direction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Type = "Gene context")

# Parte B: CpG island context
island_data <- dados_sig %>%
  filter(!is.na(IslandContext)) %>%
  group_by(Context = IslandContext, Direction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Type = "CpG island context")

# 6. Combinar os dois
dados_comb <- bind_rows(gene_data, island_data)

# 7. Plotar
ggplot(dados_comb, aes(x = reorder(Context, Count), y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Type ~ ., scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = c("Hypermethylated" = "#2C4A77", "Hypomethylated" = "#E63946")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Genomic and CpG Context of Significant Sites Time",
    x = "Genomic Feature",
    y = "Count",
    fill = "Direction"
  )


ggsave("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\sgGenomic.jpg", plot = last_plot(), width = 8, height = 6, dpi = 1200)

facet_wrap(~Type)