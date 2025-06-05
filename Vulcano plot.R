# Pacotes necessários
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

# === 1. Ler o arquivo com todos os CpGs ===
dados <- read_csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\sg.csv") %>%
  as.data.frame()

# === 2. Classificar para cor (mesmo logFC pequeno será considerado) ===
dados <- dados %>%
  mutate(
    grupo = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

# === 3. Selecionar top 5 hipo e hiper CpGs com menor FDR ===
top_hipo <- dados %>%
  filter(FDR < 0.05 & logFC < 0) %>%
  arrange(FDR) %>%
  dplyr::slice(1:5)

top_hiper <- dados %>%
  filter(FDR < 0.05 & logFC > 0) %>%
  arrange(FDR) %>%
  dplyr::slice(1:5)

top10 <- bind_rows(top_hipo, top_hiper)

# === 4. Volcano plot com logFC vs -log10(FDR) ===
ggplot(dados, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = grupo), size = 0.8, alpha = 0.6, shape = 16, show.legend = TRUE) +
  scale_color_manual(
    values = c("Up" = "firebrick", "Down" = "cornflowerblue", "Not Sig" = "grey80"),
    name = "Regulation"
  ) +
  geom_text_repel(
    data = top10,
    aes(label = MarkerName, color = grupo),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = NA,
    show.legend = FALSE  # impede letra na legenda
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  ) +
  labs(
    title = "PostM: Pre vs Post",
    x = "logFC",
    y = "-log10(P.Value)"
  )


ggsave("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\sg.jpg", width = 6.5, height = 6.5, dpi = 1200)
