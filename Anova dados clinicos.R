# Instalar pacotes (caso ainda não tenha)
install.packages("readxl")
install.packages("tidyverse")
install.packages("afex")
install.packages("emmeans")

# Carregar pacotes
library(readxl)
library(tidyverse)
library(afex)
library(emmeans)

dados_largos <- read_excel("C:\\Users\\guihi\\Downloads\\VO2xlsx.xlsx")

# 2. Separar por grupo
n <- nrow(dados_largos)

# Criar ID para cada linha (duplicado para cada condição: PreM e PostM)
dados_largos$ID_Part <- 1:n

# 3. Transformar para formato longo com ID exclusivo por grupo
dados_long <- dados_largos %>%
  pivot_longer(
    cols = c(PreM_Pre, PreM_Post, PostM_Pre, PostM_Post),
    names_to = "Condicao",
    values_to = "VO2pico" #mudar para o nome da variavel que estou usando
  ) %>%
  mutate(
    Grupo = case_when(
      str_detect(Condicao, "PreM") ~ "PreM",
      str_detect(Condicao, "PostM") ~ "PostM"
    ),
    Tempo = case_when(
      str_detect(Condicao, "Pre$") ~ "Pre",
      str_detect(Condicao, "Post$") ~ "Post"
    ),
    # ID agora inclui o grupo para ser exclusivo entre participantes
    ID = paste(ID_Part, Grupo, sep = "_")
  ) %>%
  select(ID, Grupo, Tempo, VO2pico)

# 4. Rodar ANOVA de duas vias com medidas repetidas
anova_resultado <- aov_ez(
  id = "ID",
  dv = "VO2pico",#mudar para o nome da variavel que estou usando
  data = dados_long,
  within = "Tempo",
  between = "Grupo",
  type = 3
)

# 5. Resultado da ANOVA
print(anova_resultado)

# 6. Post hoc
library(emmeans)
posthoc <- emmeans(anova_resultado, ~ Grupo * Tempo)
#ajustando por FDR
summary(contrast(posthoc, interaction = "pairwise"), infer = TRUE, adjust = "fdr")
