# Carregar os quatro arquivos CSV
arquivo1CamilaPre <- read.csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\Normalização dos dados ChAMP\\pré menopausa\\Dados Camila\\Decomposição celular\\ProporcoesCelulares_Pre - Camila.csv")
arquivo2CamilaPost <- read.csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\Normalização dos dados ChAMP\\pré menopausa\\Dados Camila\\Decomposição celular\\ProporcoesCelulares_Post - Camila.csv")
arquivo3GuilhermePre <- read.csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\Normalização dos dados ChAMP\\pos menopausa\\Dados Guilherme\\Decomposição celular\\ProporcoesCelulares_Pre_CORRETO.csv")
arquivo4GuilhermePost <- read.csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\Normalização dos dados ChAMP\\pos menopausa\\Dados Guilherme\\Decomposição celular\\ProporcoesCelulares_Post.csv")

# Escolher as colunas que quer comparar (mude "Nome_da_Coluna" para o nome real)
col1 <- arquivo1CamilaPre$Neutro
col2 <- arquivo2CamilaPost$Neutro
col3 <- arquivo3GuilhermePre$Neutro
col4 <- arquivo4GuilhermePost$Neutro

# Função para fazer o teste t e calcular o effect size
fazer_t_test <- function(grupo1, grupo2, nome_grupo1, nome_grupo2) {
  # Teste F para comparar variâncias
  teste_variancia <- var.test(grupo1, grupo2)
  
  # Escolher o tipo de t-test baseado no F-test
  if (teste_variancia$p.value < 0.05) {
    resultado_ttest <- t.test(grupo1, grupo2, var.equal = FALSE)
  } else {
    resultado_ttest <- t.test(grupo1, grupo2, var.equal = TRUE)
  }
  
  # Calcular o Effect Size (Cohen's d)
  n1 <- length(grupo1)
  n2 <- length(grupo2)
  mean1 <- mean(grupo1, na.rm = TRUE)
  mean2 <- mean(grupo2, na.rm = TRUE)
  sd1 <- sd(grupo1, na.rm = TRUE)
  sd2 <- sd(grupo2, na.rm = TRUE)
  
  s_pooled <- sqrt(((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / (n1 + n2 - 2))
  cohen_d <- (mean1 - mean2) / s_pooled
  
  # Mostrar resultados
  cat("\nComparação:", nome_grupo1, "vs", nome_grupo2, "\n")
  print(resultado_ttest)
  cat("Effect size (Cohen's d) =", round(cohen_d, 3), "\n")
}

# Fazer as comparações 
fazer_t_test(col1, col2, "CamilaPre", "CamilaPost")
fazer_t_test(col1, col3, "CamilaPre", "GuilhermePre")
fazer_t_test(col1, col4, "CamilaPre", "GuilhermePost")
fazer_t_test(col2, col3, "CamilaPost", "GuilhermePre")
fazer_t_test(col2, col4, "CamilaPost", "GuilhermePost")
fazer_t_test(col3, col4, "GuilhermePre", "GuilhermePost")



