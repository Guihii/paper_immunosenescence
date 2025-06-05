"CÓDIGO BEPE GUILHERME MONASH"

#Carregando o pacote
library("ChAMP")
library("sva")



#Subindo a sample sheet e os chips
testDir <- "C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\Area de trabalho\\Camila"



#Realizando os filtros das pobres
myLoad <- champ.load(testDir,arraytype="EPIC")


champ.QC() 

#Normalizando os dados
myNorm <- champ.norm(beta=myLoad$beta, arraytype="EPIC") #defaut é 450k
champ.QC(myNorm) 


View(myNorm)


# Corrigir os efeitos de batch usando ComBat (corrigindo Slide e Array)
pheno = myLoad$pd

# Corrigir os efeitos de batch para o Slide
myCombat <- ComBat(dat = myNorm, batch = pheno$Slide, mod = NULL)

# Corrigir os efeitos de batch para o Array
myCombat <- ComBat(dat = myCombat, batch = pheno$Array, mod = NULL)

# Verificar a correção dos dados
champ.QC(beta = myCombat, pheno = pheno$Sample_Group, dendrogram = FALSE)


# Rodar o SVD após a correção de batch
champ.SVD(beta = as.data.frame(myCombat), pd = pheno)




#comandos para fazer o SVD rodar
champ.SVD(beta=myNorm %>% as.data.frame(), pd=myLoad$pd)



library(limma)
library(minfi)
library(dplyr)

#estou trocando as linhas do Sample_Group  = 1; II = 2 e III = 3,
#I e adicionando em uma nova coluna chamada SG (será meu novo sample_group)


#criando uma planilha com o pheno
pheno=myLoad$pd


pheno$SG = recode(pheno$Sample_Group, "I" = "1", "II" = "2", "III" = "3")


pheno$SG = recode(pheno$Sample_Group, "I" = "1", "II" = "3")

View(pheno)

#como o Sample_Group estava como char precisei transformar em númerico
pheno$SG=as.numeric(pheno$SG)


## Ajustando os m valor e beta valores

#transformando os beta valores em m valor 
M=logit2(myCombat)
#estou voltando os m value para valores betas 
B=ilogit2(myCombat)




# Criar a matriz de design com as variáveis fenotípicas/ troquei pheno por myLoad$pd
design=model.matrix(~VO2Pico+SG+Idade+Gordura,pheno)


design=model.matrix(~VO2Pico+SG+Idade+Gordura,pheno)


# Calcula a correlação entre amostras duplicadas para ajustar variações técnicas, levando em consideração o design experimental e os blocos de amostras (identificadas por Sample_Name).
corfit=duplicateCorrelation(as.matrix(M), design, block=pheno$Sample_Name)







# Criar a matriz de design com as células/ troquei pheno por myLoad$pd


design=model.matrix(~VO2Pico*B+VO2Pico*NK+VO2Pico*CD4T+VO2Pico*CD8T+VO2Pico*Mono+VO2Pico*Neutro,pheno)


#Outro modelo

design=model.matrix(~VO2Pico + B+ NK+ CD4T+ CD8T+ Mono+ Neutro,pheno)


# Criar coluna de ID removendo _pre ou _post
pheno$ID_participante <- gsub("_pre|_post", "", pheno$Sample_Name)




# Calcula a correlação entre amostras duplicadas para ajustar variações técnicas, levando em consideração o design experimental e os blocos de amostras (identificadas por Sample_Name).
corfit <- duplicateCorrelation(as.matrix(M), design, block = pheno$ID_participante)



# Obtém a correlação média entre amostras duplicadas ajustada.
cor=corfit$consensus

# Ajusta um modelo linear usando os valores de M valor com correção para blocos e correlação.

fit1_M=lmFit(M, design, block=pheno$Sample_Name, correlation=cor)
"rodar aqui agora"
# Aplica o modelo bayesiano com M valor ao ajuste linear.
fit2_M=eBayes(fit1_M)


# Ajusta um modelo linear usando os valores de Beta valor com correção para blocos e correlação.

fit1_B=lmFit(B, design, block=pheno$Sample_Name, correlation=cor)

# Aplica o modelo bayesiano com Beta valor ao ajuste linear.

fit2_B=eBayes(fit1_B)

colnames(design)

coef="SG"
#VO2pico:B
#VO2pico:NK
#VO2pico:CD4T
#VO2pico:CD8T
#VO2pico:Mono
#VO2pico:Neutro



# Extrai os resultados do modelo com todos os cgs e p-valor máximo de 1.




results=topTable(fit2_M, coef=coef,number=Inf, p.value=1)
results_B=topTable(fit2_B, coef=coef,number=Inf, p.value=1)

#fazendo a coluna das cgs ter um nome
results$MarkerName=rownames(results)

# Atribui o valor de logFC de 'results_B' para a linha correspondente em 'results'
results$logFC=results_B[rownames(results),"logFC"]

# Calcula o erro padrão (SE) com base no sigma e no desvio padrão não escalado do modelo ajustado
SE=fit2_B$sigma*fit2_B$stdev.unscaled

# Atribui os valores de SE calculados ao 'results' para o coeficiente especificado
results$SE=SE[rownames(results),coef]


View(results)





library(ggplot2)

# Gerar histograma com eixo X = P.Value, eixo Y = logFC
PreMint <- ggplot(results, aes(x = P.Value)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "P-Values – PreM Pre vs Post intervention",
    x = "P-Value",
    y = "Frequency"
  )

print(PreMint)

ggsave("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\PreMint.jpg", 
       plot = PreMint, dpi = 1200, width = 9, height = 7)



# Calcular FDR para o modelo final (por garantia)
results$FDR <- p.adjust(results$P.Value, method = "fdr")

# Adicionar a coluna com o nome dos CpGs
results$MarkerName <- rownames(results)

View(results)

# Exportar TODOS os CpGs, com logFC, P.Value, FDR etc.
write.csv(results,
          file = "C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\sg.csv",
          row.names = FALSE)








"essa parte para baixo exporta só as cpgs significativas"

table$FDR=p.adjust(table$P.Value,method = "fdr")



hist(table$P.Value, main = "VO2pico")


View(table)

#os que deram siginificativos
TABLE2=table%>%filter(FDR<0.05)


#fazendo a coluna das cgs ter um nome
TABLE2$MarkerName=rownames(TABLE2)

# Exportar a tabela como um arquivo CSV
write.csv(TABLE2, file = paste0("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\cells\\CpGs_FDR_CD48T.csv"), row.names = FALSE)








#sUBINDO OS DADOS DA ASSOCIAÇÃO PARA FILTRAR AS CGS QUE ESTÃO ASSOCIADAS 


directory="C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Dados Camila\\Associação celulas\\Vo2 + cells"
setwd(directory)



table <- results


#filtrar os valores do arquivo results 

cpgs_fdr <- read.csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\CpGs_FDR_SG.csv")


View(cpgs_fdr)
head(results)
head(cpgs_fdr)

#fazendo a coluna das cgs ter um nome
results$MarkerName=rownames(results)

View(results)

# Filtrar os CpGs presentes na lista de CpGs significativos
cpg_list <- cpgs_fdr$CPG  # Ajuste o nome da coluna se necessário, aqui estamos assumindo que a coluna se chama 'cpg'


View(cpg_list)

filtered_results <- results[results$CPG%in% cpg_list, ] 

#dados camila
filtered_results <- results[results$MarkerName %in% cpg_list, ]

View(filtered_results)

# Usando a função match para associar a coluna FDR do arquivo cpgs_fdr com os CpGs filtrados
filtered_results$FDR <- cpgs_fdr$FDR[match(filtered_results$MarkerName, cpgs_fdr$CPG)]


# 6. Exportar os dados filtrados com a coluna FDR para um novo arquivo CSV
write.csv(filtered_results, "C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\CGS_Significativas_Camila_SG.csv", row.names = FALSE)



#FAZENDO O GRAFICO

# Pacotes
library(ggplot2)
library(dplyr)
library(readr)
library(ggtext)



#vo2
# grupo
#  Substitua pelos caminhos dos seus dois arquivos CSV:
dados_post <- read_csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\sg.csv")
dados_pre  <- read_csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Dados Camila\\CGS_Significativas_Camila_VO2.csv")




View(dados_post)
#  Selecionar top 15 CpGs com menor FDR de cada grupo
top_post <- dados_post %>%
  arrange(FDR) %>%
  slice(1:15) %>%
  mutate(Grupo = "PostM")

top_pre <- dados_pre %>%
  arrange(FDR) %>%
  slice(1:15) %>%
  mutate(Grupo = "PreM")

#  Unir os dois
dados_plot <- bind_rows(top_post, top_pre)

#  Selecionar apenas colunas relevantes
dados_plot <- dados_plot %>%
  select(CpG = MarkerName, logFC, Grupo)  # mude "rownames" para o nome da coluna com os IDs CpG se for diferente

#  Definir cores com base em direção e grupo
dados_plot <- dados_plot %>%
  mutate(
    Direcao = case_when(
      Grupo == "PostM" & logFC < 0 ~ "PostM ↓",
      Grupo == "PostM" & logFC >= 0 ~ "PostM ↑",
      Grupo == "PreM"  & logFC < 0 ~ "PreM ↓",
      Grupo == "PreM"  & logFC >= 0 ~ "PreM ↑"
    ),
    Cor = case_when(
      Direcao == "PostM ↓" ~ "darkred",
      Direcao == "PostM ↑" ~ "darkblue",
      Direcao == "PreM ↓"  ~ "lightcoral",
      Direcao == "PreM ↑"  ~ "lightblue"
    )
  )

#  Remover CpGs duplicados (caso apareçam nos dois grupos)
dados_plot <- dados_plot %>% distinct(CpG, Grupo, .keep_all = TRUE)

#  Ordenar CpGs pela média de logFC (para eixo Y mais legível)
ordem_cpg <- dados_plot %>%
  group_by(CpG) %>%
  summarise(m = mean(logFC)) %>%
  arrange(m) %>%
  pull(CpG)

#mudando a forma que a legenda vai aparecer
dados_plot$Direcao <- factor(
  dados_plot$Direcao,
  levels = c("PostM ↓", "PreM ↓", "PostM ↑", "PreM ↑")
)




grafico1 <- ggplot(dados_plot, aes(x = logFC, y = factor(CpG, levels = ordem_cpg), color = Direcao)) +
  geom_point(size = 4) +
  geom_text(aes(label = round(logFC, 4)), hjust = -0.4, size = 3.2, color = "black") +  # valor ao lado do ponto
  scale_color_manual(
    values = c(
      "PostM ↓" = "darkred", "PostM ↑" = "darkblue",
      "PreM ↓"  = "lightcoral", "PreM ↑"  = "lightblue"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(
    x = "logFC",
    y = "CpG",
    color = "Group and Direction",
    title = "VO"["2"]~"peak"
  ) +
  xlim(min(dados_plot$logFC) - 0.002, max(dados_plot$logFC) + 0.005)  # espaço pro texto não colar na borda

print(grafico1)

# grupo

# Carregar os dados
dados_postt <- read_csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\sg.csv")
dados_pret  <- read_csv("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Dados Camila\\CGS_Significativas_Camila_SG.csv")

# Selecionar top 15 CpGs com menor FDR do grupo PostM
top_postt <- dados_postt %>%
  arrange(FDR) %>%
  slice(1:15) %>%
  mutate(Grupo = "PostM")

# Selecionar os mesmos CpGs no grupo PreM
top_pret <- dados_pret %>%
  filter(MarkerName %in% top_postt$MarkerName) %>%
  mutate(Grupo = "PreM")

# Unir os dois
dados_plott <- bind_rows(top_postt, top_pret)

# Selecionar apenas colunas relevantes
dados_plott <- dados_plott %>%
  select(CpG = MarkerName, logFC, Grupo)

# Classificar direção e cor
dados_plott <- dados_plott %>%
  mutate(
    Direcao = case_when(
      Grupo == "PostM" & logFC < 0 ~ "PostM ↓",
      Grupo == "PostM" & logFC >= 0 ~ "PostM ↑",
      Grupo == "PreM"  & logFC < 0 ~ "PreM ↓",
      Grupo == "PreM"  & logFC >= 0 ~ "PreM ↑"
    ),
    Cor = case_when(
      Direcao == "PostM ↓" ~ "darkred",
      Direcao == "PostM ↑" ~ "darkblue",
      Direcao == "PreM ↓"  ~ "lightcoral",
      Direcao == "PreM ↑"  ~ "lightblue"
    )
  )

# Remover duplicatas
dados_plott <- dados_plott %>% distinct(CpG, Grupo, .keep_all = TRUE)

# Ordenar CpGs pela média de logFC
ordem_cpgt <- dados_plott %>%
  group_by(CpG) %>%
  summarise(m = mean(logFC)) %>%
  arrange(m) %>%
  pull(CpG)

#grafico

grafico2 <- ggplot(dados_plott, aes(x = logFC, y = factor(CpG, levels = ordem_cpgt), color = Direcao)) +
  geom_point(size = 4) +
  geom_text(aes(label = round(logFC, 4)), hjust = ifelse(dados_plott$logFC >= 0, -0.6, 1.1), 
            size = 3.5, color = "black") +
  scale_color_manual(
    values = c(
      "PostM ↓" = "darkred", "PostM ↑" = "darkblue",
      "PreM ↓"  = "lightcoral", "PreM ↑"  = "lightblue"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 11, margin = margin(r = 10)),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    x = "logFC",
    y = "CpG",
    color = "Group and Direction",
    title = "Pre and post-exercise intervention"
  ) +
  coord_cartesian(xlim = c(-0.003, 0.005), clip = "off")

print(grafico2)

ggsave("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\vo2 + cells 1505\\grafico_sg_top15.jpg", 
       plot = grafico2, dpi = 1200, width = 9, height = 7)

library(patchwork)

# Combinar lado a lado, bem ajustado
painel <- grafico1 + grafico2 + 
  plot_layout(ncol = 2, widths = c(1, 1)) & 
  theme(legend.position = "right")  # legenda só à direita

# Exibir
print(painel)

# Criar o painel com 2 linhas e 2 colunas
panel <- grid.arrange(grafico1, grafico2, nrow = 2, ncol = 2)

# Salvar o painel com resolução de 1200 dpi
ggsave("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\grafico_painel3004.png", plot = panel, dpi = 1200, width = 8, height = 6)









