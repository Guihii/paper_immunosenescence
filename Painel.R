library(magick)     # para ler e manipular imagens
library(cowplot)    # para montar o painel


img1 <- image_read("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\cells 1505\\B cells.jpg")
img2 <- image_read("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\cells 1505\\NK cells.jpg")
img3 <- image_read("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\cells 1505\\Monocytes.jpg")


p1 <- ggdraw() + draw_image(img1)
p2 <- ggdraw() + draw_image(img2)
p3 <- ggdraw() + draw_image(img3)


painel <- plot_grid(p1, p2, p3,
          labels = c("A", "B", "C"),
          ncol = 2)  # ou nrow = 2


print(painel)

ggsave("C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\DOUTORADO FAPESP\\Artigo DOC - Guilherme\\Artigo BEPE\\ANALISES 2804\\Meus dados\\Associações celulas\\cells 1505\\painel_volcano.jpg", width = 12, height = 8, dpi = 1200)
