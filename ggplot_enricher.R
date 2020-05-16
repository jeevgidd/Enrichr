###################################
#     Vanessa Escolano Maso       #
#     14/05/2020                  #
#     ggplot para Enrichr         #
###################################

# carrega pacotes necessarios
# faz grafico
library(ggplot2)
# tem as funcoes para manipulacao/limpeza de dados
library(tidyverse)

# abre o arquivo de vias enriquecidas de genes subexpressos
reactome_down <- read_delim("C:/Users/vanes/CHIKV_Natalia/Enrichr/down_regulated_pathways/Reactome_2016_table.txt", delim = "\t", col_names = T)

# abre o arquivo de vias enriquecidas com genes superexpressos
reactome_up <- read_delim("C:/Users/vanes/CHIKV_Natalia/Enrichr/up_regulated_pathways/Reactome_2016_table (1).txt", delim = "\t", col_names = T)

# filtra apenas vias com p-adjusted < 0.05 para os down
down_significant <- reactome_down %>% 
  filter(`Adjusted P-value` < 0.05)

# cria uma coluna que vai ajudar a estabelecer a cor da barra no grafico 
down_significant <- down_significant %>%
  add_column(status = "Down")

# filtra apenas vias com p-adjusted < 0.05 para os up
up_significant <- reactome_up %>% 
  filter(`Adjusted P-value` < 0.05)

# cria uma coluna que ajuda a estabelecer a cor das barras do grafico
up_significant <- up_significant %>%
  add_column(status = "Up")

# junta as vias enriquecidas tanto para gene up como down
all_p_adjust_0.05_pathways <- rbind(down_significant, up_significant)

# cria uma coluna com o valor de -log10 p-adjusted
all_p_adjust_0.05_pathways <- all_p_adjust_0.05_pathways %>% 
  mutate(log_adjusted = -log10(`Adjusted P-value`))

# retirando Homo sapiens do nome das vias
all_p_adjust_0.05_pathways$Term <- str_remove(all_p_adjust_0.05_pathways$Term, "Homo sapiens")

# organiza as vias na ordem que eu quero plotar o grafico
# ordem crescente de acordo com status ele vai organizar por ordem alfabetica, primeiro Down (D) depois Up (U) e depois de acordo com o valor de log em ordem crescente tambem. Quanto maior o log menor o p
all_p_adjust_0.05_path_arranged <- all_p_adjust_0.05_pathways %>% 
  arrange(status, log_adjusted)

#transformando em caracter todas as colunas
caracter <- apply(all_p_adjust_0.05_path_arranged, 2, as.character)

#atribuindo os characteres da coluna 1 que eh a que contem o nome das vias a um objeto, que sera um vetor de caracteres
nome_vias_vetor <- caracter[,1]

# transforma em fator ordenado os termos
# como a ordem esta correta pq deu arrange e depois transformou cada nome em character e depois colocou em um vetor temos um vetor com a ordem correta dos meus levels dentro de nome_vias_vetor, entao ? s? colocar isso dentro de c de levels
# o mutate transformara a coluna Term em vetor com os niveis ordenados de acordo com a ordem contida de nomes_vias_vetor que ja foi ordenado com o arrange la em cima
all_p_adjust_0.05_path_arranged <- all_p_adjust_0.05_path_arranged %>%
  mutate(Term = factor(Term, levels = c(nome_vias_vetor), ordered = TRUE))

# plota o grafico
grafico <- all_p_adjust_0.05_path_arranged %>% 
  # abre a janela grafica
  ggplot() +
  # estabelece como eixo x o Term que e o nome das vias
  # estabelece como eixo y o log_adjusted que ? o -10log do p valor ajustado
  # como a cor de preenchimento da barra depende do status fill vai dentro do aes e ? igual status, mesma coisa para o color, queremos que o contorno da barra seja da mesma cor que a barra
  geom_col(aes(x = Term, y = log_adjusted, fill = status, color = status)) +
  # scale_fill_brewer determina a palheta de cores para o preenchimento da barra. direction = -1 para os genes down ficarem azuis e up vermelhos, se tirar esse argumento as cores ficam invertidas
  scale_fill_brewer(palette = "Set1", direction = -1) +
  # scale_colour determina a palheta para contorno das barras, como queremos o preenchimento e o contorno de cores iguais, os argumentos de scale fill e scale colour sao os mesmos
  scale_colour_brewer(palette = "Set1", direction = -1) +
  # adiciona nome nos eixos x, y e o titulo
labs(x = "Reactome pathway name",
     y = "-log10 adjusted p-value",
     title = "Enriched pathways for Chikungunya infection") +
  # rota o grafico 90?
  coord_flip() +
  # retira a grade de tr?s, deixando fundo branco
  theme_classic()

# salva grafico no formato tiff
grafico %>% 
  ggsave(filename = "Vias_Reactome.tiff", device = "tiff", path = "C:/Users/vanes/CHIKV_Natalia/Enrichr/All_pathways_ggplot_R/")

# salva grafico no formato png
grafico %>% 
  ggsave(filename = "Vias_Reactome.png", device = "png", path = "C:/Users/vanes/CHIKV_Natalia/Enrichr/All_pathways_ggplot_R/")

# salva grafico no formato jpeg
grafico %>% 
  ggsave(filename = "Vias_Reactome.jpeg", device = "jpeg", path = "C:/Users/vanes/CHIKV_Natalia/Enrichr/All_pathways_ggplot_R/")

# # salva grafico no formato svg
grafico %>% 
  ggsave(filename = "Vias_Reactome.svg", device = "svg", path = "C:/Users/vanes/CHIKV_Natalia/Enrichr/All_pathways_ggplot_R/")

#####################################################################

# Para verificar as cores das paletas que podem ser usadas no ggplot2

# carrega o pacote de cores que ggplot usa
# library(RColorBrewer)

# verifica o nome e as op?oes da paleta de cores
# display.brewer.all()