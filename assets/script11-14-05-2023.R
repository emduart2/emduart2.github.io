####### Exercícios. #######################

# significancia alpha  = 0.05
#crear os dados 
data <- data.frame(person = rep(1:5, each=4),
                   drug = rep(c(1, 2, 3, 4), times=5),
                   score = c(30, 28, 16, 34, 14, 18, 10, 22, 24, 20,
                             18, 30, 38, 34, 20, 44, 26, 28, 14, 30))

#view data
data
#
friedman.test(y=data$score, groups=data$drug, blocks=data$person)

#perform post-hoc tests
pairwise.wilcox.test(data$score, data$drug, p.adj = "bonf")
?pairwise.wilcox.test
# Teste do Friedman

# Para um nível de significância de 0.05, rejeita-se H0,
# dado que o valor-p=0.0009<0.05.

?friedman.test()

#Ex. 1
# --- apla= 0.05
y <- c(4000,3210,6120,1600,1040,2410,1600,647,2210,
       1200,570,2060,840,445,1400,352,156,249,
       224,155,224,200,99,208,184,70,227)
Paciente<- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9)
method<-c(rep(c("A","B","C"),9))
friedman.test(y,method,Paciente)
pairwise.wilcox.test(y, method, p.adj = "bon")

MetA<-c(4000,1600,1600,1200,840,352,224,200,184)
MetB<-c()
getwd()
library(readxl)
MetodosPacientes <- read_excel("/Users/ElianaDuarte/Documents/Teaching/2023-FundamentosDeEstatistica/EA/ex1.xlsx")
boxplot(MetodosPacientes$`Metodo A`,MetodosPacientes$`Metodo B`,
        MetodosPacientes$`Metodo C`,names=c("A","B","C"))
hist(MetodosPacientes$`Metodo A`)
hist(MetodosPacientes$`Metodo B`)
hist(MetodosPacientes$`Metodo C`)
MetodosPacientes
stack(MetodosPacientes)
