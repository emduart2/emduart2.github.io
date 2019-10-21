## 20.10.2019
## Introduction to R
## For a more detailed introduction go to : https://www.kahle.io/r-crash
## Introduction to algstat package: https://github.com/dkahle/algstat

## Instaling packages  ##
##  install.packages() or install_github()
install.packages("devtools") # install packege
library("devtools")          # load the installed package 

## Help in R ##
?lm
?fisher.test
?summary


## Variables and assignments
## use " <- " for assignmet
x<-5
x
x+1
ls()

## Data structures
v<- c(0,0,1)  # This is a column vector
v
names(v)<- c("a","b","c") # Give names to the entries of the vector
v

tab<- as.table(matrix(c(43,16,3,
                        6,11,10,
                        9,18,16), nrow = 3, byrow= TRUE))
dimnames(tab)<- list( X = c("X1","X2","X3"), Y = c("Y1","Y2","Y3"))
?as.table
?dimnames

## Install packaged relevant to Algebraic Statistics
