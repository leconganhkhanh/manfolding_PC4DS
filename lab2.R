
# Read a txt file, named "data.txt"
minimnist <- as.matrix(read.table("data.txt"))
labels <- as.matrix(read.table("labels.txt"))
#head(my_data)
p <- sqrt(784)
image(matrix(minimnist[3,],p)[,p:1], col = grey(256:1/256))


# faire un ACP sur chiffre 3

# install package
#install.packages(c("FactoMineR", "factoextra"))
# loading the package
library("FactoMineR")
library("factoextra")

all3 <- minimnist[labels == 3,]

cte3 <- which(apply(all3,2, sd) == 0)

cte3.pca <- PCA(all3[, -cte3], graph = FALSE)

pca3 <- prcomp(scale(all3[,-cte3]))
biplot(pca3)
plot(pca3$rotation[,1:2])

# identify the pointe interesse
ids <- identify(pca3$rotation[,1:2]) # return the coordone of point

## Exercise 3 ####

matsimu <- function (n, D) matrix(runif(n*D), n, D)
mat <- matsimu(100,2)
# norme de chaque observation
norm2 <- sqrt(rowSums(mat^2))
# visu norm2
hist(norm2)
plot(density(norm2))