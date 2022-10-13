

#*******************-Network Theory Project Script****************************#






#install_github("insilico/npdro")
#install.packages('npdro')
#install.packages("devtools")

library(igraph)
library(devtools)
library(npdro)
library(ggplot2)




#---------------------------- SIMULATE GENE EXPRESSION DATA -----------------------------------------
num.samples <- 200
num.variables <- 100
dataset <- npdro::createSimulation2(num.samples=num.samples,
                                    num.variables=num.variables,
                                    pct.imbalance=0.5,
                                    pct.signals=0.2,
                                    main.bias=0.4,
                                    interaction.bias=0.8,
                                    hi.cor=0.95,
                                    lo.cor=0.2,
                                    mix.type="main-interactionScalefree",
                                    label="class",
                                    sim.type="mixed",
                                    pct.mixed=0.5,
                                    pct.train=0.5,
                                    pct.holdout=0.5,
                                    pct.validation=0,
                                    plot.graph=F,
                                    verbose=T)
dats <- rbind(dataset$train, dataset$holdout, dataset$validation)
colnames(dats)
dim(dats)

# > dim(dats)
# [1] 200 101 # there are 200 samples, 100 genes, and 1 class variable.


# functional variables: main effects and interaction effects
dataset$signal.names
dataset$main.vars # 10 main variables
dataset$int.vars# 10 interaction variables... the remainder are var
functional_mask <- colnames(dats[,-ncol(dats)]) %in% dataset$signal.names
# mainvars and intvars are the signal vars
# underlying variable network
A.graph <- graph_from_adjacency_matrix(dataset$A.mat, mode="undirected")
isSymmetric(dataset$A.mat)
plot(A.graph, vertex.size = 8)
deg.vec <- degree(A.graph)
hist(deg.vec, breaks = 30,
     main = "Simulated Data Degree",
     xlab = "degree",
     ylab = "frequency")



#------------------------------------------------------------------------------------------


# set working directory
setwd("~/Dropbox/Ph.D/Classes/Spring_2022/Network_Theory/Final Project")

source('npdro_functions.R')

#---------> reGAIN fit

test.regain <- regain(dats, 
                      numVars=dim(dats)[2]-1, 
                      indVarNames=colnames(dats)[-ncol(dats)], 
                      depVarName="class", 
                      reg.fam="binomial", 
                      nCovs=0, 
                      excludeMainEffects=F)

end.time <- Sys.time()
end.time - start.time
#test.regain
#x<-pmax(test.regain$stdbetas)
#dim(x)
betas <- as.matrix(test.regain$stdbetas)
regain.nodiag <- betas
diag(regain.nodiag) <- 0 # not going to consider main effects
regain.nodiag.adj <- as.matrix(regain.nodiag>1.75)+0  # define adjacency matrix based on 2.0 beta coeff filter
library(igraph)
A.graph <- graph_from_adjacency_matrix(regain.nodiag.adj, mode="undirected")
#isSymmetric(dataset$A.mat)
layout.fr.iterate <- layout_with_fr(A.graph, niter = 10000)
# node radius by main effect std beta
plot.igraph(A.graph,layout=layout.fr.iterate, 
            #vertex.color=id1_colors,               # color by identifier1
            #vertex.color=louvain$membership,      # color by cluster
            #vertex.color=knn_clust_colors,        # also cluster colors
            vertex.size=abs(diag(betas))
)
plot(A.graph)
my.ranks <- EpistasisRank(test.regain$stdbetas)
my.ranks[1:20,]

# now we have






#---------> dcGAIN fit


dats$class <- as.numeric(dats$class) - 1 # class must be numeric 1s and 0s


test.dcgain <- dcgain(dats)




#--------------------> Define CINC 



cinc <- function(A,
                 gamma = 0.85) {
  
  
  k <- rowSums(A)
  R <- matrix(0,1, ncol(A))
  for (i in seq(1:ncol(A))){
    term3 = c()
    for (j in seq(1:ncol(A))){
      
      if (k[j] != 0){
        
        denom <- sum(diag((A)))
        columns <- ncol(A)
        A_ii <- A[i,i]
        term1 <- A[i,i]/(ncol(A)*sum(diag((A))))
        term2 <- (1.0 - gamma)/ncol(A)

        if (j!=i){
          term3_ij <- A[i,j]*R[j]/k[j]
        } else {
          term3_ij <- 0.0
        }
        
        term3 <- c(term3, term3_ij)
      }
      
      else{
        term1 <- A[i,i]/(ncol(A)*sum(diag((A))))
        term2 <- 1/ncol(A)
        term3 <- 0.0
      }
      
    }
    
    R[i] <- term1 + term2 + sum(term3)
    
  }
  
  return(R)
}



# get the gene importance score 

R <- cinc(test.regain$stdbetas)

hist(R, breaks = 30,
     main = "CINC Gene Importance Score",
     xlab = "Score",
     ylab = "frequency")

plot()


# gaussian mixture model assuming two clusters

#install.packages('ClusterR')
# library(ClusterR)
# 
# gmm <- GMM(R, 2, dist_mode = 'maha_dist', seed_mode = 'random_subset', km_iter = 10, 
#            em_iter = 10, verbose = F)

library(mclust)

mcl.model <- Mclust(R, 2)

mcl.model$classification

pred <- dats[,which(mcl.model$classification==2)]



pred.idx <- which(mcl.model$classification==2)

act.idx <- dataset$int.vars





TPR <- function(pred.idx, act.idx){
  
  TP <- length(intersect(pred.idx, act.idx))
  
  FN <- length(act.idx) - length(intersect(pred.idx, act.idx))
  
  TPR <- TP/(FN+TP) # something must be going wrong
  
  return(TPR)
  
}




# > TPR
# [1] 0.2




plot(mcl.model, what = 'classification', main = 'Mclust Classification')

#------------------------------------- SPECTRAL CLUSTERING ---------------------------




# Have to get the A matrix from the B matrix since the A matrix in simulation2 doesn't
# include the perturbance


betas <- as.matrix(test.regain$stdbetas)
regain.nodiag <- betas
diag(regain.nodiag) <- 0
regain.nodiag.adj <- as.matrix(regain.nodiag>1.75)+0  # this is where the adjacency matrix is defined
library(igraph)
A.graph <- graph_from_adjacency_matrix(regain.nodiag.adj, mode="undirected")
layout.fr.iterate <- layout_with_fr(A.graph, niter = 10000)
# node radius by main effect std beta
plot.igraph(A.graph,layout=layout.fr.iterate, 
            #vertex.color=id1_colors,               # color by identifier1
            #vertex.color=louvain$membership,      # color by cluster
            #vertex.color=knn_clust_colors,        # also cluster colors
            vertex.size=abs(diag(betas))
)

#----------> get the Laplacian

# G has K disconnected components. This means that K eigenvalues
# of the laplacian matrix L are zero: Lap * y.i = lambda.i * y.i


# we'll try on the immuno network first since we know how many defined subclusters
# are present

#g <- immuno

g <- A.graph

#g.adj <- get.adjacency(g)

Lap <- function(adj) {
  
  g.adj <- as.matrix(adj)
  
  g.deg <- rowSums(as.matrix(g.adj))
  
  g.Lap <- diag(g.deg) - g.adj
  
  return(g.Lap)
  
}



g.Lap <- Lap(regain.nodiag.adj)

# choose the largest non-zero eigenvec

n <- length(V(g))
x <- eigen(g.Lap)$vectors[,n-1]
x_val <- eigen(g.Lap)$values[n-1]

plot(x)


# Find step functions based on derivative

deriv <- function(x){
  dy.list <- c()
  for (i in seq(0,length(x))){
    dy <- x[i+1]-x[i]
    dy.list <- c(dy.list, dy)
  }
  return(as.array(dy.list))
}


dy.dx <- deriv(x)
df <- data.frame(x, dy.dx)
threshold <- 0.05 * max(abs(dy.dx), na.rm = TRUE)
df['index'] <- seq(1, nrow(df))
df['bound'] <- ifelse(abs(df$dy.dx) > threshold, 1, 0)


ggplot(df, aes(index, x)) + 
  geom_point(col = ifelse(df$bound ==1, 'red','black')) # maybe mention here that
# as you decrease the threshold, more points deviate from the baseline

#plot(x, xlim=c(0,500), ylim = c(-0.05, -0.04)) # all the same eigenvalues

library(ggplot)

layout.fr.iterate <- layout_with_fr(A.graph, niter = 10000)
# node radius by main effect std beta
plot.igraph(A.graph,layout=layout.fr.iterate, 
            #vertex.color=id1_colors,               # color by identifier1
            #vertex.color=louvain$membership,      # color by cluster
            #vertex.color=knn_clust_colors,        # also cluster colors
            vertex.size=abs(diag(betas))
)

plot(x) # from this plot we can easily index the disconnected components

ggplot(data.frame(x))


#----------> these bounds index the original genes

g.eigs <- eigen_centrality(g, scale=FALSE)$vector

#V(g)$size <- 40*g.eigs

V(g)$size <- 5 #*eigen(g.Lap)$values
V(g)$color <- df$bound
V(g)$color <- ifelse(V(g)$color == 1,"green","yellow")

plot(g)


plot.igraph(A.graph,layout=layout.fr.iterate, 
              vertex.color= id1_colors,               # color by identifier1
              #vertex.color=louvain$membership,      # color by cluster
              #vertex.color=knn_clust_colors,        # also cluster colors
              vertex.size=abs(diag(betas))
  )

  # doesnt really help with clustering
  
  
  
  #------> lets try the NG-Jordan-Weiss algorithm
  
  
#---------------------------------------------- NG-JORDAN-WEISS (NGW) ALGORITHM -------------------------
  
  # 1.) Form the affinity matrix A defined by 
  
  affinity <- function(s, sigma){
    
    A <- matrix(data = 0, nrow = length(s), ncol = length(s))
  
    for (j in seq(1,ncol(A))){
      
      for (i in seq(1,nrow(A))){
        
        if (i!=j){
          
          term <- -((abs(s[i] - s[j]))^2)/(2*(sigma^2))
          A[i,j] <- exp(term)
        } else {
          A[i,j] <- 0
        }
      }
    }
    return(A)
  }
  
  
# for the s matrix, we will need the mean expression for the health controls


cntrl.genes <- dataset$train[,which(dataset$train$class == -1)]
case.genes <- dataset$train[,which(dataset$train$class == 1)]

cntrl.means <- rowMeans(cntrl.genes)

case.means <- rowMeans(case.genes)


NJW <- function(expr.means, sigma, k){
  
  A.mat <- affinity(case.means, sigma)
  D.mat <- diag(rowSums(A.mat), nrow = nrow(A.mat), ncol = ncol(A.mat))
  
  # take the inverse sqrt of the D matrix
  
  #install.packages('expm')
  
  library(expm)
  
  sqrt.D <- sqrtm(D.mat)
  
  D.negsqrt <- solve(sqrt.D, power = -1)
  
  # calculate the graph laplacian
  
  L.mat <- D.negsqrt%*%A.mat%*%D.negsqrt
  
  # solve eigenvalue decomposition and construct X
  
  g.eig.vecs <- eigen(L.mat)$vectors
  g.eig.vals <- eigen(L.mat)$values
  
  #k <- 2
  
  
  X <- t(g.eig.vecs[seq(1,k),])
  
  # renormalize each of X's rows to have unit length
  
  #install.packages('wordspace')
  library(wordspace)
  X.norm <- normalize.rows(X)
  
  
  # Perform K-Means clustering
  
  km <- kmeans(X.norm, centers = 2, iter.max = 10, nstart = 1)
  
  return(km$cluster)
  
}



clusters <- NJW(cntrl.means, 1, 2)

# 
# act.clust <- data.frame(seq(1,100,1))
# act.clust <- ifelse(data.frame==idx)

g <- A.graph

V(g)$size <- 5 #*eigen(L.mat)$values
V(g)$color <- clusters


plot(g)


# compare the clustering between case and control

case.clusters <- km$cluster
  
#cntrl.clusters <- km$cluster

pred.idx <- which(case.clusters != cntrl.clusters)

tpr <- TPR(pred.idx, dataset$int.vars)


# this predicts far too many.. hardly better than random

pred.idx

dataset$int.vars


length(pred.idx)
  
#------------------------------------------ TRY ON REAL GENE EXPRESSION DATA ----------------



# load gene expression data
# set the current working directory first (setwd())
load("sense.filtered.cpm.Rdata")
dim(sense.filtered.cpm) 

sense.filtered.cpm

x<-nrow(sense.filtered.cpm)*ncol(sense.filtered.cpm)
x<-ncol(sense.filtered.cpm)
print(x)
colnames(sense.filtered.cpm)


# phenotype (mdd/hc) is in a separate file
# match phenotypes to expression data subject ids
subject.attrs <- read.csv("Demographic_symptom.csv", stringsAsFactors = FALSE)
dim(subject.attrs)  # 160 subjects x 40 attributes
colnames(subject.attrs)  # interested in X (sample ids) and Diag (diagnosis)
subject.attrs$X
subject.attrs$Diag


# There are more samples in the demographic data than in the gene expression data, so we need to match gene expression samples with thier diagnosis


#* 2) How many cases and controls are there? 
#*  There are 5692 healthy controls.. after you filter there are only 79. 

hc_vec <- which(subject.attrs$Diag == "HC")
hcs<-sum(hc_vec)



library(dplyr) # install.packages("dplyr")
# create a phenotype vector
# grab X (subject ids) and Diag (Diagnosis) from subject.attrs that 
# intersect %in% with the RNA-Seq data
phenos.df <- subject.attrs %>% 
  filter(X.1 %in% colnames(sense.filtered.cpm)) %>%
  dplyr::select(X.1, Diag)  
colnames(phenos.df) # $Diag is mdd diagnosis
# grab Diag column and convert character to factor
mddPheno <- as.factor(phenos.df$Diag)  # this is our phenotype/class vector 

summary(mddPheno) # MDD -- major depressive disorder, HC -- healthy control


# Normalized and transform
library(preprocessCore)
mddExprData_quantile <- normalize.quantiles(sense.filtered.cpm)
mddExprData_quantileLog2 <- log2(mddExprData_quantile)

# coefficient of variation filter sd(x)/abs(mean(x))
CoV_values <- apply(mddExprData_quantileLog2,1,
                    function(x) {sd(x)/abs(mean(x))})
# smaller threshold, the higher the experimental effect relative to the 
# measurement precision
thresh <- .02
sum(CoV_values<thresh)
# there is one gene that has 0 variation -- remove
sd_values <- apply(mddExprData_quantileLog2,1, function(x) {sd(x)})
rownames(mddExprData_quantileLog2)[sd_values==0]  
# filter the data matrix 
GxS.covfilter <- mddExprData_quantileLog2[CoV_values< thresh & sd_values>0,]
dim(GxS.covfilter)

# now we have only 781 genes

#----> instead of getting ajacency matrix here lets seperate the dataframe into HC and MDD

hc.genes <- GxS.covfilter[,which(phenos.df$Diag == 'HC')]
mdd.genes <- GxS.covfilter[,which(phenos.df$Diag == 'MDD')]

hc.means <- rowMeans(hc.genes)

mdd.means <- rowMeans(mdd.genes)

#------> lets try the NG-Jordan-Weiss algorithm

# 1.) Form the affinity matrix A defined by 

affinity <- function(s, sigma){
  
  A <- matrix(data = 0, nrow = length(s), ncol = length(s))
  
  for (j in seq(1,ncol(A))){
    
    for (i in seq(1,nrow(A))){
      
      if (i!=j){
        
        term <- -((abs(s[i] - s[j]))^2)/(2*(sigma^2))
        A[i,j] <- exp(term)
        } else {
        A[i,j] <- 0
      }
    }
  }
  return(A)
}


A.mat <- affinity(hc.means, 1)

D.mat <- diag(rowSums(A.mat), nrow = nrow(A.mat), ncol = ncol(A.mat))

# take the inverse sqrt of the D matrix

install.packages('expm')

library(expm)

sqrt.D <- sqrtm(D.mat)

D.negsqrt <- solve(sqrt.D, power = -1)

# calculate the graph laplacian

L.mat <- D.negsqrt%*%A.mat%*%D.negsqrt

# solve eigenvalue decomposition and construct X

g.eig.vecs <- eigen(L.mat)$vectors
g.eig.vals <- eigen(L.mat)$values

k <- 8


X <- t(g.eig.vecs[seq(1,k),])

# renormalize each of X's rows to have unit length

#install.packages('wordspace')
library(wordspace)
X.norm <- normalize.rows(X)


# Perform K-Means clustering

km <- kmeans(X.norm, centers = 5, iter.max = 10, nstart = 1)

km$cluster

df <- data.frame(X1 = X.norm[,1],X2 = X.norm[,2],X3 = X.norm[,3],X4 = X.norm[,4],X5 = X.norm[,5],X6 = X.norm[,6],
                 X7 = X.norm[,7],X8 = X.norm[,8], cluster=km$cluster)

plot(x = df$X1, y = df$X2,  col = df$cluster)


# lets get the adjacency matrix

mddCorr<-cor(t(GxS.covfilter))  # correlation between genes

thresh<-.70 # controls sparsity of network
# threshold and turn T/F to 1/0
adjMat <- (abs(mddCorr)>thresh)+0
diag(adjMat) <- 0  # remove self-connections
rownames(adjMat) <- row.names(GxS.covfilter)
colnames(adjMat) <- row.names(GxS.covfilter)

adjMat

GxS.covfilter

g.adj <- adjMat


g <- graph_from_adjacency_matrix(g.adj, mode='undirected')

layout.fr.iterate <- layout_with_fr(g, niter = 10000)
# node radius by main effect std beta
plot.igraph(g,layout=layout.fr.iterate, 
            #vertex.color=id1_colors,               # color by identifier1
            #vertex.color=louvain$membership,      # color by cluster
            #vertex.color=knn_clust_colors,        # also cluster colors
            vertex.size=abs(diag(betas))
)


g.deg <- rowSums(as.matrix(g.adj))

g.Lap <- diag(g.deg) - g.adj

#-----> calculate the normalized laplacian

D.Lap <- diag(g.deg)^-0.5

g.Lap.norm <- (diag(g.deg)^-0.5)%*%g.Lap%*%(diag(g.deg)^-0.5)


n <- length(V(g))

g.eig.vecs <- eigen(g.Lap)$vectors
g.eig.vals <- eigen(g.Lap)$values

lb <- min(which(g.eig.vals == 0)) # values below this are positive
ub <- max(which(g.eig.vals == 0)) # values above this are negative


# x <- eigen(g.Lap)$vectors[,lb-6]
# x_val <- eigen(g.Lap)$values[lb-6]

x <- eigen(g.Lap)$vectors[,n-5]
x_val <- eigen(g.Lap)$values[n-5]

plot(x)

plot(x, xlim=c(0,1000), ylim = c(-0.005, 0.005)) # all the same eigenvalues



}
