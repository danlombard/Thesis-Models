#Contraints matrix construction
t1 <- Sys.time()
library(lpSolve)
# Constraint 1
n <- 8
m <- 2
u <- 4

# Matrix generation
items <- n
m.vect <- rep(0, items*items)
for(i in 1:(items*items)-1){
  m.vect[i] <- as.integer(runif(1,1,10))
}
m.vect
x1<-1
for(j in 1:items){
  m.vect[x1] <- 0
  x1 <- x1+items+1
}
m.vect
matrix <- matrix(m.vect,nrow = items, byrow = TRUE)
x1 <- 1
y <- 1
for(j in 1:(items-1)){
  for(i in 1:(items-(items-y))){
    matrix[x1+(items-y),(items-y)] = matrix[(items-y),x1+(items-y)]
    x1 <- x1+1
  }
  x1 <- 1
  y <- y+1
}
W <- matrix
W

A <- matrix(0,3*n*m*(n-1)/2 +2*n+m,(m*n*(n+1))/2)
b <- matrix(0,3*n*m*(n-1)/2+2*n+m,1)
x <- 1
# Constraint 1
for(i in 1:(n-1)){
  for(j in (i+1):n){
    for(k in 1:m){
      A[x,m*n+(k-1)*n*(n-1)/2+(i-1)*(2*n-i)/2+j-i] = 1
      A[x,(k-1)*n+i] = -1
      
      b[x] = 0
      x <- x+1
    }
  }
}

# Constraint 2
for(i in 1:(n-1)){
  for(j in (i+1):n){
    for(k in 1:m){
      A[x,(k-1)*n+j] = -1
      A[x,m*n+(k-1)*n*(n-1)/2+(i-1)*(2*n-i)/2+j-i] = 1
      b[x] = 0
      x <- x+1
    }
  }
}
# Constraint 3

for(i in 1:(n-1)){
  for(j in (i+1):n){
    for(k in 1:m){
      A[x,(k-1)*n+i] = 1
      A[x,(k-1)*n+j] = 1
      A[x,m*n+(k-1)*n*(n-1)/2+(i-1)*(2*n-i)/2+j-i] = -1
      b[x] = 1
      x <- x+1
    }
  }
}
#Contraint 4

for(i in 1:n){                          
  for(k in 1:m){
    A[x,(k-1)*n+i] = 1
  }
  b[x] = 1
  x <- x+1
}

# Contraint 5

for(i in 1:n){                          
  for(k in 1:m){
    A[x,(k-1)*n+i] = -1
  }
  b[x] = -1
  x <- x+1
}


# Contraint 6

for(k in 1:m){
  for(i in 1:n){
    A[x,(k-1)*n+i] = 1
    
  }
  b[x] = u
  x <- x+1
}

# Direction matrix

dir <- matrix("<=",3*n*m*(n-1)/2+2*n+m, 1)

# Objective matrix

obj <- matrix(0,m*n*(n+1)/2, nrow = 1)

for(k in 1:m){
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      obj[m*n+(k-1)*n*(n-1)/2+(i-1)*(2*n-i)/2+j-i] = W[i,j]
    }
  }
}
direction <- t(dir)
rhs <- t(b)
print(lp("max",obj, A, direction, rhs, int.vec = 1:(m*n*(n+1))/2))
xvector <- lp("max",obj, A, direction, rhs, int.vec = 1:(m*n*(n+1))/2)$solution
t2 <- Sys.time()
t2-t1

WMatrix <- data.frame(W)

columnNamesW <- c(1:n)

colnames(W) <- columnNamesW

W

# Manipulation of xvector 


solution <- c()
for (i in 1:((m*n*(n+1))/2)){
  if(xvector[i] == 1){
    solution[i] <- i
  }
}

for(i in 1:((m*n*(n+1))/2)){
  if(is.na(solution[i]) == TRUE){
    solution[i] <- 0
  }
}


solutionFinal <- c()
x <- 1
y <- 1
for(j in 1:m){
  for( i in 1:n){
    if((solution[x] > 0) && (solution[x] <= j*n)){
      solutionFinal[y] <- paste0("Item ", i, " is in container ", j)
      y <- y+1
      x <- x +1
    } else {
      x <- x + 1
    }
  }
}

solutiondf <- data.frame(solutionFinal)
solutiondf



