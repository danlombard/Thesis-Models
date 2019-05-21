library(igraph)
library(Surrogate)
library(ggplot2)
cap <- c()
table <- c()
numTables <- 20
maxk <- 15
ParetoMatrix <- matrix(0, numTables, maxk)

for(k in 1:maxk){
  t1 <- Sys.time()
  
  numFamilies <- 26
  
  f <- c(35,7,6,5,5,5,5,5,4,3,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
  problemsize <- 99
  
  #Vector of family sizes 
  
  for(t in 1:numTables){
    
    TableCapacity <- ceiling(problemsize/t)  # Table capacity
    x <- 1 
    MaxTableNode <- numFamilies + t + 1  #Plus one is for the sink node
    familyNodes <- c()
    for(i in 1:numFamilies+1) { 
      familyNodes <- c(familyNodes,1,x+1,f[x])
      x<-x+1
    }
    familyNodes;
    famArcs <- t(matrix(familyNodes, nrow = 3))
    famArcs;
    
    #Diversity index
    z <- 2       #First Family node (this value will always be 2)
    TableNodes <- c()
    for(j in 1:numFamilies+2){
      LastFamNode <- numFamilies + 1       # Number of families plus 1
      while(LastFamNode < MaxTableNode){ 
        TableNodes <- c(TableNodes,z,LastFamNode+1,k)
        LastFamNode = LastFamNode+1
      }
      z = z+1
    }
    TableNodes ;
    TbNodes <- t(matrix(TableNodes, nrow = 3))
    TbNodes;
    
    LastFamNode <- numFamilies + 1
    FirstTableNode <- LastFamNode + 1       #First table node
    sink <- MaxTableNode + 1     #Sink node
    SinkNode <- c()
    while(FirstTableNode < MaxTableNode+1){
      SinkNode <- c(SinkNode,FirstTableNode,sink,TableCapacity)
      FirstTableNode = FirstTableNode + 1
    }
    SinkNode;
    SkNodes <- t(matrix(SinkNode, nrow = 3))
    SkNodes;
    
    
    E <- rbind(
      famArcs,TbNodes,SkNodes
    )
    E;
    colnames(E) <- c("from", "to", "capacity")
    g1 <- graph_from_data_frame(as.data.frame(E))
    f1 <- max_flow(g1, source=V(g1)["1"], target=V(g1)[sink])$flow
    f2 <- max_flow(g1, source=V(g1)["1"], target=V(g1)[sink])$value
    
    ParetoMatrix[t,k] <- f2
    cap[t] <- TableCapacity
    table[t] <- t
  }
}

columnNames <- c()
TableCap <- data.frame(table, cap, ParetoMatrix) 

colnames(TableCap) <- c("Table", "Capacity", paste("k=",1:maxk,sep = ""))

TableCap

yvalues <- c()
for(i in 1:numTables){
  for(j in 1:maxk){
    if(ParetoMatrix[i,j] == sum(f)){
      yvalues[i] <- j
      break
    }
  }
}

xvalues <- c(1:numTables)

for(i in 1:numTables){
  if(is.na(yvalues[i])=="TRUE"){
    xvalues <- xvalues[-i]
  }
}

yvalues <- na.omit(yvalues)


xvalues
yvalues <- yvalues[ !is.na( yvalues ) ]

df <- data.frame(xvalues, yvalues)


ggplot(df, aes(x=xvalues, y=yvalues)) +
  geom_line()+
  geom_point() +
  labs(y = "Similarity threshold k") +
  labs(x = "Number of containers")
  
t2 <- Sys.time()
comptimeMaxflow <- t2-t1
comptimeMaxflow



