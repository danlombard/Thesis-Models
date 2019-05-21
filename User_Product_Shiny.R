library(shiny)
library(DT)
ui <- fluidPage(
  titlePanel("Two approaches to the SPP of items to containers"),
  fluidRow(
    column(3,
           helpText("IP modelling approach"),
           
              textInput("Items",
                        label = "Number of items"),
           
             textInput("IPContainers",
                      label = "Number of containers"),
           
             textInput("IPContainerCap",
                      label = "Capacity of containers"),
           
             actionButton("runIP", "Run IP")
           
             ),
    column(4, offset = 1,
           
           helpText("Max Flow Approach"),
           
           textInput("Classes",
                     label = "Number of families"),
           
           textInput("Containers",
                     label = "Maxmimum number of containers"),
           
           textInput("ContainerCap",
                     label = "Capacity of containers"),
           
           textInput('Pallets',
                     label = "Number of pallets in various families (seperate quantities with a comma)"),
           textInput('maxk',
                     label = "Maxmimum similarity threshold"),
           
           actionButton("run", "Run")
           )
  ),

  
    
  
      mainPanel(
      h1("The network flows"),
      textOutput("selected_var"),
      plotOutput("Plot"),
      textOutput("comptime"),
      textOutput("IP"),
      DT::dataTableOutput("DS")
    )
  )


server <- function(input, output) {
  
  Classes <- eventReactive(input$run, {input$Classes})
  
  Containers <- eventReactive(input$run, {input$Containers})
  
  ContainerCap <- eventReactive(input$run, {input$ContainerCap})
  
  Pallets <- eventReactive(input$run, {input$Pallets})
  
  Items <- eventReactive(input$runIP, {input$Items})
  
  maximumk <- eventReactive(input$run, {input$mk})
  
  IPContainers <- eventReactive(input$runIP, {input$IPContainers})
  
  IPContainerCap <- eventReactive(input$runIP, {input$IPContainerCap})
  
  output$DS = DT::renderDataTable({
   
  })
  
  
  output$selected_var <- renderDataTable({ 
   
  library(igraph)
  library(Surrogate)
    x <- 1 
    
    numFamilies <- as.numeric(Classes())
    f <- as.numeric(unlist(strsplit(input$Pallets,",")))         # Vector of family sizes 
    familyCapacity <- max(f) # take note this value became 15 at 1010.
    numTables <- as.numeric(Containers())
    TableCapacity <- as.numeric(ContainerCap())  # Table capacity
    
    
    MaxTableNode <- numFamilies + numTables + 1  #Plus one is for the sink node
    familyNodes <- c()
    for(i in 1:numFamilies+1) { 
      familyNodes <- c(familyNodes,1,x+1,f[x])
      x<-x+1
    }
    familyNodes;
    famArcs <- t(matrix(familyNodes, nrow = 3))
    famArcs;
    
    k <- as.numeric(k())      #Diversity index
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
    E
    
    colnames(E) <- c("from", "to", "capacity")
    g1 <- graph_from_data_frame(as.data.frame(E))
    f1 <- max_flow(g1, source=V(g1)["1"], target=V(g1)[sink])$flow
    f2 <- max_flow(g1, source=V(g1)["1"], target=V(g1)[sink])$value
    
    t2 <- Sys.time()
    comptime <- t2-t1
    print(f2)
    
    
   # for(i in 1:as.numeric(Classes())){
    #  results <- f1[(as.numeric(Classes())+1+(i-1)*as.numeric(Containers())):(as.numeric(Classes())+i*as.numeric(Containers()))]
    #  print(results)
    #  i <- i + 1
     # }
 
})
  output$IP <- renderPrint({
    t1IP <- Sys.time()
    library(lpSolve)
    # Constraint 1
    n <- as.numeric(Items())
    m <- as.numeric(IPContainers())
    u <- as.numeric(IPContainerCap())
    
    # Matrix generation
    items <- as.numeric(Items())
    m.vect <- rep(0, items*items)
    for(i in 1:(items*items)-1){
      m.vect[i] <- as.integer(runif(1,1,10))
    }
    m.vect;
    x1<-1
    for(j in 1:items){
      m.vect[x1] <- 0
      x1 <- x1+items+1
    }
    m.vect;
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
    lp("max",obj, A, direction, rhs, int.vec = 1:(m*n*(n+1))/2)
    
    t2IP <- Sys.time()
    print(lp("max",obj, A, direction, rhs, int.vec = 1:(m*n*(n+1))/2))
    
    print(lp("max",obj, A, direction, rhs, int.vec = 1:(m*n*(n+1))/2)$solution)
    
  })
  output$Plot <- renderPlot({
    library(igraph)
    library(Surrogate)
    x <- 1 
    
    numFamilies <- as.numeric(Classes())
    familyCapacity <- max(f) # take note this value became 15 at 1010.
    f <- as.numeric(unlist(strsplit(input$Pallets,",")))         # Vector of family sizes 
    numTables <- as.numeric(Containers())
    TableCapacity <- as.numeric(ContainerCap())  # Table capacity
    
    
    MaxTableNode <- numFamilies + numTables + 1  #Plus one is for the sink node
    familyNodes <- c()
    for(i in 1:numFamilies+1) { 
      familyNodes <- c(familyNodes,1,x+1,f[x])
      x<-x+1
    }
    familyNodes;
    famArcs <- t(matrix(familyNodes, nrow = 3))
    famArcs;
    
    k <- max(f)      #Diversity index
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
    result <- max_flow(g1, source=V(g1)["1"], target=V(g1)[sink])$flow
    
    t2 <- Sys.time()
    comptime <- t2-t1
   plot(g1,layout=layout_as_tree)
   
  })
  output$comptime <- renderPrint({
    library(igraph)
    library(Surrogate)
    t1 <- Sys.time()
    x <- 1 
    
    numFamilies <- as.numeric(Classes())
    familyCapacity <- max(f) # take note this value became 15 at 1010.
    f <- as.numeric(unlist(strsplit(input$Pallets,",")))         # Vector of family sizes 
    numTables <- as.numeric(Containers())
    TableCapacity <- as.numeric(ContainerCap())  # Table capacity
    
    
    MaxTableNode <- numFamilies + numTables + 1  #Plus one is for the sink node
    familyNodes <- c()
    for(i in 1:numFamilies+1) { 
      familyNodes <- c(familyNodes,1,x+1,f[x])
      x<-x+1
    }
    familyNodes;
    famArcs <- t(matrix(familyNodes, nrow = 3))
    famArcs;
    
    k <- max(f)      #Diversity index
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
    result <- max_flow(g1, source=V(g1)["1"], target=V(g1)[sink])$flow
    
    t2 <- Sys.time()
    comptime <- t2-t1
    comptime
   
  })
}

shinyApp(ui = ui, server = server)