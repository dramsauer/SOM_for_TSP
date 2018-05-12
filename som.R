library(ggplot2)
library(readr)

#**********************************
# Functions
#**********************************
{
  
  neighborhood <- function(n_nodes, winner, k, sigma) {
    #'
    #' @param n_nodes amount of nodes used for SOM
    #' @param winner winning neuron in this iteration step
    #' @param k iteration parameter over the nodes
    #' @param sigma radius of the neighborhood function
    #' @return neighborhood-factor influencing the SOM learning rule
    nodes <- c(1:n_nodes)
    distance  <- sqrt((nodes[winner] - nodes[k]) ** 2)
    distance  <- min(distance, n_nodes-distance)
    h <- 1 / exp((distance ** 2) / (2 * sigma ** 2))
    return(h)
  }
  determine_winner <- function(n_nodes, coordinates, w, current_coordinates) {
    #' Function to determine the winning neuron 
    #' @param n_nodes amount of nodes used for SOM
    #' @param coordinates matrix containing coordinations of places in travling salesman problem
    #' @param w matrix which will get updated in SOM and later contain the winning neurons with their coordinates
    #' @param current_coordinates currently considered place in travling salesman problem
    #' @return the current winning neuron
    differences <- c()
    for (i in c(1:n_nodes)) {
      differences[i] <- 
        sqrt(
          (w[i, 1] - coordinates[current_coordinates, 1]) ** 2 + 
            (w[i, 2] - coordinates[current_coordinates, 2]) ** 2)
    }
    minimum <- min(differences)
    winner <- as.integer(which(differences == minimum))
    return(winner)
  }
  
  learning_rate <- function(rate, t) {
    #' @param rate learning rate value
    #' @param t iteration step
    #' @return the learning rate or step width
    return(0.01 + rate * t ** (-0.5))
  }
  sigma <- function(sigma, t) {
    #' Halbwertsbreite der neighborhood-function
    #' @param sigma radius value
    #' @param t iteration step
    #' @return radius at step t
    return(sigma / (t ** 0.5))
  }
  
  som <- function(n_nodes, coordinates, w, sigma, rate, iterations=8000, save_plots=FALSE){
    #' Main algorithm for a Self-Organizing-Map for the Travelling Salesman Problem
    #' @param n_nodes amount of nodes used for SOM
    #' @param coordinates matrix containing coordinations of places in travling salesman problem
    #' @param w matrix which will get updated in SOM and later contain the winning neurons with their coordinates
    #' @param sigma radius value
    #' @param rate learning rate value
    #' @param iterations value how often SOM algorithm shall be executed
    #' @param save_plots if set to TRUE, the algorithm will save a plot for each iteration
    #' @return matrix w which contains the winning neurons with their coordinates
    
    # Iteration for each learning cycle
    for (t in c(1:iterations)) {
      # Iteration over all places in travelling salesman problem
      for (j in c(1:nrow(coordinates))) {
        
        # Determine the current winning neuron
        winner <- determine_winner(n_nodes,coordinates, w, j)
        
        # Update neurons
        sig <- sigma(sigma, t)
        for (k in c(1:n_nodes)) {
          # Update equatation:  w + (learning rate * neighborhood * coordinates)
          w[k,] <-
            w[k,] + learning_rate(rate, t) * neighborhood(n_nodes, winner, k, sig) * (coordinates[j,] - w[k,])
        }
        
      }
      
      # Plot every graph
      if (save_plots) {
        if (t %% 15 == 0) {
          gname <- paste(c("N", n_nodes, "-", t, ".png"), collapse = "")
          coord.df  <- as.data.frame(coordinates)
          w.df      <- as.data.frame(w)
          colnames(coord.df) <- c("Grad O", "Grad N")
          colnames(w.df)     <- c("Grad O", "Grad N")
          
          g <- ggplot() +
            geom_point(data = coord.df, aes(`Grad O`, `Grad N`), color = "red", shape = 16, size = 4)+
            geom_point(data = w.df, aes(`Grad O`, `Grad N`), color = "green", shape = 4, size = 5)+
            geom_path(data = w.df, aes(`Grad O`, `Grad N`), color = "green")+
            ggtitle(paste(c("Iteration: ", t, " Nodes: ", n_nodes), collapse = ""))
          ggsave(filename = gname, plot = g, path = "Exercise II - ANN/TSP/plots")
        }
      }
    }
    return(w)
  }
  
  initialize_w_random <- function(n_nodes, grad_o_low, grad_o_high, grad_n_low, grad_n_high) {
    #' Function for initializing matrix w with random values between certain ranges
    #' w which will contain the winning neurons later    
    #' @param n_nodes amount of nodes used for SOM - matrix w will contain as many entries
    #' @param grad_o_low min value for Grad O  
    #' @param grad_o_high max value for Grad O  
    #' @param grad_n_low min value for Grad N
    #' @param grad_n_high max value for Grad N
    #' @return matrix w with start points
    w <- matrix(nrow = n_nodes, ncol = 2)
    
    for (i in c(1:n_nodes)) {
      w[i, 1] <- runif(1, grad_o_low, grad_o_high)
      w[i, 2] <- runif(1, grad_n_low, grad_n_high)
    }
    colnames(w, c("Grad O","Grad N"), do.NULL = FALSE)
    return(w)
  }
  initialize_w_from_df <- function(n_nodes, coordinates.df) {
    #' Function for initializing matrix w with samples from given dateframe
    #' w which will contain the winning neurons later
    #' @param n_nodes amount of nodes used for SOM - matrix w will contain as many entries
    #' @param coordinates.df data.frame containing coordinations of places in travling salesman problem
    #' @return matrix w with start points
    w <- as.matrix(dplyr::sample_n(coordinates.df, n_nodes, replace = TRUE)) * 1000
    return(w)
  }
  
  
  
  plot_data <- function(w,coordinates,title){
    coordinates.df  <- as.data.frame(coordinates)
    w.df      <- as.data.frame(w)
    colnames(coordinates.df) <- c("Grad O", "Grad N")
    colnames(w.df)     <- c("Grad O", "Grad N")
    
    plot_result <- ggplot() +
      geom_point(data = coordinates.df, aes(`Grad O`, `Grad N`), color = "red", shape = 16, size = 4 ) +
      geom_point(data = w.df, aes(`Grad O`, `Grad N`), color = "green", shape = 4, size = 5)+
      geom_path(data = w.df, aes(`Grad O`, `Grad N`), color = "green")+
      ggtitle(title)
    
    return(plot_result)
  }
  
}



#**********************************
# Data preperation
#**********************************
{
  
  # Read data with coordinates and create data.frame 
  beer_gardens <- read_csv("Exercise II - ANN/TSP/beer_gardens.csv")
  beer_gardens <- as.data.frame(beer_gardens)
  
  ggplot(beer_gardens, aes(x=`Grad O`, y=`Grad N`))+
    geom_point(aes(),color="red")+
    ggtitle("Coordinates of beer gardens in Regensburg")
  
  
  # Create matrices for easier calculation
  coordinates <- as.matrix(beer_gardens)
  coordinates <- coordinates * 1000
  
  # Learning rate value
  rate <- 0.2
  # Sigma - radius value for neighborhood function
  s <- 18
  # Amount of nodes
  n_nodes <- 18
  # Iterations for SOM executions
  iterations <- 8000
}



#**********************************
# Exercise
#**********************************

# SOM with start points from existing data
w1 <- initialize_w_from_df(n_nodes, beer_gardens)
w1_som <- som(n_nodes, coordinates, w1, s, rate, iterations, save_plots = TRUE)
#p1<-plot_data(w1_som,coordinates, "Start points from data.frame with 18 nodes - 8000 iterations")
#p1

# SOM with random start points from existing data
#w2 <- initialize_w_random(n_nodes, grad_o_low=12092, grad_o_high=12180, 
#                                   grad_n_low=49009, grad_n_high=49028)
#w2_som <- som(n_nodes, coordinates, w2, s, rate, iterations, save_plots = FALSE)
#p2<-plot_data(w2_som,coordinates, "Start points random with 25 nodes - 8000 iterations")
#p2

